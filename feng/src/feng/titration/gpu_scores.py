"""Arm E axis H -- Pangolin deep-splice scorer (candidate).

Two readouts, both ``candidate`` / pre-registration-eligible:

* **H-on (on-target potency):** replace each ASO footprint with ``N`` and ask
  Pangolin how much SMN2 exon-7 5' splice-site usage changes.  Pangolin encodes
  ``N`` as an all-zero one-hot column, so this readout is explicitly an
  **N-ablation sensitivity proxy**: it measures sensitivity to missing sequence
  recognition information, not a validated prediction of physical steric
  occlusion by a bound ASO.  A positive gain is the expected sign given the
  known ISS-N1 mechanism, but the proxy does not itself model that mechanism.

* **H-off (off-target splicing burden):** the same delta-splice prediction run
  at each near-complement off-target site and aggregated into a per-mer burden.
  It retains the same N-ablation-proxy interpretation.
  The per-site function is built + wired here; the genome-wide sweep is a
  separate GPU compute step (via ``scripts/arm_e_hoff.sbatch``).

Runtime split: this module imports cleanly under the CPU ``feng-env`` (torch and
Pangolin are imported only inside the worker, which runs under the GPU env
``~/orcd/pool/feng-gpu-env``).  Every scorer is gated on
:func:`feng.gpu.env.check_gpu_env`; when the GPU env / Pangolin is absent it
returns ``value=None`` with a reason and never fabricates a number.

The N-ablation proxy and the exon-7 readout position are deliberate MODELING
CHOICES flagged ``candidate`` -- they need review before any prereg seal.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Optional

from feng import source_auth

STATE = "candidate"
PANGOLIN_MASKING_METHOD = "N-ablation sensitivity proxy"

# --- SMN2 anchor (GRCh38 v50, + strand, 1-based; from the feng phase2 recon
#     README:167-177, the same verified frame the ladder is built on). ----------
CHROM = "chr5"
EXON7_LAST_BASE_1B = 70076574          # exon 7 ends here (+ strand)
EXON7_DONOR_1B = 70076575              # intron-7 first base == exon-7 5' donor boundary
ISS_N1_START_1B = 70076584             # nusinersen target window (== the 18-mer)
ISS_N1_END_1B = 70076601

PANGOLIN_CONTEXT = 5000                # flank Pangolin needs on each side
DEFAULT_DISTANCE = 50                  # Pangolin -d; output covers +/- d around center
# Pangolin labels 5' donor-site usage at the last exonic base, not at the first
# intronic base of the genomic splice boundary.
DEFAULT_READOUT_1B = EXON7_LAST_BASE_1B

PANGOLIN_WEIGHT_FILENAMES = tuple(
    f"final.{j}.{i}.3.v2" for i in (0, 2, 4, 6) for j in (1, 2, 3)
)

DEFAULT_GPU_PYTHON = Path(
    os.environ.get("FENG_GPU_PYTHON", os.path.expanduser("~/orcd/pool/feng-gpu-env/bin/python"))
)
DEFAULT_GENOME_FASTA = Path(
    os.environ.get(
        "FENG_ARMC_FASTA",
        os.path.expanduser("~/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa"),
    )
)


# ======================================================================
# Pure helpers -- no torch/pangolin; importable + testable under feng-env
# ======================================================================
def window_bounds_1b(center_1b: int, d: int = DEFAULT_DISTANCE, context: int = PANGOLIN_CONTEXT) -> tuple[int, int]:
    """1-based inclusive [start, end] of the context window fed to Pangolin.

    Length is ``2*(context+d)+1`` so that ``compute_score`` (which trims the
    model's fixed context loss) returns exactly ``2*d+1`` positions centered on
    ``center_1b``.
    """
    return center_1b - context - d, center_1b + context + d


def footprint_mask(
    ref_window: str,
    footprint_start_1b: int,
    footprint_end_1b: int,
    window_start_1b: int,
    mask_char: str = "N",
) -> str:
    """Return ``ref_window`` with the ASO footprint bases replaced by ``mask_char``.

    Length-preserving. ``N`` one-hot-encodes to an all-zero column in Pangolin
    (``IN_MAP[0]``), representing missing/unknown sequence information. This is
    an ablation proxy, not a validated physical-occlusion model. Raises if the
    footprint falls outside the window (fail-closed -- never silently mask
    nothing).
    """
    if footprint_end_1b < footprint_start_1b:
        raise ValueError("footprint end precedes start")
    off = footprint_start_1b - window_start_1b
    end = footprint_end_1b - window_start_1b
    if off < 0 or end >= len(ref_window):
        raise ValueError(
            f"footprint {footprint_start_1b}-{footprint_end_1b} outside window "
            f"start={window_start_1b} len={len(ref_window)}"
        )
    n = end - off + 1
    return ref_window[:off] + mask_char * n + ref_window[end + 1 :]


def readout_index(center_1b: int, readout_1b: int, d: int = DEFAULT_DISTANCE) -> int:
    """Array index into a ``2*d+1`` compute_score output for a genomic position."""
    idx = d + (readout_1b - center_1b)
    if idx < 0 or idx > 2 * d:
        raise ValueError(f"readout {readout_1b} outside +/-{d} of center {center_1b}")
    return idx


# ======================================================================
# GPU env gate
# ======================================================================
def _gpu_ready(gpu_python: Path) -> tuple[bool, Optional[str]]:
    """True iff the GPU env exists and reports torch+pangolin present."""
    try:
        from feng.gpu.env import check_gpu_env
    except Exception as exc:  # pragma: no cover - feng.gpu must be importable
        return False, f"feng.gpu.env import failed: {exc}"
    prefix = Path(gpu_python).resolve().parent.parent
    # Probe only the packages the H scorer needs (torch+pangolin; NOT the slow
    # openmm / absent spliceai), with a generous timeout: torch's cold import is
    # ~25s and exceeds check_gpu_env's 30s default under CPU load, false-negating.
    report = check_gpu_env(prefix=str(prefix), packages=("torch", "pangolin"), timeout_s=180.0)
    pkgs = report.get("packages", {})
    if not report.get("env_exists"):
        return False, f"GPU env absent at {prefix}"
    if not (pkgs.get("torch") and pkgs.get("pangolin")):
        return False, f"GPU env missing torch/pangolin: {pkgs} ({report.get('errors')})"
    return True, None


# ======================================================================
# Parent-side scorers (feng-env): fetch context, mask, dispatch to worker
# ======================================================================
def _fetch_window(center_1b: int, d: int, genome_fasta: Path) -> tuple[str, int]:
    """Fetch the + strand context window; return (sequence, window_start_1b)."""
    from feng.bench.scores import WindowSpec, _fetch_transcript_sequences  # noqa: WPS433
    from feng.titration.evaluate import _resolve_samtools  # reuse the proven resolver

    start_1b, end_1b = window_bounds_1b(center_1b, d)
    spec = WindowSpec(
        exon_key="armE_smn2_ctx",
        chrom=CHROM,
        start_0base=start_1b - 1,
        end=end_1b,
        strand="+",
    )
    seqs = _fetch_transcript_sequences([spec], Path(genome_fasta), samtools=_resolve_samtools(None))
    seq = seqs["armE_smn2_ctx"].upper()
    expected = end_1b - start_1b + 1
    if len(seq) != expected:
        raise ValueError(f"context fetch returned {len(seq)} nt, expected {expected}")
    return seq, start_1b


def _run_worker(jobs: list[dict], d: int, gpu_python: Path) -> list[dict]:
    """Dispatch scoring jobs to the GPU-env Pangolin worker. Never fabricates."""
    payload = json.dumps({"d": d, "jobs": jobs})
    env = dict(os.environ)
    src = str(Path(__file__).resolve().parents[2])  # feng/src
    env["PYTHONPATH"] = src
    env["PYTHONSAFEPATH"] = "1"
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    # Cap BLAS/OpenMP thread pools: torch's OpenBLAS otherwise spawns one thread
    # per core and blows past a shared node's RLIMIT_NPROC (768) under load,
    # making `import torch` itself fail (`blas_thread_init: RLIMIT_NPROC`). A
    # small cap fixes that AND is considerate on a shared login node. Override
    # via FENG_GPU_THREADS for a dedicated GPU-node run.
    _threads = os.environ.get("FENG_GPU_THREADS", "4")
    for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[_var] = _threads
    source_root = env.get("FENG_ARME_SOURCE_ROOT")
    if source_root:
        bootstrap = Path(source_root) / "feng" / "src" / "feng" / "source_auth.py"
        worker_command = [
            str(gpu_python), str(bootstrap), "run-module",
            "feng.titration.gpu_scores", "--", "_worker",
        ]
    elif env.get("FENG_REQUIRE_SOURCE_AUTH") == "1":
        raise RuntimeError("source authentication environment missing for Pangolin worker")
    else:
        worker_command = [str(gpu_python), "-m", "feng.titration.gpu_scores", "_worker"]
    proc = subprocess.run(
        worker_command,
        input=payload,
        capture_output=True,
        text=True,
        env=env,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"pangolin worker failed (rc={proc.returncode}): {proc.stderr[-500:]}")
    try:
        out = json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception as exc:
        raise RuntimeError(f"could not parse worker output: {exc}; stdout tail={proc.stdout[-300:]}")
    if "error" in out:
        raise RuntimeError(f"worker error: {out['error']}")
    return out["results"]


def pangolin_exon7_delta(
    rungs: list[dict],
    *,
    genome_fasta: Path = DEFAULT_GENOME_FASTA,
    gpu_python: Path = DEFAULT_GPU_PYTHON,
    d: int = DEFAULT_DISTANCE,
    readout_1b: int = DEFAULT_READOUT_1B,
) -> list[dict]:
    """H-on: predicted SMN2 exon-7 inclusion shift per ladder rung (candidate).

    Returns one dict per input rung: ``pangolin_exon7_delta`` (splice GAIN at the
    exon-7 donor under the N-ablation sensitivity proxy),
    ``pangolin_max_gain``/``pangolin_max_gain_offset`` (strongest induced gain
    in +/- d and where), ``masking_method='N-ablation sensitivity proxy'``,
    ``state='candidate'``, and ``pangolin_incomplete`` + ``reason`` when the GPU
    env is unavailable (value None -- never fabricated).
    """
    def _incomplete(reason: str) -> list[dict]:
        return [
            {
                "length": r.get("length"),
                "pangolin_exon7_delta": None,
                "pangolin_max_gain": None,
                "pangolin_max_gain_offset": None,
                "pangolin_incomplete": True,
                "reason": reason,
                "masking_method": PANGOLIN_MASKING_METHOD,
                "state": STATE,
            }
            for r in rungs
        ]

    # Cheap gate only: does the GPU env python exist? The AUTHORITATIVE
    # torch/pangolin availability check is the worker's own import (below) -- we
    # do NOT hard-gate on the load-sensitive check_gpu_env probe, which can
    # false-negative under CPU load. A genuine import failure surfaces as a
    # worker error -> _incomplete(reason), never a fabricated score.
    if not Path(gpu_python).exists():
        return _incomplete(f"GPU env python absent at {gpu_python}")

    center_1b = readout_1b
    idx = readout_index(center_1b, readout_1b, d)
    try:
        ref_window, win_start_1b = _fetch_window(center_1b, d, genome_fasta)
        jobs = []
        for r in rungs:
            alt = footprint_mask(ref_window, int(r["start"]), int(r["end"]), win_start_1b)
            jobs.append(
                {"id": int(r["length"]), "ref_seq": ref_window, "alt_seq": alt, "readout_idx": idx}
            )
        results = _run_worker(jobs, d, gpu_python)
    except (KeyError, RuntimeError, ValueError, OSError) as exc:
        return _incomplete(f"pangolin worker/fetch failed: {exc}")
    by_id = {res["id"]: res for res in results}

    out = []
    for r in rungs:
        res = by_id.get(int(r["length"]), {})
        out.append(
            {
                "length": r.get("length"),
                "pangolin_exon7_delta": res.get("delta_at_readout"),
                "pangolin_max_gain": res.get("max_gain"),
                "pangolin_max_gain_offset": res.get("max_gain_offset"),
                "pangolin_incomplete": res.get("delta_at_readout") is None,
                "reason": None if res.get("delta_at_readout") is not None else "no worker result",
                "masking_method": PANGOLIN_MASKING_METHOD,
                "state": STATE,
            }
        )
    return out


def pangolin_off_target_burden(
    aso_5to3: str,
    off_target_sites: list[dict],
    *,
    genome_fasta: Path = DEFAULT_GENOME_FASTA,
    gpu_python: Path = DEFAULT_GPU_PYTHON,
    d: int = DEFAULT_DISTANCE,
) -> dict:
    """H-off: aggregate delta-splice burden across an ASO's off-target sites (candidate).

    Each site dict needs ``chrom``, ``start`` and ``end`` (1-based inclusive),
    ``strand`` (target transcript strand), and ``target_5to3`` for the
    near-complement footprint. Burden = sum across sites of
    ``max(0, max_gain) + max(0, -max_loss)`` under the N-ablation sensitivity
    proxy. Returns ``value=None`` (incomplete) when the GPU env is absent or no
    sites were given -- never a fabricated zero.
    """
    def _incomplete(reason: str) -> dict:
        return {"pangolin_off_target_burden": None, "pangolin_off_target_incomplete": True,
                "reason": reason, "n_sites": len(off_target_sites),
                "masking_method": PANGOLIN_MASKING_METHOD, "state": STATE}

    if not off_target_sites:
        return _incomplete("no off-target sites supplied (not scanned)")  # != zero burden
    if not Path(gpu_python).exists():
        return _incomplete(f"GPU env python absent at {gpu_python}")

    try:
        aso = aso_5to3.upper().replace("U", "T")
        if not aso or set(aso) - set("ACGT"):
            raise ValueError("aso_5to3 must be a nonempty A/C/G/T sequence")
        from feng.bench.scores import reverse_complement
        from feng.scan.seed_extend import _alignment_cost
        hybridization_query = reverse_complement(aso)
        jobs = []
        for i, site in enumerate(off_target_sites):
            chrom = str(site.get("chrom", CHROM))
            strand = str(site["strand"])
            if strand not in {"+", "-"}:
                raise ValueError(f"off-target site {i} has invalid strand {strand!r}")
            center = (int(site["start"]) + int(site["end"])) // 2
            ref_window, win_start_1b = _fetch_window_generic(chrom, center, d, genome_fasta)
            lo = int(site["start"]) - win_start_1b
            hi = int(site["end"]) - win_start_1b + 1
            genomic_target = ref_window[lo:hi]
            oriented_target = (
                reverse_complement(genomic_target) if strand == "-" else genomic_target
            )
            declared_target = str(site["target_5to3"]).upper().replace("U", "T")
            if declared_target != oriented_target:
                raise ValueError(
                    f"off-target site {i} target_5to3 does not match reference on {strand} strand"
                )
            rc_distance = _alignment_cost(hybridization_query, declared_target)[0]
            direct_distance = _alignment_cost(aso, declared_target)[0]
            # Palindromic or near-symmetric targets can tie.  Their orientation
            # is established by the strand-aware reference check above, so only
            # a strictly closer direct-ASO match is invalid.
            if rc_distance > direct_distance:
                raise ValueError(
                    f"off-target site {i} has direct-ASO orientation "
                    f"(hybridization edit distance={rc_distance}, direct={direct_distance})"
                )
            alt = footprint_mask(ref_window, int(site["start"]), int(site["end"]), win_start_1b)
            jobs.append({
                "id": i,
                "ref_seq": ref_window,
                "alt_seq": alt,
                "readout_idx": d,
                "strand": strand,
                # Provenance metadata only; worker scoring uses ref_seq/alt_seq.
                "hybridization_query_5to3_provenance": hybridization_query,
            })
        results = _run_worker(jobs, d, gpu_python)
    except (KeyError, RuntimeError, ValueError, OSError) as exc:
        return _incomplete(f"pangolin worker/fetch failed: {exc}")
    by_id = {res.get("id"): res for res in results}
    if len(by_id) != len(results) or set(by_id) != set(range(len(off_target_sites))):
        return _incomplete("pangolin worker returned duplicate or missing site results")
    burden = 0.0
    for i in range(len(off_target_sites)):
        gain, loss = by_id[i].get("max_gain"), by_id[i].get("max_loss")
        if gain is None or loss is None:
            return _incomplete(f"pangolin worker returned missing value for site {i}")
        if not math.isfinite(float(gain)) or not math.isfinite(float(loss)):
            return _incomplete(f"pangolin worker returned non-finite value for site {i}")
        burden += max(0.0, float(gain)) + max(0.0, -float(loss))
    return {"pangolin_off_target_burden": burden, "pangolin_off_target_incomplete": False,
            "reason": None, "n_sites": len(off_target_sites),
            "masking_method": PANGOLIN_MASKING_METHOD, "state": STATE}


def _fetch_window_generic(chrom: str, center_1b: int, d: int, genome_fasta: Path) -> tuple[str, int]:
    from feng.bench.scores import WindowSpec, _fetch_transcript_sequences
    from feng.titration.evaluate import _resolve_samtools

    start_1b, end_1b = window_bounds_1b(center_1b, d)
    spec = WindowSpec(exon_key="armE_ot_ctx", chrom=chrom, start_0base=start_1b - 1, end=end_1b, strand="+")
    seqs = _fetch_transcript_sequences([spec], Path(genome_fasta), samtools=_resolve_samtools(None))
    return seqs["armE_ot_ctx"].upper(), start_1b


# ======================================================================
# Worker (runs under the GPU env python) -- torch + Pangolin here only
# ======================================================================
# NOTE: ``one_hot_encode`` and ``compute_score`` below are reproduced verbatim
# (bar the pyfastx/vcf-dependent variant plumbing) from Pangolin
# (github.com/tkzeng/Pangolin, pangolin/pangolin.py, GPL-3.0) so we can score a
# supplied ref/alt sequence pair directly without importing that module's
# vcf/gffutils/pyfastx top-level dependencies (absent in the GPU env). Only
# ``pangolin.model`` (numpy+torch) is imported.
def _worker_main() -> int:
    source_auth.verify_pangolin_package()
    # The worker's real import is the AUTHORITATIVE availability check -- a
    # genuine ImportError here (torch/pangolin actually missing) is reported as a
    # structured error, never a fabricated score. This is why the parent does not
    # hard-gate on the (load-sensitive) check_gpu_env probe.
    try:
        import numpy as np
        import torch
        import pangolin.model as pm
        from pangolin.model import Pangolin, L, W, AR
    except ImportError as exc:
        print(json.dumps({"error": f"gpu import failed: {type(exc).__name__}: {exc}"}))
        return 1

    IN_MAP = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

    def one_hot_encode(seq, strand):
        seq = seq.upper().replace("A", "1").replace("C", "2").replace("G", "3").replace("T", "4").replace("N", "0")
        if strand == "+":
            seq = np.asarray(list(map(int, list(seq))))
        else:
            seq = np.asarray(list(map(int, list(seq[::-1]))))
            seq = (5 - seq) % 5
        return IN_MAP[seq.astype("int8")]

    def compute_score(ref_seq, alt_seq, strand, d, models):
        ref_seq = one_hot_encode(ref_seq, strand).T
        ref_seq = torch.from_numpy(np.expand_dims(ref_seq, axis=0)).float()
        alt_seq = one_hot_encode(alt_seq, strand).T
        alt_seq = torch.from_numpy(np.expand_dims(alt_seq, axis=0)).float()
        if torch.cuda.is_available():
            ref_seq = ref_seq.to(torch.device("cuda"))
            alt_seq = alt_seq.to(torch.device("cuda"))
        pangolin = []
        for j in range(4):
            score = []
            for model in models[3 * j : 3 * j + 3]:
                with torch.no_grad():
                    ref = model(ref_seq)[0][[1, 4, 7, 10][j], :].cpu().numpy()
                    alt = model(alt_seq)[0][[1, 4, 7, 10][j], :].cpu().numpy()
                    if strand == "-":
                        ref = ref[::-1]
                        alt = alt[::-1]
                    l = 2 * d + 1
                    ndiff = np.abs(len(ref) - len(alt))
                    if len(ref) > len(alt):
                        alt = np.concatenate([alt[0 : l // 2 + 1], np.zeros(ndiff), alt[l // 2 + 1 :]])
                    elif len(ref) < len(alt):
                        alt = np.concatenate(
                            [alt[0 : l // 2], np.max(alt[l // 2 : l // 2 + ndiff + 1], keepdims=True), alt[l // 2 + ndiff + 1 :]]
                        )
                    score.append(alt - ref)
            pangolin.append(np.mean(score, axis=0))
        pangolin = np.array(pangolin)
        loss = pangolin[np.argmin(pangolin, axis=0), np.arange(pangolin.shape[1])]
        gain = pangolin[np.argmax(pangolin, axis=0), np.arange(pangolin.shape[1])]
        return loss, gain

    req = json.loads(sys.stdin.read())
    d = int(req["d"])

    mdir = os.path.join(os.path.dirname(pm.__file__), "models")
    models = []
    for weight_filename in PANGOLIN_WEIGHT_FILENAMES:
        model = Pangolin(L, W, AR)
        weights = torch.load(
            os.path.join(mdir, weight_filename), map_location=torch.device("cpu")
        )
        if torch.cuda.is_available():
            model.cuda()
        model.load_state_dict(weights)
        model.eval()
        models.append(model)

    results = []
    for job in req["jobs"]:
        strand = job.get("strand", "+")
        loss, gain = compute_score(job["ref_seq"], job["alt_seq"], strand, d, models)
        ridx = int(job["readout_idx"])
        max_i = int(gain.argmax())
        results.append(
            {
                "id": job["id"],
                "delta_at_readout": float(gain[ridx]),
                "loss_at_readout": float(loss[ridx]),
                "max_gain": float(gain[max_i]),
                "max_gain_offset": max_i - d,
                "max_loss": float(loss[int(loss.argmin())]),
            }
        )
    print(json.dumps({"results": results}))
    return 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_worker":
        raise SystemExit(_worker_main())
    raise SystemExit("usage: python -m feng.titration.gpu_scores _worker  (reads JSON on stdin)")
