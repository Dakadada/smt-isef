"""Deterministic J4 scale-out config builder (GRCh38 / GENCODE v50).

Turns the J2 loci manifest (cols: locus_id chrom start1 end1 strand, one row per
canonical protein-coding exon padded by +-EXON_PAD nt with start clipped at 1)
into per-site engine configs for conference.model_v2.gpu_perturbation_atlas.

Site convention (matches feng.titration.gpu_scores and the engine's Site docstring:
"position is exonic base"):
  donor    = LAST exonic base in TRANSCRIPT orientation  (SMN2 exon 7: 70076574 +)
  acceptor = FIRST exonic base in TRANSCRIPT orientation
  + strand exon [s, e]: acceptor = s, donor = e
  - strand exon [s, e]: acceptor = e, donor = s
Terminal sites that do not exist are never emitted: exon 1 has no acceptor and
exon N has no donor (N = max exon index of that transcript). Exon indices must be
exactly 1..N; a transcript with gaps is skipped whole, never guessed.

Interval rule (one absolute genomic tiling interval per config, per the engine):
  grid = [site - HALF_SPAN, site + HALF_SPAN], 1-based inclusive (2*HALF_SPAN+1 nt).
  Footprint starts are range(grid.start, grid.end - length + 2, stride), i.e. every
  footprint [start, start+length-1] lies FULLY inside the grid interval - identical
  arithmetic to gpu_perturbation_atlas.tasks().
  Window center = the site itself, distance = HALF_SPAN, so the readout index for the
  site is exact and the grid sits inside the Pangolin context window.
  Donor and acceptor of one exon are generally >HALF_SPAN apart, so they CANNOT share
  a tiling interval: each exon yields up to TWO configs (one per site), never one
  config with two grids.

Exon-bound recovery: exon_start = start1 + EXON_PAD, exon_end = end1 - EXON_PAD.
If start1 <= EXON_PAD the left pad was clipped at 1 and the true exon start is
unrecoverable -> explicit skipped manifest row, never guessed.

Evidence class for N-masked rows stays sequence_ablation_proxy (set by the engine);
this builder only prespecifies sites and grids, independently of any measured outcome.

Outputs (all atomic, no timestamps -> byte-identical across runs):
  shards/shard_%04d.jsonl     one canonical config JSON per line
  j4_configs_manifest.tsv     one row per (locus, site), emitted or skipped, + sha256
  j4_config_shards.tsv        per-shard file sha256 and counts
  j4_build_provenance.json    input sha256, parameters, totals

CLI: python -m conference.model_v2.j4_build_configs [--loci TSV] [--out-dir DIR]
     [--configs-per-shard N] [--half-span N] [--lengths 10,14,18] [--stride N]
     [--no-curated-smn2]
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re

from conference.model_v2 import gpu_perturbation_atlas as atlas
from conference.model_v2 import trans_intervention as ti
from feng.titration import gpu_scores as gs

VERSION = 1
EXON_PAD = 100                       # J2 manifest pad, see jobs/j2_build_loci.sbatch
DEFAULT_HALF_SPAN = 50
DEFAULT_LENGTHS = (10, 14, 18)
DEFAULT_STRIDE = 5
DEFAULT_LOCI = Path("conference/model_v2/manifests/j2_loci_phase1.tsv")
DEFAULT_OUT_DIR = Path("~/orcd/scratch/feng/model_v2/j4_configs")
DEFAULT_CONFIGS_PER_SHARD = 64
TISSUE_CHANNELS = len(atlas.TISSUES)

# Curated SMN2 exon-7 fine grid: stride 1 over the ISS-N1 region +- 30 nt.
CURATED_ID = "SMN2_exon7_issn1_fine"
CURATED_LENGTHS = (8, 10, 14, 18)
CURATED_FLANK = 30

LOCUS_RE = re.compile(r"^(?P<tx>[^\s:]+):ex(?P<idx>\d+)$")

MANIFEST_COLUMNS = (
    "config_id", "status", "reason", "locus_id", "transcript_id", "exon_index",
    "chrom", "strand", "exon_start", "exon_end", "site_position", "site_kind",
    "center", "distance", "grid_start", "grid_end", "stride", "lengths",
    "footprints", "rows", "inferences", "shard", "shard_file", "config_sha256",
)
SHARD_COLUMNS = ("shard", "shard_file", "configs", "footprints", "rows", "inferences", "file_sha256")


@dataclass(frozen=True)
class Locus:
    locus_id: str
    transcript_id: str
    exon_index: int
    chrom: str
    start1: int
    end1: int
    strand: str


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temp.write_text(text)
    os.replace(temp, path)


def read_loci(path) -> list[Locus]:
    """Parse the J2 manifest. Comment/blank lines are ignored; anything else malformed
    or duplicated is a hard error (fail closed - never a silently dropped locus)."""
    rows: list[Locus] = []
    seen: set[str] = set()
    for number, line in enumerate(Path(path).expanduser().read_text().splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(f"{path}:{number}: expected 5 tab-separated fields, got {len(fields)}")
        locus_id, chrom, start1, end1, strand = fields
        match = LOCUS_RE.match(locus_id)
        if not match:
            raise ValueError(f"{path}:{number}: locus_id {locus_id!r} is not '<transcript>:ex<N>'")
        if locus_id in seen:
            raise ValueError(f"{path}:{number}: duplicate locus_id {locus_id!r}")
        seen.add(locus_id)
        try:
            start, end = int(start1), int(end1)
        except ValueError as exc:
            raise ValueError(f"{path}:{number}: non-integer coordinates") from exc
        index = int(match.group("idx"))
        if index < 1:
            raise ValueError(f"{path}:{number}: exon index must be >= 1")
        if not 1 <= start <= end:
            raise ValueError(f"{path}:{number}: require 1 <= start1 <= end1")
        if strand not in ("+", "-"):
            raise ValueError(f"{path}:{number}: strand must be + or -")
        rows.append(Locus(locus_id, match.group("tx"), index, chrom, start, end, strand))
    if not rows:
        raise ValueError(f"{path}: no loci")
    return sorted(rows, key=lambda r: (r.transcript_id, r.exon_index))


def exon_bounds(locus: Locus) -> tuple[int | None, int | None, str | None]:
    """Recover the unpadded exon interval, or (None, None, reason) when ambiguous."""
    if locus.start1 <= EXON_PAD:
        return None, None, "padding_clipped_ambiguous"
    start, end = locus.start1 + EXON_PAD, locus.end1 - EXON_PAD
    if start > end:
        return None, None, "padding_trim_empty"
    return start, end, None


def transcript_site_plan(loci: list[Locus]) -> dict[str, tuple[int, str | None]]:
    """Per transcript: (max exon index, skip reason or None) from exon-index contiguity."""
    indices: dict[str, list[int]] = {}
    for locus in loci:
        indices.setdefault(locus.transcript_id, []).append(locus.exon_index)
    plan = {}
    for transcript, values in indices.items():
        top = max(values)
        reason = None if sorted(values) == list(range(1, top + 1)) else "exon_index_noncontiguous"
        plan[transcript] = (top, reason)
    return plan


def site_kinds(locus: Locus, exon_count: int) -> tuple[str, ...]:
    kinds = []
    if locus.exon_index > 1:
        kinds.append("acceptor")
    if locus.exon_index < exon_count:
        kinds.append("donor")
    return tuple(kinds)


def site_position(exon_start: int, exon_end: int, strand: str, kind: str) -> int:
    if kind == "donor":
        return exon_end if strand == "+" else exon_start
    if kind == "acceptor":
        return exon_start if strand == "+" else exon_end
    raise ValueError(f"unknown site kind {kind!r}")


def footprint_starts(grid_start: int, grid_end: int, length: int, stride: int) -> range:
    """Identical arithmetic to gpu_perturbation_atlas.tasks(): fully-contained footprints."""
    return range(grid_start, grid_end - length + 2, stride)


def footprint_count(grid_start: int, grid_end: int, lengths, stride: int) -> int:
    return sum(len(footprint_starts(grid_start, grid_end, n, stride)) for n in sorted(lengths))


def context_fits(chrom: str, center: int, distance: int) -> bool:
    """True iff the full Pangolin input window lies inside the GRCh38 contig."""
    length = atlas.GRCH38_PRIMARY_LENGTHS.get(chrom)
    if length is None:
        return False
    low, high = gs.window_bounds_1b(center, distance, gs.PANGOLIN_CONTEXT)
    return low >= 1 and high <= length


def config_identifier(locus_id: str, kind: str) -> str:
    return f"{locus_id.replace(':', '_')}__{kind}"


def make_config(config_id: str, locus: str, chrom: str, strand: str, center: int,
                distance: int, site: int, kind: str, grid_start: int, grid_end: int,
                lengths, stride: int) -> dict:
    return {
        "assembly": ti.ASSEMBLY,
        "annotation": ti.ANNOTATION,
        "builder_version": VERSION,
        "config_id": config_id,
        "windows": [{
            "locus": locus, "chrom": chrom, "strand": strand,
            "center": center, "distance": distance,
            "sites": [{"position": site, "kind": kind}],
        }],
        "grid": {"lengths": [int(n) for n in sorted(lengths)],
                 "start": grid_start, "end": grid_end, "stride": stride},
    }


def curated_smn2_config(stride: int = 1) -> dict:
    """SMN2 exon-7 donor, fine grid over ISS-N1 +- CURATED_FLANK nt."""
    center = gs.EXON7_LAST_BASE_1B
    return make_config(
        CURATED_ID, "SMN2_exon7", "chr5", "+", center, DEFAULT_HALF_SPAN,
        center, "donor",
        gs.ISS_N1_START_1B - CURATED_FLANK, gs.ISS_N1_END_1B + CURATED_FLANK,
        CURATED_LENGTHS, stride,
    )


def plan_configs(loci: list[Locus], *, half_span: int = DEFAULT_HALF_SPAN,
                 lengths=DEFAULT_LENGTHS, stride: int = DEFAULT_STRIDE,
                 curated: bool = True) -> tuple[list[dict], list[dict]]:
    """Return (emitted config records, skipped manifest records). Deterministic order."""
    plan = transcript_site_plan(loci)
    emitted: list[dict] = []
    skipped: list[dict] = []

    def skip(locus: Locus, kind: str | None, reason: str, exon=(None, None), site=None):
        skipped.append({
            "config_id": config_identifier(locus.locus_id, kind or "site"),
            "status": "skipped", "reason": reason, "locus_id": locus.locus_id,
            "transcript_id": locus.transcript_id, "exon_index": locus.exon_index,
            "chrom": locus.chrom, "strand": locus.strand,
            "exon_start": exon[0], "exon_end": exon[1],
            "site_position": site, "site_kind": kind,
        })

    for locus in loci:
        exon_count, transcript_reason = plan[locus.transcript_id]
        if transcript_reason:
            skip(locus, None, transcript_reason)
            continue
        start, end, reason = exon_bounds(locus)
        if reason:
            skip(locus, None, reason)
            continue
        kinds = site_kinds(locus, exon_count)
        if not kinds:
            skip(locus, None, "single_exon_transcript_no_splice_site", (start, end))
            continue
        for kind in kinds:
            site = site_position(start, end, locus.strand, kind)
            if not context_fits(locus.chrom, site, half_span):
                skip(locus, kind, "context_out_of_chrom", (start, end), site)
                continue
            grid_start, grid_end = site - half_span, site + half_span
            config_id = config_identifier(locus.locus_id, kind)
            config = make_config(config_id, locus.locus_id, locus.chrom, locus.strand,
                                 site, half_span, site, kind, grid_start, grid_end,
                                 lengths, stride)
            emitted.append({
                "config_id": config_id, "config": config, "locus_id": locus.locus_id,
                "transcript_id": locus.transcript_id, "exon_index": locus.exon_index,
                "exon_start": start, "exon_end": end,
                "site_position": site, "site_kind": kind,
            })
    if curated:
        config = curated_smn2_config()
        emitted.append({
            "config_id": CURATED_ID, "config": config, "locus_id": "SMN2_exon7",
            "transcript_id": "SMN2_curated", "exon_index": 7,
            "exon_start": None, "exon_end": None,
            "site_position": gs.EXON7_LAST_BASE_1B, "site_kind": "donor",
        })
    return emitted, skipped


def _record_accounting(config: dict) -> tuple[int, int, int]:
    grid = config["grid"]
    footprints = footprint_count(grid["start"], grid["end"], grid["lengths"], grid["stride"])
    sites = sum(len(w["sites"]) for w in config["windows"])
    rows = footprints * sites * TISSUE_CHANNELS
    inferences = footprints + len(config["windows"])
    return footprints, rows, inferences


def build(loci_path=DEFAULT_LOCI, out_dir=DEFAULT_OUT_DIR, *,
          configs_per_shard: int = DEFAULT_CONFIGS_PER_SHARD,
          half_span: int = DEFAULT_HALF_SPAN, lengths=DEFAULT_LENGTHS,
          stride: int = DEFAULT_STRIDE, curated: bool = True) -> dict:
    if configs_per_shard < 1:
        raise ValueError("configs_per_shard must be >= 1")
    if half_span < 1 or stride < 1:
        raise ValueError("half_span and stride must be >= 1")
    lengths = tuple(sorted({int(n) for n in lengths}))
    span = 2 * half_span + 1
    if not lengths or lengths[0] < 1 or lengths[-1] > span:
        raise ValueError(f"lengths must fit the {span} nt tiling interval")

    loci_path = Path(loci_path).expanduser()
    out = Path(out_dir).expanduser()
    loci = read_loci(loci_path)
    emitted, skipped = plan_configs(loci, half_span=half_span, lengths=lengths,
                                    stride=stride, curated=curated)

    shards: list[dict] = []
    manifest_rows: list[dict] = []
    for position, record in enumerate(emitted):
        shard = position // configs_per_shard
        if shard == len(shards):
            shards.append({"shard": shard, "name": f"shard_{shard:04d}.jsonl", "lines": [],
                           "configs": 0, "footprints": 0, "rows": 0, "inferences": 0})
        bucket = shards[shard]
        config = record["config"]
        text = canonical(config)
        footprints, rows, inferences = _record_accounting(config)
        bucket["lines"].append(text)
        bucket["configs"] += 1
        bucket["footprints"] += footprints
        bucket["rows"] += rows
        bucket["inferences"] += inferences
        window, grid = config["windows"][0], config["grid"]
        manifest_rows.append({
            "config_id": record["config_id"], "status": "emitted", "reason": None,
            "locus_id": record["locus_id"], "transcript_id": record["transcript_id"],
            "exon_index": record["exon_index"], "chrom": window["chrom"],
            "strand": window["strand"], "exon_start": record["exon_start"],
            "exon_end": record["exon_end"], "site_position": record["site_position"],
            "site_kind": record["site_kind"], "center": window["center"],
            "distance": window["distance"], "grid_start": grid["start"],
            "grid_end": grid["end"], "stride": grid["stride"],
            "lengths": ",".join(str(n) for n in grid["lengths"]),
            "footprints": footprints, "rows": rows, "inferences": inferences,
            "shard": shard, "shard_file": f"shards/{bucket['name']}",
            "config_sha256": sha256_text(text),
        })
    manifest_rows.extend(skipped)

    for bucket in shards:
        write_atomic(out / "shards" / bucket["name"], "".join(f"{line}\n" for line in bucket["lines"]))

    def cell(row, column):
        value = row.get(column)
        return "NA" if value is None else str(value)

    header = [
        "# J4 scale-out configs for conference.model_v2.gpu_perturbation_atlas",
        f"# builder_version={VERSION} assembly={ti.ASSEMBLY} annotation={ti.ANNOTATION}",
        f"# loci={loci_path} loci_sha256={sha256_file(loci_path)}",
        f"# exon_pad={EXON_PAD} half_span={half_span} lengths={','.join(map(str, lengths))} stride={stride}",
        "# interval rule: grid=[site-half_span, site+half_span] 1-based inclusive; footprint"
        " starts=range(grid_start, grid_end-length+2, stride); every footprint fully inside the grid",
        "# site convention: donor=last exonic base, acceptor=first exonic base, transcript orientation",
        "# config_sha256 = sha256 of the canonical config JSON text (the shard line, no newline)",
        "# evidence_class for N-masked rows stays sequence_ablation_proxy (set by the engine)",
        "#" + "\t".join(MANIFEST_COLUMNS),
    ]
    manifest_text = "\n".join(header + ["\t".join(cell(r, c) for c in MANIFEST_COLUMNS)
                                        for r in manifest_rows]) + "\n"
    write_atomic(out / "j4_configs_manifest.tsv", manifest_text)

    shard_rows = []
    for bucket in shards:
        path = out / "shards" / bucket["name"]
        shard_rows.append({"shard": bucket["shard"], "shard_file": f"shards/{bucket['name']}",
                           "configs": bucket["configs"], "footprints": bucket["footprints"],
                           "rows": bucket["rows"], "inferences": bucket["inferences"],
                           "file_sha256": sha256_file(path)})
    shard_text = "\n".join(["#" + "\t".join(SHARD_COLUMNS)]
                           + ["\t".join(str(row[c]) for c in SHARD_COLUMNS) for row in shard_rows]) + "\n"
    write_atomic(out / "j4_config_shards.tsv", shard_text)

    reasons: dict[str, int] = {}
    for row in skipped:
        reasons[row["reason"]] = reasons.get(row["reason"], 0) + 1
    summary = {
        "builder_version": VERSION,
        "assembly": ti.ASSEMBLY, "annotation": ti.ANNOTATION,
        "loci_path": str(loci_path), "loci_sha256": sha256_file(loci_path),
        "loci_rows": len(loci),
        "parameters": {"exon_pad": EXON_PAD, "half_span": half_span,
                       "lengths": list(lengths), "stride": stride,
                       "configs_per_shard": configs_per_shard, "curated_smn2": bool(curated),
                       "curated_lengths": list(CURATED_LENGTHS), "curated_flank": CURATED_FLANK,
                       "tissue_channels": TISSUE_CHANNELS,
                       "pangolin_context": gs.PANGOLIN_CONTEXT},
        "configs": len(emitted), "shards": len(shards),
        "skipped": len(skipped), "skipped_reasons": dict(sorted(reasons.items())),
        "total_footprints": sum(b["footprints"] for b in shards),
        "total_rows": sum(b["rows"] for b in shards),
        "total_inferences": sum(b["inferences"] for b in shards),
        "shard_files": [row["shard_file"] for row in shard_rows],
        "shard_sha256": {row["shard_file"]: row["file_sha256"] for row in shard_rows},
        "manifest_sha256": sha256_file(out / "j4_configs_manifest.tsv"),
        "shard_table_sha256": sha256_file(out / "j4_config_shards.tsv"),
    }
    write_atomic(out / "j4_build_provenance.json", canonical(summary) + "\n")
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--loci", type=Path, default=DEFAULT_LOCI)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--configs-per-shard", type=int, default=DEFAULT_CONFIGS_PER_SHARD)
    parser.add_argument("--half-span", type=int, default=DEFAULT_HALF_SPAN)
    parser.add_argument("--lengths", default=",".join(str(n) for n in DEFAULT_LENGTHS))
    parser.add_argument("--stride", type=int, default=DEFAULT_STRIDE)
    parser.add_argument("--no-curated-smn2", action="store_true")
    args = parser.parse_args(argv)
    summary = build(args.loci, args.out_dir, configs_per_shard=args.configs_per_shard,
                    half_span=args.half_span,
                    lengths=tuple(int(n) for n in args.lengths.split(",") if n.strip()),
                    stride=args.stride, curated=not args.no_curated_smn2)
    print(canonical({k: summary[k] for k in (
        "configs", "shards", "skipped", "skipped_reasons", "total_footprints",
        "total_rows", "total_inferences", "loci_sha256")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
