"""Checkpointed Arm E GENCODE-v50 off-target transcriptome scan.

The reference is the deterministic union of annotated exonic intervals on each
strand.  Each FASTA record is in transcript (sense) orientation and retains a
direct genomic coordinate in its header.  This prevents transcript-isoform
double counting and makes every seed-and-extend hit consumable by H-off.

Scientific universe limitation: this reference covers only contiguous genomic
pre-mRNA exonic loci.  It deliberately does not contain exon-exon junction
sequences, so neither the scan nor H-off measures mature-transcript junction
off-targets.  Every terminal artifact records that limitation.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import signal
from typing import Iterator, Mapping

from feng.gpu.checkpoint import checkpoint_on_signal, load_checkpoint, save_checkpoint
from feng.scan.seed_extend import HybridizationHit, seed_and_extend
from feng.titration.ladder import build_ladder


SCRATCH_ROOT = Path(os.environ.get("SCRATCH", "/home/dashii/orcd/scratch"))
DEFAULT_GENOME_FASTA = SCRATCH_ROOT / "feng/gencode/v50/GRCh38.primary_assembly.genome.fa"
DEFAULT_GTF = SCRATCH_ROOT / "feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf"
DEFAULT_ARM_E_ROOT = SCRATCH_ROOT / "feng/arm_e"
PINNED_TRANSCRIPTOME_REFERENCE = DEFAULT_ARM_E_ROOT / "reference/gencode_v50_exonic_union.fa"
PINNED_TRANSCRIPTOME_CONFIG = DEFAULT_ARM_E_ROOT / "reference/gencode_v50_exonic_union.pin.json"
DEFAULT_MAX_EDIT = 5
MINIMUM_INFORMATIVE_LENGTH = 8
SCAN_ALGORITHM_VERSION = "armE-site-scan-v3"
SCAN_ARTIFACT_SCHEMA_VERSION = 2
REFERENCE_BUILD_ALGORITHM_VERSION = "armE-exonic-union-v2"
SCIENTIFIC_UNIVERSE = (
    "contiguous GENCODE-v50 pre-mRNA exonic loci only; exon-exon junction "
    "sequences and mature-transcript junction off-targets are excluded"
)
DEDUP_RULE = (
    "one best alignment per genomic site (chrom, start, end, strand); rank "
    "alternatives within a site by edit_distance, mismatches, target sequence, "
    "reference name"
)
EDIT_POLICY = (
    "requested max_edit capped at floor(query_length/5) (20% edit fraction) "
    "and 5; lengths below 8 are terminal-uninformative with burden=None"
)

SITE_COLUMNS = (
    "length", "aso_5to3", "chrom", "start", "end", "strand",
    "target_5to3", "reference_name", "edit_distance", "mismatches",
)


@dataclass(frozen=True, order=True)
class Exon:
    chrom: str
    start: int
    end: int
    strand: str

    @property
    def name(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}:{self.strand}"


def _sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _digest(payload: Mapping) -> str:
    canonical = json.dumps(dict(payload), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def effective_edit_policy(length: int, requested_max_edit: int) -> dict:
    """Return the preregistration-candidate, length-aware edit policy.

    A constant edit count changes its biological meaning across the ladder and
    becomes exhaustive when it reaches query length.  A fixed 20% fraction
    preserves comparable stringency, while the existing five-edit ceiling
    limits long-rung work.  Below 8 nt, even exact matches are too ubiquitous to
    identify useful loci; those rungs are explicitly incomplete, never zero.
    """
    if type(length) is not int or length < 1:
        raise ValueError("length must be a positive integer")
    if type(requested_max_edit) is not int or requested_max_edit < 0:
        raise ValueError("max_edit must be a nonnegative integer")
    informative = length >= MINIMUM_INFORMATIVE_LENGTH
    effective = min(requested_max_edit, DEFAULT_MAX_EDIT, length // 5) if informative else None
    return {
        "policy": EDIT_POLICY,
        "informative": informative,
        "minimum_informative_length": MINIMUM_INFORMATIVE_LENGTH,
        "requested_max_edit": requested_max_edit,
        "effective_max_edit": effective,
    }


def _validated_state(checkpoint: Path, input_digest: str) -> dict:
    state = load_checkpoint(checkpoint) or {}
    if state and state.get("input_digest") != input_digest:
        raise ValueError("checkpoint input/config digest mismatch; refusing stale resume")
    return state


def load_unique_exons(gtf: Path) -> list[Exon]:
    """Return sorted, merged exonic intervals, separately for each strand."""
    if not gtf.is_file() or gtf.stat().st_size == 0:
        raise FileNotFoundError(f"GENCODE GTF missing or empty: {gtf}")
    raw: set[Exon] = set()
    with gtf.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"malformed GTF line {line_number}: expected 9 fields")
            if fields[2] != "exon":
                continue
            strand = fields[6]
            if strand not in {"+", "-"}:
                raise ValueError(f"invalid exon strand at GTF line {line_number}: {strand!r}")
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"invalid exon interval at GTF line {line_number}")
            raw.add(Exon(fields[0], start, end, strand))

    merged: list[Exon] = []
    for exon in sorted(raw, key=lambda item: (item.chrom, item.strand, item.start, item.end)):
        if (merged and exon.chrom == merged[-1].chrom and exon.strand == merged[-1].strand
                and exon.start <= merged[-1].end + 1):
            previous = merged[-1]
            merged[-1] = Exon(previous.chrom, previous.start, max(previous.end, exon.end), previous.strand)
        else:
            merged.append(exon)
    if not merged:
        raise ValueError(f"GENCODE GTF contains no stranded exon records: {gtf}")
    return merged


def write_fai(fasta: Path) -> Path:
    """Write a minimal .fai for ordinary fixed-width FASTA (test/helper use)."""
    entries = []
    with fasta.open("rb") as handle:
        name = None
        length = offset = line_bases = line_width = 0
        while True:
            line = handle.readline()
            if not line:
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_width))
                break
            if line.startswith(b">"):
                if name is not None:
                    entries.append((name, length, offset, line_bases, line_width))
                name = line[1:].split()[0].decode("ascii")
                length = 0
                offset = handle.tell()
                line_bases = line_width = 0
            else:
                bases = len(line.rstrip(b"\r\n"))
                if line_bases == 0:
                    line_bases, line_width = bases, len(line)
                length += bases
    out = Path(str(fasta) + ".fai")
    out.write_text("".join(f"{n}\t{l}\t{o}\t{b}\t{w}\n" for n, l, o, b, w in entries))
    return out


def _load_fai(fasta: Path) -> dict[str, tuple[int, int, int, int]]:
    fai = Path(str(fasta) + ".fai")
    if not fasta.is_file() or fasta.stat().st_size == 0:
        raise FileNotFoundError(f"genome FASTA missing or empty: {fasta}")
    if not fai.is_file() or fai.stat().st_size == 0:
        raise FileNotFoundError(f"genome FASTA index missing or empty: {fai}")
    index = {}
    for line in fai.read_text().splitlines():
        fields = line.split("\t")
        index[fields[0]] = tuple(map(int, fields[1:5]))
    return index


def _fetch(handle, index: Mapping[str, tuple[int, int, int, int]], exon: Exon) -> str:
    try:
        contig_len, offset, line_bases, line_width = index[exon.chrom]
    except KeyError as exc:
        raise ValueError(f"GTF contig absent from FASTA index: {exon.chrom}") from exc
    if exon.end > contig_len:
        raise ValueError(f"exon {exon.name} exceeds contig length {contig_len}")
    start0, remaining = exon.start - 1, exon.end - exon.start + 1
    handle.seek(offset + (start0 // line_bases) * line_width + start0 % line_bases)
    pieces = []
    while remaining:
        chunk = handle.readline().strip()
        take = min(remaining, len(chunk))
        pieces.append(chunk[:take])
        remaining -= take
    seq = b"".join(pieces).decode("ascii").upper()
    if len(seq) != exon.end - exon.start + 1 or set(seq) - set("ACGTRYSWKMBDHVN"):
        raise ValueError(f"invalid sequence fetched for {exon.name}")
    if exon.strand == "-":
        seq = _reverse_complement(seq)
    return seq


def _reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
    return seq.translate(table)[::-1]


def _resume_file(path: Path, state: Mapping | None, header: str | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    offset = int((state or {}).get("output_offset", 0))
    mode = "r+b" if path.exists() else "w+b"
    handle = path.open(mode)
    if offset:
        if handle.seek(0, os.SEEK_END) < offset:
            raise ValueError(f"partial output shorter than checkpoint offset: {path}")
        handle.truncate(offset)
        handle.seek(offset)
    else:
        handle.seek(0)
        handle.truncate()
        if header is not None:
            handle.write(header.encode())
    return handle


def build_reference(*, gtf: Path, genome_fasta: Path, out: Path, pin: Path,
                    checkpoint: Path | None = None) -> dict:
    """Build and SHA-pin the deterministic coordinate-bearing exonic reference."""
    exons = load_unique_exons(gtf)
    index = _load_fai(genome_fasta)
    checkpoint = checkpoint or out.with_suffix(out.suffix + ".ckpt.json")
    partial = out.with_suffix(out.suffix + ".partial")
    build_inputs = {
        "algorithm_version": REFERENCE_BUILD_ALGORITHM_VERSION,
        "source_gtf_sha256": _sha256(gtf),
        "source_genome_fasta_sha256": _sha256(genome_fasta),
        "output_path": str(out.resolve()),
    }
    input_digest = _digest(build_inputs)
    state = _validated_state(checkpoint, input_digest)
    start_index = int(state.get("record_index", 0))
    output = _resume_file(partial, state)

    committed = [{
        "record_index": start_index,
        "output_offset": int(state.get("output_offset", output.tell())),
        "input_digest": input_digest,
        "inputs": build_inputs,
    }]

    def commit(record_index: int):
        output.flush()
        os.fsync(output.fileno())
        snapshot = {
            "record_index": record_index,
            "output_offset": output.tell(),
            "input_digest": input_digest,
            "inputs": build_inputs,
        }
        save_checkpoint(snapshot, checkpoint)
        committed[0] = snapshot

    def checkpoint_committed(_signum=None):
        save_checkpoint(dict(committed[0]), checkpoint)

    current = [start_index]
    try:
        with genome_fasta.open("rb") as genome, checkpoint_on_signal(
            checkpoint_committed, signals=(signal.SIGTERM,), requeue=False
        ):
            for i, exon in enumerate(exons[start_index:], start_index):
                seq = _fetch(genome, index, exon)
                output.write(f">{exon.name}\n{seq}\n".encode())
                current[0] = i + 1
                if current[0] % 1000 == 0:
                    commit(current[0])
            commit(len(exons))
    finally:
        output.close()
    os.replace(partial, out)
    payload = {
        "schema_version": 1,
        "reference_kind": "GENCODE-v50 stranded exonic union",
        "reference_path": str(out.resolve()),
        "reference_sha256": _sha256(out),
        "record_count": len(exons),
        "source_gtf": str(gtf.resolve()),
        "source_gtf_sha256": _sha256(gtf),
        "source_genome_fasta": str(genome_fasta.resolve()),
        "source_genome_fasta_sha256": _sha256(genome_fasta),
        "algorithm_version": REFERENCE_BUILD_ALGORITHM_VERSION,
        "scientific_universe": SCIENTIFIC_UNIVERSE,
    }
    save_checkpoint(payload, pin)
    return payload


def validate_pinned_reference(reference: Path, pin: Path) -> dict:
    if not reference.is_file() or reference.stat().st_size == 0:
        raise FileNotFoundError(f"pinned transcriptome reference missing or empty: {reference}")
    if not pin.is_file() or pin.stat().st_size == 0:
        raise FileNotFoundError(f"transcriptome pin config missing or empty: {pin}")
    payload = json.loads(pin.read_text())
    if Path(payload.get("reference_path", "")) != reference.resolve():
        raise ValueError("pinned transcriptome reference path mismatch")
    actual = _sha256(reference)
    if payload.get("reference_sha256") != actual:
        raise ValueError(f"pinned transcriptome SHA-256 mismatch: expected {payload.get('reference_sha256')}, got {actual}")
    return payload


def iter_reference(reference: Path) -> Iterator[tuple[Exon, str]]:
    name = None
    pieces: list[str] = []
    with reference.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield _parse_name(name), "".join(pieces)
                name, pieces = line[1:].split()[0], []
            else:
                if name is None:
                    raise ValueError("sequence before FASTA header")
                pieces.append(line.upper())
    if name is not None:
        yield _parse_name(name), "".join(pieces)


def _parse_name(name: str) -> Exon:
    chrom, interval, strand = name.rsplit(":", 2)
    start, end = map(int, interval.split("-"))
    return Exon(chrom, start, end, strand)


def _atomic_json(payload: Mapping, path: Path) -> None:
    tmp = path.with_suffix(path.suffix + ".partial")
    save_checkpoint(dict(payload), tmp)
    os.replace(tmp, path)


def _scan_inputs(rung: Mapping, reference_sha256: str, requested_max_edit: int,
                 policy: Mapping) -> dict:
    return {
        "algorithm_version": SCAN_ALGORITHM_VERSION,
        "reference_sha256": reference_sha256,
        "rung": {key: rung.get(key) for key in (
            "length", "aso_5to3", "chrom", "start", "end", "strand"
        )},
        "requested_max_edit": requested_max_edit,
        "effective_max_edit": policy["effective_max_edit"],
        "edit_policy": EDIT_POLICY,
        "dedup_rule": DEDUP_RULE,
        "scientific_universe": SCIENTIFIC_UNIVERSE,
    }


def validate_burden_artifact(burden_out: Path, manifest: Path) -> dict:
    """Authenticate a terminal scan result before a wrapper creates ``.done``."""
    if not burden_out.is_file() or not manifest.is_file():
        raise FileNotFoundError("scan burden or canonical manifest missing")
    payload = json.loads(burden_out.read_text())
    if payload.get("schema_version") != SCAN_ARTIFACT_SCHEMA_VERSION:
        raise ValueError("scan artifact schema mismatch")
    if payload.get("algorithm_version") != SCAN_ALGORITHM_VERSION:
        raise ValueError("scan artifact algorithm version mismatch")
    if payload.get("terminal") is not True:
        raise ValueError("scan artifact is not terminal")
    for key in ("input_digest", "reference_sha256"):
        value = payload.get(key)
        if not isinstance(value, str) or len(value) != 64 or any(
            char not in "0123456789abcdef" for char in value
        ):
            raise ValueError(f"scan artifact {key} missing or malformed")
    actual = _sha256(manifest)
    if payload.get("manifest_sha256") != actual:
        raise ValueError("scan manifest SHA-256 mismatch")
    try:
        rows = list(read_site_manifest(manifest))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("scan manifest schema invalid") from exc
    length, aso = payload.get("length"), payload.get("aso_5to3")
    if any(row["length"] != length or row["aso_5to3"] != aso for row in rows):
        raise ValueError("scan manifest rung/ASO mismatch")
    if payload.get("complete") is True:
        burden = payload.get("off_target_hits")
        if type(burden) is not int or burden < 0 or burden != len(rows):
            raise ValueError("scan artifact burden/manifest count mismatch")
    elif not (
        payload.get("complete") is False
        and payload.get("reason") == "length_below_minimum_informative_length"
        and payload.get("off_target_hits") is None
        and not rows
    ):
        raise ValueError("scan artifact has invalid incomplete state")
    return payload


def scan_rung(rung: Mapping, *, reference: Path, pin: Path, manifest: Path,
              burden_out: Path, checkpoint: Path, max_edit: int = DEFAULT_MAX_EDIT) -> dict:
    """Scan one ladder rung, resuming at a FASTA-record checkpoint."""
    pin_payload = validate_pinned_reference(reference, pin)
    length = int(rung["length"])
    policy = effective_edit_policy(length, max_edit)
    inputs = _scan_inputs(rung, pin_payload["reference_sha256"], max_edit, policy)
    input_digest = _digest(inputs)
    if manifest.is_file() and burden_out.is_file():
        prior = validate_burden_artifact(burden_out, manifest)
        if prior.get("input_digest") != input_digest:
            raise ValueError("completed scan input/config digest mismatch")
        return prior
    state = _validated_state(checkpoint, input_digest)
    start_index = int(state.get("record_index", 0))
    burden = int(state.get("off_target_hits", 0))
    partial = manifest.with_suffix(manifest.suffix + ".partial")
    header = "\t".join(SITE_COLUMNS) + "\n"
    output = _resume_file(partial, state, header)
    output.flush()
    os.fsync(output.fileno())
    committed = [{
        "record_index": start_index,
        "output_offset": int(state.get("output_offset", output.tell())),
        "off_target_hits": burden,
        "input_digest": input_digest,
        "inputs": inputs,
    }]

    def commit(record_index: int):
        output.flush()
        os.fsync(output.fileno())
        snapshot = {"record_index": record_index, "output_offset": output.tell(),
                    "off_target_hits": burden, "input_digest": input_digest,
                    "inputs": inputs}
        save_checkpoint(snapshot, checkpoint)
        committed[0] = snapshot

    def checkpoint_committed(_signum=None):
        # Never read output.tell(), burden, or current record here: a signal may
        # arrive after some rows from an exon were written but before its commit.
        save_checkpoint(dict(committed[0]), checkpoint)

    if not policy["informative"]:
        commit(start_index)
        output.close()
        os.replace(partial, manifest)
        payload = {
            "schema_version": SCAN_ARTIFACT_SCHEMA_VERSION,
            "algorithm_version": SCAN_ALGORITHM_VERSION,
            "length": length,
            "aso_5to3": rung["aso_5to3"],
            "off_target_hits": None,
            "off_target_burden": None,
            "complete": False,
            "terminal": True,
            "reason": "length_below_minimum_informative_length",
            **policy,
            "dedup_rule": DEDUP_RULE,
            "scientific_universe": SCIENTIFIC_UNIVERSE,
            "reference_sha256": pin_payload["reference_sha256"],
            "manifest_sha256": _sha256(manifest),
            "input_digest": input_digest,
            "state": "candidate",
        }
        _atomic_json(payload, burden_out)
        return payload

    query = _reverse_complement(str(rung["aso_5to3"]).upper().replace("U", "T"))
    try:
        current = [start_index]
        with checkpoint_on_signal(checkpoint_committed, signals=(signal.SIGTERM,), requeue=False):
            for i, (exon, sequence) in enumerate(iter_reference(reference)):
                if i < start_index:
                    continue
                candidates: dict[tuple[str, int, int, str], list[tuple]] = {}
                for hit in seed_and_extend(
                    query, sequence, max_edit=policy["effective_max_edit"],
                    reference_name=exon.name,
                ):
                    if exon.strand == "+":
                        start, end = exon.start + hit.start, exon.start + hit.end - 1
                    else:
                        start, end = exon.end - hit.end + 1, exon.end - hit.start
                    intended = (exon.strand == "+" and exon.chrom == str(rung["chrom"])
                                and start <= int(rung["end"]) and end >= int(rung["start"]))
                    if intended:
                        continue
                    target = sequence[hit.start:hit.end]
                    site = (exon.chrom, start, end, exon.strand)
                    candidates.setdefault(site, []).append((hit, target))
                for (chrom, start, end, strand), alternatives in sorted(candidates.items()):
                    hit, target = min(
                        alternatives,
                        key=lambda item: (
                            item[0].edit_distance, item[0].mismatches,
                            item[1], exon.name,
                        ),
                    )
                    row = (length, rung["aso_5to3"], exon.chrom, start, end,
                           exon.strand, target, exon.name, hit.edit_distance,
                           hit.mismatches)
                    output.write(("\t".join(map(str, row)) + "\n").encode())
                    burden += 1
                current[0] = i + 1
                if current[0] % 250 == 0:
                    commit(current[0])
            commit(current[0])
    finally:
        output.close()
    os.replace(partial, manifest)
    payload = {"schema_version": SCAN_ARTIFACT_SCHEMA_VERSION,
               "algorithm_version": SCAN_ALGORITHM_VERSION,
               "length": length, "aso_5to3": rung["aso_5to3"],
               "off_target_hits": burden, "off_target_burden": burden,
               "complete": True, "terminal": True, "reason": None,
               **policy, "dedup_rule": DEDUP_RULE,
               "scientific_universe": SCIENTIFIC_UNIVERSE,
               "reference_sha256": pin_payload["reference_sha256"],
               "manifest_sha256": _sha256(manifest),
               "input_digest": input_digest, "state": "candidate"}
    _atomic_json(payload, burden_out)
    return payload


def read_site_manifest(path: Path) -> Iterator[dict]:
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle, delimiter="\t"):
            row = dict(raw)
            for key in ("length", "start", "end", "edit_distance", "mismatches"):
                row[key] = int(row[key])
            yield row


def _main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("build-reference")
    build.add_argument("--gtf", type=Path, default=DEFAULT_GTF)
    build.add_argument("--genome-fasta", type=Path, default=DEFAULT_GENOME_FASTA)
    build.add_argument("--out", type=Path, default=PINNED_TRANSCRIPTOME_REFERENCE)
    build.add_argument("--pin", type=Path, default=PINNED_TRANSCRIPTOME_CONFIG)
    build.add_argument("--checkpoint", type=Path)
    validate = sub.add_parser("validate-reference")
    validate.add_argument("--reference", type=Path, default=PINNED_TRANSCRIPTOME_REFERENCE)
    validate.add_argument("--pin", type=Path, default=PINNED_TRANSCRIPTOME_CONFIG)
    scan = sub.add_parser("scan")
    scan.add_argument("--task", type=int, required=True)
    scan.add_argument("--reference", type=Path, default=PINNED_TRANSCRIPTOME_REFERENCE)
    scan.add_argument("--pin", type=Path, default=PINNED_TRANSCRIPTOME_CONFIG)
    scan.add_argument("--genome-fasta", type=Path, default=DEFAULT_GENOME_FASTA)
    scan.add_argument("--manifest", type=Path, required=True)
    scan.add_argument("--burden-out", type=Path, required=True)
    scan.add_argument("--checkpoint", type=Path, required=True)
    scan.add_argument("--max-edit", type=int, default=DEFAULT_MAX_EDIT)
    check = sub.add_parser("validate-burden")
    check.add_argument("--manifest", type=Path, required=True)
    check.add_argument("--burden-out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "build-reference":
        build_reference(gtf=args.gtf, genome_fasta=args.genome_fasta, out=args.out,
                        pin=args.pin, checkpoint=args.checkpoint)
    elif args.command == "validate-reference":
        validate_pinned_reference(args.reference, args.pin)
    elif args.command == "validate-burden":
        validate_burden_artifact(args.burden_out, args.manifest)
    else:
        rungs = build_ladder(genome_fasta=args.genome_fasta)
        if not 1 <= args.task <= len(rungs):
            parser.error(f"--task must be in 1..{len(rungs)}")
        scan_rung(rungs[args.task - 1], reference=args.reference, pin=args.pin,
                  manifest=args.manifest, burden_out=args.burden_out,
                  checkpoint=args.checkpoint, max_edit=args.max_edit)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
