"""Typed, versioned, immutable representation of an ASO intervention on a pre-mRNA locus.

Scope: bookkeeping only (identity, genomic placement, splice-junction geometry,
evidence labels). Nothing here predicts potency, occupancy or PSI; those are
carried as explicit ``unidentifiable`` predictions with no value.

Coordinates: 1-based inclusive, GRCh38, GENCODE v50, genomic (+ strand) numbering.
``target_sequence`` is always transcript-sense DNA == reverse complement of the ASO.

CLI: ``python -m conference.model_v2.trans_intervention {validate,place,schema}``
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from dataclasses import dataclass, fields, is_dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

SCHEMA_VERSION = "1.0.0"
ASSEMBLY = "GRCh38"
ANNOTATION = "GENCODE v50"

EVIDENCE_CLASSES = (
    "observed",
    "derived",
    "proxy",
    "model_assumption",
    "unlabeled",
    "not_available",
)
STATUSES = ("placed", "needs_human_confirmation", "no_match", "mismatch")
ROLES = ("reference", "alternate")
STRANDS = ("+", "-")
IDENTIFIABILITY = ("identifiable", "unidentifiable")
REGIONS = ("exonic", "intronic", "exon_intron_spanning")
END_REGIONS = ("exonic", "intronic_upstream", "intronic_downstream")

DEFAULT_FASTA = "~/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa"
FASTA_ENV = "TRANS_INTERVENTION_FASTA"
MAX_SLICE_BASES = 100_000
SCHEMA_DIR = Path(__file__).resolve().parent / "schema"
SCHEMA_FILENAME = "trans_intervention.schema.json"

_COMP = str.maketrans("ACGTN", "TGCAN")


class SchemaError(ValueError):
    """Structural or invariant violation in an intervention record."""


class ParalogJoinError(ValueError):
    """Attempt to join placements that belong to different loci."""


class GenomeAccessError(RuntimeError):
    """Fasta / index unavailable or a slice request is out of bounds."""


def _req(cond: bool, msg: str) -> None:
    if not cond:
        raise SchemaError(msg)


def normalize_sequence(seq: str) -> str:
    """Upper-case DNA alphabet (U -> T); only A/C/G/T accepted."""
    out = str(seq).strip().upper().replace("U", "T")
    _req(len(out) > 0, "empty sequence")
    _req(set(out) <= set("ACGT"), f"non-ACGT characters in sequence: {seq!r}")
    return out


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMP)[::-1]


# --------------------------------------------------------------------------- genome access
class InMemoryGenome:
    """Tiny sequence source for tests: {chrom: plus-strand sequence}."""

    def __init__(self, sequences: Mapping[str, str]):
        self._seqs = {k: v.upper() for k, v in sequences.items()}

    def length(self, chrom: str) -> int:
        if chrom not in self._seqs:
            raise GenomeAccessError(f"unknown contig {chrom!r}")
        return len(self._seqs[chrom])

    def fetch(self, chrom: str, start: int, end: int) -> str:
        n = self.length(chrom)
        if not (1 <= start <= end <= n):
            raise GenomeAccessError(f"slice {chrom}:{start}-{end} outside 1..{n}")
        return self._seqs[chrom][start - 1 : end]


class FaiFastaReader:
    """Pure-python random access via the ``.fai`` index (read-only, small slices)."""

    def __init__(self, fasta_path: Optional[str] = None):
        raw = fasta_path or os.environ.get(FASTA_ENV) or DEFAULT_FASTA
        self.path = Path(os.path.expanduser(str(raw)))
        self.fai_path = Path(str(self.path) + ".fai")
        if not self.path.is_file():
            raise GenomeAccessError(f"genome fasta not found: {self.path}")
        if not self.fai_path.is_file():
            raise GenomeAccessError(f"fasta index not found: {self.fai_path}")
        self._index: Dict[str, Tuple[int, int, int, int]] = {}
        with open(self.fai_path, "r", encoding="ascii") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                self._index[parts[0]] = (int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]))

    def length(self, chrom: str) -> int:
        if chrom not in self._index:
            raise GenomeAccessError(f"unknown contig {chrom!r}")
        return self._index[chrom][0]

    def _offset(self, chrom: str, pos0: int) -> int:
        _, off, line_bases, line_width = self._index[chrom]
        return off + (pos0 // line_bases) * line_width + (pos0 % line_bases)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        n = self.length(chrom)
        if not (1 <= start <= end <= n):
            raise GenomeAccessError(f"slice {chrom}:{start}-{end} outside 1..{n}")
        if end - start + 1 > MAX_SLICE_BASES:
            raise GenomeAccessError(f"slice larger than {MAX_SLICE_BASES} bases refused")
        first = self._offset(chrom, start - 1)
        last = self._offset(chrom, end - 1)
        with open(self.path, "rb") as fh:
            fh.seek(first)
            raw = fh.read(last - first + 1)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(seq) != end - start + 1:
            raise GenomeAccessError(f"short read for {chrom}:{start}-{end}")
        return seq


# --------------------------------------------------------------------------- data model
@dataclass(frozen=True)
class Modification:
    position: int  # 1-based from the ASO 5' end
    sugar: str = "unknown"
    backbone: str = "unknown"
    base: str = "unknown"

    def __post_init__(self) -> None:
        _req(isinstance(self.position, int) and self.position >= 1, "modification position must be >= 1")


@dataclass(frozen=True)
class Chemistry:
    family: str = "unknown"
    modifications: Tuple[Modification, ...] = ()
    evidence: str = "not_available"

    def __post_init__(self) -> None:
        _req(self.evidence in EVIDENCE_CLASSES, f"bad evidence class {self.evidence!r}")
        _req(bool(self.family), "chemistry family must be non-empty ('unknown' allowed)")
        if self.family == "unknown":
            _req(self.evidence in ("not_available", "unlabeled"), "unknown chemistry must be not_available/unlabeled")
            _req(len(self.modifications) == 0, "unknown chemistry cannot list modifications")
        else:
            _req(self.evidence != "not_available", "known chemistry needs an evidence class")
        pos = [m.position for m in self.modifications]
        _req(len(pos) == len(set(pos)), "duplicate modification positions")
        _req(pos == sorted(pos), "modifications must be sorted by position")


@dataclass(frozen=True)
class ASO:
    aso_id: str
    sequence_5to3: str
    chemistry: Chemistry = Chemistry()
    requires_human_confirmation: bool = False
    sequence_evidence: str = "unlabeled"
    note: str = ""

    def __post_init__(self) -> None:
        _req(bool(self.aso_id), "aso_id required")
        object.__setattr__(self, "sequence_5to3", normalize_sequence(self.sequence_5to3))
        _req(self.sequence_evidence in EVIDENCE_CLASSES, f"bad evidence class {self.sequence_evidence!r}")
        for m in self.chemistry.modifications:
            _req(m.position <= self.length, f"modification position {m.position} beyond ASO length {self.length}")

    @property
    def length(self) -> int:
        return len(self.sequence_5to3)

    @property
    def target_sequence(self) -> str:
        """Transcript-sense DNA the ASO base-pairs with."""
        return reverse_complement(self.sequence_5to3)


@dataclass(frozen=True)
class Locus:
    """One exon-centred pre-mRNA locus. Paralogs are distinct loci, never merged."""

    locus_id: str
    gene: str
    paralog_group: str
    role: str
    chrom: str
    strand: str
    search_start: int
    search_end: int
    exon_label: str
    exon_start: int  # genomic, 1-based inclusive
    exon_end: int
    coords_evidence: str = "derived"
    assembly: str = ASSEMBLY
    annotation: str = ANNOTATION
    note: str = ""

    def __post_init__(self) -> None:
        _req(bool(self.locus_id), "locus_id required")
        _req(self.role in ROLES, f"bad role {self.role!r}")
        _req(self.strand in STRANDS, f"bad strand {self.strand!r}")
        _req(self.coords_evidence in EVIDENCE_CLASSES, f"bad evidence class {self.coords_evidence!r}")
        _req(1 <= self.search_start <= self.search_end, "bad search window")
        _req(self.search_start <= self.exon_start <= self.exon_end <= self.search_end, "exon outside search window")

    @property
    def sign(self) -> int:
        return 1 if self.strand == "+" else -1

    @property
    def donor_last_exonic_base(self) -> int:
        return self.exon_end if self.strand == "+" else self.exon_start

    @property
    def donor_first_intronic_base(self) -> int:
        return self.donor_last_exonic_base + self.sign

    @property
    def acceptor_first_exonic_base(self) -> int:
        return self.exon_start if self.strand == "+" else self.exon_end

    @property
    def acceptor_last_intronic_base(self) -> int:
        return self.acceptor_first_exonic_base - self.sign


@dataclass(frozen=True)
class SpliceGeometry:
    """Signed, strand-aware offsets (transcript orientation; downstream is positive).

    Both donor conventions are carried: origin at the last exonic base (first
    intronic base == +1, HGVS-like) and origin at the first intronic base
    (first intronic base == 0). Likewise for the acceptor: origin at the first
    exonic base, and origin at the last intronic base.
    """

    donor_last_exonic_base: int
    donor_first_intronic_base: int
    acceptor_first_exonic_base: int
    acceptor_last_intronic_base: int
    target_5p_pos: int
    target_3p_pos: int
    donor_dist_5p_from_last_exonic: int
    donor_dist_3p_from_last_exonic: int
    donor_dist_5p_from_first_intronic: int
    donor_dist_3p_from_first_intronic: int
    acceptor_dist_5p_from_first_exonic: int
    acceptor_dist_3p_from_first_exonic: int
    acceptor_dist_5p_from_last_intronic: int
    acceptor_dist_3p_from_last_intronic: int
    region_5p: str
    region_3p: str
    region: str
    evidence: str = "derived"

    def __post_init__(self) -> None:
        _req(self.evidence in EVIDENCE_CLASSES, f"bad evidence class {self.evidence!r}")
        _req(self.region in REGIONS, f"bad region {self.region!r}")
        _req(self.region_5p in END_REGIONS and self.region_3p in END_REGIONS, "bad end region")


@dataclass(frozen=True)
class Placement:
    site_id: str
    locus_id: str
    assembly: str
    annotation: str
    chrom: str
    start: int
    end: int
    strand: str
    target_sequence: str
    geometry: SpliceGeometry
    evidence: str = "derived"

    def __post_init__(self) -> None:
        _req(self.strand in STRANDS, f"bad strand {self.strand!r}")
        _req(self.evidence in EVIDENCE_CLASSES, f"bad evidence class {self.evidence!r}")
        _req(1 <= self.start <= self.end, "bad placement interval")
        _req(self.end - self.start + 1 == len(self.target_sequence), "interval length != target length")
        expect = target_site_id(self.locus_id, self.assembly, self.chrom, self.start, self.end, self.strand, self.target_sequence)
        _req(self.site_id == expect, "site_id does not match its inputs")


@dataclass(frozen=True)
class PlacementResult:
    locus_id: str
    role: str
    status: str
    reason: str
    window_start: int
    window_end: int
    placement: Optional[Placement] = None
    candidates: Tuple[Placement, ...] = ()

    def __post_init__(self) -> None:
        _req(self.role in ROLES, f"bad role {self.role!r}")
        _req(self.status in STATUSES, f"bad status {self.status!r}")
        if self.status == "placed":
            _req(self.placement is not None and not self.candidates, "placed requires exactly one placement")
        else:
            _req(self.placement is None, f"status {self.status} must not carry a chosen placement")
        if self.status == "needs_human_confirmation":
            _req(len(self.candidates) >= 1, "needs_human_confirmation must enumerate candidates")
        if self.status in ("no_match", "mismatch"):
            _req(not self.candidates, f"status {self.status} must not carry candidates")
        for p in ((self.placement,) if self.placement else ()) + tuple(self.candidates):
            _req(p.locus_id == self.locus_id, "placement locus_id differs from result locus_id (paralog join)")
        ids = [p.site_id for p in self.candidates]
        _req(len(ids) == len(set(ids)), "duplicate candidate site_id")


@dataclass(frozen=True)
class Prediction:
    """A quantity the record refuses to assign by fiat."""

    name: str
    value: Optional[float] = None
    evidence: str = "not_available"
    identifiability: str = "unidentifiable"
    note: str = ""

    def __post_init__(self) -> None:
        _req(bool(self.name), "prediction name required")
        _req(self.evidence in EVIDENCE_CLASSES, f"bad evidence class {self.evidence!r}")
        _req(self.identifiability in IDENTIFIABILITY, f"bad identifiability {self.identifiability!r}")
        if self.identifiability == "unidentifiable":
            _req(self.value is None, f"{self.name}: unidentifiable prediction cannot carry a value")
            _req(self.evidence == "not_available", f"{self.name}: unidentifiable prediction must be not_available")
        if self.value is None:
            _req(self.evidence == "not_available", f"{self.name}: missing value must be not_available")
        else:
            _req(self.evidence != "not_available", f"{self.name}: value present but evidence not_available")
        if self.name == "observed_psi" and self.value is not None:
            _req(self.evidence == "observed", "observed_psi may only be set from an observation")
        if self.name == "molecular_occupancy":
            _req(self.value is None, "molecular_occupancy is unidentifiable here; never assigned")


def default_predictions() -> Tuple[Prediction, ...]:
    return (
        Prediction("molecular_occupancy", note="not identifiable from sequence/placement; dG/occupancy != potency"),
        Prediction("observed_psi", note="requires a measurement; never assigned by fiat"),
    )


@dataclass(frozen=True)
class Intervention:
    schema_version: str
    aso: ASO
    reference: PlacementResult
    alternates: Tuple[PlacementResult, ...] = ()
    predictions: Tuple[Prediction, ...] = ()

    def __post_init__(self) -> None:
        errs = validate(self)
        _req(not errs, "; ".join(errs))


# --------------------------------------------------------------------------- placement logic
def target_site_id(locus_id: str, assembly: str, chrom: str, start: int, end: int, strand: str, target_sequence: str) -> str:
    """Deterministic ID; includes locus_id so identical paralog sequences never collide."""
    key = "|".join([locus_id, assembly, chrom, str(start), str(end), strand, target_sequence])
    return "ts_" + hashlib.sha256(key.encode("ascii")).hexdigest()[:20]


def _end_region(locus: Locus, pos: int) -> str:
    if locus.exon_start <= pos <= locus.exon_end:
        return "exonic"
    return "intronic_downstream" if (pos - locus.donor_last_exonic_base) * locus.sign > 0 else "intronic_upstream"


def splice_geometry(locus: Locus, start: int, end: int) -> SpliceGeometry:
    s = locus.sign
    p5, p3 = (start, end) if s == 1 else (end, start)
    d5 = (p5 - locus.donor_last_exonic_base) * s
    d3 = (p3 - locus.donor_last_exonic_base) * s
    a5 = (p5 - locus.acceptor_first_exonic_base) * s
    a3 = (p3 - locus.acceptor_first_exonic_base) * s
    r5, r3 = _end_region(locus, p5), _end_region(locus, p3)
    if r5 == r3 == "exonic":
        region = "exonic"
    elif r5 == r3:
        region = "intronic"
    else:
        region = "exon_intron_spanning"  # also covers a target that swallows the whole exon
    return SpliceGeometry(
        donor_last_exonic_base=locus.donor_last_exonic_base,
        donor_first_intronic_base=locus.donor_first_intronic_base,
        acceptor_first_exonic_base=locus.acceptor_first_exonic_base,
        acceptor_last_intronic_base=locus.acceptor_last_intronic_base,
        target_5p_pos=p5,
        target_3p_pos=p3,
        donor_dist_5p_from_last_exonic=d5,
        donor_dist_3p_from_last_exonic=d3,
        donor_dist_5p_from_first_intronic=d5 - 1,
        donor_dist_3p_from_first_intronic=d3 - 1,
        acceptor_dist_5p_from_first_exonic=a5,
        acceptor_dist_3p_from_first_exonic=a3,
        acceptor_dist_5p_from_last_intronic=a5 + 1,
        acceptor_dist_3p_from_last_intronic=a3 + 1,
        region_5p=r5,
        region_3p=r3,
        region=region,
    )


def _make_placement(locus: Locus, start: int, end: int, target: str) -> Placement:
    return Placement(
        site_id=target_site_id(locus.locus_id, locus.assembly, locus.chrom, start, end, locus.strand, target),
        locus_id=locus.locus_id,
        assembly=locus.assembly,
        annotation=locus.annotation,
        chrom=locus.chrom,
        start=start,
        end=end,
        strand=locus.strand,
        target_sequence=target,
        geometry=splice_geometry(locus, start, end),
    )


def _transcript_sense(locus: Locus, genomic_plus: str) -> str:
    return genomic_plus if locus.strand == "+" else reverse_complement(genomic_plus)


def place(aso: ASO, locus: Locus, genome: Any, window: Optional[Tuple[int, int]] = None) -> PlacementResult:
    """Exact reverse-complement placement of ``aso`` inside one locus window.

    Never picks among several matches and never resolves an ASO flagged
    ``requires_human_confirmation``: those return ``needs_human_confirmation``
    with every candidate enumerated. No fuzzy matching.
    """
    ws, we = window if window is not None else (locus.search_start, locus.search_end)
    if not (locus.search_start <= ws <= we <= locus.search_end):
        raise ValueError(f"window {ws}-{we} outside locus {locus.locus_id} search bounds")
    base = dict(locus_id=locus.locus_id, role=locus.role, window_start=ws, window_end=we)
    if aso.length > we - ws + 1:
        return PlacementResult(status="no_match", reason="aso_longer_than_window", **base)
    tx = _transcript_sense(locus, genome.fetch(locus.chrom, ws, we))
    target = aso.target_sequence
    hits: List[Placement] = []
    i = tx.find(target)
    while i >= 0:
        start = ws + i if locus.strand == "+" else we - i - aso.length + 1
        hits.append(_make_placement(locus, start, start + aso.length - 1, target))
        i = tx.find(target, i + 1)
    hits.sort(key=lambda p: (p.start, p.end))
    if not hits:
        return PlacementResult(status="no_match", reason="no_exact_match", **base)
    if aso.requires_human_confirmation:
        return PlacementResult(status="needs_human_confirmation", reason="aso_flagged_requires_human_confirmation", candidates=tuple(hits), **base)
    if len(hits) > 1:
        return PlacementResult(status="needs_human_confirmation", reason="multiple_exact_matches", candidates=tuple(hits), **base)
    return PlacementResult(status="placed", reason="unique_exact_match", placement=hits[0], **base)


def place_at(aso: ASO, locus: Locus, start: int, end: int, genome: Any) -> PlacementResult:
    """Verify a claimed placement against the genome (off-by-one guard)."""
    if not (locus.search_start <= start <= end <= locus.search_end):
        raise ValueError(f"interval {start}-{end} outside locus {locus.locus_id} search bounds")
    base = dict(locus_id=locus.locus_id, role=locus.role, window_start=start, window_end=end)
    if end - start + 1 != aso.length:
        return PlacementResult(status="mismatch", reason="length_mismatch", **base)
    tx = _transcript_sense(locus, genome.fetch(locus.chrom, start, end))
    target = aso.target_sequence
    n_mm = sum(1 for a, b in zip(tx, target) if a != b)
    if n_mm:
        return PlacementResult(status="mismatch", reason=f"sequence_mismatch:n={n_mm}", **base)
    hit = _make_placement(locus, start, end, target)
    if aso.requires_human_confirmation:
        return PlacementResult(status="needs_human_confirmation", reason="aso_flagged_requires_human_confirmation", candidates=(hit,), **base)
    return PlacementResult(status="placed", reason="verified_exact_match", placement=hit, **base)


def join_placements(left: Placement, right: Placement) -> bool:
    """True iff both are the same target site. Cross-locus joins raise, even when
    the sequences are identical (SMN1 vs SMN2)."""
    if left.locus_id != right.locus_id:
        raise ParalogJoinError(f"refusing to join across loci: {left.locus_id} vs {right.locus_id}")
    return left.site_id == right.site_id


def build_intervention(
    aso: ASO,
    reference_locus: Locus,
    genome: Any,
    alternate_loci: Sequence[Locus] = (),
    window: Optional[Tuple[int, int]] = None,
) -> Intervention:
    _req(reference_locus.role == "reference", "reference locus must have role 'reference'")
    for loc in alternate_loci:
        _req(loc.role == "alternate", f"{loc.locus_id}: alternate locus must have role 'alternate'")
    return Intervention(
        schema_version=SCHEMA_VERSION,
        aso=aso,
        reference=place(aso, reference_locus, genome, window=window),
        alternates=tuple(place(aso, loc, genome) for loc in alternate_loci),
        predictions=default_predictions(),
    )


def validate(iv: Intervention) -> List[str]:
    """Cross-field invariants. Returns a list of problems (empty == valid)."""
    errs: List[str] = []
    if iv.schema_version != SCHEMA_VERSION:
        errs.append(f"schema_version {iv.schema_version!r} != {SCHEMA_VERSION!r}")
    if iv.reference.role != "reference":
        errs.append("reference placement result must have role 'reference'")
    seen_loci = {iv.reference.locus_id}
    for alt in iv.alternates:
        if alt.role != "alternate":
            errs.append(f"{alt.locus_id}: alternate placement result must have role 'alternate'")
        if alt.locus_id in seen_loci:
            errs.append(f"{alt.locus_id}: duplicate locus_id across reference/alternates")
        seen_loci.add(alt.locus_id)
    site_ids: List[str] = []
    for res in (iv.reference,) + tuple(iv.alternates):
        for p in ((res.placement,) if res.placement else ()) + tuple(res.candidates):
            site_ids.append(p.site_id)
            if p.target_sequence != iv.aso.target_sequence:
                errs.append(f"{p.site_id}: target_sequence is not the reverse complement of the ASO")
            if p.assembly != ASSEMBLY or p.annotation != ANNOTATION:
                errs.append(f"{p.site_id}: assembly/annotation must be {ASSEMBLY}/{ANNOTATION}")
        if iv.aso.requires_human_confirmation and res.status == "placed":
            errs.append(f"{res.locus_id}: ASO requires human confirmation but status is 'placed'")
    if len(site_ids) != len(set(site_ids)):
        errs.append("site_id shared between placements (cross-locus join)")
    names = [p.name for p in iv.predictions]
    if len(names) != len(set(names)):
        errs.append("duplicate prediction names")
    for required in ("molecular_occupancy", "observed_psi"):
        if required not in names:
            errs.append(f"missing explicit prediction entry: {required}")
    return errs


# --------------------------------------------------------------------------- serialization
_NESTED: Dict[type, Dict[str, Tuple[str, type]]] = {
    Chemistry: {"modifications": ("list", Modification)},
    ASO: {"chemistry": ("one", Chemistry)},
    Placement: {"geometry": ("one", SpliceGeometry)},
    PlacementResult: {"placement": ("optional", Placement), "candidates": ("list", Placement)},
    Intervention: {
        "aso": ("one", ASO),
        "reference": ("one", PlacementResult),
        "alternates": ("list", PlacementResult),
        "predictions": ("list", Prediction),
    },
}
_ENUMS: Dict[Tuple[type, str], Tuple[str, ...]] = {
    (Chemistry, "evidence"): EVIDENCE_CLASSES,
    (ASO, "sequence_evidence"): EVIDENCE_CLASSES,
    (SpliceGeometry, "evidence"): EVIDENCE_CLASSES,
    (SpliceGeometry, "region"): REGIONS,
    (SpliceGeometry, "region_5p"): END_REGIONS,
    (SpliceGeometry, "region_3p"): END_REGIONS,
    (Placement, "evidence"): EVIDENCE_CLASSES,
    (Placement, "strand"): STRANDS,
    (Placement, "assembly"): (ASSEMBLY,),
    (Placement, "annotation"): (ANNOTATION,),
    (PlacementResult, "role"): ROLES,
    (PlacementResult, "status"): STATUSES,
    (Prediction, "evidence"): EVIDENCE_CLASSES,
    (Prediction, "identifiability"): IDENTIFIABILITY,
    (Intervention, "schema_version"): (SCHEMA_VERSION,),
}
_PRIMITIVES = {
    "int": ("integer", int),
    "str": ("string", str),
    "bool": ("boolean", bool),
    "Optional[float]": (["number", "null"], (int, float, type(None))),
}
_SCHEMA_CLASSES = (Modification, Chemistry, ASO, SpliceGeometry, Placement, PlacementResult, Prediction, Intervention)


def to_dict(obj: Any) -> Any:
    """Plain JSON-able structure; keys in dataclass field order."""
    if is_dataclass(obj) and not isinstance(obj, type):
        return {f.name: to_dict(getattr(obj, f.name)) for f in fields(obj)}
    if isinstance(obj, (tuple, list)):
        return [to_dict(x) for x in obj]
    return obj


def _decode(cls: type, data: Any, path: str) -> Any:
    _req(isinstance(data, dict), f"{path}: expected object")
    names = [f.name for f in fields(cls)]
    missing = [n for n in names if n not in data]
    extra = sorted(k for k in data if k not in names)
    _req(not missing, f"{path}: missing keys {missing}")
    _req(not extra, f"{path}: unknown keys {extra}")
    kwargs: Dict[str, Any] = {}
    for f in fields(cls):
        v, sub = data[f.name], f"{path}.{f.name}"
        nested = _NESTED.get(cls, {}).get(f.name)
        if nested is None:
            _, py = _PRIMITIVES[str(f.type)]
            ok = isinstance(v, py) and not (isinstance(v, bool) and str(f.type) != "bool")
            _req(ok, f"{sub}: expected {f.type}, got {type(v).__name__}")
            kwargs[f.name] = v
        elif nested[0] == "one":
            kwargs[f.name] = _decode(nested[1], v, sub)
        elif nested[0] == "optional":
            kwargs[f.name] = None if v is None else _decode(nested[1], v, sub)
        else:
            _req(isinstance(v, list), f"{sub}: expected array")
            kwargs[f.name] = tuple(_decode(nested[1], x, f"{sub}[{i}]") for i, x in enumerate(v))
    return cls(**kwargs)


def from_dict(data: Any) -> Intervention:
    return _decode(Intervention, data, "$")


def to_json(iv: Intervention) -> str:
    """Canonical JSON: sorted keys, fixed separators, trailing newline."""
    return json.dumps(to_dict(iv), sort_keys=True, indent=2, ensure_ascii=True) + "\n"


def from_json(text: str) -> Intervention:
    return from_dict(json.loads(text))


def json_schema() -> Dict[str, Any]:
    defs: Dict[str, Any] = {}
    for cls in _SCHEMA_CLASSES:
        props: Dict[str, Any] = {}
        for f in fields(cls):
            nested = _NESTED.get(cls, {}).get(f.name)
            if nested is None:
                spec: Dict[str, Any] = {"type": _PRIMITIVES[str(f.type)][0]}
                if (cls, f.name) in _ENUMS:
                    spec["enum"] = list(_ENUMS[(cls, f.name)])
            elif nested[0] == "one":
                spec = {"$ref": f"#/$defs/{nested[1].__name__}"}
            elif nested[0] == "optional":
                spec = {"anyOf": [{"$ref": f"#/$defs/{nested[1].__name__}"}, {"type": "null"}]}
            else:
                spec = {"type": "array", "items": {"$ref": f"#/$defs/{nested[1].__name__}"}}
            props[f.name] = spec
        defs[cls.__name__] = {
            "type": "object",
            "additionalProperties": False,
            "required": [f.name for f in fields(cls)],
            "properties": props,
        }
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": f"conference/model_v2/schema/{SCHEMA_FILENAME}",
        "title": "TransIntervention",
        "description": "ASO intervention on a pre-mRNA locus; 1-based inclusive GRCh38 / GENCODE v50 coordinates.",
        "x-schema-version": SCHEMA_VERSION,
        "x-evidence-classes": list(EVIDENCE_CLASSES),
        "$ref": "#/$defs/Intervention",
        "$defs": defs,
    }


def schema_text() -> str:
    return json.dumps(json_schema(), sort_keys=True, indent=2, ensure_ascii=True) + "\n"


def write_schema(out_dir: Optional[Path] = None) -> Path:
    out = Path(out_dir) if out_dir is not None else SCHEMA_DIR
    out.mkdir(parents=True, exist_ok=True)
    path = out / SCHEMA_FILENAME
    path.write_text(schema_text(), encoding="ascii")
    return path


_PARQUET_COLUMNS = (
    "schema_version", "aso_id", "aso_sequence_5to3", "reference_locus_id", "status",
    "site_id", "chrom", "start", "end", "strand", "payload_json",
)


def to_rows(interventions: Sequence[Intervention]) -> List[Dict[str, Any]]:
    rows = []
    for iv in interventions:
        p = iv.reference.placement
        rows.append({
            "schema_version": iv.schema_version,
            "aso_id": iv.aso.aso_id,
            "aso_sequence_5to3": iv.aso.sequence_5to3,
            "reference_locus_id": iv.reference.locus_id,
            "status": iv.reference.status,
            "site_id": p.site_id if p else None,
            "chrom": p.chrom if p else None,
            "start": p.start if p else None,
            "end": p.end if p else None,
            "strand": p.strand if p else None,
            "payload_json": to_json(iv),
        })
    return rows


def write_parquet(interventions: Sequence[Intervention], path: Any) -> None:
    """One row per intervention; flat reference columns + canonical JSON payload."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    str_t, int_t = pa.string(), pa.int64()
    types = {"start": int_t, "end": int_t}
    schema = pa.schema([(c, types.get(c, str_t)) for c in _PARQUET_COLUMNS])
    rows = to_rows(interventions)
    table = pa.table({c: [r[c] for r in rows] for c in _PARQUET_COLUMNS}, schema=schema)
    pq.write_table(table, str(path), compression="none")


def read_parquet(path: Any) -> List[Intervention]:
    import pyarrow.parquet as pq

    table = pq.read_table(str(path))
    return [from_json(s) for s in table.column("payload_json").to_pylist()]


# --------------------------------------------------------------------------- registry
# ASO sequences as recorded in conference/aso_pipeline/sequences.py. Chemistry is not
# recorded there -> "unknown"/not_available. 3UP8 placement needs human confirmation.
ASOS: Dict[str, ASO] = {
    "F18": ASO("F18", "TCACTTTCATAATGCTGG", note="ISS-N1 core"),
    "F14": ASO("F14", "TTTCATAATGCTGG", note="ISS-N1 3-prime anchored"),
    "F10": ASO("F10", "ATAATGCTGG", note="ISS-N1 3-prime anchored"),
    "3UP8": ASO("3UP8", "GCTGGCAG", requires_human_confirmation=True, note="three-nt-upstream candidate window"),
}

# Anchor window: feng/src/feng/titration/ladder.py (chr5:70076584-70076601, +).
ANCHOR_WINDOW = (70076584, 70076601)

# Exon 7 last base 70076574 / first intron-7 base 70076575: feng/titration/gpu_scores.py.
# Exon 7 first base 70076521 and all SMN1 coordinates are DERIVED here (54-nt exon,
# flanking AG/GT, exact-match offset +875420), not read from the GENCODE GTF.
LOCI: Dict[str, Locus] = {
    "SMN2_exon7": Locus(
        locus_id="SMN2_exon7", gene="SMN2", paralog_group="SMN", role="reference",
        chrom="chr5", strand="+", search_start=70076421, search_end=70076674,
        exon_label="exon7", exon_start=70076521, exon_end=70076574,
        note="exon end from repo; exon start derived (54 nt)",
    ),
    "SMN1_exon7": Locus(
        locus_id="SMN1_exon7", gene="SMN1", paralog_group="SMN", role="alternate",
        chrom="chr5", strand="+", search_start=70951841, search_end=70952094,
        exon_label="exon7", exon_start=70951941, exon_end=70951994,
        note="derived by exact-sequence offset from SMN2_exon7; never joined to it",
    ),
}


# --------------------------------------------------------------------------- CLI
def _parse_window(text: Optional[str]) -> Optional[Tuple[int, int]]:
    if not text:
        return None
    a, _, b = text.partition("-")
    return int(a), int(b)


def _cmd_validate(args: argparse.Namespace) -> int:
    try:
        iv = from_json(Path(args.path).read_text(encoding="utf-8"))
    except (SchemaError, ValueError, OSError) as exc:
        print(f"INVALID: {exc}")
        return 1
    print(f"OK schema_version={iv.schema_version} aso={iv.aso.aso_id} status={iv.reference.status}")
    return 0


def _cmd_place(args: argparse.Namespace) -> int:
    if args.seq:
        aso = ASO(args.aso or "custom", args.seq, requires_human_confirmation=args.require_confirmation)
    elif args.aso in ASOS:
        aso = ASOS[args.aso]
    else:
        print(f"unknown ASO {args.aso!r}; known: {sorted(ASOS)} (or pass --seq)", file=sys.stderr)
        return 2
    for name in [args.locus] + list(args.alt_locus):
        if name not in LOCI:
            print(f"unknown locus {name!r}; known: {sorted(LOCI)}", file=sys.stderr)
            return 2
    try:
        genome = FaiFastaReader(args.fasta)
    except GenomeAccessError as exc:
        print(f"genome unavailable: {exc}", file=sys.stderr)
        return 2
    iv = build_intervention(
        aso, LOCI[args.locus], genome,
        alternate_loci=[LOCI[n] for n in args.alt_locus],
        window=_parse_window(args.window),
    )
    text = to_json(iv)
    if args.out:
        Path(args.out).write_text(text, encoding="ascii")
    else:
        sys.stdout.write(text)
    return {"placed": 0, "needs_human_confirmation": 3}.get(iv.reference.status, 4)


def _cmd_schema(args: argparse.Namespace) -> int:
    if args.stdout:
        sys.stdout.write(schema_text())
    else:
        print(write_schema(Path(args.out_dir) if args.out_dir else None))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m conference.model_v2.trans_intervention")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="validate an intervention JSON file")
    v.add_argument("path")
    v.set_defaults(fn=_cmd_validate)
    p = sub.add_parser("place", help="place an ASO on a locus (exit 0 placed, 3 needs confirmation, 4 no match)")
    p.add_argument("--aso", help=f"registry id {sorted(ASOS)} or a label for --seq")
    p.add_argument("--seq", help="ASO sequence 5'->3' (overrides the registry)")
    p.add_argument("--locus", default="SMN2_exon7")
    p.add_argument("--alt-locus", action="append", default=[])
    p.add_argument("--window", help="START-END (1-based inclusive) inside the locus search window")
    p.add_argument("--fasta", help=f"genome fasta with .fai (default ${FASTA_ENV} or {DEFAULT_FASTA})")
    p.add_argument("--require-confirmation", action="store_true")
    p.add_argument("--out")
    p.set_defaults(fn=_cmd_place)
    s = sub.add_parser("schema", help="emit the JSON schema to schema/")
    s.add_argument("--out-dir")
    s.add_argument("--stdout", action="store_true")
    s.set_defaults(fn=_cmd_schema)
    args = ap.parse_args(argv)
    if args.cmd == "place" and not (args.aso or args.seq):
        ap.error("place needs --aso or --seq")
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
