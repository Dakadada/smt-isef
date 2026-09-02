"""Arm B ENCODE RBP-layer ingest for the clean-room Feng arm.

Claim boundary: this module only *reads and exposes* processed ENCODE data. It
computes no concordance statistic and asserts no interaction; those belong to a
separate, post-freeze step. Every derived summary a caller reports from these
frames must carry an explicit :class:`~feng.honesty.states.ResultState` and go
through the claim linter. The two data roles handled here are:

* ``eclip_peaks`` -- ENCODE eCLIP CLIPper peak calls (BED narrowPeak) that record
  where an RBP crosslinks, with fold enrichment over size-matched input and
  significance. Source publication-data collection ENCSR456FVU.
* ``kd_splicing_rmats`` -- batch-corrected rMATS differential-splicing tables from
  shRNA/CRISPR RBP knockdown RNA-seq, exposing the knockdown delta-PSI direction
  per splicing event. Source publication-data collection ENCSR870OLK.

Provenance state for both sources is ``paper_supported``: they are the published
ENCODE processed outputs, not a clean-room re-derivation. No file here is labelled
``clean_reanalysis_reproduced``.

Assembly note: every ENCODE file staged for this layer is on ``hg19``. The rest of
the arm works in the GENCODE v50 (GRCh38) frame, so coordinate-level use downstream
requires a lift to GRCh38. Liftover is deliberately deferred and is NOT applied by
this module; :func:`load_manifest` surfaces the ``needs_liftover`` flag so callers
cannot silently mix frames.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
import io
from pathlib import Path
import re
import tarfile

import pandas as pd

from ..honesty.states import ResultState

SOURCE_ASSEMBLY = "hg19"
TARGET_ASSEMBLY = "GRCh38"  # GENCODE v50 frame used elsewhere in the arm
PROVENANCE_STATE = ResultState.PAPER_SUPPORTED

# ENCODE eCLIP "peaks" are distributed as BED narrowPeak. For eCLIP the numeric
# columns carry CLIPper/enrichment semantics rather than generic peak scores.
ECLIP_PEAK_COLUMNS = (
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signal_value",  # log2 fold enrichment over size-matched input
    "neg_log10_p",   # -log10 p-value of enrichment
    "neg_log10_q",   # -log10 q-value (often -1 when not computed)
    "peak_offset",   # point-source offset from start (-1 if absent)
)

# rMATS alternative-splicing event types, one MATS table each inside the tar.
RMATS_EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")

# The ENCODE RBP-map tars use the older rMATS "JunctionCountOnly" naming; newer
# rMATS-turbo writes "JC". Accept both so the loader spans both conventions.
_RMATS_TABLE_SUFFIXES = (".MATS.JunctionCountOnly.txt", ".MATS.JC.txt")

# Columns every rMATS MATS.JC table exposes that this layer treats as canonical.
_RMATS_REQUIRED_COLUMNS = (
    "GeneID",
    "geneSymbol",
    "IncLevelDifference",  # knockdown - control delta-PSI; its SIGN is the direction
    "PValue",
    "FDR",
)

MANIFEST_COLUMNS = (
    "collection",
    "accession",
    "role",
    "target_rbp",
    "cell_line",
    "file_format",
    "output_type",
    "assembly",
    "needs_liftover",
    "md5sum",
    "sha256",
    "file_size",
    "experiment",
    "url",
    "curated_subset",
)

_VALID_ROLES = frozenset({"eclip_peaks", "kd_splicing_rmats"})


class AssemblyDataFrame(pd.DataFrame):
    """Coordinate frame whose merge operations enforce assembly provenance."""

    @property
    def _constructor(self):  # type: ignore[override]
        return AssemblyDataFrame

    def merge(self, right, *args, **kwargs):  # type: ignore[override]
        assembly = assert_same_assembly(self, right)
        result = super().merge(right, *args, **kwargs)
        result.attrs["assembly"] = assembly
        return result

    def join(self, other, *args, **kwargs):  # type: ignore[override]
        others = list(other) if isinstance(other, list) else [other]
        assembly = assert_same_assembly(self, *others)
        result = super().join(other, *args, **kwargs)
        result.attrs["assembly"] = assembly
        return result


@dataclass(frozen=True, slots=True)
class EncodeFile:
    """One processed ENCODE file as recorded in the RBP-layer manifest."""

    collection: str
    accession: str
    role: str
    target_rbp: str
    cell_line: str
    file_format: str
    output_type: str
    assembly: str
    needs_liftover: bool
    md5sum: str
    sha256: str | None
    file_size: int | None
    experiment: str
    url: str
    curated_subset: bool

    def __post_init__(self) -> None:
        if self.role not in _VALID_ROLES:
            raise ValueError(f"unknown RBP-layer role: {self.role!r}")
        if not self.accession.strip():
            raise ValueError("manifest row missing ENCODE accession")
        if self.assembly != SOURCE_ASSEMBLY:
            raise ValueError(
                f"{self.accession}: expected source assembly {SOURCE_ASSEMBLY!r}, "
                f"got {self.assembly!r}"
            )
        if not self.needs_liftover:
            raise ValueError(
                f"{self.accession}: hg19 source must be marked needs_liftover"
            )
        if not re.fullmatch(r"[0-9a-f]{32}", self.md5sum):
            raise ValueError(f"{self.accession}: invalid md5sum")
        if self.sha256 is not None and not re.fullmatch(r"[0-9a-f]{64}", self.sha256):
            raise ValueError(f"{self.accession}: invalid sha256")
        if self.curated_subset and self.sha256 is None:
            raise ValueError(f"{self.accession}: curated file is missing sha256")


def _as_bool(value: str) -> bool:
    token = str(value).strip().casefold()
    if token in {"yes", "true", "1"}:
        return True
    if token in {"no", "false", "0", ""}:
        return False
    raise ValueError(f"expected a yes/no flag, got {value!r}")


def load_manifest(path: str | Path) -> tuple[EncodeFile, ...]:
    """Read the RBP-layer manifest TSV into typed, validated rows."""
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = [c for c in MANIFEST_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"manifest is missing columns: {missing}")
    rows: list[EncodeFile] = []
    for record in frame.to_dict("records"):
        size = record["file_size"].strip()
        rows.append(
            EncodeFile(
                collection=record["collection"],
                accession=record["accession"],
                role=record["role"],
                target_rbp=record["target_rbp"],
                cell_line=record["cell_line"],
                file_format=record["file_format"],
                output_type=record["output_type"],
                assembly=record["assembly"],
                needs_liftover=_as_bool(record["needs_liftover"]),
                md5sum=record["md5sum"],
                sha256=record["sha256"].strip() or None,
                file_size=int(size) if size else None,
                experiment=record["experiment"],
                url=record["url"],
                curated_subset=_as_bool(record["curated_subset"]),
            )
        )
    accessions = [row.accession for row in rows]
    if len(accessions) != len(set(accessions)):
        raise ValueError("manifest accessions must be unique")
    return tuple(rows)


def _attach_source_provenance(
    frame: pd.DataFrame, source: EncodeFile, *, expected_role: str
) -> pd.DataFrame:
    """Attach the hard assembly boundary to every row returned by a reader."""
    if source.role != expected_role:
        raise ValueError(
            f"{source.accession}: expected role {expected_role!r}, got {source.role!r}"
        )
    # EncodeFile validates hg19 + needs_liftover at construction time. Repeat
    # the assertions here so a future relaxed manifest model cannot weaken the
    # coordinate-reader boundary.
    if source.assembly != SOURCE_ASSEMBLY or not source.needs_liftover:
        raise ValueError(
            f"{source.accession}: refusing coordinate data without enforced "
            f"{SOURCE_ASSEMBLY}->{TARGET_ASSEMBLY} liftover provenance"
        )
    frame = AssemblyDataFrame(frame.copy())
    frame["source_accession"] = source.accession
    frame["source_assembly"] = source.assembly
    frame["assembly"] = source.assembly
    frame["needs_liftover"] = source.needs_liftover
    frame["source_sha256"] = source.sha256
    frame.attrs["assembly"] = source.assembly
    frame.attrs["needs_liftover"] = source.needs_liftover
    return frame


def assert_same_assembly(*frames: pd.DataFrame) -> str:
    """Hard-fail a coordinate join unless every input declares one assembly."""
    assemblies: set[object] = set()
    for frame in frames:
        declared = frame.attrs.get("assembly")
        if "assembly" in frame.columns:
            column_values = set(frame["assembly"].dropna().unique())
            if len(column_values) != 1:
                raise ValueError(
                    f"coordinate frame has mixed assembly values: {column_values}"
                )
            column_assembly = next(iter(column_values), None)
            if declared is not None and column_assembly != declared:
                raise ValueError(
                    "coordinate frame column/attribute assembly provenance disagrees"
                )
            declared = column_assembly
        assemblies.add(declared)
    if None in assemblies:
        raise ValueError("coordinate frame lacks assembly provenance")
    if len(assemblies) != 1:
        raise ValueError(f"cannot join coordinate frames across assemblies: {assemblies}")
    return str(next(iter(assemblies)))


def coordinate_frame(frame: pd.DataFrame, *, assembly: str) -> AssemblyDataFrame:
    """Declare a coordinate frame so every subsequent ``merge``/``join`` is fenced."""
    if not assembly.strip():
        raise ValueError("coordinate assembly must be non-empty")
    result = AssemblyDataFrame(frame.copy())
    result["assembly"] = assembly
    result.attrs["assembly"] = assembly
    return result


def _verify_source_sha256(path: Path, source: EncodeFile) -> str:
    """Hash bytes on disk and reject absent or mismatched manifest provenance."""
    if source.sha256 is None:
        raise ValueError(
            f"{source.accession}: cannot read an unverified file without sha256"
        )
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
    except OSError as exc:
        raise ValueError(f"{source.accession}: cannot read source file {path}") from exc
    actual = digest.hexdigest()
    if actual != source.sha256:
        raise ValueError(
            f"{source.accession}: SHA-256 mismatch for {path}: "
            f"expected {source.sha256}, got {actual}"
        )
    return actual


def read_eclip_peaks(path: str | Path, *, source: EncodeFile) -> pd.DataFrame:
    """Load one eCLIP narrowPeak BED (optionally gzipped) into a tidy frame.

    Coordinates are raw ``hg19`` half-open BED intervals. Source assembly,
    accession, true source assembly, hash, and required-liftover provenance are
    attached to every row; no lift to GRCh38 is applied here.
    """
    path = Path(path)
    _verify_source_sha256(path, source)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        frame = pd.read_csv(
            handle,
            sep="\t",
            header=None,
            comment="#",
            names=list(ECLIP_PEAK_COLUMNS),
            usecols=range(len(ECLIP_PEAK_COLUMNS)),
        )
    for col in ("start", "end", "peak_offset"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce").astype("Int64")
    for col in ("signal_value", "neg_log10_p", "neg_log10_q"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return _attach_source_provenance(frame, source, expected_role="eclip_peaks")


def _read_mats_member(handle: io.BufferedReader) -> pd.DataFrame:
    frame = pd.read_csv(handle, sep="\t")
    missing = [c for c in _RMATS_REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"rMATS table missing required columns: {missing}")
    frame["IncLevelDifference"] = pd.to_numeric(
        frame["IncLevelDifference"], errors="coerce"
    )
    frame["FDR"] = pd.to_numeric(frame["FDR"], errors="coerce")
    return frame


def read_kd_splicing_rmats(
    path: str | Path, *, source: EncodeFile
) -> dict[str, pd.DataFrame]:
    """Extract per-event-type rMATS MATS.JC tables from a knockdown tar(.gz).

    Returns a mapping from event type in :data:`RMATS_EVENT_TYPES` to its tidy
    frame. ``IncLevelDifference`` is the knockdown-minus-control delta-PSI; only
    its sign (the direction) is used downstream. This function parses and exposes;
    it computes no concordance and makes no interaction claim.
    """
    tables: dict[str, pd.DataFrame] = {}
    path = Path(path)
    _verify_source_sha256(path, source)
    with tarfile.open(path, "r:*") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            base = Path(member.name).name
            suffix = next((s for s in _RMATS_TABLE_SUFFIXES if base.endswith(s)), None)
            if suffix is None:
                continue
            event = base[: -len(suffix)]
            if event not in RMATS_EVENT_TYPES:
                continue
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            tables[event] = _attach_source_provenance(
                _read_mats_member(extracted),
                source,
                expected_role="kd_splicing_rmats",
            )
    if not tables:
        raise ValueError(f"no rMATS MATS event tables found in {path}")
    return tables
