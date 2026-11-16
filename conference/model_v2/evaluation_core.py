"""F3 - evaluation core: label states, leakage-free splits, guarded metrics, clustered uncertainty.

Statistics / software-engineering module only. It fits nothing and reads no outcome
data. Five pieces, each refusing loudly instead of degrading silently:

1. ``validate_table``  - exhaustive label states; a ``verified_negative`` needs explicit
   ``negative_evidence`` (assayed_and_powered is True + finite nonnegative detection limit) and can
   never come from a responder-only table (absence is not a negative).
2. ``make_splits`` / ``make_nested_splits`` - purged group splits. Five leakage axes:
   near-identical ASO sequence, chemistry family, study, window, independent unit.
   The test set is one unit of the mode's primary axis; every other row sharing ANY
   active axis group with a test row is purged (in neither train nor test).
3. ``compute_metric`` - AUROC/AUPRC/FDR/accuracy/precision/recall/calibration only with
   enough verified negatives and positives (counted as independent units), otherwise
   ``MetricBlocked``. ``pu_report`` is the positive-unlabeled alternative and never
   emits a supervised metric. Unlabeled / unassayed / censored are never negatives.
4. ``cluster_bootstrap_ci`` - resamples effective independent clusters
   (``independent_unit_id`` merged with sequence / window overlap); refuses below
   ``min_clusters``.
5. ``BASELINES`` (M0..M6) and ``run_ablation`` (leave-feature-family-out interface).

Identifiers (``assay_id``, ``study_id``, ``independent_unit_id``, ``observation_id``)
are opaque, CASE-SENSITIVE strings of 1..128 ASCII characters from
``[A-Za-z0-9._:+-]``. Whitespace, control/format characters and all non-ASCII
characters raise ``SchemaError``; IDs are never stripped, casefolded or Unicode
normalized. Public table consumers (including ``read_table``) validate this contract.
Sequence, chromosome and chemistry aliases remain canonicalized for assay identity.
Callers must assign distinct assay IDs to genuine replicate assays within a study.

``check_source_not_forbidden`` refuses the July exon-7 outcome workbook, any ``psi25``
directory and anything under ``feng/terminal`` (SCIENTIFIC_LOCK).
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import random
import re
import unicodedata
from dataclasses import asdict, dataclass
from enum import Enum
from types import MappingProxyType
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import unquote

from conference.model_v2.trans_intervention import EVIDENCE_CLASSES

__all__ = [
    "LabelState", "ObservationTable", "validate_table", "read_table",
    "GroupingConfig", "Fold", "SplitPlan", "NestedSplitPlan", "axis_groups",
    "make_splits", "make_nested_splits", "verify_no_leakage",
    "MetricPolicy", "MetricResult", "compute_metric", "precision_recall_curve",
    "PUQuantity", "PUReport", "pu_report",
    "CIResult", "effective_clusters", "cluster_bootstrap_ci", "bootstrap_metric_ci",
    "FeatureFamily", "BaselineSpec", "BASELINES", "baseline_feature_columns",
    "permuted_reference_scores", "AblationCell", "AblationResult", "run_ablation",
    "check_source_not_forbidden",
]


# ------------------------------------------------------------------ errors
class SchemaError(ValueError):
    """A row or argument violates the observation-table contract."""


class NegativeEvidenceError(SchemaError):
    """A verified_negative lacks explicit, powered negative evidence."""


class LabelStateError(SchemaError):
    """A non-negative label state was used (or encoded) as a negative."""


class ForbiddenSourceError(RuntimeError):
    """An input path matches a forbidden outcome source."""


class SplitRefused(RuntimeError):
    """A split was refused (independent N too small, or nothing left to train on)."""


class LeakageError(RuntimeError):
    """A group straddles train and test. Indicates a bug; never expected."""


class MetricBlocked(RuntimeError):
    """A metric was refused; ``reason`` says why."""

    def __init__(self, metric: str, reason: str):
        super().__init__(f"{metric}: BLOCKED - {reason}")
        self.metric = metric
        self.reason = reason


class UncertaintyRefused(RuntimeError):
    """A confidence interval was refused (too few independent clusters)."""


class BaselineUndefined(ValueError):
    """The baseline has no agreed feature set yet."""


# ------------------------------------------------------------------ forbidden sources
_FORBIDDEN_PATTERNS: Tuple[Tuple[str, "re.Pattern[str]"], ...] = (
    ("july_exon7_outcome_workbook", re.compile(
        r"smn2[\W_]*exon[\W_]*7[\W_]*7[\W_]*27[\W_]*28[\W_]*complete[\W_]*data", re.I)),
    ("psi25", re.compile(r"(^|/)psi25(?![A-Za-z0-9])", re.I)),
    ("feng_terminal", re.compile(r"(^|/)feng/+terminal(?![A-Za-z0-9])", re.I)),
)


def _normalized_source_forms(raw: str) -> set[str]:
    forms = {raw}
    current = raw
    for _ in range(16):
        normalized = unicodedata.normalize("NFKC", unquote(current)).replace("\\", "/")
        forms.add(normalized)
        if normalized == current:
            break
        current = normalized
    else:
        raise ForbiddenSourceError("source normalization did not converge within 16 steps")
    for item in tuple(forms):
        expanded = os.path.expanduser(item.replace("\\", "/"))
        forms.update((os.path.normpath(expanded), os.path.abspath(expanded)))
        try:
            resolved = os.path.realpath(expanded)
            # Symlink targets may themselves use encoded or Unicode names.
            for _ in range(16):
                normalized = unicodedata.normalize("NFKC", unquote(resolved)).replace("\\", "/")
                forms.update((normalized, os.path.normpath(normalized)))
                if normalized == resolved:
                    break
                resolved = normalized
            else:
                raise ForbiddenSourceError("resolved source normalization did not converge")
        except (OSError, ValueError):
            pass
    return forms


def check_source_not_forbidden(path: Any) -> str:
    """Return ``str(path)`` or raise ``ForbiddenSourceError``.

    Checks the literal string, its URL-decoded form, the normalised absolute path
    and the symlink-resolved path. Only path names are inspected; nothing is opened.
    """
    raw = os.fspath(path) if isinstance(path, os.PathLike) else str(path)
    for form in _normalized_source_forms(raw):
        text = form.replace("\\", "/")
        for label, pattern in _FORBIDDEN_PATTERNS:
            if pattern.search(text):
                raise ForbiddenSourceError(
                    f"refused forbidden outcome source ({label}): {raw!r}")
    return raw


def _strings(value: Any) -> List[str]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, os.PathLike):
        return [os.fspath(value)]
    if isinstance(value, Mapping):
        return [s for v in value.values() for s in _strings(v)]
    if isinstance(value, (list, tuple, set, frozenset)):
        return [s for v in value for s in _strings(v)]
    return []


# ------------------------------------------------------------------ label states / table
class LabelState(str, Enum):
    MEASURED_POSITIVE = "measured_positive"
    VERIFIED_NEGATIVE = "verified_negative"
    UNLABELED = "unlabeled"
    UNASSAYED = "unassayed"
    CENSORED = "censored"

    @classmethod
    def parse(cls, value: Any) -> "LabelState":
        if isinstance(value, cls):
            return value
        try:
            return cls(value)
        except ValueError:
            raise SchemaError(
                f"label_state {value!r} not in {[s.value for s in cls]}") from None


LABEL_STATES = tuple(s.value for s in LabelState)
POSITIVE = LabelState.MEASURED_POSITIVE.value
NEGATIVE = LabelState.VERIFIED_NEGATIVE.value
IDENTIFIER_FIELDS = ("assay_id", "study_id", "independent_unit_id", "observation_id")
_IDENTIFIER_RE = re.compile(r"[A-Za-z0-9._:+-]{1,128}")
REQUIRED_FIELDS = ("observation_id", "label_state", "assay_id", "study_id",
                   "independent_unit_id", "provenance")
RESPONDER_ONLY_KINDS = frozenset({"responder_only", "responders_only", "hits_only",
                                  "positives_only"})
ABSENCE_BASES = frozenset({"absence", "absent", "absent_from_table", "not_listed",
                           "not_reported", "inferred", "assumed", "absence_from_responder_list"})
NEGATIVE_BASES = frozenset({"dose_response", "powered_assay", "direct_measurement",
                            "assayed_nonresponder", "measured_negative"})
SOURCE_TABLE_KINDS = RESPONDER_ONLY_KINDS | {"full_screen", "all_assayed", "complete",
                                           "dose_response", "assay_results"}
RESPONDER_ONLY_SOURCES = frozenset({"holgersen_s4", "holgersen_table_s4",
                                   "holgersen_supplementary_table_s4"})
BINARY_LABEL_FIELDS = ("binary_label", "label", "y")
# Closed identity contract: assay_id within a study must distinguish experimental
# conditions and genuine replicate assays. Row/unit names, labels, evidence and
# arbitrary annotations cannot turn one assay into multiple observations.
ASSAY_IDENTITY_FIELDS = ("study_id", "assay_id", "aso_sequence", "chemistry_family",
                         "chrom", "start", "end", "strand")


def _category(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", unicodedata.normalize("NFKC", value).casefold()).strip("_")


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({k: _freeze(v) for k, v in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(v) for v in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(v) for v in value)
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise SchemaError(f"unsupported mutable metadata value: {type(value).__name__}")


@dataclass(frozen=True)
class ObservationTable:
    """Deeply immutable validated rows; public consumers revalidate even this wrapper.

    ``duplicate_assay_records`` lists rows sharing the canonical assay-defining
    fields, regardless of metadata. They count as one independent unit for grouping;
    statistics refuse duplicate contributors rather than giving copies extra weight.
    """
    rows: Tuple[Mapping[str, Any], ...]
    duplicate_assay_records: Tuple[Tuple[str, ...], ...] = ()

    def __len__(self) -> int:
        return len(self.rows)

    @property
    def ids(self) -> Tuple[str, ...]:
        return tuple(r["observation_id"] for r in self.rows)

    def subset(self, indices: Sequence[int]) -> "ObservationTable":
        return ObservationTable(tuple(self.rows[i] for i in indices))


def _nonempty(row: Mapping[str, Any], name: str, where: str) -> str:
    value = row.get(name)
    if not isinstance(value, str) or not value.strip():
        raise SchemaError(f"{where}: field {name!r} must be a non-empty string")
    return value


def _identifier(value: Any, name: str, where: str) -> str:
    if not isinstance(value, str) or _IDENTIFIER_RE.fullmatch(value) is None:
        raise SchemaError(
            f"{where}: field {name!r} must be a case-sensitive ASCII identifier "
            "of 1..128 characters from [A-Za-z0-9._:+-]")
    return value


def _canonicalize_assay_identity(row: Dict[str, Any], where: str) -> None:
    """Canonicalize supplied identity fields once, before hashing or consumption.

    Optional fields may be absent, but supplied values must be valid. Coordinates
    are already canonical integers; alternate scalar types fail closed.
    """
    if "aso_sequence" in row:
        row["aso_sequence"] = normalize_sequence(row["aso_sequence"])
    if "chemistry_family" in row:
        row["chemistry_family"] = normalize_chemistry(_nonempty(row, "chemistry_family", where))
    if "chrom" in row:
        row["chrom"] = _canonical_chromosome(_nonempty(row, "chrom", f"{where}: window axis"))
    for name in ("start", "end"):
        if name in row:
            _positive_int(row[name], f"{where}: window axis {name}")
    if "start" in row and "end" in row and row["end"] < row["start"]:
        raise SchemaError(f"{where}: window axis needs start<=end")
    if "strand" in row:
        strand = unicodedata.normalize("NFKC", _nonempty(row, "strand", f"{where}: window axis")).strip()
        if strand not in ("+", "-"):
            raise SchemaError(f"{where}: window axis strand must be +/-")
        row["strand"] = strand


def validate_table(rows: Any, *, responder_only_studies: Sequence[str] = ()) -> ObservationTable:
    """Validate observation rows (mappings) and return an ``ObservationTable``.

    ``responder_only_studies``: study ids whose tables list responders only (e.g. a
    hits-only supplement). No ``verified_negative`` may originate from them, and none
    from a row declaring ``source_table_kind`` in ``RESPONDER_ONLY_KINDS``.
    """
    if isinstance(rows, ObservationTable):
        rows = rows.rows
    responder_only = frozenset(
        _identifier(s, "study_id", "responder_only_studies") for s in responder_only_studies)
    out: List[Dict[str, Any]] = []
    seen = set()
    for index, row in enumerate(rows):
        where = f"row {index}"
        if not isinstance(row, Mapping):
            raise SchemaError(f"{where}: not a mapping")
        r = dict(row)
        for name in IDENTIFIER_FIELDS:
            _identifier(r.get(name), name, where)
        oid = r["observation_id"]
        where = f"row {index} ({oid})"
        if oid in seen:
            raise SchemaError(f"{where}: duplicate observation_id")
        seen.add(oid)
        if "label_state" not in r:
            raise SchemaError(f"{where}: missing label_state")
        state = LabelState.parse(r["label_state"])
        r["label_state"] = state.value
        _canonicalize_assay_identity(r, where)
        prov = r.get("provenance")
        if isinstance(prov, Mapping):
            if not isinstance(prov.get("source"), str) or not prov["source"].strip():
                raise SchemaError(f"{where}: provenance mapping needs a non-empty 'source'")
        elif not isinstance(prov, str) or not prov.strip():
            raise SchemaError(f"{where}: provenance must be a non-empty string or mapping")
        for text in _strings(prov) + _strings(r.get("source_path")):
            check_source_not_forbidden(text)
        if r.get("evidence_class") is not None:
            r["evidence_class"] = _category(str(r["evidence_class"]))
            if r["evidence_class"] not in EVIDENCE_CLASSES:
                raise SchemaError(f"{where}: bad evidence_class {r['evidence_class']!r}")
        kind = _category(str(r.get("source_table_kind", "")))
        if kind and kind not in SOURCE_TABLE_KINDS:
            raise NegativeEvidenceError(f"{where}: unknown source_table_kind {kind!r}")
        if "source_table_kind" in r:
            r["source_table_kind"] = kind

        if state is LabelState.VERIFIED_NEGATIVE:
            sources = [r["study_id"]] + _strings(prov) + _strings(r.get("source_path"))
            source_tokens = [_category(form) for s in sources for form in _normalized_source_forms(s)]
            known_responder = any(
                f"_{known}_" in f"_{token}_" for known in RESPONDER_ONLY_SOURCES for token in source_tokens)
            if kind in RESPONDER_ONLY_KINDS or r["study_id"] in responder_only or known_responder:
                raise NegativeEvidenceError(
                    f"{where}: verified_negative from a responder-only table; absence "
                    "from a responder list is not a measured negative")
            ne = r.get("negative_evidence")
            if not isinstance(ne, Mapping):
                raise NegativeEvidenceError(
                    f"{where}: verified_negative requires an explicit negative_evidence mapping")
            if ne.get("assayed_and_powered") is not True:
                raise NegativeEvidenceError(
                    f"{where}: negative_evidence.assayed_and_powered must be exactly True")
            limit = ne.get("detection_limit")
            if isinstance(limit, bool) or not isinstance(limit, (int, float)) \
                    or not math.isfinite(limit) or limit < 0:
                raise NegativeEvidenceError(
                    f"{where}: negative_evidence.detection_limit must be a finite nonnegative "
                    "number; signed/log-scale conventions are not supported")
            basis = _category(str(ne.get("basis", "")))
            if basis not in NEGATIVE_BASES:
                raise NegativeEvidenceError(
                    f"{where}: negative_evidence.basis={ne.get('basis')!r} is absence-inferred "
                    "or not a recognized measured-evidence category")
            r["negative_evidence"] = {**ne, "basis": basis}
        for name in BINARY_LABEL_FIELDS:
            value = r.get(name)
            if value is None:
                continue
            if value not in (0, 1):
                raise SchemaError(f"{where}: {name}={value!r} is not binary")
            needed = LabelState.MEASURED_POSITIVE if value == 1 else LabelState.VERIFIED_NEGATIVE
            if state is not needed:
                raise LabelStateError(
                    f"{where}: {name}={int(value)} encodes {needed.value} but label_state is "
                    f"{state.value}; {state.value} is never a "
                    f"{'negative' if value == 0 else 'positive'}")
        out.append(_freeze(r))
    records: Dict[str, List[str]] = {}
    for r in out:
        identity = {k: r.get(k) for k in ASSAY_IDENTITY_FIELDS}
        records.setdefault(_fingerprint(identity), []).append(r["observation_id"])
    duplicates = tuple(tuple(ids) for ids in records.values() if len(ids) > 1)
    return ObservationTable(tuple(out), duplicates)


def _as_table(rows: Any) -> ObservationTable:
    return validate_table(rows)


def read_table(path: Any, **kwargs: Any) -> ObservationTable:
    """Read a ``.jsonl`` / ``.json`` observation table, after the forbidden-source guard."""
    name = check_source_not_forbidden(path)
    with open(os.path.expanduser(name), "r", encoding="utf-8") as handle:
        if name.endswith(".jsonl"):
            rows = [json.loads(line) for line in handle if line.strip()]
        elif name.endswith(".json"):
            rows = json.load(handle)
        else:
            raise SchemaError(f"unsupported table format: {name!r} (use .jsonl or .json)")
    return validate_table(rows, **kwargs)


# ------------------------------------------------------------------ grouping
AXES = ("sequence", "chemistry", "study", "window", "independent_unit")
MODES = {"leave_one_aso_out": "sequence", "leave_study_out": "study",
         "region_blocked": "window"}
_SEQ_RE = re.compile(r"^[ACGTN]+$")

# Closed vocabulary. PS alone is a backbone, not an unambiguous chemistry family.
CHEMISTRY_ALIASES = {
    "moe": "moe", "2moe": "moe", "2omoe": "moe", "2omethoxyethyl": "moe",
    "methoxyethyl": "moe",
    "ome": "ome", "2ome": "ome", "2omethyl": "ome",
    "lna": "lna", "lockednucleicacid": "lna",
    "pmo": "pmo", "morpholino": "pmo", "phosphorodiamidatemorpholino": "pmo",
    "pna": "pna", "peptidenucleicacid": "pna", "dna": "dna", "rna": "rna",
    "cet": "cet", "constrainedethyl": "cet",
}


def normalize_chemistry(value: str) -> str:
    text = unicodedata.normalize("NFKC", value).casefold().strip()
    # Remove only known orthographic separators; do not turn MOE/LNA or MOE? into aliases.
    key = re.sub(r"[\s\-'′’]", "", text)
    if key not in CHEMISTRY_ALIASES:
        raise SchemaError(f"unknown or ambiguous chemistry_family {value!r}")
    return CHEMISTRY_ALIASES[key]


def _positive_int(value: Any, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise SchemaError(f"{name} must be a finite positive integer")


@dataclass(frozen=True)
class GroupingConfig:
    """Leakage axes and their thresholds.

    ``relax`` removes axes from the keep-together set; it is recorded on every plan
    as ``relaxed_axes`` (a declared leakage risk), never applied silently.
    ``window_padding`` pads each window on both sides before clustering (1-based
    closed coordinates, as in J6 ``start1``/``end1``).
    """
    sequence_metric: str = "edit"            # "edit" (Levenshtein) or "hamming"
    sequence_threshold: int = 2
    containment_is_identical: bool = True    # truncations (e.g. 14-mer inside 18-mer)
    min_containment_length: int = 8
    window_padding: int = 0
    strand_aware: bool = True
    relax: Tuple[str, ...] = ()
    max_unique_sequences: int = 4000         # O(n^2) pairwise guard (login-node safe)

    def __post_init__(self) -> None:
        if self.sequence_metric not in ("edit", "hamming"):
            raise SchemaError(f"bad sequence_metric {self.sequence_metric!r}")
        for name in ("sequence_threshold", "window_padding"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise SchemaError(f"{name} must be a finite nonnegative integer")
        for name in ("min_containment_length", "max_unique_sequences"):
            _positive_int(getattr(self, name), name)
        for name in ("containment_is_identical", "strand_aware"):
            if not isinstance(getattr(self, name), bool):
                raise SchemaError(f"{name} must be a boolean")
        if not isinstance(self.relax, (tuple, list)) or any(
                not isinstance(axis, str) for axis in self.relax):
            raise SchemaError("relax must be a sequence of axis names")
        if len(set(self.relax)) != len(self.relax):
            raise SchemaError("relax must not repeat axes")
        bad = set(self.relax) - set(AXES)
        if bad:
            raise SchemaError(f"unknown axes in relax: {sorted(bad)}")
        if "independent_unit" in self.relax:
            raise SchemaError("independent_unit cannot be relaxed")
        object.__setattr__(self, "relax", tuple(self.relax))

    @property
    def axes(self) -> Tuple[str, ...]:
        return tuple(a for a in AXES if a not in self.relax)


class _UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def normalize_sequence(seq: Any) -> str:
    if not isinstance(seq, str):
        raise SchemaError(f"aso_sequence must be a string, got {seq!r}")
    out = "".join(seq.split()).upper().replace("U", "T")
    if not _SEQ_RE.match(out):
        raise SchemaError(f"aso_sequence {seq!r} is not ACGTUN")
    return out


def _edit_within(a: str, b: str, limit: int) -> bool:
    if abs(len(a) - len(b)) > limit:
        return False
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        if min(cur) > limit:
            return False
        prev = cur
    return prev[-1] <= limit


def _near_identical(a: str, b: str, cfg: GroupingConfig) -> bool:
    if a == b:
        return True
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    if cfg.containment_is_identical and len(short) >= cfg.min_containment_length \
            and short in long_:
        return True
    if cfg.sequence_metric == "hamming":
        return len(a) == len(b) and sum(x != y for x, y in zip(a, b)) <= cfg.sequence_threshold
    return _edit_within(a, b, cfg.sequence_threshold)


def _sequence_groups(table: ObservationTable, cfg: GroupingConfig) -> List[str]:
    seqs = [normalize_sequence(r.get("aso_sequence")) for r in table.rows]
    uniq = sorted(set(seqs))
    if len(uniq) > cfg.max_unique_sequences:
        raise SplitRefused(
            f"{len(uniq)} unique sequences exceeds max_unique_sequences="
            f"{cfg.max_unique_sequences}; pairwise clustering refused on this node")
    uf = _UnionFind(len(uniq))
    for i in range(len(uniq)):
        for j in range(i + 1, len(uniq)):
            if _near_identical(uniq[i], uniq[j], cfg):
                uf.union(i, j)
    label = {s: f"seq:{uniq[uf.find(i)]}" for i, s in enumerate(uniq)}  # root = smallest
    return [label[s] for s in seqs]


def _canonical_chromosome(chrom: str) -> str:
    """GRCh38 primary chromosomes only; unknown contig aliases fail closed."""
    name = chrom.strip().upper()
    if name.startswith("CHR"):
        name = name[3:]
    if name in ("M", "MT"):
        return "chrM"
    if name in {str(i) for i in range(1, 23)} | {"X", "Y"}:
        return "chr" + name
    raise SchemaError(f"unrecognized GRCh38 primary chromosome {chrom!r}")


def _window_groups(table: ObservationTable, cfg: GroupingConfig) -> List[str]:
    keyed: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = {}
    for i, r in enumerate(table.rows):
        chrom, start, end, strand = r.get("chrom"), r.get("start"), r.get("end"), r.get("strand")
        ok = (isinstance(chrom, str) and chrom and strand in ("+", "-")
              and all(isinstance(v, int) and not isinstance(v, bool) for v in (start, end))
              and 1 <= start <= end)
        if not ok:
            raise SchemaError(
                f"row {i} ({r['observation_id']}): window axis needs chrom, 1-based closed "
                "int start<=end, strand in +/-")
        chrom = _canonical_chromosome(chrom)
        keyed.setdefault((chrom, strand if cfg.strand_aware else "*"), []).append((start, end, i))
    out = [""] * len(table)
    for (chrom, strand), items in sorted(keyed.items()):
        items.sort()
        cluster: List[int] = []
        lo = hi = 0
        for start, end, i in items + [(None, None, -1)]:  # sentinel flushes the last cluster
            if cluster and (start is None or start > hi + 2 * cfg.window_padding):
                for member in cluster:
                    out[member] = f"win:{chrom}:{strand}:{lo}-{hi}"
                cluster = []
            if start is None:
                break
            if not cluster:
                lo, hi = start, end
            cluster.append(i)
            hi = max(hi, end)
    return out


def axis_groups(rows: Any, config: GroupingConfig = GroupingConfig(),
                axes: Optional[Sequence[str]] = None) -> Dict[str, Tuple[str, ...]]:
    """Per-row group id on each requested leakage axis."""
    config.__post_init__()
    table = _as_table(rows)
    out: Dict[str, Tuple[str, ...]] = {}
    for axis in (config.axes if axes is None else axes):
        if axis == "sequence":
            out[axis] = tuple(_sequence_groups(table, config))
        elif axis == "chemistry":
            out[axis] = tuple(
                "chem:" + normalize_chemistry(_nonempty(r, "chemistry_family", f"row {i}"))
                for i, r in enumerate(table.rows))
        elif axis == "study":
            out[axis] = tuple(f"study:{r['study_id']}" for r in table.rows)
        elif axis == "window":
            out[axis] = tuple(_window_groups(table, config))
        elif axis == "independent_unit":
            out[axis] = tuple(f"unit:{r['independent_unit_id']}" for r in table.rows)
        else:
            raise SchemaError(f"unknown axis {axis!r}")
    return out


def _components(table: ObservationTable, group_lists: Sequence[Sequence[str]]) -> List[str]:
    """Connected components of rows linked through any of the given group ids."""
    uf = _UnionFind(len(table))
    for groups in group_lists:
        first: Dict[str, int] = {}
        for i, g in enumerate(groups):
            uf.union(first.setdefault(g, i), i)
    ids = table.ids
    smallest: Dict[int, str] = {}
    for i in range(len(table)):
        root = uf.find(i)
        smallest[root] = min(smallest.get(root, ids[i]), ids[i])
    return [f"component:{smallest[uf.find(i)]}" for i in range(len(table))]


def _stable_seed(*parts: Any) -> int:
    blob = json.dumps(parts, sort_keys=True, default=str).encode()
    return int(hashlib.sha256(blob).hexdigest()[:16], 16)


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_plain(v) for v in value), key=repr)
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    return value


def _fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(_plain(value), sort_keys=True,
                                    separators=(",", ":")).encode()).hexdigest()


# ------------------------------------------------------------------ splits
@dataclass(frozen=True)
class Fold:
    fold_id: int
    test_idx: Tuple[int, ...]
    train_idx: Tuple[int, ...]
    purged_idx: Tuple[int, ...]              # share a group with test; in neither side
    test_group_ids: Tuple[str, ...]          # held-out primary-axis units
    held_out_groups: Mapping[str, Tuple[str, ...]]  # per axis, every group touching test


@dataclass(frozen=True)
class SplitPlan:
    mode: str
    primary_axis: str
    axes: Tuple[str, ...]
    relaxed_axes: Tuple[str, ...]
    seed: int
    n_independent_groups: int                # connected components over all active axes
    observation_ids: Tuple[str, ...]
    group_ids: Mapping[str, Tuple[str, ...]]
    folds: Tuple[Fold, ...]
    config: GroupingConfig = GroupingConfig()
    n_folds: Optional[int] = None
    min_independent_groups: int = 3
    data_fingerprint: str = ""
    config_fingerprint: str = ""


@dataclass(frozen=True)
class NestedSplitPlan:
    outer: SplitPlan
    inner: Tuple[Tuple[Fold, ...], ...]      # inner[k] folds live inside outer.folds[k].train_idx
    inner_mode: str


def verify_no_leakage(plan: SplitPlan, rows: Any = None) -> None:
    """Raise ``LeakageError`` if any active-axis group appears in both train and test."""
    if rows is not None:
        table = _as_table(rows)
        if plan.observation_ids != table.ids or plan.data_fingerprint != _fingerprint(table.rows):
            raise SchemaError("split plan was built for a different table (data fingerprint differs)")
        expected = make_splits(table, plan.mode, config=plan.config, seed=plan.seed,
                               n_folds=plan.n_folds,
                               min_independent_groups=plan.min_independent_groups)
        if plan != expected:
            raise LeakageError("split plan differs from recomputed groups/config/folds")
    for fold in plan.folds:
        if set(fold.test_idx) & set(fold.train_idx):
            raise LeakageError(f"fold {fold.fold_id}: a row is in both train and test")
        for axis in plan.axes:
            groups = plan.group_ids[axis]
            both = {groups[i] for i in fold.train_idx} & {groups[i] for i in fold.test_idx}
            if both:
                raise LeakageError(f"fold {fold.fold_id}: {axis} groups straddle: {sorted(both)}")


def make_splits(rows: Any, mode: str, *, config: GroupingConfig = GroupingConfig(),
                seed: int = 0, n_folds: Optional[int] = None,
                min_independent_groups: int = 3) -> SplitPlan:
    """Purged group split.

    ``n_folds=None``: leave-one-primary-unit-out (one ASO cluster / study / region
    block per fold, sorted order). ``n_folds=k``: seeded k-fold over independent
    components. Refuses (``SplitRefused``) when the number of independent components
    over all active axes is below ``min_independent_groups``.
    """
    _positive_int(min_independent_groups, "min_independent_groups")
    if mode == "nested":
        raise SchemaError("use make_nested_splits for mode='nested'")
    if mode not in MODES:
        raise SchemaError(f"unknown mode {mode!r}; expected {sorted(MODES)} or 'nested'")
    primary = MODES[mode]
    axes = config.axes
    if primary not in axes:
        raise SchemaError(f"mode {mode!r} cannot relax its own primary axis {primary!r}")
    table = _as_table(rows)
    groups = axis_groups(table, config, axes)
    return _splits_from_groups(table, mode, config, seed, n_folds, min_independent_groups, groups)


def _splits_from_groups(table: ObservationTable, mode: str, config: GroupingConfig,
                        seed: int, n_folds: Optional[int], min_independent_groups: int,
                        groups: Mapping[str, Tuple[str, ...]]) -> SplitPlan:
    _positive_int(min_independent_groups, "min_independent_groups")
    if mode not in MODES or MODES[mode] not in config.axes:
        raise SchemaError(f"invalid split mode or relaxed primary axis: {mode!r}")
    primary, axes = MODES[mode], config.axes
    comp = _components(table, [groups[a] for a in axes])
    n_comp = len(set(comp))
    if n_comp < min_independent_groups:
        raise SplitRefused(
            f"{mode}: independent N too small - {n_comp} independent group(s) over axes "
            f"{list(axes)} from {len(table)} rows, need >= {min_independent_groups}. "
            "Rows sharing a sequence cluster, chemistry family, study or overlapping "
            "window are one unit; relaxing an axis must be explicit (GroupingConfig.relax).")
    n = len(table)
    if n_folds is None:
        units = sorted(set(groups[primary]))
        tests = [[i for i in range(n) if groups[primary][i] == u] for u in units]
    else:
        if not 2 <= n_folds <= n_comp:
            raise SplitRefused(
                f"{mode}: n_folds={n_folds} impossible with {n_comp} independent group(s)")
        order = sorted(set(comp))
        random.Random(_stable_seed("kfold", mode, seed)).shuffle(order)
        bucket = {c: k % n_folds for k, c in enumerate(order)}
        tests = [[i for i in range(n) if bucket[comp[i]] == k] for k in range(n_folds)]
    folds = []
    for fold_id, test in enumerate(tests):
        test_set = set(test)
        held = {a: frozenset(groups[a][i] for i in test) for a in axes}
        train, purged = [], []
        for i in range(n):
            if i in test_set:
                continue
            (purged if any(groups[a][i] in held[a] for a in axes) else train).append(i)
        if not train:
            raise SplitRefused(
                f"{mode}: fold {fold_id} has no training rows after purging groups shared "
                "with the test unit")
        folds.append(Fold(fold_id, tuple(test), tuple(train), tuple(purged),
                          tuple(sorted({groups[primary][i] for i in test})),
                          _freeze({a: tuple(sorted(held[a])) for a in axes})))
    plan = SplitPlan(mode, primary, axes, tuple(config.relax), seed, n_comp, table.ids,
                     _freeze(groups), tuple(folds), config, n_folds, min_independent_groups,
                     _fingerprint(table.rows), _fingerprint({
                         "grouping": asdict(config), "mode": mode, "seed": seed,
                         "n_folds": n_folds, "min_independent_groups": min_independent_groups}))
    verify_no_leakage(plan)
    return plan


def make_nested_splits(rows: Any, *, outer_mode: str, inner_mode: Optional[str] = None,
                       config: GroupingConfig = GroupingConfig(), seed: int = 0,
                       outer_n_folds: Optional[int] = None, inner_n_folds: Optional[int] = None,
                       min_outer_groups: int = 3, min_inner_groups: int = 3) -> NestedSplitPlan:
    """Nested purged splits. Inner folds are built only from each outer fold's train
    rows, so they never contain an outer test (or purged) row or any of its groups.
    Refuses, naming the outer fold, when either level has too few independent groups."""
    _positive_int(min_outer_groups, "min_outer_groups")
    _positive_int(min_inner_groups, "min_inner_groups")
    table = _as_table(rows)
    inner_mode = inner_mode or outer_mode
    outer = make_splits(table, outer_mode, config=config, seed=seed, n_folds=outer_n_folds,
                        min_independent_groups=min_outer_groups)
    inner_all = []
    for fold in outer.folds:
        index = fold.train_idx
        try:
            # Keep full-data group identities: removing a bridge cannot undo a cluster.
            inherited = {a: tuple(g[i] for i in index) for a, g in outer.group_ids.items()}
            sub = _splits_from_groups(
                table.subset(index), inner_mode, config,
                _stable_seed("inner", seed, fold.fold_id) % (2 ** 31),
                inner_n_folds, min_inner_groups, inherited)
        except SplitRefused as exc:
            raise SplitRefused(
                f"nested refused: outer fold {fold.fold_id} (held out "
                f"{list(fold.test_group_ids)}) leaves too little for inner folds - {exc}"
            ) from None
        inner = tuple(
            Fold(f.fold_id, tuple(index[i] for i in f.test_idx),
                 tuple(index[i] for i in f.train_idx), tuple(index[i] for i in f.purged_idx),
                 f.test_group_ids, f.held_out_groups)
            for f in sub.folds)
        excluded = set(fold.test_idx) | set(fold.purged_idx)
        for f in inner:
            if (set(f.test_idx) | set(f.train_idx)) & excluded:
                raise LeakageError(f"inner fold of outer {fold.fold_id} sees outer test rows")
            for axis, g in outer.group_ids.items():
                if {g[i] for i in f.train_idx} & {g[i] for i in f.test_idx}:
                    raise LeakageError(f"inner fold splits original {axis} groups")
                if {g[i] for i in f.train_idx + f.test_idx} & {g[i] for i in fold.test_idx}:
                    raise LeakageError(f"inner fold sees outer test {axis} groups")
        inner_all.append(inner)
    return NestedSplitPlan(outer, tuple(inner_all), inner_mode)


# ------------------------------------------------------------------ metrics
@dataclass(frozen=True)
class MetricPolicy:
    """Minimum counts of independent connected components contributing to each class."""
    min_verified_negatives: int = 10
    min_positives: int = 10

    def __post_init__(self) -> None:
        _positive_int(self.min_verified_negatives, "min_verified_negatives")
        _positive_int(self.min_positives, "min_positives")


@dataclass(frozen=True)
class MetricResult:
    name: str
    value: float
    evidence_class: str
    mode: str
    n_positive_rows: int
    n_negative_rows: int
    n_positive_units: int
    n_negative_units: int
    excluded_by_state: Mapping[str, int]
    n_tied_rows: int
    threshold: Optional[float] = None


_ALIASES = {"auc": "auroc", "roc_auc": "auroc", "roc": "auroc", "average_precision": "auprc",
            "pr_auc": "auprc", "ppv": "precision", "sensitivity": "recall", "tpr": "recall",
            "ece": "calibration_error", "calibration": "calibration_error",
            "false_discovery_rate": "fdr"}
SUPERVISED_METRICS = ("auroc", "auprc", "accuracy", "precision", "recall", "fdr", "brier",
                      "calibration_error")
_THRESHOLD_METRICS = ("accuracy", "precision", "recall", "fdr")


def metric_key(name: str) -> str:
    key = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    for prefix in ("sklearn.metrics.", "metrics."):
        if key.startswith(prefix):
            key = key[len(prefix):]
    key = key.replace("_score", "").strip("_")
    return _ALIASES.get(key, key)


def _resolve_scores(table: ObservationTable, scores: Any) -> List[float]:
    if isinstance(scores, Mapping):
        missing = [oid for oid in table.ids if oid not in scores]
        if missing:
            raise SchemaError(f"scores missing for {len(missing)} rows, e.g. {missing[:3]}")
        out = [scores[oid] for oid in table.ids]
    else:
        out = list(scores)
        if len(out) != len(table):
            raise SchemaError(f"{len(out)} scores for {len(table)} rows")
    for value in out:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise SchemaError(f"score {value!r} is not a finite number")
    return [float(v) for v in out]


def _midranks(values: Sequence[float]) -> List[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return ranks


def _rank_probability(pos: Sequence[float], other: Sequence[float]) -> float:
    """P(pos > other) + 0.5 P(tie), tie-exact (Mann-Whitney U / (n_pos * n_other))."""
    ranks = _midranks(list(pos) + list(other))
    u = sum(ranks[:len(pos)]) - len(pos) * (len(pos) + 1) / 2
    return u / (len(pos) * len(other))


def _pr_points(y: Sequence[int], s: Sequence[float]) -> List[Tuple[float, float, float]]:
    """(threshold, precision, recall) at each distinct score, descending; ties grouped."""
    n_pos = sum(y)
    order = sorted(range(len(s)), key=lambda i: -s[i])
    points, tp, fp, i = [], 0, 0, 0
    while i < len(order):
        j = i
        while j < len(order) and s[order[j]] == s[order[i]]:
            tp += y[order[j]]
            fp += 1 - y[order[j]]
            j += 1
        points.append((s[order[i]], tp / (tp + fp), tp / n_pos))
        i = j
    return points


def _metric_value(key: str, y: Sequence[int], s: Sequence[float],
                  threshold: Optional[float], n_bins: int) -> float:
    n_pos = sum(y)
    if key == "auroc":
        return _rank_probability([v for v, t in zip(s, y) if t], [v for v, t in zip(s, y) if not t])
    if key == "auprc":  # average precision; a tie group contributes one step
        total, prev = 0.0, 0.0
        for _, precision, recall in _pr_points(y, s):
            total += (recall - prev) * precision
            prev = recall
        return total
    if key in _THRESHOLD_METRICS:
        if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) \
                or not math.isfinite(threshold):
            raise SchemaError(f"{key} needs an explicit finite non-boolean numeric threshold "
                              "(score >= threshold is positive)")
        pred = [v >= threshold for v in s]
        tp = sum(1 for p, t in zip(pred, y) if p and t)
        fp = sum(1 for p, t in zip(pred, y) if p and not t)
        if key == "accuracy":
            return (tp + sum(1 for p, t in zip(pred, y) if not p and not t)) / len(y)
        if key == "recall":
            return tp / n_pos
        if tp + fp == 0:
            raise MetricBlocked(key, f"no predictions at threshold {threshold}; {key} undefined")
        return tp / (tp + fp) if key == "precision" else fp / (tp + fp)
    if any(v < 0.0 or v > 1.0 for v in s):
        raise SchemaError(f"{key} needs probabilities in [0, 1]")
    if key == "brier":
        return sum((v - t) ** 2 for v, t in zip(s, y)) / len(y)
    if n_bins < 1:
        raise SchemaError("n_bins must be >= 1")
    bins: Dict[int, List[Tuple[float, int]]] = {}
    for v, t in zip(s, y):
        bins.setdefault(min(int(v * n_bins), n_bins - 1), []).append((v, t))
    return sum(abs(sum(v for v, _ in b) - sum(t for _, t in b)) for b in bins.values()) / len(y)


def _supervised(key: str, display: str, states: Sequence[str], units: Sequence[str],
                scores: Sequence[float], policy: MetricPolicy, threshold: Optional[float],
                n_bins: int) -> MetricResult:
    policy.__post_init__()  # public policy instances are not a validation bypass
    keep = [i for i, st in enumerate(states) if st in (POSITIVE, NEGATIVE)]
    y = [1 if states[i] == POSITIVE else 0 for i in keep]
    s = [scores[i] for i in keep]
    pos_units = {units[i] for i in keep if states[i] == POSITIVE}
    neg_units = {units[i] for i in keep if states[i] == NEGATIVE}
    excluded: Dict[str, int] = {}
    for st in states:
        if st not in (POSITIVE, NEGATIVE):
            excluded[st] = excluded.get(st, 0) + 1
    problems = []
    if not neg_units or len(neg_units) < policy.min_verified_negatives:
        problems.append(f"{len(neg_units)} independent verified_negative unit(s) "
                        f"({len(y) - sum(y)} rows) < minimum {policy.min_verified_negatives}")
    if not pos_units or len(pos_units) < policy.min_positives:
        problems.append(f"{len(pos_units)} independent measured_positive unit(s) "
                        f"({sum(y)} rows) < minimum {policy.min_positives}")
    if problems:
        raise MetricBlocked(display, "; ".join(problems) + (
            f"; excluded (never negatives): {excluded}" if excluded else "")
            + ". Use pu_report() for positive-unlabeled quantities.")
    value = _metric_value(key, y, s, threshold, n_bins)
    return MetricResult(key, value, "derived", "supervised", sum(y), len(y) - sum(y),
                        len(pos_units), len(neg_units), excluded, len(s) - len(set(s)),
                        threshold if key in _THRESHOLD_METRICS else None)


def _check_supervised_request(name: str, mode: str, negative_states: Sequence[str]) -> str:
    key = metric_key(name)
    if mode == "pu":
        raise MetricBlocked(name, "positive-unlabeled mode never emits supervised metrics; "
                                  "use pu_report()")
    if mode != "supervised":
        raise SchemaError(f"unknown mode {mode!r}")
    extra = set(negative_states) - {NEGATIVE}
    if extra or NEGATIVE not in negative_states:
        raise LabelStateError(
            f"negative_states={list(negative_states)}: only verified_negative is a negative; "
            f"{sorted(extra)} can never be treated as negative")
    if key not in SUPERVISED_METRICS:
        raise SchemaError(f"unknown metric {name!r}; supported: {list(SUPERVISED_METRICS)}")
    return key


def _has_duplicate_contributors(table: ObservationTable, indices: Sequence[int]) -> bool:
    ids = {table.rows[i]["observation_id"] for i in indices}
    return any(sum(oid in ids for oid in copies) > 1 for copies in table.duplicate_assay_records)


_DUPLICATE_REASON = ("duplicate assay records: renamed copies provide only 1 independent "
                     "assay per group; resolve copies before computing statistics")


def compute_metric(name: str, rows: Any, scores: Any, *, policy: MetricPolicy = MetricPolicy(),
                   threshold: Optional[float] = None, n_bins: int = 10,
                   mode: str = "supervised",
                   negative_states: Sequence[str] = (NEGATIVE,)) -> MetricResult:
    """Guarded supervised metric over measured_positive vs verified_negative rows only.

    Higher score = more positive. Raises ``MetricBlocked`` (with ``.reason``) when
    either class has fewer independent units than ``policy`` requires.
    """
    key = _check_supervised_request(name, mode, negative_states)
    table = _as_table(rows)
    keep = [i for i, r in enumerate(table.rows) if r["label_state"] in (POSITIVE, NEGATIVE)]
    if _has_duplicate_contributors(table, keep):
        raise MetricBlocked(name, _DUPLICATE_REASON)
    s = _resolve_scores(table, scores)
    return _supervised(key, name, [r["label_state"] for r in table.rows],
                       effective_clusters(table), s, policy, threshold,
                       n_bins)


def precision_recall_curve(rows: Any, scores: Any, *, policy: MetricPolicy = MetricPolicy()
                           ) -> List[Tuple[float, float, float]]:
    """Guarded (threshold, precision, recall) points; same blocking rule as AUPRC."""
    table = _as_table(rows)
    s = _resolve_scores(table, scores)
    compute_metric("auprc", table, s, policy=policy)  # raises MetricBlocked if not allowed
    keep = [i for i, r in enumerate(table.rows) if r["label_state"] in (POSITIVE, NEGATIVE)]
    return _pr_points([int(table.rows[i]["label_state"] == POSITIVE) for i in keep],
                      [s[i] for i in keep])


# ------------------------------------------------------------------ positive-unlabeled mode
PU_QUANTITIES = ("pu_recall_at_top_k", "pu_positive_mean_percentile_rank",
                 "pu_positive_median_percentile_rank",
                 "pu_positive_vs_unlabeled_rank_probability",
                 "pu_prior_adjusted_rank_probability")
_PU_BASE = ("positive-unlabeled: unlabeled candidates are a mixture of undiscovered "
            "positives and negatives and are never counted as negatives")


@dataclass(frozen=True)
class PUQuantity:
    name: str
    value: float
    evidence_class: str
    pu_assumption: str
    k: Optional[int] = None


@dataclass(frozen=True)
class PUReport:
    mode: str
    quantities: Tuple[PUQuantity, ...]
    n_positive_rows: int
    n_positive_units: int
    n_universe: int
    universe_states: Tuple[str, ...]
    excluded_by_state: Mapping[str, int]
    class_prior: Optional[float]
    blocked: Mapping[str, str]

    def get(self, name: str, k: Optional[int] = None) -> PUQuantity:
        for q in self.quantities:
            if q.name == name and q.k == k:
                return q
        raise KeyError((name, k))


def pu_report(rows: Any, scores: Any, *, top_k: Sequence[int] = (), class_prior: Optional[float] = None,
              min_positives: int = 10, include_unassayed: bool = False) -> PUReport:
    """PU-valid quantities only. The candidate universe is measured_positive + unlabeled
    (+ unassayed if ``include_unassayed``); verified_negative and censored rows are excluded.

    ``class_prior`` is a STATED assumption (fraction of true positives among the
    non-positive universe rows), never estimated here; the prior-adjusted statistic is
    only emitted when it is given and is tagged ``model_assumption``.
    """
    _positive_int(min_positives, "min_positives")
    table = _as_table(rows)
    s_all = _resolve_scores(table, scores)
    states = [LabelState.UNLABELED.value] + ([LabelState.UNASSAYED.value] if include_unassayed else [])
    pos = [i for i, r in enumerate(table.rows) if r["label_state"] == POSITIVE]
    unl = [i for i, r in enumerate(table.rows) if r["label_state"] in states]
    if _has_duplicate_contributors(table, pos + unl):
        raise MetricBlocked("pu_report", _DUPLICATE_REASON)
    excluded: Dict[str, int] = {}
    for r in table.rows:
        if r["label_state"] != POSITIVE and r["label_state"] not in states:
            excluded[r["label_state"]] = excluded.get(r["label_state"], 0) + 1
    clusters = effective_clusters(table)
    pos_units = {clusters[i] for i in pos}
    if len(pos_units) < min_positives:
        raise MetricBlocked("pu_report", f"{len(pos_units)} independent measured_positive "
                            f"unit(s) ({len(pos)} rows) < minimum {min_positives}")
    if class_prior is not None and not 0.0 <= class_prior < 1.0:
        raise SchemaError("class_prior must be in [0, 1)")
    universe = pos + unl
    s = [s_all[i] for i in universe]
    is_pos = [1] * len(pos) + [0] * len(unl)
    n = len(universe)
    out: List[PUQuantity] = []
    blocked: Dict[str, str] = {}
    universe_note = f"{_PU_BASE}; universe = {n} rows ({' + '.join([POSITIVE] + states)})"

    order = sorted(range(n), key=lambda i: -s[i])
    for k in top_k:
        if isinstance(k, bool) or not isinstance(k, int) or not 1 <= k <= n:
            raise SchemaError(f"top_k entry {k!r} must be an int in [1, {n}]")
        hits, taken, i = 0.0, 0, 0
        while i < n and taken < k:  # a tie group straddling k gets fractional credit
            j = i
            while j < n and s[order[j]] == s[order[i]]:
                j += 1
            group_pos = sum(is_pos[order[m]] for m in range(i, j))
            hits += group_pos * min(1.0, (k - taken) / (j - i))
            taken += j - i
            i = j
        out.append(PUQuantity(
            "pu_recall_at_top_k", hits / len(pos), "derived",
            f"{universe_note}; recall of KNOWN positives only, ties at the cut-off get "
            "expected credit under random tie-breaking; says nothing about precision", k))

    ranks = _midranks(s)
    pct = sorted((ranks[i] - 0.5) / n for i in range(len(pos)))
    mid = len(pct) // 2
    median = pct[mid] if len(pct) % 2 else (pct[mid - 1] + pct[mid]) / 2
    note = f"{universe_note}; percentile of mid-rank within the universe (1 = top)"
    out.append(PUQuantity("pu_positive_mean_percentile_rank", sum(pct) / len(pct), "derived", note))
    out.append(PUQuantity("pu_positive_median_percentile_rank", median, "derived", note))

    if not unl:
        blocked["pu_positive_vs_unlabeled_rank_probability"] = "no unlabeled rows in the universe"
    else:
        prob = _rank_probability(s[:len(pos)], s[len(pos):])
        out.append(PUQuantity(
            "pu_positive_vs_unlabeled_rank_probability", prob, "derived",
            f"{universe_note}; P(positive outranks an unlabeled row) with ties = 1/2; NOT a "
            "ranking statistic against negatives - hidden positives pull it toward 0.5"))
        if class_prior is None:
            blocked["pu_prior_adjusted_rank_probability"] = "no class_prior stated"
        else:
            out.append(PUQuantity(
                "pu_prior_adjusted_rank_probability",
                (prob - class_prior / 2) / (1 - class_prior), "model_assumption",
                f"{universe_note}; assumes labelled positives are selected completely at "
                f"random (SCAR) and a STATED, not estimated, class prior of {class_prior} "
                "among unlabeled rows; may leave [0, 1] if the prior is wrong"))
    return PUReport("positive_unlabeled", tuple(out), len(pos), len(pos_units), n,
                    tuple([POSITIVE] + states), excluded, class_prior, blocked)


# ------------------------------------------------------------------ clustered uncertainty
@dataclass(frozen=True)
class CIResult:
    point: float
    lower: float
    upper: float
    alpha: float
    n_clusters: int
    n_boot: int
    n_failed: int
    seed: int
    cluster_axes: Tuple[str, ...]
    method: str
    evidence_class: str

    def __post_init__(self) -> None:
        if not all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)
                   for v in (self.point, self.lower, self.upper)):
            raise UncertaintyRefused("CI point and bounds must be finite")


def effective_clusters(rows: Any, *, merge_axes: Sequence[str] = ("sequence", "window"),
                       config: GroupingConfig = GroupingConfig()) -> Tuple[str, ...]:
    """Independent cluster per row: ``independent_unit_id`` merged with every listed
    leakage axis, so overlapping windows / near-identical ASOs filed under different
    unit ids still count as ONE cluster. Exact renamed assay copies always collapse,
    including when optional merge_axes are disabled."""
    table = _as_table(rows)
    lists: List[Sequence[str]] = [[f"unit:{r['independent_unit_id']}" for r in table.rows]]
    duplicate_id = {oid: ids[0] for ids in table.duplicate_assay_records for oid in ids}
    lists.append([duplicate_id.get(oid, oid) for oid in table.ids])
    lists += list(axis_groups(table, config, tuple(merge_axes)).values())
    return tuple(_components(table, lists))


def _quantile(sorted_values: Sequence[float], q: float) -> float:
    pos = q * (len(sorted_values) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(sorted_values) - 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (pos - lo)


def cluster_bootstrap_ci(rows: Any, statistic: Callable[[Sequence[int]], float], *,
                         n_boot: int = 1000, seed: int = 0, min_clusters: int = 10,
                         alpha: float = 0.05, merge_axes: Sequence[str] = ("sequence", "window"),
                         config: GroupingConfig = GroupingConfig(),
                         max_failed_fraction: float = 0.1,
                         eligible_indices: Optional[Sequence[int]] = None) -> CIResult:
    """Percentile cluster bootstrap. ``statistic(indices)`` receives row indices (with
    repeats) of whole resampled clusters. Refuses with ``UncertaintyRefused`` when
    fewer than ``min_clusters`` effective clusters exist, or when too many replicates
    are degenerate (``MetricBlocked`` / non-finite). ``eligible_indices`` restricts
    the point/resamples and cluster count to contributing observations, while keeping
    the original clustering (excluded bridge rows cannot split an established unit).
    Statistic indices always refer to the original input rows.
    Sequence and window merging are mandatory; extra merge axes may be requested.
    Exact assay copies among contributing rows are refused before calling statistic.
    """
    merge_axes = tuple(merge_axes)
    if not {"sequence", "window"}.issubset(merge_axes):
        raise UncertaintyRefused("CI requires sequence and window merging; axes cannot be weakened")
    table = _as_table(rows)
    _positive_int(n_boot, "n_boot")
    _positive_int(min_clusters, "min_clusters")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) \
            or not math.isfinite(alpha) or not 0.0 < alpha < 1.0:
        raise SchemaError("need finite 0 < alpha < 1")
    if isinstance(max_failed_fraction, bool) or not isinstance(max_failed_fraction, (int, float)) \
            or not math.isfinite(max_failed_fraction) or not 0 <= max_failed_fraction <= 1:
        raise SchemaError("max_failed_fraction must be finite in [0, 1]")
    clusters = effective_clusters(table, merge_axes=merge_axes, config=config)
    eligible = tuple(range(len(table))) if eligible_indices is None else tuple(eligible_indices)
    if any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(table)
           for i in eligible) or len(set(eligible)) != len(eligible):
        raise SchemaError("eligible_indices must be unique valid row indices")
    if _has_duplicate_contributors(table, eligible):
        raise UncertaintyRefused(_DUPLICATE_REASON)
    members: Dict[str, List[int]] = {}
    for i in eligible:
        members.setdefault(clusters[i], []).append(i)
    ids = sorted(members)
    if len(ids) < min_clusters:
        raise UncertaintyRefused(
            f"CI refused: {len(ids)} independent cluster(s) from {len(eligible)} rows "
            f"(units merged over {list(merge_axes)}) < minimum {min_clusters}; correlated "
            "rows are not independent replicates")
    point = float(statistic(eligible))
    if not math.isfinite(point):
        raise UncertaintyRefused("CI refused: point estimate is not finite")
    rng = random.Random(_stable_seed("cluster_bootstrap", seed))
    values, failed = [], 0
    for _ in range(n_boot):
        draw = [ids[rng.randrange(len(ids))] for _ in ids]
        try:
            value = float(statistic(tuple(i for c in draw for i in members[c])))
        except MetricBlocked:
            failed += 1
            continue
        if math.isfinite(value):
            values.append(value)
        else:
            failed += 1
    if not values or failed > max_failed_fraction * n_boot:
        raise UncertaintyRefused(
            f"CI refused: {failed}/{n_boot} bootstrap replicates were degenerate "
            f"(> {max_failed_fraction:.0%}); too few clusters per class for a stable interval")
    values.sort()
    return CIResult(point, _quantile(values, alpha / 2), _quantile(values, 1 - alpha / 2), alpha,
                    len(ids), n_boot, failed, seed, tuple(merge_axes),
                    "percentile_cluster_bootstrap", "derived")


def bootstrap_metric_ci(name: str, rows: Any, scores: Any, *,
                        policy: MetricPolicy = MetricPolicy(), threshold: Optional[float] = None,
                        n_bins: int = 10, **kwargs: Any) -> CIResult:
    """Cluster-bootstrap CI for a guarded metric. The point estimate must pass the full
    ``policy``; replicates only need both classes present (else counted as failed)."""
    key = _check_supervised_request(name, "supervised", (NEGATIVE,))
    table = _as_table(rows)
    s = _resolve_scores(table, scores)
    compute_metric(name, table, s, policy=policy, threshold=threshold, n_bins=n_bins)
    keep = [i for i, r in enumerate(table.rows) if r["label_state"] in (POSITIVE, NEGATIVE)]
    states = [r["label_state"] for r in table.rows]
    units = effective_clusters(table)
    relaxed = MetricPolicy(1, 1)

    def statistic(indices: Sequence[int]) -> float:
        return _supervised(key, name, [states[i] for i in indices], [units[i] for i in indices],
                           [s[i] for i in indices], relaxed, threshold, n_bins).value

    return cluster_bootstrap_ci(table, statistic, eligible_indices=keep, **kwargs)


# ------------------------------------------------------------------ baselines / ablation
@dataclass(frozen=True)
class FeatureFamily:
    name: str
    columns: Tuple[str, ...]
    evidence_class: str
    introduced_by: str


@dataclass(frozen=True)
class BaselineSpec:
    name: str
    kind: str                                 # permuted_reference | feature_set | undefined
    families: Tuple[FeatureFamily, ...]
    description: str

    @property
    def feature_columns(self) -> Tuple[str, ...]:
        return tuple(c for fam in self.families for c in fam.columns)


# Numeric J9 ``feature_name`` values (long-format ``value_double`` / ``value_int``).
_FAMILIES = (
    FeatureFamily("complementarity", (
        "m1_aso_length", "m1_wc_pairs", "m1_gu_pairs", "m1_mismatch_count", "m1_wc_fraction",
        "m1_n_sites_in_locus_window", "m1_n_exact_sites_in_locus_window"), "derived", "M1"),
    FeatureFamily("duplex", ("m2_duplex_energy_kcal_mol",), "model_assumption", "M2"),
    FeatureFamily("accessibility", ("m3_accessibility_p_unpaired",), "proxy", "M3"),
    FeatureFamily("splice_geometry", (
        "m4_abs_donor_dist_5p_from_last_exonic", "m4_abs_acceptor_dist_3p_from_first_exonic",
        "m4_target_span_nt"), "derived", "M4"),
    FeatureFamily("sequence_ablation_proxy", (
        "m5_raw_model_delta_mean", "m5_raw_model_delta_max"), "proxy", "M5"),
)


def _build_registry() -> Dict[str, BaselineSpec]:
    registry = {"M0": BaselineSpec("M0", "permuted_reference", (),
                                   "seeded permuted / random-order reference; uses no features")}
    for k, fam in enumerate(_FAMILIES, 1):
        registry[fam.introduced_by] = BaselineSpec(
            fam.introduced_by, "feature_set", _FAMILIES[:k],
            f"cumulative J9 baseline: families {[f.name for f in _FAMILIES[:k]]}")
    registry["M6"] = BaselineSpec(
        "M6", "undefined", (), "reserved: no M6 feature family is defined upstream (J9 stops "
        "at M5); register one with register_baseline() before use")
    return registry


BASELINES: Dict[str, BaselineSpec] = _build_registry()


def register_baseline(spec: BaselineSpec) -> None:
    """Define / replace a baseline (e.g. M6 once its feature family is agreed)."""
    if spec.kind not in ("permuted_reference", "feature_set"):
        raise SchemaError(f"bad baseline kind {spec.kind!r}")
    for fam in spec.families:
        if fam.evidence_class not in EVIDENCE_CLASSES:
            raise SchemaError(f"bad evidence class {fam.evidence_class!r}")
    BASELINES[spec.name] = spec


def baseline_feature_columns(name: str) -> Tuple[str, ...]:
    if name not in BASELINES:
        raise SchemaError(f"unknown baseline {name!r}; known: {sorted(BASELINES)}")
    spec = BASELINES[name]
    if spec.kind == "undefined":
        raise BaselineUndefined(f"{name}: {spec.description}")
    return spec.feature_columns


def permuted_reference_scores(observation_ids: Sequence[str], seed: int = 0) -> Dict[str, float]:
    """M0: seeded random scores, independent of input order."""
    observation_ids = tuple(
        _identifier(oid, "observation_id", "permuted_reference_scores") for oid in observation_ids)
    rng = random.Random(_stable_seed("m0_permuted_reference", seed))
    return {oid: rng.random() for oid in sorted(observation_ids)}


@dataclass(frozen=True)
class AblationCell:
    variant: str                              # "full", "drop:<family>", "permuted_reference"
    dropped_family: Optional[str]
    feature_columns: Tuple[str, ...]
    fold_id: Optional[int]                    # None = pooled out-of-fold scores
    test_group_ids: Tuple[str, ...]
    held_out_groups: Mapping[str, Tuple[str, ...]]
    n_train: int
    n_test: int
    n_purged: int
    metrics: Tuple[MetricResult, ...]
    pu: Optional[PUReport]
    blocked: Mapping[str, str]                # metric name -> reason


@dataclass(frozen=True)
class AblationResult:
    baseline: str
    mode: str                                 # "supervised" | "positive_unlabeled"
    split_mode: str
    relaxed_axes: Tuple[str, ...]
    seed: int
    cells: Tuple[AblationCell, ...]           # per variant x fold
    pooled: Tuple[AblationCell, ...]          # per variant, out-of-fold scores pooled
    pooled_note: str = ("pooled cells combine scores from different per-fold models; "
                        "they are comparable only if the model's score scale is")

    def blocked_reasons(self) -> Dict[Tuple[str, Optional[int], str], str]:
        return {(c.variant, c.fold_id, m): r for c in self.cells + self.pooled
                for m, r in c.blocked.items()}


FitPredict = Callable[[List[Dict[str, Any]], List[Dict[str, Any]], Tuple[str, ...]],
                      Mapping[str, float]]


def _evaluate(table: ObservationTable, scores: Mapping[str, float], pu: bool,
              metrics: Sequence[str], policy: MetricPolicy, threshold: Optional[float],
              top_k: Sequence[int], class_prior: Optional[float]
              ) -> Tuple[Tuple[MetricResult, ...], Optional[PUReport], Dict[str, str]]:
    blocked: Dict[str, str] = {}
    if pu:
        n_universe = sum(r["label_state"] in (POSITIVE, LabelState.UNLABELED.value)
                         for r in table.rows)
        for k in top_k:
            _positive_int(k, "top_k entry")
            if k > n_universe:
                blocked[f"pu_recall_at_top_k[{k}]"] = f"k={k} exceeds universe of {n_universe}"
        try:
            report = pu_report(table, scores, top_k=[k for k in top_k if k <= n_universe],
                               class_prior=class_prior, min_positives=policy.min_positives)
        except MetricBlocked as exc:
            return (), None, {**blocked, exc.metric: exc.reason}
        blocked.update(report.blocked)
        return (), report, blocked
    results = []
    for name in metrics:
        try:
            results.append(compute_metric(name, table, scores, policy=policy, threshold=threshold))
        except MetricBlocked as exc:
            blocked[name] = exc.reason
    return tuple(results), None, blocked


def _validate_feature_sources(features: Mapping[str, Mapping[str, Any]],
                              ids: Sequence[str], columns: Sequence[str],
                              source_paths: Sequence[Any],
                              feature_sources: Any) -> Mapping[str, Mapping[str, Any]]:
    """Require declared lineage for every feature cell before passing it to a fitter.

    Manifest shape: {observation_id: {column: (source_path, ...)}}. This validates
    coverage and forbidden names; it cannot authenticate caller-declared lineage.
    """
    if not source_paths or not isinstance(feature_sources, Mapping):
        raise SchemaError("feature inputs require nonempty source_paths and feature_sources manifest")
    allowed = set()
    for path in source_paths:
        if not isinstance(path, (str, os.PathLike)) or not os.fspath(path).strip():
            raise SchemaError("source_paths entries must be nonempty paths")
        allowed.add(check_source_not_forbidden(path))
    # Check ALL declared sources, including entries not consumed by this baseline.
    for oid, manifest in feature_sources.items():
        if not isinstance(manifest, Mapping):
            raise SchemaError(f"feature_sources[{oid!r}] must map columns to source paths")
        for col, paths in manifest.items():
            if not isinstance(paths, (list, tuple)) or not paths:
                raise SchemaError(f"feature source manifest for {oid}/{col} needs nonempty paths")
            for path in paths:
                if not isinstance(path, (str, os.PathLike)) or not os.fspath(path).strip():
                    raise SchemaError("feature source paths must be nonempty")
                if check_source_not_forbidden(path) not in allowed:
                    raise SchemaError(f"feature source {path!r} absent from source_paths")
    for oid, values in features.items():
        if not isinstance(values, Mapping):
            raise SchemaError(f"features[{oid!r}] must be a mapping")
        required = set(values) | (set(columns) if oid in ids else set())
        if required - set(feature_sources.get(oid, {})):
            raise SchemaError(f"feature source manifest missing coverage for {oid}")
    return _freeze(features)


def run_ablation(rows: Any, features: Mapping[str, Mapping[str, Any]], baseline: str,
                 plan: SplitPlan, fit_predict: Optional[FitPredict] = None, *,
                 metrics: Sequence[str] = ("auroc", "auprc"), pu: bool = False,
                 policy: MetricPolicy = MetricPolicy(), threshold: Optional[float] = None,
                 top_k: Sequence[int] = (), class_prior: Optional[float] = None,
                 seed: int = 0, source_paths: Sequence[Any] = (),
                 feature_sources: Optional[Mapping[str, Mapping[str, Sequence[Any]]]] = None
                 ) -> AblationResult:
    """Leave-feature-family-out runner (interface only: the caller supplies the model).

    ``fit_predict(train, test, columns) -> {observation_id: score}``. Train records
    carry ``label_state``; test records do NOT (labels cannot leak into prediction).
    M0 ignores ``fit_predict`` and scores with a seeded permutation. Every blocked
    metric is recorded with its reason instead of being dropped.
    ``source_paths``: every input file behind ``rows`` / ``features``; each is checked
    against the forbidden outcome sources before anything runs.
    Feature baselines also require ``feature_sources[observation_id][column]`` to
    list nonempty source paths from that inventory for every feature input.
    """
    for path in source_paths:
        check_source_not_forbidden(path)
    table = _as_table(rows)
    if plan.observation_ids != table.ids:
        raise SchemaError("split plan was built for a different table (observation ids differ)")
    verify_no_leakage(plan, table)
    policy.__post_init__()
    for name in metrics:
        _check_supervised_request(name, "supervised", (NEGATIVE,))
    spec = BASELINES.get(baseline)
    columns = baseline_feature_columns(baseline)  # raises for unknown / undefined
    if spec.kind == "permuted_reference":
        variants: List[Tuple[str, Optional[str], Tuple[str, ...]]] = [("permuted_reference", None, ())]
    else:
        if fit_predict is None:
            raise SchemaError(f"{baseline} needs a fit_predict callable")
        missing = [oid for oid in table.ids if oid not in features]
        if missing:
            raise SchemaError(f"features missing for {len(missing)} rows, e.g. {missing[:3]}")
        features = _validate_feature_sources(features, table.ids, columns,
                                              source_paths, feature_sources)
        variants = [("full", None, columns)] + [
            (f"drop:{fam.name}", fam.name,
             tuple(c for other in spec.families if other is not fam for c in other.columns))
            for fam in spec.families]

    def record(i: int, cols: Tuple[str, ...], with_label: bool) -> Dict[str, Any]:
        row = table.rows[i]
        rec: Dict[str, Any] = {"observation_id": row["observation_id"],
                               "features": {c: features[row["observation_id"]].get(c) for c in cols}}
        if with_label:
            rec["label_state"] = row["label_state"]
        return rec

    cells, pooled = [], []
    for variant, dropped, cols in variants:
        oof: Dict[str, float] = {}
        for fold in plan.folds:
            test_ids = [table.ids[i] for i in fold.test_idx]
            if spec.kind == "permuted_reference":
                scores = permuted_reference_scores(test_ids, _stable_seed(seed, fold.fold_id) % 2 ** 31)
            else:
                scores = dict(fit_predict([record(i, cols, True) for i in fold.train_idx],
                                          [record(i, cols, False) for i in fold.test_idx], cols))
                if set(scores) != set(test_ids):
                    raise SchemaError(f"fit_predict must score exactly the test rows of fold "
                                      f"{fold.fold_id}")
            oof.update(scores)
            sub = table.subset(fold.test_idx)
            res, report, blocked = _evaluate(sub, scores, pu, metrics, policy, threshold, top_k,
                                             class_prior)
            cells.append(AblationCell(variant, dropped, cols, fold.fold_id, fold.test_group_ids,
                                      fold.held_out_groups, len(fold.train_idx),
                                      len(fold.test_idx), len(fold.purged_idx), res, report,
                                      blocked))
        index = [i for i, oid in enumerate(table.ids) if oid in oof]
        res, report, blocked = _evaluate(table.subset(index), oof, pu, metrics, policy, threshold,
                                         top_k, class_prior)
        pooled.append(AblationCell(
            variant, dropped, cols, None,
            tuple(sorted({g for f in plan.folds for g in f.test_group_ids})), {}, 0, len(index),
            0, res, report, blocked))
    return AblationResult(baseline, "positive_unlabeled" if pu else "supervised", plan.mode,
                          plan.relaxed_axes, seed, tuple(cells), tuple(pooled))
