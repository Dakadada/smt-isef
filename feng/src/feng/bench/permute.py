"""Arm D Phase-9 matched-stratum permutation test -- GUARDED SHUT until freeze.

The locked benchmark's confirmatory readout is a one-sided, matched-stratum label
permutation over 10^6 shuffles with a frozen seed, comparing the baseline score
(``min 18-nt Hamming distance``) against the test score (accessibility-aware
IntaRNA interaction energy). The effect size is::

    delta_r = mean over positives of ( percentile_test - percentile_baseline )

where each percentile is the positive's rank within its own matched stratum.

FREEZE DISCIPLINE. Computing delta_r or its permutation p-value inspects the
benchmark OUTCOME and must not happen before preregistration is written and the
git freeze tag is cut (Phase 8), after which the test runs ONCE (Phase 9). This
module therefore refuses to run unless an explicit unlock is present:

    FENG_PHASE9_UNLOCK=1  (environment)   -- or --
    a sentinel file whose path is passed to ``assert_phase9_unlocked``.

Every outcome entry point calls :func:`assert_phase9_unlocked` first. The pure
building block :func:`percentile_within_stratum` is safe to unit-test on
synthetic data; the aggregate / permutation functions are never run pre-freeze
and are intentionally left inert.
"""

from __future__ import annotations

from pathlib import Path
import hashlib
import json
import os
import posixpath
import re
import subprocess

import numpy as np
import pandas as pd

from feng.bench.strata import FreezeGuardError

PHASE9_UNLOCK_ENV = "FENG_PHASE9_UNLOCK"
DEFAULT_N_PERMUTATIONS = 1_000_000
DEFAULT_SEED = 20260901
FENG_REPO_ROOT = Path(__file__).resolve().parents[4]
FREEZE_MANIFEST_PATH = FENG_REPO_ROOT / "feng/artifacts/bench/freeze_manifest.json"
EXPECTED_PREREGISTRATION_PATH = "feng/PREREGISTRATION.md"
EXPECTED_SCORED_STRATA_PATH = "feng/artifacts/bench/scored_strata.parquet"
EXPECTED_ARTIFACT_PATHS = frozenset({
    "feng/artifacts/bench/detectable_universe_v50.parquet",
    "feng/artifacts/bench/detectable_universe_v50_README.md",
    "feng/artifacts/bench/positives_v50.tsv",
    "feng/artifacts/bench/positives_v50_README.md",
    "feng/artifacts/bench/strata_v50.parquet",
    EXPECTED_SCORED_STRATA_PATH,
})
UNIVERSE_POINTER_PATH = Path(
    "feng/artifacts/bench/detectable_universe_v50.current"
)
UNIVERSE_BUNDLE_ROOT_PATH = Path(
    "feng/artifacts/bench/detectable_universe_v50.bundles"
)
POINTER_BACKED_ARTIFACTS = frozenset({
    Path("feng/artifacts/bench/detectable_universe_v50.parquet"),
    Path("feng/artifacts/bench/detectable_universe_v50_README.md"),
})
FREEZE_MANIFEST_SCHEMA = "feng-phase9-freeze-manifest-v1"
APPROVED_FREEZE_TAG = "feng-freeze"
# External trust root. Phase 9 remains impossible until the authorized freeze
# ceremony publishes the annotated tag object's SHA-256 through the immutable
# freeze record and this sentinel is replaced with that approved digest.
APPROVED_FREEZE_TAG_OBJECT_SHA256 = "0" * 64


def _git(args: list[str], *, text: bool = True):
    try:
        return subprocess.run(
            ["git", "-C", str(FENG_REPO_ROOT), *args],
            check=True,
            capture_output=True,
            text=text,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise FreezeGuardError(f"could not verify frozen git state: {exc}") from exc


def _tagged_and_worktree_bytes(tag_commit: str, relative: Path) -> tuple[bytes, bytes]:
    """Read the same logical artifact bytes from the worktree and frozen tag."""
    logical_path = FENG_REPO_ROOT / relative
    if relative not in POINTER_BACKED_ARTIFACTS:
        current_path = logical_path.resolve()
        if not current_path.is_relative_to(FENG_REPO_ROOT.resolve()):
            raise FreezeGuardError("manifest path escapes the repository")
        return (
            current_path.read_bytes(),
            _git(["show", f"{tag_commit}:{relative.as_posix()}"], text=False).stdout,
        )

    expected_logical_target = Path(UNIVERSE_POINTER_PATH.name) / relative.name
    if not logical_path.is_symlink():
        raise FreezeGuardError(f"pointer-backed artifact is not a symlink: {relative}")
    logical_target = os.readlink(logical_path)
    tagged_logical_target = _git(
        ["show", f"{tag_commit}:{relative.as_posix()}"], text=False
    ).stdout
    if (
        logical_target != expected_logical_target.as_posix()
        or tagged_logical_target != logical_target.encode("utf-8")
    ):
        raise FreezeGuardError(f"tag/pointer mismatch for {relative}")

    pointer_path = FENG_REPO_ROOT / UNIVERSE_POINTER_PATH
    if not pointer_path.is_symlink():
        raise FreezeGuardError("detectable-universe current pointer is not a symlink")
    pointer_target_text = os.readlink(pointer_path)
    tagged_pointer_target = _git(
        ["show", f"{tag_commit}:{UNIVERSE_POINTER_PATH.as_posix()}"], text=False
    ).stdout
    if tagged_pointer_target != pointer_target_text.encode("utf-8"):
        raise FreezeGuardError("detectable-universe pointer differs from frozen tag")

    if Path(pointer_target_text).is_absolute():
        raise FreezeGuardError("detectable-universe pointer target is absolute")
    resolved_pointer_text = posixpath.normpath(
        (UNIVERSE_POINTER_PATH.parent / pointer_target_text).as_posix()
    )
    resolved_pointer = Path(resolved_pointer_text)
    try:
        resolved_pointer.relative_to(UNIVERSE_BUNDLE_ROOT_PATH)
    except ValueError as exc:
        raise FreezeGuardError(
            "detectable-universe pointer target escapes its bundle root"
        ) from exc

    resolved_artifact = resolved_pointer / relative.name
    current_path = FENG_REPO_ROOT / resolved_artifact
    bundle_root = FENG_REPO_ROOT / UNIVERSE_BUNDLE_ROOT_PATH
    if (
        bundle_root.is_symlink()
        or current_path.is_symlink()
        or not current_path.resolve().is_relative_to(bundle_root.resolve())
        or logical_path.resolve() != current_path.resolve()
    ):
        raise FreezeGuardError(f"invalid versioned-bundle path for {relative}")
    return (
        current_path.read_bytes(),
        _git(
            ["show", f"{tag_commit}:{resolved_artifact.as_posix()}"], text=False
        ).stdout,
    )


def _verified_freeze_manifest(path: Path | None) -> dict:
    """Return the exact tagged Phase-9 manifest, or fail closed."""
    if path is None or Path(path).resolve() != FREEZE_MANIFEST_PATH.resolve():
        raise FreezeGuardError("freeze manifest is not at its fixed repository path")
    if not Path(path).is_file():
        raise FreezeGuardError("freeze manifest is missing")
    try:
        record = json.loads(Path(path).read_text(encoding="utf-8"))
        if set(record) != {"schema", "tag", "preregistration", "artifacts"}:
            raise FreezeGuardError("freeze manifest has an unexpected schema")
        if record["schema"] != FREEZE_MANIFEST_SCHEMA:
            raise FreezeGuardError("freeze manifest schema version is not accepted")
        tag = record["tag"]
        if tag != APPROVED_FREEZE_TAG:
            raise FreezeGuardError("manifest does not name the exact approved freeze tag")
        if APPROVED_FREEZE_TAG_OBJECT_SHA256 == "0" * 64:
            raise FreezeGuardError("approved external freeze tag digest is not pinned")
        repo = Path(_git(["rev-parse", "--show-toplevel"]).stdout.strip()).resolve()
        if repo != FENG_REPO_ROOT.resolve():
            raise FreezeGuardError("git command resolved a different repository")
        tag_commit = _git(["rev-parse", f"{tag}^{{commit}}"]).stdout.strip()
        if re.fullmatch(r"[0-9a-f]{40}", tag_commit) is None:
            raise FreezeGuardError("freeze tag did not resolve to a commit")
        tag_object = _git(["cat-file", "tag", f"refs/tags/{tag}"], text=False).stdout
        if hashlib.sha256(tag_object).hexdigest() != APPROVED_FREEZE_TAG_OBJECT_SHA256:
            raise FreezeGuardError("freeze tag object does not match the external trust root")

        manifest_rel = FREEZE_MANIFEST_PATH.resolve().relative_to(FENG_REPO_ROOT.resolve())
        tagged_manifest = _git(
            ["show", f"{tag_commit}:{manifest_rel.as_posix()}"], text=False
        ).stdout
        if tagged_manifest != Path(path).read_bytes():
            raise FreezeGuardError("manifest bytes differ from the tagged commit")

        prereg = record["preregistration"]
        if set(prereg) != {"path", "sha256", "locked"}:
            raise FreezeGuardError("preregistration manifest entry has unexpected fields")
        if prereg["path"] != EXPECTED_PREREGISTRATION_PATH or prereg["locked"] is not True:
            raise FreezeGuardError("manifest does not bind the expected locked preregistration")
        artifacts = record["artifacts"]
        if not isinstance(artifacts, list) or any(
            not isinstance(item, dict) or set(item) != {"path", "sha256"}
            for item in artifacts
        ):
            raise FreezeGuardError("artifact manifest entries have unexpected fields")
        artifact_paths = [item["path"] for item in artifacts]
        if len(artifact_paths) != len(set(artifact_paths)) or set(artifact_paths) != EXPECTED_ARTIFACT_PATHS:
            raise FreezeGuardError("manifest does not contain the exact Phase-9 artifact set")

        bound = [prereg, *artifacts]
        for item in bound:
            relative = Path(item["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise FreezeGuardError("manifest contains an unsafe path")
            current, tagged = _tagged_and_worktree_bytes(tag_commit, relative)
            digest = hashlib.sha256(current).hexdigest()
            if current != tagged or digest != item["sha256"]:
                raise FreezeGuardError(f"tag/hash mismatch for {relative}")
        record["tag_commit"] = tag_commit
        return record
    except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError) as exc:
        raise FreezeGuardError(f"invalid freeze manifest: {exc}") from exc


def verify_freeze_manifest(path: Path | None) -> bool:
    """Return whether the fixed manifest and every bound file verify."""
    try:
        _verified_freeze_manifest(path)
        return True
    except FreezeGuardError:
        return False


def verify_scored_strata(path: Path, freeze_manifest: Path | None) -> bool:
    """Bind the exact scored dataframe bytes to the verified frozen manifest."""
    try:
        record = _verified_freeze_manifest(freeze_manifest)
        expected = (FENG_REPO_ROOT / EXPECTED_SCORED_STRATA_PATH).resolve()
        if Path(path).resolve() != expected:
            return False
        item = next(
            entry for entry in record["artifacts"]
            if entry["path"] == EXPECTED_SCORED_STRATA_PATH
        )
        return hashlib.sha256(Path(path).read_bytes()).hexdigest() == item["sha256"]
    except (FreezeGuardError, OSError, StopIteration):
        return False


def phase9_is_unlocked(
    sentinel: Path | None = None, *, freeze_manifest: Path | None = None
) -> bool:
    """True iff explicit unlock and a commit-bound freeze manifest are present."""
    explicit_unlock = os.environ.get(PHASE9_UNLOCK_ENV) == "1" or (
        sentinel is not None and Path(sentinel).is_file()
    )
    return explicit_unlock and verify_freeze_manifest(freeze_manifest)


def assert_phase9_unlocked(
    sentinel: Path | None = None, *, freeze_manifest: Path | None = None
) -> None:
    """Refuse to run any outcome computation before the Phase-9 unlock."""
    if not phase9_is_unlocked(sentinel, freeze_manifest=freeze_manifest):
        raise FreezeGuardError(
            "Arm D permutation test is sealed until Phase 9. A verified freeze "
            "manifest must bind the intended feng-freeze* commit and hashed "
            "preregistration/artifacts, and an explicit unlock must be supplied. Set "
            f"{PHASE9_UNLOCK_ENV}=1 (or pass an existing sentinel file) only "
            "after the preregistration freeze tag is cut. Refusing to compute "
            "delta_r or a permutation p-value pre-freeze."
        )


def percentile_within_stratum(
    positive_value: float, background_values: np.ndarray, *, lower_is_better: bool
) -> float:
    """Percentile rank of a positive against its matched background.

    Returns the fraction of the (positive + background) set that the positive
    ranks at or better than, in [0, 1], where "better" means lower when
    ``lower_is_better`` (Hamming distance, IntaRNA energy) and higher otherwise.
    A pure ranking primitive -- safe to test on synthetic data. It does not, by
    itself, read any real score.
    """
    background = np.asarray(background_values, dtype=float)
    if background.size == 0:
        raise ValueError("cannot rank a positive against an empty stratum")
    total = background.size + 1
    if lower_is_better:
        rank_better = int(np.sum(background > positive_value))
    else:
        rank_better = int(np.sum(background < positive_value))
    ties = int(np.sum(background == positive_value))
    # Mid-rank convention for ties, counting the positive itself.
    return (rank_better + (ties + 1) / 2) / total


def run_phase9_test(
    scored_strata: Path,
    *,
    baseline_col: str = "baseline_min_hamming",
    test_col: str = "intarna_min_energy",
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    seed: int = DEFAULT_SEED,
    sentinel: Path | None = None,
    freeze_manifest: Path | None = None,
) -> dict:  # pragma: no cover - Phase 9 only, never run pre-freeze
    """Run the one-shot matched-stratum permutation test.

    SEALED until Phase 9. Calls :func:`assert_phase9_unlocked` before touching
    any score, then would compute delta_r and its one-sided permutation p-value
    over ``n_permutations`` shuffles with ``seed``. The body is intentionally
    left unimplemented so it cannot run pre-freeze even if the guard were
    bypassed.
    """
    assert_phase9_unlocked(sentinel, freeze_manifest=freeze_manifest)
    if not verify_scored_strata(Path(scored_strata), freeze_manifest):
        raise FreezeGuardError(
            "scored_strata is not the exact tagged dataframe bound by the freeze manifest"
        )
    raise NotImplementedError(
        "Phase 9 execution: implement against the frozen preregistered score "
        "columns after the freeze tag is cut."
    )
