"""Verify load-bearing claim boundaries for the clean-room Feng arm."""

import json

import pytest

from feng.honesty.claims import Claim
from feng.honesty.linter import (
    FORBIDDEN_FIELD_NAMES,
    FORBIDDEN_STRINGS,
    ClaimLintError,
    lint_claim,
)
from feng.honesty.states import ALLOWED_RESULT_STATES, ResultState


def test_allowed_states_are_exactly_frozen_vocabulary() -> None:
    assert ALLOWED_RESULT_STATES == {
        "candidate",
        "paper_supported",
        "clean_reanalysis_reproduced",
    }
    assert {state.value for state in ResultState} == ALLOWED_RESULT_STATES


@pytest.mark.parametrize("phrase", FORBIDDEN_STRINGS)
def test_forbidden_claim_strings_are_blocked(phrase: str) -> None:
    with pytest.raises(ClaimLintError):
        lint_claim(f"This output is {phrase} across the evaluation.")


@pytest.mark.parametrize("field_name", sorted(FORBIDDEN_FIELD_NAMES))
def test_forbidden_field_names_are_blocked(field_name: str) -> None:
    with pytest.raises(ClaimLintError):
        lint_claim("A bounded candidate result.", [field_name])


def test_phrase_matching_is_case_and_whitespace_insensitive() -> None:
    with pytest.raises(ClaimLintError):
        lint_claim("A GENOME-WIDE\n safety map is shown.")


def test_override_requires_reason_and_log(tmp_path) -> None:
    with pytest.raises(ClaimLintError):
        lint_claim("This is generalizable.", override_reason="reviewed")
    with pytest.raises(ClaimLintError):
        lint_claim("This is generalizable.", override_log=tmp_path / "override.jsonl")


def test_explicit_override_is_logged(tmp_path) -> None:
    log_path = tmp_path / "audit" / "claim-overrides.jsonl"
    violations = lint_claim(
        "This is the first bounded demonstration.",
        ["safe"],
        override_reason="Quoted historical wording required by the figure caption.",
        override_log=log_path,
    )

    record = json.loads(log_path.read_text(encoding="utf-8"))
    assert {(item.kind, item.value) for item in violations} == {("string", "first"), ("field", "safe")}
    assert record["reason"].startswith("Quoted historical wording")
    assert len(record["claim_sha256"]) == 64
    assert record["violations"] == [
        {"kind": "string", "value": "first"},
        {"kind": "field", "value": "safe"},
    ]


def test_claim_constructor_enforces_state_and_linter() -> None:
    claim, violations = Claim.create("A candidate under this specification.", "candidate")
    assert claim.state is ResultState.CANDIDATE
    assert violations == ()
    with pytest.raises(ValueError):
        Claim.create("A bounded candidate.", "confirmed")
    with pytest.raises(ClaimLintError):
        Claim.create("A clinical response claim.", ResultState.CANDIDATE)
    with pytest.raises(ClaimLintError):
        Claim("This is generalizable.", ResultState.CANDIDATE)
