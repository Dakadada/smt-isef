"""Hard project constraint: no catalysis, no degradation, anywhere in the model.

This project models a NON-CATALYTIC steric blockade only. dCas13 is
catalytically dead (the "d" is deliberate: its HEPN nuclease is inactivated by
point mutation) and this codebase never gives it a k_cat. Likewise no
RNase-H-style degradation pathway is modelled for the ASO comparison case --
an ASO here is a passive steric occupant, not a recruiter of RNase H. Both
mechanisms are out of scope by design, not by omission, so their vocabulary
(k_cat, kcat, cleavage, RNase, degradation rate, turnover) must never appear
as a live identifier or an undisclaimed data value in the executable source
of ``smt.steric`` or ``smt.thermo``.

How the check works
--------------------
We parse (never just grep) every ``.py`` file under ``src/smt/steric`` and
``src/smt/thermo`` with :mod:`ast` and apply two different rules:

1. IDENTIFIERS (variable/argument/function/class names, attribute access,
   keyword-argument names, import names) matching a forbidden root FAIL
   unconditionally. An identifier is never documentation -- ``k_cat = 0``
   is still modelling catalytic turnover as a concept, just setting it to
   zero, and that is exactly the thing this project must never do. Because
   comments are stripped by the tokenizer before ``ast`` ever sees them,
   this rule is structurally incapable of matching a comment.

2. STRING LITERALS (this also covers docstrings, since a docstring is just a
   string literal in statement position) MAY mention a forbidden term only
   if the same string also carries an explicit disclaiming cue ("no",
   "not", "non-", "never", "out of scope", "cannot", "zero", "excluded",
   "disclaim", ...). This is what lets
   ``"No turnover chemistry ... is modelled"`` (the real text in
   ``shadow.py``) pass, while a bare, undisclaimed
   ``"apply degradation rate"`` would fail.

Rule (1) is what stops trivial bypass: hiding the concept behind a
zero-valued variable, a differently-cased identifier, or a dynamically
dispatched attribute name still trips it, because it is a structural
property of the AST, not a textual pattern that string formatting or
whitespace could dodge. Rule (2) is what stops trivial over-breaking: it
would otherwise be impossible to even document, in the modules' own
docstrings, that these mechanisms are intentionally absent.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCANNED_DIRS = [REPO_ROOT / "src" / "smt" / "steric", REPO_ROOT / "src" / "smt" / "thermo"]

# Forbidden roots, case-insensitive. Substring match is deliberate: it also
# catches "kcat_placeholder", "cleavage_rate", "RNaseH_site", etc.
_FORBIDDEN_PATTERNS: dict[str, re.Pattern[str]] = {
    "k_cat / kcat": re.compile(r"k_?cat", re.IGNORECASE),
    "cleavage": re.compile(r"cleav", re.IGNORECASE),
    "RNase": re.compile(r"rnase", re.IGNORECASE),
    "degradation": re.compile(r"degrad", re.IGNORECASE),
    "turnover": re.compile(r"turnover", re.IGNORECASE),
}

# A string literal containing a forbidden term is allowed only if it also
# contains one of these disclaiming cues. Checked case-insensitively over the
# whole string (a docstring/message is short enough that "same string" is a
# meaningful, easily-explained scope -- no need for a token-distance window).
_DISCLAIM_CUE = re.compile(
    r"\b(no|not|non-|never|zero|cannot|can't|doesn't|does not|"
    r"out of scope|excluded|disclaim|deliberately)\b",
    re.IGNORECASE,
)


def _forbidden_hits(text: str) -> list[str]:
    return [label for label, pat in _FORBIDDEN_PATTERNS.items() if pat.search(text)]


def _iter_source_files() -> list[Path]:
    files: list[Path] = []
    for d in SCANNED_DIRS:
        assert d.is_dir(), f"expected directory {d} to exist"
        files.extend(sorted(p for p in d.rglob("*.py") if "__pycache__" not in p.parts))
    assert files, f"no .py files found under {SCANNED_DIRS}; check the scan paths"
    return files


class _Violation(str):
    """A human-readable violation report line; subclassing str for cheap repr."""


def _check_file(path: Path) -> list[str]:
    source = path.read_text()
    tree = ast.parse(source, filename=str(path))
    rel = path.relative_to(REPO_ROOT)
    violations: list[str] = []

    # --- Rule 1: identifiers -------------------------------------------------
    identifier_nodes: list[tuple[int, str, str]] = []  # (lineno, kind, name)
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            identifier_nodes.append((node.lineno, "name", node.id))
        elif isinstance(node, ast.Attribute):
            identifier_nodes.append((node.lineno, "attribute", node.attr))
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            identifier_nodes.append((node.lineno, "function name", node.name))
            for arg in (
                list(node.args.posonlyargs)
                + list(node.args.args)
                + list(node.args.kwonlyargs)
                + ([node.args.vararg] if node.args.vararg else [])
                + ([node.args.kwarg] if node.args.kwarg else [])
            ):
                identifier_nodes.append((arg.lineno, "argument", arg.arg))
        elif isinstance(node, ast.ClassDef):
            identifier_nodes.append((node.lineno, "class name", node.name))
        elif isinstance(node, ast.keyword) and node.arg is not None:
            identifier_nodes.append((node.lineno, "keyword argument", node.arg))
        elif isinstance(node, ast.alias):
            for nm in (node.name, node.asname):
                if nm:
                    identifier_nodes.append((getattr(node, "lineno", 0), "import name", nm))

    for lineno, kind, name in identifier_nodes:
        hits = _forbidden_hits(name)
        if hits:
            violations.append(
                f"{rel}:{lineno}: {kind} {name!r} implies catalysis/degradation "
                f"({', '.join(hits)}). dCas13 is catalytically dead and this project "
                "never models turnover or degradation kinetics -- remove the "
                "identifier rather than zeroing it out."
            )

    # --- Rule 2: string literals (docstrings included) -----------------------
    docstring_ids = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            doc_node = ast.get_docstring(node, clean=False)
            if doc_node is not None and node.body:
                first = node.body[0]
                if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
                    docstring_ids.add(id(first.value))

    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            hits = _forbidden_hits(node.value)
            if not hits:
                continue
            if _DISCLAIM_CUE.search(node.value):
                continue  # explicitly disclaiming -- allowed, docstring or not
            kind = "docstring" if id(node) in docstring_ids else "string literal"
            snippet = node.value.strip().replace("\n", " ")[:80]
            violations.append(
                f"{rel}:{node.lineno}: undisclaimed {kind} mentions catalysis/degradation "
                f"({', '.join(hits)}) with no disclaiming cue: {snippet!r}. "
                "Either state explicitly that the mechanism is absent/out of scope, "
                "or remove the reference."
            )

    return violations


@pytest.mark.parametrize("path", _iter_source_files(), ids=lambda p: str(p.relative_to(REPO_ROOT)))
def test_no_catalysis_vocabulary_in_source(path: Path) -> None:
    violations = _check_file(path)
    assert not violations, "\n" + "\n".join(violations)


def test_scanner_actually_flags_identifiers() -> None:
    """Meta-test: the identifier rule must fire on a synthetic k_cat use.

    Guards against the checker being silently broken (e.g. by an import
    error swallowed somewhere) so this file can't pass by accidentally
    checking nothing.
    """
    src = "def rate(k_cat, substrate):\n    return k_cat * substrate\n"
    tree = ast.parse(src)
    names = [n.id for n in ast.walk(tree) if isinstance(n, ast.Name)]
    args = [
        a.arg
        for n in ast.walk(tree)
        if isinstance(n, ast.FunctionDef)
        for a in n.args.args
    ]
    assert any(_forbidden_hits(n) for n in names + args)


def test_scanner_allows_explicit_disclaimer() -> None:
    """Meta-test: a disclaiming docstring must NOT be flagged."""
    text = "No turnover chemistry and no RNase-H-style degradation is modelled here."
    assert _forbidden_hits(text)  # the words are present...
    assert _DISCLAIM_CUE.search(text)  # ...but so is the disclaiming cue

    src = f'"""{text}"""\ndef f():\n    pass\n'
    tmp_tree = ast.parse(src)
    # exercise the same Constant-node path _check_file uses
    for node in ast.walk(tmp_tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            assert _forbidden_hits(node.value)
            assert _DISCLAIM_CUE.search(node.value)


def test_scanner_rejects_undisclaimed_string() -> None:
    """Meta-test: a bare, undisclaimed mention must be flagged."""
    text = "apply degradation rate to the bound complex"
    assert _forbidden_hits(text)
    assert not _DISCLAIM_CUE.search(text)
