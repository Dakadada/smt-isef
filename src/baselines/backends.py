"""Fail-closed contracts for frozen cis-splicing and accessibility backends."""
from abc import ABC, abstractmethod
from dataclasses import dataclass
import hashlib
import subprocess


@dataclass(frozen=True)
class FrozenBackendContract:
    name: str
    package_version: str
    weight_sha256: str
    context_bases: int
    padding: str
    coordinate_convention: str
    strand_handling: str
    reference_convention: str
    channels: tuple
    golden_outputs: tuple = ()

    def validate(self):
        if self.name not in ("spliceai", "pangolin"):
            raise ValueError("only frozen SpliceAI/Pangolin contracts are scientific backends")
        if not self.package_version or len(self.weight_sha256) != 64:
            raise ValueError("backend requires version and 64-hex weight checksum")
        try: int(self.weight_sha256, 16)
        except ValueError as exc: raise ValueError("invalid weight SHA-256") from exc
        if self.context_bases <= 0 or self.padding not in ("N", "reference"):
            raise ValueError("context and padding must be explicit")
        if self.coordinate_convention != "0-based-sequence-position" or self.strand_handling not in ("transcript-5to3", "genomic-with-explicit-strand"):
            raise ValueError("unsupported coordinate/strand convention")
        if set(self.channels) != {"donor", "acceptor"} or not self.reference_convention:
            raise ValueError("typed donor/acceptor channels and reference convention required")


class CisSplicingBackend(ABC):
    name, scientific = "abstract", False

    @abstractmethod
    def score_channels(self, sequence):
        """Return donor and acceptor vectors on transcript 5'->3'."""


class StubCisBackend(CisSplicingBackend):
    """Typed deterministic test backend; explicitly not biological."""
    def __init__(self, name="stub-cis"): self.name = name

    def score_channels(self, sequence):
        seq = sequence.upper().replace("U", "T"); donor, acceptor = [], []
        for i, base in enumerate(seq):
            digest = hashlib.sha256((str(i) + seq[max(0, i - 1):i + 2]).encode()).digest()
            donor.append((digest[0] + {"A": 0, "C": 13, "G": 29, "T": 47}.get(base, 0)) / 302.0)
            acceptor.append((digest[1] + {"A": 41, "C": 23, "G": 7, "T": 0}.get(base, 0)) / 296.0)
        return {"donor": donor, "acceptor": acceptor}


class VerifiedCisBackend(CisSplicingBackend):
    """Adapter accepted only after manifest, reference and official golden checks."""
    scientific = True
    def __init__(self, contract, scorer, *, loaded_weight_sha256, reference_convention, golden_observations):
        contract.validate()
        if loaded_weight_sha256 != contract.weight_sha256: raise ValueError("loaded weights do not match frozen checksum")
        if reference_convention != contract.reference_convention: raise ValueError("reference convention mismatch")
        observed = tuple((case, tuple(float(x) for x in values)) for case, values in golden_observations)
        if not observed or observed != contract.golden_outputs: raise ValueError("official golden-output verification missing or mismatched")
        self.name, self.contract, self._scorer = contract.name, contract, scorer

    def score_channels(self, sequence):
        seq = sequence.upper().replace("U", "T")
        raw = self._scorer(seq, context=self.contract.context_bases, padding=self.contract.padding, strand=self.contract.strand_handling)
        if not isinstance(raw, dict) or set(raw) != {"donor", "acceptor"}: raise ValueError("backend must return donor and acceptor channels")
        result = {key: list(map(float, raw[key])) for key in raw}
        if any(len(v) != len(seq) for v in result.values()): raise ValueError("backend channel length does not match sequence")
        if any(not 0 <= x <= 1 for values in result.values() for x in values): raise ValueError("backend propensity outside [0,1]")
        return result


class CallableCisBackend(VerifiedCisBackend):
    """Compatibility name retaining the strict VerifiedCisBackend signature."""


class AccessibilityBackend(ABC):
    @abstractmethod
    def accessibility(self, sequence, start, end): pass


class StubAccessibilityBackend(AccessibilityBackend):
    def accessibility(self, sequence, start, end):
        target = sequence[start:end].upper()
        return 0.0 if not target else sum(b in "AUAT" for b in target) / len(target)


class RNAplfoldBackend(AccessibilityBackend):
    def __init__(self, executable="RNAplfold", window=80, span=40): self.executable, self.window, self.span = executable, int(window), int(span)
    def accessibility(self, sequence, start, end):
        subprocess.run([self.executable, "-W", str(self.window), "-L", str(self.span), "-u", str(end-start)], input=sequence+"\n", text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        try:
            with open("plfold_lunp") as handle: rows = [x.split() for x in handle if x and not x.startswith("#")]
            value = float(next(r for r in rows if int(r[0]) == end)[end-start])
        except (OSError, StopIteration, IndexError, ValueError) as exc: raise RuntimeError("could not parse RNAplfold output") from exc
        if not 0 <= value <= 1: raise RuntimeError("RNAplfold accessibility outside [0,1]")
        return value
