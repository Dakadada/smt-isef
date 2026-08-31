"""Exhaustive, recall-oriented full-pre-mRNA candidate enumeration."""

import math
from dataclasses import replace
from typing import Iterable, List, Sequence, Tuple

from .alignment import global_edit_alignment, normalize, reverse_complement
from .chemistry import estimate_delta_g, parameters_for
from .models import Candidate, EnumeratorConfig, Evidence, Transcript


def _overlapping_evidence(evidence: Iterable[Evidence], transcript_id: str,
                          start: int, end: int) -> Tuple[Evidence, ...]:
    return tuple(sorted((item for item in evidence
                         if item.transcript_id == transcript_id
                         and item.start < end and item.end > start),
                        key=lambda x: (x.start, x.end, x.source, x.rbp)))


def _effect_sign(items: Sequence[Evidence]) -> str:
    """A sign is impossible without explicit signed RBP/perturbation evidence."""
    supported = {item.effect_sign for item in items
                 if item.effect_sign in {"enhancer", "silencer"}
                 and (item.source.lower() == "eclip" or item.perturbation_supported)}
    return supported.pop() if len(supported) == 1 else "unknown"


def _expression_localization_weight(tpm: float, nuclear: float,
                                    config: EnumeratorConfig) -> float:
    if tpm < 0 or not 0 <= nuclear <= 1:
        raise ValueError("expression must be nonnegative and nuclear_fraction in [0,1]")
    expression = tpm / (tpm + config.expression_half_saturation_tpm) if tpm else 0.0
    localization = config.localization_floor + (1 - config.localization_floor) * nuclear
    return expression * localization


def enumerate_offtargets(effector_sequence: str, transcripts: Iterable[Transcript],
                         config: EnumeratorConfig = EnumeratorConfig(),
                         evidence: Iterable[Evidence] = ()) -> List[Candidate]:
    """Scan every offset of every eligible complete pre-mRNA; return ranked hits."""
    query = reverse_complement(effector_sequence)
    params = parameters_for(config.chemistry, config.allow_unknown_chemistry_heuristic)
    evidence = tuple(evidence)
    candidates: List[Candidate] = []
    for transcript in sorted(transcripts, key=lambda x: x.transcript_id):
        sequence = normalize(transcript.sequence)
        weight = _expression_localization_weight(
            transcript.expression_tpm, transcript.nuclear_fraction, config)
        minimum = (config.min_alignment_length if config.min_alignment_length is not None
                   else max(1, len(query) - config.max_indel_bases))
        maximum = (config.max_alignment_length if config.max_alignment_length is not None
                   else len(query) + config.max_indel_bases)
        for start in range(len(sequence)):
            for length in range(minimum, maximum + 1):
                end = start + length
                if end > len(sequence):
                    continue
                alignment = global_edit_alignment(query, sequence[start:end])
                indels = alignment.insertions + alignment.deletions
                if alignment.distance > config.max_edits or indels > config.max_indel_bases:
                    continue
                items = _overlapping_evidence(evidence, transcript.transcript_id, start, end)
                dg = (estimate_delta_g(len(query), alignment.mismatches, indels, params)
                      if params is not None else None)
                edit_affinity = math.exp(-alignment.distance / max(1.0, len(query) * 0.22))
                energy_affinity = (1.0 / (1.0 + math.exp((dg + len(query)) / 3.0))
                                   if dg is not None else 0.0)
                rbp_support = min(1.0, sum(max(0.0, x.occupancy) for x in items))
                # Scores rank/annotate the exhaustive bounded universe; they never gate it.
                score = (0.55 * edit_affinity + 0.25 * energy_affinity
                         + 0.20 * weight + config.rbp_weight * rbp_support)
                if params is None or not params.known:
                    score *= config.unknown_chemistry_penalty
                score = min(1.0, score)
                confidence = min(0.99, max(0.01,
                    0.15 + 0.55 * edit_affinity + 0.15 * weight + 0.15 * rbp_support))
                candidates.append(Candidate(
                    0, transcript.transcript_id, transcript.gene_id, start, end,
                    sequence[start:end], alignment.mismatches, alignment.insertions,
                    alignment.deletions, alignment.distance, alignment.cigar,
                    config.chemistry if params is None else params.name,
                    bool(params and params.known), dg,
                    ("none" if params is None else
                     "heuristic-non-calibrated" if params.heuristic else "chemistry-specific"),
                    transcript.expression_tpm, transcript.nuclear_fraction, weight,
                    items, _effect_sign(items), score, confidence,
                    transcript.expression_tpm >= config.min_expression_tpm,
                    score >= config.min_score,
                ))
    candidates.sort(key=lambda x: (-x.recall_score, x.edit_distance,
                                   x.transcript_id, x.start, x.end))
    return [replace(candidate, rank=index) for index, candidate in enumerate(candidates, 1)]
