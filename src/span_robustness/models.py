"""Honesty-gated types for the closed-world span method fixture."""
from __future__ import annotations
from dataclasses import dataclass,fields,is_dataclass
from enum import Enum
from typing import Any
from design_engine import PROSPECTIVE_BANNER,PROSPECTIVE_STATUS
METHOD_LABEL="method/engine-wiring fixture only"
CLAIM_BOUNDARY=("Metrics describe only the frozen non-locus method fixture and its explicit event geometry; "
 "NOT STMN2/UNC13A or any locus inference, NOT the P4 data-conditional benchmark, NOT experimentally validated.")
def _plain(v:Any)->Any:
    if is_dataclass(v): return {f.name:_plain(getattr(v,f.name)) for f in fields(v)}
    if isinstance(v,Enum): return v.value
    if isinstance(v,(tuple,list)): return [_plain(x) for x in v]
    return v
class PredictedDirection(str,Enum):
    INCREASES_INCLUSION="increases_inclusion"; DECREASES_INCLUSION="decreases_inclusion"; ABSTAIN="abstain"
@dataclass(frozen=True)
class PerturbationSpec:
    window_size:int; stride:int; coordinate_convention:str="zero-based half-open"; operation:str="slide fixed-width ASO footprint"
    def __post_init__(self):
        if self.window_size<1 or self.stride<1 or self.stride>self.window_size: raise ValueError("stride must be <= window_size")
        if self.coordinate_convention!="zero-based half-open" or self.operation!="slide fixed-width ASO footprint": raise ValueError("unsupported perturbation")
@dataclass(frozen=True)
class Lead:
    locus:str; role:str="lead-only-not-evaluated"; limitation:str="No committed cryptic-exon context or qualifying matched reads."
    def __post_init__(self):
        if self.locus not in {"STMN2","UNC13A"} or self.role!="lead-only-not-evaluated": raise ValueError("invalid lead")
@dataclass(frozen=True)
class WindowPrediction:
    index:int; start:int; end:int; direction:PredictedDirection
    def __post_init__(self):
        if not isinstance(self.direction,PredictedDirection) or self.index<0 or self.start<0 or self.end<=self.start: raise ValueError("invalid window")
@dataclass(frozen=True)
class FlipPair:
    left_index:int; right_index:int; left_direction:PredictedDirection; right_direction:PredictedDirection
    def __post_init__(self):
        if self.right_index!=self.left_index+1 or self.left_direction==self.right_direction or PredictedDirection.ABSTAIN in (self.left_direction,self.right_direction): raise ValueError("invalid flip")
@dataclass(frozen=True)
class SpanRobustnessReport:
    schema_version:str; banner:str; prediction_status:str; method_label:str; claim_boundary:str
    fixture_id:str; sequence:str; sequence_sha256:str; accessible_elements:tuple[tuple[int,int],...]
    silencers:tuple[tuple[int,int],...]; perturbation:PerturbationSpec; engine_provenance:tuple[tuple[str,str],...]
    windows:tuple[WindowPrediction,...]; evaluable_window_count:int; abstention_count:int
    modal_direction:PredictedDirection|None; modal_direction_fraction:float|None
    adjacent_evaluable_pair_count:int; flip_pairs:tuple[FlipPair,...]; adjacent_flip_fraction:float|None
    leads:tuple[Lead,...]
    def __post_init__(self): self.validate_for_rendering()
    def validate_for_rendering(self):
        from . import method
        if self.schema_version!="prospective-span-robustness-method/3": raise ValueError("unsupported schema")
        if self.banner!=PROSPECTIVE_BANNER or self.prediction_status!=PROSPECTIVE_STATUS: raise ValueError("prospective status required")
        if self.method_label!=METHOD_LABEL or self.claim_boundary!=CLAIM_BOUNDARY: raise ValueError("claim boundary changed")
        if (self.fixture_id,self.sequence,self.sequence_sha256,self.accessible_elements,self.silencers)!=(method.FIXTURE_ID,method.FIXTURE_SEQUENCE,method.FIXTURE_SEQUENCE_SHA256,method.FIXTURE_ACCESSIBLE,method.FIXTURE_SILENCERS): raise ValueError("fixture identity/context is not allowlisted")
        if self.engine_provenance!=method._approved_provenance(): raise ValueError("engine provenance/code/config mismatch")
        if self.perturbation.window_size>len(self.sequence): raise ValueError("window exceeds fixture")
        expected=method._windows(self.perturbation)
        if self.windows!=expected: raise ValueError("directions are not deterministic approved-engine output")
        ev,modal,fraction,pairs,flips=method._summaries(expected)
        actual=(self.evaluable_window_count,self.abstention_count,self.modal_direction,self.modal_direction_fraction,self.adjacent_evaluable_pair_count,self.flip_pairs,self.adjacent_flip_fraction)
        derived=(len(ev),len(expected)-len(ev),modal,fraction,len(pairs),flips,len(flips)/len(pairs) if pairs else None)
        if actual!=derived: raise ValueError("derived span summary mismatch")
        if self.leads!=(Lead("STMN2"),Lead("UNC13A")): raise ValueError("lead declaration changed")
    def to_dict(self): self.validate_for_rendering(); return _plain(self)
