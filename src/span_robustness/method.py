"""Closed-world method fixture for sliding-window engine-direction robustness."""
from __future__ import annotations
from dataclasses import asdict
import hashlib, inspect, json
from pathlib import Path
from design_engine import PROSPECTIVE_BANNER, PROSPECTIVE_STATUS
from design_engine.backends import DeterministicStubFolder
from offtarget.alignment import reverse_complement
from physics_adapter import (CoordinateConvention, DuplexOrientation, EffectorPhysicsAdapter,
    EffectorSpec, EventDefinition, IntermolecularComplexSpec, IntermolecularPair, Interval,
    ModelConfig, TranscriptContext, protocol_hash)
from .models import (CLAIM_BOUNDARY, METHOD_LABEL, FlipPair, Lead, PerturbationSpec,
    PredictedDirection, SpanRobustnessReport, WindowPrediction)

FIXTURE_ID = "METHOD_FIXTURE_NON_LOCUS_001"
FIXTURE_SEQUENCE = "ACGTTGCATGCAACGT"
FIXTURE_SEQUENCE_SHA256 = hashlib.sha256(FIXTURE_SEQUENCE.encode("ascii")).hexdigest()
FIXTURE_ACCESSIBLE = ((2, 5),)
FIXTURE_SILENCERS = ((9, 12),)
APPROVED_ENGINE_IDENTITY = "physics_adapter.EffectorPhysicsAdapter+design_engine.DeterministicStubFolder/1"

def _sha(value: str) -> str: return hashlib.sha256(value.encode()).hexdigest()
def _source_hash(obj) -> str: return _sha(inspect.getsource(obj))
def _canonical(value) -> str: return json.dumps(value,sort_keys=True,separators=(",",":"))

def _approved_provenance() -> tuple[tuple[str,str],...]:
    config=ModelConfig()
    package_sources={}
    package_dirs={Path(inspect.getsourcefile(obj)).resolve().parent
                  for obj in (EffectorPhysicsAdapter,DeterministicStubFolder)}
    for directory in sorted(package_dirs):
        for path in sorted(directory.glob("*.py")):
            package_sources[f"{directory.name}/{path.name}"]=_sha(path.read_text(encoding="utf-8"))
    return tuple(sorted({
        "engine_identity":APPROVED_ENGINE_IDENTITY,
        "backend_version":DeterministicStubFolder.backend_version,
        "model_config_json":_canonical(asdict(config)),
        "model_protocol_sha256":protocol_hash(config),
        "adapter_source_sha256":_source_hash(EffectorPhysicsAdapter),
        "backend_source_sha256":_source_hash(DeterministicStubFolder),
        "span_predictor_source_sha256":_APPROVED_PREDICTOR_SOURCE_SHA256,
        "package_source_manifest_json":_canonical(package_sources),
        "package_source_manifest_sha256":_sha(_canonical(package_sources)),
        "outcomes_consumed":"false",
        "fixture_scope":"method-engine-wiring-only",
    }.items()))

def _event() -> EventDefinition:
    return EventDefinition("METHOD_FIXTURE_EVENT",
        tuple(Interval(*x) for x in FIXTURE_ACCESSIBLE), (),
        tuple(Interval(*x) for x in FIXTURE_SILENCERS))

def _predict_direction(start:int,end:int) -> PredictedDirection:
    context=TranscriptContext(FIXTURE_ID,FIXTURE_SEQUENCE); target=Interval(start,end)
    effector=reverse_complement(FIXTURE_SEQUENCE[start:end])
    pairs=tuple(IntermolecularPair(start+i,len(effector)-1-i) for i in range(len(effector)))
    complex_spec=IntermolecularComplexSpec(target,CoordinateConvention.ZERO_BASED_HALF_OPEN,
        effector,FIXTURE_ID,f"fixture:{start}-{end}",DuplexOrientation.ANTIPARALLEL,start,"PMO",pairs)
    spec=EffectorSpec(f"fixture:{start}-{end}",target,complex_spec,False,())
    prediction=EffectorPhysicsAdapter(DeterministicStubFolder(),ModelConfig()).predict(context,spec,_event())
    return {1:PredictedDirection.INCREASES_INCLUSION,-1:PredictedDirection.DECREASES_INCLUSION,
            None:PredictedDirection.ABSTAIN}[prediction.primary_direction]

_APPROVED_PREDICTOR_SOURCE_SHA256 = _source_hash(_predict_direction)

def _windows(spec:PerturbationSpec, _approved_predictor=_predict_direction) -> tuple[WindowPrediction,...]:
    """Use the definition-time approved predictor even if the module name is rebound."""
    return tuple(WindowPrediction(i,s,s+spec.window_size,_approved_predictor(s,s+spec.window_size))
        for i,s in enumerate(range(0,len(FIXTURE_SEQUENCE)-spec.window_size+1,spec.stride)))

def _summaries(windows):
    evaluable=tuple(x for x in windows if x.direction is not PredictedDirection.ABSTAIN)
    counts={d:sum(x.direction is d for x in evaluable) for d in
        (PredictedDirection.INCREASES_INCLUSION,PredictedDirection.DECREASES_INCLUSION)}
    modal=fraction=None
    if evaluable and len(set(counts.values()))>1:
        modal=max(counts,key=counts.get); fraction=counts[modal]/len(evaluable)
    pairs=[(a,b) for a,b in zip(windows,windows[1:]) if PredictedDirection.ABSTAIN not in (a.direction,b.direction)]
    flips=tuple(FlipPair(a.index,b.index,a.direction,b.direction) for a,b in pairs if a.direction!=b.direction)
    return evaluable,modal,fraction,pairs,flips

def evaluate_sliding_windows(*,fixture_id:str=FIXTURE_ID,window_size:int=3,stride:int=1)->SpanRobustnessReport:
    """Run only the frozen non-locus wiring fixture; arbitrary contexts are refused."""
    if fixture_id!=FIXTURE_ID: raise ValueError("fixture_id is not in the frozen method-fixture allowlist")
    spec=PerturbationSpec(window_size,stride)
    if window_size>len(FIXTURE_SEQUENCE): raise ValueError("window exceeds fixture")
    windows=_windows(spec); ev,modal,fraction,pairs,flips=_summaries(windows)
    return SpanRobustnessReport("prospective-span-robustness-method/3",PROSPECTIVE_BANNER,
        PROSPECTIVE_STATUS,METHOD_LABEL,CLAIM_BOUNDARY,FIXTURE_ID,FIXTURE_SEQUENCE,
        FIXTURE_SEQUENCE_SHA256,FIXTURE_ACCESSIBLE,FIXTURE_SILENCERS,spec,
        _approved_provenance(),windows,len(ev),len(windows)-len(ev),modal,fraction,len(pairs),
        flips,len(flips)/len(pairs) if pairs else None,(Lead("STMN2"),Lead("UNC13A")))
