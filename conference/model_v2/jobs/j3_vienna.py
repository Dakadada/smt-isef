"""Private JSONL worker: run under feng-env, never fold a genome on the login node."""
import json
import math
import sys

import RNA


def evaluate(request):
    context = request["context"]
    if request["method"] == "duplex":
        duplex = RNA.duplexfold(context, request["aso"])
        if duplex.energy >= 99999 or "(" not in duplex.structure:
            return {"status": "not_applicable", "reason": "no RNA duplex found"}
        energy = float(duplex.energy)
        if not math.isfinite(energy):
            raise ValueError("nonfinite duplex energy")
        return {"status": "ok", "result": {
            "energy_kcal_mol": energy, "structure": duplex.structure}}
    if request["method"] != "accessibility":
        raise ValueError("unknown method")
    start, length = request["start1"], request["length"]
    end = start + length - 1
    window = min(request["window"], len(context))
    span = min(request["span"], window)
    if length > window:
        return {"status": "not_applicable", "reason": "footprint exceeds RNAplfold window"}
    if not 1 <= start <= end <= len(context):
        raise ValueError("footprint outside context")
    probabilities = RNA.pfl_fold_up(context, length, window, span)
    # Like J2 and RNAplfold's _lunp file: segments END at index end (1-based).
    probability = float(probabilities[end][length])
    if not math.isfinite(probability) or not -1e-12 <= probability <= 1 + 1e-12:
        raise ValueError("invalid RNAplfold probability")
    # ViennaRNA summation roundoff can produce 1.0000000000000002.
    return {"status": "ok", "result": {"p_unpaired": min(1.0, max(0.0, probability))}}


def main():
    print(json.dumps({"tool": "ViennaRNA", "version": RNA.__version__, "temperature_c": 37}), flush=True)
    for line in sys.stdin:
        try:
            result = evaluate(json.loads(line))
        except Exception as exc:
            result = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
