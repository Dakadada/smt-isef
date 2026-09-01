# Engine-internal effector score characterization

This package characterizes the existing design engine's `manufacturability`
heuristic over the authenticated, committed SMN2 tiling. It re-runs the frozen
engine and derives candidate, chemistry/length, and class distributions from that
run. Rendering repeats the derivation from committed files and rejects summaries
that differ.

The output is a prospective description of software behavior. It is not an
experimental estimate of ASO, U7, or dCas13 performance, efficacy, potency, or a
biological effector-class effect. Architecture P1 remains off because its
multi-locus repeated-intervention data gate is not met.

The current committed context is 21 nt. Consequently, the engine emits ASO and U7
18/20-mers but emits no configured dCas13 22/24-mers. The report retains the
dCas13 base score and `grna_shadow=4` configuration while explicitly recording
that neither a dCas13 score distribution nor a shadow interval was instantiated.
