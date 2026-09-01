# feng

`feng` is a clean-room prospective research arm for testing whether an inferred
effector-specific steric operator can transfer splice-modulation predictions between
modalities. Its claim boundary is deliberately narrow: this package will benchmark
candidate hybridization-associated splice effects in specified public experiments. It
does not establish clinical safety, CNS exposure, patient response, therapeutic
efficacy, or general performance beyond the preregistered datasets and endpoints.

The arm is implemented independently. Code in this package must not import from or
path-reference `src/`, `pre/src/`, `design_engine`, `offtarget`, `physics_adapter`, or
`design_panel`. Published methods may be independently re-derived, with provenance,
but legacy implementation code must not be copied or renamed into this package.

The `terminal/` directory remains sealed until the month-8 preregistration freeze. No
feature, model, benchmark, or arm implementation may inspect terminal results before a
freeze tag exists.

