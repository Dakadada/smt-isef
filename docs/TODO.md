# Compute-native plan execution queue

Updated 2026-08-31 UTC. The detailed scientific constraints and evidence remain in
`docs/HANDOFF.md` and `docs/architecture.md`.

**OWNER: CODEX ONLY.** This is Codex's execution queue. Do not hand these items to
Claude or another orchestrator unless the user explicitly changes ownership.

- [x] Task #8 — prospective superior-molecule design set. Committed at `0cf3603`;
  focused tests 12/12 and checkpoint smoke suite 147/147. The honest superior set
  is empty because the committed off-target corpus is empty.
- [x] Task #9 — prospective in-silico cryptic-exon span-robustness method.
  STMN2/UNC13A are leads only; do not represent this as the conditional-on-data P4
  benchmark or emit uncommitted target designs. Implemented as a closed-world,
  non-locus engine-wiring fixture; focused tests 10/10 and independent gate clean.
- [x] Task #10 — in-silico characterization of design-engine effector-class score
  structure on the committed SMN2 tiling. Label all results engine-internal; do not
  infer a biological or data-backed class effect. Independently reviewed CLEAN;
  focused 12/12; artifacts reproduce byte-for-byte.
- [x] Task #11 — minigene context compression/representation tooling. Preserve an
  explicit minigene domain label; this is tooling, not biological validation.
  Independently reviewed CLEAN; focused 13/13; full target/configuration re-derivation
  and artificial-junction gates are enforced.
- [x] Run the standing independent review and fix/log findings after each of tasks
  #9, #10, and #11; run focused tests plus `make smoke-test`; checkpoint each clean
  increment. Final evidence: #9 10/10, #10 12/12, #11 13/13, unified 182/182;
  checkpoints `5cbb864` and `b6840c6`.
- [ ] Task #16 additional adversarial lock gate — deferred by prior user decision,
  not silently treated as completed.

Out of scope for this queue: archived `pre/` work, run logs, and `work/` data.
