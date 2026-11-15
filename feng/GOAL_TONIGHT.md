# /goal — Arm D + Arm E compute out tonight (STAGED 2026-09-04)

A keep-going directive, not a full spec. Fix details live in the HANDOFF 2026-09-04 audit block.

```
/goal Get Arm D + Arm E heavy compute launched tonight, fixed and re-reviewed first, with live
blocking monitors. Delegate substantive work to Codex sol (medium). Keep feng/HANDOFF.md current.

1. FIX the D + E HIGH blockers from the 5-arm audit (HANDOFF 2026-09-04 top block). Arm D core is
   sound — these are pre-freeze data-hygiene fixes (rebuild universe from 3 ScrMOE controls, cohort
   + score-cache guards). Arm E: de-double-count off-target burden, indeterminate on incomplete
   curve, durable seal, fix Pangolin H-off (aso+strand, no zero-fabrication), checkpoint fails loud.
2. GATE: feng tests green + a Codex sol adversarial re-review of the fixes. No launch until clean.
3. LAUNCH on preemptable + checkpoint: Arm D IntaRNA score sweep over the 8x101 strata; Arm E
   off-target seed-extend scan + Pangolin H-off sweep to populate therapeutic_window.
4. MONITORS: one lean blocking watcher per job, fires to MAIN at terminal. Verify job PID before
   trusting "running." Update HANDOFF with live job + monitor IDs after every action.

RULES:
- Seal stays UNOPENED / freeze UNCUT without explicit user sign-off. Sealed Phase-9 is one-shot,
  never autonomous. Do NOT manufacture a GO.
- COMPUTE IS NOT HARD-RESTRICTED. The core-hour estimates are rough — spend what's needed to reach
  honest verdicts. Respect HPC guard rails + quota (read /etc/motd.d/*, ~/orcd/.quota), but don't
  under-provision or stall a run just to stay under an estimate.
- If Codex hits its 5-HR rolling limit mid-work: set a monitor for 1 min after the reset and
  continue then — don't checkpoint-and-stop for the 5-hr case.
```

## Pending user inputs before launch
- (a) go-ahead to dispatch the fixes.
- (b) confirm Arm-C-scale compute on Arm E off-target scan despite its against-hypothesis signal.
- (c) whether Arm D seal fires tonight (default NO).

## Sources
- Arm C sig-test + review: `feng/artifacts/phase0_arm_c_probe/NUMBERS_AND_HEADLINES.md`
- 5-arm audit + ranking + fix file:lines: HANDOFF top block (2026-09-04).
