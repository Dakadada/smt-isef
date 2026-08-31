# Working style: Claude + Codex

- Outsource almost ALL substantive work to Codex — implementation, investigation,
  analysis, research, data/reanalysis, plan drafting, and review. Claude
  orchestrates and relays concise results; Claude does not do heavy tool work or
  long analysis inline.
- Keep this main Claude chat CLEAN. Do NOT start substantive work unless the user
  explicitly says to get started. Default posture: wait for the go-ahead.
- Default Codex model is gpt-5.6-sol. Use luna only when the user explicitly asks.
- Delegate to Codex via the `codex:codex-rescue` subagent (`task`). Continue the
  last resumable thread for a given piece of work instead of starting fresh,
  unless the user asks for a clean slate.
- Delegate review to Codex, preferring `/codex:adversarial-review` over plain
  `/codex:review` unless a quick read-only pass is all that's needed.
- STANDING GATE: after every main deliverable (each completed task / artifact),
  run a Codex sol adversarial review for STRUCTURAL issues before advancing.
  Fix (or explicitly log) the findings; do not move to the next task until the
  review is clean or the open issues are recorded.
- Only truly trivial, obvious fixes may be done directly by Claude without
  delegating.
- Forking: as work scales, use Claude fork subagents for parallel Claude-side
  work; until then prefer Codex. (Fork = guaranteed Opus 4.8 + full context.)
