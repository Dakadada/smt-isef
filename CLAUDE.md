# Working style: Claude + Codex

- Outsource almost ALL substantive work to Codex — implementation, investigation,
  analysis, research, data/reanalysis, plan drafting, and review. Claude
  orchestrates and relays concise results; Claude does not do heavy tool work or
  long analysis inline.
- Keep this main Claude chat CLEAN. Do NOT start substantive work unless the user
  explicitly says to get started. Default posture: wait for the go-ahead.
- Model/effort tiering (pass model+effort EXPLICITLY on every dispatch; default is sol):

  | Task | Codex | Claude fallback (when Codex exhausted) |
  |---|---|---|
  | Development / fixes | sol · high | opus 5 · deep (think-hard prompt) |
  | Trivial mechanical edits (renames, format-preserving tweaks, single-fn rewrites, test scaffolding) | luna · max | sonnet 5 |
  | First review (per arm) | sol · high | opus 5 · deep |
  | Re-reviews (later rounds) | terra · lower effort | opus 5 · light (targeted re-verify) |
  | Escalated re-review (significant structural) | sol · high | opus 5 · deep |
  | Monitoring / status | luna | haiku 4.5 |
  | Orchestrator / main chat | — | opus 4.8 (this session) |

  Reviews must NOT all be sol: first review sol, re-reviews terra, re-escalate to sol only on
  significant structural findings. Do NOT cut the NUMBER of rounds to save cost — tier the model down,
  keep the rounds. Claude has no effort knob → tune effort via prompt depth. opus 5 is a
  fire-and-return WORKER only; never hand it the driver seat / never converse with it (that is why the
  main chat is opus 4.8).
- Codex route — DEFAULT MCP, not the subagent. Call `mcp__codex__codex` (new
  task) / `mcp__codex__codex-reply` (continue thread) directly from the main
  chat; sandbox=workspace-write writes fine (old write-broken blocker CLEARED).
  MCP has no Claude driver layer, so it's leaner than the subagent — no second
  Claude instance spawned just to babysit Codex (which stalls waiting on
  notifications anyway). MCP fires NO harness notification → poll
  `jobs/<task>.json` status (arm a Monitor on it). Model+effort explicit every
  dispatch; bare `sol`/`luna`/`terra` REJECTED → use `gpt-5.6-sol`/`-terra`/`-luna`.
  Continue the last resumable thread for a FIX (cheap context reuse); only
  REVIEWS start fresh (reviewer independence).
- Use the `codex:codex-rescue` subagent ONLY when Codex output needs isolating
  from main context (noisy multi-file investigation) — it compresses the result
  before it hits main. Otherwise MCP.
- Token discipline: scope review prompts to the CHANGED files + the specific open
  findings, not "the whole diff" every round. Keep feng/HANDOFF.md trimmed to LIVE
  state + evergreen orientation — archive superseded checkpoints to
  feng/HANDOFF_ARCHIVE.md (a bloated HANDOFF is re-read in full every session).
- Delegate review to Codex, preferring `/codex:adversarial-review` over plain
  `/codex:review` unless a quick read-only pass is all that's needed.
- STANDING GATE: after every main deliverable (each completed task / artifact),
  run a Codex sol adversarial review for STRUCTURAL issues before advancing.
  Fix (or explicitly log) the findings; do not move to the next task until the
  review is clean or the open issues are recorded.
- Only truly trivial, obvious fixes may be done directly by Claude without
  delegating.
- Subagent delegation (Codex AND Claude fallback): give workers LIMITED, SCOPED context — a
  self-contained prompt, like a Codex dispatch — never the full conversation. NEVER use forks
  (subagent_type: fork inherits full context + is locked to opus 4.8 = wrong tier + expensive); use
  fresh Agent subagents with an explicit model override per the tiering table.
- Checkpoint rule: anything running ~1 hour, or planned to run that long, MUST be checkpointed to
  feng/HANDOFF.md before/while waiting, so a plain /clear is always safe. Keep HANDOFF current with
  live job/monitor/task IDs.
