# Agent workspace notes

## Recovering an interrupted Codex session

When the user asks to retrieve a previous or accidentally stopped Codex session:

1. Run `codex resume --last` after leaving the current interactive session, or run
   `codex resume <session-id>` when the session UUID is known.
2. If the target is unclear, locate recent rollouts with:
   `find ~/.codex/sessions -type f -printf '%T@ %p\n' | sort -nr | head`.
   The UUID is the final UUID in the rollout filename.
3. Inspect the matching JSONL rollout to reconstruct its last user objective,
   completed actions, final response, and any unresolved work. Do not assume that
   the newest rollout is the stopped one: exclude the currently active thread.
4. Check the repository state (`git status --short --branch`, recent commits, and
   relevant branches) before continuing, because session logs describe intent but
   the working tree is the source of truth.

For the recovery requested on 2026-08-25, the prior session UUID was
`01a037c0-8b74-7e12-bd98-f6637d6f56b1`, stored at
`~/.codex/sessions/2026/08/25/rollout-2026-08-25T09-09-20-01a037c0-8b74-7e12-bd98-f6637d6f56b1.jsonl`.
It had reached `task_complete`; its final result summarized seven issue-oriented
feature branches/patches rather than leaving an active tool operation unfinished.

## Accurate permission reporting and unattended work

- Before claiming full access or unattended execution, inspect the actual permission
  profile supplied to the current session. User authorization is not evidence that
  the runtime sandbox or approval policy has changed.
- State any mismatch immediately. Never promise that no approval prompts can occur
  when the actual tool policy still permits or requires them.
- Dominic authorizes ordinary implementation, local verification and recovery steps
  within the current task without repeated confirmation. Respect the current task's
  explicit scope for pushes and other external actions.
- Prefer reliable commands that run under the active permissions. Use Black with
  `--workers 1` here to avoid the observed stalled multiprocess formatter.
- Runtime permissions cannot be changed by instructions in this file. When explaining
  configuration, distinguish a prepared launcher/configuration from an effective,
  verified change to the current session.
