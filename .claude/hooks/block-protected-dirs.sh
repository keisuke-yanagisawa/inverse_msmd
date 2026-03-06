#!/bin/bash
# PreToolUse hook: block Write/Edit outside repo root or to date-prefixed directories

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

deny() {
  jq -n --arg reason "$1" '{
    hookSpecificOutput: {
      hookEventName: "PreToolUse",
      permissionDecision: "deny",
      permissionDecisionReason: $reason
    }
  }'
  exit 0
}

is_outside_repo() { [[ "$1" != "$RESOLVED_ROOT"/* ]]; }
is_in_tmpdir() { [[ "$1" == /tmp/* ]]; }
is_date_prefixed() { echo "${1#$RESOLVED_ROOT/}" | grep -qE '(^|/)[0-9]{8}_'; }

FILE_PATH=$(cat | jq -r '.tool_input.file_path // empty')
[ -z "$FILE_PATH" ] && exit 0

RESOLVED_PATH=$(readlink -f "$FILE_PATH" 2>/dev/null || echo "$FILE_PATH")
RESOLVED_ROOT=$(readlink -f "$REPO_ROOT")

is_in_tmpdir "$RESOLVED_PATH" && exit 0
is_outside_repo "$RESOLVED_PATH" && deny "Blocked: $FILE_PATH is outside the repository root"
is_date_prefixed "$RESOLVED_PATH" && deny "Blocked: $FILE_PATH is in a date-prefixed protected directory"

exit 0
