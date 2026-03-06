#!/bin/bash
# PreToolUse hook: block Write/Edit outside repo root or to protected directories

REPO_ROOT="/home/4/ud02114/workspace/9999_git_repositories/inverse_msmd"
PROTECTED_DIR="20251024_tyk2liang"

# Read JSON input from stdin
INPUT=$(cat)

# Extract file_path from tool_input
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty')

if [ -z "$FILE_PATH" ]; then
  exit 0
fi

# Resolve symlinks for comparison
RESOLVED_PATH=$(readlink -f "$FILE_PATH" 2>/dev/null || echo "$FILE_PATH")
RESOLVED_ROOT=$(readlink -f "$REPO_ROOT")

# Allow /tmp
if [[ "$RESOLVED_PATH" == /tmp/* ]]; then
  exit 0
fi

# Block if outside repo root
if [[ "$RESOLVED_PATH" != "$RESOLVED_ROOT"/* ]]; then
  jq -n --arg reason "Blocked: $FILE_PATH is outside the repository root" '{
    hookSpecificOutput: {
      hookEventName: "PreToolUse",
      permissionDecision: "deny",
      permissionDecisionReason: $reason
    }
  }'
  exit 0
fi

# Block if in protected directory
if echo "$RESOLVED_PATH" | grep -q "$PROTECTED_DIR"; then
  jq -n --arg reason "Blocked: $FILE_PATH is in protected directory '$PROTECTED_DIR'" '{
    hookSpecificOutput: {
      hookEventName: "PreToolUse",
      permissionDecision: "deny",
      permissionDecisionReason: $reason
    }
  }'
  exit 0
fi

exit 0
