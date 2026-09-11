#!/usr/bin/env bash
# Submit the Arm E H-off launcher from an already-created CODE commit.
# Usage: SCRATCH=/orcd/scratch/orcd/011/dashii \
#   feng/scripts/submit_arm_e_hoff.sh <full-code-commit-sha> <UTC-run-id>

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $0 <full-code-commit-sha> <YYYYMMDDTHHMMSSZ>" >&2
    exit 2
fi

SEALED_COMMIT="$1"
RUN_ID="$2"
: "${SCRATCH:?SCRATCH must name the approved scratch root}"
[[ "${SEALED_COMMIT}" =~ ^([0-9a-f]{40}|[0-9a-f]{64})$ ]] || {
    echo "ERROR: code commit must be a full lowercase SHA" >&2
    exit 2
}
[[ "${RUN_ID}" =~ ^[0-9]{8}T[0-9]{6}Z$ ]] || {
    echo "ERROR: run id must be a UTC basic timestamp" >&2
    exit 2
}

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "${REPO_ROOT}"
RESOLVED_COMMIT="$(git -C "${REPO_ROOT}" rev-parse --verify "${SEALED_COMMIT}^{commit}")"
[[ "${RESOLVED_COMMIT}" == "${SEALED_COMMIT}" ]] || {
    echo "ERROR: code commit does not resolve exactly" >&2
    exit 2
}
LAUNCHER_PATH="feng/scripts/arm_e_hoff.sbatch"
LAUNCHER_BLOB_SHA256="$(
    git -C "${REPO_ROOT}" show "${SEALED_COMMIT}:${LAUNCHER_PATH}" |
        sha256sum | awk '{print $1}'
)"

# stdin is the committed blob itself; no working-tree launcher bytes are sent.
git -C "${REPO_ROOT}" show "${SEALED_COMMIT}:${LAUNCHER_PATH}" |
    sbatch --export="ALL,SCRATCH=${SCRATCH},FENG_ARME_SEALED_COMMIT=${SEALED_COMMIT},FENG_ARME_LAUNCHER_BLOB_SHA256=${LAUNCHER_BLOB_SHA256},FENG_ARME_RUN_ID=${RUN_ID}"
