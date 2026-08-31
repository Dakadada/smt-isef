#!/usr/bin/env bash
set -uo pipefail

# Operator-tunable policy and site paths.
POLL_SECONDS=300
CODEX_POLL_SECONDS=30
CODEX_TIMEOUT_SECONDS=720
MAX_RESUBMITS=6
MAX_FIXES=2
MAX_IDENTICAL_FAILURES=3
REPO_ROOT=/home/dashii/smt-isef
SCRATCH_ROOT=/orcd/scratch/orcd/011/dashii
SBATCH_SCRIPT=scripts/real_run/run_stageA.sbatch
SBATCH_PARTITION=mit_data_transfer
CODEX_COMPANION=/home/dashii/.claude/plugins/cache/openai-codex/codex/1.0.6/scripts/codex-companion.mjs
PYTHON_BIN=/usr/bin/python3.12
WORK_DIR_REL=data/reanalysis

usage() {
    printf 'usage: %s JOB_ID RUN_DIR\n' "${0##*/}" >&2
    exit 2
}

[[ $# -eq 2 ]] || usage
cur=$1
RUN_DIR=$2
[[ -n "$cur" && -n "$RUN_DIR" ]] || usage

mkdir -p "$RUN_DIR" 2>/dev/null || {
    printf 'cannot create/access RUN_DIR: %s\n' "$RUN_DIR" >&2
    exit 2
}
RUN_DIR=$(cd "$RUN_DIR" 2>/dev/null && pwd -P) || {
    printf 'cannot resolve RUN_DIR: %s\n' "$RUN_DIR" >&2
    exit 2
}

LOG_FILE=$RUN_DIR/babysitter.log
ESCALATE_FILE=$RUN_DIR/BABYSITTER_ESCALATE
DONE_FILE=$RUN_DIR/BABYSITTER_DONE
PREVFAIL_FILE=$RUN_DIR/.babysitter_prevfail
FAILURES_FILE=$REPO_ROOT/$WORK_DIR_REL/fetch_failures.tsv
COUNTERS_FILE=$RUN_DIR/.babysitter_counters
SIGNATURE_FILE=$RUN_DIR/.babysitter_signature

log() {
    printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >>"$LOG_FILE" 2>/dev/null ||
        printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" >&2
}

resubmit_count=0
fix_count=0
last_signature=
same_signature_count=0

if [[ -r "$COUNTERS_FILE" ]]; then
    read -r resubmit_count fix_count <"$COUNTERS_FILE" || {
        resubmit_count=0
        fix_count=0
    }
fi
[[ $resubmit_count =~ ^[0-9]+$ ]] || resubmit_count=0
[[ $fix_count =~ ^[0-9]+$ ]] || fix_count=0
if [[ -r "$SIGNATURE_FILE" ]]; then
    IFS=$'\t' read -r same_signature_count last_signature <"$SIGNATURE_FILE" || {
        same_signature_count=0
        last_signature=
    }
fi
[[ $same_signature_count =~ ^[0-9]+$ ]] || same_signature_count=0

save_state() {
    printf '%s %s\n' "$resubmit_count" "$fix_count" >"$COUNTERS_FILE.tmp" 2>/dev/null &&
        mv -f "$COUNTERS_FILE.tmp" "$COUNTERS_FILE" 2>/dev/null
    printf '%s\t%s\n' "$same_signature_count" "$last_signature" >"$SIGNATURE_FILE.tmp" 2>/dev/null &&
        mv -f "$SIGNATURE_FILE.tmp" "$SIGNATURE_FILE" 2>/dev/null
    return 0
}

newest_slurm_log() {
    local job=$1 candidate
    candidate=$(find "$REPO_ROOT" "$RUN_DIR" -maxdepth 1 -type f -name "stageA-${job}.out" \
        -printf '%T@ %p\n' 2>/dev/null | sort -rn | head -n 1 | cut -d' ' -f2-)
    printf '%s' "$candidate"
}

sacct_report() {
    sacct -j "$1" --format=JobID,State,ExitCode,Elapsed 2>&1 || true
}

write_escalation() {
    local reason=$1 report slurm_log tmp
    if [[ -e "$ESCALATE_FILE" ]]; then
        log "ESCALATE already exists; preserving it and exiting: reason=$reason job=$cur"
        exit 1
    fi
    report=$(sacct_report "$cur")
    slurm_log=$(newest_slurm_log "$cur")
    tmp=$ESCALATE_FILE.tmp.$$
    {
        printf 'reason: %s\n' "$reason"
        printf 'current job id: %s\n' "$cur"
        printf '\nsacct -j %s --format=JobID,State,ExitCode,Elapsed:\n%s\n' "$cur" "$report"
        printf '\nlast 40 lines of newest matching slurm log:\n'
        if [[ -n "$slurm_log" && -r "$slurm_log" ]]; then
            printf 'log: %s\n' "$slurm_log"
            tail -n 40 "$slurm_log" 2>&1 || true
        else
            printf '(no matching stageA-%s.out found)\n' "$cur"
        fi
        if [[ -f "$FAILURES_FILE" ]]; then
            printf '\ncontents of %s:\n' "$FAILURES_FILE"
            cat "$FAILURES_FILE" 2>&1 || true
        fi
    } >"$tmp" 2>/dev/null
    if [[ ! -e "$ESCALATE_FILE" ]]; then
        mv "$tmp" "$ESCALATE_FILE" 2>/dev/null || true
    else
        rm -f "$tmp" 2>/dev/null || true
    fi
    log "ESCALATE reason=$reason job=$cur resubmits=$resubmit_count fixes=$fix_count"
    exit 1
}

write_done() {
    local kind=$1 summary=$2 tmp
    if [[ -e "$DONE_FILE" ]]; then
        log "DONE already exists; exiting successfully job=$cur"
        exit 0
    fi
    tmp=$DONE_FILE.tmp.$$
    {
        printf 'final job id: %s\n' "$cur"
        printf 'result: %s\n' "$kind"
        printf 'summary: %s\n' "$summary"
    } >"$tmp" 2>/dev/null
    if [[ ! -e "$DONE_FILE" ]]; then
        mv "$tmp" "$DONE_FILE" 2>/dev/null || true
    else
        rm -f "$tmp" 2>/dev/null || true
    fi
    log "DONE result=$kind job=$cur summary=$summary"
    exit 0
}

check_thresholds() {
    save_state
    if (( resubmit_count > MAX_RESUBMITS )); then
        write_escalation "resubmit count exceeded ${MAX_RESUBMITS}"
    fi
    if (( fix_count > MAX_FIXES )); then
        write_escalation "auto-fix count exceeded ${MAX_FIXES}"
    fi
}

record_signature() {
    local signature=$1
    if [[ "$signature" == "$last_signature" ]]; then
        ((same_signature_count += 1))
    else
        last_signature=$signature
        same_signature_count=1
    fi
    save_state
    log "failure signature count=$same_signature_count signature=$signature"
    if (( same_signature_count >= MAX_IDENTICAL_FAILURES )); then
        write_escalation "same failure ${MAX_IDENTICAL_FAILURES}x"
    fi
}

reset_signature() {
    last_signature=
    same_signature_count=0
    save_state
}

resubmit() {
    local why=$1 output new_id
    ((resubmit_count += 1))
    check_thresholds
    log "resubmit attempt=$resubmit_count reason=$why old_job=$cur"
    output=$(
        cd "$REPO_ROOT" 2>/dev/null || exit 1
        export SCRATCH="$SCRATCH_ROOT"
        export RUN_DIR
        sbatch --parsable --partition="$SBATCH_PARTITION" \
            --export=ALL,REPO="$PWD",RUN_DIR="$RUN_DIR" "$SBATCH_SCRIPT" 2>&1
    )
    new_id=$(printf '%s\n' "$output" | tail -n 1 | sed 's/;.*//' | tr -d '[:space:]')
    if [[ $new_id =~ ^[0-9]+([_.][A-Za-z0-9_.-]+)?$ ]]; then
        log "resubmit accepted old_job=$cur new_job=$new_id output=$output"
        cur=$new_id
        save_state
        return 0
    fi
    log "resubmit command failed old_job=$cur output=$output"
    write_escalation "sbatch resubmission failed"
}

failure_accessions() {
    if [[ -r "$FAILURES_FILE" ]]; then
        awk -F '\t' 'NR == 1 && tolower($1) ~ /accession|run/ { next }
            NF && $1 !~ /^[[:space:]]*$/ { gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1); print $1 }' \
            "$FAILURES_FILE" 2>/dev/null | LC_ALL=C sort -u
    fi
}

handle_partial() {
    local current_set
    current_set=$RUN_DIR/.babysitter_current_failures.$$
    failure_accessions >"$current_set" 2>/dev/null || true
    log "PARTIAL marker found; failure_accessions=$(paste -sd, "$current_set" 2>/dev/null || true)"
    if [[ -f "$PREVFAIL_FILE" ]] && cmp -s "$current_set" "$PREVFAIL_FILE"; then
        rm -f "$current_set" 2>/dev/null || true
        write_escalation "repeated fetch failures unresolved after retry"
    fi
    mv -f "$current_set" "$PREVFAIL_FILE" 2>/dev/null ||
        write_escalation "could not save partial failure accession set"
    reset_signature
    resubmit "Stage-A partial; retrying uncheckpointed accessions"
}

run_codex_fix() {
    local traceback=$1 launch task_id status_output elapsed status_log prompt
    status_log=$RUN_DIR/codex-fix-status.$$.log
    prompt="A long SLURM Stage-A run failed with this traceback tail:
$traceback

Fix the minimal code bug in src/rawread_reanalysis. Keep make smoke-test green. Do not touch pre/."
    log "dispatching Codex minimal-fix task"
    launch=$(node "$CODEX_COMPANION" task --background --write "$prompt" 2>&1)
    task_id=$(printf '%s\n' "$launch" | sed -nE \
        's/.*([Tt]ask[ _-]?[Ii][Dd]|taskId)[^A-Za-z0-9_-]*([A-Za-z0-9_-]+).*/\2/p' | tail -n 1)
    if [[ -z "$task_id" ]]; then
        task_id=$(printf '%s\n' "$launch" | tail -n 1 | awk '{print $NF}')
    fi
    if [[ -z "$task_id" || "$task_id" == *[!A-Za-z0-9_-]* ]]; then
        log "Codex task launch failed output=$launch"
        return 1
    fi
    log "Codex task launched id=$task_id output=$launch"
    elapsed=0
    while (( elapsed < CODEX_TIMEOUT_SECONDS )); do
        status_output=$(node "$CODEX_COMPANION" status "$task_id" 2>&1)
        printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$status_output" >>"$status_log" 2>/dev/null || true
        if printf '%s\n' "$status_output" | grep -Eiq '^[[:space:]]*completed[[:space:]]*$|(^|[[:space:]"'"'])status["'"']?[[:space:]]*[:=][[:space:]]*["'"']?completed(["'"']|[[:space:]}]|$)'; then
            log "Codex task completed id=$task_id"
            return 0
        fi
        if printf '%s\n' "$status_output" | grep -Eiq 'failed|cancelled|canceled|error'; then
            log "Codex task failed id=$task_id status=$status_output"
            return 1
        fi
        sleep "$CODEX_POLL_SECONDS"
        ((elapsed += CODEX_POLL_SECONDS))
    done
    log "Codex task timed out id=$task_id after=${CODEX_TIMEOUT_SECONDS}s"
    return 1
}

handle_dead_job() {
    local report state slurm_log tail_text lower class detail signature traceback
    report=$(sacct_report "$cur")
    state=$(printf '%s\n' "$report" | awk 'NR > 2 && $1 !~ /[.+]/ && $2 != "" {print $2; exit}')
    [[ -n "$state" ]] || state=UNKNOWN
    state=${state%%+*}
    slurm_log=$(newest_slurm_log "$cur")
    if [[ -n "$slurm_log" && -r "$slurm_log" ]]; then
        tail_text=$(tail -n 120 "$slurm_log" 2>&1 || true)
    else
        tail_text=
    fi
    lower=$(printf '%s\n%s\n' "$state" "$tail_text" | tr '[:upper:]' '[:lower:]')
    class=TRANSIENT
    detail="state=$state"

    if [[ "$state" == OUT_OF_MEMORY* ]] || grep -Eqi 'out of memory|oom-kill|oom kill' <<<"$tail_text"; then
        class=LARGE; detail=OUT_OF_MEMORY
    elif grep -Eqi 'no space left|disk quota exceeded|quota exceeded' <<<"$tail_text"; then
        class=LARGE; detail=$(grep -Eim1 'no space left|disk quota exceeded|quota exceeded' <<<"$tail_text")
    elif grep -Eqi 'sha256(sum)?:.*(failed|mismatch)|checksum.*(failed|mismatch)|computed checksum.*(did not match|expected)|: FAILED$' <<<"$tail_text"; then
        class=LARGE; detail=$(grep -Eim1 'sha256(sum)?:.*(failed|mismatch)|checksum.*(failed|mismatch)|computed checksum.*(did not match|expected)|: FAILED$' <<<"$tail_text")
    elif { grep -Eqi 'assert(ionerror)?.*(config|validation)|config(uration)? validation.*(assert|failed)|validate_config.*assert' <<<"$tail_text" ||
           { grep -Eqi 'validate_config|validation.*config|stageA\.targeted\.json' <<<"$tail_text" &&
             grep -Eqi 'assertionerror|assert ' <<<"$tail_text"; }; }; then
        class=LARGE; detail="config validation assertion"
    elif [[ "$state" == TIMEOUT* ]] && ! find "$REPO_ROOT/$WORK_DIR_REL/checkpoints" -type f -print -quit 2>/dev/null | grep -q .; then
        class=LARGE; detail="TIMEOUT after full walltime with no progress"
    elif grep -q 'Traceback (most recent call last)' <<<"$tail_text"; then
        if grep -Eqi 'curl|fetch|download|network|connection|timed out|temporary failure|http[^ ]* (429|5[0-9][0-9])' <<<"$tail_text"; then
            class=TRANSIENT; detail="network/fetch traceback"
        else
            class=CODE_ERROR; detail="python traceback"
        fi
    elif [[ "$state" == NODE_FAIL* || "$state" == PREEMPTED* ]] ||
         grep -Eqi 'curl: \([0-9]+\)|connection (reset|refused)|temporary failure|network is unreachable|timed out|preempt|node_fail' <<<"$lower"; then
        class=TRANSIENT; detail="state/network transient: $state"
    fi

    if [[ "$detail" == "state=$state" ]]; then
        signature="$class|$state|$(printf '%s' "$tail_text" | sed -E 's/[0-9]+/N/g' | tail -n 20 | sha256sum | awk '{print $1}')"
    else
        signature="$class|$state|$(printf '%s' "$detail" | tr '\n\t' '  ')"
    fi
    log "job left squeue job=$cur state=$state class=$class detail=$(printf '%s' "$detail" | tr '\n' ' ') log=${slurm_log:-none}"
    record_signature "$signature"

    case $class in
        LARGE)
            write_escalation "$detail"
            ;;
        CODE_ERROR)
            traceback=$(printf '%s\n' "$tail_text" | tail -n 80)
            if run_codex_fix "$traceback"; then
                log "running smoke test with python=$PYTHON_BIN"
                if cd "$REPO_ROOT" 2>/dev/null && PATH="$(dirname "$PYTHON_BIN"):$PATH" make smoke-test >>"$LOG_FILE" 2>&1; then
                    ((fix_count += 1))
                    check_thresholds
                    log "smoke test passed auto_fix_count=$fix_count"
                    resubmit "Codex minimal fix passed smoke-test"
                else
                    log "smoke test failed after Codex task"
                    write_escalation "codex auto-fix failed"
                fi
            else
                write_escalation "codex auto-fix failed"
            fi
            ;;
        TRANSIENT)
            resubmit "$detail"
            ;;
    esac
}

log "START initial_job=$cur run_dir=$RUN_DIR repo=$REPO_ROOT work_dir=$WORK_DIR_REL poll=${POLL_SECONDS}s"
while true; do
    if [[ -e "$ESCALATE_FILE" ]]; then
        log "existing BABYSITTER_ESCALATE detected; exiting non-zero"
        exit 1
    fi
    if [[ -e "$DONE_FILE" ]]; then
        log "existing BABYSITTER_DONE detected; exiting zero"
        exit 0
    fi
    if [[ -e "$RUN_DIR/STAGE_A_COMPLETE" ]]; then
        write_done COMPLETE "Stage-A completion marker found; all accessions completed."
    fi

    queue_status=$(squeue -h -j "$cur" -o '%i %T %M %R' 2>&1)
    queue_rc=$?
    log "squeue job=$cur rc=$queue_rc status=$(printf '%s' "$queue_status" | tr '\n' ';')"
    if (( queue_rc != 0 )); then
        log "squeue query failed transiently; retaining current job and retrying"
        sleep "$POLL_SECONDS"
        continue
    fi
    if [[ -n "$queue_status" ]]; then
        sleep "$POLL_SECONDS"
        continue
    fi

    log "job no longer present in squeue job=$cur"
    if [[ -e "$RUN_DIR/STAGE_A_COMPLETE" ]]; then
        write_done COMPLETE "Stage-A completed cleanly after leaving the queue."
    elif [[ -e "$RUN_DIR/STAGE_A_PARTIAL" ]]; then
        handle_partial
    else
        handle_dead_job
    fi
    sleep "$POLL_SECONDS"
done
