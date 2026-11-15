#!/bin/bash
# non-LLM heartbeat: writes status only. expires T0+48h.
D=/home/dashii/smt-isef/conference/model_v2/ops
[ "$(date -u +%s)" -gt "$(date -u -d 2026-11-17T06:52:49Z +%s)" ] && exit 0
mkdir "$D/.hb.lock" 2>/dev/null || exit 0; trap 'rmdir "$D/.hb.lock"' EXIT
J=$(squeue -u $USER -h -o '%i %j %T' 2>&1 | head -50 | tr '\n' ';')
printf '{"ts":"%s","host":"%s","squeue":"%s"}\n' "$(date -u +%FT%TZ)" "$(hostname)" "$J" > "$D/LAST_HEARTBEAT.json.tmp" && mv "$D/LAST_HEARTBEAT.json.tmp" "$D/LAST_HEARTBEAT.json"
# self-reschedule via at (no daemon); stops at expiry check above
echo "$D/heartbeat.sh >> $D/heartbeat.log 2>&1" | at now + 55 minutes >/dev/null 2>&1
