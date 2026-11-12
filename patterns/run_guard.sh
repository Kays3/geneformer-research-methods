#!/bin/bash
# Skeleton of a resumable run driver. Start detached:
#   setsid nohup bash run_guard.sh > run.log 2>&1 < /dev/null &
# Every path comes from a required environment variable; nothing absolute is committed.
# Exit codes: 3 cap reached, 6 already running, 7 host not ready, 9 outside the time window,
#             10 an input differs from its pin, 12 the weights differ from their pin.
set -euo pipefail
: "${RUN_PYTHON:?not set: the interpreter of the analysis environment}"
: "${RUN_WORK:?not set: the working folder of this run}"
: "${RUN_MODEL_DIR:?not set: the model folder}"
: "${RUN_WEIGHTS_SHA256:?not set: the pinned sha256 of the weights file}"
: "${RUN_STOP_EPOCH:?not set: the end of the granted window, in seconds since the epoch}"
: "${RUN_CAP_GPU_S:?not set: the cap on charged GPU seconds}"
W="${RUN_WORK%/}"; HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$W"

# one run at a time; the exit code is kept in a file for whoever polls
[ -f "$W/run.pid" ] && kill -0 "$(cat "$W/run.pid")" 2>/dev/null && { echo "already running"; exit 6; }
echo $$ > "$W/run.pid"; rm -f "$W/run.exit"; trap 'echo $? > "$W/run.exit"' EXIT

# start checks: code, weights, inputs, host, window
(cd "$HERE" && sha256sum -c --quiet MANIFEST.sha256) || { echo "code differs from its manifest"; exit 10; }
got="$(sha256sum "$RUN_MODEL_DIR/model.safetensors" | cut -c1-64)"
[ "$got" = "$RUN_WEIGHTS_SHA256" ] || { echo "weights differ (a pointer file in place of the weights?)"; exit 12; }
(cd "$W/inputs" && sha256sum -c --quiet INPUTS.sha256) || { echo "an input table differs from its pin"; exit 10; }
[ -z "$(nvidia-smi --query-compute-apps=pid --format=csv,noheader)" ] || { echo "GPU not idle"; exit 7; }
[ "$(date +%s)" -lt "$RUN_STOP_EPOCH" ] || { echo "past the stop time"; exit 9; }

charged() { cat "$W/charged_s" 2>/dev/null || echo 0; }

# step NAME COMMAND...: skip if done; refuse past the window or the cap; charge the seconds used
step() {
  local name="$1"; shift
  [ -f "$W/done_$name" ] && return 0
  [ "$(date +%s)" -lt "$RUN_STOP_EPOCH" ] || { echo "time stop before $name"; exit 9; }
  [ "$(charged)" -lt "$RUN_CAP_GPU_S" ] || { echo "cap reached before $name"; exit 3; }
  local t0; t0="$(date +%s)"
  "$@"
  echo $(( $(charged) + $(date +%s) - t0 )) > "$W/charged_s"
  date -u +%FT%TZ > "$W/done_$name"
}

# gates first; a failed gate leaves no done file and stops the run (set -e)
step noop           "$RUN_PYTHON" "$HERE/your_noop_gate.py"   "$W"
step resolvability  "$RUN_PYTHON" "$HERE/your_resolve_gate.py" "$W"
step screen         "$RUN_PYTHON" "$HERE/your_screen.py"       "$W"   # checks the cap between chunks itself
step output_check   "$RUN_PYTHON" "$HERE/your_output_check.py" "$W"
step readout        "$RUN_PYTHON" "$HERE/your_readout.py"      "$W"   # CPU; written by the driver, not by hand
echo "run done"
