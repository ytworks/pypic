#!/usr/bin/env bash
# After-run analyses of a campaign run, then the pruning of its tracks and particle dumps (campaign-a.md).
#   research/topics/spherical-shock-injection/post_run.sh NAME          (a run with debris)
#   research/topics/spherical-shock-injection/post_run.sh NAME --null   (a null run)
# From the repository root. Every step logs to output/.../logs/NAME-post.log; a failed step stops before pruning.
set -u
n=$1; null=${2:-}
ROOT=$(pwd)
R=$ROOT/output/research/spherical-shock-injection/$n
L=$ROOT/output/research/spherical-shock-injection/logs/$n-post.log
PY=$ROOT/.venv/bin/python
cd research/topics/spherical-shock-injection || exit 1
step() { echo "== $(date +%H:%M:%S) $*" >> "$L"; "$PY" analyze.py "$@" "$R" >> "$L" 2>&1 || { echo "FAILED $*" >> "$L"; exit 1; }; }
if [ "$null" = "--null" ]; then
  step null
  "$PY" analyze.py prune --null "$R" >> "$L" 2>&1
else
  step front
  step drive
  step sites
  step counts
  step campaign
  "$PY" analyze.py prune "$R" >> "$L" 2>&1
fi
echo "== $(date +%H:%M:%S) post done" >> "$L"
