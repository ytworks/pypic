#!/usr/bin/env bash
# Runs campaign runs one after another on 8 GPUs (campaign-a.md). Usage, from the repository root:
#   research/topics/spherical-shock-injection/run_queue.sh NAME [NAME ...]
# Before each run it waits until the working tree is clean (a recorded run starts from a committed tree), and
# it skips a run whose output directory already holds a finished log. Outputs go to the NVMe through the link
# output/research/spherical-shock-injection.
set -u
S=output/research/spherical-shock-injection
C=research/topics/spherical-shock-injection/configs
for n in "$@"; do
  if [ -f "$S/logs/$n.log" ] && grep -q "^run exit 0" "$S/logs/$n.log"; then echo "skip $n (done)"; continue; fi
  while [ -n "$(git status --porcelain)" ]; do sleep 30; done
  echo "$(date +%H:%M:%S) start $n at $(git rev-parse --short HEAD)"
  .venv/bin/python -m pypic.launch -n 8 research/topics/spherical-shock-injection/run_piston.py --config "$C/$n.yaml" \
    --out "$S/$n" --progress 144 > "$S/logs/$n.log" 2>&1
  rc=$?
  echo "run exit $rc" >> "$S/logs/$n.log"
  echo "$(date +%H:%M:%S) end $n exit $rc; disk $(df -h /mnt/nvme | tail -1 | awk '{print $4}') free"
done
