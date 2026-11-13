#!/usr/bin/env bash
# Runs post_run.sh on each named run as soon as its log ends with "run exit 0" (in the given order).
#   research/topics/spherical-shock-injection/post_watch.sh NAME[:null] ...
set -u
S=output/research/spherical-shock-injection
for item in "$@"; do
  n=${item%%:*}; kind=${item#*:}
  until [ -f "$S/logs/$n.log" ] && grep -q "^run exit" "$S/logs/$n.log"; do sleep 60; done
  if grep -q "^run exit 0" "$S/logs/$n.log"; then
    if [ "$kind" = "null" ]; then research/topics/spherical-shock-injection/post_run.sh "$n" --null
    else research/topics/spherical-shock-injection/post_run.sh "$n"; fi
    echo "$(date +%H:%M:%S) post $n: $(tail -1 "$S/logs/$n-post.log")"
  else
    echo "$(date +%H:%M:%S) $n failed; not analysed"
  fi
done
