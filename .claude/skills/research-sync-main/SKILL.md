---
name: research-sync-main
description: Take main into research/trunk after an engine change (merge --no-ff at a gated point), update the engine baseline in research/RESEARCH.md, then take trunk into the topic branches. Use when main has new engine commits, when the SessionStart report says origin/main is ahead, or when the user asks to update the research branches from main.
---

# Take main into the research branches

Main flows only one way: main → `research/trunk` → `research/<topic>`. Never merge main directly into a topic.

1. **Fetch, then look at what is new:**
   ```bash
   git fetch origin main research/trunk
   git switch research/trunk && git merge --ff-only origin/research/trunk
   git log --oneline research/trunk..origin/main
   git diff --stat research/trunk...origin/main -- pypic
   ```
2. **Choose the commit to take.** By default this is `origin/main`. Its `pypic/` tree must have an all-passed gate
   record:
   ```bash
   python -c "import sys; sys.path.insert(0, 'research/tools'); import check_branch as c; \
   print(c.gate_record(c.ROOT, 'origin/main'))"
   ```
   `None` means the engine at that commit has not been gated. Then take the newest gated commit of
   main instead, or ask the user.
3. **Merge:** `git merge --no-ff <commit> -m "Merge main <sha9> into research/trunk: <what changed in the engine>"`.
   - A conflict should not happen (`research/README.md`, "Why merges of main never conflict").
   - If one does: for main-owned files, take main's side (`git checkout --theirs <file>`). In
     `.claude/settings.json`, keep the hooks of both sides. Then report the cause to the user.
4. **Check:**
   - `python research/tools/check_branch.py`
   - `python -m pytest -q research/tests`
   - `ruff check research`
5. **Record it in `research/RESEARCH.md`:**
   - Rewrite the Engine baseline section: the commit and its gate record.
   - Add a Change Log entry `[progress] main <sha9> taken`, with the engine changes that matter to the topics.
   - Commit, then `git push origin research/trunk`.
6. **Take trunk into the topics.** List the topics that do not yet contain trunk:
   ```bash
   for b in $(git for-each-ref --format='%(refname:short)' 'refs/remotes/origin/research/*'); do
     git merge-base --is-ancestor research/trunk "$b" || echo "$b"; done
   ```
   For each topic the user wants updated:
   ```bash
   git switch research/<topic> && git merge --ff-only origin/research/<topic>
   git merge --no-ff research/trunk
   python research/tools/check_branch.py
   git push origin research/<topic>
   ```
   Runs recorded before the merge keep their commits. A re-run on the new engine is a new record.
