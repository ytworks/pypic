---
name: research-new-topic
description: Start a new research topic (theme) on its own branch research/<topic>, cut from research/trunk, with research/topics/<topic>/ scaffolded from research/templates/topic and its README filled in with the user. Use when the user wants to start, open or set up a new research topic or theme on the research branches.
---

# Start a research topic

1. **Choose the slug** with the user: `[a-z0-9-]+`, short. It must not be `trunk` and must not already exist
   (`git branch -a --list 'research/*'`, `ls research/topics`).
2. **Cut the branch from an up-to-date trunk.** Never cut it from main or from another topic.
   ```bash
   git fetch origin research/trunk
   git switch -c research/<slug> origin/research/trunk
   git branch --unset-upstream   # never track trunk or main: push with -u to its own name
   ```
3. **Scaffold:** `cp -r research/templates/topic research/topics/<slug>`, then replace the placeholders.
4. **Fill in the README with the user.** Write each section in their words, following
   `.claude/rules/research/stance.md`:
   - Background;
   - Research question;
   - Prior work;
   - Idea;
   - Related work and novelty;
   - Approach, with its criteria and null tests fixed before any run.

   Leave a section marked open rather than inventing content.
   Add the works found to `research/literature/` and relate them in the topic's `literature.md`.
   Then have the `research-reviewer` agent review the README, and show the review to the user.
5. **Check:** `python research/tools/check_branch.py --worktree`.
6. **Commit:** `research(<slug>): start topic`. Push with `git push -u origin research/<slug>`.
7. **Register on trunk.** Add the topic's row to the Topics table of `research/RESEARCH.md`, on
   `research/trunk` only. Do it now on trunk, or at the topic's first merge into trunk. Topic
   branches do not edit `RESEARCH.md`.
