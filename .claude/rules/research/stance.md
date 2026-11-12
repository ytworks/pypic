# Research stance

This stance applies to all research work on research branches. It does not apply to engine or
code work; that follows main's rules.

Source: A. Morihata, "研究とは何か" (What is research), a memo for students at the University of Tokyo,
https://www.graco.c.u-tokyo.ac.jp/labs/morihata/research_memo.htm. It is summarised here in our own
words and applied to computational plasma physics. The author notes that practice differs between
fields; the principles below are the ones we adopt.

## What research is
- Research makes a novel result publicly known. Without novelty it is not research. Without
  publication it is not research either.
  - Runs, figures and code are not results until a write-up states them.
  - A topic is finished when the write-up is out, not when the runs are done.
  - Code without a document is only weakly public.
- Novelty means new to the field, not new to us. "We did not know it" or "our search found nothing"
  is not novelty.
  - Novelty cannot be proven. It is supported by a systematic survey that shows how far the
    related works go.
  - Ask every time: would this claim survive in front of the field's leading expert?
- The value of a result lies in how certain its novelty is, not how large it is. A small, certain
  advance beats a large, uncertain one.
  - When extending a result, do not ask "how do we make it more impressive?". Ask "where can we not
    yet say it is certain?".
- Knowledge is valuable in itself. Do not justify a topic by short-term applications.
- Studying is a means, not research. A summary of what we learned is not a result.

## The elements, kept in view at all times
Every topic must be able to answer, at any moment, what each element is or that it is still missing:
- **Background:** what every intended reader already knows and accepts.
- **Research question:** one question, small enough to answer within the topic, worth answering,
  and with an answer that is not already public.
- **Prior work:** the direct predecessors that the question improves on or combines.
  - Prefer questions that have prior work. If no one has worked on a question, it is likely
    unimportant, already solved, or too hard.
  - If a predecessor listed this as future work, ask why its authors left it (too hard? not
    interesting?) and what lets us do better.
- **Idea:** what makes the question answerable now. Question and idea come together: choose a
  question for which an idea exists. The ideal is a problem that is 80 % solved, with one gap
  left. If the difficulty is effort or computing time, say so instead of claiming an idea.
- **Evidence:** surveys, runs, analysis and theory that raise the certainty of the claim, with
  their criteria and null tests fixed before the runs.
- **Related work:** works that support the novelty and the value of the claim. This is not a list
  of what we read along the way.
- **Publication:** where and how the result becomes public. Peer review is preferred.

## Reading and searching
- **Find the question first.** Read a paper's introduction until its question and idea can be
  stated. If they cannot be stated, re-read the introduction before going further.
  - Read the evidence sections in detail only for direct prior work.
  - Note for every paper: whether it is peer reviewed, its venue (core or edge of the venue's scope),
    who sets the direction among the authors, and whether a later journal version exists.
- **Search for related work** by these routes:
  - citations, both backward and forward;
  - the field's own keywords;
  - 5–10 years of the venues where the topic sits at the centre of the scope (this also shows how
    much work exists);
  - later papers by the same authors;
  - asking experts.

## Writing
- The introduction states every element, including concrete results with numbers, and how they
  stand against prior work. Nothing is held back for later sections; a fuller introduction is better
  than a sparse one.
- A research plan is unfinished research: background, question, idea, prior work, some related
  work, and a novelty that can be named.

## In this repository
- **Topic README:** its sections are these elements; `research/tools/check_branch.py` checks
  their headings. Mark a section as open rather than filling it with guesses.
- **Claims:** each one names its records in `runs/` and its limits. A claim without a record is
  a hypothesis.
- **Review:** use the `research-reviewer` agent on the question, plan, results and drafts of a
  topic, and before every milestone merge into trunk. Never use it for code or engineering.
