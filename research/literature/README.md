# Literature

One file per work: `<citekey>.md`. The citekey is the first author's family name, the year and the first
significant word of the title, all in lowercase: `spreiter1966hydromagnetic`. Add a letter on a
clash.

A topic can add a work, so these files may change on topic branches. Before adding one, check
trunk and the other topics for the same work. Each topic relates the works to its own question in its `literature.md`.

```markdown
# <Authors> (<year>): <title>

- Ref: <venue, volume, pages>; DOI or arXiv id
- Read: YYYY-MM-DD, <abstract | skimmed | full>; topics: <topic>, ...

## Claim
<!-- What the work claims, in one or two sentences. -->

## Method
<!-- Model, set-up, numerical method, parameters. -->

## Numbers to compare against
<!-- Values with their equation, figure or table number. -->

## Assessment
<!-- What the work shows against what it states; its limits; what it leaves open. -->
```
