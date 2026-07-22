# FAQ — using jansky-research, and how the research in it works

Common questions from researchers who want to use the toolkit or understand the papers in this
repository. This is the project's stance plus general open-science norms — **not legal advice**; if
a licensing question matters to you, take formal advice.

## Using the toolkit

### How do I use jansky-research in my own work?

Install it as a dependency and import the slice modules or call the CLIs
(`python -m jansky_research.<slice>`), and reuse the shared helpers (`data.py`, `pipeline.py`,
`report.py`) as building blocks. See [`usage.md`](usage.md). You **don't fork it to use it** — depend
on it like any Python library; your own analysis lives in your own project.

### Should I fork the repository?

Only if you're developing a change to the *toolkit itself* — a new tool/slice or a bug fix (see
[`CONTRIBUTING.md`](../CONTRIBUTING.md)). To *use* the toolkit in your research, depend on it; don't
fork.

### How do I cite it?

Cite the archived release via its Zenodo DOI (`10.5281/zenodo.21482378`, all-versions) — see
[`CITATION.cff`](../CITATION.cff) and the README badge. Once the software paper appears in JOSS, cite
that too.

## The papers in this repository

### Why does a software repository contain research papers?

This repo is two things at once: a **reusable toolkit** and the **maintainer's own research that
demonstrates it**. Each `papers/<slice>/` is an analysis authored by the maintainer (Joseph
Barbere), run on real public data, showing the toolkit end to end. They are the worked examples the
tool was built to produce — not a shared publication venue.

### I did research with the toolkit. Where does my paper go?

In **your own** repository/space, authored by **you**, citing jansky-research. The `papers/` here are
the maintainer's work. If you want to contribute a new *tool or slice* back to the toolkit, that's a
pull request (your code contribution is credited) — but your **science write-ups and their authorship
stay yours**, in your space.

### Could someone else publish the papers that are already in this repo?

Practically, no — and the distinction that matters is *copying* versus *misconduct*:

- **Licensing.** The repository carries the MIT license. MIT is a *software* license; applied to a
  repo it is generally read to permit reuse of the contents *with attribution*, though its fit for
  prose is imperfect (see [Licensing](#licensing)).
- **Misconduct.** Regardless of any license, submitting someone else's analysis to a journal as your
  own is **academic misconduct** — journals require original authorship, and it would be caught.
- **Provenance protects priority.** Every paper is authored, dated, and public in the git history,
  and the releases are archived with **timestamped Zenodo DOIs**. That public record establishes who
  did the work first far more strongly than keeping it private would.

So the papers are not "unprotected." Others are welcome to **build on** the work — cite it, extend it,
reproduce it — which is the point of open science; they just can't **re-present it as their own**.

### If I leave a paper unpublished (in a journal) forever, is it "up for grabs"?

No. A publicly posted, DOI-archived, dated analysis is **already disseminated and attributed** —
"not submitted to a journal" is not the same as "unauthored" or "public domain." Priority rests on
the public record. And leaving it public does **not** stop the author from submitting it to a journal
later: for astronomy venues, prior public posting (a repo, preprint, or Zenodo archive) is normal and
accepted — RNAAS and the arXiv-friendly journals expect it. (Check the specific journal's policy, but
astronomy has no "prior public posting disqualifies you" rule.)

## Licensing

### Is there one license for the code and the papers?

The repository carries the **MIT** license — ideal and permissive for the **code**. Because MIT is a
software license and the papers are creative works, how the *papers* may be reused is a softer
question than for the code. If that matters to you as a reader, the safe default is **cite and link
rather than copy**, and open an issue to ask. (The maintainer may add a separate license for
`papers/` — e.g. Creative Commons — or a NOTICE clarifying this.) *Not legal advice.*

## Contributing

### How do I contribute?

See [`CONTRIBUTING.md`](../CONTRIBUTING.md): tools, slices, and fixes via pull request; issues and
questions via the tracker. Your own research papers belong in your own project, citing the toolkit.
