---
name: traditional-style
description: Convert or draft a paper in traditional pre-LLM radio-astronomy journal style — empirically derived from a 2,000-paper 1933–2021 corpus — stripping AI prose tells (em-dash pivots, \emph openers, editorializing titles, bloated abstracts) while mechanically protecting every number, macro, and citation. Use when writing or revising anything under papers/, when prose "reads like AI", or before a submission pass.
---

# Traditional style — pre-LLM journal prose

Rewrites prose only. The science is untouchable: every generated-macro invocation,
`\cite*` key, numeric literal, and the `\software{}` block must survive the edit
exactly (multiset-equal), and `--diff-guard` enforces that mechanically.

This skill carries a `references/` subdirectory (a first for this repo's skills —
the evidence base is too large for one file):

- `references/style-guide.md` — the conversion rules, every one carrying corpus
  evidence (a fingerprint percentile or a verbatim quote with its bibcode/arXiv id).
- `references/before-after.md` — worked conversions of this repo's own passages.
- Quantitative baselines live in `results/stylecorpus_fingerprints.json` (1,137
  pre-LLM papers; per-era + overall percentiles) — the linter reads it directly.

## Converting an existing paper

1. Read `references/style-guide.md` in full, then the paper's `main.tex`.
2. Baseline: `uv run python .claude/skills/traditional-style/prose_lint.py papers/<slice>`
   — the findings are the work list, worst first.
3. Rewrite prose only, per the guide: abstract toward the corpus ~160-word median and
   result-first; em-dash parentheticals → parentheses or plain sentences; `\emph`
   rhetorical openers → said in words; editorializing section titles → nouns;
   **the paper title too** (guide §S4b: noun phrase ≤~18 words; 0.9% of 113,829 corpus
   titles are sentences, RNAAS 6.0%) — but a title is the author's most visible wording,
   so propose the new title in the PR body for sign-off rather than silently rewording;
   self-referential epistemics and reader address → deleted or restated as fact;
   limitations stay but move in-line with magnitudes. Do NOT touch `generated/`,
   `refs.bib`, macros, numbers, `\software{}`, or any claim's strength.
   **Run-in labels often carry ranking.** A `\textbf{Label.}` or `\emph{Label.}` opener
   frequently encodes how to weight what follows ("(validation)", "(a near-null)", "An
   honest null:"). Dissolving it into a topic sentence, or promoting it to a
   `\subsection{}`, changes that weighting: a heading that names an attempted measurement
   reads as a claimed one, and a one-sentence aside given heading rank acquires the
   visual weight of a result. Carry the label's scope into the prose, or keep the label.
   **Never relocate a numeric literal across sections to satisfy the diff-guard.** The
   guard compares multisets, so a moved number passes by construction — and in batch 2 a
   relocated `$\pm$230\,s` landed beside a ratio it contradicted, while its own abstract
   lost the clause that made the step mandatory. If shortening requires dropping a number,
   drop the sentence, not the number's home.
4. Gates, in order — all must pass:
   - `... prose_lint.py papers/<slice> --diff-guard` (prose-only edit, mechanical)
   - `... prose_lint.py papers/<slice>` (zero HIGH findings)
   - `uv run python scripts/triage_papers.py --paper <slice>` (zero HIGH/MED)
5. Show the user `git diff --word-diff papers/<slice>/main.tex` for approval.

For a whole-paper conversion, prefer delegating to the **style-editor** agent, which
runs this procedure end-to-end; check `git status` afterwards regardless (a reviewer
here once gutted a results JSON without editing a file).

## RNAAS notes

Notes are a different genre, with their own baseline (`results/stylecorpus_rnaas.json`,
391 pre-LLM notes) and guide section (`references/style-guide.md` §8). Lint and guard with
`prose_lint.py papers/<slice> --file rnaas.tex --genre rnaas` (and `--diff-guard`). Key
deltas from papers: ~830-word median (journal cap 1,500), a REQUIRED short abstract
since 2020-05-01 (~100 words; the baseline's abstract p50 = 0 averages across that
policy break — see guide §8), title carries the verdict, no closing recap, voice runs
hotter.

## Drafting new prose

Apply the guide from the start: open with the object/instrument ("We present <data>
of <object> with <instrument>"), state results with numbers by the third abstract
sentence, hedge only interpretations (once, causally, then commit), weave limitations
in with magnitudes, and let plain connectives repeat. Run the linter before review.

## Honesty is not negotiable

This repo's framing rules (validations, limits, negatives reported plainly) are fully
compatible with the corpus: pre-LLM papers state their limits in-line, with numbers,
in declarative sentences. Conversion changes the register, never the claim strength —
if a sentence gets stronger during rewriting, that is a defect, not style.
