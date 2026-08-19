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
   self-referential epistemics and reader address → deleted or restated as fact;
   limitations stay but move in-line with magnitudes. Do NOT touch `generated/`,
   `refs.bib`, macros, numbers, `\software{}`, or any claim's strength.
4. Gates, in order — all must pass:
   - `... prose_lint.py papers/<slice> --diff-guard` (prose-only edit, mechanical)
   - `... prose_lint.py papers/<slice>` (zero HIGH findings)
   - `uv run python scripts/triage_papers.py --paper <slice>` (zero HIGH/MED)
5. Show the user `git diff --word-diff papers/<slice>/main.tex` for approval.

For a whole-paper conversion, prefer delegating to the **style-editor** agent, which
runs this procedure end-to-end; check `git status` afterwards regardless (a reviewer
here once gutted a results JSON without editing a file).

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
