# 89 — `stylecorpus` / `traditional-style`: empirical pre-LLM prose for the papers

Status: 🚧 Stage 1 in progress (2026-08-17). Approved session plan; user decisions recorded below.

## Context

All 44 papers pass the science gates (triage, presenter/referee) but read stylistically
modern-AI: em-dash parentheticals, editorializing subsection titles ("…and what it is worth"),
rule-of-three constructions, `\emph{}`-led rhetorical mini-headings, self-referential epistemics,
reader-address, 25+-line methodology-narrating abstracts, inconsistent I/we. `paper-referee`
explicitly declines style commentary, so nothing gates this today.

Goal: derive traditional style **empirically** from the pre-LLM radio-astronomy literature
(1933–2021) — structure, tone, formatting — and ship a `traditional-style` **skill +
style-editor agent pair** that (a) converts existing papers and (b) guides future ones.
Meta-slice: plan → tested module → real corpus run → committed evidence → deliverable.

**User decisions (2026-08-17):** full-text sample ~2,000 papers (~5 GB, gitignored `data/`);
Joe supplies a free NASA ADS API token (`ADS_API_TOKEN` env var, never committed) for the
1933–1990 literature; deliverable = skill + agent pair. Defaults: pre-LLM cutoff = submitted
before 2022-01-01 (safely pre-ChatGPT); no OCR initially (qualitative pass reads scanned PDFs
directly; add tesseract only if the pre-1991 era needs numbers); "done" = one validated pilot
conversion, remaining 43 papers are a follow-on campaign in PRs of 4–6.

## Stages

1. **Scoping + total-PDF-size ballpark** (branch `style-corpus-scoping`): enumerate "all radio
   astronomy papers pre-2022" via ADS decade×journal count queries + arXiv per-year counts;
   stratified ~50-papers/era size sample (HTTP HEAD/Range for arXiv-era, sample-PDF download for
   ADS-scan era); estimate = Σ count × mean size with bootstrap CI. Code: pure logic in
   `src/jansky_research/stylecorpus.py` (tested), CLI `scripts/style_corpus_scoping.py`,
   network `# pragma: no cover`. Evidence: `results/stylecorpus_scoping.json` +
   `survey/stylecorpus-findings.md`.
2. **Acquisition** (branch `style-corpus-acquire`): maximal ADS/arXiv metadata+abstract harvest →
   `data/style_corpus/metadata/`; ~2,000 full texts stratified by era×journal×doctype
   (1933–59: 60; 60s/70s/80s: 150 each; 90s: 250; 2000s: 400; 2010–15: 400; 2016–21: 440),
   preferring arXiv LaTeX source; manifest committed as `results/stylecorpus_manifest.json`.
3. **Analysis**: quantitative fingerprints (sentence length, passive rate, I/we, section-title
   tables, abstract length, em-dash/`\emph`/list rates, hedging, rule-of-three,
   self-reference/reader-address) over corpus **and** our own 44 papers →
   `results/stylecorpus_fingerprints.json` + `results/stylecorpus_selfscan.json`; parallel
   subagent qualitative reads (~10 papers per decade×journal cell) → synthesized style guide.
4. **Deliverables** (branch `traditional-style-skill`):
   `.claude/skills/traditional-style/{SKILL.md, prose_lint.py, references/{style-guide,
   fingerprints, before-after}.md}` (lint = thin CLI over stylecorpus; `--diff-guard` fails if
   macros/`\cite*` keys/numeric literals/`\software{}`/`generated/` change) +
   `.claude/agents/style-editor.md` (sonnet, Edit-capable) + `.claude/README.md` bullet.
5. **Validation** (branch `style-pilot-<slice>`): pilot on atlas3i or fashienv; gates in order:
   diff-guard clean → triage zero HIGH/MED → fingerprints inside corpus IQR →
   presenter/referee round with no new science findings → blind A/B (fresh subagent picks the
   pre-2022-sounding version, ≤2 residual markers) → Joe approves the word-diff. Verdict in
   `survey/stylecorpus-findings.md`.

Full approved plan text: session plan 2026-08-17 (declarative-bouncing-raven).
