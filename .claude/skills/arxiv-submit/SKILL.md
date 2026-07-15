---
name: arxiv-submit
description: Prepare, validate, and guide an arXiv paper submission — assemble the source/PDF upload package, capture every arXiv submission property (categories, license, comments, MSC/ACM, journal-ref, DOI, ORCID), and walk through the (manual) web upload. Use when the user wants to submit a paper to arXiv.
---

# arXiv submission helper

**Reality first (state this to the user):** arXiv has **no public API for *submitting* papers** — the
arXiv API is read-only (search/retrieve), and SWORD is for institutional repositories only. So this
skill cannot upload for the user. What it *does*: assemble a clean, valid submission package, capture
**all** submission-form properties so nothing is missed, validate against arXiv's rules, and give a
precise step-by-step for the web upload. The user performs the final upload at
<https://arxiv.org/submit>.

## When invoked

1. **Find the paper.** Ask for the paper directory if not given (default `paper/`). Identify the
   main `.tex` (the one with `\documentclass`) or, if PDF-only, the `.pdf`.
2. **Run the assembler:** `uv run python .claude/skills/arxiv-submit/assemble_arxiv.py --paper <dir>
   --out arxiv-submission`. It produces, under `arxiv-submission/`:
   - `arxiv-source.tar.gz` — the upload (LaTeX source + figures + `.bbl` + `\input`s), or the PDF.
   - `metadata.yaml` — every submission property, auto-filled where possible (including a
     **suggested primary category + cross-lists** — see *Choosing the categories* below),
     `TODO` otherwise.
   - `CHECKLIST.md` — the web-upload steps + validation results.
3. **Fill the metadata with the user.** Walk through `metadata.yaml` and confirm/complete each field
   (below). Extract title/authors/abstract from the `.tex` automatically; ask for the rest.
4. **Validate** (the assembler checks these; fix any failures before submitting):
   - LaTeX **source** is strongly preferred over PDF-only (PDF-only is allowed but flagged/penalised
     and can't be reprocessed). For AASTeX/BibTeX, **include the compiled `.bbl`** — arXiv runs
     `latex` but not always your `.bst`, so ship the `.bbl` and don't ship `.bib`-only.
   - Abstract ≤ ~1920 characters; no unescaped LaTeX-only macros in the abstract field (arXiv stores
     plain text).
   - All `\includegraphics`/`\plotone` figures are present and referenced by relative path (no
     absolute paths, no `..` escaping the source dir).
   - No `\input{generated/...}` pointing outside the tarball; include generated macros.
   - A `00README.json` (or `00README.XXX`) is included only if you need to set the build order / mark
     the toplevel file (usually unnecessary for a single `main.tex`).

## Every arXiv submission property (capture all in `metadata.yaml`)

- **title** — plain text (LaTeX math allowed sparingly).
- **authors** — full names + affiliations, in submission order; one author has the submitting
  account. (An AI/LLM is **not** an eligible author — credit it in acknowledgements + a software
  citation, never the author list.)
- **abstract** — plain text, ≤ ~1920 chars.
- **primary_category** / **cross_lists** — which arXiv groups the paper is submitted to. The
  assembler pre-fills a suggestion in `metadata.yaml`; confirm it using **Choosing the
  categories** below.
- **comments** — free text; conventionally page/figure count + status, e.g.
  `"4 pages, 3 figures. Code and data: github.com/joebarbere/jansky-research"`.
- **license** — the user must choose one at submission:
  - `CC BY 4.0` (recommended for open work) · `CC BY-SA 4.0` · `CC BY-NC-SA 4.0` · `CC0` (public
    domain) · `arXiv non-exclusive license to distribute` (most restrictive; default).
- **report_number** — institutional report id, if any (usually none for independent authors).
- **journal_ref** — only once published in a journal (added later via "journal-ref").
- **doi** — a journal/Zenodo DOI, if any (added later).
- **msc_class** / **acm_class** — only for math (`math.*`) / CS (`cs.*`) submissions.
- **orcid** — link the submitting author's ORCID (e.g. `0009-0008-3289-4447`) — strongly
  recommended; it's set on the arXiv account/profile, and `\author[orcid]{...}` is in the source.
- **endorsement** — note whether the submitting category needs endorsement (see Account setup).

## Choosing the categories (which groups to submit to)

arXiv routes a paper by **one primary** category (the single best fit — it sets the listing) plus
up to ~3 **cross-lists** (secondary categories whose readers should also see it). The assembler
pre-fills a suggestion in `metadata.yaml` — curated per slice, or keyword-inferred for a new paper
(flagged *VERIFY*). Confirm it with the user using the rules below. Radio astronomy lives almost
entirely under **astro-ph**:

| Code | Scope |
|---|---|
| `astro-ph.GA` | galaxies, the Milky Way, ISM, AGN/jets, HI 21 cm, Faraday rotation |
| `astro-ph.CO` | cosmology, large-scale structure, isotropy/dipole, source counts, lensing |
| `astro-ph.HE` | FRBs, pulsars/magnetars, transients, accretion/jets, explosive phenomena |
| `astro-ph.SR` | the Sun & solar radio bursts, stars, M/white/brown dwarfs, coherent emission |
| `astro-ph.EP` | planets & the solar system: Jovian/Saturnian/ice-giant radio, magnetospheres |
| `astro-ph.IM` | instrumentation, methods, software, pipelines, data analysis, statistics, RFI |

**Rules:**

1. **Primary = where the headline result lives.** A science result → its science domain (an HI
   rotation curve is `astro-ph.GA`; an FRB census is `astro-ph.HE`). A paper whose contribution
   *is* the tool/benchmark/DSP kernel/pipeline → `astro-ph.IM`.
2. **These papers are reproducible tools, so `astro-ph.IM` is almost always a cross-list** — but
   the *primary* stays the science domain unless the method itself is the point (the IM-primary
   ones here are `driftsearch`, `torchfdmt`, `torchdsp`, `rfitrend`, `ecallisto_pipeline`).
3. **Cross-list only for a genuinely different readership**, ≤3, most-relevant first. Add a second
   *science* domain when the result spans them: AGN variability = `GA` + `HE`; an isotropy/dipole
   test = `CO` + `GA`; a WD-pulsar = `SR` + `HE`.
4. **Non-astro cross-lists** only when another community should see it: `eess.SP` (a pure
   DSP/algorithm kernel), `physics.space-ph` (solar-wind/magnetosphere heliophysics), `physics.ins-det`
   (a hardware/receiver design), `stat.ML` (an SBI/ML method). Not routine.
5. **Endorsement is per-category** (see Account setup) — a first submission to a category you're not
   endorsed in needs an endorser there, so prefer a primary you can actually submit to.

**Per-slice quick reference** (the assembler's `SLICE_CATEGORIES` holds the full primary+cross map):

| Primary | Slices |
|---|---|
| `astro-ph.HE` | frbstats, frbperiod, frbwait, frblens, lpt, lptv, pulsarspec, ppdot, pte2, glitchpop |
| `astro-ph.GA` | hi, peaked, southern, offsets, vlbi, vlass, stacking, spectra, fashienv, rmsky, rmstructure |
| `astro-ph.CO` | rmdipole, sourcecounts |
| `astro-ph.SR` | solarbursts, windwaves, swaves, triangulate, type3synthesis, typeii, ecallisto_census, stokesv, stokesv_discovery, svsbi, wdpulsar |
| `astro-ph.EP` | junodam, skr, vgpra |
| `astro-ph.IM` | driftsearch, torchfdmt, torchdsp, rfitrend, ecallisto_pipeline |

## Web-upload steps (what to tell the user to do)

1. Sign in at <https://arxiv.org/login> → **Start a new submission** (<https://arxiv.org/submit>).
2. **License:** pick the one from `metadata.yaml`.
3. **Upload:** the `arxiv-source.tar.gz` (preferred) — arXiv compiles it; check the **AutoTeX
   preview PDF** matches `paper/main.pdf`. (If PDF-only, upload the `.pdf`; expect a "no source"
   flag.)
4. **Metadata:** paste title, authors, abstract, comments; set **primary** + **cross-list**
   categories; add report-no / MSC / ACM only if applicable.
5. **Preview & submit.** Note: new submitters often need **endorsement** in the chosen category, and
   there's a ~daily announcement cadence (submit before the deadline to appear next business day).

## Account setup (guide the user)

- Register at <https://arxiv.org/user/register> with a valid email; **verify** it.
- Add an **affiliation** (independent researchers use *"Independent researcher"*) and **link ORCID**
  in the account profile — ORCID + a consistent name is your durable identity.
- **Endorsement:** arXiv requires a first-time submitter to be *endorsed* for the category they
  submit to (e.g. `astro-ph.*`), unless auto-endorsed (auto-endorsement comes from affiliation/email
  domain history, which independent authors usually lack). To get endorsed: submit the paper to
  trigger an endorsement request, then ask an established arXiv author in that category to endorse
  you via the link arXiv emails/shows; or email the moderators. Plan for this — it can take days.
- There is **no fee**. Submissions are moderated; off-topic or low-quality (e.g. unverified
  AI-generated) content can be reclassified, held, or rejected.

Return to the user: the path to `arxiv-submission/`, the filled `metadata.yaml`, the validation
verdict, and the next concrete action (fix X, or "ready — upload at arxiv.org/submit").
