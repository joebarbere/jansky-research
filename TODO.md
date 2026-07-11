# Publishing TODO — jansky-research

The manual, account-bound steps to get the papers and software out. Nothing here can be
automated by the repo; each needs you signed in to an external service.

Identity to use everywhere: **Joe / Joseph Barbere**, ORCID **0009-0008-3289-4447**
(the papers and `CITATION.cff` use "Joseph"; keep the name consistent on each account).

---

## Recommended order (and why) — read this first

The steps are ordered by **dependency** (each unblocks the next) and by **effort/payoff** (start with
the one permanent, low-friction win). Do them top to bottom:

| # | Step | Why it's in this slot | Blocks |
|---|------|----------------------|--------|
| 0 | One-time setup (ORCID, affiliation) | Every downstream account wants these | everything |
| 1 | **Zenodo** DOI for the repo | A permanent citable DOI. **JOSS requires it.** Quick: flip a switch + cut a release | step 2 |
| 2 | **JOSS** software paper | Reviews the *software* (not novelty) — the natural home for the toolkit. Needs the step-1 DOI | — |
| 3 | **RNAAS** note (frbstats) | A short citable note; editorial screening only (no peer review). Independent of 1–2 | — |
| 4 | **arXiv** (lptv, junodam, fashienv; then frblens, rmstructure) | Preprints for the genuine-novelty papers (triage 2026-07-11); top 3 **packaged & ready**, `lptv` leads. Needs endorsement (lead time) | — |

**Plain-English walkthrough:**

1. **Do Zenodo first.** It's the only hard dependency in the chain (JOSS asks for the DOI at
   acceptance) and it's the fastest: connect Zenodo to the GitHub repo, cut a `v0.0.1` release, copy
   the DOI. One permanent, citable artifact for the whole toolkit.
2. **Then JOSS**, because you now have the DOI it needs. JOSS judges whether the *software* is a real
   scholarly tool — so lead with the validated analyses + the Airflow/Podman reproducibility layer,
   not a novelty claim.
3. **RNAAS and arXiv are independent of 1–2** and of each other — do them whenever. RNAAS is the
   quickest (no peer review); start arXiv early because the **endorsement** step has multi-day lead
   time on a first `astro-ph` submission.
4. **Order among the arXiv papers doesn't matter**, but if you submit only one first, make it
   `vlass` or `peaked` (the strongest methodology + recover-a-known stories).

**Scope as of now:** the repo has **twenty-six** completed slices, each with a tested tool + an
AASTeX paper under `papers/<slice>/` (the full list is in `README.md`). The earlier "reliable
no-auth data sources are largely used up" assumption is **retired** — the 2026-07 re-survey
(`survey/opportunity-scan-2026-07.md`) found major new/upgraded open sources (LoTSS DR3, CHIME/FRB
Cat 2, SPICE-RACS DR2, MeerKAT image products, VLASS SE images, OVRO-LWA solar, Juno/Waves, Parkes
PSRDA), and its ranked shortlist has since been **executed**: plans 33–37 are merged slices
(`stokesv_discovery`, `torchfdmt` — absorbing plan 28 — `lpt`, `rmstructure`, `junodam`). The
current opportunity file is **`fable-ideas.md`** (2026-07-05 deep re-scan; supersedes the scan's
shortlist, records corrections/closed doors). `plans/29-lotss-deep-144mhz-counts.md` remains
scoped (should migrate to LoTSS DR3 — `fable-ideas.md` F33). The publishing steps below operate
on the completed slices.

**`stokesv` is now complete** (CASDA recovered; the forced-photometry leg is wired). It is a
methods + honest-limits paper (I recovered from real RACS-low V images; single-epoch V is
variability-limited), so it belongs in the **repo + Zenodo** bucket, not the arXiv shortlist — its
result is an honest limitation, not a fresh discovery. Real numbers come from `make reproduce` with
`CASDA_USERNAME` set (the leg needs an OPAL login + `~/.casda_pw`).

**`ecallisto_census` is now complete** — the coverage-corrected type III occurrence census
(`rate = N_events / C` vs the SILSO sunspot number). It builds and *recover-a-known* validates the
census statistic (Pearson r=0.97, recovers the injected k=0.03), but the event stream is synthetic and
the SILSO series is the only real data; the real multi-cycle e-Callisto ingest is coverage-limited
future work (`survey/ecallisto-census-findings.md`). So it is a **repo + Zenodo** method paper, not an
arXiv discovery. GATE-2 passed with no blockers.

---

## 0. One-time setup

- [ ] Confirm your **ORCID** (0009-0008-3289-4447) is verified and public.
- [ ] Decide an affiliation string — the papers use **"Independent researcher"**; use the
      same on every service.
- [ ] `git tag v0.0.1 && git push --tags` is **not** done yet — hold it until Zenodo is
      wired (step 1), so the first tagged release is the one Zenodo archives.

---

## 1. Zenodo — archive + DOI (do this first)

Gives a permanent, citable DOI for the repo. JOSS requires this DOI at acceptance.

- [ ] Sign in at <https://zenodo.org> with GitHub; **log in with ORCID** linked.
- [ ] **GitHub** tab → flip the switch **ON** for `joebarbere/jansky-research`.
      (Zenodo only archives releases created *after* the switch is on.)
- [ ] On GitHub, **Releases → Draft a new release**: tag `v0.0.1`, title
      "jansky-research v0.0.1", publish.
- [ ] Back on Zenodo, confirm the new record appeared and that it picked up
      **`.zenodo.json`** (title, description, creator = Joseph Barbere, MIT, keywords).
      Fix anything in the Zenodo UI if needed.
- [ ] Copy the **concept DOI** (the version-agnostic one).
- [ ] Add the DOI badge to `README.md` (top): `[![DOI](https://zenodo.org/badge/DOI/<doi>.svg)](https://doi.org/<doi>)`.

---

## 2. JOSS — the software paper (`joss/paper.md`)

JOSS reviews the **software**, not novelty — the right home for this toolkit. Free.

Pre-submission checklist (JOSS will check these):
- [x] OSI-approved license in the repo (MIT — `LICENSE`). ✓
- [x] `joss/paper.md` + `joss/paper.bib` present, with an ORCID'd author. ✓
- [x] Tests + CI (`make test`, 85% coverage floor). ✓
- [ ] **Archive the software** (the Zenodo DOI from step 1) — JOSS asks for it.
- [ ] Skim the JOSS "substantial scholarly effort" bar
      (<https://joss.readthedocs.io/en/latest/submitting.html>) — a thin wrapper can be
      desk-rejected; lead the submission notes with the twenty-six slices (recover-a-known
      validations + methodology contributions) + the Airflow/Podman reproducibility layer.
- [ ] Submit at <https://joss.theoj.org/papers/new>: give the repo URL, the branch/tag,
      `joss/paper.md`'s path, and the Zenodo DOI.
- [ ] Respond to the editor/reviewers in the public review issue; expect requests for
      docs (a `CONTRIBUTING.md` and a short usage doc help).

---

## 3. RNAAS — the frbstats note (`papers/frbstats/rnaas.tex`)

A short, citable Research Note of the AAS. Editorial screening, not peer review.

- [ ] `make paper` (or build just the note) and check `papers/frbstats/rnaas.pdf`
      (≤1000 words, 1 figure — it currently is).
- [ ] Create/sign in to an AAS Journals author account at
      <https://journals.aas.org/research-notes/>; link ORCID.
- [ ] Confirm the **current RNAAS publication charge** in the author guidelines
      (there is a small fee; the amount changes — don't assume).
- [ ] Submit via the AAS peer-review system (Editorial Manager): upload `rnaas.tex`
      with `refs.bib`, the figure (`figures/dm_by_class.pdf`), and the generated
      `generated/macros.tex` (run `make figures` first so they exist), or inline the
      numbers if their system won't take `\input`.
- [ ] Author list: **Joseph Barbere** only; AI disclosure stays in the acknowledgements.
- [ ] This is **not** posted to arXiv.

---

## 4. arXiv — lptv, junodam, fashienv (packaged & ready) + frblens, rmstructure next

**Shortlist refreshed 2026-07-11** after a full-repo arXiv-worthiness triage (five parallel reviewers
over all 42 papers). The three genuine-novelty, real-data papers below are **packaged and validated**
— `papers/{lptv,junodam,fashienv}/arxiv-submission/` built with **0 assembler errors**, abstracts under
the 1920-char limit, figures + compiled `.bbl` included, `metadata.yaml` pre-filled with categories &
page counts:

- **`lptv`** (astro-ph.HE × astro-ph.SR) — the first uniform multi-epoch forced Stokes-V survey of the
  long-period-transient class: one secure single-epoch circular detection (ASKAP J1745−5051, 15%,
  21.6σ) + one candidate + a uniform V-limit table. *Frame as "first uniform V survey + serendipitous
  burst catches," not a discovery; the secure source is a known magnetic CV. 15/16 measured (one CASDA
  access-dropout — footnote it).*
- **`junodam`** (astro-ph.EP × astro-ph.IM) — a from-orbit Juno/Waves DAM census showing the apparent
  ~180× range-occurrence "law" is a 1/r² detection effect (→2.2× residual) and the Earth-canonical Io
  boxes do not organise the orbital vantage. *An honest debunking; the residual is entangled with
  beaming geometry (conceded, framed as an upper bound).*
- **`fashienv`** (astro-ph.GA × astro-ph.CO) — the first environment-split FASHI DR1 HI mass function;
  void knee suppressed 2.9σ, an independent FAST confirmation of the ALFALFA void trend. *Referee risk:
  the significance falls to ~1σ under an EdS distance frame (disclosed, defended via matched cosmology).*

Next after these: **`frblens`** (first catalogue-level lensed-repeater search + a real fraction limit)
and **`rmstructure`** (first structure function of SPICE-RACS DR2). The **old shortlist**
(`type3synthesis`, `vlass`, `peaked`, `triangulate`) is **demoted** — the triage rated them hold /
not-a-paper (recover-a-knowns or method demos with no discovery); keep them repo + Zenodo, not arXiv.
**`frbstats`, `torchdsp`, `torchfdmt`** are JOSS software candidates, not science preprints. Before any
upload, ensure the paper's real-data macros are populated (`make reproduce`, or the slice's real
driver — already done for these three).

The steps below are ordered for **acceptance**, not just submission — the two things that actually gate
a first independent-researcher arXiv paper are the **endorsement** (long lead time) and **moderation**
(your paper must read as a genuine contribution and be correctly categorised). Do them in this order:

- [ ] **(do first — long lead) Get an endorser.** A first submission to `astro-ph.*` from an account
      with no institutional email/history needs an **endorsement**. Identify a recently-published
      `astro-ph` author and ask (most will, for a sound paper); or register and let arXiv issue an
      endorsement code to pass on. This can take days–weeks, so start it before the paper is even final.
- [ ] Register + verify email at <https://arxiv.org/user/register>; affiliation **"Independent
      researcher"**, link ORCID (0009-0008-3289-4447).
- [ ] **Paper-ready gate** (do before uploading): run **`make reproduce`** so every macro is a
      real-data number (not the offline-synthetic CI value); confirm `make paper` compiles the chosen
      paper with no undefined refs; confirm the AI-use disclosure is present and the author is a real
      person. Pick the **single strongest paper first** (`lptv` — a positive detection; or `junodam`
      for the cleanest debunking) — one clean submission builds the track record that makes later ones
      routine.
- [ ] `make arxiv` → use the `arxiv-submission/` package under
      `papers/{lptv,junodam,fashienv}/` (already assembled; `metadata.yaml` pre-filled 2026-07-11):
  - [ ] Open the `metadata.yaml`; fill the remaining TODOs — **`primary_category`** and **`comments`**
        ("N pages, M figures") — from the category table below.
  - [ ] Upload `arxiv-source.tar.gz` (LaTeX **source** + `.bbl` included — arXiv re-runs AutoTeX, do
        not upload only a PDF); check the AutoTeX preview matches the local `main.pdf` exactly.
  - [ ] Set the license (default CC BY 4.0), the primary category, and cross-lists; paste
        title/authors/abstract/comments from `metadata.yaml`.
- [ ] Author list is **Joseph Barbere** only; an AI/LLM is not an author (it lives in the disclosure +
      `\software{}`). A real human author is required for endorsement and moderation to pass.
- [ ] **Submit, then clear moderation.** Submissions are held for moderator review; expect a possible
      **category reclassification** (accept it or reply with a one-line justification) and, for anything
      that reads as a reproduction, a request to sharpen the novelty in the abstract/comments. Correct
      categorisation (below) and a clear "what's new here" sentence in the abstract are what get a paper
      released rather than held — which is exactly why only the fresh-angle papers are on this list.

### Fill these `metadata.yaml` TODOs before any arXiv upload
- `primary_category` and `comments` (page count) are left as TODO on purpose.
- Recommended primary categories if you ever do submit a given paper:

| Paper | primary | reasonable cross-list |
|-------|---------|-----------------------|
| `frbstats` | `astro-ph.HE` | `astro-ph.IM` |
| `frbperiod` | `astro-ph.HE` | `astro-ph.IM` |
| `hi` | `astro-ph.GA` | `astro-ph.IM` |
| `driftsearch` | `astro-ph.IM` | `astro-ph.EP` |
| `spectra` | `astro-ph.GA` | `astro-ph.IM` |
| `offsets` | `astro-ph.GA` | `astro-ph.IM` |
| `vlass` | `astro-ph.HE` | `astro-ph.SR`, `astro-ph.IM` |
| `peaked` | `astro-ph.GA` | `astro-ph.HE`, `astro-ph.IM` |
| `triangulate` | `astro-ph.SR` | `astro-ph.IM`, `physics.space-ph` |
| `type3synthesis` | `astro-ph.SR` | `astro-ph.IM`, `physics.space-ph` |
| `ecallisto_pipeline` | `astro-ph.IM` | `astro-ph.SR` |
| `ecallisto_census` | `astro-ph.SR` | `astro-ph.IM` |

---

## Summary

| Artifact | Venue | Status |
|----------|-------|--------|
| Toolkit (the repo) | Zenodo DOI → JOSS | not submitted |
| `papers/frbstats/rnaas.tex` | RNAAS | not submitted |
| `papers/type3synthesis/main.tex` | arXiv (optional, lead) | not submitted |
| `papers/vlass/main.tex` | arXiv (optional) | not submitted |
| `papers/peaked/main.tex` | arXiv (optional) | not submitted |
| `papers/triangulate/main.tex` | arXiv (optional) | not submitted |
| `papers/ecallisto_pipeline/main.tex` | arXiv (optional, astro-ph.IM) | not submitted |
| `papers/frbstats/main.tex` | repo + Zenodo only (Airflow re-homed → ecallisto_pipeline) | — keep in repo |
| `papers/stokesv/main.tex` | repo + Zenodo only (methods + honest single-epoch limits) | — keep in repo |
| `papers/ecallisto_census/main.tex` | repo + Zenodo (method + recover-a-known; real census is future work) | — keep in repo |
| the other 14 `papers/<slice>/` (reproductions/negatives) | repo + Zenodo only | — keep in repo |
