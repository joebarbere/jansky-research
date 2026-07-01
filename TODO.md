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
| 4 | **arXiv** (type3synthesis, vlass, peaked, triangulate) | Preprints for the papers with a fresh angle; `type3synthesis` leads. Needs endorsement (lead time) | — |

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

**Scope as of now:** the repo has **nineteen** completed slices, each with a tested tool + an AASTeX
paper under `papers/<slice>/` (the full list is in `README.md`). Slice-building is **paused** — the
reliable no-auth data sources are largely used up; two larger efforts are scoped in
`plans/28-breakthrough-listen-singlepulse.md` and `plans/29-lotss-deep-144mhz-counts.md`. The
publishing steps below operate on the current nineteen; revisit the arXiv shortlist if #28/#29 land.

**`stokesv` is now complete** (CASDA recovered; the forced-photometry leg is wired). It is a
methods + honest-limits paper (I recovered from real RACS-low V images; single-epoch V is
variability-limited), so it belongs in the **repo + Zenodo** bucket, not the arXiv shortlist — its
result is an honest limitation, not a fresh discovery. Real numbers come from `make reproduce` with
`CASDA_USERNAME` set (the leg needs an OPAL login + `~/.casda_pw`).

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
      desk-rejected; lead the submission notes with the nineteen slices (recover-a-known
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

## 4. arXiv — the type3synthesis, vlass, peaked, and triangulate papers (optional)

Reserve arXiv for the papers with a genuinely fresh angle. **`type3synthesis` is the lead candidate**: a
solar type III electron beam tracked across four public instruments from the corona to 0.4 AU in one
reproducible drift-to-distance framework, with the genuinely new element being an *independent geometric*
check on the density-model distance (STEREO/WAVES and a STEREO-A+B triangulation analyse the **same
2013-05-15 event**; the plasma-frequency and geometric distances track in shape at r=0.989). The others:
**`vlass`** (a 703 deg² multi-epoch variability census with a real recovery — FK Comae Berenices — plus
the Quick-Look-systematics methodology), **`peaked`** (a three-frequency curvature selector with the
TGSS-upper-limit + resolution-floor method and two recover-a-known validations), and **`triangulate`**
(the 3D triangulation slice the synthesis builds on, standalone). **`frbstats` is removed** from the
shortlist — its only fresh claim was the Airflow pattern, now **re-homed** on a streaming archive as
`ecallisto_pipeline` (`plans/31`: a daily-schedule, catchup-backfill, per-station-fan-out Airflow-on-Podman
ingest of e-Callisto — a legitimate astro-ph.IM automation note); frbstats stays a JOSS/RNAAS tool. Do **not** post the pure reproductions/negatives as a
preprint batch. **Before any upload, run `make reproduce`** so the paper macros hold real-data numbers,
not the offline-synthetic CI values.

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
      person. Pick the **single strongest paper first** (`type3synthesis`) — one clean submission
      builds the track record that makes later ones routine.
- [ ] `make arxiv` → use the `arxiv-submission/` package under
      `papers/{type3synthesis,vlass,peaked,triangulate}/`:
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
| the other 14 `papers/<slice>/` (reproductions/negatives) | repo + Zenodo only | — keep in repo |
