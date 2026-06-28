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
| 4 | **arXiv** (frbstats, vlass, peaked) | Preprints for the 3 papers with a fresh angle. Needs endorsement (lead time) | — |

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

**Not ready yet:** the `stokesv` (RACS Stokes-V) slice is still in progress — its live forced-photometry
run is blocked on a CASDA outage (see `survey/stokesv-findings.md`). It is **not** in the table above;
add it once the slice merges a paper.

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
      desk-rejected; lead the submission notes with the five validated analyses + the
      Airflow/Podman reproducibility layer.
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

## 4. arXiv — the frbstats, vlass, and peaked papers (optional)

Reserve arXiv for the **three** papers with a genuinely fresh angle: `frbstats` (the Airflow-on-Podman
reproducibility pattern), **`vlass`** (a 703 deg² multi-epoch variability census with a real
recovery — FK Comae Berenices — plus the Quick-Look-systematics methodology), and **`peaked`** (a
three-frequency curvature selector with the TGSS-upper-limit + resolution-floor method and two
recover-a-known validations: high purity vs MHz-peaked, 100% recovery of a known HFP sample). Do
**not** post the pure reproductions/negatives as a preprint batch.

- [ ] Register + verify email at <https://arxiv.org/user/register>; add affiliation,
      link ORCID.
- [ ] **Endorsement:** a first submission to `astro-ph.*` usually needs an endorser.
      Submitting triggers the request; line up an established astro-ph author, or email
      the moderators. Plan for a few days.
- [ ] `make arxiv` → use `papers/frbstats/`, `papers/vlass/`, and `papers/peaked/arxiv-submission/`:
  - [ ] Open each `metadata.yaml`; fill the remaining TODOs (see the table below).
  - [ ] Upload `arxiv-source.tar.gz` (LaTeX source + `.bbl` included); check the AutoTeX
        preview matches the local `main.pdf`.
  - [ ] Set license (default in the file: CC BY 4.0) and categories; paste
        title/authors/abstract/comments from `metadata.yaml`.
- [ ] Author list is **Joseph Barbere**; an AI/LLM is not an author (it's in the
      disclosure + `\software{}`).

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
| `vlass` | `astro-ph.HE` | `astro-ph.SR`, `astro-ph.IM` |
| `peaked` | `astro-ph.GA` | `astro-ph.HE`, `astro-ph.IM` |

---

## Summary

| Artifact | Venue | Status |
|----------|-------|--------|
| Toolkit (the repo) | Zenodo DOI → JOSS | not submitted |
| `papers/frbstats/rnaas.tex` | RNAAS | not submitted |
| `papers/frbstats/main.tex` | arXiv (optional) | not submitted |
| `papers/vlass/main.tex` | arXiv (optional) | not submitted |
| `papers/peaked/main.tex` | arXiv (optional) | not submitted |
| `papers/{frbperiod,hi,driftsearch,spectra}/` | repo + Zenodo only | — keep in repo |
