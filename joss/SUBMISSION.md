# JOSS submission — paste-ready notes

Helper for submitting `jansky-research` to the Journal of Open Source Software. The submission
itself is manual and account-bound (sign in with GitHub/ORCID at
<https://joss.theoj.org/papers/new>); this file just holds the fields and text to paste. The
software paper is `joss/paper.md`; the archived release is on Zenodo.

## Form fields

| Field | Value |
|---|---|
| **Repository URL** | `https://github.com/joebarbere/jansky-research` |
| **Software version** | `v1.0.0` |
| **Git branch** | `main` (leave default) |
| **Archive DOI** | `10.5281/zenodo.21482378` (concept / all-versions). At acceptance JOSS may want the version-specific archive of the released tag: `10.5281/zenodo.21482379` — give that if it asks for one exact release. |
| **Submitting author** | Joseph Barbere · ORCID `0009-0008-3289-4447` · Independent researcher |
| **Suggested keywords/topics** | Python; astronomy; radio astronomy; fast radio bursts; pulsars; reproducibility |

## Comments to the editor (paste into the message box)

> **jansky-research** is a tested, CPU-first Python toolkit (with optional opt-in GPU acceleration)
> for reproducible radio-astronomy analyses of public survey data. It is **not a thin wrapper or
> single-function utility**: it comprises **more than forty self-contained analysis "slices"**, each
> built to one discipline — a tested tool (pure NumPy/SciPy/Astropy, or pure PyTorch where device
> portability pays) → real public survey data → an adversarial science-review gate → an honest
> write-up whose every number regenerates from the pipeline.
>
> Two things mark the substantial scholarly effort:
>
> 1. **Recover-a-known validations + methodology, with honest nulls at scale.** Each tool must
>    recover an injected or historically established signal (e.g. FRB 20180916B's 16.35-day period,
>    the flat inner Milky Way rotation curve, the Crab's DM, published planetary rotation periods)
>    before being trusted on new data; a large fraction of slices then report bounded non-detections
>    or "the apparent signal is really a systematic" results rather than manufactured discoveries.
>    The toolkit doubles as a reproducible reference for what careful, gate-checked amateur analysis
>    of public data can and cannot establish.
> 2. **A dual reproducibility layer.** Static-input slices build through a server-less **Snakemake**
>    file-target DAG; the frequently-updated e-Callisto archive is ingested by an **Apache Airflow
>    DAG on rootless Podman** (daily schedule, backfill, per-station fan-out). One set of tested
>    entry points feeds the CLI, the notebooks, and both orchestrators, so `make reproduce`
>    regenerates data → analyses → papers → arXiv packages from a clean checkout.
>
> **Related publications (disclosure):** Two Research Notes of the AAS derived from this software are
> in preparation — a blind reanalysis of the Voyager 2 PRA ice-giant rotation periods
> (`papers/vgpra/`) and a cautionary note on TGSS×NVSS ultra-steep-spectrum selection
> (`papers/spectra/`). The repository additionally contains per-slice reproducible AASTeX write-ups
> under `papers/`; these are in-repository research records, not separate journal submissions. The
> JOSS contribution is the toolkit and its methodology, distinct from any individual science result.
>
> **AI-use disclosure:** The software was developed collaboratively with Anthropic's Claude (a large
> language model). The author directed the work and independently reviewed and verified all code,
> results, and citations; an AI/LLM is not an author. Full disclosure is in the paper's
> Acknowledgements.
>
> **Install note for reviewers:** the sole dependency `jansky` (the sibling teaching course) is
> resolved as a local path source, so please clone both repositories side by side (`jansky` and
> `jansky-research`) before `uv sync` — documented in `docs/usage.md`.

## Suggested reviewers

**Access note.** The JOSS reviewer pool (<https://reviewers.joss.theoj.org/>) is **editor-gated** —
GitHub-login required, no public browse or export. Two ways forward: suggest the verified handles
below, and/or ask your handling editor to search the pool (area **"Astronomy, Astrophysics, and
Space Sciences"** plus keywords like `pulsar`, `FRB`, `radio`, `reproducibility`, `workflow`).

Every handle below was **verified** by reading the linked `openjournals/joss-reviews` issue — each
is a confirmed reviewer of an accepted, closely-related JOSS paper (so demonstrably in the pool and
relevant).

**Recommended (suggest ~4–5):**

| Handle | Reviewed for JOSS | Fit for `jansky-research` |
|---|---|---|
| `@pravirkr` | *Your* — FRB/pulsar reader ([#2750](https://github.com/openjournals/joss-reviews/issues/2750)) | `sigpyproc` author; FRB/pulsar single-pulse & dedispersion — matches `frbstats`, `frbperiod`, `torchfdmt`/`singlepulse`, `pte2` |
| `@paulray` | *Your* ([#2750](https://github.com/openjournals/joss-reviews/issues/2750)) | Pulsar timing (PINT co-author, NICER) — matches `ppdot`, `pulsarspec`, `glitchpop`, `vgpra` |
| `@matteobachetti` | *Hasasia* — PTA sensitivity ([#1775](https://github.com/openjournals/joss-reviews/issues/1775)) | X-ray/time-series timing (Stingray/HENDRICS) **and** a prominent reproducible-research / Astropy contributor — spans the analysis **and** reproducibility pillars |
| `@ygrange` | *Virgo* — radio spectrometer ([#3067](https://github.com/openjournals/joss-reviews/issues/3067)) | ASTRON/LOFAR radio-data infrastructure & pipelines |
| `@garrettj403` | *Blimpy* — Breakthrough Listen I/O ([#1554](https://github.com/openjournals/joss-reviews/issues/1554)) | SETI/filterbank I/O — matches `driftsearch` + the Breakthrough Listen Voyager slices |

**Further verified options** (for specific sub-domain coverage): `@ptiede` or `@David-McKenna`
(VLBI/interferometry — *pyuvdata v3*, [#7482](https://github.com/openjournals/joss-reviews/issues/7482));
`@astrom-tom` / `@zhampel` / `@cmbiwer`. For a reproducibility-tooling-focused reviewer, the editor
can also search the pool for `snakemake`/`workflow`/`containers`.

**Conflict of interest.** Before suggesting anyone, confirm none is a recent co-author (~4 years),
same institution, or a current collaborator/advisor/advisee (JOSS COI rules). Availability isn't
guaranteed — suggesting 4–5 lets the editor choose.

## RNAAS submission (AAS Editorial Manager)

The two RNAAS notes (`papers/vgpra/rnaas.tex`, `papers/spectra/rnaas.tex`) are submitted separately
from JOSS, via the AAS **Editorial Manager** (article type **"Research Note"**). RNAAS is
editorially screened, not peer-reviewed; confirm the current AAS publication fee first. Build the PDF
(`make paper` or the per-note `tectonic` build) and upload it. Cover text to paste into the
comments/cover field:

**`vgpra` note — cover text**

> Dear Editors,
>
> Please consider the attached Research Note, *"A Blind Modern Reanalysis of the Voyager 2 PRA Radio
> Data Recovers Neither the Uranus nor the Neptune Rotation Period,"* for publication in RNAAS.
>
> The note reports a controlled null: a blind Lomb–Scargle reanalysis of the open PDS-PPI Voyager 2
> PRA data recovers a clean injected rotation in synthetic tests but recovers neither historical
> ice-giant rotation period — showing that the original determinations' emission-beaming and
> magnetic-longitude modelling was essential, and that a naive periodogram of these brief flybys
> constrains little. It is a concise methods-and-limits result well suited to RNAAS.
>
> The note is under 1000 words with a single figure, is not under consideration elsewhere, and uses
> only public data (NASA PDS-PPI). Every quantitative result is generated by open, tested software
> (`jansky-research`; archived at Zenodo, DOI 10.5281/zenodo.21482378) and regenerates from a clean
> checkout. Development used an AI assistant, disclosed in the acknowledgements; I directed and
> verified all work, am the sole author, and have no conflicts of interest.
>
> Thank you for your consideration.
> Joseph Barbere (ORCID 0009-0008-3289-4447), Independent researcher

**`spectra` note — cover text** (same structure; swap the first two paragraphs)

> …Please consider the attached Research Note, *"Ultra-Steep-Spectrum Selection from Raw TGSS ADR1 ×
> NVSS Is Dominated by the TGSS Flux-Scale Systematic,"* for publication in RNAAS.
>
> The note is a reproducible cautionary result: a small tested tool recovers the textbook mean radio
> spectral index (validation), but its ultra-steep-spectrum candidates do not survive a cross-check
> against the flux-scale-corrected de Gasperin et al. catalogue — the raw candidates are the known
> TGSS ADR1 flux-scale systematic, not a genuine ultra-steep population. …
>
> (Same closing paragraph: <1000 words, one figure, public data (TGSS/NVSS), pipeline-generated
> numbers, `jansky-research` provenance + Zenodo DOI, AI-use disclosed, sole author, no COI.)

## Before you click submit
- Confirm the Zenodo release metadata is correct and the DOI resolves.
- Sanity-check the repo against the JOSS review criteria — the `research-publish` skill's readiness
  check (`uv run python .claude/skills/research-publish/check_publish_readiness.py`) is all `[x]`
  except the inherently-manual items.
- Weigh the repo's youth (public since 2026-06-27) against JOSS's preference for a development
  history; the in-preparation RNAAS notes are evidence of active research use.
