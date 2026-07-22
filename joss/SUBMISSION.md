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

JOSS lets you suggest reviewers by GitHub handle from its pool
(<https://reviewers.joss.theoj.org/>). Suggest 3–5, each **verified to be in the pool** and **free
of conflicts of interest** — per JOSS COI rules, not a recent co-author (last ~4 years), not the
same institution, not a current collaborator, advisor, or advisee.

> _Candidate shortlist is being compiled from the JOSS reviewer pool (radio-astronomy and
> reproducible-research backgrounds) and will be filled in here once each handle is verified._

## Before you click submit
- Confirm the Zenodo release metadata is correct and the DOI resolves.
- Sanity-check the repo against the JOSS review criteria — the `research-publish` skill's readiness
  check (`uv run python .claude/skills/research-publish/check_publish_readiness.py`) is all `[x]`
  except the inherently-manual items.
- Weigh the repo's youth (public since 2026-06-27) against JOSS's preference for a development
  history; the in-preparation RNAAS notes are evidence of active research use.
