# jansky-research — guide for Claude

**What this is.** Amateur radio-astronomy *research*, end to end. A public sibling of the
[`jansky`](https://github.com/joebarbere/jansky) teaching course (installed from its pinned git tag;
clone it at `../jansky` + `eval "$(make -s dev-env)"` for cross-repo work): where
jansky *teaches*, this repo *does original work* as reproducible "slices" — one gap → one tested tool
→ real public data → honest write-up. It **depends on `jansky` as a library** (`from jansky import …`)
and mirrors its `uv`/ruff/mypy/pytest conventions.

## The two repos (each should point at the other)

- **`../jansky`** — the course (library we depend on). Its `.claude/skills/` hold the general
  research helpers: `find-radio-papers`, `radio-source-lookup`, `dataset-watch`, `radio-mastodon`.
- **this repo** — the research. Its `.claude/skills/` hold the publishing/data helpers: `arxiv-submit`,
  `casda-cutout-fetch`, `pull-station-data` (pulls the `jansky-observe` station's codified
  observation bundles — averaged HI spectra + provenance — into `data/station/` for plan 78's
  `hline.read_capture`). (Ports of `find-radio-papers` + `radio-source-lookup` live here too so they
  work without the course checked out.)

A new session in *either* repo: read this file (or jansky's `CLAUDE.md`) to learn the other exists.

## Working rules (non-negotiable)

- **Branch before committing — never commit on `main`.** Squash-merge PRs; delete the branch.
- **85% coverage floor**, `ruff` (line-length 100) + `mypy` clean, before every PR.
- **Real-run outputs are committed evidence** (changed 2026-08-04 after the synthetic-clobber
  incident): `results/*.json|csv`, `papers/*/figures/*`, and `papers/*/generated/*` are tracked
  and reviewed in PRs. The offline Snakemake DAG (`make figures`) is CI smoke ONLY — never
  commit its synthetic outputs; `make guard-real` gates packaging. Compiled paper PDFs stay
  gitignored.
- **Honest framing**: validations, limits, and negatives reported plainly; no overclaiming. The
  `science-reviewer` agent gates this.
- **Versioning**: SemVer per [`VERSIONING.md`](VERSIONING.md) (version lives in `pyproject.toml`
  + `CITATION.cff`; Zenodo takes it from the tag). Every PR adds a `CHANGELOG.md` `Unreleased`
  entry; `python scripts/next_version.py` recommends the next bump from it.
- Commit footer: `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.
  PR footer: `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

## Where numbers go missing (lessons, 2026-08-10)

Every one of these produced prose that **still read as a finished claim with the number gone**
— never a crash, never a wrong value, just a hole that parses. Check for them explicitly; the
test suite cannot see any of them.

**A slice with two run modes has two sets of metrics, and one `macros.tex`.** The offline
synthetic validation and the real census (and CPU vs GPU) each emit the *other* mode's macros
as `--`, because their own metrics dict has no such keys. Whichever ran last silently blanked
the other's numbers. Four papers shipped abstracts citing blanked macros:
*"the detector separates type II from type III and RFI at purity --"*. The abstracts cite
**both** namespaces, so no single run can populate them — the file has to accumulate.
`report.preserve_live_macros` now enforces the invariant the two-namespace design assumed:
**a run may only add information; a real value always beats a placeholder, never the reverse.**
Call it from every `_write_macros`.

**Running an offline mode in the repo root destroys the real results JSON.** `run(".",
offline=True)` overwrites `results/<slice>_metrics.json` with synthetic output — verified live:
`typeii` lost 3429 lines and `is_real` flipped True→False. `make guard-real` catches this
*at packaging time*, which is late. **Run offline modes with `out=<tmpdir>`** and, if you need
the macros, call `_write_macros` on the real path afterwards. This is the mirror of the
2026-08-04 synthetic-clobber incident: same root cause (one artifact, two writers), opposite
direction.

**Merging is not enough — mode-dependent macros must be NAMESPACED.** `preserve_live_macros`
only arbitrates real-vs-placeholder. When one macro *name* means different things in the two
modes, both runs write a real value and there is nothing to arbitrate: `\tiiNEvents` meant
768 real observing days in the paper's prose and 48 synthetic events offline, so a rebuild
turned *"768 days, zero failures"* into *"48 days, zero failures"* — a wrong published number,
worse than the blank it replaced. Every mode-dependent macro is now `<slice>Syn*`/`<slice>Real*`.
Caught by the science reviewer, which reproduced the clobber on itself while checking
determinism.

**A bootstrap SE is not the uncertainty when the input is one random-field realization.**
`rmstructure`'s synthetic validation quotes 4.64 +/- 0.35 for an injected boost of 5. The
bootstrap resamples sources *within one fixed field*, so it misses the realization variance:
across 30 seeds the recovered ratio has mean **3.15** and std **1.11** — 3.2x the quoted SE —
with seed 0 (the published number) a high outlier at 4.64. Any injection test on a correlated
field must vary the field seed, not just resample within it.

**`--` in a macro means "the results JSON had null here".** It reads as an en-dash in a table
and as a hole in a sentence. The arXiv assembler now blocks on it; do not "fix" it with an
`arxiv.yaml` override, which hides a paper problem behind a packaging file.

**Regex traps that silently drop numbers**, all found in `assemble_arxiv.py` in one session:
`\newcommand{\x}{...}` matched with `[^{}]*` **skips every value containing braces** (any
exponent, subscript or `\mathrm{}`), and the macro then vanishes entirely — `p=\rfRealHUMAINP`
became `p=,`. `re.sub`'s replacement is a **template**, so a value containing a backslash
raises `bad escape` or expands a group reference; pass a callable. And `re.sub` with the `s`
flag makes `.` match newlines, which ate a whole config file.

**`\citet` is not deletable.** A textual citation is the grammatical subject of its sentence;
dropping it leaves an abstract starting *"derived the definitive modern..."*. It is resolved
from `refs.bib` now.

## The slice pattern (how every result is built)

tested helper (pure NumPy/SciPy/astropy + synthetic offline fixture) → real-data run (network,
`# pragma: no cover`) → GATE-2 science review → AASTeX paper (`papers/<slice>/`) → arXiv package
(`make arxiv`). Plan in `plans/NN-*.md`; findings in `survey/<slice>-findings.md`.

Every new slice paper / RNAAS note `\software{}`-cites **`jansky-research`** (the toolkit, Zenodo
concept DOI `10.5281/zenodo.21482378`) alongside `jansky` — copy the `@misc{janskyresearch}`
`refs.bib` entry from `papers/vgpra/` or `papers/spectra/`. See the `research-publish` skill.

**Merged slices:** forty plus the type III synthesis — the authoritative per-slice table
(tool + outcome) is in `README.md`. Each has a paper under `papers/<slice>/`. Publishing steps
(Zenodo → JOSS → RNAAS → arXiv) are tracked in Joe's personal notes, outside this repo
(Obsidian vault: `efforts/radio_astronomy/research_paper_todo.md`).

## Active direction (2026-07)

- **Pick the next slice from `fable-ideas.md`** (2026-07-05, a 12-agent deep re-scan; supersedes
  the shortlist in `survey/opportunity-scan-2026-07.md`, whose Tier-1 items are now merged:
  `stokesv_discovery`, `lpt`, `rmstructure`, `torchfdmt`, `junodam`). Suggested first moves
  there: F4 (WD-pulsar sweep), F8 (FASHI DR2). Executed so far: F1 → `rmdipole` (plan 38);
  F2+F5 → `frbwait`+`frblens` (plans 39+42; Cat 2 mirrored from CANFAR DOI 10.11570/25.0066 —
  chime-frb.ca itself is still 503); F6 → `torchdsp` (plan 43; CHIME baseband + ROCm GPU legs
  done). All fable-ideas have plans (`plans/38`–`84`).
- **Standing GATE-0 for anything from that file:** the scan session couldn't fetch primary
  sources (egress-blocked), so do a full-text novelty pass + a data-URL check before writing the
  plan. Its "Corrections & closed doors" section lists ideas already killed — check it first.
- Earlier blockers are resolved: CASDA recovered (`stokesv` complete, both legs; the discovery
  census is merged as `stokesv_discovery`); the southern GLEAM-X×RACS curvature catalogue is
  merged (`southern`).

## Layout

`src/jansky_research/` (slice modules + `data.py`/`pipeline.py`/`report.py`) · `tests/` · `plans/` ·
`fable-ideas.md` (current plan-ready idea list) · `survey/` (committed findings) ·
`papers/<slice>/` · `station/` (build guides for the physical rooftop
station — self-collected data, WIP; owner's working notes live in an Obsidian vault, not here) ·
`workflow/Snakefile` (static-slice file-DAG, drives `make figures`) · `airflow/` (streaming e-Callisto
ingest, Podman DAG) · `.claude/` (agents + skills) ·
`Makefile` (`setup`/`test`/`cov`/`lint`/`typecheck`/`figures`/`paper`/`arxiv`/`reproduce`).
