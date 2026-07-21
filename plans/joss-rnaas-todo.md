# joss-rnaas-todo — v1.0.0 release, versioning policy, JOSS-paper refresh, RNAAS notes

Status: 📋 planned (2026-07-21). A multi-task publishing plan, not a single slice. Tracks the work to
get `jansky-research` release-ready for JOSS, codify versioning, and package the first RNAAS notes.
Publishing sequencing is coordinated with other in-progress slices in the personal publishing todo
(`efforts/radio_astronomy/research_paper_todo.md`, outside this repo).

## Context / why now

- The JOSS submission needs a **tagged, Zenodo-archived release**. Version decision: **`v1.0.0`**,
  not the `v0.0.1` the publishing todo currently prescribes. Rationale: a JOSS paper asserts "this is
  a real, citable scholarly tool"; `0.0.1` (SemVer for "pre-alpha, first commit") contradicts a
  feature-complete toolkit with ~40 tested slices, CI, an 85% coverage floor, and three
  reproducibility paths. Publishing is exactly the SemVer 1.0.0 moment (a public API you commit to),
  and it is the JOSS norm. `v1.0.0` also lets the paper describe the real scope.
- **The JOSS paper (`joss/paper.md`) is a stale snapshot** — dated 27 June 2026, it describes *six*
  modules (`frbstats`, `frbperiod`, `hi`, `driftsearch`, `spectra`, `vlass`) and calls the toolkit
  "CPU-only". The repo now has ~40 merged slices (46 slice modules in `src/jansky_research/`),
  **including GPU/ROCm slices** (`torchfdmt`, `torchdsp`, `svsbi`) that break the "CPU-only" claim. It
  must be refreshed before submission. Detailed findings in **T2**.
- **No versioning infrastructure exists** — no git tags, no CHANGELOG, no policy. Codify it now (T1)
  so the 1.0.0 bump is principled and every future bump is deterministic and Claude-decidable.
- **RNAAS notes:** package **`vgpra`** first — its recover-a-known→honest-null arc is the purest
  compact showcase of the house method. **`spectra`** (a "the apparent signal is really a systematic"
  result) is a natural follow-up. Both need a short-form `rnaas.tex` rewrite (T4, T5). Which to
  actually submit, and when, is coordinated with the broader publishing sequence (personal todo).

**Dependency order:** T1 (versioning policy) → T2 (JOSS paper refresh) → T3 (cut v1.0.0). T4/T5
(RNAAS notes) are independent of T1–T3 and can proceed in parallel; **`vgpra` (T4) is packaged
first**.

---

## T1 — Codify release-versioning guidance (do first)

Make the next version number a deterministic function of "what changed since the last tag", so any
Claude session can decide it without judgement calls.

**Deliverables**
- `VERSIONING.md` (repo root) — the SemVer policy for this repo (table below).
- `CHANGELOG.md` (repo root, [Keep a Changelog] format) with a `## [Unreleased]` section.
- `scripts/next_version.py` — reads the last git tag + the `Unreleased` changelog section and prints
  the recommended next version **with its reasoning** (the "easy for Claude" mechanism).
- A pointer bullet in `CLAUDE.md` → `VERSIONING.md`, plus a working-rule that every PR appends to
  `CHANGELOG.md`'s `Unreleased`.

**The policy (put in `VERSIONING.md`).** Public API = the `python -m jansky_research.<slice>` CLIs,
each module's `__all__`, the `run()` signatures, and the public functions of
`data.py`/`pipeline.py`/`report.py`.

| Bump | When | Examples |
|------|------|----------|
| **MAJOR** `X.0.0` | Backward-incompatible public-API change | remove/rename a slice module or CLI; change a `run()`/public-function signature incompatibly; drop a Python version; change a results-file schema consumed downstream |
| **MINOR** `1.X.0` | Backward-compatible addition | a new slice/module; a new public function or `--flag`; a new optional extra; additive behavior |
| **PATCH** `1.0.X` | No public-API change | bug fix; paper/doc/test-only change; dependency pin with no API effect; figure/macro regeneration |

**The "next version" recipe (codify in `VERSIONING.md` and `scripts/next_version.py`):**
1. `git describe --tags --abbrev=0` → last released tag (`v1.0.0` after T3).
2. Read `CHANGELOG.md`'s `## [Unreleased]` section (and cross-check `git log <last-tag>..HEAD`).
3. Classify entries as Added / Changed / Fixed / Removed. **Highest severity present wins:** any
   Removed or breaking Changed → MAJOR; else any Added → MINOR; else → PATCH.
4. `scripts/next_version.py` prints e.g. `next: 1.1.0 (MINOR — Unreleased has "Added: lineconf slice", no breaking changes)`.

**Steps**
- [ ] Write `VERSIONING.md` (policy table + recipe + the single-source-of-truth note: version lives
      in `pyproject.toml` and `CITATION.cff`; Zenodo takes it from the release tag).
- [ ] Write `CHANGELOG.md` with a `[Unreleased]` section seeded from `git log 0e2ad04..HEAD` (i.e.
      everything since the last `main` commit, grouped Added/Changed/Fixed) — this becomes the
      `[1.0.0]` section in T3.
- [ ] Add `scripts/next_version.py` (pure stdlib: parse `git tag`/`git log` + the changelog; print
      recommended bump + reason; exit non-zero if `Unreleased` is empty). A unit test on a fixture
      changelog keeps it in the 85% floor.
- [ ] `CLAUDE.md`: add a "Versioning" bullet under **Working rules** pointing at `VERSIONING.md`, and
      "every PR adds a `CHANGELOG.md` `Unreleased` entry".
- **Done when:** `python scripts/next_version.py` prints a correct recommendation on the seeded
      `Unreleased`; `VERSIONING.md` + `CHANGELOG.md` committed; CLAUDE.md points at them.

---

## T2 — Refresh the JOSS paper `joss/paper.md` (from the review)

**Review findings (why it needs updating):** the paper is a June-2026, six-slice snapshot. Concrete
staleness:
1. **Scope:** "six self-contained analysis modules" + the six-item bullet list → now ~40 merged
   slices (README slice table is authoritative). Don't list all 40; **group by domain** with a few
   exemplars each.
2. **"CPU-only" is now false** (title + Summary + Statement of need). The GPU/ROCm slices
   (`torchfdmt`, `torchdsp`, `svsbi`) are a real part of the toolkit. Reframe as **"CPU-first, with
   optional GPU (ROCm/CUDA-portable, pure-PyTorch) acceleration"** — the core install and CI stay
   CPU-only; GPU is opt-in extras. Fix the title accordingly.
3. **"Two of the six results are negative"** → the honest-null/limit body is now large (e.g.
   `vgpra`, `typeii`, `rfitrend`, `pte2`, `frbwait`, `frblens`, `skr`, `wdpulsar`, `svsbi` coverage).
   Recast the methodological-contribution paragraph around the *pattern* (recover-a-known validation +
   honest nulls at scale), not a 2-of-6 count.
4. **"six AASTeX papers under `papers/`"** → current count.
5. **Missing dimensions to add:** neural SBI (`svsbi`), the streaming + static dual-orchestration is
   already there (keep), and the **`jansky-observe` station bridge** (`pull-station-data` → self
   collected HI joins the public-archive slices) as a data-source note.
6. **Date** → release date; **version** → 1.0.0 scope.
7. **`paper.bib`:** add citations for any newly-referenced datasets/tools (PyTorch/ROCm, SBI, CHIME
   Cat 2, SPICE-RACS, FASHI, etc. as cited).

**Steps**
- [ ] Rewrite **Summary**: one-para intro + a **domain-grouped** capability list (FRB/time-domain;
      pulsars; HI & spectral-line; solar/heliospheric; planetary radio; RM/Faraday & cosmology;
      continuum variability; GPU/ML — SBI + PyTorch DSP), each with 1–2 named exemplar slices; point
      to the README table as the full inventory.
- [ ] Fix the **CPU-only** framing in the title and body → CPU-first + optional GPU extras.
- [ ] Rewrite the **Statement of need** methodological paragraph around recover-a-known + honest-null
      at scale; keep the three-audience framing.
- [ ] Update **Functionality**: current counts, the GPU extras, the `arxiv-submit` +
      `pull-station-data` helpers, `make reproduce`.
- [ ] Update `paper.bib`; refresh the `date`.
- [ ] Rebuild/validate the JOSS paper (JOSS renders `paper.md`; check it compiles via the JOSS
      Docker preview or `make` target if present).
- **Done when:** `joss/paper.md` accurately describes the v1.0.0 scope, makes no CPU-only overclaim,
      and lists no stale counts; `science-reviewer` (or a careful read) confirms no overclaiming.

---

## T3 — Cut the `v1.0.0` release (version bump + publishing-todo update)

Depends on T1 (policy) and T2 (paper describes the released scope).

**Steps**
- [ ] Bump version to `1.0.0` in **`pyproject.toml`** (`version = "1.0.0"`) and **`CITATION.cff`**
      (`version: "1.0.0"`).
- [ ] Check **`.zenodo.json`** for a hardcoded version (it currently has none — Zenodo takes it from
      the tag; leave as-is unless a version field appears).
- [ ] Move `CHANGELOG.md`'s `Unreleased` → a dated `## [1.0.0] — 2026-…` section.
- [ ] Update the Obsidian publishing todo (`efforts/radio_astronomy/research_paper_todo.md`): change
      every `v0.0.1` → `v1.0.0` in §0 (the "hold the tag" note) and §1 (the Zenodo release steps and
      title "jansky-research v1.0.0"); note the version decision + rationale.
- [ ] Per the todo §1: wire Zenodo (GitHub switch ON) **first**, then draft the GitHub release
      `v1.0.0`, publish, confirm the Zenodo record + `.zenodo.json` metadata, copy the concept DOI,
      add the DOI badge to `README.md`.
- **Done when:** `pyproject.toml` + `CITATION.cff` say `1.0.0`, `CHANGELOG.md` has a `[1.0.0]`
      section, the todo no longer says `v0.0.1`, and the tag/release/Zenodo checklist is ready to
      execute (the account-bound steps stay Joe's to click).

---

## T4 — `vgpra` RNAAS short-form `papers/vgpra/rnaas.tex` (package first)

A compact showcase of the house method (recover-a-known → honest null). Template =
`papers/frbstats/rnaas.tex` (`\documentclass[RNAAS]{aastex631}`, `\input{generated/macros}`, title +
author, one `\begin{abstract}`, a single `\section{Result}`, one figure, acknowledgements +
`\software{}` + bibliography).

**Source material (already in the repo):** `papers/vgpra/main.tex` (the full paper — title *"Why a
Blind Modern Reanalysis of the Voyager 2 PRA Data Does Not Recover the Uranus and Neptune Radio
Rotation Periods"*), `papers/vgpra/generated/macros.tex`, and the figure `papers/vgpra/figures/vgpra.pdf`.
Key macros to reuse: `\vgSynInjected`, `\vgSynRecovered`, `\vgRealUPeriod`, `\vgRealUSpread`,
`\vgRealURecovers`, `\vgRealNPeriod`, `\vgRealNRecovers`.

**Steps**
- [ ] Write `papers/vgpra/rnaas.tex` (≤1000 words, 1 figure): abstract states the controlled null;
      the Result section makes the recover-a-known→null arc explicit (synthetic recovers the injected
      period to ~1 min; neither real ice-giant period recovered; blind total-power can't do it —
      historical geometric modelling was essential). Every number `\input` from `generated/macros`.
- [ ] Reuse `figures/vgpra.pdf` (the periodogram / recovery figure) as the single figure.
- [ ] Keep the AI-use acknowledgement + `\software{}` block from the frbstats template; author =
      Joseph Barbere only.
- [ ] `make figures`/`make paper` so macros exist; build and check the RNAAS PDF fits the 1-page
      / ≤1000-word limit.
- [ ] Add a Summary-table row and check-box to the publishing todo §3 (already listed as a
      candidate; mark it "packaged" once the `rnaas.tex` builds).
- **Done when:** `papers/vgpra/rnaas.tex` builds a ≤1000-word, 1-figure note with all numbers from
      macros; ready for the §3 submission steps.

---

## T5 — `spectra` RNAAS short-form `papers/spectra/rnaas.tex`

Same pattern; the thematically-matched follow-up ("the apparent USS population is really the TGSS
flux-scale systematic").

**Steps**
- [ ] Confirm `papers/spectra/` has `generated/macros.tex` + a suitable figure (regenerate with
      `make figures` if needed); note its title/result macros.
- [ ] Write `papers/spectra/rnaas.tex` (≤1000 words, 1 figure) from `papers/spectra/main.tex`: the
      cautionary result that raw TGSS×NVSS ultra-steep-spectrum selection is dominated by the TGSS
      flux-scale systematic and the candidates don't survive the de Gasperin cross-check. Numbers from
      macros; AI-use ack + `\software{}`; author Joseph Barbere.
- [ ] Build and check the ≤1000-word / 1-figure limit.
- [ ] Mark it "packaged" in the publishing todo §3.
- **Done when:** `papers/spectra/rnaas.tex` builds within RNAAS limits, numbers from macros.

---

## Verification / gates (whole plan)

- Each new/edited paper: `make figures` first so macros exist; build the PDF; confirm **no
  hand-typed numbers** (every result is `\input`).
- T1: `scripts/next_version.py` passes its unit test and prints a correct recommendation; repo stays
  ruff/mypy/pytest-green with the 85% floor.
- T2: `joss/paper.md` re-read against the README slice table — no stale count, no CPU-only overclaim;
  `science-reviewer` sign-off on framing.
- All code changes go on a branch, not `main`; squash-merge (repo working rules).
- **Sequencing gate:** don't tag `v1.0.0` (T3) until T1 + T2 land, and don't push or submit anything
  outward until Joe green-lights it (other in-progress slices have their own separate gates,
  tracked in the personal todo).

## Open decisions for Joe

- **Version = `1.0.0`** (recommended here) vs a hedged `0.1.0` — plan assumes 1.0.0.
- **RNAAS order:** `vgpra` first (methods showcase), `spectra` second (follow-up) — package both, or
  just `vgpra`? Submission timing is coordinated in the personal publishing todo.
- **JOSS paper:** full rewrite of Summary/Statement-of-need in one pass (T2 as written), or a lighter
  "update counts + fix CPU-only" touch if you'd rather keep it close to the current text.
