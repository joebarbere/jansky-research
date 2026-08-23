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

**Crossref settles a citation in one call, with no search budget.**
`curl -s https://api.crossref.org/works/<doi>` returns the authoritative author list, article
number, volume and title from the DOI already in the `.bib` entry. It caught `lawrance2024`
carrying an author who was never on the paper and a page number that was an article id.
Prefer it to any search — this is the tractable form of *a search summary is not a source*.

**A read-only reviewer can still mutate the repo, by running the code.** The science reviewer
made no edits and nonetheless left `results/typeii_metrics.json` gutted, because *checking
determinism* meant re-running an offline mode that writes to the repo root. Check
`git status` after any agent that executes slice code, and treat "no file edits" as a claim
about intent, not about effect.

**Overclaiming can live in a verb.** `rmstructure` said the pipeline "recovers" an injected
enhancement; the honest verb was "responds to". Nothing numeric had to change for the sentence
to become true — and no test can see the difference. When a validation sentence and its number
disagree in strength, the sentence is usually the thing that is wrong.

**Do not manufacture an expectation to make a measurement look predicted.** When my
closed-form estimate of `rmstructure`'s statistic gave 12.8 against a measured 3.15, the
tempting move was to keep hunting for a model that reproduced 3.15 and present it as "matches
expectation". The honest output is the empirical ensemble plus an explicit statement of what
is *not* modelled.

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

**A GPU benchmark is not a "real vs synthetic" property.** `torchdsp` used one `device` field
to label both the science leg and the timing run, so recording a real GPU benchmark meant
mislabelling CPU-run science as GPU. Separate them (`benchmark_device`, `benchmark_hardware`).
The ROCm venv is `~/.venvs/rocm-test` (Python 3.14, outside the repo's `<3.13` pin) — run with
`PYTHONPATH=src:../jansky/src ~/.venvs/rocm-test/bin/python`, not `uv run`.

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

**A robustness check that cannot vary what it claims to test is worse than none.**
`dr20radio` reported its north/south contrast as "robust" because it survived a conservative
5 mJy RACS limit. That check was vacuous by construction: the common luminosity limit is the
*RACS* one in **both** legs, so raising the flux floor rescales north and south identically and
the ratio cannot move (1.4391 -> 1.4434). The parameter the contrast actually depends on is the
assumed spectral index in the cross-frequency K-correction, which is not measured for the
sample: sweeping alpha from 0 to -1 moves the gap 0.23 -> 1.66 pp, and at alpha = 0 it nearly
vanishes. **Ask what a robustness check is free to change.** If the perturbed quantity enters
both arms of a comparison through the same term, it tests precision, not the claim — and it
lends the claim a credibility it never earned.

**A fitted parameter sitting on its bound is not a converged fit.** `innerrc`'s sensitivity
scan reported `n_converged: 8` and a range 0.19--0.32 GeV/cm^3 because `curve_fit` did not
raise. Two of the eight variants have `v_bulge` at **exactly** 800.0000 km/s — the upper bound
in `decompose_rc` — and the variant supplying the quoted maximum (`\irRhoMax` = 0.3223, the
number carrying "fully compatible with the consensus density") is one of them. SciPy reports
success for a solution glued to a wall. **Record per-variant whether any parameter is at or
near a bound, and exclude those from any quoted range**; a bound is a modelling choice, so a
result that rests on one is reporting the choice, not the data.

**When a systematic acts on a large ensemble, the number to report is a bias, not a variance.**
`dr20radio` applies one spectral index to ~10^5 quasars. Asking "how much does the answer move
if I draw each source's index from the measured distribution?" and quoting the *spread across
realizations* gives 0.016 pp — because a fraction over 10^5 objects averages the draw away as
1/sqrt(N), so any breadth of distribution looks harmless. The scatter's real effect is a
**shift**: the detection criterion is a threshold, thresholds are not linear in alpha, and the
steep and flat halves do not cancel. Measured, it moves each fraction by ~0.3 pp (7-10%
relative) while moving their *difference* by only 0.03. Both facts had to go in the paper —
the contrast survives the scatter, the absolute fractions do not. This is the same shape as
the `rmstructure` bootstrap-SE error: an uncertainty estimated by a procedure that cannot see
the effect it is meant to bound. Before quoting a spread, ask what would make it small
regardless of the truth.

**A .bib entry's coordinates outlive its title and DOI — adjudicate, don't search.** 19
citation defects were unfixable by fuzzy title search but all fell to one procedure: pull the
recorded DOI's full Crossref metadata and compare it field-by-field against the entry's own
author/journal/volume/pages. Five entries had the RIGHT DOI under a wrong title (a repair pass
that trusted titles would have "corrected" them into error); the wrong DOIs were mostly the
*adjacent* identifier (e085 for e084, 169090 for 169089, pages 120 for 112) — constructed by
pattern instead of looked up. **Never construct an identifier; look it up, and when repairing
one, decide which field is broken before changing any.** Also: `export.arxiv.org` 301s HTTP to
HTTPS — an API that returns empty may just need `-L`.

**Run the test even when you expect it to fail the claim — it may not.** The referee
expected a void jackknife to destroy `fashienv`'s 2.9-sigma HIMF knee offset, and I expected
it too (the quoted error was Poisson-within-one-realisation, the `rmstructure` shape exactly).
Measured: jackknife 0.039 dex against a 0.087 fit error, so the offset *strengthened*. The
value of running it was not the verdict but the relocation — it proved the problem is **bias,
not variance**, and no amount of resampling would ever have found the two biases that do limit
that measurement. **A resampling test bounds variance; it is silent on bias, so passing one
is not a clean bill of health.**

**A null result divides by sensitivity, not by sample size.** `frblens` searched 33 FRB
repeaters, found none lensed, and quoted f < 2.996/33 = 0.091. That is right only if every
source was fully sensitive. Measured, the mean injection efficiency was **0.24** and four
sources had eps = **0** -- too few bursts for any injected signal to beat the null, so they
constrained nothing while inflating the denominator. The honest limit is 2.996/sum(eps) =
0.368, **four times weaker**. The paper had also argued that measuring efficiency on only the
deepest source made the limit "conservative", which is backwards: assuming eps = 1 where it is
smaller tightens the limit. **Before quoting "we searched N and saw none", ask how many of the
N the search could actually have seen into.**

**A "forced" measurement that searches is a noise maximum.** `stokesv_discovery`'s method
section said "forced peak flux at the propagated pixel"; the code took the brightest Stokes-I
pixel within 12" and read V there. The census exposed it without any new data: **I > 0 for 54
of 54** quiescent targets, p = 2^-54 for a genuine fixed-pixel measurement. On blank sky a
peak search returns the largest of several independent beams, so it is positive essentially
always and biased high, and the companion quantity is read wherever that maximum fell. For an
upper-limit census at known positions this invalidates the limits. **Test a photometry routine
on pure noise: forced should go negative about half the time.**

**When two quantities fall by the same fraction, suspect the calibration.** `stokesv_discovery`
reported a 10-sigma inter-epoch decline in Stokes V. Stokes I fell 26% and V 23% over the same
pair, leaving V/I constant to 3.4% -- which is what a flux-scale difference between two
independently calibrated observations does, and is not what a change of emission state need
do. The quoted significance carried image noise only; a 5% per-epoch scale term takes it to
3.5 sigma. **Before quoting a significance on a difference of two epochs, ask what else moved
by the same factor.**

**A test can lock a defect in.** `test_run_offline_writes_artifacts` asserted
`\svbSynNTargets` must *not* exist — i.e. it required the un-namespaced behaviour that let an
offline rebuild write the synthetic parent size (400 stars) into the macro the abstract used
for the real census (38). The test passed for months *because* the bug was present. When a
test encodes "X must not exist", ask whether X is the fix.

**To tell a measurement from a prior, move the prior.** `svsbi` reported a break luminosity of
14.5 from a posterior whose mass piled against a box closed at 15. Widening the box to 16.5
moved the median to 15.5 — a 0.96 dex shift against a 0.02–0.14 seed scatter. The number was
reporting the wall. A posterior/prior width ratio does not catch this (the ratio was 0.41,
i.e. "informative"); only re-inference under a different box does. Applies to any bounded
parameter: **quote a bounded fit only after showing the bound does not set it.**

**A rule stated in this file is not a rule the repo follows — audit it.** Three
CLAUDE.md requirements were checked mechanically on 2026-08-12 and each was met by a small
minority of slices: `preserve_live_macros` ("call it from every `_write_macros`") was called
by **5 of 42**; the `\software{}` block citing `jansky-research` was present in **4 of 46**
papers; 41 `refs.bib` files lacked the `janskyresearch` entry entirely. The
`preserve_live_macros` gap was live, not cosmetic: `make figures` runs all 33 offline slices
with `--out .` in the repo root, so one invocation would have blanked every real macro in 33
papers. **Periodically grep for each stated invariant rather than assuming new code inherited
it** — the cost is one command per rule.

**Accumulated hedging is its own distortion.** After four review rounds `dr20radio` spent
46% of its length on the one derived quantity it concludes is a survey artefact, while the two
clean census results shared 22 lines. Every qualification was individually true and the whole
was misleading: length signals importance, so a heavily-caveated section reads as the paper's
headline. When a result survives review by being qualified rather than by being confirmed,
**shorten it to the conclusion and move the workings to the committed JSON.** Cutting 499 to
361 lines made the paper more honest, not less complete.

**A difference is not scale-free; a ratio is.** `dr20radio` compared two survey legs above a
common limit and reported the contrast in percentage points. Raising *either* survey's flux
limit deepens *both* cuts, shrinking both fractions and therefore the pp gap, while leaving
their ratio alone. So a pp gap that moves under a deeper cut is measuring normalisation, not
contrast — and I read a 28% reduction as "the contrast is materially definitional" when the
ratio had moved 2%. **Before interpreting a difference as sensitive to a parameter, check
whether the ratio moves too.** Note this is the same error the 5 mJy variant made, arrived at
from the opposite direction: I built the "mirror" check, and it landed on the same axis
(RACS 5 mJy gives cuts (1.994, 5.000); VLASS 2 mJy gives (2.000, 5.014)).

**Look for the algebraic identity before building machinery around a quantity.** Four rounds
of sweeps, censored estimators and systematic budgets went into `dr20radio`'s
luminosity-matched comparison before anyone wrote down that the K-correction cancels between
a source and its own survey's limit — making the whole thing a redshift-independent flux cut
with alpha entering at exactly one place. That one line answers the double-counting question,
explains why a per-source index moves both fractions but not their ratio, and collapses the
systematic budget to a single number. **Simplify the estimator algebraically first; the
machinery you then need is much smaller.**

**Fixing censoring on one side can bias you as hard as the side you fixed.** Kaplan-Meier
over RACS-detected quasars correctly keeps the steep sources a completeness cut discards —
and silently drops the 2,179 detected by VLASS but not RACS, which are the flat ones. Same
mechanism, opposite sign. When data are censored from both sides no median is identified;
bound it (`censored_median_bounds`) rather than picking whichever one-sided estimate is
convenient. And KM is unbiased under *independent* censoring: here the censoring limit is a
function of the same flux the estimand depends on, so "unbiased" was the wrong word.

**When you kill a robustness check, run its mirror before claiming the axis is closed.**
`dr20radio`'s 5 mJy RACS variant was shown vacuous in one round; the VLASS-side variant --- the
one that *can* move the ratio, because VLASS quotes a per-epoch reliability threshold where
RACS quotes a 95% completeness limit --- had still never been run a round later. It moves the
gap 1.27 -> 0.91 pp. Removing a bad check leaves a hole exactly where a real check belongs.

**A measured systematic is only an improvement if its own systematics are measured too.**
Replacing an assumed alpha with a measured one felt like closing the issue, and the quoted
+/-0.015 bootstrap SE made it look settled. Three checks that could each have failed --- vary
the assumed completeness floor, use one epoch instead of the max, bin by flux --- moved the
median by 0.06, 0.08 and 0.34. The statistical error was the smallest term by an order of
magnitude. **Before quoting a measurement's precision, list what you assumed to obtain it and
vary each one.** If the checks are chosen so they cannot fail, the precision is decoration.

**Write the number after you run it, not before.** Having decided to add the VLASS-limit
sweep, I wrote "most of the contrast is attributable to the limit asymmetry" into the draft
and then ran it: the answer was 28%, and the gap survived. The sentence was in the file
before any evidence existed for it. Order the work so the prose cannot precede the number.

**A flux-limited sample cannot measure a spectral-index distribution without a completeness
cut.** Requiring detection in *both* of two surveys truncates asymmetrically — at the
shallower survey's limit a steep source has already dropped below the deeper survey's limit at
the higher frequency — so the joint-detection median comes out too flat. Measured: -0.62 for
all 5571 joint detections against -0.72 for the 4190 above the flux where the truncation
cannot operate for any alpha >= -1.5. The bias is 0.1 in alpha, the same size as the effect
being measured. `alpha_complete_limit_mjy` computes the threshold; quote the complete sample.

**A committed results file can omit the numbers its own headline is computed from.**
`innerrc` filtered its sensitivity variants down to a curated key list before writing them,
dropping `v_halo` and `h_halo` — the only two parameters `rho_dm_gev` is derived from. The
quoted maximum was therefore unauditable from the evidence, and the railed halo bound behind
three variants was invisible. Commit the whole fit; a curated subset of a fit is not evidence.

**Do not assume a derived quantity is monotonic in its inputs.** Propagating `innerrc`'s halo
uncertainties by pairing `(v+dv, h-dh)` gave a 1-sigma maximum of 0.24; the true maximum is
0.31, because rho_DM *rises* with the NFW scale radius at R0. Scan all corners (or sample)
rather than reasoning about which one is extremal — a wrong corner understates an interval
without ever looking wrong.

**The results JSON needed the same merge rule the macros got, and a guard is only as good
as its marker.** `preserve_live_macros` fixed the cross-run clobber for `generated/macros.tex`
in 2026-08; the results JSON kept the hole for another two weeks and bit three more slices
(`typeii`'s gutted census, `southern`'s near-miss, `torchfdmt`'s hand-patched benchmark row).
`report.preserve_live_results` now applies the same invariant to structured data: a run may only
add information, and a synthetic run never overwrites real evidence. Two things it taught:
**a forced `make figures` is the only honest test** (a plain `make figures` says "nothing to be
done" and proves nothing), and **the first version of the guard protected 3 files out of 25**
because it keyed on `is_real`, which most slices never set, instead of the `source` field that
`guard_real_results.py` already used. Match the marker the rest of the repo uses, then re-run the
destructive test and count. A mixed source (`"synthetic recover-a-known + real RACS-mid epoch
pair"`) must count as real, or the file is unprotected by the very string that documents it.
Note the guard still cannot see 76 results files that carry no marker at all; `make guard-real`
now prints that number rather than skipping them silently.

**A `make` variable that is not the one the target reads fails silently and slowly.**
`make paper SLICE=dr20radio` looked like a single-paper build and rebuilt all 44 (the variable
is `SLICES`). Make does not warn about unused command-line variables. `SLICE` now narrows
`SLICES` and errors on an unknown name.

## The slice pattern (how every result is built)

tested helper (pure NumPy/SciPy/astropy + synthetic offline fixture) → real-data run (network,
`# pragma: no cover`) → GATE-2 science review → AASTeX paper (`papers/<slice>/`) → arXiv package
(`make arxiv`). Plan in `plans/NN-*.md`; findings in `survey/<slice>-findings.md`.

Every new slice paper / RNAAS note `\software{}`-cites **`jansky-research`** (the toolkit, Zenodo
concept DOI `10.5281/zenodo.21482378`) alongside `jansky` — copy the `@misc{janskyresearch}`
`refs.bib` entry from `papers/vgpra/` or `papers/spectra/`. See the `research-publish` skill.

**Merged slices:** forty-four — the authoritative per-slice table
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

- **Review campaign (2026-08):** every paper now passes the mechanical triage
  (`uv run python scripts/triage_papers.py` — zero HIGH/MED findings; run it before any
  submission and after any paper edit). Seven papers have had the full presenter/referee
  round (`.claude/agents/paper-{presenter,referee}.md`, both read-only): atlas3i, dr20radio
  (5 rounds), innerrc, svsbi, stokesv_discovery, frblens, fashienv — each verdict and every
  retraction is in its `survey/<slice>-findings.md`, and the README papers table tracks the
  Reviewed column. Deep-review the remaining papers before submitting any of them; the triage
  finding history says the referee pays for itself.

## Layout

`src/jansky_research/` (slice modules + `data.py`/`pipeline.py`/`report.py`) · `tests/` · `plans/` ·
`fable-ideas.md` (current plan-ready idea list) · `survey/` (committed findings) ·
`papers/<slice>/` · `station/` (build guides for the physical rooftop
station — self-collected data, WIP; owner's working notes live in an Obsidian vault, not here) ·
`workflow/Snakefile` (static-slice file-DAG, drives `make figures`) · `airflow/` (streaming e-Callisto
ingest, Podman DAG) · `.claude/` (agents + skills) ·
`Makefile` (`setup`/`test`/`cov`/`lint`/`typecheck`/`figures`/`paper`/`arxiv`/`reproduce`).
