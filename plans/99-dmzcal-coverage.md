# Plan 99 — `dmzcal`: do per-burst DM→z credible intervals cover the hosts' true redshifts?

Status: **plan, controls frozen 2026-10-10 before any posterior was computed on real data.** The
localized-FRB tables have been downloaded (GATE-0 data check) and their `z` column exists on
disk, but no p(z|DM) has been evaluated for any real burst. Changes after the first real run go
in `survey/dmzcal-findings.md` with the reason (CLAUDE.md, slice pattern).

## Context

FRB 20240304B (Caleb et al. 2026, *Science*, doi:10.1126/science.adz2675; arXiv:2508.01648)
has a JWST/NIRSpec host redshift z = 2.148 ± 0.001 (Table 2), twice the previous record for a
localized FRB. Its DM-only prediction was z_Macquart = 2.8 (+0.6, −1.2), and its 95% DM range
is [1.628, 3.397] ("Clues" section). The truth sits inside that range but well below the
centre. Per the "Foreground contributions to the DM" section, the burst sits about 550 pc cm⁻³
above the median DM_cosmic(z), and about 235 of that is Virgo (Supplementary S1.4).

Most FRBs have no host. Their redshifts, energies and any population inference built on them
come from p(z|DM). The question here is whether the published per-burst intervals are
**calibrated**: do 68% intervals contain the true z 68% of the time on bursts whose z is known?

## GATE-0 (run 2026-10-10)

**Novelty.** Building a new estimator is **not** novel.
- `zdm` (github.com/FRBs/zdm; James et al. 2022, MNRAS 509, 4775) already computes p(z|DM) per
  survey, including Macquart scatter (`logF`) and a log-normal DM_host. The relevant script is
  `zdm/scripts/Plotting_zDM/pz_given_dm.py`.
- fruitbat 2.x gives a p(z) from EAGLE sightline scatter. It has no host-DM distribution and was
  last pushed 2022-01-11 (read in its source: `_frb.py:334`, `methods.py:185`).
- `FRBs/FRB` `frb.dm.igm.z_from_DM` gives a point estimate with no scatter.
- Cordes, Ocker & Chatterjee 2022 (arXiv:2108.01172; abstract only) report the RMS error of a
  DM-only estimator on 14–32 localized FRBs. That is a point-estimate accuracy check, not an
  interval-coverage test.

No paper or repo test was found that measures **empirical coverage / PIT calibration** of
per-burst p(z|DM) on the localized sample. That is the open gap. The search was one pass plus
code greps of `zdm`, `FRB` and fruitbat, not an exhaustive sweep, so the novelty check is
re-run before the note is drafted.

**Data.**
- `github.com/FRBs/FRB`, `frb/data/FRBs/FRB*.json` (BSD-3, last commit 2026-05-06): 187
  files, 111 with non-null `z`, all with `DM`, `DMISM` (NE2001), `ra`, `dec`, `refs`. Max z < 1.2.
- `frb/data/Galaxies/public_hosts.csv` has a `Spectrum` column (70 of 94 hosts flagged), used
  for the spectroscopic-z rule below.
- 10 bursts appear under two names (e.g. `FRB20240210` / `FRB20240210A`) and must be deduplicated.
- FRB 20240304B is **not** in the repo. It is added by hand from arXiv:2508.01648:
  - DM_obs 2458.20 ± 0.01 (scattering-corrected, Table 1)
  - l, b = 269.8676°, 72.1035° (Table 1)
  - DM_ISM 28.1 (NE2001) / 20.8 (YMW16) (Table 1)
  - z 2.148 (Table 2)
- Dead or blocked sources: FRBSTATS (domain parked, catalogue built from a private sheet), TNS
  (403 without a bot-ID), and VizieR (Gordon+2023 not hosted).
- `pygedm` does not build (`f2c.h` missing), so the per-burst NE2001 `DMISM` from the JSONs is
  used.
- `zdm` installs from git into a throwaway Python 3.12 environment (114 packages incl. lalsuite).
  It is **not** added to the repo's dependencies. The real-data leg runs it in a separate env
  and commits its outputs (the ROCm-venv precedent).

**Citation correction.** `survey/candidate-gaps.md` cited arXiv:2511.01195 for F ≈ 0.32. That
paper (Sang & Lin) is about the DM_host distribution: log-normal μ = 5.03, σ = 0.96 on 117 FRBs
(abstract only). F = 0.32 is `zdm`'s default (`parameters.py`, `IGMParams.logF`). Its source
paper must be found and given a locator before it is quoted. Baptista et al. 2024
(arXiv:2305.07022, abstract) gives only log10 F > −0.89 (3σ).

## Estimators under audit (frozen)

- **E1 — `zdm` p(z|DM), per survey.** Default parameter state of the pinned `zdm` commit
  (recorded in the results JSON). Each burst is evaluated under the survey file that matches its
  discovery instrument: CRAFT ICS / CRACO, DSA, CHIME, MeerTRAP, Parkes, FAST. Bursts with no
  matching survey model are excluded and counted.
- **E2 — minimal Macquart + log-normal host posterior** (new, in `src/jansky_research/dmzcal.py`):
  - the Macquart et al. 2020 p(DM_cosmic|z) form (Nature 581, 391, eq. 4 — locator to be
    confirmed);
  - the same F and host parameters as E1;
  - a uniform-in-comoving-volume prior in z;
  - no survey selection.

  E2 exists so the audit can say how much of E1's calibration comes from the survey selection
  function. It is a comparison arm, not a proposed replacement for `zdm`.
- **E3 — fruitbat 2.0.1 `Batten2021`** (added 2026-10-10, before any real posterior). Settings
  are the package defaults: `subtract_host=False`, `prior="uniform"`, and `dm_galaxy` = the
  per-burst `DMISM`. PIT is computed from `calc_redshift_pdf`. This is the widely-installed
  simulation-based (EAGLE) alternative. Its scatter was **not** fitted to localized FRBs, so
  every burst is on the certification side for E3. It is still scored on the same
  certification sample as E1 so the estimators are compared like for like, and the full-sample
  E3 number is reported alongside.

  **Environment.** fruitbat no longer runs on a current stack:
  - it fails to import `astropy.cosmology.core` on astropy ≥ 6;
  - it needs `np.long` / `np.asscalar`, so numpy < 1.23;
  - it imports `pygedm` at module load, and `pygedm` does not build here.

  The recipe that works is Python 3.9, numpy 1.22.4, astropy 4.3.1, scipy < 1.10, h5py < 3.8,
  setuptools < 70, `pip install --no-deps fruitbat`, and a `pygedm` stub that **raises** if
  called (so it cannot silently supply a Milky Way DM). The recipe is recorded in the driver
  script. The decay itself is a finding for the note, stated as fact and not as criticism.
- **E4 — `FRBs/FRB` `frb.dm.igm.z_from_DM`, point estimate only** (`corr_nuisance=True`,
  the default). It returns no interval, so it **cannot enter the coverage test**. It enters only
  the secondary point-accuracy table below.
- Milky Way: per-burst `DMISM` (NE2001), plus a halo term of 50 pc cm⁻³ that is folded into the
  host distribution's floor. No other foreground (Virgo, intervening groups) is modelled for any
  burst. That is what a DM-only user gets.

## Provenance (*an oracle that fed a fit may never certify it*)

`zdm`'s default parameters were fitted on localized FRBs. Any burst used in that fit is on the
**production side** and may not certify E1.

**Step 0 (before any posterior):**
1. Identify the paper the default parameter state comes from, with a locator.
2. List the localized bursts used in that fit, and also those used for the 2511.01195 host
   prior if E2 adopts it.
3. Commit the list to `results/dmzcal_provenance.json`.

The **certification sample** is the deduplicated, spectroscopic-z localized bursts **not** on
that list. The full sample's numbers are reported alongside, labelled in-sample. If the
certification sample has fewer than 25 bursts, the outcome is UNDERPOWERED by rule, and the
in-sample numbers are still reported, labelled as such.

## Sample rules (frozen)

1. **Dedupe.** Group JSON files on the UTC date plus the burst name with any trailing letter
   removed, and keep the file with more non-null fields. Every merge is listed in the results JSON.
2. **Spectroscopic z only.** Keep a burst if `public_hosts.csv` flags a spectrum, or if its
   reference names a spectroscopic redshift. Unflagged hosts form a separate sensitivity arm
   and are never pooled silently.
3. **Repeaters** stay in, flagged, and their coverage is also reported separately. Some repeater
   hosts are known to have large DM_host, so this split is pre-stated, not post hoc.
4. **FRB 20240304B** is excluded from every coverage statistic, because it was chosen *after*
   its z was known. It is reported as a single recover-a-known case:
   - its PIT value under E1 and E2;
   - whether its 68% and 95% intervals contain 2.148;
   - E1 compared with the paper's own [1.628, 3.397].

## Statistics and controls (frozen)

For each burst b and estimator E: **PIT_b** = ∫₀^{z_true} p_E(z|DM_b) dz.
- *Coverage at q* is the fraction of bursts whose central q-interval contains z_true, for
  q = 0.68 and 0.95.
- *Bias* is the median of PIT − 0.5.

**C0 — null false-fail rate, computed before freezing the thresholds** (the hiblend lesson).
- Draw M = 2000 synthetic samples of the certification-sample size N from E2's own generative
  model, so the model is exactly calibrated by construction. Each draw takes z from the prior,
  then DM from p(DM|z), then computes PIT.
- The **calibration rule** is:
  - |coverage_q − q| within the central 95% band of its null distribution, for both q;
  - a KS test of PIT against U(0,1) with p > 0.01.
- **Gate:** the rule's measured false-fail rate on the null must be ≤ 0.07. If it exceeds 0.07,
  the band is widened to the null's 97.5% quantiles and the change is logged *before* the real
  run.
- The final thresholds and the measured false-fail rate are written to
  `results/dmzcal_controls.json` and committed **before** step 3.

**C1 — power against a planted misspecification.**
- Generate data with DM_host drawn from a distribution with twice the mean, or with
  F × 2 (two separate cases), and score with the unmodified model.
- At the certification N, the rule must flag mis-calibration in ≥ 80% of draws for at least the
  host case.
- If not, the real outcome is reported as "the sample cannot distinguish calibrated from
  mis-calibrated at this N". Coverage numbers are still quoted, but no calibration verdict is
  given.

**C2 — planted truth through the real pipeline.**
- Replace each certification burst's z_true with a z drawn from that burst's own E2 posterior,
  then run the full E2 pipeline (dedupe, MW subtraction, posterior, PIT) on the result.
- PIT must pass the calibration rule. This tests the code path, not the physics.

**Pre-stated selection caveat (not a control).** Bursts with spectroscopic hosts are biased
toward low z at fixed DM, because faint high-z hosts drop out. That pushes PIT values low
independently of any mis-calibration of p(z|DM). Neither estimator models host
identifiability. A low-PIT bias therefore **cannot by itself** be read as an IGM or host
mis-specification, and the note must say so. If the `P(O|x)` column permits, a post-hoc
stratification by host magnitude is run and labelled post hoc.

## Outcomes (frozen wording)

- **CALIBRATED** — the rule passes on the certification sample and C1 had power.
- **OVERCONFIDENT** — coverage falls below the band at either q.
- **UNDERCONFIDENT** — coverage lies above the band.
- **BIASED** — coverage is in band but the KS test fails with |median PIT − 0.5| > 0.1. The
  direction is reported, and the selection caveat is attached.
- **UNDERPOWERED** — N < 25, or C1 failed.

Each outcome is reported per estimator. No outcome is upgraded by the in-sample arm.

**Secondary point-accuracy table.** This is descriptive only, with no verdict.
- For E1–E4 it reports the median and the normalized median absolute deviation (NMAD) of
  Δz / (1 + z_true). For E1–E3 the point estimate is the posterior median; E4 has only its point
  value.
- This puts the interval-free tool (E4) on the same footing as the others. Point accuracy is
  Cordes et al. 2022 territory and is not the claim.

**Already seen.** The install checks on 2026-10-10 evaluated FRB 20240304B under E3 (median 2.40,
68% interval [2.11, 2.62]) and under E4 (2.52; 2.47 with the nuisance correction against
DM − DM_ISM − 40). That burst is already excluded from every statistic. No other burst has been
evaluated.

## Steps

0. Branch `slice/dmzcal`. Provenance list (above), with the source paper of the `zdm` default
   parameter state and the source of F = 0.32 located.
1. `dmzcal.py`:
   - sample assembly: dedupe, spectroscopic-z rule, 20240304B row;
   - E2 posterior, PIT, coverage, the calibration rule.

   All of it runs on a synthetic offline fixture, with ≥ 85% coverage, and ruff and mypy clean.
2. C0 and C1 on the synthetic generator. Commit `results/dmzcal_controls.json`.
3. Real run (`# pragma: no cover`). E1 runs in the throwaway `zdm` env through a small driver
   script whose output JSON is committed. E2 runs in-repo. C2 runs here as well.
4. `survey/dmzcal-findings.md`, then GATE-2 science review, then the RNAAS note
   `papers/dmzcal/`, with `\software{}` citing `jansky-research`, `jansky` and `zdm`.
   Macros are namespaced `dmzSyn*` / `dmzReal*`, and `preserve_live_macros` is called from
   `_write_macros`.
5. Upstream: `jansky.transients.macquart_redshift` is a linear 900·z. On a quick check it agrees
   with the Planck18 integral to 2% at z = 2.148 (1933 vs 1901 pc cm⁻³, f_d = 0.84), so it
   needs no change. The docstring could point to the full p(z|DM) for anything quantitative;
   that is a course-repo PR, optional.

## Risks

- **Small N after provenance.** Many of the 111 bursts are CRAFT, which is likely `zdm`'s
  fitting sample. The DSA-110 (Sharma24, Law2023) and CHIME Outrigger (Leung+2025) bursts are
  the probable certification sample, and step 0 decides that. If N < 25, the note becomes "an
  audit protocol, and the sample size it needs", which is still a reportable result.
- **Survey-model mismatch.** The CHIME Outrigger bursts may not match `zdm`'s CHIME survey file,
  which models the main CHIME/FRB pipeline. Exclusions are counted, not hidden.
- **Scoop.** `zdm`'s authors are the natural people to do this. Move at slice speed and re-run
  the novelty check before drafting.
