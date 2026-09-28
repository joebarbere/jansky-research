# Findings — blind VLASS proper motions across four epochs (plan 64)

`jansky_research.vlasspm` + `scripts/vlasspm_real.py`. Status 2026-09-27 (run 5, shape-aware
errors + compactness cut): **UV Ceti recovered blind and the only candidate left; zero new
movers; limit < 9.2e-5 per deg^2 (>= 3 mJy, point sources)** -- see the last two sections. (Earlier status lines below
are kept as the record of runs 1-3.)

## GATE 0 refresh (2026-09-26)

- ADS re-search since 2026-07-01: still OPEN (nearest hits are flux-only VLASS transient work and
  known-object studies; arXiv:2606.27414 "VLASS Data Products" is context, not a scoop).
- **Four epochs, not three:** a QL4.1 catalogue exists (2025-08-26 to 2026-02-08, 1.09M
  components). E1-E4 baseline ~8 yr.
- **Data URLs moved:** the NRAO `vlass-dl` webcache is offline after a disk failure, so every URL
  in `vlass.py` 404s. CIRADA E1/E2 components + subtile tables come from CANFAR
  (`https://www.canfar.net/storage/vault/file/cirada/continuum/vlass_data/`); QL3.x/QL4.1 from the
  temporary public SharePoint `VLASS_catalogs/QL/` folder NRAO links. Per-component epochs:
  E1/E2 via the subtile table's `DATEOBS`; E3/E4 have per-component `MJD`.
- **The plan's recover-a-known was unreachable as written:** UV Ceti is not detected in E1 (subtile
  observed 2018-02-09, rms ~0.14 mJy), so E1xE2 linkage can never find it. It is present in
  E2/E3/E4. Its PM is 3.23"/yr (the plan's 3.4 is BL Cet's).
- Clean counts: E1 1,874,643 / E2 1,873,524 / E3 2,360,730 / E4 1,081,153 (Gordon+2021 flags for
  E1/E2; Memo-22 `Flag==0` for E3/E4).

## Method as built (and one simplification that sized it)

Orphans (no counterpart within 2.5" in the other epoch) -> E1xE2 linkage in the rate window ->
E3 collinearity at each component's own epoch -> scramble null -> injection completeness, over
**every epoch triple** (E1-E2-E3, E1-E2-E4, E1-E3-E4, E2-E3-E4). The plan's "multi-day GPU
all-pairs" framing was wrong: linkage is a fixed-radius KD-tree query (< 20" at 5"/yr), and the
whole real leg runs in ~6 minutes on a CPU. Astrometric floors are **measured** from bright
(>= 10 mJy) static sources, not assumed: 0.10 / 0.09 / 0.10 / 0.11" per epoch (pair scatters
0.13-0.15").

Two design findings from the synthetic fixture:
1. **The static-match radius sets a rate floor** (~2.5"/baseline), and E2->E3 is the shorter
   baseline. The plan's 0.3"/yr lower bound is not reachable; completeness must be reported vs
   rate (below).
2. **E2's positional error enters the E3 residual twice** (fitted rate + anchor):
   residual = -d1*r + d2*(1+r) - d3 with r = dt23/dt12. Treating the terms as independent
   rejected real movers at ~3 sigma. Fixed; residuals are Rayleigh-calibrated (median 1.167 vs
   1.177), locked by a test.

## Run 1 — failed, kept as a record (`data/vlass/work/run1/`, not in results/)

746 candidates against a scramble null of ~0.3, and **UV Ceti not recovered**. Diagnosed rather
than tuned:
- **86% of candidates had a same-epoch neighbour within 30"** vs 25% for a random component.
  Extended sources decomposed differently in each epoch produce spatially *correlated* orphans
  that line up; an arcminute RA-scramble destroys that correlation, so the null was blind to the
  dominant false-positive population by construction (CLAUDE.md: a resampling test is silent on
  what it cannot vary).
- **The flux-consistency cut rejected the target population:** UV Cet is 1.9 / 11.0 / 1.2 mJy
  across E2-E4 (ratio 9.5 > 3). Radio stars flare.

Fix: orphans must be isolated (no other component within 30" in their own epoch); the flux cut
is off by default.

## Run 2 — 2026-09-26 (`results/vlasspm_metrics.json`, `results/vlasspm_candidates.csv`)

| triple | orphans (a, b, c) | pairs | triplets | candidates | null (chance/scramble, 50 reps) |
|---|---|---|---|---|---|
| E1-E2-E3 | 256k, 264k, 253k | 1,336 | 4 | 4 | 0.02 |
| E1-E2-E4 | 256k, 264k, 108k | 1,336 | 1 | 1 | 0.00 |
| E1-E3-E4 | 205k, 384k, 76k | 2,826 | 0 | 0 | 0.00 |
| E2-E3-E4 | 211k, 382k, 86k | 1,819 | 1 | 1 | 0.00 |

- **Recover-a-known PASSES:** UV Ceti found blind via E2-E3-E4 at **3.45"/yr** (Gaia 3.23"/yr;
  the 7% excess is plausibly the unresolved UV/BL Cet pair — BL Cet sits ~2" away). E3
  residual 2.3 sigma.
- **Completeness** (20,000 injections, union over triples): 0.00 at 0.30-0.53"/yr, 0.08 at
  0.53-0.92, **0.92-0.99 above 0.92"/yr**; identical at 1.5 and 3 mJy. The search is a
  >~1"/yr search; that is the honest statement of its range.
- **Five other candidates** against a null expectation of ~0.02: rates 1.2-1.6"/yr, peak
  1-6 mJy; three with no Gaia/CatWISE counterpart, one with a CatWISE source, two whose Gaia
  query failed (network) and must be retried.

### Why the five are NOT discoveries yet

- **Four of five sit at Dec -15 to -39 deg**, where E1 astrometry is documented to degrade
  (~1" south of Dec -20 vs ~0.5" north; Memo 22). The floors above are a single all-sky number
  from bright sources, so southern faint-source errors are likely underestimated, inflating the
  significance of E1 offsets — and the scramble null preserves the error model, so it cannot
  see this either. Three of the five have E3 residuals of 2.2-3.0 sigma, near the cut.
- The 0.02 null therefore bounds chance alignments **under the assumed errors**, not the
  candidates' reality. Same shape as run 1's failure, one level down.

### Next

1. Declination-dependent floors (fit per Dec band from bright statics) and re-run.
2. Retry the two failed Gaia queries; image-level inspection of each survivor (VLASS cutouts via
   the `radio-cutout` skill) — a real mover shows a point source walking across three images.
3. Only then compute a surface-density limit; the one in the metrics file
   (`real_density_limit_per_deg2_3mJy`) treats unvetted candidates as detections and is **not**
   a result.

## Vetting (2026-09-26): zero new movers; the southern clustering is beam elongation

Three steps, in order of cost.

1. **Gaia retry.** The two candidates whose Gaia query had failed have no Gaia DR3 or CatWISE
   source either. So 5 of the 6 non-UV-Cet candidates were "optically dark", which is what an
   astrometric artefact looks like as much as a dark mover.
2. **Declination-banded floors (run 3).** Floors fitted per Dec band from bright statics:
   ~0.17-0.20" at Dec -40..-20 against ~0.08-0.10" further north. The candidates did not go
   away (7, including UV Cet), and **5 of the 6 non-UV-Cet candidates sit south of Dec -15**, a
   band holding ~23% of the area (binomial p ~ 0.003). A systematic, not a population. The
   bright-source floor cannot see what is wrong with these faint, partly resolved sources.
3. **Image inspection** (`scripts/vlasspm_vet_images.py`; CADC SODA cutouts at every epoch;
   `results/vlasspm_vetting.json`, `results/vlasspm_vetting_montage.png`). Two tests per
   epoch: S/N at the predicted track position, and S/N at the *first* detection's position --
   a mover leaves that spot empty.
   - **UV Ceti is textbook:** a point source walks the track, its first position is empty
     afterwards (S/N 2.0, 1.6), and a marginal S/N 4.6 source sits on-track in 2018 -- just below
     the E1 catalogue threshold, which is why E1 has no entry.
   - **All six others are static extended sources.** Emission persists at the first position in
     later epochs (S/N 5-16). c3 and c5 are compact knots fixed on diffuse structure; c0 and c1
     are faint extended blobs. **c2 and c4 (Dec -39) are compact in 2019/2022 but stretched N-S in
     2024; the image headers give beams of 4.8-4.9" x 1.7-1.8" in 2024 against ~3.5" x 2.2"
     otherwise.** At low declination the VLA observes near the horizon, the beam elongates, and
     the elongation differs between epochs; the fitted centroid of a partly resolved source
     moves along it. That, plus extended structure, is the southern systematic -- and the
     scramble null holds each epoch's beam and error model fixed, so it could not see it.

**Result: UV Ceti recovered blind (3.45"/yr vs Gaia 3.23"/yr); zero new movers.** The 95%
upper limit on optically dark (no Gaia DR3 / CatWISE2020) radio movers at 0.92-5"/yr, for
isolated compact sources above ~3 mJy at 3 GHz, is **< 9.2 x 10^-5 per deg^2** over the
33,838 deg^2 of the E1-E2-E3 sky: 2.996 / (area x the 0.961 mean completeness in those rate
bins). That is fewer than ~3.1 such movers in the whole VLASS sky. The
`real_density_limit_per_deg2_3mJy` in `vlasspm_metrics.json` treated unvetted candidates as
detections and is superseded by `vlasspm_vetting.json`.

**Caveats that belong in any write-up:** the completeness comes from injections that assume the
catalogue's positional errors plus the measured floors, and the vetting shows faint partly
resolved sources have larger, beam-dependent position errors -- so completeness for a real
faint mover near other emission is overstated. The isolation cut (no neighbour within 30")
also removes real movers that pass near other sources; the injections, placed 60-120" from
real components, do not pay that cost. A shape-aware (beam-deconvolved) astrometric error
model is the fix before a paper.

## Shape-aware astrometric errors (run 4, 2026-09-27)

The fix the caveats above called for: per-component error ellipses, an empirically calibrated
systematic for resolved sources, and injections that pay the isolation cost. Evidence:
`results/vlasspm_metrics.json` (`real_calibration`, `real_ablation_measurement_only`,
`real_run3_candidates_rescored`, `real_previous_run3`), `results/vlasspm_candidates.csv`
(`vlasspm_candidates_run3.csv` keeps the old list), `results/vlasspm_vetting.json`,
`results/vlasspm_vetting_montage.png`.

### What the catalogues provide (per epoch)

All four tables carry the same PyBDSF shape set per component: fitted `Maj`/`Min`/`PA`,
deconvolved `DC_Maj`/`DC_Min`/`DC_PA`, the restoring beam `BMAJ`/`BMIN`/`BPA` (arcsec, deg),
`Peak_flux` and `Isl_rms`. E1/E2 are the CIRADA CSVs (no `SNR` column; S/N = Peak/Isl_rms),
E3/E4 the NRAO QL FITS (`SNR` equals Peak/Isl_rms). No image headers are needed. Deconvolved
size is exactly 0 for 26% (E1), 25% (E2), 23% (E3), 20% (E4) of clean components.

### The thermal term was already in the catalogue -- and runs 1-3 read it wrong

Condon (1997, PASP 109, 166), eq. 21, with the effective S/N of Condon et al. (1998, AJ 115,
1693), eq. 26: `8 ln2 sigma^2(x0)/theta_M^2 = 2/rho^2` with (aM, am) = (5/2, 1/2) along the fitted
major axis and (1/2, 5/2) along the minor; theta = fitted (convolved) FWHM, theta_N^2 = bmaj*bmin.
Exponents checked against Prandoni et al. (2000, A&AS 146, 41, App. A) and PyBDSF's
`get_errors`, which produced these catalogues. Recomputed from the columns, it reproduces the
catalogue `E_DEC` to 1 part in 10^5 (median) in every epoch, and `E_RA` likewise **only without a
cos(dec) factor**: `E_RA` is already an on-sky angle. The run 1-3 loader multiplied it by
cos(dec), shrinking RA errors by 1.15x at Dec -30 and up to 2.7x at Dec +68. So the "Condon
model" adds no new thermal information; what changed is (a) the RA bug, (b) the ellipse and the
RA-Dec correlation are kept instead of a circularised error, (c) the floor bands.

The run-3 floor bands (-40,-20,0,30,90) left Dec < -40 on the all-sky floor and averaged
0.26" (Dec < -35) with 0.12" (-25..-20). Fine southern bands (-90,-35,-30,-25,-20,0,30,90),
still from >= 10 mJy statics, flatten the statics' Dec dependence (table below).

### The systematic, measured on statics

Statics: components isolated (no other within 30") in both epochs, matched within 5" (wider
than the 2.5" orphan radius, so the tail that *creates* orphans is in the sample), all six
epoch pairs, 4.87 M pairs. Fit on even 1-degree RA strips, validated on odd ones. Model per
detection: `Condon ellipse + floor^2 I + k^2 * (deconvolved shape FWHM^2 covariance) + q^2 *
|beam - mean beam|`. Exploration showed the excess scatter is along the *source's* major axis
(0.40" vs 0.21" across it for 4-8" sources, S/N > 15) and grows with deconvolved size far more
than with beam change -- hence the shape covariance, not an isotropic term. Fit target: the 3-sigma
tail (P(z>3) = 0.0111 for Rayleigh), per (size x S/N) and (size x beam-change) bin, because the
search's decisions are 3-sigma cuts.

**Result: k_struct = 0.10 (a resolved source's centroid wanders by 10% of its deconvolved FWHM
along each axis), q_beam = 0.0** -- once size is modelled, the beam-change term does not reduce
the loss (a random-subsample fit during exploration gave 0.02-0.03; it is poorly constrained).

Validation half (2.43 M pairs), fraction beyond 3 sigma (Rayleigh 0.0111):

| error model | all | size 0-1" | 4-6" | 6-10" | Dec < -35 | -35..-30 | P(z>5) all |
|---|---|---|---|---|---|---|---|
| runs 1-3 (cos-dec bug, circular, old bands) | 0.041 | 0.025 | 0.082 | 0.126 | 0.146 | 0.069 | 0.0084 |
| catalogue ellipse + fine floors | 0.032 | 0.020 | 0.063 | 0.101 | 0.058 | 0.039 | 0.0069 |
| + k_struct = 0.10 | **0.0103** | 0.018 | 0.011 | 0.019 | 0.025 | 0.015 | 0.0012 |

Where the calibration **fails** (reported, not tuned away):
- **The core and the far tail cannot both be matched by a Gaussian.** With the 3-sigma tail right,
  the median z is 0.81 (Rayleigh 1.18): errors are ~30% too large for the typical static. And
  P(z>5) is 0.0012 against 3.7e-6 -- a heavy tail ~300x Gaussian in every bin. The model is
  conservative in the core and still optimistic beyond ~4 sigma.
- **Dec < -35: 2.3x** too many 3-sigma outliers (0.025); -35..-30: 1.4x. The southern
  systematic is not fully captured by size + floor.
- **Compact sources (deconvolved size < 1"): 1.6x** (0.018); these have the smallest errors, so
  the floor's own non-Gaussian tail dominates.
- **Bright, large sources (S/N 15-50, size 6-10"): 3.8x** (0.042).
- Beam change > 0.8 (fractional): 1.2x (0.014) even with q = 0 preferred overall.

### Run 4: search, null, candidates

| | run 3 | measurement-only (ablation) | run 4 (k = 0.10) |
|---|---|---|---|
| E1-E2 pairs | 1,336 | 1,334 | 1,142 |
| candidates (all triples) | 7 | 7 | **11** |
| scramble null (E1-E2-E3, per scramble) | 0.02 | -- | 0.02 |

Larger errors cut both ways: they make a static source's apparent shift less significant
(fewer pairs), but they also widen the E3 tolerance for resolved sources, which admits more
static-source triplets. The net was *more* candidates, not fewer. The null did not move, so
the excess over it is still the correlated static population the scramble cannot see.

**The run-3 false positives became less significant, as they should** (E1->E2 shift, sigma,
measurement-only -> run 4): c0 10.4 -> 9.1, c1 4.35 -> **1.99**, c2 9.1 -> 5.1, c3 11.5 -> 4.9,
c4 6.5 -> 3.4, c5 15.7 -> **2.59**. c1 and c5 now fall below the 3-sigma linkage cut; c0, c2, c3,
c4 still pass. Their E3 residuals shrink too (e.g. c2 2.32 -> 1.14), which is the loosening
described above -- the E3 test does not reject them under either model.

**UV Ceti: still recovered** (E2-E3-E4, 3.45"/yr), E3 residual **1.79 sigma** (run 3: 2.2).

**Six new candidates** (J170.9129+76.8236, J198.1384-39.1432, J202.1059-35.0560,
J38.2031-30.7428, J275.6343+64.9963, J102.1962-37.9361), every one with a deconvolved size of
7-15" in at least one epoch, four of six south of Dec -30. **Image vetting: all six are static
extended sources** -- emission persists at the first-detection position in the later epochs
(S/N 4.3-47); J198.1384-39.1432 is the c2/c4 pattern again (2024 beam 4.89" x 1.77", source
stretched N-S). J202.1059-35.0560 also has Gaia and CatWISE sources within 5".
**Result unchanged: UV Ceti, and zero new movers among 11 candidates.**

### Completeness and the limit

Injections now (i) draw their positional scatter from the same model -- the Condon ellipse the
catalogue would quote at the local beam and noise for the assumed flux and size, plus the local
floor, plus k^2 x shape -- and carry that covariance into the search; (ii) half are placed at a
random sky position (uniform over a 15' disc round a host) instead of 60-120" from one.

**The isolation cost turned out to be ~1%, and the old injections were already paying most of
it.** Only 1.3% of random sky positions have a component within 30" in E1: the 25-36%
non-isolated fraction of *real components* reflects multi-component sources, and a mover is not
one. The 60-120" placement only kept injections away from their own host, not from other
sources, so its neighbour fraction was the same 1.3%. Realistic and isolated classes agree to
within the binomial noise.

Completeness per rate bin (0.30, 0.53, 0.92, 1.62, 2.85, 5.0"/yr), 3 mJy, 20,000 injections:

| | 0.30-0.53 | 0.53-0.92 | 0.92-1.62 | 1.62-2.85 | 2.85-5.0 |
|---|---|---|---|---|---|
| run 3 (isolated, old errors) | 0.000 | 0.092 | 0.915 | 0.986 | 0.983 |
| run 4, realistic placement, point source | 0.000 | 0.065 | 0.917 | 0.986 | 0.981 |
| run 4, realistic, 2" deconvolved size | 0.000 | 0.143 | 0.890 | 0.985 | 0.985 |

The 2" variant exists because a faint point source does not come out of PyBDSF with size 0
(UV Ceti: 1.2", 0, 2.9"); size 0 gives an injection the smallest structure term the model
allows. It raises the 0.53-0.92 bin (larger tolerance) and lowers 0.92-1.62 (weaker linkage
significance). 1.5 mJy is indistinguishable from 3 mJy, as in run 3.

**Limit (95%, 0 new movers, 0.92-5"/yr, E1-E2-E3 area 33,838.5 deg^2): < 9.3 x 10^-5 per deg^2**
(the weaker of point: 9.21e-5, mean completeness 0.961; 2": 9.29e-5, 0.953), i.e. < 3.8 across
the whole sky. Run 3 quoted 9.2e-5 -- **the limit is effectively unchanged** (1% weaker). The
completeness above ~1"/yr is insensitive to the error model because fast movers are displaced
many sigma under any of them.

### What this does and does not fix

- Fixed: the RA-error bug; circularised errors; the southern floor bands; resolved-source
  optimism in the 3-sigma tail (4x -> 1x overall); injections that did not share the error
  model; the isolation cost is now measured rather than assumed.
- Not fixed: the heavy >4-sigma tail (300x Gaussian); the Dec < -35 excess (2.3x); the E3 test's
  loss of power for resolved sources, which is why candidates went *up*. Every candidate in runs
  3 and 4 other than UV Ceti has been a resolved, often southern, static source. A **compactness
  requirement** (deconvolved size consistent with zero given `E_DC_Maj`, calibrated on
  injections with the measured faint-source size noise) is the natural next cut, and would
  also let the structure term be dropped for the survivors -- not done here.
- The scramble null remains blind to the static-source false positives by construction; the
  vetting, not the null, is what bounds them.

## Compactness cut (run 5, 2026-09-27)

Every run 3-4 candidate except UV Ceti was a resolved static source, so run 5 asks the obvious
question: does a mover have to look like a point source? Evidence: `real_run5_compactness` in
`results/vlasspm_metrics.json`, `run5_compactness` in `results/vlasspm_vetting.json`.

**Metric.** Per detection, `c = DC_Maj / BMAJ`: the catalogue's deconvolved major axis over the
restoring beam's major axis (0 where PyBDSF could not deconvolve). Threshold and cross-epoch rule
are parameters (`search(..., compact_max, compact_rule)`), applied after the E3 test, so the cut
only removes.

**Pre-stated criterion** (in `scripts/vlasspm_real.py` before the cut touched any real
candidate): threshold = the *smallest* value on the grid 0.3-3.0 (step 0.1) keeping >= 95% of
recovered injected point movers overall **and** in every S/N bin with >= 200 recoveries; rule =
"all three detections compact" vs "at least two of three", whichever rejects the larger fraction
of isolated static sources seen in E1, E2 and E3 (the population false candidates come from, not
the candidates), ties within 0.02 going to "all". ("Loosest threshold that keeps >= 95%" read as
"loosest *requirement*": the most aggressive cut that still meets the keep target -- any larger
threshold trivially keeps more.)

**Realistic size noise.** Faint point sources do not come out of PyBDSF at size 0. Reference
sample: GCNS stars (Gaia EDR3, < 100 pc; VizieR J/A+A/649/A6), proper-motion-propagated to each
component's epoch, matched within 1.0" -- **303 detections** (E1 73, E2 83, E3 102, E4 45); the
same match with stars shifted 2' gives **16** chance matches (5%, likely background AGN, which
if anything widens the size distribution and loosens the threshold). Their `c` has
95th percentiles of 1.60 (S/N < 7), 1.15 (7-10), 1.09 (10-15), 0.83 (15-25), 0.54 (25-50),
0.32 (> 50); 33-39% are exactly 0 in every bin. Injections draw (DC_Maj, DC_Min)/BMAJ from the
same epoch and S/N bin (epochs pooled below 30 detections -- which, with 303 in total, is most
bins: the epoch dependence is not resolved by this sample).

**Calibration** (40,000 injections at 1, 1.5, 3, 6 mJy, random-sky placement, 23,932 recovered).
Kept fraction of recovered point movers:

| threshold | all: overall | all: S/N < 7 | 2of3: overall | 2of3: S/N < 7 |
|---|---|---|---|---|
| 0.7 | 0.724 | 0.403 | 0.921 | 0.785 |
| 0.9 | 0.890 | 0.714 | **0.984** | **0.951** |
| 1.0 | 0.915 | 0.785 | 0.990 | 0.971 |
| 1.2 | 0.960 | 0.881 | 0.997 | 0.990 |
| 1.6 | **0.984** | **0.953** | 0.999 | 0.998 |
| 2.0 | 0.997 | 0.987 | 1.000 | 1.000 |

Criterion thresholds: **all -> 1.6; 2of3 -> 0.9.** Rejection of 953,850 E1-E2-E3 statics:
all@1.6 6.6%, 2of3@0.9 **14.1%** (of statics larger than 2 beams in some epoch: 100% vs 89%).
**Chosen: 2 of 3 detections with c <= 0.9.** One noisy faint detection is tolerated, which is
what lets the threshold be tight; "all" has to stay loose to survive the single worst of three
noisy sizes.

Per rate bin (2of3 @ 0.9): 0.931 at 0.53-0.92"/yr, 0.983 / 0.987 / 0.989 above 0.92"/yr (the
0.30-0.53 bin has 10 recoveries). **The 95% target is missed in the 0.53-0.92 bin (93%)** -- the
criterion was stated per S/N bin, not per rate bin, so this is reported rather than re-tuned.

**UV Ceti survives.** Its c values are 0.41 (E2), 0 (E3), 1.46 (E4) -- 1.2", 0 and 2.9" against
beams of 2.84", 2.14" and 2.0". The E4 detection alone fails 0.9; two of three pass. (Under
"all" it would also have passed at 1.6.) It is a real point mover with a 2.9"-deconvolved
detection, which is exactly the size noise the calibration was built for.

**Real search.** Run 4: 11 candidates -> run 5: **1 (UV Ceti)**. No new candidates (the cut only
removes). All 10 image-vetted static sources are removed. Their c triples, for the record:
several sit close to the threshold -- J198.1384-39.1432 (0.92, 0, 1.77) would survive at 1.0, and
run-3 c0 = J314.1346-32.6145 (1.01, 0, 1.12) at 1.1 -- so the clean sweep depends on the grid step
the criterion landed on, and the images, not the cut, are what established these as static.

**Null.** Scramble chance candidates with the cut: 0.02 / 0 / 0 / 0 per scramble (unchanged).

**Completeness and limit** (20,000 realistic-placement point-source injections with size noise;
mean over the 0.92-5"/yr bins; E1-E2-E3 area 33,838.5 deg^2; 0 new movers):

| | no cut | with cut |
|---|---|---|
| 3 mJy | 0.9648 -> 9.18e-5 / deg^2 | 0.9645 -> **9.18e-5** |
| 1.5 mJy | 0.9610 -> 9.21e-5 | 0.9520 -> **9.30e-5** |

The cut costs ~0 at 3 mJy (S/N ~ 20, where stars are compact) and 1% at 1.5 mJy; the limit is
**< 9.2e-5 per deg^2 (>= 3 mJy)** and < 9.3e-5 at 1.5 mJy. Run 4 quoted 9.3e-5, using the weaker
of a size-0 and a size-2" injection; the size-noise draw supersedes both.

**What the cut excludes, physically.** The limit now applies to movers unresolved at 2.5":
- moving objects with resolved radio emission -- a pulsar-wind-nebula bow shock or trail around a
  high-velocity pulsar, a star with extended wind/jet emission;
- radio-emitting binaries whose separation is a sizeable fraction of the beam (~1-3"): blended,
  they fit as one elongated Gaussian. UV/BL Cet (~2") passed, but only because the cut tolerates
  one resolved detection; a pair emitting comparably in two epochs would fail;
- a mover blended with unrelated background emission in two of its three epochs (also removed by
  the isolation cut in most cases);
- faint movers: the keep fraction is 95% by construction at S/N < 7 and 93% at 0.53-0.92"/yr.

**Unresolved.** The reference sample is small (303) and pools epochs; the size noise is
calibrated on stars, and a population with different S/N or position-in-tile distribution could
differ. The 0.53-0.92"/yr bin misses the 95% target.

## Referee round 1 on the RNAAS note (2026-09-27): MINOR revision

The limit is verified: 2.996 / (33,838.5 deg^2 x 0.9645) = 9.18e-5, matching the JSON; every
macro checked, the figure matches, and all 10 DOIs match Crossref. There is also a property the
paper does not claim: with the cut, no by-eye step lies between the data and "0 new movers".
Claim strength and definitions need work:

1. **(major) The UV Ceti rate is a two-epoch number with parallax ignored.** The triplet rate is
   the E2->E3 pair rate; E2->E4 gives 3.364, and 3.41 once parallax (374 mas) is removed.
   Parallax alone moves two-epoch rates by 0.14-0.23"/yr, the size of the whole "7%
   discrepancy". UV Cet and BL Cet are 2.26" apart (confirmed, Gaia DR3 2016.0); both have RUWE
   10-12 and orbital motion, and the system value is 3.344. The radio positions drift toward BL
   Cet's prediction, so the component cannot be assigned. **Fix:** a three-epoch,
   parallax-aware rate with an uncertainty, compared with the system motion; call it a detection
   of the unresolved Luyten 726-8 system.
2. **(major) "Optically dark" is positional coincidence within 5".** A 5" circle contains
   >= 1 Gaia or CatWISE source about 42% of the time; 5 of the 10 static candidates had one. A
   truly dark mover would be called "not dark" about 40% of the time. It cost nothing here (the
   counterpart step removed no candidate), but the definition must say "no co-moving
   counterpart", or fold the coincidence rate in (limit ~1.5e-4).
3. Rate domain and weighting: the completeness is an unweighted mean over log bins (a log-uniform
   rate prior), and the worst bin gives 9.6e-5. The figure's line joins imply ~0.55 at the
   0.92"/yr edge; use steps and state the weighting.
4. The cut's threshold was set on injections, but the idea came after the vetted candidates, and
   one commit holds everything. Say so, give the 0.02 margin (J198.1 survives at 1.0, J314.1 at
   1.1), and claim the real strength: no by-eye step in the limit's chain.
5. Error tails overstate completeness by at most ~1% (3 mJy) or ~2% (1.5 mJy); one sentence.
6. Straight-line injections ignore parallax. State the distance below which completeness does
   not apply (~5 pc; UV Cet's parallax adds ~1.2" to its E3 residual).
7. The scramble null underpredicted false positives about 500-fold; say what it can and cannot
   see.
8. E4 covers only half the sky (~17,300 deg^2), so "four epochs" needs qualifying.
9. Prior work: Atri et al. 2022 (MNRAS 517, 5810; VLBA proper motions of compact variable radio
   sources without optical counterparts; targeted, not blind). Narrow "so far" to survey-scale
   searches and cite it. ADS finds no earlier blind survey-scale search.
10-12 (nits):
   - flux domain: steady 2.4-3.75 mJy per epoch, or state >= 1.5 mJy;
   - "whole sky" assumes isotropy;
   - explain "17 of 25";
   - the size-noise sample is thin at S/N < 7;
   - the 0.53-0.92 sentence could be cut.

## Response to referee round 1 (2026-09-27)

New evidence: `results/vlasspm_referee1.json` (from `scripts/vlasspm_referee1.py --out .`), the
re-vetted counterparts in `results/vlasspm_metrics.json` / `vlasspm_candidates.csv`, and
`vpmReal*` / `vpmSyn*` macros for every number below that the note quotes.

1. **UV Ceti rate (major).** `vlasspm.fit_track` fits a straight line plus parallax to the three
   committed E2/E3/E4 detections, each weighted by its full covariance (catalogue ellipse + floor
   + structure term); parallax factors from astropy's Earth ephemeris (`parallax_factors`, tested
   against the geometric displacement; a synthetic 8-epoch track with a 0.25" parallax is
   recovered). Results: **3.413 +/- 0.042"/yr with the Gaia DR3 parallax (374 mas) fixed**
   (chi^2 2.9, 2 dof); 3.363 +/- 0.043 without parallax; free parallax 192 +/- 113 mas (1 dof),
   3.389. Gaia DR3 (queried by TAP): UV Cet 3.232 (RUWE 10.5), BL Cet 3.429 (RUWE 12.4), 2.3"
   apart; SIMBAD's system motion (GJ 65 / G 272-61, UCAC4) 3.344. The fixed-parallax rate is
   1.6 sigma from the system and 0.4 sigma from BL Cet's Gaia value, 4.3 sigma from UV Cet's.
   The note now calls it a detection of the unresolved Luyten 726-8 system and lists parallax,
   orbital motion and blending as the effects that can move a three-epoch rate at this level;
   the component is not assigned. Parallax alone moves UV Cet's E3 point by 1.2".
2. **"Optically dark" (major).** Measured at 500 random footprint positions (seed 20260927, a
   random E1 component displaced 2-5'): P(>= 1 source within 5") = **0.28** (Gaia 0.14,
   CatWISE 0.24), lower than the review's ~42% estimate; within 1.5" 0.032, within 1" 0.010.
   Adopted definition, now in `vet_counterparts` and `vlasspm.comoving`: **dark = no co-moving
   counterpart** (a Gaia/CatWISE source within 3 sigma + 0.5" of the three-epoch track's
   prediction at the catalogue epoch, or with proper motion within 30% of the radio rate). The
   run-4 list re-vetted: 6 of 11 dark under the co-moving rule vs 5 under the 5" rule (one
   candidate had only an off-track source inside 5"); the four static candidates with a
   CatWISE source *at* the track position are most likely the hosts of those static radio
   sources. After the compactness cut only UV Ceti reaches this step (three co-moving Gaia sources), so the step
   removed no candidate and chance coincidence cost nothing.
3. **Rate weighting and the edge.** Stated: the average is an unweighted mean over the log bins
   (log-uniform prior). Worst coarse bin: 0.921 -> 9.6e-5. Fine bins (3 mJy, cut, 20,000
   injections over 0.8-5"/yr): 0.26 (0.80-0.92), **0.65 (0.92-1.0)**, 0.87 (1.0-1.1), >= 0.97
   above 1.1. A mover just above 0.92"/yr is constrained only to < 1.4e-4. The figure is now
   steps per bin with the fine bins as points.
4. **Provenance of the cut.** The note now says the idea followed the vetted candidates, the
   threshold and rule were fixed on injections, and gives the margins (J198.1-39.1 passes at
   0.92, J314.1-32.6 at 1.01, against 0.9), then claims the property the review identified: no
   by-eye step lies between the catalogues and the limit.
5. **Tails.** `tail_bound` takes the largest excess of the compact statics' 3-sigma tail over
   Rayleigh among the S/N bins a detection at that flux can occupy (per-epoch peak 0.8-1.25 x
   flux over the 10th-90th percentile rms): **0.6% at 3 mJy, 1.9% at 1.5 mJy**. An upper bound,
   quoted as such.
6. **Parallax domain.** `parallax_floor`: over 2,000 E1-E2-E3 positions and dates, the median
   E3 residual per arcsec of parallax is 2.35 and the median E3 standard error for a 3 mJy point
   source 0.26"; parallax exceeds one standard error inside **8 pc** (the review estimated ~5).
7. **Scramble null.** One sentence on what it measures (chance alignments of unrelated orphans)
   and cannot (a static source whose centroid shifts).
8. **Four epochs.** E4 area from `real_per_triple`: 17,301 deg^2.
9. **Prior work.** Atri et al. 2022 cited (Crossref-verified); the sentence is narrowed to
   survey-scale searches.
10-12. Flux domain stated per epoch (2.4-3.75 mJy for the 3 mJy injections, 1.2-1.9 mJy for
   1.5 mJy; from `INJ_FLUX_SPREAD`); "whole sky" now says "if movers are isotropic"; the
   synthetic 17 of 25 is explained from the offline run (all 8 missed movers are slower than
   1.1"/yr, one of them also not isolated; `vpmSynNMissed`, `vpmSynMissedMuMax`); the
   0.53-0.92 sentence is cut (that bin is below the limit's rate range).

Side fix: `vlasspm_real.run5` now writes the metrics file as one payload, so it no longer
carries a `_merge` block listing the whole file as "retained from a previous run".

Length: 111-word abstract, ~1,270 words of text in the two sections (pdftotext, including the
figure caption), over the ~1,000 target and under the RNAAS 1,500 limit.

## Referee round 2 on the RNAAS note (2026-09-27): MINOR revision

All 12 round-1 findings are addressed. All 86 macros resolve with no `--`; the ~60 traced match
the JSON; the figure matches; the three new DOIs match Crossref. RNAAS guideline (fetched): at
most 1,500 words and a single figure or table, so the note is within the limit.

1. **(major) "1.6 sigma from the system" compares magnitudes only.** The vectors are radio
   (3.391, 0.385) vs SIMBAD/UCAC4 system (3.296, 0.564): |d| = 0.20"/yr, chi2 = 13 for 2 dof.
   The orbital scale, UV - BL in Gaia, is 0.21"/yr. Parallax cannot explain it: it is fixed in the
   fit, and the UV/BL parallax difference is about 0.014". Drop SysSigma; give the vector offset
   next to the orbital scale, and name orbital motion and blending. Fixing the parallax is right;
   +/-0.04 is a fair statistical error (0.05 scaled); component switching would show at >= 6 sigma.
   The findings' scalar sigmas vs BL and UV should become vector values (chi2 6.2 and 35.5).
2. **(major, framing) Lead with 1.1-5"/yr.** Every fine bin there is >= 0.967 complete, giving a
   worst-bin limit of 9.15e-5, the same 9.2e-5 with no rate prior. The averaged, worst-coarse and
   edge limits can leave the text. At 1.5 mJy there are no fine bins (worst coarse bin 0.896 ->
   9.9e-5): run fine bins at 1.5 mJy, or quote it as averaged.
3. The co-moving rule text should say "0.5" floor added in quadrature" (radius >= 1.5"). The 0.03
   chance rate is a lower bound (0.048 at 2", 0.10 at 3"). The PM clause never fired. The
   positional clause flags 4 static candidates as "co-moving". Not applying it to injections is
   correct.
4. The 8 pc domain is computed correctly. Say what it means physically (v_t > 36 km/s at
   >= 0.92"/yr) and put "beyond ~8 pc" in the abstract. Check "fast pulsars" as an example (the
   rate domain far exceeds typical pulsar PMs).
5. The tail bound is an estimate, not "at most": pair offsets vs a three-position residual gives
   ~1.5x, so about 1% and 2-3%.
6. "No by-eye step" depends on the 2-of-3 rule winning: under "all at 1.6", J314.1 would need
   images. Legitimate (pre-stated, decided on statics); fold it into one clause.
7. Add "persistent (three-epoch)" to the abstract's domain (no flux-ratio cut; 0.8-1.25x spread).
8-10 (nits):
   - units on "system motion 3.34";
   - "at least 0.97" vs 0.967;
   - "rises across 0.92-1.1" wording;
   - Atri+2022 selected compact flat-spectrum variable Galactic-plane sources (no counterpart was
     a result, not the selection).

**Cuts proposed** (~-230 words, to ~950-1,000):
- the synthetic-field sentence;
- the chance-coincidence numbers;
- the worst-bin and edge limits;
- the fine-bin number chain;
- "16 expected by chance" and the "all at 1.6" numbers;
- "3.36 without parallax" and the component values;
- the "reproduces catalogue errors" and 5 sigma clauses;
- the scramble 0.02 sentence;
- the without-cut completeness pairs;
- the long Reproducibility data list.

## Response to referee round 2 (2026-09-28)

New evidence: `results/vlasspm_referee2.json` (`scripts/vlasspm_referee2.py --out .`), reaching
the note through new `vpmReal*` macros; `\vpmRealSysSigma` is removed.

**Run 1: 1.5 mJy fine bins** (cut, 20,000 injections over 0.8-5"/yr): 0.331, 0.640, 0.830,
**0.937** (1.1-1.25), 0.972, 0.979, 0.976. Worst bin over 1.1-5"/yr 0.937, so the limit for
1.2-1.9 mJy is **9.4e-5 per deg^2** (no parallax).

**Run 2: parallax-aware injections.** `_inject(..., distance_pc=)` displaces every detection by
(1/d) x `parallax_factors` at its own date (tested: the displacement equals 0.5 x the factors at
2 pc; a field whose epochs repeat the season shows no loss, one half a year out of phase does).
3 mJy, cut, fine bins, same seed; worst fine bin over 1.1-5"/yr:

| distance | none | 16 pc | 8 pc | 4 pc |
|---|---|---|---|---|
| completeness | 0.967 | 0.963 | 0.949 | 0.888 |
| relative loss | -- | 0.4% | 1.9% | 8.2% |

Adopted domain: the smallest tested distance with a loss below 5%, **8 pc** (the 4-8 pc crossing
is not resolved by this grid; the round-1 one-sigma proxy also gave 8 pc). Physically,
v_t > 4.74 x 1.1 x 8 = **42 km/s**. The headline limit uses the completeness AT 8 pc (0.949):
**< 9.3e-5 per deg^2** for persistent compact sources of 2.4-3.75 mJy per epoch at 1.1-5"/yr
beyond 8 pc, with no rate prior. Without parallax the same range gives 9.15e-5; the averaged,
worst-coarse and edge limits of round 1 stay in `vlasspm_referee1.json` and leave the text.

**F1, vectors.** Radio (3.391, 0.385) "/yr (fixed Gaia parallax) against:
- the SIMBAD/UCAC4 system (3.296, 0.564): |d| = **0.20**"/yr, chi^2 = **13.0** (2 dof);
- BL Cet (Gaia DR3): |d| = 0.16, chi^2 = **6.2**;
- UV Cet (Gaia DR3): |d| = 0.29, chi^2 = **35.5**.

These replace round 1's scalar 1.6 / 0.4 / 4.3 sigma, which compared magnitudes only. The
orbital scale |UV - BL| in Gaia is **0.21**"/yr, and the components' parallaxes differ by 6 mas,
so parallax cannot supply the offset. The note names orbital motion and blending.

**F2-F9.**
- F2, F7: the abstract leads with the 1.1-5"/yr no-prior limit for persistent (three-epoch)
  compact sources of 2.4-3.75 mJy beyond 8 pc.
- F3: the co-moving rule says the 0.5" floor is added in quadrature. The chance-coincidence
  numbers are dropped ("only UV Ceti reaches this step, so chance coincidence cannot have
  removed a mover"), and so is the word "co-moving" for the positional clause.
- F4: the distance domain is expressed as a tangential velocity. "Fast pulsars" is dropped, since
  the fastest known pulsars move at ~0.1-0.4"/yr, below the range; "nearby cool brown dwarfs,
  invisible to Gaia" is kept.
- F5: the tail effect is now "an estimated 0.6% / 1.9%". This is still the pair-offset excess;
  the review's ~1.5x amplification for a three-position residual would give ~1% and ~3%. That
  factor is not computed by the pipeline, so it is not quoted.
- F6: the 2-of-3 dependence is folded into one clause.
- F8: units on the system motion; 0.967; "climbs steeply between 0.92 and 1.1".
- F9: Atri+2022's targets described as compact, flat-spectrum, variable Galactic-plane sources.

**Cuts.** Every item in the review's list is applied, and the south-of-Dec -35 sentence on the
vetted candidates is also dropped. Length (pdftotext): abstract **117** words; Search + Results
+ Reproducibility **~1,120** words including the 57-word figure caption; the two sections without
Reproducibility are ~1,080.
