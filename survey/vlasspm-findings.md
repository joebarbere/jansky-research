# Findings — blind VLASS proper motions across four epochs (plan 64)

`jansky_research.vlasspm` + `scripts/vlasspm_real.py`. Status 2026-09-26: **pipeline validated
on its recover-a-known; five unexplained candidates, not yet vetted — no discovery claim.**

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
