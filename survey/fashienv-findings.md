# Findings — the environment-split FASHI HI mass function (plan 45, DR1 first leg)

`jansky_research.fashienv` splits the FASHI HI mass function (HIMF) by large-scale environment
for the first time — void/wall (Douglass+2023 VoidFinder) and group/field (Tempel+2017) —
recovering the ALFALFA void-HIMF suppression from FAST data.

## GATE 0 (2026-07-08)

- **DR2 is NOT public yet.** arXiv:2606.31539 (156,411 sources) is out, but the catalogue links
  on zcp521.github.io/fashi 404 ("available upon publication"; release ~Aug 2026). Verified the
  404 directly. → **DR1 first leg**, the same DR1-while-DR2-embargoed pattern `rmstructure` used;
  DR2 swap is a one-line source change.
- **FASHI DR1 is on VizieR: `J/other/SCPMA/67.19511/table2`, 41,741 sources** (the ".19511" is
  the SCPMA article number, NOT the row count — the agent's 41,741 was right; my first pass
  mis-read it as 19,511 and I corrected it). Columns: RA/Dec, cz, z, W50, Ssum (flux), Dist,
  logMass — everything the HIMF needs, precomputed.
- **Novelty confirmed**: the DR2 paper and DR1 paper both publish only a GLOBAL HIMF; the sole
  environment study (arXiv:2510.22902) used 230 group galaxies. No DR2 environment-split HIMF /
  deficiency / void study exists.
- **Cross-match catalogues resolve** (all VizieR): Tempel+2017 groups `J/A+A/602/A100`
  (table2 = groups with R200/M200), Douglass+2023 voids `J/ApJS/265/7` (table1 = VoidFinder
  spheres, Planck2018 cosmology, Mpc/h coords — verified the frame matches standard RA/Dec).
- **Plan citation corrected**: Lim+2017 is MNRAS 470, 2982 (not ApJ 854, 62); not used in the
  end (Tempel groups suffice).

## Scope corrections vs plan 45 (both forced by what FASHI lacks)

- **"Gas fraction at fixed M*" DROPPED**: FASHI carries no stellar masses. The literature void
  gas-fraction excess is also weak/dwarf-only (Kreckel+2012), so void-vs-wall HIMF is the
  cleaner statement — exactly what Moorman+2014 / Jones+2018 did.
- **"HI-deficiency vs clustercentric radius" DROPPED**: classical deficiency needs optical
  diameters/types (absent), and the raw median-HI-vs-R/R200 of DETECTED sources is
  selection-biased (stripped galaxies drop out of a flux-limited sample). Replaced with the
  cleaner group-member vs field HIMF split.

## Recover-a-known (offline, in CI)

- Injected two Schechter HIMFs (void: logM*=9.70, α=−1.45; wall: 9.95, −1.25) into a
  flux-limited mock; the 1/Vmax + Schechter fit recovers both knees/slopes within 0.25 dex and
  the correct knee-offset sign. Single-Schechter round-trip also passes.

## Result (FASHI DR1, SDSS-cap overlap)

| environment | logM* | α | N |
|---|---|---|---|
| global (all DR1) | 9.94 | −1.73 | 41,741 |
| void | 9.70 ± 0.08 | −1.6 | 2,538 |
| wall | 9.95 ± 0.05 | −1.7 | (cap) |
| group member | 10.08 ± 0.05 | −1.74 | 6,119 |
| field | 9.89 ± 0.05 | −1.68 | (cap) |

(void/wall use the Planck2018-matched void geometry; group/field use Tempel R200 membership)

- **Void knee offset = −0.256 ± 0.087 dex = 2.93σ** (void suppressed vs wall) — consistent in
  sign, somewhat larger, than Moorman+2014's ALFALFA void HIMF (−0.14 dex, ~2σ). **The headline:
  an independent FAST-based measurement.** Void faint-end marginally flatter (Moorman-like), not
  steeper (Jones-like); within errors at DR1 size — report the sign only.
- **Group knee offset = +0.19 ± 0.073 dex = 2.64σ** (members *higher* than field) — statistically
  the STRONGER offset, but NOT the headline because it's survivor-biased: a flux-limited HI
  survey detects the group galaxies that RETAINED HI (the gas-rich survivors, group knee 10.08 >
  global 9.94); stripped members drop out. So the group knee is an UPPER ENVELOPE and the true
  stripping effect is more negative. The void–wall comparison avoids this (voids don't strip).
- **Distance-model sensitivity (GATE-2 catch):** the void offset depends on the comoving-distance
  frame. Using the EdS (q0=+0.5) relation dilutes it to ~1σ; using the Douglass catalogue's own
  Planck2018 cosmology (q0=−0.527, so the frames match) gives the −0.256 dex / 2.9σ above. The
  matched-cosmology geometry is the correct choice and is what we report.

## Honest caveats (GATE-2 material)

- **Absolute α (~−1.73) is steeper than the published FASHI global (~−1.3)**: the simple 1/Vmax
  uses a single flux limit, not FASHI's full completeness function. The RELATIVE knee offsets
  (same estimator both bins) are robust; the absolute slope is not, and we do NOT claim to
  reproduce the global HIMF.
- **Footprint**: void/group cross-match is the SDSS-cap subset (~27k classifiable of 41,741),
  not the full survey — FASHI's southern/high-z sky has no SDSS optical catalogue.
- **Single void-finder**: VoidFinder (Douglass table1 spheres) only; the V²/VIDE/REVOLVER
  galaxy-membership tables (keyed by NSAID) need an NSA cross-match — a stated follow-on for the
  "report per void-finder" robustness the plan wants.

## Reproduce

`uv run python -m jansky_research.fashienv --out .` (VizieR fetches, ~min). Offline CI leg:
`--offline`. DR2 swap: point `fetch_fashi_dr1` at the DR2 table when it publishes (~Aug 2026).

## Referee round (2026-08-12) — the variance was fine; the biases are not

**A test that could have failed, and didn't.** The referee expected a void jackknife to kill
the 2.93σ offset, on the grounds that the quoted error is Poisson counting noise within one
realisation and cannot see void-to-void variance. Measured by deleting each of the 186
occupied voids in turn and refitting: **jackknife error 0.039 dex**, *smaller* than the fit
error of 0.087, so the offset strengthens to 6.5σ on that axis. A bootstrap over voids agrees
(0.031). Sample variance across voids is **not** what limits this measurement. Both numbers
are now committed; the expectation was wrong and the paper says so.

**What does limit it is bias, which no resampling can see.** Two, both acting in the direction
of the reported signal, and both now stated where the reader meets them:

- *1/V_max assumes uniform density within the accessible volume* — exactly false for a sample
  selected *by* density. For the void bin the accessible volume is the union of void spheres,
  whose fraction falls to zero at the near and far edges of the box while staying finite for
  the wall bin; since V_max is mass-dependent, so is the distortion, and it suppresses φ at the
  high-mass end of the void bin. The offset is now framed as an upper bound on any true
  suppression, with SWML named as the standard fix and explicitly not attempted.
- *The "wall" bin is a Cartesian bounding box of a wedge-shaped survey.* Now quantified in the
  evidence file, which never carried it: **n_wall = 24,175 — 58% of all of DR1** — with a knee
  of 9.954 against the all-sky 9.944, i.e. **0.010 dex apart, a fifth of its own error**. The
  comparison is closer to *void versus everything* than to void versus wall.

**Numbers that were hard-typed and wrong.** The paper quoted the void knee error as ±0.08
against a committed 0.071, and the wall as ±0.05; combined those give 0.094, which does not
reproduce the quoted 0.087 (√(0.071²+0.05²) = 0.0868 does). Both are macros now. The claim
that "every number above is pipeline-generated" is narrowed to the data-derived ones, since
literature values, the flux limit and catalogue sizes are hard-typed.

**Independence overstated.** FAST against Arecibo is an independent HI measurement, but the
*environment* half of both comes from the same SDSS DR7 parent and largely the same voids —
a deeper HI sample over much the same volume, not an independent realisation of large-scale
structure.

Open and reported, not fixed: the single void-finder headline (plan 45 says never to quote
one); beam confusion, which plan 45 required and which reproduces the sign of both offsets;
the injection validation, which places galaxies uniformly in comoving volume and so switches
off the one systematic that matters; the α–M* covariance, discarded by keeping only
`sqrt(diag(pcov))`; `curve_fit` without `absolute_sigma=True`; `n_bins` silently dropped from
every metrics dict by an `isinstance(v, float)` filter; the EdS distance-frame variant quoted
as "~1σ" with no committed number; and `FASHI_FLUX_LIMIT = 0.30`, which sets every V_max and
appears nowhere in the paper.
