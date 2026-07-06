# Findings — lensed-repeater delay-pattern search in Cat 2 (plan 42)

`jansky_research.frblens` runs the first catalogue-level search for the Dai & Lu (2017) / Li+
(2018) lensed-repeater signature — a recurring, fixed burst-pair delay — across the CHIME/FRB
Catalog 2 repeaters, and converts the (expected) null into the first empirical upper limit on
the lensed-repeater fraction.

## GATE 0 (2026-07-06; data mirror shared with plan 39 — see `survey/frbwait-findings.md`)

- **Kill-condition read of arXiv:2605.19653** ("IMBH microlensing in Cat 2"): millisecond-scale
  INTRA-burst echo search via dynamic-spectrum autocorrelation, non-repeaters only (2
  candidates, ~10²–10³ M☉ lenses). Different regime entirely; cited as fence. Same for the
  Cat-1 autocorrelation candidate (arXiv:2406.19654) and the baseband work (Leung/Kader 2022).
- **No catalogue-level burst-to-burst delay-pattern search exists anywhere** — the review
  (arXiv:2412.01536) states the idea as a proposal only.

## The design lessons (worth the findings file by themselves)

**Lesson 1 — first contact with the real catalogue falsified our planned null.** The plan (and
first implementation) used the frbwait day-scramble (keep sidereal phase, redraw day) with a
300-s match tolerance — and "detected" exactly the three most intrinsically clustered
repeaters (20220912A M_max=82, 20180916B, 20201124A, all p=0.005): dense activity epochs put
many pairs at near-integer sidereal-day delays within a loose tolerance, and a day-scramble
destroys that intrinsic structure, so beating it only proves clustering — frbwait's result,
not lensing. Fix: **phase permutation, not day scramble** — keep each burst's sidereal DAY,
permute the within-transit-window phases. Per-day burst counts survive exactly; only
sub-window fixed-delay coherence is destroyed. A synthetic clustered control (40-day epoch, 3
bursts/day) is null under this scramble (regression-tested); injected lensed trains are
recovered at the scramble floor. Tolerance tightened 300 s → 5 s (fixed-delay coherence is a
sub-window statement).

**Lesson 2 — GATE-2 caught a frame error: mjd_400 is TOPOCENTRIC.** A lens delay is fixed at
the barycenter; the annual Roemer term drifts a fixed barycentric delay by ~Δ×9 s per day of
delay (±150 s at Δ=26 d) in topocentric pair delays — ~30× the 5-s tolerance. The search would
have been blind to real signals. Fix: barycenter all TOAs (astropy, CHIME location, per-source
positions; offline ephemeris; regression-tested — a barycentre-fixed injection is recovered in
the bary frame and missed in the topo frame). The null re-applies day-keyed offsets (Roemer
varies <0.01 s across a transit window).

**Lesson 3 — GATE-2 caught the DM cut rejecting genuine pairs.** The planned 1 pc/cc
"measurement floor" is far below the per-source fitted-DM scatter (structure-driven, 2–7 pc/cc
for active repeaters; 67–95% of intra-source pairs differ by >1) — genuine image pairs, whose
DMs are fitted independently, would fail it at the same rate. Fix: per-source tolerance
3√2·σ_DM (robust scatter), the SAME cut applied in the injections with image-DM fit scatter
drawn N(0, √2·σ_DM) — the cut now costs injections what it costs real pairs.
(Consistency arithmetic: ΔDM = 1 pc/cc ↔ 26 ms at 400 MHz ≪ 5 s.)

## Result (200 scrambles/source, barycentric; `results/frblens_metrics.json`)

- **0 detections among 33 searched repeaters** (≥5 bursts, span > 2 d); cleanest non-detection
  p=0.81 (all sources fully consistent with their phase-permutation nulls).
- **First empirical lensed-repeater fraction limit: < 0.091 per searched repeater (95% CL)**,
  scoped to the injection-mapped sensitivity region.
- **Transit-survey selection function, made explicit by injections**: both images are jointly
  detectable only when the delay is within the ~15-min transit window of an integer number of
  sidereal days — off-comb injection cells are dark (0 images detected). The deepest train
  (20220912A, 373 bursts) is sensitive over 77.8% of the injection grid (delays 2–30 sidereal
  days + off-comb controls, magnification ratios 0.1–1.0).

## Interpretation limits (GATE-2 material)

- The limit is per-searched-repeater within the sensitivity region — NOT an absolute
  optical-depth measurement; delays off the sidereal comb and magnification ratios below each
  source's empirical fluence floor are unconstrained (stated in the paper).
- The injection map is computed on the deepest train only (20220912A; 77.8% of its grid
  sensitive — the dark cells are the deliberately off-comb delays and the lowest magnification
  ratio); the other 32 searched trains are shallower — stated in the paper.
- Sub-day delays (< 1 transit spacing) are excluded: image pairs and intrinsic same-window
  clustering are indistinguishable there.
- The transit window is treated as 15 min uniformly; the real window width is
  declination-dependent — a second-order effect on the injection map.

## Reproduce

`uv run python -m jansky_research.frblens --n-scramble 200 --out .` (~2 min CPU with the local
mirror). Offline CI leg: `--offline`.
