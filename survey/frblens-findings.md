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

## Referee round (2026-08-12) — the limit was four times too tight

**The denominator was the source count, not the sensitivity.** A null bounds the lensed
fraction at −ln(1−C)/N only if every searched source was fully sensitive. The efficiency had
been measured for exactly one train (FRB20220912A, the deepest of 33), and the paper argued
the other 32 being shallower made the limit "conservative in depth" — which is backwards.
Assuming ε = 1 for sources where it is smaller makes the denominator too *large* and the limit
too *tight*.

Measured properly — injecting into every searched train, at the census's own detection
threshold rather than the fivefold looser one the map had been using:

| | published | corrected |
|---|---|---|
| denominator | 33 (source count) | **8.27** (Σεᵢ) |
| mean efficiency | assumed 1 | **0.24** |
| sources with ε = 0 | — | **4** |
| **limit** | **f < 0.091** | **f < 0.368** |

Four sources carry too few bursts for any injected image train to beat the phase-scramble
null, so they constrain nothing while inflating a count-based denominator. Every per-source
row in the results file is bit-identical after the change; only the limit moved.

**The map's threshold was looser than the census's.** `run()` called `sensitivity_map` with
`n_scramble//5`, so `detection_p` defaulted to 2/41 = 0.049 while a real source needed
2/201 = 0.0099 to count as a detection. Cells were therefore called "sensitive" that the
census could not have detected in — again tightening the limit.

**The injection grid could not go dark where it mattered.** Seven of nine delays sat exactly
on the sidereal comb and two were deliberately off it, so the only dark cells were the
controls and the magnification axis contributed nothing (every ratio down to 0.1 was
recovered). The grid now runs to 0.02 in magnification and includes delays at ~3 and ~9
minutes off the comb, which measures the transit window's width instead of asserting it.

**The abstract targeted an object the pipeline excludes.** It opened on a lensed *one-off*
masquerading as a repeater. The ≥5-burst cut and the M_max ≥ 2 requirement exclude that
channel entirely: a galaxy-lens one-off gives 2–4 images, and **30 of 33 sources have
M_max = 1, of which 29 have p = 1.000 exactly**. Both cited theory papers are about lensed
*repeaters*, which is also what the title says. Reframed.

Also corrected: σ_DM was quoted as "2–7 pc cm⁻³ for the most active repeaters" against a
committed median of 0.4, a range of 0.2–8.4, and 9 of 33 sources sitting at the 1 pc cm⁻³
floor — only 2 of 33 fall in the quoted band. "All p ≥ 0.81" against a minimum of 0.806. And
the Roemer drift, quoted as ±150 s at Δ = 26 d, is ±230 s by the paper's own 9 s/day rule
(and 221 s from 2 × 499 s × sin(π·26/365.25)); the derived "two orders of magnitude beyond
our tolerance" is 46×.

Open: committing the per-cell sensitivity grid rather than one scalar; including the figure
(the paper has none); a lens-mass range; the top-32 truncation in the candidate-delay scan.

## Correctness pass (2026-08-22) — the retracted limit was still the first conclusion

Two defects found by the post-conversion referee round, both predating the style campaign.

**The conclusion quoted the retracted limit.** Discussion item (i) read "repeater samples are
not lensing-contaminated at the $\gtrsim$9\% level", hard-coded as a literal. That 9% is
`\flRealLimitNaive` = 0.091 — the count-based limit this slice's own referee round retracted
in favour of 0.368, and which the Results paragraph three paragraphs earlier calls "wrong by a
factor of four, in the optimistic direction". The paper argued against its own headline. Now
macro-backed (`\flRealLimit`), so it cannot drift again.

**The printed equation was not reproducible.** `\flRealEpsSum` was 8.3 and `\flRealEpsMean`
0.24, both computed over **34** sources, while the search ran on **33** and the limit divides
by the restricted sum (8.1333). So the page showed `f < 2.996/8.3 = 0.3683` where 2.996/8.3 is
0.3624. The extra source is **FRB20210601A** (eps = 0.1333): `injection_efficiency` accepts any
train with `n_bursts >= 5`, while the search additionally requires a span > 2 d.

The limit itself was always correct — `run()` already restricted to the searched names — so
**only the reported macros were wrong**, and the fix needed no re-run: `_write_macros` now
restricts `per_source` to `m["rows"]` before summing, and the macros were regenerated from the
committed real metrics JSON. `\flRealEpsSum` 8.3 → 8.1, `\flRealEpsMean` 0.24 → 0.25;
`\flRealEpsZero` (4) and `\flRealEpsMax` (0.69) are unchanged. 2.996/8.1333 = 0.3684 against
the quoted 0.3683.

Pinned by `test_macros_restrict_efficiency_to_the_searched_sources` and
`test_committed_frblens_equation_is_reproducible`, the latter checking the committed evidence
directly so the page and the JSON cannot drift apart again.
