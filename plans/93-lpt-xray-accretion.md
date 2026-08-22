# 93 — `lptxray`: does accretion predict radio loudness in white-dwarf binaries?

Status: 📋 planned 2026-08-20. Sourced from Rose et al. 2026 (Nature Astron. 10, 1166).

## Context

Rose et al. classify ASKAP J174508.9−505149 as an **accreting** cataclysmic variable, with
orbitally modulated X-rays and an ongoing X-ray outburst, and argue this strengthens the
link between LPTs and white-dwarf binaries. Meanwhile this repo's `wdpulsar` census searched
49 RACS-covered AR Sco-like candidates from Pelisoli et al. and found **zero** persistent
radio emission (f < 6.1% at 95%, persistent emission only, control efficiency 2/5).

Those two results sit next to each other without being connected. AR Sco itself is
famously *non*-accreting (a propeller/synchrotron system), while the Rose source is
accreting. The obvious question the pairing raises: **does accretion state predict whether a
white-dwarf binary is a radio emitter?** If the radio-loud LPTs are preferentially the
accreting systems, then the `wdpulsar` null is not just a limit — it is evidence that the
photometrically-selected, non-accreting candidates are the wrong place to look.

## The measurement

Cross-match three position lists against public X-ray catalogues:

- the 16-member LPT catalogue (`data/lpt_sample.csv`),
- the 56 Pelisoli AR Sco-like candidates (`data/wdpulsar_candidates.csv`),
- a control sample matched in magnitude/colour/sky, to calibrate the chance-match rate.

against 2RXS (ROSAT), the XMM-Newton serendipitous source catalogue, Chandra CSC, and
eROSITA-DE DR1 where the footprint allows. Then ask whether X-ray detection fraction differs
between radio-loud LPTs and radio-quiet candidates.

## The things that will make this wrong

- **Measure the chance-coincidence rate; do not assume it.** Use rigid position-shift trials
  against the real catalogues, exactly as `dr20radio` does. X-ray positional errors are
  large (2RXS ~10–30″), so the false-match rate is not negligible and is the whole ballgame
  for a fraction comparison.
- **The comparison is confounded by depth and footprint, not just by astrophysics.** X-ray
  catalogue sensitivity varies by orders of magnitude across the sky, and eROSITA-DE covers
  only the German half (l = 180–360°) — a constraint `dr20radio` already scoped out for
  rights reasons and which must be re-checked, not assumed, before use. Compare fractions
  only above a common flux limit, and expect the same lesson `dr20radio` learned the hard
  way: a **difference** in percentage points inherits the normalisation, a **ratio** does
  not.
- **Selection is not physics.** The Pelisoli candidates were selected photometrically from
  Gaia+WISE; the LPTs were selected by being radio-detected. Any X-ray difference between
  them may simply reflect that one list was built from radio detections. State this as the
  primary limitation, and treat the control sample as the thing that makes the comparison
  interpretable at all.
- **Small numbers.** Sixteen LPTs is not a sample that supports a strong claim. Pre-register
  what N would be needed for the comparison to mean anything, and if the answer exceeds what
  exists, report the fractions with intervals and no inference.

## Deliverables

- `src/jansky_research/lptxray.py` — cross-match with measured chance rate, common-limit
  fraction comparison, ratio-and-difference reporting, offline fixture.
- `results/lptxray_metrics.json`; `survey/lptxray-findings.md`.
- Write-up only if the chance-corrected fractions separate; otherwise a recorded negative.

## GATE 0

1. Novelty: has anyone cross-matched the LPT class against X-ray catalogues? (Rose et al.
   and several discovery papers report per-source X-ray follow-up; a class-wide archival
   match is the gap being claimed.)
2. Rights/footprint check on eROSITA-DE DR1 before any use.
3. Confirm a defensible control sample can be built; without it, do not start.

### GATE 0 outcome, 2026-08-22 — passed for the candidates, **CUT for the LPT leg**

1. **Novelty holds.** No class-wide, chance-corrected X-ray cross-match exists for either
   list. Per-source follow-up is well established, and the LPT review (arXiv:2601.10393)
   tabulates it, but a uniform re-search with a measured false-match rate does not exist.
   Pelisoli et al.'s own `xray` column is explicitly *not* a cross-match: the paper states
   it was compiled from prior literature and "not used in the candidate selection, but
   defined follow-up priority".
2. **eROSITA-DE is public, half-sky.** DR1 (Merloni et al. 2024, A&A 682, A34) and the
   deeper **DR2** (released 2026-07-31, ~2M sources) both cover 180 ≤ l ≤ 360 only; the
   eastern half is Russian-consortium and not public. Split: **9 of 16 LPTs and 38 of 56
   candidates** are in the public half.
3. **The LPT leg cannot be done with source catalogues, and this is the slice's main
   methodological result.** All three LPTs with a published X-ray detection
   (ASKAP J1832−0911, ASKAP J144834−685644, ASKAP J174508.9−505149) are **absent from
   every serendipitous catalogue** — 5XMM-DR15, 2RXS, eRASS1, CSC 2.1.1 — while each has
   *pointed* archival coverage: XMM 0953011101 ("VAST J1448-6856", 2024-08),
   XMM 0973390301 ("ASKAP J1745", 2025-10), Chandra 26681/26682/29265/29266 (2024). Those
   are dedicated follow-up taken *after* the catalogues were built, so the cross-match is
   structurally blind to exactly the data that produced the detections. **Recall on the LPT
   sample is 1/3.** On the candidates the same machinery recovers **16 of 16** of Pelisoli's
   independent X-ray identifications, so it has power there and none on the LPTs.

**Consequent scope change.** The two legs are reported separately and never differenced.
The headline is also softer than the plan assumed: only **one** LPT is a confirmed accretor,
so "accretion predicts radio loudness" is a one-versus-fifteen contrast, not a population
split. What the slice can honestly deliver is (a) the accretion split *within* the candidate
list, where the catalogues work, against `wdpulsar`'s radio limits, and (b) a pointed-coverage
census identifying which LPTs have archival X-ray data nobody has analysed.

## Related

`wdpulsar` (the null this contextualises), `lptv`, `dr20radio` (chance-rate measurement,
common-limit comparison, ratio-vs-difference), plans 90–92.
