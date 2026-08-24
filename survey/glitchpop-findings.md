# Findings — JBO glitch waiting-time classification + post-2018 delta (plan 48, F11)

`jansky_research.glitchpop`: a monitoring-gap-robust per-pulsar waiting-time classification of the
Jodrell Bank glitch catalogue (exponential / quasi-periodic / clustered), plus the diff of which
classifications change when the post-2018 glitches are added. Howitt+2018 already classified
per-pulsar waiting times (Poisson vs quasi-periodic) for the few most-glitching pulsars pre-2018; the
wedge here is extending that to the full ≥5-glitch sample, adding monitoring-gap robustness, and the
post-2018 change table. Differentiated from Basu+2022 (glitch-SIZE mixture) and Zhu&Zheng 2025
(glitch-CLUSTER period vs age).

## GATE 0 (2026-07-11) — PASS (conditional)

- **Data LIVE + scrapable**, no auth: `jb.man.ac.uk/pulsar/glitches/gTable.html`, **727 glitches /
  222 pulsars** (up from Basu+2022's 543/178 → ~184 new post-2018). 11 columns: idx, JNAME, Bname,
  glitch#, MJD+err, dF/F (Δν/ν ×1e-9)+err, dF1/F1 (Δν̇/ν̇ ×1e-3)+err, References. The JNAME cell mixes
  J- and B-names (Vela = **B0833-45**, Crab = B0531+21); carried forward when blank on repeat rows.
- **Basu+2022** (MNRAS 510, 4049; arXiv:2111.06835): catalogue + two-component-Gaussian SIZE mixture;
  end-2018; NO per-pulsar waiting-time model selection.
- **IN-PRESS check**: no JBO catalogue-update/population paper post-2022. **Zhu & Zheng 2025**
  (arXiv:2501.01862) use the same ≥5-glitch JBO+ATNF sample but for glitch-CLUSTER period vs age — a
  different deliverable, cited and differentiated (we do not do clustering).
- **Recover-a-known**: most pulsars exponential (Poisson, Howitt+2018); J0537-6910 (~100 d) and Vela
  quasi-periodic, NOT exponential.

## Method — and a monitoring-gap insight (the load-bearing step)

- **Waiting times**: scrape → group by name → successive-MJD diffs (years) for pulsars with ≥5
  glitches.
- **First attempt failed the recover-a-known**: the naive classifier called **J0537-6910
  EXPONENTIAL** — because its recorded intervals include a **2264-day monitoring gap** (the
  RXTE→NICER hiatus) among ~100-day intervals, which inflates its raw CV to ~2. The archival
  catalogue conflates real inter-glitch intervals with monitoring gaps.
- **Fix — monitoring-gap excision**: drop waits > 6× the pulsar's median (almost certainly
  unobserved intervals for a regular glitcher; a genuine exponential loses only ~1.5% far tail). This
  recovered J0537 (cv 0.53 → quasi-periodic) AND Vela (cv 0.57 → quasi-periodic). **Cost**: a
  genuinely *clustered* pulsar (bursts + real long gaps) has its long gaps excised too, so long-gap
  clustering is under-detected — stated plainly (and the committed census has 1 nominal clustered call (chance-consistent; see the revision note below)).
- **Classification — parametric-bootstrap CV test** (replaced a gamma-vs-exponential BIC that
  underweighted Vela's mild-but-real regularity): CV = std/mean of excised waits, calibrated against
  the exponential null by bootstrapping the CV of n exponential waits. Quasi-periodic if CV
  significantly low (p<0.05), clustered if significantly high, else exponential.
- **Post-2018 delta** (`classification_delta`): reclassify on the pre-2019 (Basu) subset vs the full
  catalogue; report `newly_classifiable` (crossed ≥5 only with new data) and `flipped` (real class in
  both epochs, changed).

## Synthetic recover-a-known (offline, CI)

- Exponential vs quasi-periodic separated with accuracy rising with glitch count (exp 0.85→0.93,
  qp completeness 1.0), FP ~0.075.
- A quasi-periodic series with an injected monitoring gap: misclassified without excision, recovered
  with it (the J0537 scenario in miniature).
- Clustered recoverable only with excision off (its long gaps look like monitoring gaps).

## The real census (ran 2026-07-11): a working classification with a passing recover-a-known

**727 glitches / 223 pulsars; 3 magnetars dropped; 31 with ≥5 glitches classified.**

- **Recover-a-known PASSES**: J0537-6910 and Vela (B0833-45) both **quasi_periodic**
  (`known_quasiperiodic_ok = yes`) — required the gap excision (without it J0537 fails). Crab
  (B0531+21) comes out exponential, consistent with it being less regular than Vela/J0537.
- **Census**: **20 exponential, 10 quasi-periodic (~32%), 1 clustered**. Only **1.55** false
  quasi-periodics are expected under an all-Poisson null, so the QP *fraction* is highly significant
  (binomial **p≈2e-6**) — but with a per-pulsar p<0.05 test ~1–2 borderline members are expected
  spurious, so the class is asserted, not every member. The 1 clustered is a lower bound (gap
  excision suppresses long-gap clustering — the signal Zhu&Zheng 2025 report on this same sample).
- **Post-2018 change**: **1 flip** among the 28 pulsars classifiable in both epochs — J2229+6114
  (exponential→quasi_periodic, on its 9th glitch, near the p=0.05 boundary) — plus **3 newly
  classifiable** (B1830-08 exp, J2021+3651 qp, B2224+65 clustered). Explicitly a data-driven
  statement, not physical: a single new glitch moves a borderline case (≤10 waits/pulsar).
- **Sensitivity**: J0537 → quasi_periodic and Vela → quasi_periodic for gap_factor ∈ [3,10]; the QP
  count is stable (8–11). Vela passes with ZERO gap excision; J0537's gap is ~22× its median (excised
  by any threshold 3–20×) — so the recover-a-known is not sensitive to the specific 6×.

**Honest limits**: incomplete/uneven glitch detection (regularity is a lower bound); the ≥5-glitch
cut selects active glitchers; gap excision assumes large waits are monitoring gaps (so it CANNOT
detect genuine long-gap clustering — a method asymmetry); flips are single-glitch-sensitive; the
recover-a-known is partly in-sample (the statistic was reselected after a gamma-BIC failed Vela — the
synthetic injection test is the independent validation).

## GATE-2 (2026-07-11) — PASS, no blockers

The reviewer confirmed the numbers are internally consistent, the recover-a-known passes, the
1.5% excision-loss is exact (2^-6), the J0537 monitoring gap is physically real (RXTE→NICER hiatus),
and all citations check out. Required should-fixes, all applied:

- **Howitt+2018 under-credited**: it already did per-pulsar waiting-time classification (and even
  reclassified B1338-62 to quasi-periodic — which we independently reproduce). Novelty restated as
  extending Howitt to the full ≥5-glitch sample + gap-robustness + the post-2018 delta.
- **Multiple testing**: added `population_significance` — 10 QP vs 1.55 expected → binomial p≈2e-6
  (aggregate significant), but ~1–2 borderline members expected spurious (stated).
- **"0 clustered" was defined away**: reframed as a method asymmetry (gap excision removes the
  long-gap clustering Zhu&Zheng detect on the same sample), not evidence against clustering.
- **Method-iteration circularity disclosed**: the CV-bootstrap replaced a gamma-BIC that failed Vela,
  so the two-object recover-a-known is partly in-sample; the synthetic test is the independent one
  (and Vela passes with zero excision, J0537 for any threshold — so not sensitive to 6×).
- **Bootstrap null now replicates the gap excision** (identical processing of null and data).
- **Magnetars dropped** (3; X-ray-outburst-driven, not rotation-powered glitches).
- **Sensitivity sweep** on gap_factor added (robust over [3,10]).
- Nits: stale 222→223 pulsar count fixed; unused ATNF citation removed (the Ė join was not
  implemented; the wedge is waiting-time class, not spin-down).

## Reproduce

Offline (classifier + synthetic recover-a-known + vendored-sample parser test):
`uv run python -m jansky_research.glitchpop --offline --out .`
Real (scrapes the live JBO table): `uv run python scripts/glitchpop_real.py`
(writes `results/glitchpop_metrics.json` + `results/glitchpop_census.json`, `is_real=True`).

## Full referee round (2026-08-24): MAJOR REVISION, 19 findings, two BLOCKERs

The gap-excision insight is real and load-bearing; the macros, the binomial arithmetic and all
four DOIs check exactly. The problems: the title result's stability, the sample accounting, and
regenerability.

**BLOCKERs:** (1) The headline flip (J2229+6114, exponential -> quasi-periodic) has margin
p = 0.0495 vs 0.05 -- TWO bootstrap replicates out of 4000, 0.15 MC-SE from its own threshold; a
different seed flips the title result about half the time, and the pre-2019 margin isn't recorded
at all. B1758-23 (0.0473) is in the same condition. (2) The abstract's first sentence counts
727 glitches on the UNFILTERED catalogue but 220 pulsars on the magnetar-FILTERED sample -- the
source of the 220/222/223 confusion across the findings file, refs.bib and the module docstring.

**MAJORs:** the census is a live scrape with no committed snapshot (three scrapes, three counts;
not regenerable even in principle -- fix: commit the 727-row parsed table + --from-csv). The
synthetic validation cannot fail: injected qp at cv 0.12 vs real detections at 0.18-0.62, hence
completeness exactly 1.0 at every grid point; the selection function at cv~0.5, n~6 (where 4 of
10 detections and the flip live) was never measured. \gpSynFP=0.0938 pools false-qp,
false-clustered AND insufficient-dropped -- not "a Poisson series called regular". The population
p=1.6e-6 uses the nominal 0.05 rate the paper's own validation contradicts (measured 0.0938;
recomputed P(>=10|31,0.0938)=3.7e-4 -- 220x weaker, still significant). "The single clustered
source is a lower bound" asserts existence where 1 detection IS the chance expectation (1.55),
resting on 4 retained waits at the min_waits floor. "31 pulsars with >=5 glitches" is really
survived-TWO-cuts (excision can drop a >=5-glitch pulsar entirely; the drop count is unrecorded
and biased against clustered candidates). The gap-factor 3-10x stability claim has no committed
sweep, and the range in the findings prose (8-11 qp) is a +/-15% swing in the headline. The
out-of-sample anchor EXISTS and is unused: Howitt+2018's independent B1338-62 reclassification,
which the committed census reproduces (cv 0.555, p 0.0) -- promote it, ideally as a full agreement
table.

**MINOR:** the null is never calibration-checked against the real data, and the committed
p_regular values hint it may be off (mean 0.28 vs 0.5 expected; z=-2.3 on the exponential
subset) -- either the population is more regular than Poisson wholesale, or the null is too
dispersed; both interesting, neither examined. All pulsars with the same wait count share ONE
4000-draw null (seed never varied); p=0.0 occurs (invalid MC p). Documentation contradicts the
JSON in four places (README says 33/23 and omits the clustered class entirely; findings says
222/223 pulsars and "0 clustered"). The paper contains neither the census table nor any figure.
Three verbs stronger than the evidence (title's "Change"; "still shifting"; "warranted").
The synthetic leg has no committed results file.

**Status: fixes pending** (needs one live scrape to commit the snapshot; everything else offline).

### Resolved 2026-08-24 (revision): the flip was Monte-Carlo noise, and the census is pinned

All 19 findings addressed. The two blockers resolved the strong way:

**The title flip dissolved exactly as predicted.** At n_boot = 2e5 with per-pulsar (crc32) null
seeds and (k+1)/(B+1) p-values carrying their MC standard errors, J2229+6114 does not flip --
nor does B1758-23, the other MC-marginal call, survive as quasi-periodic. n_flipped = 0, n_qp
10 -> 8. The aggregate excess is untouched: nominal binomial p = 1.1e-4 and Poisson-binomial
p = 7.9e-5 under the empirically measured count-dependent false rates. The title became "the
Post-2018 Increment" and the paper states the finding plainly: no classification change is
secure at current per-pulsar sample sizes.

**The census is pinned to a committed snapshot** (`results/glitchpop_glitches.csv`, 728 rows /
224 pulsars, retrieved 2026-08-24 -- a fourth distinct count in four retrievals, proving the
point) with a `from_csv` analysis path; the paper states the retrieval date and cites jbocat.
The abstract's mismatched sample counts are fixed with matched pairs (728/224 raw, 714/221
analysed).

**"~184 post-2018 glitches" was catalogue-vintage arithmetic**: measured against the Basu epoch
in this snapshot, the post-2018 increment is **89**, and 21 pre-2019 glitches are retroactive
additions (the is_new flag, previously parsed and discarded, now counts them).

**The validation can now fail, and does**: the injection surface (cv 0.1-0.7 x n 6-40) shows
completeness 1.0 only for cv <= 0.3 at n >= 10, falling to 0.44 at (cv 0.5, n 6) and 0.16 at
(0.7, 6) -- the regime where the borderline calls live; the census's QP count is a floor. The
pooled "false-positive rate" 0.0938 is disaggregated (false-qp 0.04-0.12 by count,
false-clustered 0.04-0.10, dropped-by-excision <= 0.02), and the count-dependent false-qp rates
feed an exact Poisson-binomial population significance.

**The clustered claim is withdrawn to chance**: B2224+65 (4 retained waits, p = 0.039) against
an expected 1.55 false clustered (binomial p = 0.80); no clustering claim in either direction,
deferring to Zhu & Zheng. The hidden second cut is stated (32 pulsars pass >=5 glitches, 1 is
dropped when excision leaves too few intervals). The gap-factor sweep is committed (QP 8-10
across 3-10x; the clustered call vanishes below 6x). The out-of-sample Howitt anchor
(B1338-62, cv 0.555) is promoted into the paper as the validation the reselected statistic
needed. And the p-uniformity diagnostic is committed and REAL: the census's p_regular values
are non-uniform (mean 0.28, KS p = 1e-4) -- either the population is mildly regular wholesale
or the null lacks a monitoring dead-time; stated as the principal open methodological question.

Doc mismatches fixed (README row, this file's "0 clustered"); the synthetic paths run at 4000
boots so the census-grade default does not slow CI (test file 2:03 -> 3.5 s).
