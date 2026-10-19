# lvdrift findings — the within-burst type III drift law from LV614 (plan 101)

Running log, newest section last. Plan: `plans/101-lvdrift-within-burst.md` (controls frozen
2026-10-18, PR #349). Numbers are in committed `results/lvdrift_*.json` unless marked.

## Step 1 — estimator and false-failure rates (2026-10-18; no LV614 array opened)

`lvdrift.burst_index` implements the frozen estimator. One addition found while testing: with no
frequency leverage SciPy reports an error of exactly 0 on a free index, so a zero or undefined
error (or a negligible span) is flagged `unconstrained` and treated as a failed fit, like a fit
on the parameter bound.

False-failure rates, simulated before either rule is used (results/lvdrift_falsefailure.json;
white noise, not yet the real noise windows):
- **P1** (median recovered a within ±0.15, 50 plants per index): 0 failures in 40 trials.
- **N1** (≥ 95% of channel-shuffled bursts rejected): 300 / 300 rejected.
Neither rule needs widening.

**A systematic low bias at high indices**, small in spread (σ of trial medians ≈ 0.01) but not
small against the comparison values: planted 1.00 → 0.95, **1.23 → 1.14**, 1.84 → 1.79; planted 0
and 0.43 are unbiased to 0.01. It is inside the P1 tolerance, so the rule passes, but it is the
same size as the difference the slice wants to measure against Zhang et al. (2018, index 1.23).
Recorded now, before any real array: the within-burst result will be reported with the P1 bias
measured in the **real** noise windows alongside it, and any comparison with 1.23 or 1.84 will
state it. Likely cause (not tested): 1 s arrival-time quantisation and the 3 × 3 smoothing
flattening t(f) at the fast, high-frequency end.

## Step 2 setup — real noise windows (2026-10-18)

Zenodo zip downloaded (MD5 755862f4…, matching the record); 290 type III arrays. Array lengths:
median 60 s (59–2820 s); the summed flux peaks at a median of 51 s, but 13% peak inside the first
20 s, so "the first 20 s" is not pre-burst for all. **Refinement before P1 runs:** noise windows
are the first 20 s of arrays whose summed flux peaks at ≥ 30 s (10 s margin) — chosen from flux
timing only, no drift information. Plants use t_hi = 3 s; even the slowest planted index crosses
the band in ≲ 12 s, so a planted burst fits inside a 20 s window.
Channel layouts differ (from the file names): 285 arrays use subbands 150–311 (29.3–60.7 MHz), 4
use 100–343 (19.5–67.0 MHz), 1 uses 150–271. Each array's frequencies are parsed from its name.
P1 uses windows from the 285 standard-layout arrays only, so every plant shares one frequency grid.

## P1 in real noise — PASS (2026-10-18; results/lvdrift_p1_real.json)

217 noise windows (first 20 s of standard-layout arrays peaking ≥ 30 s); 50 plants per index at
S/N 5–50, FWHM 2–6 s. Median recovered a (bias): 0 → −0.001 (−0.001); 0.43 → 0.448 (+0.018);
1.00 → 0.983 (−0.017); **1.23 → 1.118 (−0.112)**; 1.84 → 1.820 (−0.020). All within ±0.15: **P1
passes**. The low bias near 1.23 seen in white noise (−0.09) persists in real noise (−0.11) while
1.0 and 1.84 are nearly unbiased here, so it is specific to that range, not a uniform compression;
any comparison with Zhang et al.'s 1.23 must carry it. Only 25–30 of 50 plants per index pass the
quality cut (the weak ones keep < 40 channels above S/N 5): the real-data result describes bursts
bright enough to track, and the fraction of the 290 that pass will be reported.

## Step 3 — N1, R0 and the within-burst index on the real arrays (2026-10-18; results/lvdrift_real.json)

- **N1 — PASS.** 290 / 290 channel-shuffled real bursts fail the quality cut (rule ≥ 95%).
- **Within-burst index (headline setting, S/N 5, smoothing 3):** **43 of 290** type III bursts
  pass the quality cut (15%); median a = **1.45**, bootstrap 95% interval **1.25 – 1.96**. Most
  failures are "non-monotonic", often by ≳ 50 s — a ridge that jumps to another burst or to
  interference in the same array rather than a slow drift. The result describes a 15% subsample
  selected by the cut, and that selection is the first caveat. Against Zhang et al. (2018, 1.23 at
  10–80 MHz): the interval's lower edge sits just above 1.23; P1 showed the estimator reads a true
  1.23 as ≈ 1.12 in this noise, so the bias runs toward agreement, not away from it.
- **R0 — FAILS its frozen rule.** Per burst, our peak-time drift (±3 s regression) over Table 3's
  drift has median **1.34** for the 201 joined bursts with both values (rule [0.8, 1.25]). Means:
  type III **4.16** MHz s⁻¹ (paper 4.383; Table 3's own mean of the listed drifts 4.24), type II
  **1.57** (paper 0.043) — our ridge method does not capture slow, long type II bursts. Join:
  281 unique, 8 ambiguous (two bursts in one minute), 1 unmatched.
- **The paper's across-burst law does not reproduce from its own Table 3:** unbinned log-log fit
  a = **0.23** (229 bursts); the binned, mean-flux-weighted recipe as we read it gives **−3.99**,
  against the published 0.43 ± 0.31. The weighting convention is not stated (Sect. III.5).

R0 failing means our reconstruction of *their* per-burst quantity does not match theirs; it does
not test the within-burst estimator, which P1 validates. Both facts go in the note. Robustness
(S/N 4/5/7, smoothing 1/3/5) follows.

## Robustness grid (2026-10-18; results/lvdrift_robustness.json)

| S/N min | smoothing | bursts passing | median a | 95% interval |
|---|---|---|---|---|
| 4 | 1 | 33 | 1.82 | 1.29 – 2.32 |
| 4 | 3 | 42 | 1.42 | 1.25 – 1.96 |
| 4 | 5 | 51 | 1.66 | 1.07 – 2.12 |
| 5 | 1 | 34 | 1.75 | 1.29 – 2.23 |
| 5 | 3 | 43 | 1.45 | 1.25 – 1.96 |
| 5 | 5 | 52 | 1.62 | 1.03 – 2.12 |
| 7 | 1 | 33 | 1.81 | 1.29 – 2.32 |
| 7 | 3 | 44 | 1.42 | 1.25 – 1.92 |
| 7 | 5 | 52 | 1.62 | 1.03 – 2.12 |

Medians span 1.42 – 1.82, always above 1.23, but the 95% interval
contains Zhang et al.'s 1.23 in 3 of 9 settings (all at smoothing 5, where more, noisier bursts pass).
**The within-burst index is not robustly different from 1.23.** What the data support: a
within-burst index of about 1.4–1.8 for the 12–18% of bursts that can be tracked, consistent
with 1.23 at some settings and well above the paper's across-burst 0.43. The interval width is
driven by the small passing sample (33–52 bursts).

