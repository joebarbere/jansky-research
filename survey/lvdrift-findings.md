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
