# Findings — Jovian DAM occurrence census from Juno/Waves (plan 37)

`jansky_research.junodam`: an occurrence census over the (CML, Io-phase) plane from the public
Juno/Waves calibrated dataset (doi:10.25935/6jg4-mk86; daily CDFs, 110 chan to 40.5 MHz with
shipped per-channel Background/Sigma), using Horizons for the **sub-Juno** System III CML (the
naive IAU W_III formula errs by up to ~40° at Juno — GATE-0 verified) and Lieske (1987) for Io.

## Recover-a-known (synthetic orbit)

Real CML/Io rates + Bernoulli activity in the canonical Io-A/B/C/D boxes: injected contrast 8.75,
recovered **7.2** (cell-edge smearing accounts for the gap); 324/324 plane cells covered in a
synthetic month.

## Multi-orbit real leg (7 months / 210 days, 2016–2019, v01+v02, 8.6 GB local)

| quantity | value |
|---|---|
| 15-s bins / active | 1,209,600 / 44,294 (3.7%) |
| **duty cycle by range quartile (raw)** | **12.9% → 1.5% → 0.20% → 0.07% — raw near/far ~196×** |
| **duty cycle by range quartile (1/r² sensitivity-corrected)** | **0.25% → 0.27% → 0.20% → 0.11% — corrected near/far only 2.2×** |
| Io-box contrast, aggregate | **1.12** |
| per-month contrasts (7 months) | 1.56, **2.22**, 0.87, 0.35, 0.93, 0.70, 0.84 — median 0.87, both sides of 1 |
| distance-resolved contrast (near→far) | 1.24 / 0.36 / 1.37 / 0.71 (far quartiles activity-starved) |
| Io-box contrast, Feb-2017 single orbit (SUPERSEDED by the multi-orbit rows above) | 2.22 |

**GATE-2 caught a real convention blocker**: Io phase was computed as Λ_Io−CML instead of the
standard Φ_Io = CML+180°−Λ_Io — displacing every canonical box (except Io-B by coincidence) and
faking a weak 1.38 contrast. With the correct convention the single Feb-2017 month gave contrast 2.22 — but the multi-orbit
extension shows that was the high tail of month-to-month scatter: **across 7 months the
canonical boxes do not coherently organise Juno-frame occurrence** (median 0.87, aggregate
1.12, no strengthening with distance). One orbit can mislead; the census needed the spread. A one-sample test of the 7 monthly
contrasts against unity does not reject (t p≈0.76; sign test 2-above/5-below p≈0.23) — stated
as consistency with no coherent organisation, not proof; and genuine month-to-month Io-DAM
variability (volcanism / plasma-torus state) is a real published alternative to pure scatter,
acknowledged, not separated here.
Lesson recorded: **a synthetic round-trip cannot validate frame conventions** (it injected and
recovered in the same wrong frame). A beaming-model fit over the full v02 archive is the
definitive follow-on.

## Follow-up (2026-07-08): 1/r² sensitivity null model — the ~180× is almost all sensitivity

Prompted by the sibling `skr` slice (which found its SKR proximity trend was dominated by 1/r²
detection sensitivity), the same null model was retrofitted here. `read_waves_cdf` now also emits a
per-bin 90th-percentile channel SNR (`snr_p90`, which reproduces the `active_frac ≥ 0.1` detection
up to linear-interpolation ties at the 10% boundary — marginally stricter there), and
`sensitivity_corrected_active` distance-corrects it (SNR → SNR·(r/rref)²) and re-thresholds. Its
validity for the null rests on `snr_p90` scaling *linearly* under the correction (a self-consistent
p90 detector applied identically at every range), not on identity with `active_frac`; the raw
near/far on the same p90 detector is reported alongside the corrected one so the collapse isolates
the distance correction. Re-running the 7-month real leg:

- **Raw near/far 196×** (`active_frac` detector; 330× on the same p90 detector as the corrected
  column; 39→111 Rj quartile medians) **→ sensitivity-corrected near/far 2.2×.** The 330→2.2
  same-detector collapse isolates the distance correction (not a detector swap). Corrected
  quartiles 0.25/0.27/0.20/0.11% are **non-monotonic** (q2>q1) — no clean intrinsic trend.
- So the celebrated "proximity dominates ~180×" is a **threshold-amplified 1/r² visibility
  effect**, NOT a ~180× intrinsic occurrence rise: the DAM SNR distribution falls steeply across
  the background+5σ floor, so the modest ~8× flux change over the quartile range spread is
  amplified into the ~196× occurrence swing. Correcting for it leaves only ~2.2×.
- The residual 2.2× is **an upper bound, not intrinsic**: Juno's polar orbit couples range to
  magnetic latitude and hence CMI beaming geometry (same confound the `skr` slice flagged).
- The offline synthetic path and all pre-existing raw numbers are byte-unchanged — the null is
  purely additive. Test: `test_sensitivity_corrected_active_flattens_pure_1r2_trend`.

## Caveats

- One orbit; uneven (CML, phase) exposure (per-cell exposure carried; <3-visit cells masked).
- 10%-of-band activity criterion is one choice (suppresses narrowband RFI, dilutes narrowband DAM).
- The 1/r² null divides out inverse-square sensitivity but not the co-varying beaming geometry —
  the corrected 2.2× residual is an upper bound on intrinsic occurrence-vs-range, not a measurement.
- Units are V²m⁻²Hz⁻¹ (occurrence needs only detection vs the shipped background).
- Reproduce: download a month of CDFs to `data/junodam/` (URL pattern in module), then
  `uv run python -m jansky_research.junodam --out .` (Horizons access needed for CML).

> Cross-slice note: the 1/r² sensitivity null added here was audited repo-wide (2026-07-09) — see [sensitivity-null-audit.md](sensitivity-null-audit.md). No other merged slice needs the fix.

## Correctness pass (2026-08-23) — a stale ratio, an unevidenced validation, and a live clobber

Three defects found by the batch-5 referee round, all predating the style campaign, plus a
fourth found while fixing them.

**A hand-typed ratio, stale against the macro four lines above it.** The abstract denied "an
~180x intrinsic rise" while rendering `\jdNearFarRaw` = **196.2** in the same passage. Now
macro-backed. Same shape as `torchfdmt`'s "~24x" against a macro-derived 29x.

**The recover-a-known was not in committed evidence, and its recovered value was wrong.** The
paper hand-typed "injected contrast 8.75, recovered 7.2"; `expected_contrast` is `null` in
`results/junodam_metrics.json` and `\jdExpContrast` rendered as `--`. Re-running the offline leg
(into a tmpdir) gives injected **8.75**, recovered **6.95**, not 7.2. Both are now macros.

**The injection cannot fail near the boundary it validates.** The injected contrast 8.75 is
several times the aggregate contrast measured on the real orbits (`\jdRealContrast` = 1.12,
per-month range 0.35--2.22), so the round trip says nothing about whether a weak contrast near
unity would be recovered. The paper now says so.

**A live clobber, found while fixing the above.** `io_contrast` and `expected_contrast` mean
different things in the two run modes -- the measured Io-region enhancement on the real leg, the
recovery of an injected one offline -- and they shared **one macro name**. An offline rebuild
would have written the synthetic 6.95 into `\jdContrast`, the macro the paper uses for the real
1.12. This is the documented `\tiiNEvents` clobber, and `preserve_live_macros` cannot arbitrate
it because both runs write real values. Now namespaced `\jdRealContrast` / `\jdSynContrast` /
`\jdSynExpContrast`, with a test that reproduces the clobber and asserts it no longer happens.

Two existing tests asserted the un-namespaced names, i.e. they were pinning the defect --- the
same pattern `southern` showed. Updated.

**Also fixed:** `louis2021jgr` sat in `refs.bib` uncited. It measures latitudinal beaming from
Juno/Waves flux densities, i.e. exactly this paper's preferred interpretation of the canonical
regions and the basis of its proposed follow-up. Now cited in both places.
