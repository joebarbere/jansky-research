# Findings — FRB burst-statistics on CHIME/FRB Catalog 1 (pre-GATE-2)

Run of `jansky_research.pipeline` on the real **CHIME/FRB Catalog 1** (600 bursts: 94 from 18
repeaters + 506 non-repeaters; CHIME/FRB Collaboration 2021, ApJS 257, 59). Numbers are the
machine output in `results/metrics.json`; this document is the **honest interpretation** that
GATE 2 (science-reviewer + human) must sign off before any paper drafting.

## Result 1 — Repeater vs non-repeater differences (robust; reproduces the literature)

Two-sample KS tests show repeaters and non-repeaters differ significantly in every measured
property:

| property | KS *D* | *p* | direction (median) |
|----------|:------:|-----|--------------------|
| width    | 0.39   | 1.9e-11 | repeaters **wider** (2.0 ms vs 1.0 ms) |
| DM       | 0.46   | 5.7e-16 | repeaters **lower** (349 vs 563 pc cm⁻³) |
| fluence  | 0.29   | 1.4e-6  | repeaters **higher** (5.7 vs 3.8 Jy ms) |

The **larger temporal width of repeaters** is the well-established CHIME Cat 1 result (Pleunis et
al. 2021; CHIME/FRB 2021). Recovering it from the public CSV with this independent tool is the
**validation finding** — it shows the pipeline reproduces a known, peer-reviewed result. *Honesty
note:* these are catalogue-level differences subject to CHIME's selection function (e.g. the
lower repeater DM partly reflects that nearer sources are easier to re-detect); we report them as
*observed catalogue differences*, not intrinsic population claims.

## Result 2 — Fluence/energy power law (robust with a proper completeness cut)

Differential fluence distribution dN/dF ∝ F^(−γ): **γ = 2.38 ± 0.10** above a Clauset-selected
lower bound **f_min = 7.7 Jy ms** (183-burst tail). The automatic x_min selection matters — fitting
from the global minimum instead gives a spuriously low γ ≈ 1.41 by including the incomplete faint
end. γ ≈ 2.4 (cumulative slope ≈ 1.4) is consistent with published CHIME source-count analyses.
*Honesty note:* a single power law above one completeness cut; not corrected for the
DM/sky-dependent selection function.

## Result 3 — Wait-time clustering (caveated; a methodological illustration, not a claim)

Pooled **within-source** inter-burst waiting times (76 waits across 18 repeaters) fit a Weibull
with shape **k = 0.14** (95% CI 0.12–0.16) — formally "highly clustered". **This is dominated by
CHIME's transit-instrument cadence, not intrinsic clustering:** each source is seen ~once per
sidereal day for minutes, so waits are bimodal (short intra-transit vs long inter-day), which
drives any single-Weibull fit to k ≪ 1. We therefore present this **not** as a clustering
measurement but as an illustration of why dense, single-source monitoring (à la FRB 121102) is
required to measure burst clustering — and as a caution the tool surfaces honestly.

## Overall assessment for GATE 2

- **Non-trivial & honest:** yes — Result 1 reproduces a known peer-reviewed result (tool
  validation); Result 2 is a defensible measurement with the right completeness treatment.
- **No overclaiming:** the wait-time result is explicitly downgraded to a cadence-bias
  illustration; population differences are framed as catalogue-level/selection-affected.
- **Framing of the paper:** a *reproducibility / tooling* contribution — "a lightweight, tested,
  CPU-only, offline-reproducible FRB burst-statistics tool, validated by recovering the CHIME
  Cat 1 repeater/non-repeater differences" — **not** a discovery claim.
