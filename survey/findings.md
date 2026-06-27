# Findings — FRB burst-statistics on CHIME/FRB Catalog 1 (post GATE-2 review)

Run of `jansky_research.pipeline` on the real **CHIME/FRB Catalog 1** (CHIME/FRB Collaboration
2021, ApJS 257, 59; arXiv:2106.04352). The public CSV has 600 rows; these are **536 independent
events** stored as multiple `sub_num` components per multi-part burst, so the pipeline keeps one
row per event (`sub_num == 0`) — **62 repeater bursts from 18 sources + 474 non-repeaters**.
Treating sub-bursts as independent would pseudo-replicate near-identical DMs and inflate
significance, so it is avoided. Numbers below are `results/metrics.json`; this is the honest
interpretation, revised after the GATE-2 science review.

## Result 1 — Burst width: repeaters are wider (reproduces the literature)

KS test on temporal width: **D = 0.45, p ≈ 2e-10**; repeater median **2.0 ms** vs non-repeater
**0.94 ms**. Repeaters being **wider** is the established CHIME Cat 1 morphology result
(**Pleunis et al. 2021**, ApJ 923, 1; arXiv:2106.04356). Recovering it independently from the
public CSV is the tool's **validation finding**.

## Result 2 — Fluence/energy power law (defensible measurement)

Differential fluence distribution dN/dF ∝ F^(−γ) via the Clauset–Shalizi–Newman / Hill MLE with
KS-based x_min selection: **γ = 2.54 ± 0.13** above **f_min = 7.8 Jy ms** (138-event tail). This
matches the Cat 1 paper's reported cumulative slope α = −1.40 ± 0.11 (differential γ = 1 + |α| =
2.40) within ~1σ. *Caveats:* the auto-selected f_min = 7.8 Jy ms sits **above** CHIME's nominal
~5 Jy ms peak-sensitivity threshold — expected, since the effective completeness of a population
spanning many DMs/widths/transit positions is higher than the peak figure; the fit is a single
power law not corrected for the DM/sky-dependent selection function.

## Result 3 — Wait-time clustering (caveated; not a headline result)

Pooled **within-source** inter-burst waits (44 waits across 18 repeaters) fit a Weibull with
shape **k = 0.41 (95% CI 0.32–0.56)** — i.e. clustered (k < 1), broadly consistent with
dense-monitoring estimates for FRB 121102 (k ≈ 0.34; Oppermann et al. 2018). **Reported with a
caution, not as a measurement:** CHIME is a transit instrument (~one look/source/day), so the
waits remain cadence-biased (bimodal short/long structure). Note that removing the sub-burst
pseudoreplication moved k from a spurious 0.14 to this plausible 0.41 — itself a lesson in why the
event-level treatment matters. The k value and CI appear **only** in the methods/pitfalls
discussion, never in an abstract or summary table.

## Additional observed differences (exploratory; selection-affected)

- **DM:** KS D = 0.47, p ≈ 1e-11; repeaters **lower** (350 vs 565 pc cm⁻³). *Not* a reproduction
  of Pleunis et al. 2021, who found the Cat 1 repeater/non-repeater DM distributions statistically
  **consistent**. A significant difference at larger samples was later reported by **CHIME/FRB
  Collaboration 2023** (ApJ 947, 83; arXiv:2301.08762), which attributes the lower repeater DM
  partly to selection (nearer sources are easier to re-detect). We present this as exploratory.
- **Fluence:** KS D = 0.23, **p ≈ 6e-3** (marginal); repeaters slightly higher (4.9 vs 3.6 Jy ms).
  Plausibly a selection effect (repeaters identified via multiple detections); not literature-
  validated. Reported as a weak difference.

## Assessment for GATE 2

- **Honest & non-trivial:** Result 1 reproduces a peer-reviewed result (tool validation); Result 2
  is a defensible measurement matching the catalogue paper; the wait-time and DM/fluence results
  are explicitly downgraded with their selection/cadence caveats and correct citations.
- **Paper framing:** a **reproducibility / tooling contribution** — "a lightweight, tested,
  CPU-only, offline-reproducible FRB burst-statistics tool, validated by recovering the CHIME
  Cat 1 width result" — **not** a discovery claim. Only **width** carries "reproduces the
  literature"; DM/fluence are "additional observed, selection-affected differences."
