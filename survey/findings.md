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

## Full referee round (2026-08-25): MAJOR REVISION, 18 findings, one BLOCKER

The headline validated result is real and survives everything the referee threw at it:
repeaters are wider at burst level, at SOURCE level (per-source-median KS D = 0.563,
p = 1.1e-5), under a source-cluster bootstrap, and after deleting the two dominant sources.
All 25 macros reproduce exactly from data/chimefrbcat1.csv; all six DOIs clean; disowning k
was the right call.

**BLOCKER: Figure 1 plots pooled-across-source waits under a within-source caption.**
report.py calls `wait_times` (not `grouped_wait_times`), so 44 of the 61 plotted intervals
cross source boundaries — the quantity Methods itself calls "meaningless" — while the overlaid
Weibull was fitted to the 44 within-source waits. The mismatch displays a good fit as a bad
one: max |empirical − curve| = 0.307 (KS p = 1.3e-5) as plotted, vs 0.110 (p = 0.39) against
the data actually fitted. It survived because the offline fixture's repeaters are a single
source, so the multi-source path is never exercised.

**MAJORs:**
- Figure 2's caption promises "the maximum-likelihood power law" — no fit is drawn and F_min
  is unmarked; the one place a reader could judge the fit has been removed from the evidence.
- The γ error conditions on an x_min chosen from the same data: joint bootstrap gives
  sd 0.29 (95% [2.07, 3.19]) vs the quoted ±0.13 — 2.2× too small — and F_min's own 95%
  interval is [2.9, 15.7]; γ runs 2.09→3.08 across x_min 3→20.
- "Matches the cumulative source-count slope" has no number, citation, or transform — and the
  comparand sits in this repo's own survey/findings.md (Cat 1's α_cum = −1.40 ± 0.11 →
  γ_diff = 2.40; the match is real at 0.84σ). Put it in the sentence.
- 62 bursts from 18 sources treated as 62 independent draws: two sources supply 48% of the
  bursts. Source-level restatement: width D = 0.563, p = 1.1e-5 (READS BETTER than the
  burst-level 0.45); DM collapses from p = 1e-11 to p = 0.044 — the marginal value the paper
  already believes.
- report.write_macros never calls preserve_live_macros (the CLAUDE.md invariant), and
  build_catalog silently falls back to synthetic on ANY fetch failure with out_dir=".": a
  documented one-liner (make pipeline ARGS=--offline, or make reproduce offline) rewrites the
  paper's macros to the synthetic values (\nBursts 600, D 0.08, p 0.4) while the JSON is
  preserved and guard-real stays green — split-brain evidence with every gate passing.
- The gap claim ("no pip-installable library") describes a property the artifact lacks (no
  PyPI release; the dependency resolves from a git tag); the tool is named frbstats while the
  paper contrasts it with the FRBSTATS web platform; the whole landscape paragraph carries
  zero citations.

**MINOR:** two undisclosed cuts (26 width upper limits, ALL non-repeaters, dropped — including
them is conservative, D→0.456; six fluence=0.0 placeholders inside the KS and \medFluenceOne);
the catalogue's own excluded_flag (39 events) ignored (γ 2.543→2.584 without them); "44 waits
across 18 sources" is 16 (two single-burst sources contribute none); the Weibull CI resamples
waits not sources (cluster CI 0.32–0.69; conclusion unchanged) and the asserted bimodality has
a clean uncomputed number (7 same-transit + 8 sidereal-day waits of 44); the one hand-typed
number (0.14, reproducible but uncommitted); "temporal width" never defined (width_fitb of the
leading sub-burst; bc_width gives 8.85/3.93 ms — same sign, factor-4 different medians — and
the validation passes under every implementation error tried: it needs Pleunis' own numbers as
comparand); the selection function named three times, used zero; the 5 Jy ms threshold
uncited; no data provenance in either manuscript (third-party mirror URL, no checksum,
\catalogSource unused); stale arXiv package; RNAAS \software omits jansky-research and
attributes the tool to the course.

**Status: fixes pending.** The single change: make the SOURCE the unit of analysis for every
repeater statistic — width quoted as D = 0.56, p = 1e-5 in both abstracts, DM restated at its
source-level 0.044, and Methods stating that 62 bursts from 18 sources are not 62 draws.
