# Findings — VLASS multi-epoch variability (image-vetted)

`jansky_research.vlass` cross-matches the three VLASS Quick-Look (QL) component catalogues (CIRADA
Epochs 1–2, NRAO Epoch 3; Gordon et al. 2021) and selects variable/transient candidates by the
η (weighted reduced χ², significance) and V (coefficient of variation, amplitude) metrics with a 3σ
2-D outlier cut (Rowlinson et al. 2019). First real run: a 2.5° cone at RA 190°, Dec +20° (b≈+80°),
Epochs 1–3 — 1148 anchor sources, 1022 detected in ≥2 epochs.

## Systematics handled (each one would manufacture variability)

- **Per-epoch flux scale.** QL images underestimate peak flux epoch-dependently (E1 ~15%, E2 ~7%,
  E3 ~3%; VLASS Memos 13/22). Corrected to the Perley–Butler scale, with the ~7% residual scatter
  folded into the errors. Without it the epoch offset alone flags every source.
- **A µJy/mJy unit trap.** The VizieR Epoch-1 `Fpeak` is micro-Jy mislabelled "mJy/beam" — 1000×
  off from the genuine-mJy Epochs 2–3. Caught empirically (the catalogue flux floor is the ~0.7 mJy
  detection threshold) and rescaled.
- **Deblending of crowded sources** (see below) — handled with an isolation filter.

## Result: every candidate is a catalogue artefact — an honest negative

At 3σ in both metrics, after the isolation filter (which removes 33 crowded sources), the field
yields 1–2 candidates. **None survives image vetting.** The decisive check is `image_lightcurve`,
which reads the actual VLASS QL image peak at the position in each epoch (the ground truth) and
compares it to the catalogue light curve:

| candidate | catalogue (mJy) | image peak (mJy/beam) | verdict |
|-----------|-----------------|------------------------|---------|
| J124205+193241 | 1.3 → 5.3 → 0.7 | 4.6 → 4.7 → 4.6 (flat) | steady ~5 mJy source, **deblending artefact** |
| J123318+192851 | 5.5 → 0.9 → — | 5.3 → 5.0 → 4.6 (flat) | steady ~5 mJy source, **cross-match/extraction artefact** |

The images show **steady sources in all three epochs**; the catalogue "variability" is spurious.

### Automatic confirmation by forced photometry

This vetting is now automatic. `run` calls `confirm_candidates`, which does **forced peak photometry**
(`measure_image_flux`) at the locked candidate position in every epoch's image, recomputes η/V on that
*image* light curve, and marks a candidate `image_confirmed` only if it is significantly variable
(p < 0.01) with real amplitude (V > 0.3). On this field both candidates are **auto-rejected**: forced
light curves [4.28, 3.54, 4.32] (V = 0.11) and [2.32, 1.69, 1.95] (V = 0.16) — far below the catalogue
V ≈ 1.0 and well under the threshold. Forced photometry at a fixed position is immune to the deblending
*and* cross-match failures, so the surviving candidate list (`n_image_confirmed`) is trustworthy
without manual inspection.

### The failure modes (why catalogue-only QL variability is unreliable)

1. **Component deblending.** A single slightly extended source is fit with a different number of
   Gaussian components in each epoch, so a *secondary* component's flux jumps around. Candidate 1 is
   the secondary component (1.3 mJy) of a steady ~5 mJy source whose primary sits 4–5″ away; image
   peak is flat at 4.6 mJy. The isolation filter (`isolated_mask`, drop sources with a neighbour
   within 5″) removes this class.
2. **Cross-match incompleteness.** A real source present in an epoch but with its centroid shifted
   >2.5″ is recorded as a *non-detection*, faking a fade. Candidate 2's "missing" Epoch 3 is a clear
   2 mJy source in the image — just unmatched. (Not fixed by isolation; needs forced photometry or a
   non-detection check.)
3. **Peak-flux vs image inconsistency.** PyBDSF peak flux for a blended/extended component can differ
   from the image peak by a factor of a few, inflating η.

## A larger census (50 deg², with completeness)

Scaling to a 4° cone (≈50 deg², RA 190°, Dec +20°): **3139** sources, 2826 detected in ≥2 epochs, 107
dropped as crowded, **2719 usable**. The selection flags **3** catalogue candidates — and forced
photometry **auto-rejects all three** (forced V = 0.10, 0.25, 0.16; one is a steady ~8.5 mJy NVSS
source with a missing-Epoch-2 cross-match, one a ~250 mJy source with QL high-flux scatter). So
**zero confirmed variables in 2719 sources.**

The data-driven completeness (`injection_recovery` — inject a single-epoch flare of known factor into
the real steady light curves and re-run the cut) explains why the yield is so low and bounds the
result honestly:

| flare factor | 1.25 | 1.5 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|---|
| recovered | 0.00 | 0.00 | 0.01 | 0.04 | 0.42 | ~0.50 |

The selection is **severely incomplete and saturates near 50%**: with only three epochs, a
single-epoch flare's coefficient of variation $V$ tops out at $1/\sqrt{3}\,/\,(1/3)=1.73$ as the
factor → ∞, and the noisy 3-epoch $V$ threshold (≈1.3) sits just below that ceiling, so even extreme
flares are recovered only ~half the time. The **0 confirmed variables** is therefore an *upper limit*,
not a measured rate: at ~50% completeness for strong flares, the 95% upper limit on the
strong-single-epoch-variable fraction is of order $\lesssim10^{-3}$ — consistent with, but far less
constraining than, the few-percent variable fractions the professional 2-epoch censuses report at
lower thresholds over thousands of deg².

The practical lesson: the standard 3σ-in-both-metrics cut is too stringent for 3-epoch QL data.
Because forced-photometry confirmation removes false positives downstream, a **looser** catalogue cut
(or more epochs) would raise completeness without sacrificing purity — the natural next step.

## Honest conclusion

This is a **cautionary / negative result**, in the spirit of the USS and SETI slices: catalogue-only
VLASS QL multi-epoch variability selection is **dominated by source-extraction artefacts**, and every
statistical candidate must be confirmed against the actual images before it can be believed. Across a
50 deg² census, **zero candidates survive forced-photometry confirmation**, and the selection's ~50%
completeness ceiling makes that an upper limit ($\lesssim10^{-3}$ strong variables) rather than a
measurement. The deliverable is the reproducible tool plus the QL-systematic-aware methodology —
per-epoch flux-scale correction, the deblending isolation filter, automatic forced-photometry
confirmation, and the data-driven completeness — not a discovery.

## Honest limitations

- **50 deg², single cone.** A pilot-scale census, not the thousands of deg² of the professional
  surveys; the upper limit is correspondingly weak.
- **Selection saturates at ~50%** for single-epoch flares (3-epoch ceiling). A looser cut — viable now
  that forced photometry cleans up false positives — or more epochs is needed for a real rate.
- **Isolation is anchor-epoch only.** A source isolated in Epoch 1 but blended in a later epoch is
  not yet caught; a full treatment checks isolation per epoch and detects ambiguous matches.
- **Completeness is candidate-driven, not all-source.** Forced photometry runs on the candidates, so a
  real variable whose *catalogue* light curve was flattened by extraction issues can still be missed
  before it ever reaches confirmation; forced-photometering every source would close this.
- A genuine variable population almost certainly exists at this depth (the literature reports a few
  percent at 2σ over larger areas); reaching it needs much more sky and a looser, completeness-aware
  cut — future work.

## Bottom line

The pipeline runs end-to-end on real VLASS data with the QL systematics handled, and it correctly
finds that this field's variability candidates are **all extraction artefacts** — demonstrated by
image vetting, which is therefore mandatory. A tooling + methodology + honest-negative contribution.
