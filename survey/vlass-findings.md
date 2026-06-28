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

## Honest conclusion

This is a **cautionary / negative result**, in the spirit of the USS and SETI slices: catalogue-only
VLASS QL multi-epoch variability selection is **dominated by source-extraction artefacts**, and every
statistical candidate must be confirmed against the actual images before it can be believed. On this
field, **zero candidates survive image vetting**. The deliverable is the reproducible tool plus the
QL-systematic-aware methodology — per-epoch flux-scale correction, the deblending isolation filter,
and the **mandatory `image_lightcurve` ground-truth check** — not a discovery.

## Honest limitations

- **One field, conservative cut.** A 2.5° cone with a strict 3σ-in-both-metrics selection; this is a
  pipeline + methodology validation, not a population census.
- **Isolation is anchor-epoch only.** A source isolated in Epoch 1 but blended in a later epoch is
  not yet caught; a full treatment checks isolation per epoch and detects ambiguous matches.
- **Non-detections are not yet forced-photometered**, so real sources missed by the cross-match
  (failure mode 2) still leak in until image-vetted.
- **`image_lightcurve` uses the cutout peak**, adequate to confirm flatness but not a calibrated
  light curve.
- A genuine variable population almost certainly exists at this depth (the literature reports a few
  percent at 2σ over larger areas); a real census needs much more sky, completeness vs amplitude, and
  forced photometry — future work.

## Bottom line

The pipeline runs end-to-end on real VLASS data with the QL systematics handled, and it correctly
finds that this field's variability candidates are **all extraction artefacts** — demonstrated by
image vetting, which is therefore mandatory. A tooling + methodology + honest-negative contribution.
