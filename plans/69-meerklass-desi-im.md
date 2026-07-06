# 69 — MeerKLASS 2019 intensity-mapping cube × DESI: an independent-tracer cross-correlation

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

The only published single-dish HI intensity-mapping × galaxy cross-correlation detection used
WiggleZ galaxies (7.7σ, arXiv:2206.01579, MeerKLASS 2019 pilot data). An independent-tracer
check against DESI spectroscopic galaxies in the overlapping footprint is open (fable-ideas
F32): same cube, different galaxy sample, published pipeline followed step by step. Foreground
cleaning is the hard part of this field — this slice does not innovate there; it follows the
published PCA-cleaning recipe and frames every result conservatively as a cross-check. The
deliverable is a confirmation (or an honestly-bounded non-confirmation) with a second tracer.
GPU-friendly (pure-PyTorch FFTs fine on ROCm); cube footprint is modest against the ~275 GB disk.

## Deliverables

- `src/jansky_research/meerklass.py`: `fetch_cube` / `fetch_desi_targets` (network, `# pragma`),
  `pca_clean` (published-recipe foreground removal, mode count as a reported parameter),
  `cross_power` (cube × galaxy-overdensity cross-power with transit-function weighting),
  `null_tests` (rotated/shuffled galaxy catalogues), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/meerklass/`; `survey/meerklass-findings.md`; wiring.

## Approach

0. **GATE 0:** verify the MeerKLASS 2019 pilot cube is publicly downloadable (the scan did not
   live-verify this — hard gate); full-text read of arXiv:2206.01579 for the exact cleaning and
   estimator recipe; confirm DESI DR1 footprint overlap and no published MeerKLASS×DESI paper.
1. Tooling + synthetic recover-a-known: a mock cube with injected correlated signal + smooth
   foregrounds round-trips through PCA cleaning and the cross-power estimator.
2. Sanity anchor: reproduce the WiggleZ-era detection pipeline behaviour on the same cube
   (their published result is the anchor for estimator + cleaning settings).
3. Real leg: cross-correlate with DESI galaxies; report significance vs PCA mode count (the
   honest sensitivity axis); null suite (shuffles/rotations) alongside.
4. GATE-2 (signal loss from over-cleaning; transfer-function assumptions; footprint edges) →
   paper, framed as an independent-tracer consistency check, not a new detection claim.

## Verification

Synthetic injected-signal recovery; the published WiggleZ-era result reproduced as the sanity
anchor before the DESI leg; null tests consistent with zero; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Foreground cleaning is the hard part** → no new cleaning method; follow the published
  pipeline verbatim and report results as a function of modes removed — conservative framing.
- **Cube may not be public** → GATE-0 is a hard stop; no plan-B dataset, park if closed.
- **Signal loss / transfer function** → quote significance with and without the published
  transfer-function correction; over-cleaning shown explicitly in the mode-count sweep.
