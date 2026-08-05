# innerrc — findings (in progress)

Slice plan: `plans/86-inner-rotation-curve.md`. Target: Sofue & Kohno 2025 (PASJ 77, 1335;
arXiv:2509.23581), the modern TVM inner rotation curve of the Milky Way.

## Increment 1 (tooling + synthetic recover-a-known + offline anchor) — done 2026-08-04

`src/jansky_research/innerrc.py` + `tests/test_innerrc.py` (10 tests, offline): table parser
over the vendored arXiv-source RC tables (`tests/data/sofue2025/`, 161 inner + 103 unified
rows, R to 25.6 kpc), Gaussian-decomposition TVM (greedy seed + joint multi-Gaussian
refinement; recovers injected terminal velocities to ±6 km/s on crowded synthetic spectra),
the `hi`-slice threshold estimator kept verbatim for the head-to-head (it reads several km/s
high on the same synthetic spectra — the documented bias, now reproduced), the paper's
dispersion correction and Gaussian-weighted binning (eq. 10–15), the BH + Plummer + halo
decomposition (eq. 19–27; NFW and Burkert), and the E/W damped-sinusoid fit (eq. 29).

### Anchor result — their arithmetic validates; their ρ_DM is one corner of a broad degeneracy

Decomposing **their own published unified RC table** (their data, their model family, their
constants):

| solution | rms vs their table (km/s) | ρ_DM(R₀) (GeV/cm³) |
|---|---|---|
| their Table 1 parameters | 13.6 | **0.107** ✓ (reproduces their published value exactly) |
| our unconstrained least-squares refit | 11.8 | **0.244** |
| 8-variant sensitivity scan (NFW/Burkert × fit range × weighting) | — | **0.19–0.32** |

Two conclusions, both honest:

1. **Their arithmetic checks out.** Evaluating their Table-1 parameters through their own
   ρ₀ = V_h²/(G h²) convention (their eq. 25 absorbs 4π into g(x); V_h converts to our
   convention by √4π — a unit-convention trap now covered by a test) reproduces
   0.107 GeV/cm³ on the nose, and their solution describes their curve at 13.6 km/s rms.
2. **The 0.107 is not demanded by their curve.** An unconstrained refit of the same table
   fits marginally better (11.8 km/s) with a more concentrated, halo-heavier solution at
   ρ_DM = 0.24, and every converging variant in the sensitivity scan lands at 0.19–0.32 —
   the consensus range. The bulge is stable across all solutions (ours: 394 km/s / 316 pc vs
   their 406 / 333); the degeneracy lives entirely in the disc–halo split. This *quantifies*
   the paper's own "lower limit" caveat: the published curve is fully compatible with
   ρ_DM ≈ 0.3, and the factor-3 gap to consensus is a fit-degeneracy width, not a tension.

Caveats: our refit does not replicate their exact fitting procedure (weighting/algorithm
unpublished beyond "least-squares"); the scan's weighted variants use their δV column. None
of the eight variants reproduces their corner, suggesting an additional constraint in their
fit (fixed scale radius or similar) not stated in the paper — a question for the paper draft
to raise neutrally. Committed evidence: `results/innerrc_anchor.json`.

## Next increments

- HI4PI leg (fetch |b|≈0 tiles from CDS, ~15–20 GB): build the LV diagram, run both TVM
  estimators on the real spectra (the `hi`-bias head-to-head with real numbers), reproduce
  their Table 2 inner RC from raw survey data.
- E/W asymmetry from the two quadrants; sawtooth residual check.
- GATE-2, paper.
