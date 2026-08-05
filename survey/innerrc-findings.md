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

## Increment 2 (HI4PI real-data leg + paper-style figures) — done 2026-08-05

All 10 galactic-plane CAR tiles from CDS (2.6 GB — far under the plan's 15–20 GB budget),
|b|<3° LV diagrams, 1,113 sightlines through **both** terminal-velocity estimators
(ℓ = 5–89.5° and 270.5–355°; two bugs found and fixed en route: a |sin ℓ| sign error that
emptied the fourth quadrant, and a CDS transfer stall now absorbed by a timeout+resume+retry
fetcher). Committed evidence: `results/innerrc_hi4pi.json`; figures in
`papers/innerrc/figures/` (LV diagram + fitted envelope, single-spectrum decomposition at
ℓ≈30°, RC vs their Table 2, E/W + eq.-29 fit, anchor decomposition — mirroring their
Figs. 1/2/8/14/11).

- **The `hi` caveat is closed with a measured number.** On 1,113 real HI4PI spectra the
  threshold estimator reads **+17.9 km/s (median)** high of the Gaussian-decomposition
  terminal velocity (p16–p84: +9 to +23) — the same sign as, and larger than, the ~7 km/s
  the literature quotes for threshold-vs-fitting. Independently confirmed by the paper's own
  calibration route: iterating each estimator's dispersion correction to meet Θ₀ at R₀ gives
  σ_v = 9.5 km/s (Gaussian) vs 26.6 km/s (threshold) — difference 17.2 ✓. Notably the
  paper's adopted σ_HI = 15 km/s sits *between* our two estimators, suggesting their
  decomposition lands between component-centre and envelope-crossing in practice.
- **Their Table 2 reproduces at the ~4% level from raw survey data**: mean |ΔV| = 8.8 km/s
  (median −8.9) over 123 bins at R > 2 kpc, after per-estimator σ calibration (the same comparison at the first-attempt fixed σ=15 gives −14.5, committed as
  `median_dv_kms_fixed_sigma15`; the −13.4 quoted in an earlier draft came from the
  pre-Q4-fix run and is superseded). The residual is a shape difference, candidates:
  their HI+CO merge (CO traces colder gas), finer |b| windows per survey, and their 2-pc
  Gaussian binning; shape (flat, non-Keplerian, bar-region peak) is fully consistent.
- **The E/W asymmetry replicates qualitatively**: fitting their eq.-29 damped sinusoid in the
  mid-disc window (2–8 kpc; inside 2 kpc the bar rails any smooth fit — reported, not
  hidden) gives amplitude 9.5 km/s, period 6.0 kpc, phase 3.8 kpc vs their directly stated ≃15 km/s mid-disc
  amplitude and eq.-29 period/phase 4.4/4.0 kpc — same phenomenon, same phase, softer amplitude
  from HI alone; damping unconstrained over our window (railed, reported). Mid-disc δV rms
  8.0 km/s.

## Next increments



- HI4PI leg (fetch |b|≈0 tiles from CDS, ~15–20 GB): build the LV diagram, run both TVM
  estimators on the real spectra (the `hi`-bias head-to-head with real numbers), reproduce
  their Table 2 inner RC from raw survey data.
- E/W asymmetry from the two quadrants; sawtooth residual check.
- GATE-2, paper.
