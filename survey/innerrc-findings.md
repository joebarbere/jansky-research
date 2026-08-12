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

## Referee round (2026-08-12) — verdict: major revision, 19 findings

Two blockers, both in the anchor leg, both invisible to the test suite:

1. **`\irRhoMax` = 0.32 was a boundary artifact.** `decompose_rc` bounds `v_bulge` at
   800 km/s; the variant supplying the quoted maximum (`nfw:R>2kpc:w`) sat at **exactly**
   800.0. `sensitivity_scan` counted it as converged because `curve_fit` did not raise.
   Mechanism: excising `R < 2 kpc` removes every point that constrains the bulge, so
   `(v_bulge, a_bulge)` slide along the Plummer degeneracy ridge until they hit walls.
   Adding a bound-contact check found it is worse than the referee could see from the
   committed subset — **6 of 8 variants are railed**, three of them on `h_halo`'s lower
   bound, which was invisible because `sensitivity_variants` filtered `v_halo`/`h_halo`
   *out of the committed evidence* (finding 5: the only two numbers `rho_dm_gev` is computed
   from were not committed). Interior range: **0.20–0.24**, from 2 variants.
2. **"Spans the full range from their lower limit to the consensus density"** asserted an
   interval [0.107, 0.3] that no committed number supports — the scan never reaches 0.107.

The fix moved the compatibility claim onto the primary fit's own covariance, which was
already committed and unused: ρ_DM = 0.24, **1σ 0.16–0.31**, which reaches the consensus
honestly and cannot be manufactured by a bound. Note the corner calculation is not a matter
of pairing the extremes — ρ is *not* monotonic in `h_halo` at R₀ (it rises with the NFW
scale radius), so the naive `(v+dv, h−dh)` pairing gives 0.24 where the true 1σ maximum is
0.31. All four corners are scanned.

**The result got stronger, not weaker.** Comparing by rms (13.6 vs 11.8) made the two
solutions look like alternative corners of one degeneracy. In χ²/N they are not (3.18 vs
1.92), and the difference is not spread over the curve: beyond 8 kpc — exactly where
ρ_DM(R₀) is set — their published halo sits **+1.33σ per point** below their own unified
curve over 32 independent points. That is a one-sided bias, not scatter. The paper now says
what it actually shows: *their published curve prefers more local dark matter than their
published decomposition reports.*

Three further overclaims, all fixed against committed evidence:

- **"the bar-region inner peak reproduce fully"** — it does not. Inside 2 kpc this HI-only
  curve runs **−36 km/s** low over 28 bins, and their 255 km/s peak at 550 pc appears here as
  214 km/s at 1450 pc. Their inner curve is CO-dominated; this is a limitation of the
  replication, now stated as one.
- **"replicates in period and phase"** — the period is **36% longer**, the damping length
  rails against its bound, the fit is *seeded at their published period and phase* (so the
  phase agreement is partly anchored by construction), and the sinusoid accounts for under
  half the variance (5.9 vs 8.0 km/s). The findings file already said "qualitatively"; the
  abstract had upgraded it. Reverted to what was measured.
- **"the calibration route confirms it independently"** — it is not independent. V is linear
  in σ_v and both estimators are calibrated on the same sightlines with the same weights, so
  the difference in required dispersions is *algebraically the same statistic* restricted to
  the solar-circle window. Also "a difference equal to the direct bias" (17.17 vs 17.94) is
  agreement within spread, not equality.

Also: the Table-2 agreement is now stated with its zero point (a one-parameter calibration
pinned on 7.0–8.15 kpc; the uncalibrated median offset is −14.5 km/s, already committed and
previously unmacroed), "unconstrained refit" → "bounded" with the bounds stated, "reproduces
exactly" → "to the precision they quote", and `chemin2015`'s title corrected against Crossref.

`tests/test_innerrc.py` now fails if the primary refit has any parameter on a bound, if the
quoted range includes a railed variant, or if a committed variant drops the halo parameters
its ρ is computed from.
