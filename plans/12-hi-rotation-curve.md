# Plan 12 — Milky Way HI rotation curve (tangent point) ✅

> Context: fifth slice, a "recover a known result" pick. Scope: small. Status: validation.

## Context

The flat rotation curve of the Galaxy is the textbook dark-matter signature. The tangent-point
method extracts it from inner-Galaxy HI 21 cm terminal velocities — a clean, offline-capable
validation target with public data (the LAB survey).

## Deliverables

- `src/jansky_research/hi.py` — `read_lab_slice` (parse a LAB $(b,v)$ FITS), `terminal_velocity`,
  `tangent_point` ($R=R_0\sin\ell$, $V=v_\mathrm{term}+V_0\sin\ell$), `rotation_curve`,
  `fetch_lab_longitude` (VizieR, cached), `run()` (curve + figure + metrics), and
  `synthetic_lv_slice` (offline fixture with a known flat curve). Tested to the 85% floor.
- A real run over $\ell=10°$–$80°$ → `results/rotation_curve.json` + figure.
- `survey/hi-findings.md` — the flat-curve result + honest caveats.

## Approach

Reuse `jansky.data`'s LAB slice machinery; Reid et al. (2019) $R_0=8.15$ kpc, $V_0=236$ km/s.
**Caveats:** the fixed-$T_B$ terminal-velocity estimator biases the absolute level slightly high
(threshold-dependent); non-circular motions add scatter; first quadrant / $b=0$ / inner Galaxy only.
The *flatness* is the robust, threshold-independent result.

## Verification

- `make cov` ≥85% on synthetic fixtures (offline run recovers the injected flat curve).
- Real run yields a flat $V(R)$ (no Keplerian decline) over $R\approx1.4$–$8$ kpc — the sanity check.
- **science-reviewer** gate on the method, the normalisation caveat, and the citations.
