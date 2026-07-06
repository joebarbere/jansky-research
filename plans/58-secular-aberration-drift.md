# 58 — Secular aberration drift from the public ICRF3×Gaia DR3 cross-match

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm the public
ICRF3 and Gaia DR3 AGN cross-match tables and check no post-2503.03389 reproduction exists

## Context

The solar system's Galactocentric acceleration shows up as a secular aberration drift — a ~5 µas/yr
glide in AGN proper motions toward the Galactic centre, measured by Gaia (Klioner+2021) with
refinements in arXiv:2503.03389. Both underlying tables are public (ICRF3 via IERS/CDS; Gaia DR3
AGN astrometry via the archive), yet no independent-pipeline reproducibility note exists from the
public cross-match alone. This slice is the cleanest quick win in fable-ideas.md (F21): a VSH/glide
fit on the cross-matched proper motions, extending the merged `offsets` slice (radio–optical
cross-match + offset statistics, ~90% reuse) with near-zero new data access. Recover-a-known *is*
the result: an independent reproduction of the glide amplitude and direction, with the Gaia AGN
proper-motion noise floor quoted honestly. Gaia DR4 (2026-12) makes the same pipeline ~10× sharper
— this slice builds and validates it now.

## Deliverables

- `src/jansky_research/aberration.py`: `fetch_icrf3` + `fetch_gaia_agn_pm` (archive queries,
  `# pragma`), `crossmatch_icrf3_gaia` (reuse the `offsets` matcher), `vsh_glide_fit` (degree-1
  vector-spherical-harmonic / pure-glide fit with bootstrap errors), `synthetic_glide_catalogue`
  (injected glide on the real source distribution + Gaia PM noise → recovery),
  `compare_published` (Klioner+2021 + arXiv:2503.03389 table), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/aberration/`; `survey/aberration-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on Klioner+2021 and arXiv:2503.03389 plus an ADS search for any
   published independent reproduction from the public cross-match; verify the ICRF3 and Gaia DR3
   AGN table URLs and column contents.
1. Tooling + synthetic recover-a-known: inject a known glide (amplitude + apex) into a mock
   catalogue with the real sky distribution and per-source PM errors; the VSH fit must recover
   both within stated errors. Pure NumPy/astropy, CPU, days.
2. Real leg: cross-match ICRF3×Gaia DR3, quality cuts (RUWE, PM error), glide fit in mean PM;
   compare amplitude and apex direction to the published values; quote any discrepancy against
   the AGN PM noise floor rather than explaining it away.
3. GATE-2 science review: source-selection sensitivity, correlated Gaia systematics (scanning
   law), the honest scope of "reproduction note, not new measurement".
4. Paper: independent-pipeline reproducibility note; DR4 upgrade path stated as the follow-on.

## Verification

Synthetic round-trip recovers injected glide amplitude and apex; the real-leg anchor is the
~5 µas/yr Klioner+2021 glide itself — recovering it IS the result; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Gaia AGN PM noise floors / correlated systematics** → quote any amplitude or apex discrepancy
  honestly; report per-cut stability instead of tuning to the published value.
- **Modest ceiling (reproduction note)** → frame as pipeline validation that becomes ~10× sharper
  when Gaia DR4 lands (2026-12); the tooling is the durable deliverable.
