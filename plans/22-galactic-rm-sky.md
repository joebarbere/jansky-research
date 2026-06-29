# 22 — The Galactic Faraday rotation sky (Taylor+2009 RM catalogue)

Status: ✅ done (tooling + real Taylor+2009 fetch + recover-a-known + GATE-2 + paper `papers/rmsky/`)

## Context

The plane of polarisation of a radio source rotates as it passes through the magnetised, ionised
interstellar medium: :math:`\chi(\lambda) = \chi_0 + \mathrm{RM}\,\lambda^2`, with the **rotation
measure** :math:`\mathrm{RM} = 0.81\int n_e B_\parallel\,\mathrm{d}l` (rad m⁻²) integrating the
line-of-sight magnetic field weighted by electron density. The RMs of tens of thousands of
*extragalactic* sources therefore form a screen-map of the **Galactic magnetic field** along every
sightline. Two large-scale signatures are textbook (Taylor, Stil & Sunstrum 2009; Oppermann et al.
2012): (i) :math:`|\mathrm{RM}|` is strongly **enhanced toward the Galactic plane** (the path length
through the disk grows as :math:`\sim\csc|b|`); and (ii) the RM sky is **sign-organised** — in the inner
Galaxy, RMs are predominantly positive above the plane and negative below it (the antisymmetric
disk/halo field).

This slice reproduces both from the **Taylor+2009** NVSS RM catalogue (37 543 sources; VizieR
`J/ApJ/702/1230`, fully public, no auth), reusing the course's `jansky.polarization` helpers for the
underlying :math:`\lambda^2` measurement. The deliverable is a tested tool + a recover-a-known (the
plane enhancement and the sign asymmetry), with the catalogue's known limitations reported honestly —
not a new astrophysical claim.

## Reuse

- `jansky.polarization.faraday_rotate` and `rotation_measure_fit` for the foundational
  "how an RM is measured from :math:`\chi(\lambda^2)`" step (offline-testable round-trip).
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, `_figure`, AASTeX
  paper, gitignored generated artefacts), as in the `vlbi` / `solarbursts` slices.

## Deliverables

- `src/jansky_research/rmsky.py`:
  - `rm_from_angles` — wrap `rotation_measure_fit`: recover RM from polarisation angle vs wavelength.
  - `latitude_profile` — median :math:`|\mathrm{RM}|` in :math:`|b|` bins (the disk enhancement).
  - `enhancement_ratio` — median :math:`|\mathrm{RM}|` (plane) / (poles).
  - `sign_asymmetry` — mean RM in the four (north/south × inner/outer) quadrants (the antisymmetry).
  - `synthetic_rm_sky` — offline fixture: a toy disk (:math:`\csc|b|`) + sign-organised field +
    extragalactic scatter, so the tests recover both signatures.
  - `fetch_taylor2009` (network, `# pragma: no cover`) — VizieR fetch + Galactic coordinates.
  - `run(offline=...)`, `_figure` (Galactic RM sky + the :math:`|b|` profile), `_write_macros`, `_main`.
- `tests/test_rmsky.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry note for the Taylor+2009 RM catalogue.

## Approach

1. **Tooling (this step).** Pure-NumPy profile/asymmetry statistics + the `jansky.polarization`
   reuse, validated on a synthetic RM sky whose injected disk enhancement and sign pattern are
   recovered.
2. **Real-data run (next).** `run(offline=False)` fetches Taylor+2009 from VizieR, computes Galactic
   coordinates, and reports the enhancement ratio and the quadrant signs. A live probe already
   confirms the signal: median :math:`|\mathrm{RM}|` falls 62 → 11.5 rad m⁻² from :math:`|b|<10°` to
   :math:`|b|>60°` (ratio ~5), and the inner Galaxy shows mean RM +9 (north) / −24 (south).
3. **Recover-a-known.** The plane enhancement and the inner-Galaxy sign antisymmetry reproduce the
   established Galactic-RM-sky result (Taylor+2009).
4. **GATE-2** before write-up — the signatures survive the catalogue's :math:`n\pi`-ambiguity and
   intrinsic-RM caveats and a literature cross-check.
5. **Write-up** as `papers/rmsky/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers the injected plane enhancement (ratio > 3) and the correct quadrant signs.
- `rm_from_angles` round-trips `faraday_rotate` (a known RM → recovered).
- (Real-data, later) the Taylor+2009 enhancement and sign asymmetry match the published Galactic RM
  sky; GATE-2 sign-off.
