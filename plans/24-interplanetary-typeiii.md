# 24 — Interplanetary type III bursts: tracking the beam into the heliosphere (Wind/WAVES)

Status: ✅ done (tooling + real Wind/WAVES fetch + recover-a-known + GATE-2 + paper `papers/windwaves/`)

## Context

A type III electron beam does not stop at the top of the corona — it streams out along the open field
into interplanetary space, exciting radio emission at the falling plasma frequency all the way toward
1 AU. Space-based receivers (no ionospheric cutoff) see this as a slow drift from a few MHz down to
tens of kHz over tens of minutes. This slice fits that drift in a **Wind/WAVES** dynamic spectrum and
inverts it via a **heliospheric** density model (Leblanc, Dulk & Bougeret 1998) to the beam's outward
speed and the heliocentric distance it reaches — the interplanetary companion to the coronal
`solarbursts` slice (#21), which used the Newkirk corona valid only to a few solar radii.

Data are the Wind/WAVES Level-2 radio CDFs on NASA SPDF (public, no auth). The slice reuses
`solarbursts`' dynamic-spectrum tools and `jansky.solar`, adding the Leblanc model and the CDF fetch.

## Reuse

- `solarbursts.background_subtract` / `find_burst_window` / `detect_burst_ridge` / `_robust_linfit`
  (the dynamic-spectrum + robust-fit pipeline) and `jansky.solar.density_from_plasma_frequency`.
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, `_figure`, AASTeX
  paper, gitignored generated artefacts).

## Deliverables

- `src/jansky_research/windwaves.py`:
  - `leblanc_density` / `leblanc_radius` — the heliospheric density model and its numerical inverse.
  - `emission_radius` — frequency → plasma frequency → density → heliocentric radius.
  - `beam_speed` — fit radius vs time → the outward speed (km/s, c), the radius range (R⊙, AU), and R².
  - `synthetic_ip_burst` — offline fixture from the same Leblanc forward model (round-trips).
  - `fetch_windwaves` (network, `# pragma: no cover`; needs the `windwaves` extra / `cdflib`).
  - `run(offline=...)`, `_figure` (dynamic spectrum + beam track), `_write_macros`, `_main`.
- `tests/test_windwaves.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry note for Wind/WAVES; the `windwaves` (`cdflib`) optional extra.

## Approach

1. **Tooling (done).** Leblanc model + inversion composing `solarbursts`/`jansky.solar`, validated on a
   synthetic IP type III whose injected speed is recovered (within 15%).
2. **Real-data run (done).** The 2003-10-28 X17 flare IP type III (Wind/WAVES RAD2): 256 channels,
   1.075–13.825 MHz; the Leblanc inversion places the emission at 2.4–10.2 R⊙ (harmonic) — interplanetary,
   beyond the corona — and the outward radial speed is 0.045 c (fundamental) / 0.083 c (harmonic),
   reaching ~0.05 AU in RAD2 (RAD1 extends toward 1 AU).
3. **Recover-a-known.** The beam track and a speed in the (low end of the) type III range.
4. **GATE-2** before write-up — the speed survives the density-model / harmonic / projection caveats and
   the moderate fit, and a literature cross-check.
5. **Write-up** as `papers/windwaves/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers the injected beam speed; `leblanc_radius` round-trips `leblanc_density`.
- (Real-data, later) the Wind/WAVES IP type III gives a heliocentric track and an in-range speed; GATE-2.
