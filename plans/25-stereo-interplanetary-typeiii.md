# 25 — Type III to true interplanetary distances (STEREO/WAVES HFR)

Status: 🚧 in progress (tooling + real STEREO/WAVES fetch + recover-a-known done; GATE-2 + paper next)

## Context

The Wind/WAVES slice (#24) tracked a type III electron beam only to the **Alfvén surface** (~10 R⊙,
0.05 AU) because RAD2 stops at 1 MHz. **STEREO/WAVES** (the SWAVES instrument; Bougeret et al. 2008) has
a High Frequency Receiver (HFR) reaching down to **0.125 MHz**, where the emission comes from genuinely
**interplanetary** (super-Alfvénic, ≳20 R⊙) plasma — out to ~0.4 AU. This slice fits a STEREO/WAVES
type III drift across the full HFR band and inverts it, via the Leblanc heliospheric density model, to
the beam's outward speed and the true interplanetary distance it reaches — the extension the Wind/WAVES
slice flagged as out of reach.

Data are the STEREO/WAVES one-minute-averaged HFR dynamic spectra on NASA SPDF, in **ASCII** (no auth,
no CDF). The slice reuses the Wind/WAVES Leblanc tooling and the `solarbursts` dynamic-spectrum pipeline
almost wholesale; the new code is the STEREO ASCII parser.

## Reuse

- `windwaves.leblanc_density` / `leblanc_radius` / `emission_radius` / `beam_speed` (the heliospheric
  model + drift-to-speed inversion) and `solarbursts.background_subtract` / `find_burst_window` /
  `detect_burst_ridge` (the dynamic-spectrum pipeline).
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, `_figure`, AASTeX
  paper, gitignored generated artefacts).

## Deliverables

- `src/jansky_research/swaves.py`:
  - `parse_swaves_ascii` — read a STEREO/WAVES one-minute HFR ASCII file → (data, freqs MHz, times s).
  - `synthetic_ip_burst` — offline fixture (reusing the Leblanc forward model over the HFR band).
  - `fetch_swaves` (network, `# pragma: no cover`) — download the daily HFR ASCII from SPDF.
  - `run(offline=...)`, `_figure`, `_write_macros`, `_main`.
- `tests/test_swaves.py` — synthetic-fixture tests (incl. the ASCII parser on an inline sample) to the
  85% floor; no network.
- `data.py` registry note for STEREO/WAVES.

## Approach

1. **Tooling (this step).** The ASCII parser + the reused Leblanc/dynamic-spectrum pipeline, validated
   on a synthetic HFR type III whose injected speed is recovered.
2. **Real-data run (next).** A documented post-2006 STEREO/WAVES type III: track the beam from ~16 MHz
   (corona) down to ~0.125 MHz (~0.4 AU, interplanetary) and report the speed and the reach.
3. **Recover-a-known.** The beam reaches genuinely interplanetary distances (≳20 R⊙) and a type III
   beam speed; the wider distance baseline tightens the height–time fit relative to Wind/WAVES RAD2.
4. **GATE-2** before write-up — the speed/reach survive the density-model / harmonic caveats and a
   literature cross-check.
5. **Write-up** as `papers/swaves/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- `parse_swaves_ascii` parses an inline sample correctly (freqs in MHz, intensity matrix).
- Offline `run` recovers the injected beam speed.
- (Real-data, later) the STEREO/WAVES type III reaches ~0.4 AU and gives an in-range speed; GATE-2.
