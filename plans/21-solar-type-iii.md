# 21 — Solar type III radio bursts: drift rate → exciter (electron-beam) speed (e-Callisto)

Status: ✅ done (tooling + real e-Callisto fetch + recover-a-known + GATE-2 + paper `papers/solarbursts/`)

## Context

A **type III** solar radio burst is the signature of a beam of ~keV electrons streaming out along
open magnetic field from a flare. As the beam climbs through the corona it excites Langmuir waves at
the local plasma frequency $f_p \propto \sqrt{n_e}$, which convert to radio at $f_p$ (fundamental) or
$2f_p$ (harmonic). Because the density — and so $f_p$ — falls with height, the burst drifts **fast**
from high to low frequency: the drift rate encodes how quickly the exciter climbs, and with a coronal
density model the frequency drift becomes a **height-versus-time track** whose slope is the beam speed
(famously a sizeable fraction of $c$, ~0.1–0.5 c).

The **e-Callisto** network (i4ds/FHNW; Benz et al. 2009) publishes ground-based 15-minute dynamic
spectra (~20–900 MHz) as gzipped FITS, fully open over HTTP with **no authentication** — an ideal,
reproducible, CASDA-free data source. This slice builds a small tested tool that fits a type III
burst's drift in an e-Callisto spectrum and converts it to an exciter speed, reusing the course's
`jansky.solar` coronal-physics helpers. The deliverable is a tool + a recover-a-known (the canonical
~0.1–0.5 c beam speed) with the density-model systematic reported honestly — not a discovery.

## Reuse

- `jansky.solar.density_from_plasma_frequency` and `newkirk_radius` (Newkirk 1961 coronal model) to map
  observed frequency → plasma frequency → density → heliocentric radius; `plasma_frequency`,
  `newkirk_density`, `R_SUN_KM` for the synthetic forward model.
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, `_figure`, AASTeX
  paper, gitignored generated artefacts), as in the `vlbi`/`stacking` slices.

## Deliverables

- `src/jansky_research/solarbursts.py`:
  - `synthetic_burst` — offline fixture: a dynamic spectrum with an injected type III ridge of known
    exciter speed (built by the *same* `jansky.solar` forward model the analysis inverts) + noise.
  - `background_subtract` — per-channel baseline removal (e-Callisto raw has strong per-channel offset).
  - `detect_burst_ridge` — the frequency of peak intensity per time sample, thresholded → the (t, f)
    drift ridge.
  - `fit_drift_rate` — a representative df/dt (MHz/s) from the ridge.
  - `exciter_speed` — map the ridge frequencies to coronal radii (harmonic + Newkirk fold as params),
    fit radius vs time, return the beam speed in km/s and in units of c.
  - `fetch_ecallisto` (network, `# pragma: no cover`) — fetch + gunzip + parse an e-Callisto FITS to
    (data, freqs MHz, times s) via the public archive directory listing.
  - `run(offline=...)`, `_figure` (dynamic spectrum + ridge + height–time), `_write_macros`, `_main`.
- `tests/test_solarbursts.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry note for the e-Callisto archive.

## Approach

1. **Tooling (this step).** Pure-NumPy/SciPy drift fit + the density-model inversion composing
   `jansky.solar`, validated on a synthetic type III whose injected exciter speed is recovered (the
   forward fixture and the inverse analysis share the Newkirk mapping, so a clean burst round-trips).
2. **Real-data fetch (wired + validated end-to-end; clean recover-a-known next).** `run(offline=False)`
   fetches an e-Callisto file, background-subtracts, finds the burst window, detects the ridge, and
   reports the drift and exciter speed. Validated functional on the 2011-08-09 X6.9 flare (BIR, 20-90
   MHz): the pipeline downloads + parses the FITS, locates the burst at the flare time (UT 08:04.8),
   and detects a ~140-channel ridge. **Caveat surfaced by that run:** the 08:05 interval is a type III
   *storm*, so a fixed window smears many bursts into an artificially slow drift -- the clean
   recover-a-known needs an *isolated* type III event and RFI-aware single-burst isolation. Honest
   systematics to report: fundamental-vs-harmonic (a factor of 2 in density → height), the Newkirk fold
   factor (1× quiet → up to ~4× over active regions), and projection.
3. **Recover-a-known (done).** Picked clean, isolated, quality-1 type III bursts from the Monstein SGD
   burst lists (a flare interval is a *storm* the fit R² correctly flags at <0.3). The 2011-09-14 11:50
   UT BIR burst is clean (R² = 0.90, drift −3.3 MHz/s, 10–79 MHz) and gives an exciter speed of
   0.09 c (fundamental) → 0.14 c (harmonic) → 0.27 c (harmonic, 4× Newkirk) — squarely in the canonical
   0.1–0.5 c range. All four isolated test events land in-band. The R² metric, the robust sigma-clipped
   fit, and the `--recover` reproducible event were added in this leg. See `survey/solarbursts-findings.md`.
4. **GATE-2 (done).** science-reviewer verdict PASS-WITH-FIXES. Fixed: Zhang et al. 2018 is A&A 618
   **A165** (not A33); the synthetic fixture is a round-trip *algebraic consistency* check, not an
   "independent" validation; dropped the hand-wavy "true speed ~0.2 c de-biasing" (the peak-time
   0.137 c is already consistent with peak-time LOFAR 0.17 c); acknowledged the **peak-time vs onset**
   bias (~15–30% low); anchored the claim on the one clean event (R²=0.90) rather than "all four"
   (events 1–2 have R²<0.55); noted the slow drift (~25th percentile), the Newkirk 1–3 R⊙ range, the
   uncalibrated single-station data, and the ~10 MHz ionospheric cutoff.
5. **Write-up (done)** as `papers/solarbursts/main.tex` + `refs.bib` — compiles in both the CI
   offline-synthetic macro build and the real `make reproduce` build; every number `\input` from macros.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers the injected exciter speed to within a few percent for a clean burst.
- `exciter_speed` round-trips the `jansky.solar` forward model (a known ridge → its injected speed).
- (Real-data, later) the recovered speed is in the canonical type III range; GATE-2 sign-off.
