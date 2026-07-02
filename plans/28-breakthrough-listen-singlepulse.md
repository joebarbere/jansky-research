# 28 — Single-pulse / pulsar dedispersion on Breakthrough Listen open data

Status: ➡️ absorbed into `plans/34-torch-fdmt-singlepulse.md` (2026-07) — its GATE 0, synthetic
discipline, and recover-the-catalogued-DM/period science are steps 0/2/4 of the torch-fdmt arc;
the GPU addendum in `survey/opportunity-scan-2026-07.md` verified PyTorch-ROCm on this machine,
which turns the "GPU optional later" note below into the slice's headline

## Context

Every slice so far has run on catalogue-level or modest space-physics data. Breakthrough Listen (BL)
open data (seti.berkeley.edu/opendata) is the best-quality **bulk time-domain** radio data an amateur
can get: professionally-calibrated GBT / Parkes / Green Bank / MeerKAT voltage-derived filterbanks for
thousands of targets, public, no auth. The headline use is SETI, but the data-source scan flagged its
**non-SETI secondary science** (pulsar single-pulse, scintillation, RFI-vs-astrophysical) as HIGH —
and it is the natural home for the one classic radio method the repo still lacks: **incoherent
dedispersion and single-pulse / period search of a broadband dispersed signal**.

This is *not* the existing `driftsearch` slice (SETI narrowband Doppler-drift tones, `jansky.seti`). A
dispersed pulse sweeps across the *whole* band with a quadratic ν-dependent delay; recovering it needs
DM trials, dedispersion, and a boxcar matched filter (single pulses) or epoch folding (a pulsar's
period) — a different physics and a different toolchain. `jansky.transients` already implements all of
it, untouched by the research repo so far.

**The gating risk is data size.** BL products are GB–TB. The slice is viable only if a *small*,
self-contained product is used — one decimated/sub-banded filterbank of a known bright pulsar — with a
mandatory synthetic offline fixture so tests/CI never touch the network or a large file. The plan's
**GATE 0** is to confirm such a file exists and downloads in seconds–minutes before any tooling is
built.

## Deliverables

- `src/jansky_research/singlepulse.py` — pure-NumPy tooling, composing `jansky.transients`:
  - reuse `transients.dispersion_delay` / `disperse_pulse` (forward model) / `dedisperse` /
    `dm_search` (`DMSearchResult`) / `boxcar_snr` / `fold_profile` / `epoch_folding_search`.
  - reuse `jansky.formats.read_filterbank` (SIGPROC `.fil`) / `read_guppi` / `Spectrogram` to ingest
    the BL product, and `jansky.rfi` to flag the worst channels/times before dedispersion.
  - `synthetic_observation` — an offline dynamic spectrum with an injected dispersed pulse train at a
    known DM and period (built from `disperse_pulse`), so `dm_search` + folding round-trip.
  - `fetch_bl_filterbank` (network, `# pragma: no cover`) — download one small BL `.fil`/`.h5`.
  - `run(offline=...)`, `_figure` (DM–time S/N "butterfly" + folded profile), `_write_macros`, `_main`.
- `tests/test_singlepulse.py` — synthetic-fixture tests to the 85% floor; no network, no large files.
- `data.py` registry entry for the chosen BL product; `papers/singlepulse/` (AASTeX);
  `survey/singlepulse-findings.md`. Optional `blimpy` extra (like `windwaves=cdflib`) only if the
  chosen product is HDF5 rather than SIGPROC `.fil` (jansky's `read_filterbank` covers `.fil` already).

## Approach

0. **GATE 0 — data reality check (do FIRST, before tooling).** Find one BL open-data product that is
   (a) a *known bright pulsar* (e.g. B0329+54, B0531+21/Crab, B0833−45/Vela, or a BL pulsar-survey
   pointing), (b) small enough to download in seconds–minutes (a single coarse-channel or
   time/freq-decimated filterbank; use `blimpy`'s sub-band/decimate if needed), and (c) readable by
   `jansky.formats.read_filterbank`. If no compact public file is found, STOP and report — do not build
   on a TB file. Record the exact URL + size in `data.py`.
1. **Tooling + synthetic recover-a-known.** Inject a dispersed pulse train at known (DM, period, width)
   into a synthetic dynamic spectrum via `disperse_pulse`; confirm `dm_search` recovers the DM (peak of
   the DM–S/N curve), `boxcar_snr` recovers the single-pulse width, and `epoch_folding_search` recovers
   the period. Offline `run` round-trips all three within tolerance.
2. **Real recover-a-known.** On the chosen BL pulsar file: RFI-flag, dedisperse over a DM grid, detect
   the pulsar's single pulses (or fold), and **recover its catalogued DM and period** (ATNF values) —
   the validation. Report the DM–time butterfly peak and the folded profile S/N.
3. **GATE-2 science review** (science-reviewer) before the write-up: the recovered DM/period must match
   the catalogue within errors; surface the honest caveats (single pointing, RFI environment, the BL
   power scale is not flux-calibrated, scattering/scintillation broadening, the boxcar S/N depends on
   the noise model).
4. **Write-up** `papers/singlepulse/` — a reproducibility/method contribution and a recover-a-known
   (a known pulsar's DM+period pulled back out of public BL data), not a discovery.

## Verification

- `make test` / `cov` green on the synthetic fixture only (no network, no large file), 85% floor;
  `ruff` + `mypy` clean.
- Offline `run` recovers the injected DM, single-pulse width, and period within tolerance.
- (Real-data) the recovered DM and period of the known pulsar match the ATNF catalogue; GATE-2 sign-off.
- `make reproduce` fetches the one small BL product and reproduces the metrics/figure/macros.

## Risks & mitigations

- **File size (highest) →** GATE 0 hard stop; one decimated/sub-banded product only; synthetic fixture
  carries all tests; `data/*` gitignored; the real file off the default test path (jansky's HI4PI-style
  large-file handling).
- **RFI at GBT/Parkes →** `jansky.rfi` flagging before dedispersion; pick a clean band/pointing.
- **No compact public pulsar file →** fall back to the jansky-vendored `your_28.fil` as the *offline*
  fixture substrate and document the real-data leg as conditional; or pick a different bright pulsar.
- **GPU not required →** dedispersion is GPU-friendly but the slice stays CPU on the small/synthetic
  data (85% floor); a GPU dedisperser is an optional later optimisation, not a dependency.
- **Overlap with `driftsearch` →** none in method (broadband dispersed pulses + folding vs narrowband
  drift); they share only the injection-recovery philosophy.
