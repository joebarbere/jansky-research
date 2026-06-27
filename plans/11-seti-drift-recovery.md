# Plan 11 — CPU-only SETI drift-search injection-recovery 🚀

> Context: fourth backlog slice. Scope: small. Status: validation + benchmark.

## Context

SETI surveys each quote an injection-recovery efficiency, but on their own data with their own
grids — there is no shared, reproducible, CPU-only benchmark (survey gap). This slice builds one on
the pure-NumPy Doppler-drift search in `jansky.seti`, and validates the detector on a real known
signal (the Voyager-1 carrier). Honest expectation: a reproducible sensitivity benchmark + a
real-signal validation, not a technosignature.

## Deliverables

- `src/jansky_research/driftsearch.py` — `injection_recovery` ($P_\mathrm{detect}$ over SNR×drift),
  `false_positive_rate`, `completeness_snr`, `run()` (metrics + recovery heatmap), and
  `validate_voyager()` (recover the Voyager-1 carrier; optional `voyager` extra = h5py + hdf5plugin).
  Built on `jansky.seti.drifting_tone`/`drift_search`. Tested to the 85% floor (offline).
- The Voyager-1 file registered in `data.py` (`voyager1-h5`, ~50 MB, opt-in `large`).
- `survey/drift-findings.md` — benchmark result + Voyager validation + honest caveats.

## Approach

Inject synthetic tones, search, measure recovery; calibrate the threshold against the noise-only
false-positive rate. **Caveats:** the SNR/drift units are the `jansky.seti` model's internal units
(a relative benchmark, not on-sky janskys/Hz s$^{-1}$); no cadence/RFI realism.

## Verification

- `make cov` ≥85% on synthetic fixtures (offline; the benchmark is fully reproducible).
- Real validation: `validate_voyager` recovers the Voyager-1 carrier (S/N $\gg$ blank) as the sanity
  check.
- **science-reviewer** gate on the statistic, threshold calibration, units honesty, and the Voyager
  frequency/Doppler claim.
