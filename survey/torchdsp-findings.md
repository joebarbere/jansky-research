# Findings — torch-dsp: the pure-PyTorch coherent-DSP suite (plan 43)

`jansky_research.torchdsp` extends the merged `torch-fdmt` arc with three kernels no
pure-PyTorch (or JAX) implementation of which existed: coherent dedispersion, SumThreshold +
spectral-kurtosis RFI excision, and a radix-2 FFA. One `device=` argument covers CPU, CUDA, and
ROCm — the canonical real run below executed end-to-end on the RX 7600 XT (gfx1102).

## GATE 0 (2026-07-06, per-kernel repo/full-text sweep)

- **Coherent dedispersion — OPEN for torch/JAX.** dspsr (C++/CUDA), CDMT (98.7% CUDA),
  PyTorchDedispersion (incoherent only — confirmed unchanged). Fences requiring careful
  wording: CoherentDedispersion.jl (Julia, active) and a SYCL prototype (fxzjshm) exist → the
  claim is "pure-PyTorch / pip-installable", NOT "first device-portable".
- **RFI kernels — OPEN.** `jess` (Kania+2026, AJ 171, 73 — NOT Agarwal as fable-ideas had it)
  is CuPy/CUDA-locked; IQRM/AOFlagger CPU. No torch SumThreshold/SK anywhere.
- **FFA — OPEN for torch, fence appeared.** riptide is C++/CPU (GPU issue closed unimplemented
  2024-03-27: "let's face it, it's not happening"); **`gaffa`**, a CUDA FFA scaffold, appeared
  on GitHub 2026-06-12 (0 stars, build-only README) → claim "no pure-PyTorch FFA", cite gaffa
  as concurrent CUDA work, and watch it.
- **Data gate verified**: smallest CHIME baseband file (`FRB20181231C_24366209_beamformed.h5`,
  150.8 MB) fetched anonymously from CANFAR (DOI 10.11570/23.0029; release paper ApJ 969, 145).
  Format: tied-beam complex64 voltages (108 chan × 2 pol × 57,490 × 2.56 µs), per-channel
  `time0` alignment, `good_channels` = GLOBAL channel ids (must be mapped via
  `index_map/freq['id']` — a real loader bug the first run caught).

## Deltas from plan 43

- **FFA oracle = brute-force folding, not riptide**: an exact in-repo oracle (same S/N metric)
  beats an external C++ build dependency; riptide is cited, not depended on.
- **Crab period re-find: honest null.** The vendored 2.1-s Parkes file cannot support it — the
  brute fold at the published 33.7 ms gives S/N 2.2 (nothing to find; the file was vendored
  for giant-pulse work, and the regular pulse is too weak at this length). The FFA's formal
  peak (32.7 ms, S/N 5.2) is below our own synthetic noise threshold and is reported as not
  significant. The algorithm's validation is carried by the synthetics (injected 233.7-sample
  period found exactly, S/N 60.5 vs fold-oracle 93.4) — stated plainly in the paper.

## Recover-a-knowns (all also pass on the GPU — the canonical run's device is `cuda`)

- **Chirp round trip**: synthetic impulse dispersed with the exact inverse filter at DM 100
  re-collapses to peak offset 0 with 99.5% of energy re-concentrated (from 2.1% dispersed);
  wrong-DM control stays smeared.
- **SK**: max |torch − jansky.rfi| = 1.5e-14 (exact modulo float).
- **SumThreshold**: sequential mode byte-identical to the `jansky.rfi` oracle (incl. threaded
  masks); the parallel tensor path (pass-start mask, cumsum window means) has Jaccard 0.949 vs
  sequential on synthetic RFI — the evaluation-order difference is documented, not hidden.
  Injected CW line + broadband burst both fully caught. NOTE: float64 medians must come from
  numpy (torch's even-length median takes the lower element — bit us once).
- **FFA**: injected 233.7-sample period found at 233.70 (err 0.0); flat-noise control quiet.

## Real legs (run on the ROCm GPU)

- **CHIME baseband re-dedispersion**: FRB 20181231C, coherent dedispersion of both pols across
  97 good channels; boxcar S/N vs trial DM peaks exactly at the Cat-2 catalogue DM 556.11
  (S/N 4.0 there vs ≤1.4 at ±5/±20 pc cm⁻³ and 2.4 at DM 0). A modest-S/N event in only 38 MHz
  of saved band — the *structure* (peak at catalogue DM) is the validation.
- **Crab filterbank RFI**: torch parallel mask vs CPU oracle on real Parkes data:
  Jaccard 0.916.
- **Crab FFA**: honest null (above).

## GATE-2 (PASS with required fixes, all applied — every one a wording/disclosure item)

- **Chirp sign convention**: the round-trip test is sign-self-consistent by construction and
  cannot pin the sign; which kernel sign dedisperses depends on the backend's
  sideband/conjugation convention (dspsr parameterises this). Ours is anchored EMPIRICALLY by
  the CHIME leg; "the dspsr convention" label dropped; other backends may need the conjugate.
- **Phase-magnitude claim fixed**: ~1e4 rad per 0.39-MHz channel at CHIME DMs (not 1e6 — that
  is the wideband regime, where the float64 requirement genuinely binds).
- **"Peaks exactly at the catalogue DM" → grid-scoped**: six trial DMs (0, ±5, ±20 around
  catalogue); a machinery validation, not a DM measurement. The DM-0 trial's 2.4 reflects
  undispersed RFI/noise structure — noted.
- **Circular-boundary semantics stated**: chirp multiply and inter-channel gather are circular
  (no overlap-save discard region) — fine for mid-buffer bursts, needs overlap-save for
  streaming; documented in docstrings and the paper.
- Also: benchmark single-run caveat in the paper; the "9× FFT applies at scale" extrapolation
  softened to an expectation; FFA S/N gap explanation extended (integer-shift quantisation +
  power-of-two row truncation); figure now shows the run's own periodogram (provenance);
  macros never mislabel GPU timings as CPU.

## Benchmarks (same code, same venv, torch 2.12.1+rocm7.1, RX 7600 XT vs Ryzen CPU)

| kernel | CPU | GPU | verdict |
|---|---|---|---|
| chirp (64 ch × 1M samples) | 1.57 s | 1.97 s | ~parity: transfer/plan-bound at this size; f64-phase→wrapped-f32 trick made CPU 2.4× faster too |
| SumThreshold2d (8192×256) | 2.96 s | 7.85 s | **GPU slower** — per-series host loop (the torch-fdmt wall); batched variant is stated future work |
| FFA (2²² samples, 64 periods) | 6.82 s | 0.65 s | **10.5× GPU** — the gather merges vectorise fully |

The suite's honest headline mirrors torch-fdmt's: portability is delivered (identical results
on ROCm with zero code changes); speed is kernel-shape-dependent, and we say which shapes win.

## Reproduce

CPU: `uv run --extra fdmt --extra voyager python -m jansky_research.torchdsp --benchmark --out .`
GPU: `PYTHONPATH=src:../jansky/src ~/.venvs/rocm-test/bin/python -m jansky_research.torchdsp
--device cuda --benchmark --out .` (pinned torch 2.12.1+rocm7.1; h5py+matplotlib in the venv).
Offline CI leg: `--offline`.
