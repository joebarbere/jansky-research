# Findings — torch-dsp: the pure-PyTorch coherent-DSP suite (plan 43)

`jansky_research.torchdsp` extends the merged `torch-fdmt` arc with three kernels no
pure-PyTorch (or JAX) implementation of which existed: coherent dedispersion, SumThreshold +
spectral-kurtosis RFI excision, and a radix-2 FFA. One `device=` argument covers CPU, CUDA, and
ROCm — the BENCHMARKS below were measured on the RX 7600 XT (gfx1102). The science legs ran
on CPU (`results/torchdsp_metrics.json`: `device: cpu`, `benchmark_device: cuda`).

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
  peak (32.7 ms, S/N 5.2) sits at the measured same-search-volume noise level (5.0, committed
  2026-08-25) and is reported as not significant. The algorithm's validation is carried by the synthetics (injected 233.7-sample
  period found to the nearest drift-grid point, S/N 60.5 vs fold-oracle 93.4) — stated plainly in the paper.

## Recover-a-knowns (run on CPU; `benchmark_device` is `cuda`, the science `device` is `cpu`)

- **Chirp round trip**: synthetic impulse dispersed with the exact inverse filter at DM 100
  re-collapses to peak offset 0 with 99.5% of energy re-concentrated (from 2.1% dispersed);
  wrong-DM control stays smeared.
- **SK**: max |torch − jansky.rfi| = 1.5e-14 (exact modulo float).
- **SumThreshold**: sequential mode byte-identical to the `jansky.rfi` oracle (incl. threaded
  masks); the parallel tensor path (pass-start mask, cumsum window means) has Jaccard 0.949 vs
  sequential on synthetic RFI — the evaluation-order difference is documented, not hidden.
  Injected CW line + broadband burst both fully caught. NOTE: float64 medians must come from
  numpy (torch's even-length median takes the lower element — bit us once).
- **FFA**: injected 233.7-sample period found at 233.70 — the nearest drift-grid point (grid
  0.0039 samples; the injected value is not on the grid, so "err 0.0" is rounding, not
  exactness); flat-noise control quiet.

## Real legs (run on CPU)

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

## Benchmarks (single session 2026-08-25, one invocation, both columns; hardware introspected)

Superseding note: the table this section carried before 2026-08-25 (chirp 1.57/1.97, ST
2.96/7.85, FFA 6.82/0.65) was one of THREE mutually inconsistent number sets the round-7
referee found (see below); none matched the committed JSON. The table below is from the
single `--benchmark-only` session on the ROCm venv, whose numbers ARE the committed JSON, with
hardware strings from torch/platform introspection (`AMD Radeon RX 7600 XT, torch
2.12.1+rocm7.1` vs `x86_64` CPU, same interpreter, same venv).

| kernel | CPU | GPU | verdict |
|---|---|---|---|
| chirp (64 ch × 1M samples) | 1.73 s | 2.12 s | ~parity (GPU 0.8×): transfer/plan-bound at this size |
| SumThreshold2d (8192×256) | 3.01 s | 7.49 s | **GPU slower** — per-series host loop (the torch-fdmt wall); batched variant is stated future work |
| FFA (2²² samples, 64 periods) | 6.88 s | 0.65 s | **10.6× GPU** (derived macro `\tdRealSpeedupFfa`) — the gather merges vectorise fully |

Portability is now measured, not asserted: the same session re-ran the full offline kernel
suite on the ROCm device and committed it as `cross_device_check` — SK max diff 1.5e-14,
ST sequential byte-equal to the oracle, parallel Jaccard 0.949, FFA 233.7, dedisp round trips
0.995 — identical to the CPU values.

## Reproduce

CPU science legs: `uv run --extra fdmt --extra voyager python -m jansky_research.torchdsp --out .`
Benchmark (single session, both columns + cross-device check, merged via
`preserve_live_results`): `PYTHONPATH=src:../jansky/src ~/.venvs/rocm-test/bin/python -m
jansky_research.torchdsp --benchmark-only --device cuda --out .`
(pinned torch 2.12.1+rocm7.1; h5py+matplotlib in the venv). Offline CI leg: `--offline`.
Note the old GPU reproduce command (`--device cuda --benchmark --out .`, no `--benchmark-only`)
would re-run the science legs on the GPU and relabel them — the trap that invited the splice.

## Full referee round (2026-08-25): MAJOR REVISION, 17 findings, two BLOCKERs

The software contribution is real, the novelty scoping unusually careful, the Crab null
honestly framed, the sequential-SumThreshold byte-identity a genuine oracle. But the second
headline — the benchmark — rests on a committed JSON block no single invocation of the code
could have produced, and the abstract's speedup matches neither that block nor any other
number in the repo.

**BLOCKER 1: the committed benchmark table is a hand-assembled splice, and three mutually
inconsistent CPU sets exist.** `git show 602e0ca` proves it: the `benchmark_cpu` block is
byte-identical to the 2026-08-04 CPU-only run merely re-keyed, while the commit message of
that very commit quotes the *fresh* CPU numbers the GPU session measured (FFA 7.27 s, ST
3.04 s) — which appear nowhere in the repo; `benchmark_device`/`benchmark_hardware` are
written by no code (grep: only `singlepulse.py` introspects, `torchdsp.py` only reads); and
this file's own table carries a third set disagreeing with the committed JSON in 5 of 6
cells. The paper's "same code, the same virtual environment, and the same consumer hardware"
is false as the numbers stand (CPU leg: repo venv, six days earlier; GPU leg: ROCm venv,
Python 3.14; CPU hardware recorded nowhere).

**BLOCKER 2: the abstract's "10×" is hand-typed and contradicts its own parenthetical** —
"gains 10× (8.42 s to 0.65 s)" is 12.95; the 10× traces to this file's uncommitted 6.82.
No derived speedup macros exist (the torchfdmt fix never ported).

**MAJORs:** the chirp sign is "anchored empirically by the CHIME leg" but the conjugate-sign
arm was never run (one `burst_snr` call with the kernel conjugated settles it); the chirp's
99.5% round trip cannot fail (built with the exact inverse chirp — it verifies float64
self-composition, not the dispersion law) and validates `coherent_dedisperse` (complex128)
while every science/benchmark result uses `dedisperse_channelized` (wrapped, complex64) — add
a group-delay-vs-cold-plasma check and run the round trip through the shipped path; the Crab
null's "data length, not the algorithm, is the limit" has no measurement behind it (no
injected-amplitude ladder, no same-volume noise threshold — the referee estimates max-S/N ≈ 5
is the *expected* noise level for that search volume, which supports the null but must be
committed); the ST Jaccard (0.949/0.9157) is one number at the least-divergent configuration
(n_iter=1) with the over/under-flagging split uncomputed; the central portability claim
("device-agnostic... no science result depends on which device") is conceded unmeasured at
line 155 — one `run(offline=True, device="cuda")` on the ROCm venv commits the cross-device
comparison and this file's "identical results on ROCm" must be fixed; the CHIME S/N is a
whole-series boxcar maximum with no argmax recorded (cannot distinguish burst sharpening
from a wandering noise maximum — one integer per trial DM).

**MINOR/NIT:** "recovers an injected period exactly (0.0 samples)" is a rounding artifact
(the injected 233.7 is not on the drift grid; true error ~0.002 samples, grid 0.004); "gaffa
... that we cite" cites nothing; the CHIME leg is *undersold* — at DM 556 the intra-channel
smear (~3–28 ms) is 1–2 orders above the 0.256 ms boxcar, so S/N 4.0 at all is a real test
of the chirp; the tdReal* namespace labels the run, not the quantity (five synthetic-kernel
oracles sit under tdReal, and the header's "offline rebuild resets to placeholders" claim
contradicts `preserve_live_macros`); the documented GPU reproduce command would rewrite the
CPU science legs — no `--benchmark-only` mode exists, which is precisely what invited the
splice; benchmark shapes and real-leg descriptors (97 ch, 38 MHz, 2.1 s) hand-typed, not in
evidence; `cat2` lacks volume/pages (ApJS 283, 34) and `kania2026` lacks its published DOI
(10.3847/1538-3881/ae0d86); the "~10⁴ rad" phase figure is band-bottom (2.6×10³ at 600 MHz);
`\tdRealSkMaxDiff` exists unformatted and unused.

**Checked and clean:** no `--` macros; every non-macro prose number matches the JSON except
the 10×; dm_fitb 556.1104 ✓; file sizes ✓; `\software{}` cites both toolkits.

**Status: fixes pending.** The single change: make the benchmark producible by one command
(`--benchmark-only` through `preserve_live_results`, torch-introspected hardware for both
columns, derived speedup macros), run it once in a single ROCm session — that alone closes
both blockers, names the CPU, and reconciles this file and the README. The ROCm offline run
and the conjugate-sign trial in the same session turn both central claims into measurements.

**Status: RESOLVED (2026-08-25).** Both blockers closed by making the benchmark producible and
producing it; every "control arm never recorded" MAJOR got its arm in the same two runs.

**Blocker 1 (spliced benchmark):** a `--benchmark-only` mode now computes both timing columns
in one interpreter session, introspects the hardware strings for BOTH devices (the fields no
code wrote are now written only by code), and merges through `preserve_live_results` so the
`_merge` record documents any cross-run composition. Run once on the ROCm venv: chirp
1.73/2.12 s, ST 3.01/7.49 s, FFA 6.88/0.65 s. The three inconsistent CPU sets are reconciled:
this file and the README now carry the committed numbers, with the old table marked superseded.

**Blocker 2 (hand-typed 10x):** `\tdRealSpeedup{Ffa,St,Chirp}` are derived in `_write_macros`
from the two committed timings and used in the abstract and §5 (FFA 10.6x). No ratio is typed.

**MAJORs:**
- Conjugate-sign trial run: S/N 1.4 with the conjugate kernel vs 4.0 with the shipped sign at
  the catalogue DM — "anchored empirically" is now a measurement (`chime_snr_conjugate_sign`).
- The chirp is checked against the cold-plasma law itself (`chirp_group_delay_residual`: max
  1.2e-8 samples across the channel — fails loudly for a wrong constant/exponent), and the
  round trip now also runs through `dedisperse_channelized`, the shipped wrapped-complex64
  path (0.995 energy re-concentration on both paths). The abstract scopes the 99.5% as a
  numerical-precision bound.
- The Crab null has a selection function: the same search volume on pure noise gives S/N 5.0
  (so the 5.2 "peak" is the measured noise level, as the referee estimated), and an injected
  pulse train at the published period into the real dedispersed series is recovered only from
  2.0 sigma per pulse (`crab_injection_ladder_sigma`). "The data length is the limit" is now
  a statement with evidence. Duration/tsamp/counts emitted (2.097 s, 512 us, 4096 samples).
- ST agreement swept (thr 3.0/3.5/4.0 x iter 1/2): Jaccard 0.88-0.97, degrading with
  occupancy and iterations as the mechanism predicts; the symmetric difference split by
  direction (synthetic: 54/51 balanced; Parkes: 4537 par-only vs 2569 cpu-only — the parallel
  mode over-flags, data loss not missed RFI).
- Portability measured: the same ROCm session re-ran the offline kernel suite on the GPU and
  committed `cross_device_check` (identical to CPU). §3 now cites it; the old "not measured
  here" concession and this file's unsupported "identical results" claim are both replaced by
  the measurement.
- CHIME argmax committed per trial DM (`chime_peak_sample_vs_dm`): off-catalogue maxima fall
  at unrelated positions (noise excursions in trials whose ~7 ms residual smear buries the
  burst); the paper says so and adds the referee's undersold point — at DM 556 the
  intra-channel smear (~3-28 ms) is 1-2 orders above the 0.256 ms boxcar, so any detection at
  all requires the coherent chirp and its sign to be right.

**MINOR/NIT:** "recovers exactly (0.0 samples)" reworded to nearest-drift-grid-point with the
grid spacing emitted (`\tdRealFfaDriftRes` = 0.0039); the S/N-gap explanation is now a check
(predicted 0.67 vs measured 0.65, both committed); `gaffa` is a real citation
(github.com/lintian233/gaffa, created 2026-06-12 — found and verified via the GitHub API);
`cat2` completed (ApJS 283, 34) and `kania2026` given its published DOI (10.3847/1538-3881/
ae0d86); the macro header now says the tdReal/tdSyn prefix labels the RUN MODE, not the
quantity's provenance, and no longer claims a rebuild resets tdReal*; SK cited via the
formatted `\tdRealSkMaxDiffSci`; the phase figure is "10^3-10^4 rad (band bottom to top)";
benchmark shapes and real-leg descriptors (97 ch / 37.9 MHz / 2.097 s) all emitted from code.
