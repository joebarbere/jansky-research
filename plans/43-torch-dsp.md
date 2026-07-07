# 43 — torch-dsp: coherent dedispersion + RFI-excision kernels + FFA in pure PyTorch

Status: ✅ done — GATE 0 passed 2026-07-06 (all three niches OPEN for pure PyTorch; wording
fences: a SYCL coherent dedisperser and a 3-week-old CUDA FFA scaffold `gaffa` exist → claim
"pure-PyTorch", not "first portable/GPU"; jess is Kania+2026 not Agarwal). All three kernels
shipped + FFA stretch goal landed: chirp round-trip exact (99.5% energy, offset 0), SK/
SumThreshold byte-identical to jansky.rfi in sequential mode (parallel variant Jaccard 0.949,
documented), FFA exact on injections vs a brute-fold oracle (delta: riptide cited not depended
on). Real legs ON THE ROCm GPU: CHIME baseband FRB 20181231C burst S/N peaks exactly at its
Cat-2 DM 556.11; Crab RFI mask Jaccard 0.916; Crab period re-find is an honest null (2.1-s
file: fold S/N 2.2 at the published period — nothing to find). Benchmarks: FFA 10.5× GPU,
chirp parity (transfer-bound), SumThreshold GPU-hostile (host loop — stated future work). See
survey/torchdsp-findings.md.

## Context

Verified per-algorithm in the 2026-07 scan (snippet-level — GATE-0 re-checks full text): no
PyTorch/JAX coherent dedispersion exists (dspsr/CDMT are CUDA-only); no PyTorch
SumThreshold/spectral-kurtosis RFI kernels exist (the freshest tool, `jess` — AJ Jan 2026 — is
CuPy/CUDA-locked; IQRM/AOFlagger are CPU); no GPU FFA of any kind exists (riptide is CPU; its
2021 GPU issue #2 is unanswered). Prior-art fence: `PyTorchDedispersion` (GitHub, 1 star) covers
incoherent dedispersion only; this repo's own merged `torch-fdmt` covers FDMT. Corrections-section
fence: LPTs are discovered by image differencing, not folding searches — the FFA targets
classical pulsar reprocessing, NOT LPT discovery. This is the flagship GPU slice, extending the
`torch-fdmt` paper arc into a coherent-DSP suite (JOSS + science leg); each kernel unlocks later
archival slices (F10 Parkes single pulses, F24 Apertif). Pure PyTorch/ROCm only (gfx1102, 16 GiB,
rocFFT-backed `torch.fft` — the 9× FFT benchmark de-risks it); multi-day benchmarks welcome.

## Deliverables

- `src/jansky_research/torchdsp.py` (+ package scaffolding): `coherent_dedisperse` (chirp-filter
  multiply via `torch.fft`), `spectral_kurtosis` + `sumthreshold` (block moment/threshold
  tensor reductions), `ffa_search` (stretch goal — budget for the batched-gather vectorization
  wall `torch-fdmt` hit), `fetch_chime_baseband` (# pragma), `benchmark_suite` (GPU vs CPU vs
  oracle tools), synthetic fixtures per kernel, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/torchdsp/`; `survey/torchdsp-findings.md`; wiring.

## Approach

0. GATE 0: full-text/repo-level re-verification that the three niches are still empty (dspsr,
   CDMT, `jess`, IQRM, AOFlagger, riptide + a fresh GitHub/ADS sweep); confirm one small CHIME
   baseband file (DOI 10.11570/23.0029; 140 known-DM FRBs, HDF5, 2.56 µs) is fetchable from
   CANFAR without auth.
1. Tooling + synthetic recover-a-known per kernel: chirp-dedisperse a synthetic dispersed pulse
   to zero residual smearing; SK/SumThreshold masks match `jansky.rfi`'s CPU implementations on
   synthetic RFI; FFA recovers an injected period (riptide-CPU as the oracle).
2. Real legs: re-dedisperse a CHIME baseband FRB to its published DM/structure; match
   `jansky.rfi` masks on the vendored Crab `.fil` real RFI; re-find the Crab period with the
   FFA. Benchmark GPU vs CPU throughout (honest framing: portability first, speed second).
3. GATE-2 science/engineering review: ROCm pinning caveats, FFA vectorization ceiling, benchmark
   fairness (oracle tools on their native platforms). 4. JOSS-style package paper + science leg.

## Verification

CHIME baseband FRB re-dedispersed to its published DM/structure; `jansky.rfi` mask match on
synthetic + real RFI; Crab period re-found by the FFA against the riptide-CPU oracle; checks
green; GATE-2 sign-off.

## Risks & mitigations

- **Engineering, not novelty, is the risk**: ROCm/gfx1102 remains officially unlisted → pin and
  re-test torch versions (documented pattern in `torchfdmt-findings.md`).
- **FFA may hit the same vectorization wall as `torch-fdmt`** → declared a stretch goal; ship
  coherent dedispersion + RFI kernels regardless.
- **Someone ships a GPU port mid-flight** → GATE-0 sweep plus per-kernel scoping; the suite +
  ROCm portability story survives a single-kernel scoop.
