# 34 — `torch-fdmt`: a pure-PyTorch Fast DM Transform + the single-pulse slice (absorbs plan 28)

Status: 🔨 in progress — `fdmt.py` tooling MERGED (torch FDMT + brute twin validated vs the
`jansky.transients` oracle; 100% module coverage; the delay-split swap bug the plan predicted was
caught by the oracle checks). GATE 0 PASSED: a 3.25 MB real Parkes/UWL **Crab** filterbank
(sigpyproc3 test tree, no auth; ATNF DM 56.77 / P 33.27 ms; BL products ruled out — 1-18 s sampling
makes the sweep invisible). `singlepulse.py` MERGED: minimal pure-NumPy SIGPROC reader + FDMT
search + folding; synthetic recovers DM 56.63 / P 33.257 ms; **REAL Crab leg recovers DM 56.59
(0.3% error), single-pulse S/N 14 (a giant pulse)**. Benchmark (8192×1024×1009 DM): brute CPU
36.1 s → GPU 1.50 s (24×); **FDMT CPU 0.41 s beats brute-GPU 3.6×** (the algorithm beats the
hardware); FDMT GPU ≈ CPU (Python delay-loop bound — kernel vectorisation is future work).
Next: paper + Makefile/Snakefile wiring + findings + GATE-2

## Context

Two facts met in the 2026-07 scan. First, the repo still lacks the one classic radio method plan 28
scoped: **incoherent dedispersion + single-pulse/period search of a broadband dispersed signal**
(`jansky.transients` implements all of it, untouched by the research repo). Second, the GPU addendum
**empirically verified PyTorch-ROCm on this machine's RX 7600 XT** (gfx1102, 16 GiB): the
torch 2.12.1+rocm7.1 wheel runs natively (no HSA override, no system ROCm), and a *naive*
dedispersion-style gather+sum measured **24× over CPU** (1024 DM × 4096 chan × 8192 samp:
156 s → 6.5 s, results matching).

The tooling survey found the gap that turns those into one slice: **no maintained PyTorch (or JAX)
FDMT exists** — the Fast DM Transform (Zackay & Ofek 2017, ApJ 835, 11; O(N log N) vs the naive
O(N·N_DM)) has only a Julia implementation (`FastDMTransform.jl`) and a Python reference
(`pyfdmt`), while every production GPU dedisperser (Heimdall/`dedisp`, astro-accelerate, FREDDA)
is **CUDA-walled**. A pure-tensor-ops FDMT runs on ROCm (and CUDA, and CPU) for free — a genuine,
verified-empty niche, JOSS/RNAAS-shaped, and it upgrades plan 28's dedispersion leg on this exact
hardware. The rule of thumb from the survey: *pure PyTorch runs; shipped CUDA kernels don't.*

**Plan 28 is absorbed, not discarded:** its GATE 0 (find one small public filterbank of a bright
pulsar), its synthetic-fixture discipline, and its recover-the-catalogued-DM/period science are
this plan's steps 0, 2, and 4 verbatim. Update plan 28's status line to point here.

## Deliverables

- `src/jansky_research/fdmt.py` — the algorithm module:
  - `fdmt(spectrum, f_lo, f_hi, max_dm, ...)` — the Zackay & Ofek recursion in vectorised tensor
    ops (`cumsum` + gather/scatter; no custom kernels, no CUDA extensions), device-agnostic
    (`cpu` / `cuda`(ROCm) via one `device=` argument).
  - `brute_dedisperse()` — the naive roll-and-sum in the same tensor style (the benchmark
    baseline and the small-input cross-check).
  - `jansky.transients` (`dispersion_delay`/`disperse_pulse`/`dedisperse`/`dm_search`) is the
    **CPU correctness oracle** — the tested reference the torch outputs must match.
  - `benchmark()` — CPU-numpy oracle vs torch-CPU vs torch-GPU wall-time table over problem sizes
    (the paper's headline table; honest target: "faster than the CPU baselines", **not** "beats
    FREDDA" — different hardware class).
- `src/jansky_research/singlepulse.py` — plan 28's science module, now with an engine switch:
  compose `jansky.formats.read_filterbank` + `jansky.rfi` flagging + FDMT (or brute) DM–time
  plane + `transients.boxcar_snr` single-pulse detection + `epoch_folding_search`;
  `synthetic_observation` (injected dispersed pulse train at known DM/period/width);
  `fetch_bl_filterbank` (`# pragma: no cover`); `run(offline=...)`, `_figure` (DM–time
  butterfly + folded profile + the benchmark bars), `_write_macros`, `_main`.
- **Torch as an optional extra** `[project.optional-dependencies] fdmt = ["torch"]`, with a
  `[tool.uv.sources]`/CI pin to the **CPU wheel index** (`download.pytorch.org/whl/cpu`) so CI
  never pulls the multi-GB accelerator wheels; tests `importorskip("torch")` — the 85% floor is
  carried by the torch-CPU path in CI; the GPU device path differs only by the `device=` string
  (assert-same-results test runs on GPU only when available, `# pragma: no cover` guarded).
- `tests/test_fdmt.py` + `tests/test_singlepulse.py`; `data.py` entry for the chosen public
  filterbank; `papers/torchfdmt/` (AASTeX, astro-ph.IM-shaped; standalone-PyPI/JOSS extraction
  deferred to a follow-up); `survey/torchfdmt-findings.md`; Makefile/Snakefile wiring.

## Approach

0. **GATE 0 — data reality check (plan 28's, unchanged; do FIRST).** One public filterbank that is
   (a) a known bright pulsar (B0329+54 / Vela / Crab) or a catalogued FRB, (b) downloads in
   seconds–minutes, (c) reads with `jansky.formats.read_filterbank`/`read_guppi`. Search order:
   BL open data, **Parkes PSRDA** (CSIRO DAP, arXiv:2511.22702), **Apertif TD DR2** (ASTRON),
   FAST-FREX samples (Science Data Bank). Hard fallbacks that de-risk the slice entirely: the
   jansky-vendored FRB 121102 `.fil` (~1.6 MB) and `your_28.fil`. Record URL + size in `data.py`.
1. **FDMT tooling + oracle validation.** Implement `fdmt()`; on synthetic dispersed pulses
   (forward model = `transients.disperse_pulse`), the FDMT DM–time plane's peak must match the
   `transients.dm_search` oracle DM and S/N within tolerance, and `brute_dedisperse` must match
   `fdmt` element-wise on small inputs. Property tests: zero-DM passthrough, delay quantisation
   error bound (≤1 sample), energy conservation along tracks.
2. **Synthetic recover-a-known (plan 28 step 1).** Injected pulse train at known (DM, P, width):
   `run(offline=True)` recovers all three via FDMT + boxcar + folding within tolerance.
3. **Benchmark leg (the GPU addendum made real).** `benchmark()` grid over (n_chan, n_samp, n_DM):
   numpy oracle vs torch-CPU vs torch-GPU (RX 7600 XT) for brute and FDMT; report the crossover
   sizes and VRAM ceiling honestly (the addendum's 64-DM-batch OOM lesson: document batching).
   Pin the working torch/rocm version in the findings (the official matrix still omits gfx1102).
4. **Real recover-a-known (plan 28 step 2).** On the GATE-0 file: RFI-flag, FDMT over a DM grid,
   recover the **catalogued DM (ATNF/TNS) and period** — cross-checked against the
   `jansky.transients` CPU oracle on the same data. This is the validation, not a discovery.
5. **GATE-2 science review** — catalogue match within errors; caveats surfaced: uncalibrated
   power scale, scattering/scintillation broadening, boxcar noise model, single pointing, and the
   benchmark's scope ("consumer-GPU, vs CPU baselines" — name what was *not* compared).
6. **Write-up** `papers/torchfdmt/` — a tool + benchmark + recover-a-known methods paper; note
   the deferred standalone-package release so the JOSS path stays open.

## Verification

- `make test` / `cov` green with the torch-CPU extra (85% floor); `ruff` + `mypy` clean; CI never
  downloads accelerator wheels or large data files.
- FDMT matches the `jansky.transients` oracle (DM, S/N) and `brute_dedisperse` (element-wise,
  small inputs); offline `run` recovers injected DM/period/width.
- GPU-vs-CPU benchmark table reproduces on this machine (`make reproduce` runs the benchmark leg
  when a GPU is present, records `device: cpu` gracefully when not).
- Real leg recovers the catalogued DM/period; GATE-2 sign-off.

## Risks & mitigations

- **FDMT indexing complexity (highest technical) →** the recursion's partial-sum indexing is
  fiddly; the three-way oracle (transients / brute / fdmt) plus property tests catch it; keep the
  brute path as the always-correct fallback engine.
- **`torch.compile`/HIP kernel quality on gfx1102 →** don't depend on it: eager-mode tensor ops
  already measured 24×; treat `torch.compile` as an optional appendix number.
- **CI weight of torch →** CPU-wheel index pin + `importorskip`; if even that is too heavy for
  CI, the fallback is a numpy FDMT twin carrying coverage with torch tests local-only — decide at
  implementation, record in the findings.
- **File size (plan 28's gate) →** unchanged hard stop + the two vendored fallbacks.
- **Niche gets filled while we build →** re-run the "no torch-FDMT exists" search at GATE-2; if
  one appeared, reframe as an independent implementation + benchmark comparison (still useful).
- **Scope creep toward a full search pipeline →** this is one file, one pulsar/FRB, one benchmark
  table; a blind survey re-search (Apertif/PSRDA at scale) is explicitly a *later* slice.
