# Findings — torch-fdmt: a pure-PyTorch Fast DM Transform + Crab recover-a-known

`jansky_research.fdmt` implements the Zackay & Ofek (2017) Fast DM Transform as **pure tensor
ops** — no CUDA kernels, so one `device=` argument runs it on CPU, CUDA, or ROCm. Every production
GPU dedisperser (Heimdall/`dedisp`, astro-accelerate, FREDDA) is CUDA-only, and no maintained
torch/JAX FDMT existed (verified 2026-07) — this fills that portability gap.
`jansky_research.singlepulse` is the science leg: a minimal pure-NumPy SIGPROC reader + FDMT
search + boxcar single-pulse + epoch folding.

## Validation (the oracle-first discipline paid for itself)

- Zero-DM row ≡ plain channel sum (exact); the torch brute twin ≡ `jansky.transients.dedisperse`
  (exact); FDMT peak at the injected DM within delay quantisation.
- The oracle checks caught **two real indexing bugs** during development: a swapped high/low
  sub-band delay split, and channel padding that stole delay budget (pad rows must carry zero
  ν⁻² span).
- Documented semantic difference: FDMT integrates the full track *including* intra-channel smear,
  so its S/N on a smeared pulse exceeds the one-sample-per-channel oracle's (tested, not hidden).

## Recover-a-known

| leg | result |
|---|---|
| synthetic pulse train (Crab-like) | DM 56.63 (true 56.77); folded P 33.379 ms (true 33.392, ATNF epoch-corrected) |
| **real Parkes/UWL Crab** (GATE-0 3.25 MB public `.fil`) | **DM 56.59 — 0.3% from ATNF**; giant pulse at boxcar S/N 14 |

GATE 0 note: Breakthrough Listen products were *ruled out* for this demo — their 1–18 s sampling
makes the DM sweep invisible; the sigpyproc3 test-tree Parkes file (832 chan, 702–4030 MHz,
0.512 ms) shows the full 463 ms sweep inside one 2.1 s file.

## Benchmark (RX 7600 XT / gfx1102, ROCm; 8192×1024, 1009 DM trials)

| engine | CPU | GPU |
|---|---|---|
| brute roll-and-sum | 36.1 s | 1.50 s (**24×**) |
| **FDMT** | **0.41 s** | 0.41 s (~1×) |

**The algorithm beats the hardware:** O(N log N) FDMT on a CPU outruns the GPU-accelerated brute
force by 3.6×. Our FDMT's own GPU gain is ~1× — the per-delay host loop dominates; a
batched-gather vectorisation of the merge is stated future work. We do NOT compare against tuned
CUDA dedispersers on datacentre GPUs (different hardware class). Reproduce:
`uv run python -m jansky_research.singlepulse --benchmark --out .` (CPU); add `--device cuda`
from a ROCm venv for the GPU column.

## Honest caveats

- The 2.1 s file limits period precision; power scale uncalibrated; the RFI guard is a
  median/MAD clip, not an excision pipeline.
- torch enters as the `fdmt` extra pinned to the **CPU wheel index** so CI stays light; the
  official ROCm support matrix still omits gfx1102 — the working torch/rocm version is pinned in
  the scan's GPU addendum; re-test on upgrades.

## Correctness pass (2026-08-23) — a hand-typed ratio, and a recovery quoted without its grid

Four defects found by the post-conversion referee round, all predating the style campaign.

**The paper quoted "~24x" for a ratio its own adjacent macros make 29x.** Both the abstract and
the benchmark section read "the GPU accelerates the brute force itself by ~24x (\spBruteCpu ->
\spBruteGpu s)", rendering as "~24x (44.12 -> 1.5 s)" — and 44.12/1.5 = **29.4**. The 24 traces
to a superseded 36.1 s CPU timing recorded in this file's own benchmark table (36.1/1.5 = 24.07),
a number that appears nowhere in `results/`. The file header claims "Every number is `\input`
from the pipeline --- none typed by hand"; this one was, and it went stale.

Fixed at the source: `_write_macros` now **derives** `\spBruteSpeedup` and `\spFdmtSpeedup` from
the committed timings, so no ratio in this paper can be typed or drift again. A CPU-only run
emits `--` for them, which `preserve_live_macros` then refuses to write over a real value.
Pinned by three tests, one of which checks the macro on disk against the JSON on disk.

**The real-data DM recovery was quoted as a bare "0.3%".** The FDMT butterfly indexes rows by
whole samples of dispersive delay, so the recovered DM is *quantised*, not fitted. At this
file's band (702–4030 MHz) and sampling (0.512 ms) one row is **0.0627 pc cm⁻³**, so the 0.18
offset from the catalogue is **2.9 trials** — a grid-limited agreement, not a precision. The
abstract's "within delay quantisation" frame belongs to the synthetic leg and was silently
inheriting to the real one. `real_dm_step_pc` is now recorded by the real leg (`_dm_step`), the
paper states the DM is grid-quantised, and it now quotes `\spRealSnr` = 6.0, which was generated
and used nowhere despite being the number that says how well the peak can be localised.

**The benchmark row is a splice of two invocations.** `git log -p` on the metrics shows an early
CPU-only run (`bench_*_gpu_s: null`, `device: "cpu"`) with GPU values patched in later while
`device` stayed `"cpu"`. The code's own `devs = ("cpu", "cuda") if device != "cpu" else ("cpu",)`
means a `--device cuda` run writes **both** columns *and* sets `device: "cuda"` — so no single
invocation can produce the committed combination. The conclusion survives on same-run numbers
(0.45 vs 1.50 = 3.3x), so this is provenance rather than correctness. The paper now says the
columns come from separate invocations and that the ratios should be read to a significant
figure; `bench_devices` is recorded going forward so it cannot recur silently.

**The evidence file labelled its synthetic block with the real file's name.** `source` was set
in the real block to the Parkes file and sat directly above `recovered_dm: 56.63` — a
**synthetic** value, the real one being `real_recovered_dm: 56.59`. An auditor reading the JSON
top-down takes the wrong number. `source` now names both legs and states the key convention
(unprefixed = synthetic, `real_*` = real).

**Not fixed, and why.** `real_dm_step_pc` and `bench_devices` will populate on the next real run;
they are absent from the committed JSON because re-running to add them would drop the GPU
benchmark columns (a CPU run writes them null), which is the same two-mode hazard
`preserve_live_macros` exists to arbitrate for macros and which the results JSON does not yet
have. The paper therefore states the quantisation qualitatively rather than quoting a trial
count it cannot yet cite from evidence.
