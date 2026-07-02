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
