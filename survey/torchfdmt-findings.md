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

**[Superseded twice — the numbers of record are in `results/singlepulse_metrics.json`, all four
columns from a single invocation (2026-08-25). Earlier versions of this section carried a
36.1 s / 24x table from a run that was never committed, and then a spliced row combining an
Aug-4 CPU invocation with an Aug-10 GPU one; see the referee round below.]**

**The algorithm beats the hardware:** O(N log N) FDMT on a CPU outruns the GPU-accelerated brute
force. Our FDMT's own GPU gain is ~1× — the per-delay host loop dominates; a batched-gather
vectorisation of the merge is stated future work. We do NOT compare against tuned CUDA
dedispersers on datacentre GPUs (different hardware class). Reproduce (one invocation, all
columns, hardware string emitted by torch introspection):
`PYTHONPATH=src:../jansky/src ~/.venvs/rocm-test/bin/python -m jansky_research.singlepulse
--out . --device cpu --benchmark --bench-devices cpu,cuda`

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
invocation can produce the committed combination. [Correction, 2026-08-25: the "same-run" claim
that followed here was itself wrong — 0.45 s was the Aug-4 CPU invocation and 1.50 s the Aug-10
GPU one. The genuinely same-run pair (37.9 s CPU brute vs 1.5 s GPU, from the Aug-10 run's
commit message) gives 25.3x, not the 29x the spliced row implied; the committed row now comes
from one invocation with `--bench-devices` decoupled from the science device.]

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

## Full referee round (2026-08-24): MAJOR REVISION, 15 findings, two BLOCKERs

**BLOCKER 1: the committed figure shows no butterfly peak at the DM the paper says it peaks at.**
The left panel plots `plane.max(dim=1)` -- the RAW track sum, which rises with delay row because
higher rows integrate more samples -- so the displayed maximum is at DM ~118, not 56.59, and the
claimed detection is invisible. `best_dm` comes from `FDMTResult.best()`, which median/MAD-
normalises per row first. The figure plots a different statistic from the one that produced the
quoted number. Fix: plot the normalised curve (or both, labelled).

**BLOCKER 2: "recovers ... at butterfly S/N = 6.0" carries no trials factor.** The plane is
~1914 x ~4102 ~ 7.9e6 cells; the expected noise maximum is ~5.2 (Gumbel p(6.0) ~ 0.01, ~2.3
sigma), and the suite's own noise-only test asserts snr < 10 -- a bound that does not exclude
the paper's detection value. The evidence that IS strong and unquoted: the peak lands 2.9 DM
trials from the catalogue value out of ~1914 (p ~ 0.003 positionally) and the boxcar S/N 14 on
~2.5e4 trials is overwhelming. Demote the butterfly height; state the coincidence + boxcar. A
~200-rep per-channel circular-shift null (into a tmpdir) should be quoted beside the 6.0.

**MAJORs:** the 29x brute speed-up divides the Aug-4 CPU time (44.12 s) by the Aug-10 GPU time
(1.5 s) while the Aug-10 invocation's own CPU leg (37.9 s, committed twice in CHANGELOG) was
never written -- same-run ratio 25.3x, and the splice retained the CPU number that flatters the
GPU ratio by 16% (the FDMT "~1x" pair is likewise cross-invocation); `benchmark_device` /
`benchmark_hardware` are written by NO code in the repo (hand-entered in 602e0ca) and
`preserve_live_results` will now perpetuate them onto future runs on different hardware -- emit
them from torch introspection and re-run; the benchmark's device set is derived from the science
leg's device (`devs = ("cpu","cuda") if device != "cpu"`), so the combination the paper reports
(CPU science + both benchmark columns) is unreachable by any single invocation -- add a
decoupled --bench-devices; "Every production GPU dedisperser ships CUDA kernels" is falsified by
Sclocco et al. 2016 (OpenCL, auto-tuned, deployed in AMBER/ARTS) -- keep the narrower "no
maintained PyTorch/JAX FDMT" claim.

**MINOR:** the findings file still displays the retracted benchmark table (36.1 s / 24x / 3.6x)
above its own retraction, and calls 0.45-vs-1.50 "same-run" when they are not; three fixes
described in the past tense are absent from the committed evidence (`source` still names only
the real leg above synthetic values; `real_dm_step_pc` / `bench_devices` still missing, and the
stated blocker to re-running -- results-clobber -- is gone now that preserve_live_results
merges); \spTrueDm (the synthetic injection default) is cited as the ATNF catalogue value -- add
\spCatDm from CRAB_DM; sp_pos/sp_width/n_time uncommitted while ~47% of the best-DM series is
boundary-affected and the figure's tallest sample is not the marked one; the brute baseline is
gather-bound (1.5 GB index tensors -- say the baseline is the naive in-framework gather;
`numpy_oracle_reduced_s` is measured and discarded); the Crab giant-pulse rate claim is
uncited.

**NIT:** \spSnr (synthetic S/N 113.5) generated and unused -- the 113.5-vs-6.0 contrast belongs
in the paper; "Every number is \input" is not literally true (instrument constants are typed,
all verified correct); the untracked arxiv-submission still carries the retracted 24x.

**Status: fixes pending** (one fdmt.benchmark re-run on the ROCm venv; one replot; one null run
into a tmpdir).

**Status: RESOLVED (2026-08-25).** All 15 findings addressed; the referee's estimates were
confirmed by measurement on both blockers.

**Blocker 1 (figure):** `search` now computes the per-row z-scored peak -- the same statistic
`best()` maximises -- and the figure plots it (the raw track sum is committed as a diagnostic
key only). A test pins curve-argmax == best_dm and curve-max == best_snr.

**Blocker 2 (trials):** a 200-rep per-channel circular-shift null is code
(`shift_null`) and committed: median 5.23, p99 6.03, max 6.26 -- the referee's Gumbel estimate
(~5.2) to the decimal. The observed butterfly 6.0 has p = 0.030 against it; the paper now leads
with the positional coincidence (2.9 of 1,903 trials, p = 0.0035) and the boxcar S/N 14, and
quotes the peak height only against the committed null.

**The benchmark, re-measured in one invocation** (`--bench-devices` decoupled from the science
device; hardware string emitted by torch introspection, never typed): brute 36.5 s CPU ->
1.49 s GPU (**24x**, not the spliced 29x -- the splice had retained the older, slower CPU
number exactly as the referee suspected), FDMT 0.41 s CPU / 0.44 s GPU (~1x), NumPy oracle
19.9 s committed. The FDMT-CPU-beats-brute-GPU headline is 3.6x. The ~15% CPU run-to-run
spread is stated in the paper as the honest wall-clock systematic.

**A guard bug found en route:** `preserve_live_macros`' synthetic-source check used the naive
"synthetic in source" rule, so torchfdmt's MIXED source string ("synthetic injection ... +
Parkes ... real") made every real rerun look like a synthetic downgrade and its macro updates
were silently discarded -- the fresh benchmark left \spBruteSpeedup at 29 until caught. The
rule now matches `preserve_live_results` (mixed counts as real), with tests for both
directions. This is the same marker-mismatch lesson as the results-guard fix of 2026-08,
arriving in the macro guard.

**The rest:** sp_pos/sp_width/n_time committed with the boundary condition stated (~46% of the
best-DM series integrates truncated tracks); \spCatDm sourced from CRAB_DM and cited at the
catalogue sentence; \spSnr (synthetic 113.5) now used as the contrast to the real 6.0; the
universal "every production GPU dedisperser ships CUDA" claim narrowed and Sclocco et al. 2016
(OpenCL, AMBER/ARTS) cited as the existing portable route; the gather-bound baseline stated;
the uncited giant-pulse rate dropped; the stale findings table and the wrong "same-run 3.3x"
claim corrected in place; provenance comment scoped to result numbers; arXiv package rebuilt
clean under the abstract limit.
