# 51 — First external FAST-FREX benchmark + an open-weights burst classifier (ROCm)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — confirm the Science
Data Bank download works from a US connection and pin the total GB (currently unknown)

## Context

Verified in the scan: every published FAST-FREX user shares authors with the dataset/DRAFTS
team, and DRAFTS itself is CUDA-locked — so a zero-author-overlap external benchmark with a
pure-PyTorch classifier and released weights is still open, and doubles as the repo's ROCm
demonstrator (the workstation rule: pure PyTorch runs, CUDA-kernel tools do not). Data:
FAST-FREX, Science Data Bank, doi:10.57760/sciencedb.15070 — 600 real FAST bursts + 1,000 RFI
FITS files. Method: plain-torchvision ResNet/U-Net classifiers trained on the released splits,
compared against the paper's RaSPDAM baseline, with weights + training recipe published. This
extends the `torch-fdmt` GPU arc (pinned ROCm wheel, gfx1102 patterns documented in
`torchfdmt-findings.md`). Honest framing: a benchmark paper with a modest ceiling — the value is
the external, reproducible, open-weights data point.

## Deliverables

- `src/jansky_research/fastfrex.py`: `fetch_fastfrex` (Science Data Bank download + manifest,
  `# pragma`), `load_burst_fits` (FITS → normalized dynamic-spectrum tensors + labels),
  `make_splits` (the paper's split reproduced, or a documented deterministic split),
  `train_classifier` (ResNet + U-Net variants, plain torchvision, ROCm), `evaluate_benchmark`
  (recall/precision/F1 matching the paper's metric definitions), `synthetic_burst_set`
  (injected dispersed pulses + RFI textures → training loop smoke-recovery offline),
  `export_weights` (checkpoint + recipe), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/fastfrex/`; `survey/fastfrex-findings.md`; wiring.

## Approach

0. **GATE 0:** confirm doi:10.57760/sciencedb.15070 downloads from a US connection and record
   the total GB (must fit the ~275 GB disk budget); full-text pass on the FAST-FREX paper and
   its citing list to **confirm no zero-author-overlap external benchmark has landed in the
   interim**; pin the paper's exact metric definitions and splits.
1. Tooling + synthetic recover-a-known: a small synthetic burst/RFI set trained end-to-end
   offline; the loop must converge and the evaluation code must reproduce hand-computed metrics.
2. Real leg: download, train ResNet and U-Net variants on the ROCm GPU (days-long runs are
   in-budget), evaluate against the RaSPDAM baseline numbers, ablate input normalization.
3. GATE-2 science review: honest benchmark framing (external reproduction, not SOTA-chasing),
   split/metric fidelity to the original paper, dataset-shift caveats (FAST-specific RFI),
   weight-release checklist.
4. Paper: the first external FAST-FREX benchmark + released weights and recipe.

## Verification

The FAST-FREX paper's own published metrics are the oracle — the evaluation harness must
reproduce the RaSPDAM baseline numbers from its stated definitions before our models are scored;
checks green; GATE-2 sign-off.

## Risks & mitigations

- **Modest ceiling (benchmark paper)** → scoped small; the open weights + ROCm recipe are the
  durable contribution; do not inflate claims beyond the benchmark.
- **An external benchmark lands first** → the GATE-0 citing-list check is mandatory; if one
  exists, drop or rescope to the open-weights/ROCm-reproducibility angle only.
- **Science Data Bank access from the US** (speed/size unknown) → GATE-0 measures both before
  any code; mirror locally once fetched.
