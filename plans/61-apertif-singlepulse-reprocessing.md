# 61 — Apertif Time-Domain DR2 single-pulse reprocessing (torch-fdmt's real-data leg 2)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — probe the ASTRON
tape-staging request turnaround (free helpdesk) before committing to a schedule

## Context

Apertif Time-Domain DR2 (1,666 pointings, ~0.48 GB/pointing, 1-bit Stokes-I PSRFITS) was
FRB-searched with AMBER/ALERT (arXiv:2406.00482) but has no published archive-wide single-pulse
RRAT / known-pulsar census (fable-ideas F24). This slice is the second real-data leg for the
merged `torchfdmt` slice: `fdmt.py` (GPU dedispersion, ROCm-verified at 24× vs CPU) +
`singlepulse.py` over staged pointings, with ATNF position-matched folding as the ground truth.
Honest framing per the idea's own risk field: the ceiling is modest — this is tooling validation
at scale (does the pure-PyTorch stack survive 1,666 heterogeneous archival pointings?), with
redetections and any candidate pulse trains as the science by-product. Disk-aware batching is
mandatory: ~275 GB free means staging, processing, and deleting in waves.

## Deliverables

- `src/jansky_research/apertifsp.py`: `stage_pointings` (ASTRON VO/tape-staging requests,
  `# pragma`), `read_psrfits_1bit` (1-bit Stokes-I PSRFITS reader), `dedisperse_pointing`
  (`fdmt.py` reuse, GPU), `search_single_pulses` (`singlepulse.py` reuse: matched-width boxcar +
  clustering), `atnf_field_check` (known pulsars/RRATs per pointing, position-matched folding),
  `candidate_trains` (repeated-DM/position event grouping), `batch_ledger` (resumable
  stage→process→delete accounting), `synthetic_pointing` (injected pulses in 1-bit noise →
  recovery), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/apertifsp/`; `survey/apertifsp-findings.md`; wiring.

## Approach

0. GATE 0: file a small ASTRON staging request and measure the turnaround (this sets the whole
   schedule); full-text pass on arXiv:2406.00482 and the DR2 release paper to confirm no
   archive-wide single-pulse census exists; verify the 1-bit PSRFITS format docs.
1. Tooling + synthetic recover-a-known: injected single pulses (known DM, width, S/N) into
   synthetic 1-bit data; end-to-end FDMT + boxcar recovery at stated completeness.
2. Real leg: staged waves of pointings, prioritizing fields containing catalogued pulsars/RRATs;
   per-field known-source check, then the blind single-pulse pass; multi-week GPU job, resumable.
3. GATE-2 science review: 1-bit quantization sensitivity caveat, per-pointing RFI environment
   heterogeneity, candidate-train vetting discipline (no single-event discovery claims).
4. Paper: the census (redetections + candidates + per-field completeness) framed as pure-PyTorch
   tooling validated at archive scale.

## Verification

Every staged field containing a catalogued pulsar/RRAT must yield ≥1 known-source redetection
before its blind candidates count; synthetic injection recovery at stated completeness; checks
green; GATE-2 sign-off.

## Risks & mitigations

- **Tape staging stalls** → GATE-0 measures turnaround first; process in resumable waves so a
  stalled request idles the queue, not the slice; the synthetic + early-wave results stand alone.
- **Modest science ceiling** → framed from day one as tooling validation at scale (the
  `torchfdmt` arc's real-data leg 2); redetection statistics are the honest headline.
- **Disk pressure** (1,666 × 0.48 GB ≫ 275 GB) → strict stage→process→delete ledger; never hold
  more than one wave on disk.
