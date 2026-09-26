# 61 — Apertif Time-Domain DR2 single-pulse reprocessing (torch-fdmt's real-data leg 2)

Status: 📋 **re-scoped after the GATE-0 access probe (2026-09-26)** — the whole-archive census as
written is infeasible (the release is **0.76 PB**, not ~0.8 TB); a pulsar-field slice is feasible
and costs $0 on the local GPU. **Blocked on one owner action:** an ASTRON helpdesk account to
request tape staging (turnaround undocumented, "best effort").

## GATE-0 access probe (2026-09-26) — M = measured, D = documented

- **Volume [M]:** VO table `arts_dr2.frb_obs` (TAP https://vo.astron.nl/tap): 1,181,886 files in
  **2,582 observations** (not 1,666), MJD 58667-59618. 1-bit Stokes-I PSRFITS search mode, one
  file per tied-array beam (40 CB x 12 TAB = 480 files/obs). 3-h file: 318.7 MB (pre-2020-05,
  384 ch x 2.048 ms) or 1.36-1.40 GB (768 ch x 0.8 ms). Sum of `filesize` = **761.7 TB**. The
  "0.48 GB/pointing" below was wrong by ~10^5 per observation (it matches one 2019 *file*).
- **Access [M]:** online files are anonymous HTTPS WebDAV
  (`https://alta.astron.nl/webdav/APERTIF_DR2_TimeDomain/<obsid>/CBxx/ARTSxxx_CBxx_TABxx.fits`).
  Only **11 of 24** FRB-detection observations actually download (3.9 TB; 13 return HTTP 500);
  everything else is on tape via the ASTRON Jira helpdesk (account required).
- **Throughput [M]:** 23.9 MB/s single stream to the workstation; ~28 MB/s with 4 streams, but
  two parallel downloads ended short **without error** — the ledger must check Content-Length.
- **Reader [M]:** astropy's TDIM (1,384,1,500) does not match the 24000-byte rows; unpack bits
  with `np.unpackbits` (5.7 s per 3-h 2019 file).
- **Processing [M]:** FDMT to DM 500 on the RX 7600 XT (ROCm): **~21 s per 3-h 2019 file (522x
  real time)**; CPU ~69 s. Post-2020 files (~4x samples, 2x channels) are estimated at ~3 min
  each [inferred], so they are GPU-bound; 2019 files are download-bound. FDMT only — the
  boxcar single-pulse search is extra and not yet timed.
- **Recover-a-known set [M]:** ATNF psrcat v2.8.1 vs 98,491 CB centres: 166 pulsars (9 RRATs)
  within 0.25 deg, in 1,077 obs / 3,489 CBs; **~486 pulsar-calibrator observations** (B0531+21 x189,
  B1933+16 x182, B0950+08 x78, B0329+54 x11), 466 of them <= 10 min. The pulsar is in CB00:
  **CB00 alone is 0.15 TB (5,796 files)**; the 5.8 TB figure is all 40 beams, 39 of which point
  at empty sky beside the pulsar (re-measured 2026-09-26 — the first draft of this plan and of
  the staging request called 5.9 TB "CB00 only", a 40x overstatement).

| subset | size | stream time at 24-28 MB/s |
|---|---|---|
| whole release | 761.7 TB | 315-367 days — infeasible |
| obs containing a pulsar (all beams) | 230.7 TB | 95-111 days |
| **only the CBs containing a pulsar** | **7.5 TB** | 3.1-3.6 days |
| **pulsar-calibrator obs, CB00 only** | **0.15 TB** | ~1.5-2 h |
| pulsar-calibrator obs, all 40 beams (not needed) | 5.8 TB | ~2.5-2.8 days |
| online now, no staging | 3.9 TB | ~45 h |

**Cost:** stream -> process -> delete on the workstation is **$0** (download is free, 67 GB
free disk holds a wave). AWS adds nothing but speed for this; it is worth it only for the
Heimdall (CUDA-only) cross-check on a sample (a few dollars on spot). The earlier $170 estimate
assumed 2 TB and is superseded.

**Re-scoped deliverable:** (1) the recover-a-known on the ~486 pulsar-calibrator observations
(CB00 only, 0.15 TB)
— does the pure-PyTorch stack redetect B0531+21/B1933+16/B0950+08/B0329+54 across ~2 years of
heterogeneous 1-bit data?); (2) the blind single-pulse census of the 3,489 pulsar-field CBs
(7.5 TB) — RRAT/known-pulsar redetection statistics + candidate trains; (3) the Heimdall
cross-check (software note). The "archive-wide census" framing below is withdrawn.

## Context

Apertif Time-Domain DR2 (1,666 pointings, ~0.48 GB/pointing, 1-bit Stokes-I PSRFITS) was
FRB-searched with AMBER/ALERT (arXiv:2406.00482) but has no published archive-wide single-pulse
RRAT / known-pulsar census (ideas.md F24). This slice is the second real-data leg for the
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
