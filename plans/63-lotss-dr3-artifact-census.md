# 63 — LoTSS DR3 artifact archaeology: exhaustive offset-vector ghost/duplicate census

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — read the DR3 QA
section first and scope this census strictly against what it already covers

## Context

Mega-catalogue QA sections sample; nobody publishes exhaustive small-separation pair statistics
over an entire 13.7M-source catalogue (fable-ideas F26). This slice runs a GPU-chunked all-pairs
search within 10′ over the full LoTSS DR3 source list, clusters the offset vectors in two frames
— relative to the nearest bright source, and in the mosaic-tile frame — and identifies ghost
families as repeated offset/PA excesses (calibration ghosts, mosaic-overlap duplicates,
deblending fragments). Deliverable shape: "DR3 contains <X% duplicates above flux S" plus a flag
list, citable by every DR3 user — including this repo's own F1 RM-dipole lineage (plans/38) and
the F15 SSL/anomaly sweep, whose morphology-novelty claims need exactly this de-artifacting.
Pure-torch all-pairs on the 16 GiB RX 7600 XT; this is the compute-for-days idea the workstation
profile explicitly invites. It exploits compute nobody bothers to spend, not new data.

## Deliverables

- `src/jansky_research/dr3ghosts.py`: `fetch_dr3_catalogue` (lofar-surveys.org, `# pragma`),
  `allpairs_within` (GPU-chunked all-pairs ≤10′ over 13.7M sources, pure torch, spilling
  partial results to disk), `offset_vectors` (separation + PA in bright-source and mosaic-tile
  frames), `cluster_ghost_families` (repeated offset/PA excess detection vs an isotropic-pair
  null), `duplicate_fraction` (per flux bin, with errors), `flag_list` (released artifact),
  `inject_duplicates` (planted duplicates + documented deblending modes → recovery),
  `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/dr3ghosts/`; `survey/dr3ghosts-findings.md`; wiring.

## Approach

0. GATE 0: read the DR3 release paper's QA section in full and scope the census against it (only
   claim what their sampling did not cover); verify the DR3 catalogue download URL + columns;
   ADS check for any independent DR3 artifact paper.
1. Tooling + synthetic recover-a-known: plant duplicates and known deblending-mode pairs into a
   realistic mock catalogue; the clustering must recover the planted families and the isotropic
   null must stay clean.
2. Real leg: chunked all-pairs over the full catalogue (multi-day GPU job, checkpointed);
   offset-vector clustering in both frames; ghost-family identification; duplicate-fraction
   table per flux bin; the flag list.
3. GATE-2 science review: physical-pair vs artifact-pair separation honesty (real doubles exist
   at all separations — the excess-over-null framing is mandatory), QA-section overlap statement.
4. Paper: the census + released flag list, framed as a QA companion to DR3, not a criticism.

## Verification

Injected duplicates and documented deblending modes recovered at stated completeness; isotropic
mock yields no false families; checks green; GATE-2 sign-off.

## Risks & mitigations

- **DR3 QA section may already cover part of this** → GATE-0 reads it first; the plan scopes to
  the exhaustive/all-pairs remainder and cites their sampling explicitly.
- **Real close pairs masquerade as artifacts** → everything is excess-over-isotropic-null; the
  flag list marks "consistent with ghost family", never "is an artifact", per source.
- **16 GiB VRAM vs 13.7M sources** → chunked pair blocks with disk spill; checkpointed so the
  multi-day run resumes cleanly.
