# 52 — SSL/anomaly sweep of a survey nobody has swept (RACS continuum or LoTSS DR3)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — specifically: read
STRADAViT (arXiv:2603.29660) + the LoTSS DR3 paper full-text for "in prep" claims, and check the
unverified HF checkpoint `ISSA-ML/stradavit-base`

## Context

SSL/SOM/anomaly sweeps exist for MGCLS, EMU-pilot, LoTSS DR1/DR2, and partial VLASS-QL — but
(per the fable-ideas F15 verification, July 2026) **none for RACS continuum and none for LoTSS
DR3** (a ~5-month window). The recipe is established and ROCm-safe: BYOL/DINO features →
Astronomaly-Protege active anomaly ranking (arXiv:2411.04188, 2602.15930); arXiv:2409.11175
shows generic DINOv2 ViTs already transfer to radio morphology at F1 0.72–0.88, so a
linear-probe leg de-risks the pretrain. Data: LoTSS DR3 cutout API (no auth) or RACS cutouts
via CASDA (verified auth; `stokesv` cutout muscle). Fence from Corrections: RFI-benchmark-style
"new architecture" angles are saturated — this is a *sweep of an unswept survey*, not an
architecture paper. This is the long-GPU-run slice (owner-approved; weeks OK, pure
PyTorch/ROCm on the 16 GiB RX 7600 XT). Pairs with the F26 artifact census for de-artifacting.

## Deliverables

- `src/jansky_research/anomsweep.py`: `fetch_cutouts` (DR3 API / CASDA, `# pragma`),
  `probe_dinov2` (linear-probe features first), `pretrain_ssl` (BYOL/DINO fallback, pure
  PyTorch), `extract_features`, `protege_rank` (active anomaly loop), `known_anchor_check`
  (ORCs/rings/GRGs in-footprint), `synthetic_morphology_set`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/anomsweep/`; `survey/anomsweep-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text read of STRADAViT (arXiv:2603.29660) and the LoTSS DR3 paper for "in
   prep" sweep claims; check the unverified HF checkpoint `ISSA-ML/stradavit-base`; live-verify
   the DR3 cutout API (or scope to RACS via CASDA); standing full-text novelty pass.
1. Tooling + synthetic recover-a-known: inject synthetic rare morphologies (rings, doubles with
   odd axis ratios) into a mock cutout set; the feature+Protege ranking must surface them.
2. Feature leg: DINOv2 linear probe first (cheap, per arXiv:2409.11175); escalate to a full
   SSL pretrain (BYOL/DINO) on ~10⁵–10⁶ cutouts only if probe features underperform — this is
   the owner-approved multi-week GPU run (~275 GB disk budget governs cutout count).
3. Protege active-anomaly loop + human vetting of the top-N; de-artifact the list (cross-check
   offset-vector ghost/duplicate signatures — the F26 pairing).
4. GATE-2 science review: teams-may-be-mid-flight scoop check, artifact contamination of the
   top-N, honest statement of what "anomalous" means for the chosen feature space.
5. Paper: ranked rare-morphology candidate list + released feature model + training recipe.

## Verification

Known ORCs/rings/GRGs in-footprint must rank highly; injected synthetic morphologies recovered
in the offline round-trip; checks green; GATE-2 sign-off.

## Risks & mitigations

- **LOFAR/ASKAP teams could be mid-flight** → GATE-0 full-text pass + same-week re-search; ship
  a bounded first sweep fast rather than the perfect one.
- **Morphology novelty requires careful de-artifacting** → artifact vetting step is mandatory;
  pairs with the F26 artifact census, whose flag list feeds this slice.
- **ROCm pretrain instability on gfx1102** → linear-probe leg first (arXiv:2409.11175 shows it
  may suffice); pin torch versions per the `torchfdmt-findings.md` pattern.
