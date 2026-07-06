# 55 — Polarized-flux dipole from SPICE-RACS DR2 (the open axis of the dipole anomaly)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — specifically: confirm
no polarized-flux dipole with a modern catalogue has landed since Tiwari & Jain 2013/15, with a
same-week ADS re-search (the dipole field moves fast)

## Context

The cosmic-dipole-anomaly literature is total-intensity source-count/flux based and now crowded
(arXiv:2509.16732 PRL joint NVSS+RACS+LoTSS-DR2 at 5.4σ; arXiv:2509.18689; RMP colloquium
arXiv:2505.23526 — the Corrections fence). The **polarized-flux dipole** is one of the two open
differentiator axes (fable-ideas F18): last measured with NVSS by Tiwari & Jain 2013/15, never
with a modern RM/polarization catalogue. SPICE-RACS DR2 (arXiv:2605.16917; ~2.5–3.4×10⁵
polarized sources over 87.5% of sky; CSIRO DAP `csiro:64891`, ~5 GB, no auth) is already
fetched and used by the merged `rmstructure` slice, so data risk ≈ 0. This is a sibling of
plan 38 (F1, the RM dipole) — sequence it after that slice; same data, same HEALPix/dipole-fit/
footprint-scramble muscles.

## Deliverables

- `src/jansky_research/poldipole.py`: `load_dr2` (reuse the `rmstructure` fetch, `# pragma`),
  `polflux_sample` (S/N + leakage-aware polarized-flux cuts), `dipole_fit` (HEALPix-binned
  likelihood dipole in polarized flux/counts), `footprint_scramble` (RA-scramble null
  distribution preserving the Dec≤+49° footprint), `nvss_anchor` (reproduce the Tiwari & Jain
  NVSS-era measurement), `synthetic_dipole_catalogue`, `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/poldipole/`; `survey/poldipole-findings.md`; wiring.

## Approach

0. **GATE 0:** standing full-text novelty pass (Tiwari & Jain 2013/15; arXiv:2509.16732,
   2509.18689, 2505.23526) + same-week ADS re-search for any polarized-flux dipole claim;
   confirm the DR2 table on disk carries the polarized-flux columns needed.
1. Tooling + synthetic recover-a-known: inject a known dipole into a mock polarized-source
   catalogue with the real DR2 footprint and noise; confirm amplitude + direction recovery
   (mirrors plan 38's mock and `rmstructure`'s synthetic-screen pattern).
2. Recover-a-known on real data: reproduce the NVSS-era Tiwari & Jain result from the public
   NVSS RM/polarization catalogue before touching DR2.
3. DR2 leg: polarized-flux sample after S/N cuts → HEALPix bins → dipole fit; nulls from
   footprint-preserving RA scrambles (load-bearing — the footprint is the dominant systematic);
   compare direction to the CMB dipole and the source-count dipole.
4. GATE-2 science review naming the fable-ideas caveats: sparse counts after S/N cuts,
   scan-pattern leakage systematics; honest isotropy-test framing (no over-interpretation of
   amplitude without a kinematic expectation).
5. Paper: the first modern-catalogue polarized-flux dipole measurement (or honest null).

## Verification

Injected mock dipole recovered in amplitude + direction on the real footprint; the NVSS-era
Tiwari & Jain result reproduced first; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Sparse counts after S/N cuts** → report dipole stability across a ladder of cuts; quote
  per-cut source counts and widen errors honestly rather than cherry-picking one threshold.
- **Scan-pattern leakage systematics** → leakage-aware sample cuts + the RA-scramble null is
  load-bearing; check the fitted direction against the RACS scan-pattern axes explicitly.
