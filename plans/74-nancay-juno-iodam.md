# 74 — Ground (Nançay JunoN) × Juno simultaneous Io-DAM census: two vantages, one emission

Status: 📋 planned (not started) — GATE 0 pending: file-level NDA/JunoN data access + full-text
novelty pass (the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

The merged `junodam` slice found that the classic ground-based Io-A/B/C/D (CML, Io-phase) boxes
do not organise Juno's sky the way they organise Nançay's — the vantage matters. The Nançay
Decameter Array's JunoN support campaign provides a public ground-based stream simultaneous
with the Juno/Waves public CDFs; a simultaneous Io-DAM occurrence census — same emission
windows, two vantages — extends that finding with an independent ground stream (fable-ideas
F37). The comparison is the point: which activity intervals are seen from the ground, from
orbit, or both, per (CML, Io-phase) cell. The Juno side is data-risk-free (verified `junodam`
CDFs); the ground side is the hard gate — file-level access to the NDA/JunoN public products
has not been live-verified.

## Deliverables

- `src/jansky_research/junon.py`: `fetch_junon` (NDA/JunoN product reader, network, `# pragma`),
  ground-side `detect_active` matched in band/threshold convention to `junodam`'s,
  `simultaneous_census` (joint occurrence in the (CML, Io-phase) plane with a
  ground-only/orbit-only/both classification per interval), visibility masking (Jupiter above
  the Nançay horizon), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/junon/`; `survey/junon-findings.md`; wiring.

## Approach

0. **GATE 0 (hard):** file-level access to the NDA/JunoN public products — locate the archive,
   fetch one real file, record format/cadence/calibration state; plus the standing full-text
   novelty pass (Nançay–Juno comparison papers exist — confirm the *occurrence-census* framing
   is unclaimed). No access → park.
1. Tooling + synthetic recover-a-known: a shared injected emission interval seen through two
   simulated vantage/visibility masks round-trips to the correct both/ground-only/orbit-only
   classification.
2. Real leg: a bounded simultaneous interval (weeks–months, sized by GATE-0 findings); joint
   census; Io-region contrast per vantage; overlap statistics with visibility-corrected exposure.
3. GATE-2 (sensitivity mismatch between a ground array and Juno/Waves; ionospheric cutoff on the
   ground side; beaming vs sensitivity degeneracy in "orbit-only" intervals) → paper.

## Verification

Synthetic two-vantage round-trip; ground side must re-find the canonical Nançay Io-A/B/C/D
organisation, and the Juno side must reproduce `junodam`'s vantage result, before any joint
claim; checks green; GATE-2 sign-off.

## Risks & mitigations

- **NDA/JunoN file access unverified** → GATE-0 is a hard stop; no synthetic-only paper.
- **Sensitivity mismatch masquerading as beaming** → report overlap statistics only above a
  matched effective threshold; degeneracy stated plainly.
- **LESIA proximity** (their instrument, their spacecraft) → occurrence-census scope, mirroring
  the `junohom` framing; drop if a joint census is already published.
