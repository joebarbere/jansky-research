# 68 — Fermi 4FGL-DR4 × ATNF: a from-scratch cross-match and vetted unassociated-source list

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

Fermi 4FGL-DR4 and the ATNF pulsar catalogue are both public and both continuously updated, but
associations propagate from pipeline lore; an independent from-scratch positional + Ė/d²
plausibility cross-match is a reproducibility exercise with a useful by-product (fable-ideas
F31): a small vetted list of known radio pulsars sitting inside *unassociated* 4FGL error
ellipses with gamma-ray-plausible spin-down flux. The hard validation anchor is built in — the
match machinery must recover all ~294 known gamma-ray pulsars before any new-candidate claim.
Catalogue plumbing and the Ė/d² axis reuse the merged `ppdot` slice directly. Data: 4FGL-DR4
FITS from FSSC, ATNF via `psrqpy`/web — both no-auth.

## Deliverables

- `src/jansky_research/fermipsr.py`: `fetch_4fgl` / `fetch_atnf` (network, `# pragma`),
  `ellipse_match` (positional match honouring per-source 95% error ellipses),
  `edot_flux_rank` (Ė/d² plausibility score, `ppdot` reuse), `vet_candidates` (known-association
  exclusion + chance-coincidence FAP from scrambled catalogues), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/fermipsr/`; `survey/fermipsr-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text pass on the 4FGL-DR4 paper and the latest Fermi pulsar-census papers to
   confirm no published independent cross-match supersedes this; verify the DR4 FITS URL and the
   current ATNF snapshot; record the current known gamma-ray-pulsar count.
1. Tooling + synthetic recover-a-known: injected pulsars inside mock ellipses recovered at the
   stated FAP; scrambled-position nulls behave.
2. Anchor leg: the cross-match must recover all ~294 known gamma-ray pulsars (misses are bugs).
3. Real leg: match against unassociated 4FGL sources, rank by Ė/d², vet, and publish the short
   candidate table with per-source chance-coincidence probabilities.
4. GATE-2 (ellipse-shape handling; ATNF distance-model systematics in d²; chance coincidences in
   the Galactic plane) → paper (note-scale; the vetted table is the deliverable).

## Verification

All ~294 known gamma-ray pulsars recovered before any candidate claim; scramble-based FAPs
quoted per candidate; checks green; GATE-2 sign-off.

## Risks & mitigations

- **The Fermi team's own association pipeline is thorough** → the value is independence and a
  ranked, reproducible table, not a claimed discovery; frame as reproducibility + shortlist.
- **DM-distance systematics dominate d²** → quote both YMW16 and NE2001 rankings.
- **Plane crowding inflates matches** → per-candidate FAP from footprint-preserving scrambles.
