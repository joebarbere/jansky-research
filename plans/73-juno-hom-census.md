# 73 — Juno/Waves HOM occurrence census + moon-induced-emission incidence table (0.3–3 MHz)

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there)

## Context

The merged `junodam` slice built the full pipeline for the public Juno/Waves Estimated Flux
Density CDFs (doi:10.25935/6jg4-mk86): reader, background-threshold activity detection,
Horizons-based Juno-frame CML, Io phase, occurrence maps. Re-banding the same CDFs to
0.3–3 MHz opens hectometric (HOM) emission — an occurrence census plus a moon-induced-emission
incidence table, with the Louis+2023 moon-footprint encounter detections as positive controls
(fable-ideas F36). This space is crowded: the LESIA team owns the physics papers on this exact
dataset — so the slice is scoped narrowly and explicitly as the *occurrence-statistics
complement* to their event-physics work, not a competitor. Data risk ≈ 0 (same verified daily
CDFs, ~37 MB/day, no auth).

## Deliverables

- `src/jansky_research/junohom.py`: HOM-band (0.3–3 MHz) selection over the `junodam` CDF reader,
  `detect_active` re-banded, occurrence maps in (CML, Io/Europa/Ganymede phase) planes,
  `moon_incidence_table` (per-moon enhancement contrast + per-cell exposure),
  `louis2023_controls` (encounter-window recovery check), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/junohom/`; `survey/junohom-findings.md`; wiring.

## Approach

0. **GATE 0:** full-text pass on Louis+2023 and the recent LESIA HOM literature to confirm no
   published occurrence census from this public dataset covers the same ground; extract the
   exact encounter windows for the positive controls; re-verify the CDF URL pattern.
1. Tooling + synthetic recover-a-known: `junodam`'s `synthetic_orbit` fixture re-banded — an
   injected phase-organised HOM box and known moon-enhancement contrast round-trip.
2. Real leg: multi-month HOM occurrence census; per-moon phase-organisation contrasts with
   per-cell exposure; the incidence table as the headline product.
3. Positive-control leg: the Louis+2023 encounter detections must appear as active intervals.
4. GATE-2 (units caveat inherited from `junodam` — occurrence only, not flux physics; Juno-frame
   vantage/proximity effects; crowded-space framing checked) → paper (note-scale complement).

## Verification

Synthetic round-trip recovers the injected moon-phase organisation; Louis+2023 encounters
recovered as positive controls before any incidence claim; checks green; GATE-2 sign-off.

## Risks & mitigations

- **Crowded LESIA space** → scoped strictly to occurrence/incidence statistics; the GATE-0
  full-text pass re-checks the boundary; drop if their census exists.
- **Proximity/vantage effects at perijove** → perijove-day flagging inherited from `junodam`;
  the census is stated as vantage-dependent.
- **Multi-moon phase confusion** → per-moon planes with per-cell exposure masks; contrasts
  quoted with scramble-based nulls.
