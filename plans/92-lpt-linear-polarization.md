# 92 — `lptlin`: is there linear polarization at LPT positions?

Status: 🚧 GATE 0 complete 2026-08-21, and it **overturned this plan's own premise**. The
scope reduction below was generalised from a query at ONE position and is wrong: 14 of 16 LPT
positions have Q/U in CASDA, VAST does produce Q/U continuum images, and 4 of 7 detected
pulses have Q and U in the same observation. One linear detection made: ASKAP J1832-0911 at
L/I = 10.9% (33 sigma). See survey/lptlin-findings.md. Novelty pass in flight.

Was: planned 2026-08-20, **scope reduced by a feasibility check before writing**.
Sourced from Rose et al. 2026 (Nature Astron. 10, 1166).

## Context

Rose et al. describe the bursts of ASKAP J174508.9−505149 as **elliptically** polarized —
i.e. carrying linear as well as circular power. The `lptv` census is Stokes-V only and says
so explicitly: its non-detections are "limits on circularly polarized flux, not on
fractional polarization". A linear-polarization counterpart would complete that census.

## The feasibility check, and what it killed

The original idea was "forced Q/U photometry at every LPT position, like the V census".
**That is not possible with RACS.** A CASDA ObsCore query at an LPT position (2026-08-20)
returned 292 products with this Stokes breakdown:

| product | i | q | u | v |
|---|---|---|---|---|
| image | 70 | **1** | **1** | 49 |
| residual | 23 | 1 | 1 | 19 |
| weights | 33 | 1 | 1 | 19 |

and every one of those q/u files is **EMU**, not RACS (`image.q.EMU_…`). RACS continuum
products at this position are I and V only. So the uniform, all-sky-south census this was
meant to be does not exist to be done.

Note also that `data/spice-racs.dr2.fits` (9.6 GB, on disk) is a **compact-source RM
catalogue**, not images: LPTs are transient and faint and will overwhelmingly be absent from
it. It is useful here as a foreground/leakage reference, not as a measurement of LPTs.

## What is left, and whether it is worth doing

A much smaller slice: **linear polarization for whichever LPTs fall inside EMU coverage**,
at whichever epochs EMU observed them. This is worth planning only if GATE 0 answers three
questions in the affirmative:

1. **Footprint.** How many of the 16 catalogued LPTs have EMU q/u coverage at all? If the
   answer is one or two, this is a per-source note at best, not a census — say so and size
   the write-up accordingly (RNAAS, if anything).
2. **Coincidence.** Does any EMU epoch coincide with a *known pulse*? A linear-polarization
   limit taken while the source is off is not a limit on the burst's polarization; it is a
   limit on nothing. This is the same trap as quoting |V|/I at a position where Stokes I is
   consistent with zero — a mistake `lptv` explicitly avoids and this slice must not
   reintroduce.
3. **Novelty.** Several LPT discovery papers report linear polarization from their own
   follow-up (MeerKAT, ATCA, MWA). Establish what is already published per source before
   claiming a first.

If (1) and (2) both come back empty — plausible — the honest outcome is a recorded negative:
"the archival imaging cannot constrain LPT linear polarization, and here is why", which is
still worth a paragraph in `survey/` and saves the next person the same search.

## Deliverables (conditional on GATE 0)

- EMU q/u forced photometry reusing `stokesv.fetch_racs_cutout` (its filename mask is
  parameterised by Stokes; EMU products follow the same convention).
- `results/lptlin_metrics.json` with per-source coverage, epoch-vs-pulse coincidence, and
  either measurements or explicitly-scoped limits.
- `survey/lptlin-findings.md` — including the negative, if that is the answer.

## Related

`lptv`, `stokesv`, `rmstructure` (RM/leakage handling), plan 90, plan 91.
