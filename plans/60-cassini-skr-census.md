# 60 — Cassini SKR: post-2013 dual-period census + Grand-Finale proximity law

Status: 📋 planned (not started) — GATE 0 pending: full-text novelty pass + data-URL verification
(the fable-ideas scan ran egress-blocked; see the standing caveat there) — pin the pre-2008 PDS
volume ID for CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0 before any code

## Context

Fischer+2015's Saturn Kilometric Radiation dual-period tracking stops in early 2013; nobody has
run the RPWS-flux dual-period census through the 2017 proximal (Grand Finale) orbits, and nobody
has asked whether the `junodam` ~180× proximity occurrence law holds for SKR (fable-ideas F23 —
a direct port of the merged `junodam` census pattern to Saturn). Note the fence: the Saturn SED
(lightning) census is closed (Fischer et al. 2025, 10.1029/2024JA033560) — this slice is SKR
occurrence/periodicity, not lightning. Data: PDS-PPI RPWS 60-s key parameters,
`CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0` (the pre-2008 volume ID is the GATE-0 pin; later
coverage volumes to be identified alongside), plus JPL Horizons for Cassini–Saturn geometry.
Tooling: `junodam` (occurrence census, background+kσ detection, Horizons batching) + `frbperiod`
(Lomb-Scargle) reuse. Validation is built in: the known ~10.8/10.6 h north/south period split
must be re-derived before anything new is claimed.

## Deliverables

- `src/jansky_research/skr.py`: `fetch_rpws_key60s` (PDS-PPI volume download, `# pragma`),
  `read_key_params` (60-s key-parameter parser → SKR-band flux series), `detect_skr`
  (background + kσ per band, `junodam` pattern), `dual_period_ls` (sliding-window Lomb-Scargle
  via `frbperiod` reuse, tracking both periods), `proximity_duty_cycle` (range-binned occurrence
  vs Saturn distance), `magnetic_latitude_weight` (dipole-field latitude correction for the
  proximity law), `fetch_geometry` (Horizons batch, `# pragma`), `synthetic_skr` (injected
  dual-period + proximity trend → recovery), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/skr/`; `survey/skr-findings.md`; wiring.

## Approach

0. GATE 0: pin the pre-2008 PDS volume ID for CO-V/E/J/S/SS-RPWS-4-SUMM-KEY60S-V1.0 and the
   volumes covering 2013–2017; confirm direct download + format docs; full-text pass on
   Fischer+2015 and the Fischer 2025 SED paper to confirm the post-2013 SKR census is unclaimed.
1. Tooling + synthetic recover-a-known: inject two close periods and a range-dependent occurrence
   trend into a synthetic key-parameter series; both must be recovered.
2. Real leg A (anchor): re-derive the known ~10.8/10.6 h north/south SKR period split on the
   pre-2013 overlap with Fischer+2015 before trusting the pipeline further.
3. Real leg B: extend the dual-period track through 2017; range-binned duty cycles across the
   proximal orbits with magnetic-latitude weighting; compare to the `junodam` proximity law.
4. GATE-2 science review: latitude-weighting model dependence, key-parameter (not full-res)
   caveat, visibility vs intrinsic occurrence. 5. Paper: the post-2013 census + proximity result.

## Verification

Pipeline re-derives the published 10.8/10.6 h period split on the Fischer+2015 overlap first;
synthetic round-trip recovers injected periods and proximity trend; checks green; GATE-2 sign-off.

## Risks & mitigations

- **SKR "proximity" needs magnetic-latitude weighting** — a real modelling step, not a rescaling:
  SKR beaming means occurrence vs range is confounded with latitude coverage; the weighting model
  and its assumptions get their own GATE-2 scrutiny and a stated-model-dependence caveat.
- **Pre-2008 volume archaeology** (PDS3 layouts, volume splits) → GATE-0 pins IDs + format docs;
  keep the parser tested against a vendored sample block, per the `vgpra` pattern.
