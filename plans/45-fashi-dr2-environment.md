# 45 — FASHI DR2 (156,411 HI sources): the environment statistics nobody has run

Status: ✅ done (DR1 first leg) — GATE 0 2026-07-08: **DR2 not public yet** (catalogue links 404;
release ~Aug 2026), so the DR1 leg on VizieR `J/other/SCPMA/67.19511/table2` (41,741 sources) —
the `rmstructure` DR1-while-DR2-embargoed pattern; DR2 swap is one line. Headline: **void HIMF
knee suppressed −0.10 dex vs walls**, an independent FAST confirmation of the ALFALFA void HIMF
(Moorman+2014); group-member knee +0.19 dex, survivor-biased (stated). Scope corrections forced
by FASHI lacking stellar masses/optical diameters: dropped "gas fraction at fixed M*" and
"deficiency vs radius" (selection-biased), replaced by the cleaner group/field HIMF split.
Absolute faint-end slope steeper than the published global (simple 1/Vmax, no completeness
function) — relative offsets robust, stated honestly. See survey/fashienv-findings.md.

## Context

FASHI DR2 dropped ~June 2026 (arXiv:2606.31539; 156,411 sources, 19,500 deg², z<0.09 — 4× DR1)
with only a global Schechter HIMF published. The single environment paper (arXiv:2510.22902)
used 230 DR1-era group galaxies (~3% of DR2 volume, no cluster-scale densities). Open questions
this slice runs: environment-split HIMF, HI-deficiency vs clustercentric radius (the
Solanes+2001 curve at 40× the sample), and void-vs-wall gas fractions (pre-FASHI ALFALFA-era
only). Data: the FASHI DR2 table (VizieR / China-VO mirror — GATE-0 confirms the ID),
SDSS/DESI group catalogues (Tempel+2017 / Lim+2017 on VizieR), and the Douglass+2023 void
catalogues (VizieR, three void-finders). Extends the `hi.py` and catalogue-cross-match muscle
memory from merged slices; the biggest fresh-data surface in fable-ideas.md. CPU, days.

## Deliverables

- `src/jansky_research/fashienv.py`: `fetch_fashi_dr2` + `fetch_group_catalogues` +
  `fetch_void_catalogues` (# pragma), `crossmatch_environment` (position+velocity →
  clustercentric radius / local-density terciles / void membership), `censored_gas_fractions`
  (survival-analysis statistics for upper limits), `schechter_per_bin`, `deficiency_profile`,
  `flag_blends` (FAST 2.9′ beam confusion at low z), `synthetic_environment_catalogue`
  (injected deficiency profile + HIMF split → recovered), `run/_figure/_write_macros/_main`.
- Tests to the 85% floor; `papers/fashienv/`; `survey/fashienv-findings.md`; wiring.

## Approach

0. GATE 0: full-text pass on arXiv:2606.31539 and 2510.22902 (plus an ADS sweep for any DR2
   environment paper in flight) to confirm the gap; confirm the FASHI DR2 catalogue ID on
   VizieR or the China-VO mirror and that the group/void catalogues resolve.
1. Tooling + synthetic recover-a-known: build a mock catalogue with an injected HI-deficiency
   profile, an injected environment-split HIMF, and censoring; the cross-match + survival
   statistics + Schechter fits recover all three.
2. Real leg: cross-match FASHI DR2 to groups/clusters and to all three void-finder catalogues;
   censored gas-fraction statistics and Schechter fits per environment bin; report results per
   void-finder separately (algorithm dependence is a known trap); blend flags carried through.
3. GATE-2 science review: beam-confusion contamination at low z, group-catalogue membership
   errors, censoring-model sensitivity, void-finder disagreement.
4. Paper: environment-split HIMF + deficiency profile + void/wall gas fractions at DR2 scale.

## Verification

Anchors from fable-ideas.md before any new claim: reproduce the declining HI-deficiency profile
inside R200 for a well-sampled cluster, and the ALFALFA-era "void galaxies are gas-richer at
fixed M*" result; synthetic round-trip passes; checks green; GATE-2 sign-off.

## Risks & mitigations

- **FAST 2.9′ beam confusion at low z** → flag blends explicitly and test result stability with
  blends excluded.
- **The FASHI team is fast** → they have signalled global-HIMF interests, not cluster work;
  GATE-0 re-checks for in-flight environment papers; ship bounded.
- **Void-finder algorithm dependence** → report per void-finder (all three Douglass+2023
  catalogues); never quote a single-finder void result as the headline.
