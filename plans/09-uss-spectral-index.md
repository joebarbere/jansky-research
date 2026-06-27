# Plan 09 — Ultra-steep-spectrum source hunt 🚀

> Context: second backlog slice, chosen for novel-finding potential. Scope: medium.

## Context

The first slice (FRB burst-statistics) validated the tooling by *reproducing* a known result. This
slice targets a genuine **discovery mode**: ultra-steep-spectrum (USS) selection is the classic way
high-redshift radio galaxies are found (De Breuck et al. 2000). Cross-matching two public radio
surveys at widely separated frequencies and flagging the steepest-spectrum sources can surface
candidates not yet catalogued as high-z radio galaxies — an honest "novel candidate" finding for an
amateur with public data.

## Deliverables

- `src/jansky_research/spectra.py` — two-point spectral index (with error propagation), source
  classification (uss / steep / flat / inverted), positional cross-match, `find_uss`, a VizieR
  `fetch_survey` (TGSS ADR1 147.5 MHz × NVSS 1.4 GHz; robust to NVSS sexagesimal coords), a NED
  `annotate_known` (novel-vs-known split), `run()` writing a candidate CSV + figures + metrics, and
  an offline `synthetic_field` fixture. Tested to the 85% floor.
- A real run over a high-latitude extragalactic field → `results/uss_candidates.csv` + figures.
- `survey/uss-findings.md` — the honest assessment: the USS candidate list, which are already
  known (NED), and which are genuinely uncatalogued candidates for follow-up.

## Approach

Reuse the established pattern (tested helper + synthetic fallback + offline tests + GATE-2 science
review). Primary spectral index from NVSS × TGSS (both recover extended flux); VLASS (3 GHz) is an
optional curvature check. **Caveats to surface:** the TGSS ADR1 flux-scale systematic
(Hurley-Walker 2017; keep the USS cut conservative at −1.3), resolution mismatch, and cross-match
confusion. The "finding" is explicitly a *candidate list requiring follow-up*, not a discovery.

## Verification

- `make cov` ≥85% on synthetic fixtures; `mypy` + `ruff` clean.
- Real run produces a sensible median α ≈ −0.7–0.8 (the radio-source norm) as a sanity check.
- **science-reviewer** gate on the method + the novel/known classification before the write-up.
