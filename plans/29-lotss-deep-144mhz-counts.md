# 29 — Deep 144 MHz source counts from LoTSS DR2 (the sub-mJy star-forming-galaxy upturn)

Status: 📋 planned (not started) — a *larger* effort: new survey + a non-trivial completeness/area gate

## Context

The `sourcecounts` slice (#27) recovered the canonical **1.4 GHz** Euclidean-normalised counts from
NVSS and reproduced the Hopkins 2003 fit over 3.5 mJy–1 Jy — but it explicitly could **not** reach the
faint regime where the normalised counts turn up again as **star-forming galaxies** (SFGs) enter the
population below ~1 mJy. The LOFAR Two-metre Sky Survey DR2 (LoTSS DR2; Shimwell et al. 2022, A&A 659,
A1) is the deepest large-area low-frequency survey (144 MHz, 6″, ~0.8 mJy over 5634 deg², 4.4M sources)
and reaches *well into* that sub-mJy regime. This slice carries the #27 pipeline to 144 MHz and
**measures the sub-mJy upturn** that the 1.4 GHz NVSS slice flagged as out of reach.

It is the natural companion to #27 and reuses its tooling near-verbatim; the genuinely-new science is
the faint-end shape (the SFG upturn), at a frequency where the source population is different
(steep-spectrum-dominated, with the SFG/AGN transition moving in flux).

**Data access is verified.** LoTSS DR2 is on VizieR (`J/A+A/659/A1`); the integrated flux is
`SpeakTot` (= `Total_flux`, **mJy**, at 144 MHz), `Speak` is the peak (mJy/beam), `Maj`/`SCode` give
compactness/structure. Cone queries with explicit columns return in <1 s (the `columns=["**"]` form is
the only slow path — avoid it). Sample fluxes are ~0.5 mJy, confirming the sub-mJy reach.

## Deliverables

- `src/jansky_research/lotsscounts.py` — pure-NumPy/astropy tooling, reusing
  `jansky_research.sourcecounts` (and thus `jansky.sourcecounts`) wholesale:
  - `fetch_lotss_region` — VizieR `J/A+A/659/A1`, `SpeakTot` (mJy → Jy), restricted to a region fully
    inside a LoTSS field with a defensible effective area (`# pragma: no cover`).
  - `lofar144_reference` — a published 144/150 MHz Euclidean-normalised count reference (Siewert et al.
    2020 LoTSS DR1; cross-check vs Mandal et al. 2021 deep fields / Intema et al. 2017 TGSS 150 MHz),
    the analogue of #27's `hopkins2003_counts`, used only over its stated validity range.
  - `compute_counts` (reuse #27's, parameterised by the reference) + `synthetic_sky` drawn from the
    144 MHz reference (round-trips to ratio 1).
  - `run(offline=...)`, `_figure` (Euclidean counts + reference + the marked sub-mJy upturn),
    `_write_macros`, `_main`.
- `tests/test_lotsscounts.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry entry `lotss-dr2`; `papers/lotsscounts/` (AASTeX); `survey/lotsscounts-findings.md`.

## Approach

1. **Tooling + synthetic recover-a-known.** Draw a flux-limited sky from the 144 MHz reference and
   confirm the (reused) binning/normalisation round-trips to ratio ≈ 1; the synthetic must include the
   faint-end upturn so the fixture exercises it.
2. **Effective-area gate (the #27 caveat, harder here).** Unlike an NVSS cone, LoTSS DR2 is a mosaic of
   pointings with position-dependent rms/completeness. Restrict to a sky box fully inside one LoTSS
   field, adopt a conservative flux cut (≳ a few × the local 0.8 mJy limit) so completeness ≈ 1, and
   use the corresponding effective solid angle. Surface this as the dominant systematic.
3. **Real recover-a-known.** Build the 144 MHz Euclidean-normalised counts over the region; **recover
   the published LoTSS/150 MHz counts** including the **sub-mJy SFG upturn** (the new feature beyond
   #27); report the ratio to the reference and the flux of the upturn.
4. **GATE-2 science review** before the write-up: completeness/area handling, the LoTSS flux-scale
   (~10–20%) tie, multi-component sources at 6″, resolution bias, and the reference-curve choice must
   all survive; the upturn must be a real population feature, not a completeness artefact.
5. **Write-up** `papers/lotsscounts/` — a reproducibility/method contribution that extends #27 to the
   faint, low-frequency regime; not a new measurement.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` round-trips the 144 MHz reference (ratio ≈ 1) and shows the injected upturn.
- (Real-data) `python -m jansky_research.lotsscounts ...` reproduces the counts, the figure with the
  reference + marked upturn, and the macros from VizieR LoTSS DR2; the recovered counts (and the
  upturn) match the published 144/150 MHz counts; GATE-2 sign-off before the write-up.

## Alternative scope (if the counts angle proves too close to #27)

A **northern measured-turnover peaked catalogue**: LoTSS 144 MHz (a *detection*, not the TGSS upper
limit that crippled the original `peaked` slice #14) + WENSS 325 MHz + NVSS 1.4 GHz + VLASS 3 GHz →
fit `southern.fit_log_parabola` to *measure* the turnover frequency for northern GPS/CSS candidates,
mirroring how `southern` (#16) fixed this in the south with GLEAM-X. Same LoTSS data access; reuses
`southern`/`spectra` instead of `sourcecounts`; recover-a-known = known northern GPS sources with a
measured (not bounded) ν_pk.

## Risks & mitigations

- **Completeness/effective area (highest) →** region fully inside one LoTSS field + conservative flux
  cut so completeness ≈ 1; report the area assumption as the dominant systematic; cross-check the count
  normalisation against the published LoTSS counts.
- **Reference-curve choice →** adopt one published 144/150 MHz count (Siewert 2020) and cross-check vs a
  second (Mandal 2021 / Intema 2017); treat as a recover-a-known, not a new measurement.
- **Flux scale + multi-component sources →** ~10–20% LoTSS scale tie folded into errors; `SCode`/`Maj`
  cut or component-merge for the faintest bins; cite the LoTSS DR2 systematics.
- **Too close to #27 →** lead on the *new* faint-end upturn (the part NVSS can't reach), or switch to
  the measured-turnover alternative scope above.
