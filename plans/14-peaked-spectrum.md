# 14 — Peaked-spectrum (GPS/CSS) radio sources via 3-frequency spectral curvature

Status: ✅ done (real-data run + two validations + paper drafted; `papers/peaked/main.tex`)

## Context

The ultra-steep-spectrum slice (`spectra`, `survey/uss-findings.md`) ended on a clear note: a
two-point index at 150 MHz vs 1.4 GHz cannot tell a genuine spectral turnover from a TGSS flux-scale
artefact — "confirming a turnover needs a third frequency (e.g. VLASS 3 GHz)." This slice is that
follow-up. With **three** frequencies — TGSS 150 MHz, NVSS 1.4 GHz, VLASS 3 GHz — we measure the
spectral **curvature** and select **peaked-spectrum** sources (GPS/CSS candidates: compact, young
radio AGN whose spectrum rises then falls, peaking in the ~0.1–3 GHz band). Curvature (comparing two
indices) is far more robust to a constant TGSS scale offset than a single steep cut, and the
high-frequency index (NVSS→VLASS) is TGSS-independent.

It is also the **maximal-reuse** slice: it composes the two most tooling-rich modules rather than
adding much new code.

## Reuse (the point of this slice)

- `spectra.spectral_index` — the two-point index + error (used twice, for $\alpha_\mathrm{low}$ and
  $\alpha_\mathrm{high}$).
- `spectra.crossmatch`, `spectra.fetch_survey` — TGSS/NVSS cone search + positional matching.
- `spectra.annotate_known` / `vlass.vet_candidates` — NED/SIMBAD cross-check (is a candidate a known
  GPS/CSS/QSO?).
- `vlass._fetch_e1_tap` — VLASS 3 GHz peak fluxes (VizieR TAP).
- `vlass.eta_metric` / `vlass.v_metric` — flag candidates that are *variable* across VLASS epochs
  (a flat/inverted **variable** source is a blazar, not a GPS source).

## Deliverables

- `src/jansky_research/peaked.py` — the new, tested logic is small:
  - `two_point_indices` — $\alpha_\mathrm{low}$ (150→1400 MHz) and $\alpha_\mathrm{high}$
    (1400→3000 MHz) via `spectra.spectral_index`.
  - `classify_sed` — peaked / steep / flat / inverted from the two indices.
  - `peak_frequency` — parabolic log-log fit → turnover frequency + concavity.
  - `find_peaked` — cross-match the three surveys, compute indices, classify, return the table.
  - `synthetic_field` — a 3-survey offline fixture with injected peaked/steep/flat SEDs.
  - `run(offline=...)` — synthetic offline; real TGSS×NVSS×VLASS cone online, with NED/SIMBAD vetting
    and a VLASS-variability blazar flag. Writes results + figures + macros.
- `tests/test_peaked.py` — synthetic-fixture tests to the 85% floor; no network.

## Approach

1. **Tooling (this step).** The curvature classification + 3-survey orchestration, validated on a
   synthetic field that recovers the injected peaked sources with low contamination.
2. **Real data (done).** A 2° cone at RA 180°, Dec +30° (the USS field) exposed two systematics and
   their fixes — TGSS depth (use it as an **upper limit**, since peaked sources are faint at 150 MHz)
   and NVSS$\to$VLASS **resolution mismatch** (a floor on $\alpha_\mathrm{high}$ rejects 110 extended
   artefacts). Result: 6 peaked candidates, all VLASS-steady (GPS-like), three uncatalogued; one is a
   known *Fermi* BL Lac that the variability flag missed but SIMBAD caught. See
   `survey/peaked-findings.md`. (`find_peaked` redesigned to the upper-limit + resolution-floor
   method; synthetic fixture updated; offline recovery ~96%.)
3. **Validation (done).** Two complementary tests pin the selection function. `validate_known`
   (Callingham et al. 2017, GLEAM-selected MHz-peaked): 0/81 false positives among sources whose
   turnover is below the 150 MHz floor — correct **purity**. `validate_hfp` (Dallacasa et al. 2000
   bright HFP sample): **98/98 = 100%** recovered as optically-thick-rising and **55%** flagged
   `ghz_peaked` (still rising at 3 GHz) — a clean **recover-a-known**. Classical GHz-GPS no longer
   "appear in the `inverted` class" as a caveat; they are the dedicated `ghz_peaked` class, and the
   three SED classes map cleanly onto the MHz-peaked / intermediate / GHz-peaked turnover regimes. See
   `survey/peaked-findings.md`.
4. **GATE-2** before write-up — candidates survive the flux-scale, resolution, and variability caveats.
4. **Write-up** as `papers/peaked/` once the real run + review are done.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers ≳80% of injected peaked sources with low false-positive contamination.
- `classify_sed` and `peak_frequency` match hand-computed values on textbook SEDs.
- (Real-data, later) recovers known GPS/CSS sources; GATE-2 sign-off before write-up.
