# 27 — 1.4 GHz Euclidean-normalised source counts from NVSS

Status: ✅ done (tooling + synthetic + real recover-a-known + GATE-2 + paper)

## Context

The differential source count dN/dS — radio sources per unit flux per unit sky — is one of the oldest
cosmological tests in radio astronomy: a static Euclidean universe gives dN/dS ∝ S^−5/2, and the
departures of the real counts from that law trace the cosmological evolution of the radio-source
population (Condon 1984 first measured the sub-Euclidean 1.4 GHz slope). This slice is the first to
exercise the untouched `jansky.sourcecounts` helpers, building the Euclidean-normalised count from a
public NVSS region and comparing it to the canonical 1.4 GHz reference (Hopkins et al. 2003). It is a
clean, deterministic recover-a-known: NVSS *should* reproduce the published counts, and it does.

## Deliverables

- `src/jansky_research/sourcecounts.py` — pure-NumPy/astropy tooling:
  - `hopkins2003_counts` — the published 6th-order polynomial reference (S^2.5 dN/dS, Jy^1.5 sr^-1),
    used only over its 0.05–1000 mJy validity range.
  - `compute_counts` — completeness cut, log binning, dN/dS per steradian (reusing
    `jansky.sourcecounts.differential_counts` / `euclidean_normalised_counts` / `count_slope`), the
    differential slope, and the Hopkins ratio — both restricted to N≥5 bins below 1 Jy.
  - `synthetic_sky` — offline fixture drawn from the Hopkins differential count (round-trips to ratio 1).
  - `fetch_nvss_region` (VizieR VIII/65, `# pragma: no cover`), `run`/`_figure`/`_write_macros`/`_main`.
- `tests/test_sourcecounts.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry entry `nvss`; `papers/sourcecounts/` AASTeX; `survey/sourcecounts-findings.md`.

## Approach

1. **Tooling + synthetic recover-a-known.** Draw a flux-limited sky from the Hopkins differential count
   and confirm the binning/normalisation round-trips to ratio ≈ 1.00 (0.03 dex scatter).
2. **Real run (8° NVSS cone at 180°, +30°, b ≈ +80°).** 7428 sources > 3.5 mJy over 0.0611 sr; the
   Euclidean-normalised counts match Hopkins 2003 over 3.5 mJy–1 Jy (10 bins) at **ratio 1.021,
   0.061 dex scatter**, with a sub-Euclidean differential slope **−1.91**.
3. **GATE-2 science review** (science-reviewer): PASS-WITH-FIXES — restricted the Hopkins ratio + slope
   to the published <1 Jy validity range, corrected the docstring validity bounds, relabelled the
   slope, added Condon 1984 / de Zotti 2010 citations and the "sub-Euclidean ≡ rising EN counts" note.
4. **Paper** `papers/sourcecounts/` — a method/reproducibility contribution: the agreement validates
   the pipeline and NVSS, not new physics.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` round-trips the Hopkins counts (ratio ≈ 1) and recovers a sub-Euclidean slope.
- (Real-data) `python -m jansky_research.sourcecounts --ra 180 --dec 30 --radius 8` reproduces the
  metrics, the Euclidean-count figure with the Hopkins reference, and the macros from VizieR NVSS;
  GATE-2 sign-off before the write-up.
