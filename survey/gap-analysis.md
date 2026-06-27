# Gap analysis & domain pick (GATE 1)

Ranking the surveyed gaps for a **single vertical slice**: one gap → one tool → one Airflow DAG →
one public dataset → one paper. Scored on five criteria (H/M/L). The automation layer
(Airflow-on-Podman) is novel for radio astronomy *regardless* of which science domain is chosen, so
the science gap is judged on its own merits.

| # | Candidate gap (tool) | Dataset (size, access) | jansky reuse | Impact | Tract. | Data | Offline | Reuse |
|---|----------------------|------------------------|--------------|:------:|:------:|:----:|:-------:|:-----:|
| **1** | **FRB burst-statistics toolkit** — Weibull wait-times, power-law energy, activity-epoch detection, repeater vs non-repeater stats | **CHIME/FRB Cat 1 CSV** (~a few hundred KB, direct download, no auth) | `jansky.transients`, `jansky.plotting`; numpy/scipy | **H** | **H** | **H** | **H** | **M** |
| 2 | Pulsar-timing quick-look — `par,tim → residuals + Lomb–Scargle red-noise PSD + small-glitch scan` | NANOGrav `.par`/`.tim` (KB each; **NGC 6440E already vendored**) | `jansky.timing`, `jansky.transients`; `jansky[pulsar]` (PINT) | H | H | H | M¹ | H |
| 3 | RFI: reproducible SK-vs-SumThreshold comparison + flag-quality scorer | synthetic + HTRU2 CSV (<1 MB) / small `.fil` | `jansky.rfi`, `jansky.transients` | M | H | H | H | H |
| 4 | HI tangent-point rotation-curve extractor (with uncertainty) | **LAB (l,v) slice already vendored** (366 KB) | `jansky.data.synthetic_hi_cube`, spectral-cube | M | H | H | H | M |
| 5 | CPU-only SETI drift-search + injection-recovery benchmark | Voyager-1 GBT `.h5` (~50 MB) / setigen synthetic | `jansky.seti`, `jansky.formats` | M | M | M | M | H |

¹ PINT needs a one-time solar-system ephemeris download (cached after); synthetic fallback via
`jansky.timing.simulate_pta_residuals`.

## Why each is a real gap (not already solved)

1. **FRB burst-stats** — confirmed: no pip-installable burst-statistics library; Weibull/energy fits
   live in private per-paper scripts; FRBSTATS is web-only. The science question (is a repeater's
   wait-time a single Weibull across epochs?) is *live* and answerable from the public CSV.
2. **Pulsar quick-look** — confirmed: la_forge is downstream of enterprise; enterprise/discovery
   assume HPC/GPU; no `par,tim → residual+PSD` one-liner; no pip small-glitch scanner (IAR 2024).
3. **RFI comparison** — confirmed: no reproducible head-to-head; the field's acknowledged
   reproducibility hole. Strong tooling value but the "finding" is methodological, not a discovery.
4. **HI tangent-point** — confirmed: no turnkey tested extractor (RotCurves targets high-z IFU; galpy
   takes a curve as input). Solid tooling, but the result reproduces a known curve.
5. **SETI** — confirmed: turboSETI heavy, hyperseti/BLISS GPU/C++; setigen has no paired detector.
   Slightly larger data + blimpy dependency.

## Recommendation

**#1 — FRB burst-statistics on CHIME/FRB Catalog 1.** It is the best single slice: the dataset is a
tiny, openly-downloadable CSV (no auth, no GPU, fully offline after one fetch), the gap is a
confirmed absence in the open-source ecosystem, and — unlike the tooling-only options — it yields a
genuine, honest, quantitative **finding** (fitted clustering/energy parameters, repeater vs
non-repeater statistics) with real scientific interest, while staying well within an amateur,
CPU-only budget. It pipelines cleanly through Airflow (fetch CSV → fit distributions → figures/tables
→ paper) and composes existing `jansky.transients` helpers.

**Runner-up: #2 (pulsar quick-look)** — the "safest" data story (NGC 6440E is already vendored) and
highest jansky reuse; choose it if a tooling/usability contribution is preferred over a
distribution-fitting finding.

**GATE 1 is the human's call** — the choice is presented for sign-off before any tooling is built.

## DECISION (GATE 1 — approved)

**Chosen: #1 — FRB burst-statistics on the CHIME/FRB Catalog 1 CSV.**
- **Tool:** `jansky_research.frbstats` — Weibull wait-time clustering (bootstrap CI on *k*),
  Clauset/Hill power-law fluence index, repeater-vs-non-repeater KS tests.
- **Dataset:** `chime-frb-catalog` (registered in `data.py`); synthetic fallback
  `frbstats.synthetic_catalog` for offline tests/CI.
- **jansky reuse:** `jansky.transients` / `jansky.plotting` (+ numpy/scipy); no extra needed.
- **Built:** P3 tool + `pipeline.py` + `report.py` + tests (94% coverage) — committed.

