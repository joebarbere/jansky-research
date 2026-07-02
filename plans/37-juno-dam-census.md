# 37 — Jovian DAM occurrence census from Juno/Waves: the Io-controlled emission zones from orbit

Status: ✅ done — tooling + real 2017-02 month + paper merged PR #87; GATE-2 caught the Io-phase convention BLOCKER (Φ_Io = CML+180−Λ_Io), re-measured contrast 1.38→2.22; multi-orbit v02 census is the scoped follow-on

## Context

The Juno/Waves Estimated Flux Density Dataset (doi:10.25935/6jg4-mk86, public Oct 2025) provides
calibrated 1-s dynamic spectra (110 channels, 1 kHz–40.5 MHz — the full decametric window) for
2016–2019 (v01; v02 extends to 2023), with per-channel `Background` and `Sigma` arrays. The classic
ground-based result (Bigg 1964; Carr+1983; Marques+2017's 26-yr Nançay catalogue) is that Jovian
DAM occurrence is organised in the (CML, Io-phase) plane into the Io-A/B/C/D source regions. An
occurrence census of that plane **from Juno's vantage** — moving, close, polar — using the public
calibrated data is amateur-tractable and not published from this dataset. Recover-a-known: the
census must re-find Io-controlled enhancement; the new-findings edge is the orbital vantage.

GATE 0 (verified live): daily CDFs at
`maser.obspm.fr/repository/juno/waves/data/l3a/data/cdf/YYYY/MM/jno_wav_cdr_lesia_YYYYMMDD_v01.cdf`
(37 MB/day, no auth, CC-BY; vars Epoch/Frequency/Data/Background/Sigma, units V²m⁻²Hz⁻¹).
CML: sub-Juno System III via JPL Horizons (`PDObsLon`, id 599, location `500@-61`) — the naive IAU
W_III formula is WRONG for Juno by up to ~40° (verified). Io phase: Lieske (1987) mean longitude
l₁ = 106.07719 + 203.488955790·(JD−2451545) minus CML. Bounded first leg: one month of CDFs (~1.1 GB).

## Deliverables

- `src/jansky_research/junodam.py`: CDF reader (reuses the `windwaves` cdflib extra) with DAM-band
  (3–40.5 MHz) selection + time binning; `detect_active` (Data > Background + kσ per bin);
  `io_mean_longitude` (Lieske, cited); `fetch_cml` (Horizons batch, pragma); occurrence map in the
  (CML, Io-phase) plane + Io-region contrast statistic; `synthetic_orbit` fixture (linear
  CML/Io-phase rates + injected Io-B-like box + non-Io bursts) → recover-a-known; run/figure/macros.
- Tests (85% floor, offline); `papers/junodam/`; `survey/junodam-findings.md`; wiring; GATE-2.

## Approach

0. GATE 0 done. 1. Tooling + synthetic round-trip (recover the injected Io box). 2. Real leg:
one month (background download), Horizons CML at 15-min sampling interpolated, occurrence map,
Io-region contrast vs the canonical zones. 3. GATE-2 (units caveat: V²m⁻²Hz⁻¹ not W-flux —
occurrence needs only detection above background; Juno-frame CML caveat; single-month coverage
of the beat plane). 4. Paper: census method + Io-control recovery from orbit.

## Verification

Offline round-trip recovers the injected box contrast; real leg shows Io-region enhancement;
checks green; SHA-verified CI; GATE-2 sign-off.

## Risks & mitigations

- Horizons batch limits → 15-min sampling (~3k epochs/month), cached to CSV, interpolated.
- (CML, phase) plane coverage in one month → beat period ~13 h fills the plane in weeks; report
  per-cell exposure and mask empty cells.
- v01 vs v02 calibration → v01 (the DOI dataset) for the census; v02 noted as follow-on.
- Proximity effects (Juno inside the emission region at perijove) → flag perijove days; the
  census is a vantage-dependent occurrence map, stated as such, not a Nançay replacement.
