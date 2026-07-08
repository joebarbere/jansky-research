# Findings — radio survey of the 56 WD-pulsar candidates (plan 41)

`jansky_research.wdpulsar` + `scripts/wdpulsar_real.py`: forced RACS I+V photometry (every
complete I+V observation per target) + VLASS cone checks of the Pelisoli+2025 candidate list.

## GATE 0 (2026-07-07)

- **Gap open**: only two papers cite Pelisoli+2025 (MNRAS 540, 821; arXiv:2505.04693) — both
  single-object. No systematic radio survey exists. The adjacent Pelisoli+2024 VLASS survey
  (MNRAS 531, 1805) covers single WDs, not this binary candidate list.
- **No VizieR/CDS deposit** ("data available upon reasonable request"); the machine-readable
  source is Table 2 of the arXiv HTML — vendored as `data/wdpulsar_candidates.csv` (56 rows,
  positions to ~0.1″, provenance-flagged; the `lpt` pattern).
- **Transcription lesson (validation caught it)**: Table 2's asterisk means "classification
  determined as part of this work" (17 rows) — NOT the abstract's "26 previously
  uncharacterised" (= rows with no Simbad ID). Both flags kept, distinctly named.
- **Plan-anchor corrections**:
  - J1912−4410 (the confirmed WD pulsar, in the list) is NOT a plausible RACS re-detection:
    MeerKAT pulses (<4 s, ~15 mJy, ~1% duty cycle) time-average to ~0.1–0.2 mJy, and it is
    absent from the RACS-low DR1 catalogue (verified cone search). Its forced limit is the
    paper's duty-cycle caveat made concrete.
  - **AR Sco is the control instead** (the list's template, deliberately not among the 56):
    RACS-low DR1 8.58±0.92 mJy; |V|/I ≈ 22–27% at 1.5 GHz (Stanway+2018).
- **Timing**: RACS-low2 (arXiv:2606.16182) just released a circular-polarization catalogue —
  our dataset AND the likely scoop vector.

## Offline validation (done, in CI)

- Vendored-table structural validation: 56 rows, type census (26 YSO/16 polar/3 CV/2 IP/
  8 unclear/1 pulsar), J1912−4410 present, both provenance flags (26 uncharacterised, 17
  classified-this-work).
- Injection round-trip: synthetic RACS-like cutout (15″ beam, 0.25 mJy rms), injected
  I=8.6/V=−2.1 mJy recovered within noise by `stokesv.measure_circular_pol`, classified
  `circular`; blank field stays below 5σ (no fake detection).
- `summarize_sweep` unit-tested: detections, deepest limits, leakage vetting (a formally
  significant V below 0.6% of I is rejected — the characteristic ASKAP on-axis leakage level,
  a weaker raw-level veto than `stokesv`'s per-region floor), uncovered targets.

## VLASS leg (done, 2026-07-07 — local bulk QL catalogues, epoch 3)

42 of 57 targets have Dec > −40. **Two matches within 7.5″:**

| target | type | peak (mJy) | sep (″) | reading |
|---|---|---|---|---|
| AR_Sco (control) | control | 5.40 | 1.4 | control re-found at 3 GHz ✓ |
| J0408+6046 | YSO (this work; previously uncharacterised) | 1.25 | 3.75 | radio counterpart **consistent with the YSO classification** (gyrosynchrotron), not an AR Sco-like signature; chance-coincidence ~2–3% across the 42 checked |

## RACS sweep (DONE — 2026-07-07, after a CASDA outage cleared)

CASDA SODA staging was erroring for ~90 min around midday (external, not our code — the
verified `stokesv` path failed too); a probe-and-sweep loop waited it out and launched at
00:56. Two driver bugs surfaced and were fixed mid-run: (1) SODA occasionally returns I and V
cutouts differing by one pixel → crop both to the common overlap; (2) recent (2025–26) RACS-mid
observations are still proprietary → 12 obs return "no access", dropped as embargoed epochs (each
affected target has 6–21 OTHER released epochs, so none is lost). One process-management
embarrassment: a `pgrep -f wdpulsar_real.py` watch-loop matched its own launcher shell (whose
args contained the string), so the chained retry pass silently never started for ~6 h; relaunched
directly once caught. Final: **636 rows, 582 good measurements across the 51 RACS-covered
candidates + AR Sco**, median I rms 0.25 mJy.

## Result — a clean null (the deliverable)

| quantity | value |
|---|---|
| candidates with RACS coverage | 51 of 56 (5 northern: J0408+6046, J0719+6557, J1125+5012, J1907+6908, J2120+6848) |
| **candidate Stokes-I detections (>5σ)** | **0** (brightest 3.3–3.9σ noise peaks, ≲1.2 mJy) |
| **candidate Stokes-V detections (leakage-vetted)** | **0** |
| median 3σ V limit | **0.41 mJy** (range 0.33–0.53) |
| AR Sco control | **detected, I=4.20 mJy (8.6σ), classified circular, negative V** ✓ |
| J1912−4410 (confirmed WD pulsar) | undetected, 9 epochs, 3σ I<0.52 / V<0.40 mJy |
| VLASS candidate match | J0408+6046 (YSO, 1.25 mJy, 3.75″) — gyrosynchrotron, not AR Sco-like |

**Framing** (honest, as planned): the survey finds no persistent radio counterparts. The
J1912−4410 non-detection calibrates the null — the archetype itself is invisible to snapshot
surveys (~1% duty cycle → ~0.1–0.2 mJy time-averaged), so the census bounds only the
*persistently* radio-loud fraction, consistent with zero at ~0.4 mJy. The AR Sco re-find (I +
VLASS) is the pipeline's recover-a-known; that it re-finds the archetype but none of the
candidates is the internal control.

## Reproduce

Offline (validation + table): `uv run python -m jansky_research.wdpulsar --offline --out .`
Real (needs the sweep CSV + local VLASS catalogues): `uv run python -m jansky_research.wdpulsar
--out .`. The sweep itself: `uv run --extra vlass python scripts/wdpulsar_real.py` (resumable,
CASDA auth via `~/.casda_pw`; ~600 SODA cutouts).
