# atlas3i — independent reproduction of the BL 3I/ATLAS GBT L-band nondetection (plan 85)

**Date:** 2026-07-30 (search + vet of all six L-band nodes ran this day, node-serial).
**Target paper:** Jacobson-Bell et al. 2025, RNAAS, arXiv:2512.19763 (GBT, 2025-12-18,
1–12 GHz, turboSETI ±4 Hz/s, ~16σ, nondetection to ~100 mW EIRP).
**Data:** public archive `https://bldata.berkeley.edu/ATLAS/GB_ATLAS/`, fine-frequency
(~2.79 Hz) rawspec products, L-band cadence t0 = 17226 s UT, nodes blc21–26
(939.0–2064.0 MHz contiguous, 187.5 MHz per node). No independent reanalysis of these files
existed before this one (checked 2026-07-30; the ATA arXiv:2512.18142 and FAST
arXiv:2603.19023 searches are original observations).

## Result

**The L-band null reproduces.** Independent pipeline (physical-units brute-force de-Doppler,
±4 Hz/s, 16σ against per-drift MAD noise; ABACAD on/off filter; then a per-candidate
drift-coherence vet on postage stamps + a satellite-allocation exclusion):

| node | band (MHz) | hits/scan (min–max) | on/off survivors | confirmed |
|---|---|---|---|---|
| blc21 | 1876.5–2064.0 | 30–779 | 0 | 0 |
| blc22 | 1689.0–1876.5 | 436–557 | 3 | 0 |
| blc23 | 1501.5–1689.0 | 21 267–30 072 | 245 | 0 |
| blc24 | 1314.0–1501.5 | 487–14 698 | 1 | 0 |
| blc25 | 1126.5–1314.0 | 96 516–234 066 | 12 | 0 |
| blc26 | 939.0–1126.5 | 459–36 187 | 0 | 0 |

261 two-position survivors total (1 124 011 raw hits — see the caveat on detector strictness
below); **0 survive vetting**. Our 16σ narrowband EIRP limit with the paper's own L-band
parameters (SEFD ≈ 10 Jy, 300 s scans, 2.79 Hz channels, d = 1.798 au — JPL Horizons,
DES=C/2025 N1, geocentric, 2025-12-18 05:00 UT) is **99.2 mW — matching the paper's ~100 mW
headline**. GATE-2 verified this match is an *algebraic identity*, not a
coincidence: the paper's Eq. 1+2 (Gajjar et al. 2021 form) reduce exactly to our
`eirp_limit_w` once the transmit-bandwidth term cancels — under the assumption β = 1 (perfect
dedrifting efficiency; the paper carries an explicit β for high-drift smearing loss that we
drop, which makes our limit marginally optimistic at the highest drift rates).

## What the two-stage filter let through, and why (the interesting part)

- **Zero-drift terrestrial carriers** (blc25): six tones at 1302–1308 MHz (inside the
  1250–1350 MHz aeronautical-radionavigation/ATC-radar region) passed the on/off filter
  because they sat *below the 16σ search threshold in the OFFs* while above it in the ONs.
  The stamp vet sees them directly in the OFF data (per-scan S/N 3–15; committed in
  `results/atlas3i_blc25_L.json`) and kills them; exactly-zero drift is additionally flagged
  terrestrial by construction.
- **Satellite downlinks** (blc23 and blc25): blc23's nine drift-coherent candidates all sit
  at 1544 MHz (Inmarsat/MSS) or 1619–1625 MHz (Iridium), with LEO-Doppler-scale drifts
  (−3.5 to +1.8 Hz/s) and duplicate frequencies at multiple drift rates (an RFI forest, not
  a point transmitter). Likewise **5 of blc25's 6 high-drift survivors are in the GPS
  L5/Galileo E5 allocation** (1200–1204 MHz, drifts 2.6–3.9 Hz/s) — GNSS-band RFI that fails
  drift-coherence, not statistical flukes; only the 1149.99 MHz one is unclassified (and
  fails tracking too). Satellites are genuinely in the sky and intermittent, so no
  two-position filter can reject them; a frequency-allocation exclusion axis
  (`SATELLITE_BANDS_MHZ`) is *necessary*, which is why BL pipelines carry the equivalent
  masks. (Band edges: Inmarsat and Iridium ranges are ITU-exact; the internal split points
  between the GNSS entries are heuristic groupings of the 1164–1215 and 1559–1610 MHz RNSS
  allocations.)

## Caveats (honest scope)

- **L band only, and wider than the paper's L band.** The S (blc22–27, t0 21817), C (23
  nodes, t0 26882) and X (25 nodes, t0 31308) cadences are pinned but not yet processed —
  the reproduction claim covers our searched span, not "1–12 GHz". Moreover we searched the
  full *recorded* node span (939.0–2064.0 MHz) while the paper analyses its L receiver as
  **1.1–1.9 GHz**: blc21 and blc26 lie mostly outside the paper's passband, in receiver
  roll-off the original team presumably trimmed deliberately. Their hit counts (up to 36 187
  per scan) are plausibly part edge-artifact, and the nominal SEFD = 10 Jy is optimistic
  there; the like-for-like reproduction statement is for 1.1–1.9 GHz (blc22–25 + edges of
  21/26), where the result is the same: 0 confirmed.
- **Symmetric search thresholds (deviation from the paper).** The paper uses 16σ for ON
  scans but a *lower* 10σ for OFF scans, so a candidate slightly weaker in the OFFs still
  registers there and is vetoed. This run used a symmetric 16σ everywhere — a stricter OFF
  threshold, which passes *more* borderline RFI to the two-position filter (part of why 261
  survivors needed a vet stage). The stamp vet subsumes the missing veto — it reads the OFF
  S/N directly at the predicted position with an 8σ bar, stricter than the paper's 10σ — so
  the 0-confirmed outcome is unaffected, but the intermediate survivor counts are not
  comparable to the paper's candidate counts. The asymmetry (`threshold_off=10`) is now
  implemented for future runs.
- **Our first-stage detector fires far more often than turboSETI's.** 1 124 011 raw hits in
  L band alone vs the paper's 471 198 across the whole 1–12 GHz survey (>2× from ~1/8 of the
  bandwidth). The per-drift-row MAD noise estimate is robust to (i.e. deflated by) dense RFI
  forests, so our "16σ" is effectively looser than turboSETI's normalisation in RFI-heavy
  sub-bands. For a null reproduction this errs conservative — we admit more false positives
  into vetting and still confirm none — but raw hit counts must not be compared
  detector-to-detector.
- **GBT L/S hardware notch filters not mapped.** The paper's Fig. 1 marks unsampled notch
  regions (L and S band); we have not verified how those gaps appear in the archived `.h5`
  data or whether gap edges seed spurious hits. Unverified, flagged for the paper pass.
- **Low-drift blind spot.** Per-channel bandpass excision suppresses tones drifting less than
  ~1 channel per scan (|drift| ≲ 0.01 Hz/s) — such tones are observationally
  indistinguishable from terrestrial carriers in any single-dish search; turboSETI-family
  searches share this blindness (they reject drift = 0 too).
- **S/N conventions differ.** Our 16σ is against per-drift-row MAD noise of the de-Doppler
  plane; the paper's "~16σ" follows the Choza et al. (2024) turboSETI convention. They are
  comparable but not identical detectors; the matching EIRP limit is parameter-level, not a
  claim of identical completeness.
- **Distance is Horizons-pinned**: 1.79801 au geocentric at 2025-12-18 05:00 UT (mid
  L-cadence; DES=C/2025 N1, range-rate −1.52 km/s), pinned 2026-07-30. The earlier 1.80 au
  nominal was accurate to 0.1%.
- **SEFD is the nominal 10 Jy** GBT L-band figure, not measured from these scans.

## Reproduce

```
uv run --extra voyager python -m jansky_research.atlas3i --sweep --threshold 16
```

Node-serial: downloads each 6×10 GB cadence, searches in bounded frequency chunks, vets
survivors from the local files, deletes the scans (peak disk ≈ 60 GB; ~2¼ h per node, ~13 h
for the band). Per-node JSONs in `results/atlas3i_blc*_L.json` are checkpoints — rerunning
skips completed nodes. Offline round-trip (synthetic cadence, no network): `run()` /
`python -m jansky_research.atlas3i`.
