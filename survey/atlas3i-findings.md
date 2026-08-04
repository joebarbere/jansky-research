# atlas3i — independent reproduction of the BL 3I/ATLAS GBT L-band nondetection (plan 85)

**Date:** 2026-07-30 (search + vet of all six L-band nodes ran this day, node-serial).
**Target paper:** Jacobson-Bell et al. 2025, RNAAS, arXiv:2512.19763 (GBT, 2025-12-18,
1–12 GHz, turboSETI ±4 Hz/s, ~16σ, nondetection to ~100 mW EIRP).
**Data:** public archive `https://bldata.berkeley.edu/ATLAS/GB_ATLAS/`, fine-frequency
(~2.79 Hz) rawspec products, L-band cadence t0 = 17226 s UT, nodes blc21–26
(939.0–2064.0 MHz contiguous, 187.5 MHz per node). No independent reanalysis of these files
existed before this one (checked 2026-07-30; the ATA arXiv:2512.18142 and FAST
arXiv:2603.19023 searches are original observations).

## Result — full 1–12 GHz coverage (completed 2026-08-04)

**The null reproduces across the paper's entire frequency range.** All four receiver cadences,
every unique recorded node (60 node-cadences; the 12 pairwise-duplicated C/X bank-boundary
recordings are skipped as redundant, documented in `DUPLICATE_NODES`):

| band | nodes | coverage (MHz) | raw hits | on/off survivors | confirmed |
|---|---|---|---|---|---|
| L | 6 | 939.0–2064.0 | 1 124 011 | 261 | 0 |
| S | 6 | 1651.5–2776.5 | 354 423 | 29 | 0 |
| C | 23 | 3939.0–8251.5 | 15 735 | **0** | 0 |
| X | 25 | 7689.0–12376.5 | 444 691 | 4 | 0 |
| **total** | **60** | — | **1 938 860** | **294** | **0** |

The hit-density contrast is itself the picture: from ~2×10⁵ hits/scan in the worst L-band
node down to individual C-band scans recording zero raw hits (C max 3,773, mean ≈114/scan —
about 1.8–3.3 decades below L depending on the comparison); separately, C produced **zero
two-position survivors across 4.3 GHz**. X band's only four survivors are zero-drift tones at 12.0–12.35 GHz —
inside the Ku FSS/DBS satellite-TV downlink allocation and *above* the paper's analysed
7.6–11.7 GHz passband (the roll-off region again); all four fail drift-coherence vetting and
are now also caught by the extended satellite-allocation table. S band's 29 survivors (all
killed) clustered in the 1.9–2.2 GHz PCS/MSS region.

**Computational cost (measured from the run logs):** ~95 h active compute — L ≈ 2.2 h/node,
S/C/X ≈ 1.5 h/node — and ~3.7 TB transferred (366 fine-resolution scans × ~10.2 GB), peak
disk ~60 GB (fetch-search-vet-delete), on one desktop (Ryzen 5 5600X, 64 GB RAM, pure-NumPy
CPU pipeline). Runs were paused for bandwidth management; the figure quoted is active
runtime, not calendar time.

## Result — L band detail (first leg, 2026-07-30)

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

- **Wider than the analysed passbands — all four bands.** Recorded vs analysed spans:
  L 939.0–2064.0 vs 1.1–1.9 GHz; S 1651.5–2776.5 vs 1.8–2.7 GHz (~150 MHz low / ~77 MHz
  high excess); C 3939.0–8251.5 vs 4.0–7.8 GHz (~61 MHz low / ~452 MHz high); X
  7689.0–12376.5 vs 7.6–11.7 GHz (~677 MHz high — where all four X survivors sit, in the
  Ku DBS allocation). Edge-region survivors failed vetting everywhere; S's survivors are
  mid-band (1.9–2.2 GHz) and C's edges produced none. On the L band specifically: we searched the
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
