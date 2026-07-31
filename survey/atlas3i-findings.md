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

261 two-position survivors total; **0 survive vetting**. Our 16σ narrowband EIRP limit with
the paper's own L-band parameters (SEFD ≈ 10 Jy, 300 s scans, 2.79 Hz channels, d = 1.80 au)
is **99.4 mW — matching the paper's ~100 mW headline**.

## What the two-stage filter let through, and why (the interesting part)

- **Zero-drift terrestrial carriers** (blc25 pass 1): six tones at 1302–1308 MHz (inside the
  1250–1350 MHz aeronautical-radionavigation/ATC-radar region) passed the on/off filter
  because they sat *below the 16σ search threshold in the OFFs* while above it in the ONs.
  The stamp vet sees them directly in the OFF data (S/N 4–14) and kills them; exactly-zero
  drift is additionally flagged terrestrial by construction.
- **Satellite downlinks** (blc23): nine candidates passed *both* the on/off filter and the
  drift-coherence vet — every one at 1544 MHz (Inmarsat/MSS) or 1619–1625 MHz (Iridium),
  with LEO-Doppler-scale drifts (−3.5 to +1.8 Hz/s) and duplicate frequencies at multiple
  drift rates (an RFI forest, not a point transmitter). Satellites are genuinely in the sky
  and intermittent, so no two-position filter can reject them; a frequency-allocation
  exclusion axis (Iridium/Inmarsat/GNSS bands, `SATELLITE_BANDS_MHZ`) is *necessary*, which
  is why BL pipelines carry the equivalent masks. Good cautionary paragraph for the paper.
- **Chance coincidences** (blc25 pass 1, high drift): with ~10⁵ hits/scan in the noisiest
  node, ±32-channel matching across 6 scans passes ~10 chance candidates; none track their
  own drift rate through the stamps.

## Caveats (honest scope)

- **L band only.** The S (blc22–27, t0 21817), C (23 nodes, t0 26882) and X (25 nodes,
  t0 31308) cadences are pinned but not yet processed — the reproduction claim covers
  939–2064 MHz, not "1–12 GHz". No silent caps: this is 1.125 GHz of the ~9 GHz total.
- **Low-drift blind spot.** Per-channel bandpass excision suppresses tones drifting less than
  ~1 channel per scan (|drift| ≲ 0.01 Hz/s) — such tones are observationally
  indistinguishable from terrestrial carriers in any single-dish search; turboSETI-family
  searches share this blindness (they reject drift = 0 too).
- **S/N conventions differ.** Our 16σ is against per-drift-row MAD noise of the de-Doppler
  plane; the paper's "~16σ" follows the Choza et al. (2024) turboSETI convention. They are
  comparable but not identical detectors; the matching EIRP limit is parameter-level, not a
  claim of identical completeness.
- **Distance is a constant** (1.80 au) pending a Horizons pin at the epoch (TODO before the
  paper; changes the EIRP limit as d²).
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
