# ucdv findings — RACS Stokes-V two-epoch census of the Gaia UCD sample (plan 59)

Running log, newest section last. Plan: `plans/59-racs-v-ucd-census.md` (controls frozen 2026-10-15,
PR #344). Every number below is in a committed `results/ucdv_*.json` unless marked otherwise.

## GATE 0 (2026-10-15, scout, PR #343)

Novelty: Perger, Seli & Vida (arXiv:2506.21169, full text) is VLASS Stokes I only; RACS-low2
(arXiv:2606.16182, full text) measures V at Stokes-I catalogue positions and names no UCDs. ADS
keyword queries found no targeted forced-V UCD census. Data: Reylé et al. (2018) via VizieR
`J/A+A/619/L8` — **14,915 rows** confirmed on download (Gaia DR2 ids, positions at J2015.5).

## Pre-data amendments to the frozen controls (both before any real pixel was evaluated)

1. **C1 null rule.** The first frozen version (95% sign interval, one null per target-epoch,
   n ≈ 80) fails a *perfect* pipeline 6.6% of the time (17% at n = 40; 20,000 simulations). The
   synthetic offline run tripped it on pure noise (std 0.83). Now 5 nulls per target-epoch
   (n ≥ 400) and a 99% interval: 0.8% false-fail. Cutouts enlarged to 0.06° so a 2′ null plus its
   90″ annulus, and the 3′ bright-source guard, fit. (CLAUDE.md: the `hiblend` gate lesson.)
2. **Provenance split.** The first version split set A/B by Gaia `source_id` parity. **All
   14,915 Reylé ids are even**, so every target would have been in calibration set A and nothing
   could certify. Now the first bit of SHA-256 of the decimal id: 50.2% B overall.

Both were caught by building the machinery against the plan before opening data — the order the
slice pattern asks for.

## Footprint (CASDA ObsCore metadata only)

RACS V products: 1,177 low2 (λ 0.338 m, MJD 59600–59899) and 1,702 low3 (λ 0.318 m, MJD
60200–60499). RACS-low1 has none. **12,751** UCDs lie in both a low2 and a low3 V footprint
(low2 stops at Dec +48°); 6,374 in set B. Median distance 87 pc; 20 within 20 pc, 1,181 within
50 pc. (A first footprint test passed the antipodal region too — every target "covered" — caught
by the count and fixed by orienting each edge from the polygon centroid.)

## S0 — throughput and scope (results/ucdv_s0_throughput.json)

40 random set-B targets × 2 epochs with a single SODA job per I+V pair (astroquery's default is
a job per product, polling every 20 s; here 3 s). **80/80 pairs ok, 0 failures, 22.0 s per
pair.** The full both-epoch sample projects to 156 h > 120 h, so the frozen rule sets the census
to **d < 107.8 pc: 9,805 targets, 119.99 h**, chosen from parallax and timing only. Option (b),
batched multi-position jobs, was not tried: (a) already reaches 77% of the sample, and the 2,946
dropped are the most distant, whose limits are weakest. Recorded as a deviation from the plan's
"(a) and (b)".

## C1, C2 on real cutouts, set B (results/ucdv_controls_c1c2.json)

First evaluation of any pixel value; target V itself not evaluated.

- **C1 null, V (400 positions): pass.** 51.0% positive (99% interval 43.5–56.5%), std(V/σ) 0.96,
  0 at ≥ 5σ. The same positions in Stokes I: 51.0% positive, std 1.01 — I is genuinely forced too.
- **C2 injection (480 per level, target pixels + null sites): pass.** Completeness 2.9 / 47.7 /
  96.5 / 100 / 100% at 3 / 5 / 7 / 10 / 20 σ; median recovered/injected 0.982 at ≥ 10σ.
  Target pixels alone (80 per level): 95.0% at 7σ, ratio 0.972. The 2–3% low bias is the
  nearest-pixel read of a beam centred off the pixel centre; inside the frozen 0.95–1.05.
- **Image noise:** median V σ 0.149 mJy/beam (0.116–0.269) over the 80 cutouts, so 3σ limits
  near 0.45 mJy — deeper than the 0.25 mJy assumed in early estimates.

Deviation: C2 asked for 200 injections per level "at the target pixel"; 80 target cutouts give 80,
so null sites in the same cutouts were added to reach 480. Both are reported.

## C3 recover-a-known (results/ucdv_controls_c3.json)

- **C3(a): pass.** All 67 `source_type == 'star'` rows of the RACS-low2 V catalogue (CASDA
  `AS110.racs_low2_v_v01`, 221 rows) fetched in their own SBID (stars flare, so the most-central-
  tile rule is overridden by `obs_id`), forced V at the catalogue's fitted position: median
  forced/catalogue **1.021** (IQR 1.00–1.09; frozen range 0.85–1.15), **0 handedness flips in 67**,
  local σ_V / catalogue noise 0.98. Largest outlier: RACS-LOW2 J125545.3−445152, 0.60 → 0.81 mJy
  with σ 0.12 (1.7σ apart).
- **C3(b) WISE J0623−0456: not run yet.** It is a T8 dwarf, too faint for Gaia and so not in the
  Reylé sample; the plan makes it "reported, not pass/fail" and requires the full text of
  arXiv:2306.15219 to be read before any flux from it is quoted. Pending, and not a gate for the
  census.

All three frozen controls pass on real data before the census starts.

## Census launch (2026-10-15)

50-slot trial first (part of the census, kept): 50/50 ok, 21.0 s per slot, the 250 per-slot nulls
52% positive with std(V/σ) 1.10 and none at ≥ 5σ, median σ_V 0.159 mJy/beam, ~11 MB of float32
cutouts per 50 slots (≈ 4.4 GB for the census). Target V not evaluated. The full run then started
as `ucdv-census.service` (systemd user unit, resumable) from a detached worktree pinned to
`38426ad`, 9,805 targets × 2 epochs = 19,610 slots, ~115–120 h expected.

## Census stall and restart (2026-10-16)

The census ran at ~165 slots/h until 07:54 EDT on Oct 16, then every CASDA login failed for
~11.5 h; the runner's fixed 60 s retry made ~285 failed OPAL logins before it was noticed and
stopped at 19:28 (1,694 of 19,610 slots done, all `ok`). A manual login then succeeded — no lockout.
Retries now back off exponentially (60 s doubling, cap 30 min; ~29 logins for an outage that
long) and log timestamped errors. Failed slots were never written, so no data was affected.
