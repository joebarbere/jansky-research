# Findings — PTE-II per-source giant-pulse census (plan 47, fable-ideas F10)

`jansky_research.pte2`: a uniform, per-source giant-pulse (heavy-tail) test across all 363 pulsars in
the Parkes Transient Events II database, correlated with spin-down luminosity. Reuses the merged
`frbstats` power-law machinery and `ppdot` ATNF Ė helper.

## GATE 0 (2026-07-10) — CONDITIONAL PASS, novelty NARROWED

- **Paper**: Yang+2025, *Parkes transient events II* (arXiv:2508.14403, **ApJS**
  10.3847/1538-4365/adfe67): 165,592 single pulses / 363 pulsars / Parkes MB 1997–2001; SQLite3 DB,
  raw segments preserved. **Not infrastructure-only** — it already fits log-normal/Gaussian *fluence*
  distributions to the 98 highest-count pulsars + a population power-law (α=−1.7±0.1) + a 25.3%
  giant-pulse-fraction stat.
- **Prior art**: HTRU-V (Burke-Spolaor+2012, MNRAS 423, 1351) is a 315-pulsar log-normal single-pulse
  energy census. So **a generic "363-pulsar energy-distribution census" is anticipated and NOT novel.**
- **The defensible wedge (this slice)**: uniform PER-SOURCE **heavy-tail model test** — a genuine
  power-law giant-pulse tail vs the source's own log-normal bulk — across ALL 363, ranked against Ė.
  HTRU-V tested log-normal/Gaussian (never a heavy-tail model); PTE-II fit only 98 and its power-law
  is on the *population* fluence-ratio distribution, not per-source tails. The per-source
  heavy-tail-vs-Ė census is unclaimed.
- **Data**: OPEN, no auth. GitHub LFS `Astroyx/Pulsar_collection` → `Pulsar_fits_database_v1.zip`
  (~1.5 GB, sha256 775dffa0…; the `raw` URL serves the real zip, GitHub resolves LFS); mirror CSIRO
  DAP 10.25919/34am-zx04. Real leg runnable.
- **Schema**: `pulsar`(jname, p0, s1400, spNumber) ← `file`(pfLinkID, pulsarID, timeStartMJD) ←
  `fileSegment`(pfLinkID, snr_max/snr_min, data BLOB). Per-pulse **fluence is NOT stored** (derived
  from blobs by the paper's tool); we use per-pulse **S/N** (`snr_max`) as the energy proxy. Join to
  ATNF Ė by `jname`=PSRJ.

## Method (and the honest observable)

- **Observable = per-pulse S/N** (`snr_max`). Within one source (fixed Tsys/gain) S/N ∝ fluence, so
  it faithfully probes the *shape* of that source's energy distribution — exactly what a tail test
  needs — but forbids *absolute* cross-source energy comparison. We compare only tail shape and the
  heavy/not-heavy classification across sources, never absolute energies.
- **The giant-pulse excess test** (`fit_energy_tail`): estimate the log-S/N median `m` and right-side
  width `σ_R = P84.13(log S/N) − m` **from the upper half only** — which the single-pulse detection
  S/N floor does not truncate — then count pulses above `exp(m + 3σ_R)` and compare to the log-normal
  expectation with a **Poisson excess test**. A pulsar is `heavy_tailed` only on a significant
  (p<0.05), ≥3-pulse excess: a positive giant-pulse detection.
  - **Why not a bare power-law-vs-lognormal likelihood ratio**: on a short upper tail the two are
    famously degenerate, and a naive Vuong test spuriously prefers power-law even for pure log-normal
    data (measured FP ~0.35 in early tests). The floor-robust excess test brought the FP to ~0.05.
    The Vuong statistic + giant-tail power-law index (Clauset MLE) are reported as *secondary*
    diagnostics only.
- **Attribution** (`tail_vs_edot`): Spearman rank-correlation of the giant-pulse excess with
  `log Ė`, plus a Mann–Whitney comparison of `log Ė` for heavy vs non-heavy sources. Giant pulses
  are physically tied to high Ė (Crab, MSPs), so a positive trend is the expectation.

## Synthetic recover-a-known (offline, CI)

Injected log-normal bulks (median S/N ≫ detection floor) ± a power-law giant tail:
- **Completeness rises with pulse count**: ~0.33 at 80 pulses → 0.88 at 150 → 0.98 at 300 → **1.0 at
  ≥600** — the honest count-dependent floor (a giant tail in a sparsely-sampled pulsar is
  unrecoverable).
- **False-positive rate ~0.05** on pure log-normal sets (at the test's 0.05 significance).
- A clear injected heavy tail is classified `power_law`; a pure log-normal is not (`recovered`).

## The real census (ran 2026-07-10): an HONEST NULL

Downloaded the 1.5 GB PTE-II SQLite DB, loaded per-pulsar S/N for all 363 pulsars, ran the
giant-pulse test on the **136** with ≥50 single pulses, cross-matched ATNF Ė (135 matched).

- **26 of 136 (19%) flagged with a significant log-normal excess** — we call these
  *log-normal-excess* sources, NOT giant-pulse emitters (the three points below show they largely
  aren't). (Not meaningfully comparable to Yang+2025's 25% giant-pulse-fraction, which is a
  different statistic.)
- **No Ė trend** (the wedge result): Spearman(excess, log Ė) = **−0.03, p=0.76**; heavy median
  log Ė 33.29 vs non-heavy 32.84, Mann-Whitney **p=0.13**. Giant pulses are physically tied to high
  Ė; the archival S/N census shows no such trend.
- **Detection-power floor** (`count_confound`, `count_limited=True`): heavy sources have **2.4×
  more pulses** (median 552 vs 228, MW **p=0.02**); flagged fraction rises **0.088 (low-count half)
  → 0.294 (high-count half)**. A significance test has more power on more pulses, so a small
  departure from log-normal reaches p<0.05 only where the count is large → the 19% is a power floor,
  not a clean incidence.
- **Steep tails, absent archetypes**: only ~5 flagged sources have ≥8 giants to fit an index; for
  those the median is **11.6** (and a Hill index on <30 points above a high threshold is very
  uncertain) — the excess does NOT continue into a shallow (~2–3) power law, i.e. evidence *against*
  a giant-pulse tail. The classic GP emitters (Crab J0534+2200, B1937+21 J1939+2134, Vela
  J0835−4510) do **not** survive the 50-pulse cut in this 1997–2001 Parkes sample, so the census has
  essentially no real-data positive control; of known tail-heavy sources only **B0950+08
  (J0953+0755)** is present, and it is correctly flagged (excess 0.86, p=1.5e-6). An instrumental
  contribution (RFI, calibration jumps, off-axis bright sources) to the excess cannot be excluded —
  which only reinforces the null.
- Top flagged sources: J1243−6423 (n=2731, 141 giants, excess 1.53), J1224−6407, J1633−5015,
  J1535−4114, J1146−6030 — the high-excess ones are also the highest-count ones (the confound).

**Honest bottom line**: a uniform per-source heavy-tail test over PTE-II flags mild high-S/N
excesses in ~19% of well-sampled pulsars, but the flagged fraction is **detection-power limited**,
the tails are **too steep to be classic giant pulses**, and the incidence **does not correlate with
Ė**. No giant-pulse population and no Ė trend are claimed. The negative result is not a tooling
failure — the synthetic recover-a-known confirms the test recovers injected heavy tails — but a real
limit of an S/N-only, count-heterogeneous archival census. The confound is pipeline-generated
(`count_limited`, count-split fractions, median tail index), so the artifact carries the conclusion.

## GATE-2 (2026-07-10) — PASS, no blockers

The reviewer verified the math (Poisson excess one-sided; floor-robust `σ_R` genuinely immune to the
low-S/N truncation; Vuong algebra correct but demoted to a secondary diagnostic), confirmed all
citations, and confirmed the honest framing holds (null stated plainly; the 19% de-emphasized as a
detection-power floor; `count_limited` load-bearing and correctly interpreted; the Ė-null uses the
count-robust continuous `excess`, not the biased binary flag). Required should-fixes, all applied:

- **Terminology**: renamed the flagged 19% "log-normal-excess sources" (not "giant-pulse emitters")
  in the prose, reserving giant-pulse language for shallow tails.
- **γ scope**: clarified the median index ≈11.6 is over only the ~5 sources with ≥8 giants, that a
  Hill index on <30 points is very uncertain, and that steep γ is evidence *against* a real power-law.
- **RFI/instrumental caveat** added (archival 1997–2001 single-pulse excess could be RFI/calibration/
  off-axis, not pulsar emission) — reinforces the null.
- Dropped the misleading "broadly consistent with Yang+2025's 25%" comparison (different statistic).
- Removed the dead `KNOWN_HEAVY_TAIL` constant (the archetype check was never runnable on real data).
- Reviewer's read on novelty: honest and appropriately thin — a validated per-source test + a
  pipeline-generated detection-power confound + a clean Ė-null = a citable negative-results/methods
  note; does not overclaim.

## Reproduce

Offline (test + synthetic recover-a-known + synthetic-schema parser test):
`uv run python -m jansky_research.pte2 --offline --out .`
Real (downloads the ~1.5 GB PTE-II SQLite DB): `uv run python scripts/pte2_real.py --cache /tmp/pte2`
(writes `results/pte2_metrics.json` + `results/pte2_census.json`, `is_real=True`).
