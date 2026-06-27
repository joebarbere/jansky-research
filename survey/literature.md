# Survey — scientific state of the art (2023–2026)

Synthesised from five domain briefings (radio-research-assistant). Citations are URL-verified;
see each domain's gap rows in `gap-analysis.md`.

## FRBs & radio transients

Repeaters vs non-repeaters is still unsettled — a Zipf-like burst-rate continuum (non-repeaters =
repeaters caught once) vs genuine two populations (CHIME/FRB Cat 1, ApJS 257, 59; arXiv:2506.09138).
Hyperactive repeaters (FRB 20240114A reached ~500 bursts/hr on FAST) are modelled with **Weibull
wait-time** distributions (shape *k* < 1 ⇒ clustering); repeaters show "memory" over minutes–hours
(ApJL 950). Sub-µs morphology, downward drift ("sad trombone"), host-galaxy diversity (a quiescent
elliptical host, ApJL 979 L22), and the Macquart relation (Nature 581, 391) round out active work.
**Amateur-accessible question:** is the wait-time distribution of a repeater a single Weibull across
activity epochs, or do *k* and rate co-evolve? — answerable from the public CHIME catalog.

## Pulsar timing & PTAs

The 2023 nanohertz-GWB evidence (NANOGrav 15yr ApJ 951 L8; PPTA L6; EPTA DR2) via the
**Hellings–Downs** correlation. Open issues (2024–26): an unexplained ~4 nHz **monopole** that, when
included, drops HD evidence by >10× (arXiv:2411.13472); per-pulsar **red-noise** modelling is the
dominant systematic; SMBHB vs cosmological-source discrimination awaits IPTA DR3. Small glitches are
routinely missed by eye and need algorithmic search (IAR 2024). **Amateur-accessible:** load any of
the 68 NANOGrav `.par`/`.tim` (KB-scale each), fit residuals, recover HD, read the red-noise slope.

## HI 21 cm, Galactic structure & continuum source counts

The **tangent-point** method still anchors the inner Milky Way rotation curve; active work targets
bar-induced non-circular bias (~20%, arXiv:2501.12760) and outer-halo profiles (cored Einasto over
NFW, MNRAS 528, 693). Continuum **log N–log S** counts are pinned over 8 decades (MeerKAT DEEP2,
arXiv:2101.07827); RACS-mid (3.1M sources) and MIGHTEE refine faint-end completeness and the AGN/SFG
radio luminosity function (A&A 685, A79). **Amateur-accessible:** tangent-point rotation curve from a
tiny LAB (l,v) slice; log N–log S slope + spectral indices from NVSS/FIRST cone searches (no auth).

## RFI mitigation & ML

Classical **SumThreshold/AOFlagger** (A&A 670, A166) and **spectral kurtosis** (Nita & Gary 2010)
still dominate production and beat ML on CPU throughput. ML (FETCH; U-Net comparisons MNRAS 530, 613)
is capable but GPU/torch-heavy and poorly benchmarked; **reproducibility is the field's weakest
point** — no canonical <100 MB labelled benchmark, papers compare against unspecified SumThreshold
implementations. **Amateur-accessible:** classical SK-vs-SumThreshold head-to-head on the same data;
CPU-only candidate features (HTRU2 CSV, <1 MB) for an sklearn classifier.

## SETI / technosignatures

The Doppler-drift narrowband matched filter (turboSETI Taylor-tree) + ON/OFF cadence veto; **BLC1**
(Nat. Astron. 6, 352) showed intermodulation products can survive a naive cadence test. Recent deep
surveys (UCLA 10-yr >70k stars, arXiv:2605.05408; TRAPPIST-1; 3I/ATLAS) are all null with
ever-tighter EIRP limits via the **CWTFM** figure of merit. **Amateur-accessible:** the Voyager-1
GBT file (~50 MB) has a *verified real* drifting tone for end-to-end pipeline validation; setigen
generates reproducible synthetic cadences.
