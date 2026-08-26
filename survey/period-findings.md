# Findings — FRB repeater activity-periodicity search (CHIME/FRB Catalog 1)

Run of `jansky_research.frbperiod.run` over the 18 repeating sources in the First CHIME/FRB Catalog
(CHIME/FRB Collaboration 2021, ApJS 257, 59; one row per event, `sub_num == 0`), using a
phase-folding **Rayleigh $Z^2_1$ periodogram** (Buccheri et al. 1983) on each source's burst arrival
times (MJDs). Full per-source table: `survey/period_results.csv`.

## Result: the known 16.35-day period is recovered; nothing new

- **FRB 20180916B — period $= 16.33$ days.** The periodogram peak coincides with the famous
  16.35-day activity period first reported by CHIME/FRB Collaboration (2020, Nature 582, 351;
  $16.35\pm0.15$ d) — recovered here **from the public catalogue alone**. The recovered value
  (16.326 d, 19 bursts) sits within the grid resolution (0.008 d) and within the published
  $1\sigma$ uncertainty ($0.16\sigma$ away). **That coincidence with the independently-published
  value is the validation** — not the analytic peak strength ($Z^2_1=33.4$; exposure-blind
  $\mathrm{FAP}=7\times10^{-6}$, which ignores survey sampling and is *not* a rigorous significance,
  see below). See `paper/figures/periodogram.pdf`.
- **FRB 20180814A** (11 bursts, the only other searchable source) shows **no significant period**:
  the best peak is $Z^2_1=9.4$ at 2.8 d with $\mathrm{FAP}\approx0.76$ — consistent with noise. No
  multi-day activity period has been reported for FRB 20180814A in the literature, so this null is
  expected. (The $\sim$157-day period in the literature belongs to a *different* source,
  FRB 20121102A — Rajwade et al. 2020 — not to FRB 20180814A.)
- **The other 16 repeaters have $\le 3$ bursts each** in Catalog 1 and are not searchable
  (`min_bursts = 8`). FRB 20121102A appears with a single Catalog-1 event only because most of its
  many hundreds of known bursts pre-date the catalogue's 2018–2019 survey window — it is not a
  sparse repeater in the wider literature.

So: **one detection, and it is the already-known one.** No new periodicity — the public catalogue is
simply too sparse (only 2 of 18 repeaters have enough bursts).

## Honest limitations

- **Exposure-blind significance.** CHIME is a transit survey: each source is seen $\sim$once per
  sidereal day with strongly non-uniform exposure. The analytic FAP treats trial frequencies as
  independent, which the aliased transit window violates, so it is an *approximate* number, not a
  rigorous detection significance. The real confirmation is the match to the published period.
- **Topocentric times.** The Catalog-1 MJDs are topocentric (after DM correction to infinite
  frequency); no barycentric correction was applied. For a 16-day period the resulting Roemer-delay
  modulation ($\sim$0.03% of the period) is negligible and the correct recovery confirms it — but
  this would matter for any search at sub-day periods.
- **Daily aliasing.** A catalogue-only periodogram aliases near 1 day and its beats; 16.33 d sits
  well above that and coincides with the known value, so the recovery is robust, but blind peaks at
  short periods (e.g. the 2.8 d for FRB 20180814A) are not trustworthy.
- **Sparse sampling.** With $\le 3$ bursts, 16 of 18 sources yield no constraint; even the two
  searchable sources are marginal by professional standards (CHIME's own detection used a dedicated
  monitoring campaign plus the exposure model).

## Bottom line

A clean, honest **validation** result: the open, tested, CPU-only tool recovers FRB 20180916B's
16.35-day activity period from public Catalog-1 data, and correctly finds nothing significant for
the rest — no over-claimed new periodicity. The contribution is the reproducible periodicity-search
tool plus this validation and the per-repeater limits; finding *new* repeater periods would need the
denser, exposure-modelled data the professional teams use, not Catalog 1.

## Referee round on the style conversion (2026-08-23)

Verdict *minor revision*. No number, macro or citation changed. One MAJOR, fixed:

**A declarative heading turned an anticipated null into a reported measurement.** `\paragraph{A
null, as expected.}` had become `\paragraph{No significant period in the second source.}` --- and
a `\paragraph` head is the most scannable text in the section. The paper quotes **no sensitivity
at all** for that source: FRB20180814A's 11 bursts, best period 2.792 d, `z2` 9.4 and `fap` 0.757
appear only in `survey/period_results.csv`, never in `main.tex` or `results/period_metrics.json`,
and the body's mitigating clause ("too few catalog bursts to constrain anything") is explicitly
scoped to the *other* sources. Restored the expectation to the heading.

Also fixed: the abstract's prospective scoping ("a deliberately modest question") had gone,
leaving an opening that reads like a paper expecting its answer to be news; `The contribution is X`
(delimiting --- X exhausts the claim) had become `This work is X` (open-ended) in the final
paragraph; "on `\fpNbursts` bursts" could attach to the agreement claim rather than to the
periodogram, so it was moved adjacent to "peak"; the four-item Limitations enumeration was restored,
since the closed set of four labelled axes is what evidences "structural rather than incidental".

**Pre-existing, out of scope:** `results/period_results.csv` holds only synthetic rows
(`SYN-PER`, `SYN-RND`) while the real per-source table lives at `survey/period_results.csv` ---
documented behaviour of `run()` (`csv_dir = op / ("results" if offline else "survey")`), but it
puts offline output in `results/` under a name a reader will take for real evidence, and
`results/period_metrics.json` carries no `is_real` field.

## Full referee round (2026-08-26): MAJOR REVISION, 14 findings, no blocker

Unusually reproducible: the referee re-derived every committed number exactly from
data/chimefrbcat1.csv without touching slice code (19 bursts, P=16.3255, Z²=33.4000,
n_indep=129, FAP=7.2088e-06; the null source's 2.7922/9.40/0.757), all three load-bearing
citations Crossref-verify, the data URL is live, and the recover-a-known framing is honest.
What blocks acceptance: both central quantitative statements are wrong by measured amounts.

**MAJORs (all measured from the committed catalogue):**
- FIVE OF THE 19 "BURSTS" ARE THE SAME CHIME TRANSIT (19 arrival times on 14 distinct sidereal
  days; within-day separations 2–19 min = duplicate phases at P=16.33 d, exactly the
  pseudo-replication pipeline.load_catalog_csv already collapses for sub-bursts). One epoch
  per transit: n=14, Z²=23.87, analytic FAP 8.4e-4 — the recovery is untouched (16.317 d) but
  the abstract's number moves two orders of magnitude.
- THE MEASURED NULL IS 10⁻³–10⁻², NOT 7×10⁻⁶ — and for a reason the paper doesn't name.
  Four Monte-Carlo nulls on the same grid: uniform, sidereal transit comb, and pooled-Cat-1
  draws all roughly CALIBRATE the analytic FAP (aliasing — the paper's stated worry — is
  measurably not the problem); what wrecks it is BURST CLUSTERING (a null bootstrapped from
  another repeater's intervals: P(max Z² ≥ 33.4) = 7e-4; combined with transit collapse,
  7.7e-3). The abstract disclaims the FAP in the same sentence it prints it.
- "CORRECTLY FINDS NOTHING ELSEWHERE" HAS 34% POWER: injecting a 20180916B twin (duty cycle
  0.244, measured) into the null source's 11 bursts, the paper's own significance rule detects
  it 34% of the time — a real analogue would be missed two times in three. Honest form: "no
  period detected, sensitive only to duty cycles ≲0.2". (At n=19 the same injection is
  detected 100% of the time — the asymmetry is purely burst count.)
- THE COINCIDENCE PROBABILITY — the paper's own stated validation — IS NEVER COMPUTED: the
  chance a null peak lands within ±0.15 d of 16.35 is 0.002–0.01 across all four nulls
  (~1/300, ~3σ-equivalent), i.e. the real evidence is ~10³ weaker than the abstract's number
  and still perfectly good. Print THAT.
- FIGURE CLOBBER HOLE: `python -m jansky_research.frbperiod --offline` (default --out .)
  overwrites the real periodogram.pdf with a synthetic one whose peak (16.334) would make the
  real caption still read correctly. JSON and macros are guarded; the figure is not.

**MINOR/NIT:** P quoted with no uncertainty — bootstrap σ(P)=0.031 d, injection σ=0.059 d
(claim survives at <1σ; "within the grid resolution" implies 0.008 d, ~5× too precise); the
metrics JSON records only detections (the searched-null source lives outside every automated
check — the innerrc lesson in mild form); results/period_results.csv is a TRACKED synthetic
CSV in the real-evidence directory; the compiled PDF renders "18 − 2 repeaters" (unevaluated
macro arithmetic — emit \fpNunsearchable); the FAP<0.01 detection threshold and min_bursts=8
are never stated (any threshold 4–11 selects the same two sources — worth saying); "the
correct recovery confirms [the barycentric term is negligible]" → "is consistent with"; grid
bounds hand-typed (currently correct; emit macros); grid adequacy CHECKED and fine (400k
uniform-in-frequency grid finds the same global max); citations all verify; a caution
recorded: the shuffle-intervals null is INVALID here (it inherits the periodicity, P=0.24 —
would wrongly appear to destroy the detection).

**Verdict: MAJOR REVISION.** The single change: replace the analytic FAP in the abstract with
the coincidence probability the paper itself calls the real validation (p ≈ 10⁻³, computable
from committed data in minutes).
