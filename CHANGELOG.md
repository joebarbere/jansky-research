# Changelog

All notable changes to `jansky-research` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project follows
[Semantic Versioning](https://semver.org/) as codified in [`VERSIONING.md`](VERSIONING.md).

Every PR adds an entry to `## [Unreleased]`. `scripts/next_version.py` reads that section to
recommend the next version number.


## [Unreleased]

### Fixed
- **`triangulate` cited the wrong one of two same-year Krupar companions.** The entry was
  internally consistent --- DOI `10.1007/s11207-014-0522-x` matched its own volume 289, issue 8,
  pages 3121--3135 --- so a coordinate check passes it. But it is the *flux-density* paper, and the
  citation supports an apparent-source-size claim. The ~60 deg figure is in the goniopolarimetric
  companion (`10.1007/s11207-014-0601-z`), whose 125 kHz--2 MHz band is also exactly this paper's
  `\triFlo`--`\triFhi`. When two same-year companions share a title stem, the discriminator is the
  subtitle and the claim, not the identifier. The same sentence also stated ~60 deg as a constant
  where the source reports it at the *lowest* analyzed frequencies; now qualified.
- `triangulate`'s Method said channels are kept when the ray miss distance is "below a threshold"
  and never gave it. It is 60 R_sun, against inferred distances of 15.3--106.1 R_sun, so at the
  high-frequency end it admits rays missing by four times the distance being measured. Stated, with
  its permissiveness noted. Sweeping it needs a real-leg re-run: the committed results file keeps
  only summary scalars, not the per-channel arrays the headline is computed from.

### Fixed
- **`ecallisto_pipeline`'s only quantitative result was drawn from the wrong run.** Found by the
  referee round on the style conversion. `_write_macros` emitted seven un-namespaced names from both
  legs and `run()` wrote both legs to one results file, so the real 2011-09-14 archive day (which
  ran last) owned the macros while **all twelve macro uses in the paper describe the synthetic
  day**. The abstract typeset as *"a real burst at **0** stations ... confirms exactly **0** event
  and rejects 8"*, and the Method attributed the real day's `8 of 1512` and `-5.941 MHz/s` to a
  synthetic run over ten stations. `preserve_live_macros` was powerless: both legs wrote real values
  under one name. Macros are now `\ecSyn*`/`\ecReal*`, each leg writes its own results file
  (`ecallisto_synthetic_metrics.json`, allowlisted like `vgpra_synthetic_*`), and the paper cites
  the synthetic namespace. Correct values: 4 stations, 1 event, 3 rejected, 7 of 10 flagged.
  The committed real day is now cited once instead of sitting unread in `results/`.
- Two lessons worth keeping from that fix. **`preserve_live_macros` accumulates values, not names**
  --- it rewrites the lines the new run emits and drops any existing macro the new text never
  mentions, so a writer emitting only its own namespace *deletes* the other leg's numbers; each leg
  must emit both namespaces, filling its own and leaving the other as `--`. And **the un-namespaced
  test was the lock**: `test_run_offline_writes_artifacts` asserted `\ecNbursts` existed, so it
  passed because the defect was present and would have failed on the fix. Third instance in this
  repo of a test encoding "X must not exist" where X was the fix.

### Changed
- Papers restyled to traditional pre-LLM register (batch 7: `driftsearch`, `ecallisto_census`,
  `ecallisto_pipeline`, `frbperiod`, `offsets`) --- 35 of 45 papers converted. Gates clean on every
  file before review; **all five referee rounds returned revisions** (seven batches for seven).
  The MAJORs were all scope and claim-strength drift that no linter can see: `ecallisto_census`
  narrowed *"all event counts in **this paper** are synthetic"* to "in this validation" and turned
  *"This paper ... does not report a measured solar-cycle correlation"* into agentless passive,
  which reads as a claim about the literature and is false on the paper's own citation;
  `driftsearch` lost an inferential connective, so one detector on one file became *"A detector
  validated on injected tones ... fails on the real Voyager-1 data"* (any such detector), and
  softened "does not generalize **as-is**" to "as configured", offering an escape hatch the same
  paragraph rules out; `frbperiod` retitled *"A null, as expected"* to a declarative heading, in a
  paper that quotes no sensitivity whatsoever for that source; `ecallisto_pipeline` lost its
  limitations head sentence, so limitation (i) read as a description of the method, and promoted a
  LaTeX comment ("none typed by hand") into a body claim the paper falsifies one paragraph up.
- `driftsearch`'s title drops "Honest". The conversion removed all six body instances, leaving it an
  orphan; in journal register it is a claim about the authors, not the result, and `\shorttitle`
  already omitted it.
- Recurring mechanical signature this batch: **collapsing an em-dash appositive to a comma
  re-attaches it to the nearest noun.** In `offsets` it made the 24.0x excess a property of the
  Rayleigh expectation rather than the ratio of the two fractions; in `frbperiod` it let "on 19
  bursts" qualify the agreement claim rather than the periodogram. Grep for it in future batches.

### Fixed
- **Three wrong span claims in the solar-burst papers, two of them in abstracts.** Each was a
  round-number restatement of a macro pair that nobody re-derived after the runs settled.
  `triangulate` called its 0.125--1.975 MHz band "two decades of frequency" in four places; it is
  **1.20** decades (its distance range, 15.3--106.1 R_sun, is 0.84). The text now cites
  `\triFlo`/`\triFhi`, so the span cannot drift from the run again. `type3synthesis` claimed
  "more than three decades" for 0.125--78.94 MHz, which is **2.80**; the same paper counted "four
  single events" and "four regimes" where its own committed table has three, and its abstract
  claimed validation of a coronal density model the triangulated band never touches --- that leg
  is entirely interplanetary, so only Leblanc is tested, now stated in the limitations too.
- `solarbursts` printed `$0.137\,c$` for the harmonic/1x grid point three lines below
  `\sbSpeedC` = **0.1347** for the same point; that value now cites the macro. The two flanking
  grid values (0.086, 0.272) come from the same superseded run and are **not** silently patched:
  `exciter_speed` needs the raw ridge and only its summary is committed, so they cannot be
  recomputed from evidence. `survey/solarbursts-findings.md` tabulates the whole stale-vs-committed
  set (R^2 0.90 vs 0.811, drift -3.3 vs -2.55 MHz/s) and records the fix as emitting the grid on
  the next real run. Inventing two numbers to match would have been the worse repair.

### Changed
- Papers restyled to traditional pre-LLM register (batch 6: `solarbursts`, `swaves`, `windwaves`,
  `triangulate`, `type3synthesis`, `rmsky`, `sourcecounts`). Gates clean on every file before
  review; **all three referee rounds still returned revisions** (six batches for six), and all
  twenty-seven conversion defects are fixed here.
  The MAJORs: `sourcecounts`'s abstract stopped predicating "not a new measurement" on *this work*
  and attached it to the published counts it agrees with; `windwaves` turned "a recover-a-known"
  (the class of exercise) into "recovering a known result" two lines after saying its speed came
  out **below** published values, collapsed an appositive so the fitted speed attached to the
  *outer point* rather than the whole track, and downgraded "The honest limitations are
  substantial" to "Several limitations apply" in the paper whose dominant systematic exceeds a
  factor of two; `solarbursts`'s limitations list lost its head sentence, so a factor-of-3
  systematic now elaborates a positive clause.
  **One mechanical transformation had a signature**: `--- but they ---` promoted to a full stop in
  both `rmsky` and `sourcecounts`, leaving "They" pointing at the positive clause in each. A
  referee flagged the *transformation* rather than the instance, which made it greppable; it was
  confined to those two.


### Fixed
- **The two slices no guard could protect are now namespaced.** The audit's own output pointed
  at them: of the eight slices emitting no `*Source` macro, `vlass` and `type3synthesis` were the
  two with real numeric hazards, so `preserve_live_macros` could not even tell their run modes
  apart. `vlass`'s census counts are now `vlassReal*`/`vlassSyn*` (the live clobber was
  `\vlassNconfirmed` 2 -> 0 and `\vlassArea` 703 -> 0) and `type3synthesis`'s ladder is
  `synReal*`/`synSyn*` (corona speed 0.1347 real against 0.3002 offline). Both now emit their
  provenance, and `type3synthesis` records a `source` in its metrics --- it had none, which is why
  it was the single results file a forced offline rebuild could still change. Verified: a direct
  offline run leaves both papers' macros byte-identical, and the audit reports 0 for each.
- **The audit had two bugs of its own, both found by using it.** `re.IGNORECASE` made its
  `Syn|Real` marker match the `syn` *prefix* of every `type3synthesis` macro, silently exempting
  the whole slice --- a detector reporting a clean bill for a slice it never looked at. And it
  flagged the `*Source` macro itself, which is *supposed* to differ between modes since it is what
  the merge guard reads. Now case-sensitive, requiring an upper-case component boundary, with
  provenance macros exempt.


### Added
- **`scripts/audit_macro_namespaces.py`** --- the repo-wide grep for mode-dependent macros that
  are not namespaced. It runs every slice's offline leg into a throwaway directory and diffs the
  macros against the committed ones, so it measures the hazard rather than pattern-matching for
  it. **34 of 42 slices carry at least one**; about twenty change a *numeric* value under a
  shared name (`\hiVflat` 257 -> 231, `\rmRatio` 5.4 -> 8.4, `\vlassNconfirmed` 2 -> 0). Exit
  code 1 when any is found.

### Fixed
- **`preserve_live_macros` was defeated by `make figures`, silently, since it was written.**
  Snakemake **deletes a rule's declared output before running the job**, so the DAG removed each
  `generated/macros.tex` and the merge found no file to merge with. A forced run in the repo root
  **overwrote 511 real macro values across 39 papers** and blanked real values to `--`
  (`\ptRealNFit` 136 -> `--`), which the guard's own docstring says cannot happen. Direct
  invocation was protected throughout, which is why every previous test of the guard passed.
  The DAG now builds into `build/figures` and writes nothing into `papers/` or `results/` --- what
  a smoke build should always have done. Verified: a forced full `make figures` now leaves every
  committed macro file, results JSON, figure and CSV untouched, while still producing all 43
  artifacts under `build/`.
- **`preserve_live_macros` now also refuses a synthetic-over-real overwrite**, not just
  synthetic-over-real-with-a-placeholder. It reads the provenance every macro file already
  carries (the `*Source` macro) and, when a synthetic run meets a real file, keeps the existing
  values while still adding any macro the file lacks. This neutralises all 34 audited slices at
  once; per-slice namespacing remains the clean fix and is still worth doing.


### Fixed
- **`pte2`'s abstract asserted the opposite of its own body.** It listed "three facts undercut a
  giant-pulse interpretation" and gave as (i) that the excess does not correlate with spin-down
  luminosity, "contrary to the giant-pulse--Edot expectation" --- quoting the *raw* Spearman
  (-0.027), which the body establishes is count-confounded. Controlling for pulse count gives
  rho = **+0.098**, in the **predicted** direction, and the flagged sources' median log Edot sits
  **0.45 dex above** the rest, also predicted; neither is significant. The abstract now states
  the confounding, the partial correlation and the median offset, and concludes that the data
  neither support nor contradict the expectation --- which is what the body always said. No
  number changed.
- **`junodam` denied "an ~180x intrinsic rise" four lines after rendering `\jdNearFarRaw` =
  196.2.** The 180 was hand-typed and stale; now macro-backed.
- **`junodam`'s recover-a-known was not in committed evidence and its recovered value was
  wrong.** The paper hand-typed "injected 8.75, recovered 7.2" while `expected_contrast` was
  `null` and `\jdExpContrast` rendered `--`. Re-running the offline leg into a tmpdir gives
  8.75 -> **6.95**. Both are macros now, and the paper states that the injected contrast is
  several times the contrast measured on the real orbits, so the round trip cannot test recovery
  near unity.
- **`junodam` had a live macro clobber.** `io_contrast`/`expected_contrast` mean different things
  in the two run modes and shared one macro name, so an offline rebuild would have written the
  synthetic recovery into the macro the paper uses for the real measurement --- the `\tiiNEvents`
  shape, which `preserve_live_macros` cannot arbitrate because both runs write real values. Now
  `\jdRealContrast`/`\jdSynContrast`/`\jdSynExpContrast`, with a test that reproduces the clobber.
  Two existing tests asserted the un-namespaced names, i.e. they were pinning the defect.
- `junodam` now cites `louis2021jgr`, which sat in `refs.bib` uncited while measuring the
  latitudinal beaming this paper invokes as its preferred interpretation.


### Changed
- Papers restyled to traditional pre-LLM register (batch 5: `pte2`, `junodam`, `rmdipole`, `hi`,
  `glitchpop`, `vlass`). Gates clean on every file before review; **all five referee rounds still
  returned revisions** (five batches for five), and all twenty-two conversion defects are fixed
  here.
  **The first defect in this campaign where a style edit changed a FACT rather than a claim's
  strength**: folding a run-in label into a topic sentence, the `glitchpop` conversion wrote
  "the full **post-2018** catalogue" where the code (`glitchpop.py:283`) reclassifies on "the
  full catalogue". That describes a disjoint-epoch comparison never run, and it makes both
  derived quantities incoherent -- `newly classifiable` is a union operation and the single flip
  is recorded as `n_pre: 8 -> n_now: 9`, i.e. cumulative. No number, macro or citation moved, so
  no mechanical gate could see it.
  Also: `pte2`'s "Honest bottom line." became `\subsection{Summary}`, promoting to the paper's
  designated summary a paragraph quoting only the naive 19% incidence -- the endpoint its own
  referee round established is not the answer; `hi`'s closing paragraph lost the word
  *validation*, so its contribution class appeared nowhere after Results; `junodam` lost the only
  heading naming a limitation, in a paper whose one positive claim is explicitly an upper bound;
  and `vlass` swapped "pipeline" for "selection" in a sentence whose predicate is false of the
  selection alone.

### Fixed
- `glitchpop`'s abstract had lost "among these" (the denominator: the 1.55 expected false
  positives are 31 x 0.05) together with the two words distinguishing the *quasi-periodic
  fraction* from *individual members*, in a single edit.
- `rmdipole`'s abstract had split its significance claim from the caveat that kills it, so the
  formally significant p = 0.001 dipole could be excerpted without the clause establishing that
  the rejection is carried entirely by the top 1% of |residual|.


### Added
- **`report.preserve_live_results` / `report.write_results`** --- the results-JSON counterpart of
  `preserve_live_macros`, wired into **all 49** `results/*.json` write sites across 46 modules.
  A results file may now only gain information: a real result is never replaced by a synthetic
  one, a partial real re-run retains the other leg's keys rather than dropping them (the
  `torchfdmt` GPU-column case), and any key carried across runs is listed under `_merge` so a
  spliced row cannot look like one invocation. Sixteen tests.
  **Verified against the documented incident**: `typeii.run(".", offline=True)` --- the command
  CLAUDE.md records as deleting 3429 lines and flipping `is_real` True->False --- now leaves
  `results/typeii_metrics.json` byte-identical, while the offline run still writes its `tiiSyn*`
  macros. A *forced* full `make figures` in the repo root changed **25** committed results files
  before and **1** after; that one carries no provenance marker at all.
- `make guard-real` now reports how many results files carry no `source`/`is_real` marker (76
  today), so the gap neither it nor `preserve_live_results` can cover is visible rather than
  silent.

### Fixed
- The first cut of the results guard protected **3 files out of 25** because it keyed on
  `is_real`, which most slices never set. It now uses the `source` field that
  `guard_real_results.py` already treats as authoritative, and counts a *mixed* source
  ("synthetic recover-a-known + real RACS-mid epoch pair") as real --- otherwise
  `stokesv_discovery` was left unprotected by the very string that documents its real leg.
- Seven `results/*.json` writers outside the `*_metrics.json` naming convention were missed by
  the first pass, including `hi.py`'s `rotation_curve.json`. A forced offline rebuild overwrote
  it with synthetic output and **`make guard-real` failed** --- caught only because the forced
  rebuild was run as a test rather than assumed to be unnecessary.


### Fixed
- **`torchfdmt` quoted "~24x" for a ratio its own adjacent macros make 29x.** The abstract and
  benchmark section rendered "~24x (44.12 -> 1.5 s)"; 44.12/1.5 = 29.4. The 24 traces to a
  superseded 36.1 s CPU timing that appears nowhere in `results/`, in a file whose header claims
  every number is `\input` from the pipeline. Fixed at the source: `_write_macros` now
  **derives** `\spBruteSpeedup` and `\spFdmtSpeedup` from the committed timings, so no ratio in
  this paper can be hand-typed or drift again; a CPU-only run emits `--`, which
  `preserve_live_macros` will not write over a real value. Three regression tests, one checking
  the macro on disk against the JSON on disk.
- **`torchfdmt`'s real-data DM recovery was quoted as a bare "0.3%".** The FDMT butterfly indexes
  rows by whole samples of dispersive delay, so the DM is quantised rather than fitted: at this
  file's band and sampling one row is 0.0627 pc/cm^3, making the 0.18 offset **2.9 trials**. The
  abstract's "within delay quantisation" frame belonged to the synthetic leg and was inheriting
  to the real one. The real leg now records `real_dm_step_pc`, the paper says the DM is
  grid-quantised, and it quotes `\spRealSnr` (6.0) -- generated and used nowhere, despite being
  the number that says how well the peak can be localised.
- **`torchfdmt`'s benchmark row is a splice of two invocations.** `git log -p` shows a CPU-only
  run with GPU values patched in later while `device` stayed `"cpu"`; the code's own device
  logic means a `--device cuda` run writes both columns *and* sets `device: "cuda"`, so no single
  run can produce the committed row. The conclusion survives on same-run numbers (3.3x), so this
  is provenance: the paper now says the columns come from separate invocations, and
  `bench_devices` is recorded going forward.
- **`results/singlepulse_metrics.json` labelled its synthetic block with the real file's name.**
  `source` sat directly above `recovered_dm: 56.63`, a synthetic value (the real one is
  `real_recovered_dm: 56.59`). `source` now names both legs and states the key convention.


### Changed
- Papers restyled to traditional pre-LLM register (batch 4: `typeii`, `vgpra` main+note,
  `spectra` main+note, `frbstats`, `peaked`, `torchfdmt` --- nine files). Gates clean on every
  file before review; **all five referee rounds still returned revisions**, and all twenty-eight
  conversion defects are fixed here.
  **The sharpest case in the campaign so far: `typeii`'s "100%" changed meaning without changing
  value.** At HEAD the literal sat inside a disavowal ("This curve --- not a headline `100%' ---
  is the honest characterization"); after conversion it was an assertion ("recovered at 100%
  completeness"). The numeric multiset is identical, so the guard passed --- and the disavowal
  was not stylistic, it was forced by GATE-2 finding R2 because the 1.0/1.0 is an easy-synthetic
  ceiling. Also: `vgpra` lost "The result is a controlled null" and the word *controlled* fell to
  zero occurrences while its own note still used it; `peaked` lost the only *comparative* form of
  a robustness claim that survives twice in unscoped form; `spectra` and `frbstats` each lost the
  single skim-level carrier of a negative result.

### Fixed
- **A style rule that was wrong.** The `peaked` conversion stripped `\emph{}` from the mission
  name *Fermi* to stay under an `\emph`-per-kw threshold, but the style guide's own reference
  says italics are *for* mission names. Fixed with `\textit{Fermi}`. When the lint and the guide
  disagree the guide wins; when the guide and a science gate disagree the science gate wins.
- `papers/peaked/arxiv.yaml` re-derived after the conversion. Its overridden arXiv metadata
  abstract still asserted "far more robust ... and alpha high is TGSS-independent" where the
  paper now says "more robust ... because alpha high is TGSS-independent"; a metadata override
  that no longer matches its paper is worse than none.


### Fixed
- **`torchdsp` claimed its science ran on the GPU; the committed evidence says CPU.** Settled
  from git history: the commit that created `benchmark_device` (602e0ca, 2026-08-10) states
  that one field labelling both legs meant "a real GPU benchmark could only be stored by
  mislabelling the CPU-run CHIME science as GPU". `results/torchdsp_metrics.json` agrees
  (`device: cpu`, `benchmark_device: cuda`). Exactly four numbers are GPU-backed, all wall-clock
  timings. The paper's "run entirely on the ROCm GPU", "validated end-to-end" and "runs
  identically everywhere" are corrected to what the evidence covers, and the same wrong
  attribution is fixed in `README.md` (2 rows) and `survey/torchdsp-findings.md` (3 places).
  Cross-device numerical agreement is not measured and the paper now says so.
- **`stokesv` printed a reproduction command that destroyed committed evidence.** `--offline`
  defaulted to False and `--out` to `"."`, so the command in the paper ran the *real* path from
  the repo root. `--out` is now required, and the synthetic figure is no longer drawn on the
  real path at all -- previously a real leg that raised partway (a missing CASDA credential is
  enough) left synthetic output standing where the committed two-panel real figure had been.
  Two regression tests.
- **`southern`'s macros were unmerged AND un-namespaced**, the combination that ships a wrong
  number rather than a hole. One `make figures` would have replaced 1545 real matches with the
  synthetic field's count and written `\soCallTried{0}` over a real 50. Counts are now
  `soSyn*`/`soReal*`, the real-only Callingham validation defaults to a `--` placeholder rather
  than 0, and the writer merges via `preserve_live_macros`. Verified empirically: an offline
  rebuild now fills `\soSynNmatched{40}` and leaves `\soRealNmatched{1545}` untouched.
  The existing `test_run_offline` had asserted the *un-namespaced* names, i.e. it was pinning
  the defect -- another instance of a test locking a bug in.
- **Audited the `preserve_live_macros` invariant across every slice, not just the ones in
  hand.** Four more were missing it: `atlas3i` (three macro writers, and it heads the
  submission queue), `driftsearch`, `peaked`, `stacking`. All fixed. **42 of 42 slices now call
  it**, against 5 of 42 at the 2026-08-12 audit. The CLAUDE.md lesson stands: a rule stated in
  that file is not a rule the repo follows until something greps for it.


### Changed
- Papers restyled to traditional pre-LLM journal register (batch 3: `rmstructure`, `torchdsp`,
  `stokesv`, `southern`, `lpt`, `skr`). Lint 3/2/3/3/3/2 HIGH -> 0 on all six; diff-guard clean;
  triage unchanged. **All six referee rounds still returned revisions**, and all twenty-nine
  conversion defects are fixed here. The severe ones were single deleted qualifiers that invert
  a headline: `skr`'s abstract lost **"raw"** from "a raw near/far ratio of 3.33" (the paper
  exists to show that ratio is a 1/r^2 artefact collapsing to 1.39); `lpt`'s abstract lost
  "own", "already" and **"per-value provenance"** -- the contribution named in its own title --
  from its precedence concession to Rea et al.; `rmstructure` had its quality-flag result
  demoted from "a lesson promoted to the main result" to a methods aside across four edits;
  `southern` *gained* an unsupported "directly" on a fitted extremum and had "The recovery"
  become "The recovery **fraction**", a quantity the pipeline cannot produce; `torchdsp` lost
  "as of this writing" from a novelty claim its own Introduction contradicts, and fused a
  provenance sentence so synthetic numbers inherited real-archive attribution.

### Added
- `traditional-style` skill rule: **run-in labels often carry ranking.** Dissolving a
  `\textbf{Label.}` opener into a topic sentence, or promoting it to a `\subsection{}`, changes
  how the reader weights what follows -- batch 3 turned "An honest null:" into a heading naming
  an attempted measurement as though it were a claimed one, and lost a "(a near-null)" marker so
  a paragraph read as a detection for three sentences.


### Fixed
- **`frblens` quoted its own retracted limit as its first conclusion.** Discussion item (i)
  read "not lensing-contaminated at the $\gtrsim$9\% level" as a hard-coded literal --- that
  is `\flRealLimitNaive` = 0.091, the count-based value this slice's referee round retracted
  in favour of 0.368, and which the Results paragraph calls "wrong by a factor of four, in the
  optimistic direction". Now macro-backed so it cannot drift again.
- **`frblens`'s printed equation was not reproducible.** `\flRealEpsSum` (8.3) and
  `\flRealEpsMean` (0.24) were computed over 34 sources while the search ran on 33 and the
  limit divides by the restricted sum, so the page showed `2.996/8.3 = 0.3683` where that
  quotient is 0.3624. The extra source is FRB20210601A: `injection_efficiency` accepts any
  train with >= 5 bursts, the search also requires a span > 2 d. Only the macros were wrong ---
  the limit was always computed on the searched subset --- so `_write_macros` now restricts to
  `m["rows"]` and the macros were regenerated from the committed metrics with no re-run.
  8.3 -> 8.1, 0.24 -> 0.25; 2.996/8.1333 = 0.3684. Two regression tests, one of which checks
  the committed evidence so page and JSON cannot drift apart again.
- **`stokesv_discovery`'s `\svdRealDet` counted neither systems nor rows.** The paper said it
  "counts systems rather than rows" while the macro is 2: the four >= 5-sigma rows are the two
  CNS5 entries of one binary in a single beam, so the chain is 4 rows -> 2 entries -> 1 system.
  The sentence now states that chain and the estimator records what it counts.
- **`stokesv_discovery`'s census is a 12-arcsec peak search and the Method claimed otherwise.**
  `search_arcsec` defaults to 12 ("brightest Stokes-I pixel within that radius"); only
  `<= 0` is forced. Re-verified in the committed evidence: **all 54 quiescent rows have I > 0**
  at median I/sigma_I 2.19, where a fixed-pixel measurement would go negative about half the
  time (p = 2^-54). The Method now describes the search, the two remaining "forced" method
  claims are gone, and Honest limits states that the quiescent limits bound a beam-scale
  maximum rather than the star. The forced re-measurement remains outstanding.
- **`wdpulsar` quoted a Stokes-V limit as the depth of the Stokes-I null** in the abstract and
  conclusion (\wdRealMedianVLimit is `median_v_limit_mjy`); section 3 always labelled it
  correctly. Both now say Stokes-V.
- **`frbwait`'s table caption said "top 20 of 15"**; the table has all 15 rows above the
  completeness cut and is not truncated.

### Changed
- Papers restyled to traditional pre-LLM journal register via the `traditional-style` skill
  (batch 2: `svsbi`, `stokesv_discovery`, `frblens`, `wdpulsar` main.tex, `frbwait` --- the
  first three previously refereed, so their science was frozen). Lint 3/3/2/2/2 HIGH -> 0 on
  all five; diff-guard clean; triage unchanged. **All five referee rounds still returned
  revisions**, and all sixteen conversion defects are fixed here: `frbwait`'s Conclusions had
  lost the clause scoping "clustering is ubiquitous" to its own census (3 of 15 sources),
  `wdpulsar`'s abstract had collapsed a two-of-five control record into an unqualified
  re-detection and had recast a VLASS catalogue cross-match as a photometry validation,
  `svsbi` had promoted a `though` to a full stop so a power caveat stopped bounding a
  calibration pass with a 0.004 margin, and `frblens` had widened a novelty claim by dropping
  "observational".

### Fixed
- **The diff-guard can be satisfied by relocation.** Shortening `frblens`'s abstract, the
  style pass moved `$\pm$230\,s` into the Discussion, preserving the numeric multiset by
  construction so the guard passed --- and the number landed in a clause stating the Roemer
  term is "~30x our tolerance ... reaching +/-230 s there", where the ratio is 46, and which
  the Method attributes to a different delay. The abstract meanwhile lost the reason
  barycentring is mandatory. Both restored; the skill now forbids relocating a numeric literal
  across sections. Multiset equality is not content equality.

## [1.8.0] - 2026-08-22

### Added
- `lptxray` (plan 93): uniform X-ray cross-match of the 16 LPTs and the 56 white-dwarf-pulsar
  candidates against 5XMM-DR15, 2RXS, eROSITA-DE eRASS1 and CSC 2.1.1, with the
  chance-coincidence rate **measured** two independent ways (field density, and rigid
  position-shift trials computed offline from one cached cone per position; both ~1%).
  Among the candidates, **accretion strongly predicts X-ray brightness -- 20/21 accreting
  (95.2%) against 3/35 others (8.6%), ratio 11.1** -- and survives a common-footprint cut
  (14/14 vs 3/24) and an all-sky 2RXS-only cut (18/21 vs 0/35). Set against `wdpulsar`'s
  radio photometry at matched depth (0.818 vs 0.825 mJy), **accretion does not predict radio
  loudness: 0/19 accreting and 0/30 others**, though this constrains persistent emission only
  and cannot exclude LPT-like pulses. Three X-ray associations absent from Pelisoli's
  compilation, including the confirmed WD pulsar J1912-4410. Evidence:
  `results/lptxray_metrics.json`, `survey/lptxray-findings.md`.
- `lptxray.classify_coverage`: splits archival pointings into targeted and serendipitous by
  aimpoint offset rather than target name -- all 24 pointings near ASKAP J1935+2148 are
  observations of SGR 1935+2154, whose designation shares the same RA digits. Census: 12 of
  16 LPTs have archival X-ray coverage, 8 have a dedicated pointing.

### Fixed
- Plan 93's LPT leg, **cut at GATE 0 on measured evidence**: all three LPTs with a published
  X-ray detection are absent from every serendipitous source catalogue while each has pointed
  coverage taken *after* those catalogues were built (XMM 0953011101, XMM 0973390301, Chandra
  26681/26682/29265/29266). Catalogue recall on the LPT sample is **1/3**, against 16/16 on
  the candidates, so a catalogue-derived LPT detection fraction would measure catalogue
  latency rather than the sources. `lptxray.catalogue_recall` enforces this as a guard, and
  the two legs are never differenced. The plan's headline was also overstated: only one of
  the 16 LPTs is a confirmed accretor.

### Added
- `lptlin` (plan 92): linear polarization at LPT pulse positions from VAST Stokes Q/U.
  One detection -- ASKAP J1832-0911 at L/I = 10.9% (33 sigma, EVPA 98 deg) -- plus three
  non-detections. **The 10.9% is depolarization-dominated, not the source's linear fraction**:
  the published RM of +89 rad/m^2 leaves only ~9% of the intrinsic linear polarization in a
  band-averaged image, so this is a strong lower limit consistent with an intrinsically ~90%
  polarized source (matching the published 92+/-3% total). Ricean debiasing, a leakage veto,
  a `total_polarization` bound check and a `faraday_depolarization` calculator. Evidence:
  `results/lptlin_metrics.json`, `survey/lptlin-findings.md`.

### Fixed
- Plan 92's premise: an earlier feasibility check generalised from a query at **one** position
  to conclude that only EMU carries Q/U at LPT positions. Querying all sixteen shows 14 have
  Q/U imaging and VAST produces it too. The plan's scope reduction was wrong and is corrected
  in place.

### Added
- `lptspec` real run and **negative result** (plan 91, now closed): Taylor-term cutouts were
  staged for the three pulses GATE 0 cleared, and the method fails. alpha comes out -1.24,
  -2.99 and -17.30; the last is the taylor.1 image's *global minimum* at the source pixel
  (-49 sigma), and a recover-a-known over the same image puts every steady source at
  |T1/T0| ~ 1-2 against our 17. MFS fits one constant-flux power-law source across the whole
  synthesis, so an LPT pulse present for only part of it cannot be represented and the
  deconvolution absorbs the mismatch into taylor.1 -- invalid rather than merely biased, and
  not fixable by choosing brighter pulses. `results/lptspec_metrics.json` is retained as
  evidence of the failure and stamped "do not quote the alpha values".

### Added
- `lptspec` GATE 0 (plan 91): `src/jansky_research/lptspec.py` + `scripts/lptspec_gate0.py`
  decide, without fetching an image, whether an ASKAP Taylor-term in-band spectral index is
  recoverable for any LPT pulse. The short in-band lever arm makes the Taylor-1 image ~10.7x
  noisier than Taylor-0, so `sigma_alpha ~ 10.7/(S/N)`; folding in the published MT-MFS
  penalty (Rashid et al. 2024, arXiv:2405.18978) gives **3 of 7 pulses usable** at
  `sigma_alpha <= 0.3` — including ASKAP J174508.9-505149, the source whose reported
  frequency drift motivated the plan. Novelty pass clear: no LPT paper reports in-band alpha
  from Taylor terms. Evidence: `results/lptspec_gate0.json`, `survey/lptspec-findings.md`.

### Added
- `papers/lptduty/rnaas.tex`: the plan-90 note (RNAAS register, macros generated from the
  committed JSONs by `lptduty.write_paper_assets`, prose-lint and triage clean). Leads with
  the ephemeris finding: the measurement is limited by what discovery papers publish, not by
  sensitivity.

### Fixed
- `scripts/triage_papers.py` checked only `paper/main.tex`, so **no RNAAS note had ever been
  triaged** -- including wdpulsar's, at the head of the submission queue. It now checks every
  file carrying a `\documentclass`. All five previously-unchecked notes come back clean.

### Fixed
- `lptduty` GATE-2: the Poisson upper limit used `2.996 + k`, correct only at k=0. The exact
  one-sided 95% limit is `0.5*chi2.ppf(0.95, 2k+2)` (2.996, 4.744, 6.296 for k=0,1,2), so the
  three sources *with* a detection had limits 19-26% too tight -- and the write-up quoted only
  the unaffected zero-detection limits, so it was invisible. Corrected to 0.0521, 0.0668,
  0.1431 (and 0.2965 in the phase leg); `poisson_upper_95` now pins the exact values.
  Also: the detection criterion now applies `lptv`'s leakage veto (|V| > 0.006|I|) as its
  docstring already claimed (no committed number changes), and GATE 0 records the family-wise
  alpha it actually tests against.

### Fixed
- `data/lpt_sample.csv`: five errors found by the plan-90 ephemeris audit and verified
  against the arXiv record before changing anything -- three wrong discovery arXiv ids
  (GCRT J1745-3009, GLEAM-X J162759.5, GPM J1839-10 all cited unrelated papers), ASKAP
  J1832-0911's period derivative (9.0e-12 -> 9.8e-10, two orders of magnitude) and its
  period (2656.2554 -> 2656.247, the value its own cited paper gives). The three arXiv ids
  had propagated into `papers/lptv/refs.bib` notes; triage's Crossref DOI check could not
  see them because the titles and DOIs were correct. No committed `lptv` number changes.

### Added
- `lptduty` increment 3 (`scripts/lptduty_phase.py`, `results/lptduty_phase.json`):
  phase-resolved separation of f_active from the in-period duty cycle. Validated by
  reproducing exactly the pulse phases (0.489, 0.951) that `lptv` published independently.
  Seven of ten sources are excluded for having no published reference epoch at all; of the
  three attempted, one supports the single-window model (f_active ~ 0.06).

### Added
- `lptduty` GATE 0: novelty pass (clear -- no class-wide LPT activity-fraction constraint
  exists; the field review says the selection function is unmodelled) and an aliasing test
  (`scripts/lptduty_gate0.py`, `results/lptduty_gate0.json`). The aliasing test invalidates
  two of ten constraints: GPM J1839-10's snapshots are not uniform in pulse phase
  (Kuiper V 0.275, Bonferroni p 0.012), and ASKAP J142431.2-612611's period is quoted too
  imprecisely for the test to run. Both verdicts are stamped per source into
  `results/lptduty_metrics.json` (`constraint_valid`, `phase_sampling.verdict`).

### Changed
- **Airflow 3.** The container base moves 2.9.3 -> 3.3.1, the `airflow` extra to
  `>=3.3,<4.0`, DAG imports to the Task SDK (`airflow.sdk`), and compose from `webserver`
  to `api-server` (removed in 3.0) with an explicit `EXECUTION_API_SERVER_URL` and
  SimpleAuthManager (`airflow users create` belonged to the removed FAB auth manager).
  Verified by running the DAG: `airflow dags test ecallisto_ingest 2011-09-14` completes
  with all 7 tasks successful under 3.3.1.

### Fixed
- DAG `ecallisto_ingest`: `scan_station` now returns JSON-safe values. Airflow 3 serialises
  XCom as **strict** JSON, so a NaN in the returned row raises "Out of range float values
  are not JSON compliant" and fails the task; Airflow 2's encoder allowed it. The DAG parses
  cleanly either way, so this only appears at runtime on stations whose scan yields a NaN --
  found by running the upgrade rather than by CI, which never exercises this stack.

### Added
- `lptduty` (plan 90, increments 1--2): constrains how often each LPT is on, from the 647
  measured VAST snapshots `lptv` already committed -- no new data. Reports
  `p = f_active x (w + T) / P` (whose factors are not separately identifiable without an
  ephemeris) with an efficiency-weighted denominator rather than an epoch count, and an
  efficiency floor so shallow epochs cannot sum into sensitivity that does not exist.
  Three sources give p = 0.011--0.046; seven give 95% limits of 0.029--0.107 at 5 mJy.
  Evidence: `results/lptduty_metrics.json`, `survey/lptduty-findings.md`. GATE 0 and the
  aliasing check are open -- not paper-ready.

### Added
- Plans 91--93, the follow-ups Rose et al. 2026 suggests, each gated on a feasibility check
  run before the plan was written rather than after: `lptspec` (in-band spectra of LPT
  pulses from the ASKAP Taylor-term images -- 98 `taylor.1` products confirmed present at an
  LPT position, so no sub-band re-imaging is needed); `lptlin` (linear polarization --
  **scope reduced**: a CASDA query showed RACS carries no q/u at these positions, only EMU
  does, so the intended all-sky census does not exist to be done); `lptxray` (does accretion
  predict radio loudness -- cross-match the LPT catalogue and the Pelisoli candidates against
  public X-ray catalogues, with a measured chance rate).

### Fixed
- `lptv` refs: `rose2026` upgraded from an arXiv e-print to the published record
  (Nature Astronomy 10, 1166--1178, `10.1038/s41550-026-02882-x`), verified against
  Crossref. It is the discovery/classification paper for ASKAP J174508.9-505149, one of
  the census's three RACS detections, and `lptv` is in the submission queue.

### Added
- `plans/90-lpt-duty-cycle.md`: quantify how often LPTs are on, from the 966 committed VAST
  snapshots — Rose et al. state the switch-off behaviour qualitatively for one source and
  nobody has measured it across the class. Plan records up front that only the product
  (active fraction x in-period duty cycle) is identifiable without phase information, and
  that the denominator is efficiency-weighted exposure rather than epoch count.

### Fixed
- `frbstats`: annotate the repeater/one-off property dicts as `dict[str, np.ndarray]`.
  A comprehension over literal keys infers `dict[Literal[...], ...]`, which dict
  invariance rejects against `compare_populations`' parameter — an error mypy >=2.3
  reports and 1.11 did not, and the only thing blocking the dev-tooling bump (#139).

### Fixed
- CI: the tectonic package bundle is now cached (`actions/cache` on `~/.cache/Tectonic`,
  rolling key) in both `paper.yml` and `release.yml`. Every run previously re-downloaded
  the whole bundle from `relay.fullyjustified.net`, which 429s under load — it failed the
  `paper` and `release` workflows on 2026-08-19 (the v1.7.0 tag push, so the release had
  to be created by hand) and flaked twice the day before. After one successful run the
  compile step needs no network. Restore and save are split so the save runs on
  failure too — a 429-interrupted run still downloads part of the bundle, and keeping
  that partial cache lets successive runs converge instead of each starting cold.
  The cached path is `~/.cache/tectonic` (lowercase; `~/.cache/Tectonic` is also listed
  for older builds) — the first attempt named only the capitalised form, which cached
  nothing while the step still reported success, because `actions/cache` logs a missing
  path as a warning rather than an error.

## [1.7.0] - 2026-08-19

### Changed
- Papers restyled to traditional pre-LLM journal register via the `traditional-style` skill
  (batch 1: atlas3i main + RNAAS note, dr20radio, lptv, innerrc — the submission queue).
  Prose only: every macro, citation, and numeric literal is multiset-identical to its
  pre-conversion state (`prose_lint.py --diff-guard`). Each conversion then went through a
  `paper-referee` round, all four returning *minor revision* for scope drift the guard
  cannot see (deleted qualifiers, an appositive that re-pointed a claim at the wrong
  subject, negotiated verb strengths, ranking markers) — all restored; see
  `survey/stylecorpus-findings.md`.
- `wdpulsar` figure: the left panel is now real data (forced I/sigma over the 442 usable
  candidate epochs against the unit normal), replacing the synthetic injection fixture —
  the referee's deferred suggestion. The synthetic panel remains as the offline fallback.

### Fixed
- `CHANGELOG.md` structure: entries from seven PRs had accumulated *inside* the intro
  sentence rather than under `## [Unreleased]`, so `scripts/next_version.py` reported
  "nothing to release". Entries relocated and regrouped; the release recipe works again.

### Added
- `traditional-style` RNAAS extension (plan 89 follow-on): genre baseline from 391
  pre-LLM arXiv-deposited notes (`results/stylecorpus_rnaas.json`; the main corpus's
  refereed filter had excluded all 1,035), style-guide §8, `prose_lint.py --file/--genre`,
  and the wdpulsar RNAAS note converted through the full gate sequence (referee round
  restored two qualifiers the multiset guard cannot see; blind A/B passed).
- `traditional-style` skill + `style-editor` agent (plan 89, Stages 2–4): 2,000-paper
  pre-LLM corpus acquired (manifest committed), fingerprint engine + corpus-vs-us delta
  (`results/stylecorpus_{fingerprints,selfscan,manifest}.json`), corpus-derived style
  guide with before/after conversions, and `prose_lint.py` (corpus-percentile lint +
  prose-only `--diff-guard`).
- `stylecorpus` (plan 89, Stage 1): scoping of the pre-LLM radio-astronomy literature
  (1933–2021) for the `traditional-style` skill — `src/jansky_research/stylecorpus.py`
  (era strata, ADS/arXiv query builders, stratified size sampler, bootstrap total-PDF-size
  estimator; network runners no-cover), `scripts/style_corpus_scoping.py` (resumable CLI),
  committed evidence in `results/stylecorpus_scoping.json` + `survey/stylecorpus-findings.md`.
- `papers/wdpulsar/rnaas.tex`: standalone RNAAS note (refereed; ~630 words) — the null
  census condensed with the control's true 2/5 efficiency attached to the f<6.1% bound.
- `lptv.summarize_vast` + `fold_phase` (tested): VAST-sweep reduction with per-detection
  pulse phases on the published J1839-0756 ephemeris; paper gains a VAST-extension section
  reporting two pulses absent from the discovery paper's observation table (one at
  interpulse phase 0.489±0.004), with unreleased-SBID bookkeeping (107 epochs never
  publicly released, verified via obscore) and the honest 2026 phase-scrambling caveat.
- `scripts/lptv_vast_real.py`: VAST-collection extension of the lptv forced-V sweep
  (12-min epochs, ~966 across 10 covered sources; full-precision epoch MJDs and per-row
  durations, fixing the round-2 timestamp-quantization finding).

### Fixed
- `wdpulsar` main.tex: two peak-search-era numbers corrected against the forced CSV
  (max candidate significance is 2.9σ not 3.3–3.9σ; brightest measurement 3.7 mJy not
  4.3 mJy) — found by the RNAAS note's referee pass.
- Triage: a deliberately-committed `*synthetic*` metrics file with `is_real: false` no
  longer flags no-evidence when the real sibling exists (the vgpra two-file design);
  16.11 h allowlisted as a Voyager literature constant. vgpra "confirms"/"proves" and
  wdpulsar hard-typed 7.4/"confirms" fixed; zero HIGH/MED across all 46 papers again.
- `lptv` round-3 referee (VAST section): 2026-epoch chain recomputed self-consistently in
  code (20 epochs, p=0.53); phase errors include exposure smearing + anchor (±0.016);
  interpulse claim downgraded to candidate on the published V/I tension; novelty pass run
  and recorded with all section evidence in `results/lptv_vast_adjudication.json`; sweep
  end date corrected (2026 July, not August); leakage veto |I| guard; boundary tests;
  three bib entries completed/corrected (lee2025, mcsweeney2025, hurleywalker2023).
- `lptv` detection provenance: all three census detections adjudicated against their
  discovery papers' e-prints — each is an epoch the discovery paper already reported
  (J1745's SBID 20398 discovery-search epoch; J1651's single archival detection of
  2024-11-21; J1839's discovery observation). Paper reframed: zero new detections, three
  validated re-measurements (`results/lptv_detection_provenance.json`).
- `lptv` round-2 referee pass (major): rescoped negative recomputed from the forced CSV
  (median I/sigma = -0.14, no |V|/I constraint on non-detections — the 1.9 / 49–114% figures
  traced to the superseded peak-search data); J165130.3 promotion re-based on locked-pixel
  flux recovery; J183950.5 bright epoch identified as the published discovery observation
  (SBID 57929, lee2025) and reframed as a guaranteed (non-blind) known-pulse recovery with
  committed adjudication evidence (`results/lptv_j183950_adjudication.json`); ASKAP
  J1832-0911 catalogue row's erroneous name corrected (coordinate always matched the VLBA
  position); stale two-detection text, figure caption, epoch-composition, limit-count, and
  phase-coverage errors fixed.

## [1.6.0] - 2026-08-13

### Added
- `scripts/triage_papers.py` — a mechanical first-pass review over all 46 papers, checking
  only defect classes found in real papers this cycle: cited `--` placeholders, un-namespaced
  mode-dependent macros, hard-typed numbers duplicating macros, DOIs whose Crossref metadata
  contradicts the entry, `author = {others}`, missing evidence, and retracted verbs. Ships
  with an evidence-alias map for pre-convention filenames and a reasoned allowlist for
  adjudicated coincidences. Final state: zero HIGH, zero MED across all 46 papers.
- `svsbi`: `prior_sensitivity` — the same census re-inferred under a widened prior box over
  multiple torch seeds (torch was previously unseeded, so NPE training was irreproducible).
  Verdict per parameter is committed: `log_Lbreak` is prior-driven (+0.96 dex when the box
  widens, against 0.02–0.14 seed scatter), `f_beam` is prior-located, only the slope median
  is data-driven. The prior box itself is now in the evidence file.
- `stokesv.measure_circular_pol`: a genuinely forced mode (`search_arcsec <= 0`) that reads
  the pixel at the propagated position. The default 12″ mode is a peak search — on blank sky
  a noise-maximum statistic (54/54 quiescent census targets had I > 0, p = 2⁻⁵⁴ for a true
  fixed-pixel measurement) — and the docstring now says which is which. A test pins the
  distinction: forced photometry on noise goes negative about half the time.
- `frblens`: `efficiency_per_source` — injection efficiency for every searched train at the
  census's own detection threshold, and `lensed_fraction_limit(..., eps_sum=)` divides by
  Σε rather than the source count.
- `fashienv`: `void_jackknife_offset` — delete-one-void jackknife on the knee offset, plus
  `n_wall`/`n_field` committed (the comparison bin held 58% of DR1 and its size was nowhere
  in the evidence).

### Changed
- **`svsbi`'s headline is retracted to a lower limit**: the one detection does not "pin"
  log L\* (the posterior piled against the prior wall and its median tracked the box); the
  census supports log L\* ≳ 13.9. `f_beam` is stated as the product of emitter fraction,
  beaming and duty cycle — the model cannot separate them. `\svbNTargets`/`\svbSource` are
  namespaced (an offline rebuild previously wrote the 400-star synthetic parent size under
  the macro the abstract used for the 38-target real census).
- **`stokesv_discovery`'s 10σ inter-epoch decline is demoted to marginal**: I and V fell
  together (V/I constant to 3.4%), the flux-scale signature; 5.2σ/3.5σ at a 3%/5% per-epoch
  scale term. The title's "Variability Recovery" is now "Recovery". One system, not two rows:
  BL and UV Cet share a beam and the photometry was bit-identical.
- **`frblens`'s limit is 4× weaker and now honest**: f < 0.37, dividing by measured Σε = 8.3
  (mean efficiency 0.24; four searched sources have ε = 0 and constrain nothing). The
  abstract's lensed-one-off framing is corrected to lensed repeaters — the ≥5-burst and
  M_max ≥ 2 cuts exclude the one-off channel entirely (30/33 sources have M_max = 1).
- **`fashienv`'s offset is reframed as an upper bound**: the void jackknife *strengthens* it
  (0.039 dex vs the 0.087 fit error — sample variance is not the limiter), so the paper now
  names the two biases that are, both acting in the signal's direction: the 1/Vmax
  uniform-density assumption, and a bounding-box comparison bin 0.010 dex from the all-sky fit.
- 33 modules now route `_write_macros` through `report.preserve_live_macros` (5 of 42 did,
  against the stated repo rule); a single `make figures` would previously have blanked every
  real macro in 33 papers. All 46 papers now carry the `\software{}` block citing
  `jansky-research` (4 did).

### Fixed
- 19 citation defects across 12 papers, each verified against Crossref/arXiv metadata and
  annotated in the .bib: five entries had the right DOI under a wrong title, most wrong DOIs
  were the *adjacent* identifier (constructed by pattern instead of looked up), one paper was
  cited under the wrong first author (Ye et al. 2016, not "Gurnett et al.") with the prose
  corrected to match, and five `author = {others}` entries are filled in from the arXiv API —
  including `fashi_groups`, which the fashienv novelty claim depends on, and `racslow2`,
  published since the entry was written and now carrying its PASA coordinates.
- Two arXiv-assembler bugs that reached a submission payload: `\rightarrow` had no symbol
  mapping so `$V=9.26\rightarrow7.11$` rendered as the nonexistent flux "9.267.11 mJy", and
  `\citealt` (textual) was deleted like parenthetical `\citep`, leaving "(catalogue, )".
- Nine hard-typed prose numbers replaced with the macros recording the same quantity;
  `rfitrend`'s empty stable-line list now renders "none" instead of the repo's
  missing-value marker "--".

### Added
- `dr20radio`: **the spectral index is now measured, not assumed** (`run_alpha`, new
  `results/dr20radio_alpha.json`). The survey-overlap band (-40 < dec <= +30) is covered by
  both VLASS and RACS, so a quasar detected twice gives a two-point alpha directly: 5,571 of
  174,171 band quasars are detected by both. That joint sample is truncation-biased flat
  (median -0.62), because at the RACS limit a steep source has already fallen below the VLASS
  limit at 3 GHz; above S_RACS > 6.2 mJy the truncation cannot operate for any alpha >= -1.5,
  and those 4,190 quasars give **alpha = -0.72 +/- 0.02**. The canonical -0.7 is right.
  The contrast evaluated at the measured index is 4.09% vs 2.82%, a gap of 1.27 pp, moving
  only over 1.25-1.29 across +/-1 sigma -- against 0.23-1.66 for the old 0 -> -1 sweep, which
  is now the range the measurement *excludes* rather than its uncertainty.
- `dr20radio`: `kaplan_meier_median` -- survival analysis replaces the completeness cut as the
  index estimator. A RACS detection with no VLASS counterpart is not missing data: it says
  S_VLASS < the VLASS limit, i.e. alpha is LEFT-CENSORED at
  log(S_lim/S_RACS)/log(nu_V/nu_R), and those are exactly the steep objects a cut discards.
  KM uses all 6,626 RACS-detected band quasars (5,571 measured + 1,055 censored) and is
  unbiased over the whole range instead of over the range a cut assumes: **alpha = -0.755
  +/- 0.012**, against -0.722 for the cut and -0.615 for the naive joint-detection median.
  It lands steeper, as it must, and inside the floor progression the cuts were converging to.
- `dr20radio`: the flux dependence re-tested with the censoring handled. Measured on
  detections only the bins give -0.46/-0.50/-0.52/-0.84; by Kaplan-Meier they give
  -0.61/-0.57/-0.55/-0.90. So up to 0.15 of the apparent trend in the faintest bin was
  truncation -- but the trend **survives**: the brightest sources are genuinely steeper. A
  flux LIMIT is converted at the faint end, where the value is -0.61 (`ALPHA_THRESHOLD_REGIME`),
  giving a gap of 1.11 pp against 1.32 at the sample median; both are reported.
- `dr20radio`: `fetch_racs_total_flux` -- RACS `total_flux_source` over the whole overlap band
  (1.56M sources), so the beam-resolution systematic is **measured rather than estimated**.
  Recomputing with integrated flux on both sides gives alpha = -0.699, a shift of +0.056
  flatward (the 2.5" beam resolves out flux the 25" one keeps). That is a sixth of the
  flux-range term, and replaces a ~0.1 estimate the previous round could only assert.
- `dr20radio`: three sensitivity checks on the measured index, each able to fail, and each
  of which did (`completeness_floor_sensitivity`, `per_epoch`, `flux_bins`). Lowering the
  assumed completeness floor to alpha >= -2 and -2.5 steepens the median monotonically
  (-0.722 -> -0.779 -> -0.856), so -0.722 is an **upper bound on flatness**; single VLASS
  epochs give -0.703/-0.781 against -0.722 for the max-of-epochs (a positively biased
  estimator for a flux ratio); and the median runs -0.50/-0.52/-0.84 across S_RACS bins, so
  it is not transferable across flux. The paper now quotes alpha = -0.722 +/- 0.015 (stat)
  +/- 0.36 (sys) -- the bootstrap SE is the smallest term in the budget by an order of
  magnitude -- and identifies a fourth, unquantified term (peak fluxes across a 2.5" and a
  25" beam, ~0.08 in alpha per 10% flux-ratio error).
- `dr20radio`: `VLASS_S_LIM_CONSERVATIVE_MJY` and `luminosity_matched_vlass_conservative` --
  the MIRROR of the RACS conservative variant, and the side that can actually move the ratio.
  VLASS's 1 mJy is a per-epoch *reliability* threshold while RACS's 3 mJy is a 95%
  *completeness* limit, so at the measured index the north is cut at 1.25 mJy against the
  south's 3.0. Equalising them takes the gap 1.27 -> 1.10 -> 0.91 pp, a 28% reduction: the
  contrast is materially but not wholly a consequence of how the two limits are defined.
- `dr20radio`: `luminosity_matched_per_source_alpha` tests the single-index approximation by
  giving every quasar its own index drawn from the measured distribution. This reports a bias
  rather than a variance on purpose: the realization spread is 0.016 pp and would make any
  breadth of distribution look harmless. The scatter shifts each fraction by ~0.3 pp (to
  3.83%/2.53%) but is nearly common-mode, moving the gap by 0.03 pp. The paper now states
  both -- the contrast survives the scatter, the absolute fractions do not.
- `dr20radio`: spectral-index sweep (`ALPHA_SWEEP = 0, -0.35, -0.7, -1.0`) on the
  luminosity-matched fractions in both legs, committed to `results/dr20radio_{north,south}.json`
  under `luminosity_matched_alpha`. This is the sensitivity test the published 5 mJy
  "robustness check" could not perform: because the common luminosity limit is the RACS one in
  both legs, raising the flux floor rescales north and south identically and leaves their ratio
  unchanged (1.4391 -> 1.4434), whereas alpha moves it 1.08 -> 1.61.
- `make release-check` — verifies a tagged release carries the hand-built papers asset. The
  CI-built papers are synthetic and must never be attached; the upload is deliberately manual.
- `papers/atlas3i/arxiv.yaml` — restores the hand-curated arXiv metadata (astro-ph.IM primary,
  astro-ph.EP cross, page counts, cleaned abstract) that a `make arxiv` run had silently
  overwritten with keyword-inferred guesses.

### Changed
- `dr20radio`: the north/south contrast is now reported as a **range over spectral index**
  (0.23-1.66 percentage points) rather than a single 4.06%/2.82% pair, in the abstract,
  Section 4.3, and the summary. The apparent hemispheric difference is principally an artefact
  of the cross-frequency K-correction, not a measured property of the two populations; at
  alpha = 0 it nearly vanishes. No astrophysical inference is drawn from it.
- `dr20radio`: two hard-typed figures replaced with pipeline-generated macros — the
  out-of-redshift-range census fraction (`\drOutZNorthPct`/`\drOutZSouthPct`) and the northern
  census fraction lying outside the RACS footprint (`\drNorthOutsideRacsPct`).
- `atlas3i`: seven referee findings applied to `main.tex` and `rnaas.tex`, including the
  drift-grid caveat on the EIRP limit; `\aiTotScans` added so the scan count is evidence-backed.

- `innerrc`: the sensitivity scan now excludes variants with a fitted parameter on a bound
  (`bound_contact`, `railed_variants`, `n_fitted`), commits the full per-variant parameter
  vector including `v_halo`/`h_halo`, and reports `chi2_per_n` plus the primary fit's 1-sigma
  rho_DM interval from its own covariance.
- `innerrc`: the anchor now scores the paper's Table-1 solution in chi2 and measures its
  outer-curve (R > 8 kpc) residual, committed as `paper_table1_chi2_per_n` /
  `paper_table1_outer` / `paper_table1_params`.

### Changed
- `innerrc`: **the compatibility-with-consensus claim moved off the variant scan and onto the
  fit's own covariance.** rho_DM = 0.24 with 1-sigma 0.16-0.31 (the scan's old 0.19-0.32 had
  its maximum set by a variant with v_bulge at exactly the 800 km/s bound; 6 of its 8 variants
  are railed, and the 2 interior ones span only 0.20-0.24). The anchor result is also restated:
  their published halo sits +1.33 sigma per point below their own curve beyond 8 kpc over 32
  points, chi2/N 3.18 vs 1.92, so the refit is not merely "another corner of a degeneracy".
- `innerrc`: three overclaims corrected against committed evidence -- the bar-region inner peak
  does *not* reproduce (-36 km/s over 28 bins inside 2 kpc; their 255 km/s peak at 550 pc reads
  214 km/s at 1450 pc), the E/W asymmetry replicates *qualitatively* (period 36% long, damping
  railed, fit seeded at their published answer), and the estimator calibration route is not an
  independent confirmation (algebraically the same statistic on a sub-sample).

### Fixed
- `dr20radio` `refs.bib`: `arnaudova2024` title was the Macfarlane et al. 2021 title. Corrected
  to "Exploring the radio loudness of SDSS quasars with spectral stacking" against Crossref
  (10.1093/mnras/stae233).
- `dr20radio` Limitations now states that the two legs are not matched-sky samples: 14.1% of the
  northern census lies north of the RACS footprint and is unobservable from the south at any depth.

## [1.5.0] - 2026-08-10

### Added
- **`report.preserve_live_macros`** — merges `generated/macros.tex` across run modes instead of
  overwriting it, so a run may only *add* information: a real value always beats a placeholder
  and never the reverse. Wired into `typeii`, `rmstructure`, `torchdsp` and `singlepulse`.
- **`papers/<slice>/arxiv.yaml`** — a tracked, hand-authored override the arXiv assembler merges
  over its auto-extracted metadata. `arxiv-submission/` stays a disposable build artifact, so
  `make arxiv` is idempotent and the human decisions are reviewable in a PR. Written for
  `innerrc`, `peaked`, `frblens`, `vgpra`, `rmdipole`, `pte2`, `rfitrend`, `glitchpop`, `typeii`.
- **`arxiv-submit` resolves `\citet` from `refs.bib`** into natbib textual form, handling
  comma-form and braced compound surnames and `and others`.
- **`rmstructure.N_SYNTHETIC_REALIZATIONS`** and the ensemble macros `\rmsSynRatioEns`,
  `\rmsSynRatioEnsSd`, `\rmsSynNReal` — the honest uncertainty on a recover-a-known over a
  correlated random field.
- **`benchmark_device` / `benchmark_hardware`** in `torchdsp` and `singlepulse` results, so a
  GPU timing run can be recorded without mislabelling CPU-run science as GPU.
- **GPU benchmarks measured** on this workstation's AMD Radeon RX 7600 XT (gfx1102,
  torch 2.12.1+rocm7.1): `torchfdmt` brute-force dedispersion **1.5 s** GPU vs 37.9 s CPU;
  `torchdsp` FFA **0.65 s** vs 7.27 s, SumThreshold **7.63 s** vs 3.04 s — confirming with
  numbers the paper's claim that its per-series loop is GPU-hostile.
- **`tests/test_report.py`** — including the collision case `preserve_live_macros` deliberately
  does *not* rescue, which is the argument for namespacing.
- New blocking validations in the arXiv assembler: unresolvable `\citet`, an abstract starting
  lowercase, macros with no committed value, non-S-parameter and wrong-port-count Touchstone,
  and a self-parse of the generated YAML.

### Changed
- **Every mode-dependent `typeii` macro is namespaced `tiiSyn*`/`tiiReal*`.** `\tiiSource`,
  `\tiiNEvents`, `\tiiCompleteness`, `\tiiPurity` and `\tiiComp*` previously meant different
  things in the two run modes while sharing one name, so merging could not arbitrate and an
  offline rebuild turned `\tiiNEvents` from 768 real observing days into 48 synthetic events —
  under prose reading *"…days, zero failures"*. `papers/typeii/main.tex` cites the namespace it
  actually means.
- **`rmstructure`'s recover-a-known no longer overclaims.** *"recovers an injected plane
  enhancement (4.64 ± 0.35 for an amplitude boost of 5)"* became *"responds to an injected
  low-latitude amplitude boost (3.15 ± 1.11 across 30 field realizations; the
  single-realization bootstrap, 4.64 ± 0.35, understates that scatter threefold)"*. The
  bootstrap resamples within one field realization and so measures sampling noise rather than
  realization variance; the default seed sat 1.3σ high. *"…for an amplitude boost of 5"* is gone
  because the statistic is a band-average over |b| < 10 deg of a profile 5 deg wide and was
  never going to equal the injected peak.
- `make arxiv` packages every paper and fails at the end with the list, instead of stopping at
  the first failure and silently skipping the rest.

### Fixed
- **Four papers' abstracts cited macros that had been blanked to `--`.** Not uncomputed —
  *clobbered*: each slice has two run modes producing different metrics (offline synthetic
  validation vs real census, CPU vs GPU), both write the same `generated/macros.tex`, and each
  emits the other mode's macros as `--`. Whichever ran last silently won. The abstracts cite
  **both** namespaces, so no single run can populate them.
  - `report.preserve_live_macros` merges instead of overwriting: a run may only *add*
    information, a real value always beats a placeholder and never the reverse. Wired into
    `typeii`, `rmstructure`, `torchdsp` and `singlepulse`.
  - **Merging alone was insufficient.** `typeii` left `\tiiSource`, `\tiiNEvents`,
    `\tiiCompleteness`, `\tiiPurity` and `\tiiComp*` un-namespaced, so both modes wrote real
    values and nothing could arbitrate: an offline rebuild turned `\tiiNEvents` from 768 real
    observing days into 48 synthetic events, under prose reading *"…days, zero failures"*.
    Every mode-dependent macro is now `tiiSyn*`/`tiiReal*`, and `papers/typeii/main.tex` cites
    the namespace it actually means. Caught by the science reviewer.
  - Recovered: `typeii` purity **1.0**, completeness **0.917**, curve 0.333/0.625/0.917/1.0 at
    SNR 2/2.5/3/4; `rmstructure` synthetic recovery **4.64 ± 0.35**.
  - `torchfdmt` and `torchdsp` still block — their GPU benchmarks need a ROCm run, which this
    machine cannot do. The merge means a future GPU run will no longer be wiped by a CPU one.
- **Running an offline mode in the repo root destroys the real results JSON.**
  `typeii.run(".", offline=True)` overwrote `results/typeii_metrics.json` with synthetic
  output — 3429 lines deleted, `is_real` True→False, `event_list` gone. Reproduced twice
  (once by me, once by the reviewer). `make guard-real` catches it only at packaging time.
  Documented in `CLAUDE.md`: run offline modes with `out=<tmpdir>`.
- `papers/peaked/arxiv.yaml`: restored the scope-limiting caveat *"a tooling and methodology
  contribution, not a discovery"*, which the first trim dropped. That is a disclaimer against
  over-reading the 6-candidate list, not connective tissue.
- `tests/test_report.py`: unit tests for `preserve_live_macros`, including the collision case
  it deliberately does **not** rescue (which is why namespacing was also required).

- **The two GPU benchmarks were measured**, on this workstation's AMD Radeon RX 7600 XT
  (gfx1102, torch 2.12.1+rocm7.1). `torchfdmt`: brute-force dedispersion **1.5 s** on GPU
  against 37.9 s on CPU. `torchdsp`: FFA **0.65 s** GPU against 7.27 s CPU, and SumThreshold
  **7.63 s** GPU against 3.04 s CPU — confirming, with numbers, the paper's claim that its
  per-series loop is GPU-hostile. **Every paper now packages clean.**
- `torchdsp`: split `benchmark_device` from `device`. One field was labelling both the science
  leg and the timing run, so a real GPU benchmark could only be recorded by mislabelling the
  CPU-run science as GPU. The results JSON now carries `benchmark_hardware` too.
- **`rmstructure` no longer overclaims its recover-a-known.** The abstract said it "recovers an
  injected plane enhancement (4.64 ± 0.35 for an amplitude boost of 5)". Two problems: the
  bootstrap resamples within *one* field realization and so measures sampling noise rather than
  the realization variance of a correlated random field, and the default seed sits high. Across
  30 realizations the recovered ratio is **3.15 ± 1.11** — 3.2× the quoted error, with seed 0's
  4.64 a high outlier. `run(offline=True)` now computes the ensemble
  (`N_SYNTHETIC_REALIZATIONS = 30`) and emits `\rmsSynRatioEns`/`\rmsSynRatioEnsSd`, and the
  paper reports it: *"responds to an injected low-latitude amplitude boost (3.15 ± 1.11 across
  30 field realizations; the single-realization bootstrap, 4.64 ± 0.35, understates that
  scatter threefold)"*. "Recovers … for an amplitude boost of 5" is gone: the statistic is a
  band-average over |b| < 10 deg of a profile 5 deg wide, so it was never going to equal the
  injected peak, and comparing them implied a target the measurement cannot reach.
- **`papers/typeii/refs.bib` `lawrance2024` corrected against Crossref** (DOI
  10.1007/s11207-024-02317-8): authors are Lawrance, Devi, Chandra & Miteva — "Moni-Bidin" was
  not an author — and it is article **75**, not page 58.


## [1.4.0] - 2026-08-07

### Added

- `dr20radio` paper (plan 88): `papers/dr20radio/` (AASTeX, 4 pp) — GATE-2 passed (no
  blockers); should-fixes applied: SDSS Collaboration byline, conservative 5 mJy RACS-limit
  variant (contrast robust: 3.45% vs 2.39%), north carton validation split by selecting
  survey (49% RACS-selected / 27% LOFAR-selected at 3 GHz), overlap-band luminosity-matched
  3.12% surfaced. Slices registered in the Makefile paper targets.
- `dr20radio` increment 3 (plan 88): luminosity-matched north/south contrast — counterpart
  fluxes through the match chain, rest-1.4 GHz common-limit fractions (north 4.06% vs deep
  south 2.82%, overlap 3.12%), `paper_assets` figure + macros from committed evidence.
- `dr20radio` southern leg (plan 88): **the categorical first** — 73,074 DR20 quasars south
  of −40° (pure LCO, confirmed in-data) × RACS-low DR1 → **3.95% radio-detected** at 5″
  (measured false-match 0.06%); overlap band 3.76% vs VLASS's 4.67%; carton validation
  split by selecting survey resolves the pooled 66%: racsradio×RACS 88.9% (the pipeline
  validation) vs lofarradio×RACS 17.2% (cross-frequency fading). Evidence committed.
- `dr20radio` increment 2 code (plan 88): the RACS southern leg — strip-wise resumable
  RACS-low DR1 fetcher (CASDA sync TAP, 1° Dec strips, verified live: 2,123,638 sources),
  `run_south` (deep-south categorical-first census + the −40..+30 overlap band for the
  VLASS cross-check + carton validation against the SELECTING survey), and the two-survey
  synthetic that models the spectral-fading blind spot increment 1 exposed. 9 tests.
- `dr20radio` northern leg (plan 88): 202,691 clean DR20 quasars × VLASS E2/E3 — radio-
  detected fraction 4.67% (any epoch, 2.5″, measured false-match 0.009%), with the
  classic high-z selection rise; the radio-carton VLASS match rate (31–35%) reframed
  honestly as a cross-frequency measurement (selection was at 144/888 MHz — the ~100%
  validation belongs to the RACS leg). E3 interim-list schema handled. Evidence committed.
- `dr20radio` increment 1 (plan 88): census tooling + synthetic recover-a-known — spAll
  quasar selection with the radio-carton circularity exclusion (verified to matter in the
  synthetic round-trip), sky cross-match, position-shift-measured false-match rates
  (density-scaling verified), Wilson-interval detection fractions, and the resumable
  spAll-lite fetcher. 7 offline tests.
- `plans/88-dr20-radio-census.md` — new slice plan (`dr20radio`): the first radio-counterpart
  census of the SDSS-V DR20 Black Hole Mapper quasars (VLASS north / RACS south). GATE 0
  done 2026-08-07: catalog pinned to the file+column level (spAll-lite, 177 MiB, `OBS`
  hemisphere flag), novelty verified (all prior SDSS×radio work is legacy-catalog based; the
  southern SDSS×RACS pairing was categorically impossible before DR20), the eROSITA-rights
  landmine scoped out (X-ray leg deferred), and the radio-targeted open-fiber carton
  circularity identified with its exclusion baked into the plan.

## [1.3.1] — 2026-08-05

### Fixed

- `atlas3i` pre-submission hardening (PRs #169–#171, entries added retroactively): the
  self-citation re-pinned to v1.3.0; pre-flight review of the arXiv package caught and fixed
  a missing figure (`main.tex` referenced `fig:survey` with no figure environment — the PDF
  printed "Figure ??"), a 404ing Choza et al. DOI, the since-minted Jacobson-Bell RNAAS DOI,
  and the FAST paper's DOI/initials; final reader-pass fixes: the cadence duration corrected
  to the paper's 30 min, the EIRP-agreement claim restated with its exact substitutions
  (δν_t = δν, β = 1) instead of "algebraically identical", and the figure title's literal
  backslash removed at the source.

## [1.3.0] — 2026-08-05

### Added

- `innerrc` GATE-2 PASS with fixes: Sofue & Kohno DOI added (10.1093/pasj/psaf114); the
  pre-calibration Table-2 offset made traceable (`median_dv_kms_fixed_sigma15` = −14.5,
  superseding a pre-Q4-fix −13.4); bulge claim scoped and macro-fed; sensitivity variants now
  commit bulge/disc params; E/W amplitude compared to their stated ≃15 km/s mid-disc figure;
  MG&D16 apples-to-apples clause; hidden-constraint and Davis+2025 phrasing softened.
- `innerrc` increment 3 (plan 86): the full-length paper (`papers/innerrc/main.tex` +
  `refs.bib`, compiles in the tectonic container) — anchor degeneracy + raw-HI4PI
  replication + estimator head-to-head, all numbers from `paper_macros` reading the two
  committed evidence JSONs (the committed-real-results pattern end-to-end).
- `innerrc` increment 2 (plan 86): the HI4PI real-data leg — 1,113 sightlines through both
  TVM estimators close the `hi` slice's caveat with a measured number (threshold reads
  +17.9 km/s high; per-estimator σ calibration gives 9.5 vs 26.6 km/s, bracketing the
  paper's adopted 15); their Table 2 reproduces at the ~4% level from raw survey data; the
  E/W asymmetry replicates qualitatively (period/phase match, softer amplitude). Fourth-
  quadrant |sin ℓ| bug and a CDS stall (now timeout+resume+retry) found and fixed. Five
  paper-style figures + evidence JSON committed.
- `innerrc` increment 1 (plan 86): tooling + synthetic recover-a-known + the **offline anchor
  result** — decomposing Sofue & Kohno 2025's own published RC tables (vendored from the arXiv
  source) validates their ρ_DM = 0.107 GeV/cm³ arithmetic exactly (through the eq.-25 4π
  convention, now covered by a test) while showing it is one corner of a broad disc–halo
  degeneracy: an unconstrained refit of the same curve fits marginally better at
  ρ_DM = 0.24, and all eight sensitivity-scan variants land in 0.19–0.32 — the consensus
  range. Gaussian-decomposition TVM (greedy + joint refinement) recovers injected terminal
  velocities on crowded synthetic spectra and measures the `hi` threshold estimator's
  documented high bias on the same spectra. Committed evidence: `results/innerrc_anchor.json`.

### Changed

- **Committed-real-results migration (the structural fix from the 2026-07-31 incident):**
  real-run outputs — `results/*.json|csv`, every paper's `figures/` and `generated/` macros —
  are now git-tracked evidence, reviewed in PRs, replacing the "regenerable artifacts are
  gitignored" policy that made the synthetic-clobber possible. The offline Snakemake DAG is
  demoted to CI smoke only; `make guard-real` continues to gate packaging. All artifacts
  committed in this change are guard-verified real (or by-design-disclosed synthetic).

## [1.2.0] — 2026-08-04

### Fixed

- **Integrity: packaged papers can no longer mislabel synthetic output as real data.** Caught
  2026-07-31 via the released `hi` paper, which quoted the offline synthetic-fixture rotation
  curve (231 km/s, ±0.3 scatter) under a caption claiming LAB data — `make papers-zip` depended
  on the offline `figures` DAG, which regenerates every slice's figures/macros from synthetic
  fixtures into the same files the papers `\input`. Now: `papers-zip` depends on a new
  `make guard-real` (`scripts/guard_real_results.py`) that fails the build if any
  `results/*.json` is synthetic-sourced (allowlist requires per-entry justification, currently
  empty); the CI release workflow no longer attaches its offline build to Releases (renamed to
  a clearly-labelled `papers-synthetic-smoke` artifact with a MANIFEST warning) — distributable
  zips are built locally from real-sourced results and uploaded manually. The v1.1.0 release
  asset's `hi/main.pdf` was regenerated from the real LAB leg (257 ± 6 km/s) and replaced, with
  the correction noted in the zip's MANIFEST.

### Added

- `atlas3i` full 1–12 GHz completion: S/C/X sweeps done (60 unique node-cadences total,
  1,938,860 raw hits → 294 on/off survivors → **0 confirmed** — the paper's full-band null
  reproduces). `survey_summary`/`survey_macros`/`survey_figure` aggregate all bands for the
  papers; satellite table extended (S-MSS, S-DARS, Ku FSS/DBS — X's only 4 survivors are
  Ku-band TV tones above the analysed passband); findings updated with the full funnel and
  measured cost (~95 h compute, ~3.7 TB transfer, one desktop); RNAAS updated to the
  full-band claim and a full-length `main.tex` added for arXiv.
- `plans/86-inner-rotation-curve.md` — new slice plan (`innerrc`): independent replication of
  Sofue & Kohno 2025 (PASJ 77, 1335; arXiv:2509.23581), the modern inner Milky Way rotation
  curve. GATE 0 done: full-text read (their 0.107 GeV/cm³ local DM density is halo-only and
  author-framed as a lower limit — the slice maps its sensitivity, not a refutation), novelty
  pass (cited, never audited), data verified (HI4PI CDS tiles ~15–20 GB for |b|≈0; their RC
  ASCII tables in the arXiv source as the anchor). Also retires the `hi` slice's documented
  2 K-threshold bias via a Gaussian-decomposition TVM head-to-head. First plan written under
  the committed-real-results integrity rule.
- `atlas3i`: S/C/X band support — the pinned node→band tables for all four receivers
  (`BAND_NODES`, with the six pairwise-duplicated bank-boundary recordings per C/X band
  deduplicated and documented in `DUPLICATE_NODES`), and `--sweep`/`sweep_summary`/
  `sweep_figure` generalised over `--band`. Extends the L-band reproduction toward the
  paper's full 1–12 GHz claim (S: 6 nodes; C: 23; X: 25 unique nodes).

### Fixed

- CI paper builds: `atlas3i` registered in the static-slice DAG (`workflow/Snakefile`, via its
  offline `--paper` entry point over the committed result JSONs), so `make figures`/CI
  regenerate `papers/atlas3i/generated/macros.tex` + the sweep figure before tectonic runs —
  the v1.1.0 `paper`/`release` workflows failed on the missing gitignored macros file.

## [1.1.0] — 2026-07-30

### Added

- `plans/85-bl3i-atlas-drift-reproduction.md` — new slice plan: independent reproduction of the
  Breakthrough Listen 3I/ATLAS GBT nondetection (RNAAS arXiv:2512.19763) from the public
  `bldata.berkeley.edu/ATLAS` archive, reusing the plan-11 `driftsearch` machinery. GATE 0 done
  in-plan (data verified public with sizes/layout; novelty pass: no independent reanalysis of
  the released data exists — ATA/FAST/MeerKAT follow-ups are all original observations).
- `src/jansky_research/atlas3i.py` + `tests/test_atlas3i.py` — the plan-85 tool: archive-index
  parser and ABACAD cadence selector (validated against a vendored real GB_ATLAS listing),
  bandpass normalisation with DC-spike excision (the artifact that defeated plan 11's teaching
  detector on real BL data), a physical-units (Hz/s) brute-force de-Doppler search with robust
  per-drift S/N and neighbour suppression, the ABACAD on/off sky-localisation filter, a
  narrowband EIRP limit that reproduces the paper's ~100 mW headline from GBT L-band
  parameters, and a fully offline synthetic-cadence round-trip (`run`) that recovers an
  injected drifting tone while rejecting always-on RFI and the DC spike. Real-data leg
  (`--node`, network + `voyager` extra) streams each 10 GB fine-resolution file in bounded
  frequency chunks and deletes scans after searching (peak disk ≈ one cadence). The node→band
  mapping (L: blc21–26 = 939–2064 MHz, etc.) was pinned by remote HDF5 header reads and is
  embedded as constants. Candidate vetting: per-survivor drift-coherence stamps
  (`vet_stamps`, inline during the search) plus a satellite-allocation exclusion axis
  (`classify_band`; Iridium/Inmarsat/GNSS downlinks defeat two-position filters by design).
- `papers/atlas3i/` — the RNAAS note (`rnaas.tex` + `refs.bib`, compiles in the tectonic
  container): the reproduction result + the filter-evasion taxonomy, with every number
  `\input` from `sweep_macros`/`sweep_figure` (new module functions reading the committed
  per-node result JSONs — nothing typed by hand). README slice tables updated. GATE-2
  review passed with fixes (threshold asymmetry implemented + honest caveats: raw-hit-count
  comparability, analysed-band scope, β=1, notch filters); distance Horizons-pinned to
  1.79801 au → final EIRP limit 99.2 mW.
- `survey/atlas3i-findings.md` — the L-band real-data result: all six nodes (939–2064 MHz)
  searched and vetted, 261 ABACAD survivors, 0 confirmed — **the paper's L-band null
  reproduces independently**, with a matching ~99 mW EIRP limit; includes the two instructive
  filter-evasion modes (sub-threshold OFF carriers, satellite downlinks) and honest caveats
  (L band only, low-drift blind spot, S/N-convention and SEFD/distance nominals).

## [1.0.1] — 2026-07-25

### Changed

- **Install: a single clone now works.** The `jansky` dependency defaults to its pinned git tag
  (`v0.2.0`) in `[tool.uv.sources]` instead of the `../jansky` path source, so `git clone && uv sync`
  succeeds with no sibling checkout — removing a JOSS review-checklist blocker and the most likely
  first-run failure for new users. Cross-repo development is now an explicit opt-in:
  `eval "$(make -s dev-env)"` puts a sibling `../jansky/src` on `PYTHONPATH` ahead of the pinned tag
  (`uv run` preserves `PYTHONPATH`; a `uv pip install -e` overlay does *not* survive, because
  `uv run` re-syncs the environment). CI/paper/release workflows drop their second
  `actions/checkout`, so CI installs exactly what a stranger gets.
  Docs updated: `README.md`, `CONTRIBUTING.md`, `REPRODUCING.md`, `docs/usage.md`,
  `joss/SUBMISSION.md`, `CLAUDE.md`.
- Documentation: added the Zenodo DOI badge (concept DOI `10.5281/zenodo.21482378`) to the README
  following the v1.0.0 release.
- Documentation: added community guidelines for JOSS review — `CONTRIBUTING.md` (how to contribute,
  report issues, and get support; dev setup; the slice pattern; the pre-PR checks),
  `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), and `docs/usage.md` (a short install-and-run
  guide). Linked from the README.
- JOSS paper (`joss/paper.md`): added a substantial-effort paragraph to the Statement of need
  (breadth across domains + the dual Snakemake/Airflow reproducibility layer + "the software is
  itself the research"). Refreshed the stale "twenty-six slices" tally → "forty slices plus a
  synthesis" in the README Results header and `CLAUDE.md` (the slice table has 41 rows).
- Papers now cite the toolkit: the `vgpra` and `spectra` RNAAS notes `\software{}`-cite
  **`jansky-research`** via its Zenodo concept DOI `10.5281/zenodo.21482378` (new `@misc{janskyresearch}`
  `refs.bib` entry). Added a related-work disclosure to `joss/paper.md` (the two notes are in
  preparation for RNAAS, distinct from the software paper), and recorded the "cite jansky-research
  going forward" convention in `CLAUDE.md` and the `research-publish` skill.
- `research-publish` skill: the readiness check now **auto-discovers every `papers/*/rnaas.tex`**
  note (was hardcoded to `frbstats`) with a per-note 1000-word sanity check, so `vgpra`/`spectra`
  and any future note are covered without editing the script.
- Added `joss/SUBMISSION.md` — a paste-ready helper for the manual JOSS submission (form fields,
  the comments-to-editor text with the substantial-effort + related-work + AI-use disclosures, a
  verified suggested-reviewers shortlist, and RNAAS Editorial-Manager cover text).
- Added `docs/faq.md` — how others use the toolkit (depend, don't fork), how the in-repo papers work
  (authorship, provenance/priority via Zenodo DOIs, what "unpublished" does and doesn't mean), and a
  licensing note. Linked from the README and `CONTRIBUTING.md`.
- **Dual-licensed the repository:** the papers in `papers/` are now **CC BY 4.0** (new
  `papers/LICENSE`); the code remains **MIT**. README, `docs/faq.md` updated to state it.

## [1.0.0] — 2026-07-21

**Initial public release.** With no prior tag this records the full toolkit as it stands rather
than a diff from an earlier version; every later section is a diff from its predecessor.

### Added

- The `jansky-research` toolkit — ~40 self-contained, reproducible research slices (a tested tool
  → real public data → adversarial science-review gate → honest AASTeX write-up), grouped by
  domain:
  - *FRB & time-domain* — `frbstats`, `frbperiod`, `frbwait`, `frblens`, `singlepulse`;
  - *Pulsars* — `pulsarspec`, `ppdot`, `pte2`, `glitchpop`, `wdpulsar`;
  - *HI & spectral line* — `hi`, `fashienv`;
  - *Solar & heliospheric* — `solarbursts`, `windwaves`, `swaves`, `triangulate`,
    `type3synthesis`, `typeii`, `ecallisto_census`, `ecallisto_catalog`, `rfitrend`;
  - *Planetary radio* — `junodam`, `skr`, `vgpra`;
  - *RM / Faraday & cosmology* — `rmsky`, `rmstructure`, `rmdipole`;
  - *Continuum & variability* — `vlass`, `vlbi`, `stacking`, `sourcecounts`, `peaked`,
    `southern`, `spectra`, `offsets`, `stokesv`, `stokesv_discovery`, `lpt`, `lptv`;
  - *SETI* — `driftsearch`;
  - *GPU / ML* — `fdmt` + `torchdsp` (device-portable pure-PyTorch DSP), `svsbi`
    (neural simulation-based inference).
- The shared slice scaffolding in `src/jansky_research/`: `data.py`, `pipeline.py`, `report.py`.
- Reproducibility paths: `make reproduce`, the Snakemake static-slice file-DAG
  (`workflow/Snakefile`, `make figures`), and the Airflow streaming e-Callisto ingest
  (`airflow/`).
- Publishing / data helpers under `.claude/skills/`: `arxiv-submit`, `casda-cutout-fetch`,
  `pull-station-data` (the `jansky-observe` rooftop-station bridge).
- **Release-versioning infrastructure** — `VERSIONING.md` (SemVer policy for this repo),
  this `CHANGELOG.md`, and `scripts/next_version.py` (recommends the next version from the
  `Unreleased` section, with its reasoning).
- `papers/vgpra/rnaas.tex` — an RNAAS short-form of the Voyager 2 PRA ice-giant rotation-period
  reanalysis (the recover-a-known → controlled-null showcase).
- `papers/spectra/rnaas.tex` — an RNAAS short-form showing raw TGSS×NVSS ultra-steep-spectrum
  selection is dominated by the TGSS flux-scale systematic (the "apparent signal is a
  systematic" cautionary note).

### Changed

- Optional GPU acceleration (`fdmt`, `sbi` extras, and the `torchdsp` slice) is pure-PyTorch and
  ROCm/CUDA-portable; the core install and CI remain CPU-only, GPU is opt-in.
- Refreshed the JOSS paper (`joss/paper.md` + `paper.bib`), `CITATION.cff`, and `.zenodo.json` to
  the current scope: "CPU-first with optional GPU" (was the false "CPU-only"), a domain-grouped
  capability list over the full >40-slice toolkit (was a stale six-module snapshot), and a
  Statement of need reframed around recover-a-known + honest-null at scale.

[Unreleased]: https://github.com/joebarbere/jansky-research/compare/v1.3.1...HEAD
[1.3.1]: https://github.com/joebarbere/jansky-research/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/joebarbere/jansky-research/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/joebarbere/jansky-research/compare/v1.1.0...v1.2.0
[1.1.0]: https://github.com/joebarbere/jansky-research/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/joebarbere/jansky-research/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/joebarbere/jansky-research/releases/tag/v1.0.0
