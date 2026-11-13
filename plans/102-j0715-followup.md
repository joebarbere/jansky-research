# Plan 102 — GPM J0715−1141: catalogue update + follow-up leads (PARKED, not frozen)

Status: **parked 2026-11-13**, waiting on full-text access. Nothing here is frozen. No real data has
been opened for this source. The controls in §3 are *candidates*: they need to be frozen
before the first real run, per the slice rules in `CLAUDE.md`.

## 0. What happened

McSweeney, Hurley-Walker, Bahramian, Horváth, Rea, Rodriguez & Pelisoli, *A long-period transient
with a Sun-like companion*, Nat. Astron. (2026), doi:10.1038/s41550-026-03016-z. It was published
2026-11-12, received 2026-01-19 and accepted 2026-10-05 (landing page). Joe found it via Hacker
News item 50465453 and the Curtin media release of 2026-11-12. Only the paywalled landing page has
been read: the abstract, the references, the data/code availability statements and the dates. The
body, tables, Methods and supplementary material have **not** been read.

Facts and their locators:

| Fact | Value | Locator |
|---|---|---|
| Period | 7.9 min | abstract |
| System | low-mass WD + G-type star (optical spectroscopy) | abstract |
| Pulse modulation | ~1.3 d, interpreted as the orbit | abstract |
| Roche lobe | timing suggests the G star does **not** fill it | abstract |
| RM | rapid changes near WD superior conjunction, consistent with coronal Faraday rotation | abstract |
| Authors' conclusion | WD-binary LPT emission "crucially depends on its interaction with the companion" | abstract |
| Class count | "at least 4 (out of 12)" LPTs are WD + M-dwarf binaries | abstract |
| Distance | "over 5000 light-years" | Curtin release, p.3 (press, not paper) |
| Discovery path | MWA GPM survey → MeerKAT localisation → Keck LRIS spectrum | Curtin release p.3–4; refs 28, 39–40 |
| Data | dynamic spectra + per-pulse properties, Zenodo doi:10.5281/zenodo.21386095 | Data availability (ref. 41) |
| Code | github.com/robotopia/j0715-1141-analysis | Code availability |

The code repository's ephemeris, read from `timing.py:15-29` at commit `5c367e6` (2026-07-16, before
acceptance). It is **provisional** because it is a code value, not a table in the paper:

| Field | Value |
|---|---|
| PEPOCH | MJD 60648.775490179614 (UTC) |
| P | 474.059212888062 s |
| Porb | 116896 s (= 1.35296 d) |
| T0, superior conjunction | MJD 60648.80 |
| Activity window | period 1.3529 d, centre MJD 60648.785, duration 6.5 h. The code comment says this is "distinct from orbital period", but the two agree to ~5 s. Read why in the paper. |
| DM | 181.77 pc cm⁻³ |
| τ_sc at 1 GHz | 2.818e-4 s |
| Position (ICRS) | 07:15:07.365 −11:41:53.745 |
| Inclination | 35° |
| Separation | 5.5 R☉ |
| Pdot | **none given** |

### Access check, 2026-11-12 (local session, not the cloud one)

- **Full text: still not readable.** nature.com responds but the article is paywalled. The landing
  page adds the figure titles: Fig. 1, localization; Fig. 2, single-pulse properties against
  orbital phase; Fig. 3, the optical spectrum against model spectra. The page names no Table 1.
  An arXiv title or name search finds no preprint. Step 1 needs Joe's institutional access or
  the PDF dropped into the repo.
- **Crossref: reachable.** It confirms the title and the seven authors exactly as in §0, and gives
  online-first 2026-11-12 with no volume or pages yet. `mcsweeney2026` was added to
  `papers/lpt/refs.bib` and `papers/lptv/refs.bib` (step 2 done; re-check volume and pages when
  the issue is assigned). It is not cited in any prose yet.
- **Zenodo 10.5281/zenodo.21386095: reachable.** Record created 2026-07-16, CC BY 4.0. It holds 76
  `dyn_spectra_N.tar.gz` (pickled dynamic spectra, ~2 GB in total), `contents.txt` and
  `toas.ecsv`. The ToA file has sha256 `613c5886…be29d18` and 227 rows: 184 MWA (154/200 MHz,
  MJD 60648.74–61050.67) and 43 MeerKAT (813/1279 MHz, MJD 60766.53–60915.44). It has per-pulse
  ToA, ToA_err, fluence and peak S/N, and no orbital phase or polarization. Only this inventory
  was made; nothing was fitted.
- **Code repo: unchanged** since `5c367e6` (3 commits, all 2026-07-16). The §0 ephemeris is
  still the latest public one.

## 1. Next steps (do these first, in order)

1. **Get the full text.** Use Joe's institutional access or the accepted manuscript, plus the
   Zenodo record. zenodo.org, nature.com, Crossref and arXiv were all blocked from the cloud
   session, but `github.com` worked. Replace every provisional value above with the paper's table
   value and its locator (table, section or equation). Record any disagreement with the code
   ephemeris; do not average the two.
2. **Settle the citation.** Run `curl -s https://api.crossref.org/works/10.1038/s41550-026-03016-z`
   to get the authoritative metadata. Add a `refs.bib` entry to the `lpt` and `lptv` papers, and to
   `lptduty` if the source enters its sample.
3. **Catalogue v4: add row 17 to `data/lpt_sample.csv`.** Set `binary_status=yes` and
   `companion=G dwarf`, with discovery DOI and flags. The schema has no field for a spin period
   that differs from the orbital period. Decide whether to add `orbital_period_s`; ILT J1101+5521
   and GLEAM-X J0704−37 have radio period = orbit, while J0715 has spin ≠ orbit.
4. **Reconcile "4 of 12 WD + M dwarf".** Our CSV labels only 2 companions as M dwarfs: ILT
   J1101+5521 and GLEAM-X J0704−37. The other WD-ish labels are "hot WD/disk", "hot WD", "possible
   WD" and "sub-stellar donor". Find which 4 objects the paper means (introduction or reference
   list; candidates include Bloot+2025, ref. 36). Then decide whether one of our labels is stale or
   the definitions differ. The paper also says 12 LPTs where we have 16, so record their counting
   rule.
5. **Re-run `lpt` and `lptv`, then update the papers.** Use `out=<tmpdir>` for any offline leg.
   Check `git status` afterwards. Update N, the split statistics and the "hard 78-min boundary
   already falsified" sentence (`survey/lpt-findings.md:97`).

### Pre-computed what-if (scratch only, nothing committed)

On 2026-11-13 the repo's `lpt.label_sensitivity` was run with J0715 appended as
`binary_status=yes` at P = 474.06 s. The baseline row reproduces the committed numbers.

| Labelling | Δlog median | p (exact) |
|---|---|---|
| unknowns excluded (primary) | 0.413 → 0.381 | 0.167 → 0.275 |
| as published | 0.176 → 0.145 | 0.526 → 0.661 |
| unknowns as binary | 0.446 → 0.413 | 0.092 → 0.172 |
| confirmed only | 0.605 → 0.509 | 0.191 → 0.324 |

A short-period WD binary weakens the "WD binaries are long-period" split under every labelling.
This strengthens the existing null and overturns no published claim. Re-run it with the paper's
period before quoting it.

## 2. Avenues for new findings (leads; GATE-0 not run on any)

Ranked by expected value. Each needs the standing GATE-0 (full-text novelty pass + data-URL
check) before it becomes a plan.

**A. Radio census of known WD + FGK post-common-envelope binaries** (the `wdpulsar` pattern with a
new parent list). This is the strongest lead. The paper shows a G star in a day-long orbit can
power an LPT. Nobody (unchecked) has asked how many *known* WD + FGK binaries at day-scale periods
are radio sources. Parent lists: the White Dwarf Binary Pathways survey (paper refs 21–22), the
SDSS PCEB series (refs 18–20), and Gaia-selected WD + MS binaries. Measurement: forced RACS I+V
photometry plus VLASS, reusing `wdpulsar`/`ucdv` machinery (blank-sky null, injection–recovery,
recover-a-known). Output: a detection or an occupation-fraction upper limit. Divide by
sensitivity (sum of ε), not by N; see the `frblens` lesson. Each target is gated by a
~20% activity window, so a single snapshot has ε ≤ ~0.2 even for an always-active system. Fold
that in before quoting any limit. Scoop vector: the RACS-low2 circular-polarization catalogue
(arXiv:2606.16182), as for `wdpulsar`.

**B. Archival ASKAP pulses phase-folded on the published orbit.** This would be an independent
certification of the orbital gating. Dec −11.7° is inside RACS and VAST coverage. Run
`lptv`/`lptduty`-style forced photometry at the position. Freeze this prediction before opening
the data: detections fall inside the 6.5-h activity window (fraction ≈ 0.20 of the orbit). Null:
the inside fraction equals the window fraction. The window was fitted on MWA/MeerKAT data, so
ASKAP epochs could not have shaped it and may certify it (oracle rule). **Precondition:**
propagate the Porb uncertainty from the paper back to each ASKAP epoch. If the accumulated phase
error exceeds about a quarter of the window, the test is void by construction, so say so and do
not run it. Unpublished pulses would also lengthen the timing baseline (no Pdot is published).

**C. A known-answer test for `lptduty`'s implied-active-fraction estimator.** J0715's
orbit-gated window gives an externally known on-fraction (~0.20, plus intrinsic intermittency).
That would be the first planted-truth check of f_active = p/((w+T)/P) on a real source. It is
only possible if B finds enough detections; otherwise record it as not possible.

**D. Spin vs beat.** The beat period is 1/(1/P − 1/Porb) = 475.99 s, against 474.06 s for the
pulse period, a 0.4% difference. In AR Sco the radio pulses appear at the beat period (from
memory; check against Marsh+2016, paper ref. 13). Whether J0715's pulses are at spin or beat
bears on the authors' interaction claim. First read what the paper says. If it is left open, ask
whether the Zenodo ToAs can discriminate. Archival ASKAP data alone almost certainly cannot.

**E. Spin–orbit diagram of WD LPTs and WD pulsars.** Plot P_spin against P_orb for AR Sco,
J1912−4410, J0715 and the period = orbit LPTs. Is J0715 a "WD pulsar with a G star" rather than
a new LPT subclass? This is mostly catalogue work, folded into the v4 schema change (step 3).
Probably a figure in the `lpt` v4 paper, not a slice.

**F. Lower priority:** `lptv`-style Stokes-V limits at J0715. RM against orbital phase in archival
ASKAP polarization (`lptlin`, plan 92), which is probably too faint. An `lptxray` point (plan 93),
but read N. Rea's X-ray section first.

## 3. Candidate controls (to freeze before any real run)

- B: phase prediction, window and null stated above; the phase-error void condition is computed
  before opening the data.
- A: blank-sky null (forced flux negative about half the time), injection–recovery per target,
  recover-a-known (AR Sco, as in `wdpulsar`; J0715 itself if it is in RACS). The limit is divided
  by Σε with the activity-window factor.
- All: the paper's ephemeris is an *input*. Nothing fitted from it may certify it.
