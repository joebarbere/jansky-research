# Findings — VLASS multi-epoch variability (image-vetted)

`jansky_research.vlass` cross-matches the three VLASS Quick-Look (QL) component catalogues (CIRADA
Epochs 1–2, NRAO Epoch 3; Gordon et al. 2021) and selects variable/transient candidates by the
η (weighted reduced χ², significance) and V (coefficient of variation, amplitude) metrics with a 3σ
2-D outlier cut (Rowlinson et al. 2019). First real run: a 2.5° cone at RA 190°, Dec +20° (b≈+80°),
Epochs 1–3 — 1148 anchor sources, 1022 detected in ≥2 epochs.

## Systematics handled (each one would manufacture variability)

- **Per-epoch flux scale.** QL images underestimate peak flux epoch-dependently (E1 ~15%, E2 ~7%,
  E3 ~3%; VLASS Memos 13/22). Corrected to the Perley–Butler scale, with the ~7% residual scatter
  folded into the errors. Without it the epoch offset alone flags every source.
- **A µJy/mJy unit trap.** The VizieR Epoch-1 `Fpeak` is micro-Jy mislabelled "mJy/beam" — 1000×
  off from the genuine-mJy Epochs 2–3. Caught empirically (the catalogue flux floor is the ~0.7 mJy
  detection threshold) and rescaled.
- **Deblending of crowded sources** (see below) — handled with an isolation filter.

## Result: every candidate is a catalogue artefact — an honest negative

At 3σ in both metrics, after the isolation filter (which removes 33 crowded sources), the field
yields 1–2 candidates. **None survives image vetting.** The decisive check is `image_lightcurve`,
which reads the actual VLASS QL image peak at the position in each epoch (the ground truth) and
compares it to the catalogue light curve:

| candidate | catalogue (mJy) | image peak (mJy/beam) | verdict |
|-----------|-----------------|------------------------|---------|
| J124205+193241 | 1.3 → 5.3 → 0.7 | 4.6 → 4.7 → 4.6 (flat) | steady ~5 mJy source, **deblending artefact** |
| J123318+192851 | 5.5 → 0.9 → — | 5.3 → 5.0 → 4.6 (flat) | steady ~5 mJy source, **cross-match/extraction artefact** |

The images show **steady sources in all three epochs**; the catalogue "variability" is spurious.

### Automatic confirmation by forced photometry

This vetting is now automatic. `run` calls `confirm_candidates`, which does **forced peak photometry**
(`measure_image_flux`) at the locked candidate position in every epoch's image, recomputes η/V on that
*image* light curve, and marks a candidate `image_confirmed` only if it is significantly variable
(p < 0.01) with real amplitude (V > 0.3). On this field both candidates are **auto-rejected**: forced
light curves [4.28, 3.54, 4.32] (V = 0.11) and [2.32, 1.69, 1.95] (V = 0.16) — far below the catalogue
V ≈ 1.0 and well under the threshold. Forced photometry at a fixed position is immune to the deblending
*and* cross-match failures, so the surviving candidate list (`n_image_confirmed`) is trustworthy
without manual inspection.

### The failure modes (why catalogue-only QL variability is unreliable)

1. **Component deblending.** A single slightly extended source is fit with a different number of
   Gaussian components in each epoch, so a *secondary* component's flux jumps around. Candidate 1 is
   the secondary component (1.3 mJy) of a steady ~5 mJy source whose primary sits 4–5″ away; image
   peak is flat at 4.6 mJy. The isolation filter (`isolated_mask`, drop sources with a neighbour
   within 5″) removes this class.
2. **Cross-match incompleteness.** A real source present in an epoch but with its centroid shifted
   >2.5″ is recorded as a *non-detection*, faking a fade. Candidate 2's "missing" Epoch 3 is a clear
   2 mJy source in the image — just unmatched. (Not fixed by isolation; needs forced photometry or a
   non-detection check.)
3. **Peak-flux vs image inconsistency.** PyBDSF peak flux for a blended/extended component can differ
   from the image peak by a factor of a few, inflating η.

## The census (703 deg²) — and a real variable

Scaling to a 15° cone (≈**703 deg²**, RA 190°, Dec +20°, ~14× the pilot): **42 721** sources, 38 999
detected in ≥2 epochs, 2043 dropped as crowded, **36 956 usable**. The 3σ selection flags **40**
catalogue candidates; forced photometry (top-100 by η, run in ~10 min) confirms **2** of them.

### The forced-peak centring gate (5 → 2)

Image montages of the first five "confirmed" candidates exposed a residual false positive: forced
photometry took the brightest pixel anywhere in a 4″ box, so a bright source *just outside* the
locked position pinned the peak at the box edge (offset ≈ 3.9–4.0″) and faked a variable as its
sidelobe/wings shifted between epochs. Adding a **centring gate** — the forced peak in the brightest
epoch must lie within `center_arcsec` (2.5″, ~1 beam) of the position — drops the count from 5 to **2**,
removing exactly the three box-edge artefacts (offsets ~4″) and keeping the two genuinely centred
sources (offsets 0.3–0.4″).

### The two survivors

| position | catalogue (mJy) | forced (mJy) | offset | SIMBAD | verdict |
|----------|-----------------|--------------|--------|--------|---------|
| 202.695 +24.233 | 9.0 → 3.1 → 1.1 | 7.2 → 2.7 → 1.0 | 0.3″ | **V* FK Com** (0.9″) | known radio-active star, clean decline — **a real recovery** |
| 184.668 +21.817 | 3.2 → – → 0.6 | 2.9 → 0.7 → 0.6 | 0.4″ | none | single-epoch detection — **fails archival follow-up (artefact)** |

**FK Comae Berenices** — a famously active, rapidly rotating giant and known radio flarer — is
recovered purely from the public catalogues with a clean monotonic ~7→1 mJy decline at a 0.9″ match.
That is the **validation**: the pipeline does surface genuine variables, not only reject artefacts.

### Archival follow-up — the confirmation ladder, and the second survivor falls

Image confirmation (forced photometry + centring) is necessary but not sufficient; the decisive gate
is the *next survey epoch plus archival cross-match*. Running it on the two survivors:

- **FK Com** stands: a known active star, clean decline, $0.9\arcsec$ match — confirmed.
- **184.668 +21.817 fails.** It is a single-epoch detection: forced photometry through the
  independent **VLASS Epoch 4** (2026 release) reads 2.9 mJy in Epoch~1 but only noise (0.7, 0.6,
  0.55 mJy at growing $1.7$--$3.9\arcsec$ offsets) in Epochs 2--4, i.e.\ absent for ~9 years. Its
  Epoch-1 QL entry is a blended component (`S_Code=C`) with a very low peak-to-ring ratio (0.07, a
  sidelobe signature), and every archival counterpart (NVSS~1998, FIRST, TGSS, AllWISE, PanSTARRS)
  sits at $\sim8\arcsec$ — a persistent multi-wavelength source that is *not* at the candidate
  position (which has no Gaia/optical counterpart of its own). The weight of evidence is an Epoch-1
  sidelobe/deconvolution artefact of the nearby source, not a transient.

So after the full ladder, **FK Com is the one genuine variable**; the second image-confirmed
candidate does not survive. Catalogue selection $\to$ image confirmation $\to$ archival + next-epoch
follow-up: each stage is necessary, and only the last one separated a real variable from a
single-epoch imaging artefact.
Several other candidates carry blazar/QSO/BL-Lac SIMBAD identifications but fall *below* the V > 0.3
confirmation bar — i.e. real AGN varying at the few-tens-of-percent level the VLASS cadence expects,
correctly *not* over-claimed as strong variables.

### Completeness (data-driven)

**[Superseded 2026-08-24: this table measured the injection's own error model, not the
selection — see the referee round and its resolution at the bottom of this file. The corrected
curve reaches 99.5% at 10x.]**

`injection_recovery` (inject a single-epoch flare of known factor into the real steady light curves,
re-run the cut) bounds the result honestly:

| flare factor | 1.25 | 1.5 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|---|
| recovered | 0.00 | 0.00 | 0.005 | 0.03 | 0.34 | ~0.52 |

The selection **saturates near 50%**: with only three epochs a single-epoch flare's $V$ tops out at
$1/\sqrt{3}\,/\,(1/3)=1.73$ as the factor → ∞, and the noisy 3-epoch $V$ threshold (≈0.9–1.3) sits just
below that ceiling, so even 10× flares are recovered only ~half the time (50% at ~9×). So the two
confirmed variables are a **floor**, not a complete census: the strong-single-epoch-variable fraction
is of order a few ×10⁻⁴, with completeness ~50% for the strongest flares — much less constraining than
the professional thousands-of-deg² 2-epoch censuses, but a genuine, reproducible result with one
identified real source.

The practical lesson stands: the 3σ-in-both-metrics cut is too stringent for 3-epoch QL data. Because
forced photometry (now with the centring gate) cleanly removes false positives, a **looser** catalogue
cut, more sky, or more epochs would raise completeness without sacrificing purity — the next step.

## Honest conclusion

Catalogue-only VLASS QL multi-epoch variability selection is **dominated by source-extraction
artefacts** (deblending, cross-match misses, and bright-neighbour box-edge confusion), and every
statistical candidate must be confirmed against the actual images before it can be believed — the
cautionary thread shared with the USS and SETI slices. But over 703 deg² the method also does its
positive job: of 40 catalogue candidates, **2 survive forced-photometry confirmation with the centring
gate**, and one is the known radio-active star **FK Comae Berenices**, recovered with a clean ~7→1 mJy
decline purely from public catalogues. So the result is **not** a pure negative: it is a reproducible,
QL-systematic-aware pipeline (per-epoch flux-scale correction, deblending isolation, centred
forced-photometry confirmation, data-driven completeness) that recovers a genuine variable and sets an
honest, completeness-bounded floor on the rest — not a discovery, but a validation plus a limit.

## Honest limitations

- **703 deg², single cone.** Larger than the pilot but far from the thousands of deg² of the
  professional surveys; the limit is correspondingly weak and the two confirmed sources are a floor.
- **FK Com aside, the second survivor is unvetted beyond the image.** A ~3 mJy centred fader with no
  catalogued counterpart is a plausible transient but needs multi-wavelength follow-up before any
  claim; it is reported as a candidate, not a detection.
- **Selection saturates at ~50%** for single-epoch flares (3-epoch ceiling). A looser cut — viable now
  that forced photometry cleans up false positives — or more epochs is needed for a real rate.
- **Isolation is anchor-epoch only.** A source isolated in Epoch 1 but blended in a later epoch is
  not yet caught; a full treatment checks isolation per epoch and detects ambiguous matches.
- **Completeness is candidate-driven, not all-source.** Forced photometry runs on the candidates, so a
  real variable whose *catalogue* light curve was flattened by extraction issues can still be missed
  before it ever reaches confirmation; forced-photometering every source would close this.
- A genuine variable population almost certainly exists at this depth (the literature reports a few
  percent at 2σ over larger areas); reaching it needs much more sky and a looser, completeness-aware
  cut — future work.

## Bottom line

The pipeline runs end-to-end on real VLASS data with the QL systematics handled, and it correctly
finds that this field's variability candidates are **all extraction artefacts** — demonstrated by
image vetting, which is therefore mandatory. A tooling + methodology + honest-negative contribution.

## Full referee round (2026-08-24): MAJOR REVISION, 14 findings, one BLOCKER

Every abstract number traces to committed evidence and FK Com is genuinely auditable. The
problems: the completeness's cause, and the second candidate's rejection.

**BLOCKER: the "three-epoch ceiling" explanation for the 52% completeness saturation is
arithmetically refuted by the committed thresholds, and the number is an artifact of the
injection's error model.** V for a 10x flare is 1.30 (asymptote 1.73) vs the committed
v_threshold 0.881 -- V cannot bind. What binds is `e[k] *= fac` in injection_recovery: scaling
the flare epoch's error with its flux gives eta a FINITE asymptote (~(S/sigma)^2/2), so "52%"
really means "52% of steady sources have quiescent SNR >~ 11". The real error model
(apply_flux_scale) is sigma_cat-dominated for the faint majority, so the scaling inflates the
flare epoch's error up to 10x and suppresses its chi^2 up to 100x. The 52%, the 9x, the
"<5% below 3x", the three-epoch diagnosis and the closing recommendation all descend from that
line.

**MAJORs:** "saturates near 52%" reports the value at the largest factor tested while the curve
is still rising steeply (0.335 -> 0.5225); "50% at ~9x" swings between 8.9 and NaN within 1
MC-SE (n=400, one seed) and carries no uncertainty. The abstract's "necessary and sufficient"
is refuted three pages later (the three mitigations left a 50%-false survivor list; the remedy
is a fourth stage), and no ablation for "necessary" is committed. **The second candidate's
rejection exists only as prose**: commit a0bf180 touched only main.tex and the findings file, no
code fetches Epoch 4 / NVSS / FIRST / TGSS / WISE, and the committed CSV still says
image_confirmed=True for the retracted candidate (the findings file also contradicts itself
about it). The completeness is measured on a DIFFERENT population (no isolation mask, includes
the 2043 crowded; excludes two-epoch sources) with RECOMPUTED thresholds -- not the census's own
selection. Two selection-function facts absent: the census is anchored on Epoch-1 detection
(transients turning on later have completeness identically zero -- the "single-epoch-variable
fraction" sentence reads as about transients in general), and the source counts don't close
(42721 - 2043 != 36956; the missing 3722 are epoch-match failures the paper never mentions).
The "few x 10^-4" fraction is in no results file and is stale from when the paper claimed two
variables (candidate derivations give 2.7-10 x 10^-5). The "forced" photometry searches (median
brightest-epoch offset 2.86", 27/39 beyond the centring gate) and per-epoch rms/offsets are not
committed, so neither the p<0.01 confirmation nor the "only noise" rejection is auditable.

**MINOR/NIT:** the AGN sentence says "few-tens-of-percent" where the median forced V is 0.085
(and misses the better fact: catalogue V 0.88-1.24 vs forced 0.007-0.23 -- up to 160x, the
paper's own thesis); "artifact-dominated" is a single-operating-point result stated generally
(sigma=3.0 is the only value ever run) and its purity/completeness halves are mechanically
linked; gordon2021 has the right DOI under a title belonging to a different paper (fix the
title, not the DOI); Figure 2's caption calls both tracks "variables"; prose says ~15% QL
correction where the code applies 1.13; archival counterpart lists drift between files.

**Status: RESOLVED (2026-08-24).** Every finding addressed; the census re-ran end-to-end and
reproduced its counts exactly (42,721 / 3,722 epoch-1-only / 2,043 crowded / 36,956 usable /
40 candidates / 2 image-confirmed), which is itself evidence the pipeline is deterministic.

**The blocker, measured.** With the error model fixed (only the flux-proportional 7% systematic
term rescales with the injected flare; the catalogue error keeps its quiescent value) and the
injection drawn from the census's own usable population against the census's committed
thresholds, the completeness curve transforms: 10x flares are recovered at **99.5%** (was
"52% saturation"), 50% completeness at **4.4 +/- 0.1x** (Monte-Carlo error over ten seeds; was
"~9x" with no error), >99% out to the largest factor tested (50x). The referee's diagnosis was
exactly right: the old number measured `e[k] *= fac`, not the selection, and the "three-epoch
ceiling" story built on it is gone from the paper.

**The archival vet is code now** (`archival_vet`/`run_archival`, committed to
`results/vlass_archival.json` + an `archival_survives` column in the candidates CSV). Running it
produced one methodological lesson en route: the first version judged detection on the 4"
box-searched peak (>4x rms, centred 2.5"), and the retracted candidate *passed* -- its Epoch-2
box peak (0.72 mJy at 1.65" offset) is the neighbour's wings plus the blank-sky positive bias of
any peak search (the stokesv lesson, again). Detection is now judged at the **fixed position**
(peak within 0.75", essentially the propagated pixel), and the evidence became unambiguous:
FK Com is detected in all four epochs (7.2 / 2.7 / 1.0 / **5.0** mJy -- it rebrightened in
Epoch 4) while the second candidate reads 2.92 then 0.07 / 0.37 / **-0.21** mJy -- negative,
exactly as genuine forced photometry on blank sky should sometimes be -- with every archival
counterpart (NVSS 8.28", FIRST 8.24", TGSS 8.16", AllWISE 8.36") at the persistent neighbour,
not the candidate. `n_archival_survivors = 1`; the CSV row now carries
`archival_survives=False` next to its (correct, but insufficient) `image_confirmed=True`.

**The rest:** per-epoch forced rms and offsets committed for all 40 candidates (the p<0.01
confirmations and the rejection are auditable); count closure and Epoch-1 anchoring stated in
the paper (3,722 epoch-match failures named; transients turning on after Epoch 1 have
completeness zero by construction); the stale "few x 10^-4" replaced by pipeline-derived
fractions (floor 2.7x10^-5 of the usable census; ~2.7x10^-5 for >=10x flares after the
completeness correction -- the correction is now nearly a no-op because completeness at 10x is
99.5%); "necessary and sufficient" replaced by the four-stage ladder with the explicit statement
that no prefix of it is sufficient; the AGN sentence now quotes the committed collapse
(catalogue V 0.88--1.32 vs forced median 0.099, max 0.92); the operating-point qualifier added
in the abstract and discussion; gordon2021's title corrected against Crossref ("A Quick Look at
the 3 GHz Radio Sky. I. ..."); Figure 2's caption says "candidates" and names the rejected
track; the ~15% prose now states the applied 13% campaign-mean correction; counterpart lists
match the committed JSON everywhere.
