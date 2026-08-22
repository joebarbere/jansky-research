# lptxray — findings

Plan 93. GATE 0 run 2026-08-22. **The LPT leg of the planned comparison is cut**, and the
reason it is cut is the more useful result: an archival X-ray catalogue cross-match is
structurally blind to the observations that produced every published LPT X-ray detection.

## GATE 0

### 1. Novelty holds

No class-wide, chance-corrected X-ray cross-match exists for either list. Per-source X-ray
follow-up is well established and expected — the LPT review (arXiv:2601.10393) tabulates it
and confirms detections are limited to ASKAP J1832−0911 and ASKAP J1448−6856 — but a uniform
re-search with a *measured* false-match rate has not been published.

Pelisoli et al. 2025's own `xray` column (vendored in `data/wdpulsar_candidates.csv`; 17 of 56
non-blank: 11 XMM-Newton, 5 ROSAT, 1 INTEGRAL) is explicitly **not** a cross-match. The paper
says those columns were *"not used in the candidate selection, but defined follow-up
priority"* — a compilation of prior literature associations, with no false-match statistic.

### 2. eROSITA-DE: public, and half the sky

DR1 (Merloni et al. 2024, A&A 682, A34, arXiv:2401.17274) and the deeper **DR2** (released
2026-07-31, ~2M sources, eRASS:3 stack) are both public and both cover **180 ≤ l ≤ 360 only**.
The eastern half is Russian-consortium and is not public. Measured split:

| sample | in eROSITA-DE (public) | in the Russian half |
|---|---|---|
| 16 LPTs | **9** | 7 |
| 56 candidates | **38** | 18 |

`dr20radio` deferred eROSITA over exactly this footprint restriction; the restriction is
unchanged, so the leg is used here only where it covers.

### 3. The cross-match has no power on LPTs — measured, not assumed

Three LPTs carry a published X-ray detection: ASKAP J1832−0911, ASKAP J144834−685644 and
ASKAP J174508.9−505149. Queried against every serendipitous X-ray source catalogue
(5XMM-DR15 `xmmssc`, 2RXS `rass2rxs`, eRASS1 `erass1main`, CSC 2.1.1 `csc`), **only one is
recovered** — the Rose et al. source, in eRASS1.

The other two are not faint. They are missing because the detections come from **dedicated
pointed observations taken after the catalogues were built**, which the HEASARC observation
logs confirm:

| source | pointed coverage | taken |
|---|---|---|
| ASKAP J144834−685644 | XMM **0953011101**, target "VAST J1448-6856" | MJD 60522 (2024-08) |
| ASKAP J174508.9−505149 | XMM **0973390301**, target "ASKAP J1745" | MJD 60952 (2025-10) |
| ASKAP J1832−0911 | Chandra **26681 / 26682 / 29265 / 29266** | MJD 60354–60533 (2024) |

So **catalogue recall on the LPT sample is 1/3**, and the failure mode is catalogue latency
plus pointed-observation processing, not source faintness. A catalogue-derived LPT X-ray
detection fraction would measure the catalogues, not the sources.

**The same machinery recovers 16 of 16 on the candidates.** Every one of Pelisoli's 11
XMM-Newton and 5 ROSAT identifications is re-found independently (15 of the 16 also appear in
2RXS). Catalogue cross-matching has power on that sample and none on the LPTs, so the two legs
are reported separately and **never differenced**.

This is the `frblens` lesson in a new place: a null divides by sensitivity, not by sample
size. `lptxray.catalogue_recall` is the guard that enforces it — it marks a leg unusable when
recall on known detections falls below 0.8.

### Consequent scope change

The plan's headline — "does accretion predict radio loudness?" — is also weaker than it was
written. Only **one** of the 16 LPTs (the Rose et al. source) is a confirmed accretor; one is
a confirmed non-accretor, two more are detached non-accreting binaries, and the remaining ~12
have unknown accretion state. That is a one-versus-fifteen contrast, not a population split.

What the slice can honestly deliver:

1. the accretion split **within** the candidate list, where the catalogues demonstrably work,
   set against `wdpulsar`'s radio limits for the same objects;
2. a **pointed-coverage census** for the 16 LPTs — which have archival X-ray data, and which
   of those have had none of it published.

## A confound that is measured rather than asserted

The two lists sit in different parts of the sky, so soft-X-ray absorption differs by orders of
magnitude between them:

| sample | median \|b\| | fraction within 5° of the plane |
|---|---|---|
| LPTs | **1.88°** | **0.62** |
| candidates | 17.03° | 0.09 |

The LPTs are mostly in-plane and at kpc distances; the candidates are nearby and
high-latitude. A soft-band non-detection means something very different in the two samples,
which is a second, independent reason not to difference them.

## The circularity check

The accretion split within the candidate list is partly circular by construction: Pelisoli's
X-ray column set spectroscopic follow-up priority, and that spectroscopy produced the `type`
classifications. So "accreting systems are X-ray bright" could in part be an artefact of which
candidates got classified at all.

`accretion_comparison.xray_no_prior_flag` re-runs the comparison on the candidates carrying
**no** literature X-ray flag, for which that path is cut. The check is free to fail, and it
should be read before the headline number.

## Method

`scripts/lptxray_fetch.py` stages one 10-arcmin HEASARC cone per (position, catalogue) into
`data/lptxray/cones.json`. Everything downstream is offline, which matters for the chance
rate: the rigid position-shift trials (3/5/7 arcmin, 8 azimuths, every trial disc required to
lie wholly inside the cached cone) cost no extra queries, and the field source density is
measured in a 2–10 arcmin annulus around each target. Association radii are per catalogue —
15″ (5XMM), 45″ (2RXS), 20″ (eRASS1), 10″ (CSC) — set by each catalogue's positional accuracy.

Fluxes are recorded per catalogue and **never combined**: the four are different bands and
2RXS is a count rate, so a single "X-ray flux" would need a spectral assumption this sample
cannot support.

## Measurement

`results/lptxray_metrics.json`. 72 positions × 4 catalogues.

### The chance-coincidence rate is ~1%, measured two independent ways

| catalogue | radius | expected from field density | rigid position-shift trials |
|---|---|---|---|
| 5XMM-DR15 | 15″ | 0.0104 | 0.0098 (17/1728) |
| eRASS1 | 20″ | 0.0030 | 0.0029 (5/1728) |
| 2RXS | 45″ | 0.0024 | 0.0012 (2/1728) |
| CSC 2.1.1 | 10″ | 0.0014 | 0.0012 (2/1728) |

The two estimates are independent — one is an area × density calculation, the other counts
what displaced positions actually find — and they agree to better than a factor of two. Chance
coincidence cannot produce the split below.

### Accretion predicts X-ray brightness, strongly

Among the 56 candidates, split by Pelisoli's own classification (accreting = polar 16 + IP 2 +
CV 3 = 21; other = YSO 26 + unclear 8 + pulsar 1 = 35):

| group | X-ray detected | fraction (95% Wilson) |
|---|---|---|
| **accreting** | **20/21** | **95.2% [77.3, 99.2]** |
| other | 3/35 | 8.6% [3.0, 22.4] |

Ratio **11.1**, difference 86.7 pp. Both are quoted because a difference in percentage points
inherits the normalisation and a ratio does not.

The single accreting non-detection, J0452+3017 (IP), lies at l = 172.3° — **in the Russian
half, with no public eROSITA coverage**. That is a sensitivity explanation, not a physical one.

### Three cuts that could each have broken it, and did not

| leg | accreting | other | ratio |
|---|---|---|---|
| all 56 | 20/21 | 3/35 | 11.1 |
| eROSITA-DE half only (common footprint) | 14/14 | 3/24 | 8.0 |
| **2RXS only (all-sky, no footprint asymmetry)** | **18/21** | **0/35** | — |
| 5XMM only | 12/21 | 1/35 | 20.0 |
| CSC only | 0/21 | 0/35 | no power |

The footprints are balanced to begin with (67% of the accreting and 69% of the other group sit
in the eROSITA-DE half), and the 2RXS leg settles the question: a single all-sky catalogue,
identical exposure treatment for both groups, 18/21 against 0/35. CSC has no power on this
sample at all — Chandra's sky coverage is a rounding error — which is worth stating rather than
reporting as a null.

The accreting detections are corroborated across three independent catalogues (18 in 2RXS, 14
in eRASS1, 12 in 5XMM); the three "other" detections are almost all single-catalogue eRASS1.

### The circularity check: same direction, but thin

Restricted to candidates carrying **no** literature X-ray flag — the ones never prioritised for
spectroscopic follow-up on X-ray grounds — the split is 3/4 accreting against 3/35 other, ratio
**8.75**. Same direction and similar size, so the effect is not purely an artefact of how the
types were assigned.

**But the accreting arm here is four objects.** Its interval is [30.1, 95.4]% and only just
clears the other group's upper bound of 22.4%. This check constrains the circularity worry; it
does not eliminate it, and it should not be quoted as if it did.

### Three X-ray associations not in Pelisoli's compilation

| candidate | type | catalogues | separation |
|---|---|---|---|
| **J1912−4410** | pulsar | eRASS1 + 5XMM | 1.4″, 0.9″ |
| J0530+0148 | YSO | eRASS1 | 3.3″ |
| J1226−2304 | unclear | eRASS1 | 4.4″ |

J1912−4410 is the *confirmed* white-dwarf pulsar of the sample, and Pelisoli et al. report it
shows no persistent accretion. Its X-ray counterpart is positionally solid (two catalogues,
~1″). These are associations, not published identifications — eRASS1 postdates the compilation
that the `xray` column was built from.

### Accretion does not predict radio loudness

Set against `wdpulsar`'s forced photometry for the same objects:

| group | radio-detected | with a usable limit | median 5σ limit |
|---|---|---|---|
| accreting | **0** | 19/21 | 0.818 mJy |
| other | **0** | 30/35 | 0.825 mJy |

The two depths are matched to under 1%, so this is a genuine common-limit comparison rather
than two different searches. Nineteen systems that the X-ray data independently confirm are
accreting show **no** persistent radio emission at ~0.8 mJy.

**The load-bearing caveat: this constrains *persistent* emission only.** LPTs are duty-cycled
pulsators, and `wdpulsar` already established the point on this very sample — J1912−4410's
MeerKAT pulses dilute to ~0.1–0.2 mJy time-averaged and it is absent from RACS despite being a
confirmed radio pulsar. So these limits do **not** exclude LPT-like pulsed emission in either
group, and the comparison cannot answer plan 93's question as posed. What it does establish is
that accretion, which is plainly visible in X-rays for these systems, is not accompanied by
detectable persistent radio emission.

Taken with the LPT side — the one confirmed accreting LPT is radio-loud, while 19 confirmed
accretors here are not persistently radio-loud — accretion looks neither necessary nor
sufficient for radio loudness. The numbers on the LPT side are far too small to call that a
result.

### Pointed X-ray coverage of the 16 LPTs

Classified by aimpoint offset rather than by target name, because names mislead exactly where
it matters: all 24 pointings near ASKAP J1935+2148 are observations of **SGR 1935+2154**, a
different source a few arcmin away that shares the same RA digits. Targeted pointings land
within 0.12′ of the LPT; the classification is clean with no ambiguous cases.

- **12 of 16** have some archival X-ray coverage;
- **8 of 16** have a *dedicated* pointing (GLEAM-X J1627 ×2, CHIME J0630+25 ×2,
  GCRT J1745−3009, GPM J1839−10, GLEAM-X J0704−37, ASKAP J1832−0911, ASKAP J1448−6856,
  ASKAP J174508.9−505149);
- **5 of those 8** have no published X-ray *detection*.

**"No published detection" does not mean "unanalysed".** At least GLEAM-X J1627 and
GPM J1839−10 have published X-ray upper limits derived from these very observations. The
census's value is that it is uniform and archive-derived, not that it has found neglected data.

Two LPTs whose coverage is purely serendipitous — ASKAP J1935+2148 (24 pointings, all of
SGR 1935+2154) and ILT J1101+5521 (32, an NGC 3079 field) — sit in deep archival data that was
taken for other reasons. That is the one genuinely actionable item here.

## Status

A recorded result, not a paper. The novel, defensible content is (a) the accretion/X-ray split
with a measured chance rate and three robustness legs, (b) the demonstration that archival
catalogue cross-matching has no power on the LPT class and why, and (c) the pointed-coverage
census. The plan's headline question is **not answered**, and the reason is duty cycle rather
than anything about accretion.
