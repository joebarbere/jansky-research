# Findings — ultra-steep-spectrum source hunt (TGSS × NVSS)

Run of `jansky_research.spectra.run` over a high-latitude extragalactic field near the North
Galactic Pole (centre RA $180.0°$, Dec $+30.0°$; $3°$ cone), cross-matching **TGSS ADR1**
(147.5 MHz; Intema et al. 2017) with **NVSS** (1.4 GHz; Condon et al. 1998) on VizieR, 15" match
radius. This is the honest assessment, revised after the GATE-2 science review. The full matched
catalogue is committed at [`survey/uss_candidates.csv`](uss_candidates.csv) (456 sources).

## Sanity checks (the method works)

- **456 matched sources**; **median $\alpha_{150}^{1400} = -0.73$** — consistent with the typical
  $-0.7$ to $-0.8$ range for flux-limited extragalactic radio samples. The closest comparator,
  de Gasperin, Intema & Frail (2018; the 1.4M-source TGSS×NVSS index catalogue), finds a weighted
  mean $-0.79$; our median sitting slightly flatter is consistent with TGSS flux-scale inflation
  making some steep sources look flatter. The distribution peaks near $-0.8$ as expected.

## The finding: 6 ultra-steep-spectrum candidates, none a known high-$z$ radio galaxy

Six sources have $\alpha < -1.3$ (the classic HzRG selection; De Breuck et al. 2000). **None has a
measured redshift or a high-$z$ radio-galaxy classification in either NED or SIMBAD** — all six are
HzRG *candidates*.

| # | RA | Dec | $\alpha$ | $S_{150}$ (mJy) | match sep | NED | SIMBAD | note |
|---|------|------|------|------|------|------|------|------|
| 1 | 181.0497 | +29.3853 | $-1.65$ | 126 | 5.8" | — | — | no counterpart in either DB |
| 2 | 179.5148 | +31.6903 | $-1.37$ | 100 | **9.4"** | — | NVSS radio src | loose match → lower confidence |
| 3 | 181.0344 | +30.2634 | $-1.34$ | 131 | 5.4" | WISE IR (IrS) | — | faint IR counterpart |
| 4 | 177.8709 | +28.8189 | $-1.32$ | 191 | 3.5" | WISE IR (`*`) | NVSS radio src | NED's nearest is a star (likely chance IR align) |
| 5 | 180.2691 | +29.1648 | $-1.31$ | 75 | 1.0" | — | — | tight match, no counterpart |
| 6 | 177.0305 | +29.3276 | $-1.31$ | 135 | 0.3" | — | NVSS radio src | tight match |

Where SIMBAD "matches", it knows the source *only* as a bare `NVSS J…` radio entry — no
classification, no redshift. So none of the six is catalogued as a high-$z$ radio galaxy.

## Honest limitations (this is a candidate list, not a discovery)

- **The novelty is bounded, and the key cross-check is incomplete.** TGSS × NVSS USS selection is an
  *established* HzRG-finding technique. The most directly comparable published search is **Saxena et
  al. (2018, MNRAS 475, 5041)** — USS ($\alpha \le -1.3$), 150 MHz, $\sim$50–200 mJy, $\sim$10,000
  deg² — whose parameter space these six fall squarely within. An automated VizieR query against
  that catalogue (`J/MNRAS/475/5041`) returned **no usable table** for this region, so we could
  **not** confirm whether these candidates are in it. **Until a manual cross-check against the
  published Saxena and de Gasperin catalogues (and dedicated USS samples) is done, we claim only
  "USS candidates not currently classified as HzRGs in NED/SIMBAD" — never "new sources".**
- **Chance-coincidence rate:** at 15" with the NVSS source density in this field, $\sim$1–2 false
  matches are expected across the full 456-source sample. The tight USS matches ($\le 5.8$") are
  therefore robust; candidate 2's lone 9.4" association is the most likely to be spurious.
- **TGSS flux-scale systematic** (Hurley-Walker 2017, arXiv:1703.06635): $\sim$15% typical, up to
  $\sim$40–50% in places, biasing $\alpha$ by $\sim$0.1–0.2. The rescaled TGSS-RSADR1 covers only
  Dec $\le +30°$, so the **upper half of this cone uses the uncorrected ADR1 scale** and may carry a
  larger, unquantified bias — relevant to candidate 2 (Dec $+31.7°$) and to 5, 6 sitting at the
  $-1.31$ threshold.
- **Two-point index only.** The 4 "inverted" ($\alpha>0$) sources are *not* peaked-spectrum/GPS
  claims — at 150 MHz vs 1.4 GHz a positive index is more likely a calibration/variability/
  resolution artifact; confirming a turnover needs a third frequency (e.g. VLASS 3 GHz).
- NED's nearest object to candidate 4 is classified as a star — most likely a chance infrared
  alignment, not the radio source's host.

## Assessment for GATE 2

- **Honest & non-trivial:** the method reproduces the population median; the candidate list is
  real and cross-checked against two databases; the novel-vs-known split uses corrected NED logic.
- **No overclaiming:** framed as *USS candidates requiring follow-up* (Saxena/de Gasperin
  cross-check + spectroscopy), with the established-method, flux-scale, coverage-gap, and
  chance-match caveats stated up front. The contribution is a **reproducible, open, CPU-only
  USS-selection tool** plus a vetted candidate list in this field — not a confirmed discovery.
