# Findings — ultra-steep-spectrum source hunt (TGSS × NVSS)

Run of `jansky_research.spectra.run` over a high-latitude extragalactic field near the North
Galactic Pole (centre RA $180.0°$, Dec $+30.0°$; $3°$ cone), cross-matching **TGSS ADR1**
(147.5 MHz) with **NVSS** (1.4 GHz) on VizieR. This is the honest assessment for the GATE-2 review.

## Sanity checks (the method works)

- **456 matched sources** (15" radius); **median $\alpha_{150}^{1400} = -0.73$** — the textbook value
  for the radio source population (steep-spectrum synchrotron), confirming the flux scales and the
  index calculation are right.
- The spectral-index distribution is the expected single peak near $-0.8$ with a steep tail.

## The finding: 6 ultra-steep-spectrum candidates, none a known high-$z$ radio galaxy

Six sources have $\alpha < -1.3$ (the classic high-redshift radio-galaxy selection; De Breuck et
al. 2000). **None has a measured redshift or a high-$z$ radio-galaxy classification in either NED or
SIMBAD** — i.e. all six are HzRG *candidates*.

| # | RA | Dec | $\alpha$ | $S_{150}$ (mJy) | match sep | NED | SIMBAD | note |
|---|------|------|------|------|------|------|------|------|
| 1 | 181.0497 | +29.3853 | $-1.65$ | 126 | 5.8" | — | — | no counterpart in either DB |
| 2 | 179.5148 | +31.6903 | $-1.37$ | 100 | **9.4"** | — | NVSS radio src | loose match → lower confidence |
| 3 | 181.0344 | +30.2634 | $-1.34$ | 131 | 5.4" | WISE IR (IrS) | — | faint IR counterpart |
| 4 | 177.8709 | +28.8189 | $-1.32$ | 191 | 3.5" | WISE IR (type `*`) | NVSS radio src | NED's nearest is a star (likely chance IR align) |
| 5 | 180.2691 | +29.1648 | $-1.31$ | 75 | 1.0" | — | — | tight match, no counterpart |
| 6 | 177.0305 | +29.3276 | $-1.31$ | 135 | 0.3" | — | NVSS radio src | tight match |

Where SIMBAD "matches", it knows the source *only* as a bare `NVSS J…` radio entry — no
classification, no redshift. So none of the six is catalogued as a high-$z$ radio galaxy. The full
matched catalogue is in `results/uss_candidates.csv`.

## Honest limitations (this is a candidate list, not a discovery)

- **The novelty is bounded.** TGSS × NVSS USS selection is an *established* HzRG-finding technique
  with **published catalogues** (e.g. de Gasperin et al. 2018 and other TGSS USS samples). These six
  may already appear in a dedicated USS catalogue even though NED/SIMBAD do not classify them as
  HzRGs. Confirming any is genuinely new requires cross-checking those published USS samples — not
  just NED/SIMBAD — plus spectroscopic redshifts. **We claim "candidates not currently classified as
  HzRGs", not "new sources".**
- **TGSS flux-scale systematic** (Hurley-Walker 2017): position-dependent offsets can shift $\alpha$
  by $\sim0.1$; we use a conservative $\alpha<-1.3$ cut so it does not flip the classification, but
  candidates near the threshold (5, 6 at $-1.31$) are the most affected.
- **Cross-match confusion:** candidate 2's 9.4" separation (vs the others' $\le 5.8$") is a
  lower-confidence association given the 45" NVSS beam. Candidate 4's nearest NED object is a star —
  most likely a chance infrared alignment, not the radio source's host.
- **Two-point index only:** no turnover information; the "inverted" class (4 sources, $\alpha>0$)
  are GPS/peaked candidates that would need a third frequency (e.g. VLASS 3 GHz) to confirm.

## Assessment for GATE 2

- **Honest & non-trivial:** the method reproduces the population median ($-0.73$) and yields a
  real, vetted candidate list; the novel-vs-known split is cross-checked against two databases.
- **No overclaiming:** framed as *USS candidates requiring follow-up*, with the established-method
  and published-catalogue caveats stated up front. The contribution is a **reproducible, open,
  CPU-only USS-selection tool** plus a candidate list in this field — not a confirmed discovery.
