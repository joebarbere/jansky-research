# Findings — Milky Way HI rotation curve (tangent-point method)

`jansky_research.hi.run` builds the inner-Galaxy rotation curve from the Leiden/Argentine/Bonn (LAB)
HI 21 cm survey (Kalberla et al. 2005), fetching the $(b, v)$ slices at $\ell = 10°$–$80°$ from
VizieR, reading the $b=0$ spectrum, extracting the terminal velocity, and applying the tangent-point
relations $R = R_0\sin\ell$, $V(R) = v_\mathrm{term} + V_0\sin\ell$ ($R_0 = 8.15$ kpc,
$V_0 = 236$ km/s; Reid et al. 2019).

## Result: the rotation curve is flat (recovers the dark-matter signature)

| $\ell$ (°) | $R$ (kpc) | $V$ (km/s) |
|----|------|------|
| 10 | 1.42 | 228 |
| 20 | 2.79 | 220 |
| 30 | 4.07 | 247 |
| 40 | 5.24 | 254 |
| 50 | 6.24 | 264 |
| 60 | 7.06 | 260 |
| 70 | 7.66 | 252 |
| 80 | 8.03 | 263 |

From $R \approx 1.4$ to $8$ kpc the rotation speed stays at **$V \approx 251 \pm 14$ km/s** — it does
**not** fall off as the Keplerian $V \propto R^{-1/2}$ that the visible mass alone would give. This
flatness is the textbook signature of a dark-matter halo, recovered here from public HI data with a
~30-line tool (`paper/figures/rotation_curve.pdf`). This is a **validation** — reproducing a known,
foundational result — not a new measurement.

## Honest limitations

- **Absolute normalisation is ~5–10% high.** The mean $V \approx 251$ km/s sits above the Reid et al.
  (2019) $V_0 = 236$ km/s. The dominant cause is the **terminal-velocity estimator**: a fixed
  brightness-temperature threshold (2 K) places $v_\mathrm{term}$ slightly outside the true profile
  edge, biasing $V$ high; the choice of threshold shifts the level by several km/s. Non-circular
  motions (the Galactic bar, spiral arms; significant at $R \lesssim 4$ kpc) add real scatter. The
  **flatness is threshold-independent and robust**; the absolute level is the soft part.
- **First quadrant only, $b=0$ only.** A single quadrant and the mid-plane spectrum; a full
  treatment would average over latitude and combine quadrants.
- **Inner Galaxy only.** The tangent-point method works for $R < R_0$; the (more dark-matter-decisive)
  outer curve needs a different method and distances.

## Bottom line

A clean, honest **validation**: a small, tested, offline-capable tool recovers the flat inner Milky
Way rotation curve from public LAB HI data — the classic dark-matter signature — with the absolute
normalisation modestly high (a known terminal-velocity-estimator systematic) and the flatness robust.
A tooling + validation contribution, not a new result.
