# Findings — Milky Way HI rotation curve (tangent-point method)

`jansky_research.hi.run` builds the inner-Galaxy rotation curve by the tangent-point method from the
Leiden/Argentine/Bonn (LAB) HI 21 cm survey (Kalberla et al. 2005), fetching the $(b, v)$ slices at
$\ell = 10°$–$80°$ from VizieR (VIII/76), reading the $b=0$ spectrum, extracting the terminal
velocity, and applying $R = R_0\sin\ell$, $V(R) = v_\mathrm{term} + V_0\sin\ell$. We use the modern
Reid et al. (2019) Galactic constants $R_0 = 8.15$ kpc and $V_0 = 236$ km/s (their best-fit circular
speed at the Sun, not the older IAU 220 km/s). LAB velocities are already in the LSR frame, so no
solar-motion correction is applied.

## Result: an approximately flat (non-Keplerian) rotation curve

| $\ell$ (°) | $R$ (kpc) | $V$ (km/s) | note |
|----|------|------|------|
| 10 | 1.42 | 228 | bar region (R<4): unreliable |
| 20 | 2.79 | 220 | bar region |
| 30 | 4.07 | 247 | |
| 40 | 5.24 | 254 | |
| 50 | 6.24 | 264 | |
| 60 | 7.06 | 260 | |
| 70 | 7.66 | 252 | |
| 80 | 8.03 | 263 | |

Beyond the bar-dominated inner region, from $R \approx 4$ to $8$ kpc the rotation speed stays at
**$V \approx 257 \pm 6$ km/s** — it does **not** fall off as the Keplerian $V \propto R^{-1/2}$ that
the visible mass alone would give. (There is a gentle rise from $R\sim3$ to $\sim6$ kpc, likely
non-circular streaming, on top of a broadly flat curve.) Recovering this **approximately flat,
non-Keplerian** curve from public HI data with a ~30-line tool is the **validation** — a known,
foundational result, not a new measurement (`paper/figures/rotation_curve.pdf`).

## Honest limitations

- **Absolute level is ~9% high** ($257$ vs $V_0 = 236$ km/s). The main cause is the
  **terminal-velocity estimator**: a fixed 2 K brightness threshold *overestimates* $v_\mathrm{term}$
  relative to inflection-point spectral fitting (McClure-Griffiths & Dickey 2016 find threshold
  crossings $\sim7$ km/s higher), which inflates $V$; the rest is the simple estimator and
  non-circular motions. The **shape (flat, non-Keplerian) is robust; the absolute level is the soft
  part.**
- **Circular-orbit assumption.** $V = v_\mathrm{term} + V_0\sin\ell$ assumes the tangent-point gas is
  on a circular orbit; spiral-arm / bar-driven streaming biases individual points (not just adds
  scatter) and cannot be removed without a model.
- **Bar region excluded.** Points at $R \lesssim 4$ kpc ($\ell = 10°, 20°$) are dominated by the
  Galactic bar's non-circular motions and are flagged unreliable and dropped from the mean.
- **Not, by itself, a dark-matter detection.** An inner-Galaxy flat curve is *consistent* with a dark
  halo but can also arise from the disk's own mass; the unambiguous dark-matter case needs the curve
  **beyond $R \approx R_0$**, where the visible-mass contribution falls — a regime requiring a
  different method and not covered here.
- **First quadrant, $b=0$, inner Galaxy only.** One quadrant and the mid-plane spectrum; the
  tangent-point method works only for $R < R_0$.

## Bottom line

A clean, honest **validation**: a small, tested, offline-capable tool recovers the approximately flat
inner Milky Way rotation curve from public LAB HI data, with the absolute normalisation ~9% high (a
known terminal-velocity-estimator overestimate) and the non-Keplerian shape robust. A tooling +
validation contribution — and explicitly *not* a standalone dark-matter detection.
