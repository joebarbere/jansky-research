# Findings — the pulsar P–Ṗ diagram (ATNF)

`jansky_research.ppdot` turns the ATNF Pulsar Catalogue's periods and spin-down rates into the
P–Ṗ diagram — the H–R diagram of pulsars — deriving surface fields, ages, and the death line, and
classifying the population. Reuses `pulsarspec.is_millisecond` and the project's ATNF VizieR fetch.
Tooling + real-data + recover-a-known done together (VizieR is a single reliable query).

## Data

ATNF Pulsar Catalogue (Manchester et al. 2005) via VizieR `B/psr` (public, no auth): periods `P0` (s)
and period derivatives `P1` (Ṗ). Of the catalogue, **2 052** pulsars have a positive measured Ṗ and
enter the diagram.

## Recover-a-known: the population structure

**Surface field** $B = 3.2\times10^{19}\sqrt{P\dot P}$ G cleanly separates the three populations:

| class | N | median $\log_{10}(B/\mathrm{G})$ | interpretation |
|---|---|---|---|
| millisecond ($P<30$ ms) | 182 | **8.42** | recycled by accretion — low field, fast spin |
| normal | 1 787 | **12.05** | the bulk rotation-powered population |
| magnetar / high-$B$ ($B>10^{13}$ G) | 83 | **13.31** | the top-right of the diagram |

The $\sim$5-orders-of-magnitude ($\sim$10$^{4.9}$, a factor of $\sim$80 000) spread in median field
between millisecond and high-$B$ pulsars is the textbook P–Ṗ structure (Lorimer & Kramer 2004): MSPs
sit at the bottom-left (short $P$, tiny $\dot P$), ordinary pulsars at $B\sim10^{12}$ G, magnetars at
the top-right. We use the standard $B = 3.2\times10^{19}\sqrt{P\dot P}$ G (orthogonal rotator,
$\sin\alpha=1$, $I=10^{45}$ g cm², $R=10$ km); the aligned/equatorial convention gives a coefficient
of $\sim$6.4$\times10^{19}$, i.e. all our fields scale up by $\sqrt2$.

**Death line.** With the constant-$B/P^2$ polar-cap criterion ($B_{12}/P^2=0.2$; Ruderman & Sutherland
1975; Bhattacharya & van den Heuvel 1991), **98.2%** of the catalogue lies *above* the death line — i.e.
almost every catalogued radio pulsar is on the radio-loud side, as expected for a radio-selected sample.
This is **model-dependent** (the line spans a "death valley", Chen & Ruderman 1993), and the $\sim$1.8%
($\sim$37 pulsars) below the line are not noise — they are real long-period pulsars (e.g. PSR
J2144$-$3933) that pushed the death-line literature. The Crab pulsar validates the derivations:
$B = 3.8\times10^{12}$ G and $\tau = P/2\dot P \approx 1.26\times10^3$ yr, both textbook.

## Honest assessment & caveats

- **A reproduction, not a discovery.** The tool recovers the well-known P–Ṗ population structure from a
  public catalogue; the contribution is a tested, reproducible pipeline.
- **"Magnetar" is really "high-$B$".** Our $B>10^{13}$ G cut (83 sources) captures both true magnetars
  ($\sim$30 known) and high-field rotation-powered pulsars; it is a position-in-diagram label, not an
  emission-mechanism classification (most true magnetars are X-ray, not radio, sources).
- **The death line is one model.** The constant-$B/P^2$ line is a representative criterion; published
  death lines span a "death valley" (Chen & Ruderman 1993) and depend on the gap model, so 98.2% should
  be read as "almost all, for a standard line," not a precise number.
- **Only pulsars with a positive measured $\dot P$ enter the diagram.** Of the $\sim$3 500 ATNF entries,
  $\sim$1 400 are excluded: many have no measured $\dot P$, and globular-cluster pulsars often show a
  *negative* apparent $\dot P$ from cluster acceleration (not real spin-down). An acceleration-corrected
  treatment would return some of those; here they are correctly dropped.
- **Selection, not volume-limited.** The ATNF catalogue is the union of many flux- and
  period-limited surveys (MSPs and faint/long-period pulsars are under-represented), so the *relative*
  population sizes are survey-shaped, not intrinsic; the per-class B-fields are robust to this.
- **$\tau=P/2\dot P$ assumes braking index 3 and $P_0\ll P$** — a characteristic age, not the true age.
- **The synthetic fixture is a round-trip code check.** Its `classify_accuracy` (only reported for the
  offline run) confirms the classifier on injected truth labels; it is *not* a validation of the real
  ATNF results, which have no ground-truth class labels.
- **Reproducible:** `python -m jansky_research.ppdot` regenerates the metrics, the P–Ṗ-diagram figure,
  and the macros from the public VizieR catalogue.
