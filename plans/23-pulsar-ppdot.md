# 23 — The pulsar P–Ṗ diagram: B-fields, ages, populations, and the death line (ATNF)

Status: 🚧 in progress (tooling + real ATNF fetch + recover-a-known done; GATE-2 + paper next)

## Context

A rotation-powered pulsar's period $P$ and period derivative $\dot P$ encode its physics: the surface
dipole field $B\approx3.2\times10^{19}\sqrt{P\dot P}$ G, the characteristic age $\tau=P/2\dot P$, and the
spin-down luminosity $\dot E=4\pi^2 I\dot P/P^3$. Plotting $\dot P$ against $P$ — the **P–Ṗ diagram** —
is the H-R diagram of pulsars: ordinary pulsars cluster at $P\sim0.1$–1 s and $B\sim10^{12}$ G,
**millisecond pulsars** sit at the bottom-left ($P\lesssim30$ ms, $B\sim10^8$ G — recycled by accretion),
**magnetars** at the top-right ($B\gtrsim10^{14}$ G), and a **death line** marks where the polar-cap
voltage can no longer sustain the pair production that powers radio emission (Bhattacharya & van den
Heuvel 1991; Lorimer & Kramer 2004).

This slice reproduces that structure from the **ATNF Pulsar Catalogue** (VizieR `B/psr`, public, no
auth), reusing the project's ATNF fetch pattern and `pulsarspec.is_millisecond`. The deliverable is a
tested tool + a recover-a-known (the three populations separate cleanly in B-field, and nearly all radio
pulsars lie above the death line), not a discovery.

## Reuse

- `jansky_research.pulsarspec.is_millisecond` (period threshold) and the `fetch_atnf` VizieR pattern.
- the slice/paper/test conventions (synthetic offline fixture, `_write_macros`, `_figure`, AASTeX
  paper, gitignored generated artefacts).

## Deliverables

- `src/jansky_research/ppdot.py`:
  - `magnetic_field`, `characteristic_age`, `spindown_luminosity` — the standard P, Ṗ derivations.
  - `death_line` — the $\dot P_\mathrm{death}(P)$ threshold from a constant $B/P^2$ polar-cap model.
  - `classify` — magnetar / millisecond / normal by B-field and period.
  - `population_stats` — per-class counts, median B, and the fraction above the death line.
  - `synthetic_population` — offline fixture: three injected clusters with known truth labels, so the
    tests recover the classes and their B distributions.
  - `fetch_atnf_ppdot` (network, `# pragma: no cover`) — VizieR `B/psr` P0, P1 → P, Ṗ.
  - `run(offline=...)`, `_figure` (the P–Ṗ diagram + death line + B/age guide lines), `_write_macros`,
    `_main`.
- `tests/test_ppdot.py` — synthetic-fixture tests to the 85% floor; no network.
- `data.py` registry note for the ATNF catalogue (if not already present).

## Approach

1. **Tooling (this step).** Pure-NumPy derivations + classification, validated on a synthetic
   three-population sky whose injected classes and B distributions are recovered.
2. **Real-data run (next).** `run(offline=False)` fetches ATNF P0/P1, drops non-positive Ṗ, derives B,
   τ, Ė, classifies, and reports the population stats and the death-line fraction.
3. **Recover-a-known.** The three populations separate in B-field (MSP $\sim10^{8.5}$, normal
   $\sim10^{12.5}$, magnetar $\sim10^{14}$ G), and $\gtrsim$95% of radio pulsars lie above the death line.
4. **GATE-2** before write-up — the structure survives the catalogue-selection and death-line-model
   caveats and a literature cross-check.
5. **Write-up** as `papers/ppdot/`.

## Verification

- `make test` / `cov` green on synthetic fixtures (no network), 85% floor; `ruff` + `mypy` clean.
- Offline `run` recovers the three injected classes and their median B-fields.
- `magnetic_field`/`characteristic_age` match textbook values (e.g. the Crab: $B\sim4\times10^{12}$ G,
  $\tau\sim1250$ yr).
- (Real-data, later) the ATNF P–Ṗ diagram reproduces the known population structure; GATE-2 sign-off.
