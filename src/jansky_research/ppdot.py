"""The pulsar P--Pdot diagram: surface fields, ages, populations, and the death line.

A rotation-powered pulsar's period :math:`P` and spin-down rate :math:`\\dot P` fix its surface dipole
field :math:`B\\approx3.2\\times10^{19}\\sqrt{P\\dot P}` G, its characteristic age
:math:`\\tau=P/2\\dot P`, and its spin-down luminosity :math:`\\dot E=4\\pi^2 I\\dot P/P^3`. The
:math:`P`--:math:`\\dot P` plane is the H--R diagram of pulsars: ordinary pulsars at
:math:`B\\sim10^{12}` G, **millisecond** pulsars recycled to the bottom-left, **magnetars** at the
top-right, and a **death line** below which radio emission ceases (Bhattacharya & van den Heuvel 1991;
Lorimer & Kramer 2004).

This module reproduces that structure from the ATNF Pulsar Catalogue (VizieR ``B/psr``, public, no
auth), reusing ``pulsarspec.is_millisecond`` and the project's ATNF fetch pattern. Pure NumPy with a
synthetic offline fixture; the real fetch is network-gated.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "B_COEFF",
    "MAGNETAR_B_GAUSS",
    "characteristic_age",
    "classify",
    "death_line",
    "fetch_atnf_ppdot",
    "magnetic_field",
    "population_stats",
    "run",
    "spindown_luminosity",
    "synthetic_population",
]

#: Coefficient in :math:`B=3.2\times10^{19}\sqrt{P\dot P}` G (orthogonal-rotator dipole, I=1e45 g cm^2,
#: R=10 km), the standard surface-field estimate (Lorimer & Kramer 2004, eq. 3.18).
B_COEFF = 3.2e19
#: A pulsar with surface field above this is taken to be a magnetar / high-B object (G).
MAGNETAR_B_GAUSS = 1.0e13
#: Polar-cap death-line threshold in B/P^2 (G s^-2): a source emits while B/P^2 exceeds it.
DEATH_B_OVER_P2 = 1.7e11
_SECONDS_PER_YEAR = 3.155693e7
_MOMENT_OF_INERTIA = 1.0e45  # g cm^2


def magnetic_field(period_s: np.ndarray, pdot: np.ndarray) -> np.ndarray:
    """Surface dipole magnetic field (Gauss): :math:`B=3.2\\times10^{19}\\sqrt{P\\dot P}`."""
    p = np.asarray(period_s, float)
    pd = np.asarray(pdot, float)
    return B_COEFF * np.sqrt(np.clip(p * pd, 0.0, None))


def characteristic_age(period_s: np.ndarray, pdot: np.ndarray) -> np.ndarray:
    """Characteristic age in **years**: :math:`\\tau=P/2\\dot P` (assumes braking index 3, P0 ≪ P)."""
    p = np.asarray(period_s, float)
    pd = np.asarray(pdot, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return p / (2.0 * pd) / _SECONDS_PER_YEAR


def spindown_luminosity(period_s: np.ndarray, pdot: np.ndarray) -> np.ndarray:
    """Spin-down luminosity (erg s⁻¹): :math:`\\dot E=4\\pi^2 I\\dot P/P^3` with :math:`I=10^{45}` g cm²."""
    p = np.asarray(period_s, float)
    pd = np.asarray(pdot, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return 4.0 * np.pi**2 * _MOMENT_OF_INERTIA * pd / p**3


def death_line(period_s: np.ndarray, *, b_over_p2: float = DEATH_B_OVER_P2) -> np.ndarray:
    """The :math:`\\dot P_\\mathrm{death}(P)` below which a pulsar falls past the polar-cap death line.

    A common death-line criterion is constant :math:`B/P^2`: emission ceases when the polar-cap voltage
    drops too low, i.e. :math:`B/P^2 < ` ``b_over_p2`` (G s⁻²). Substituting
    :math:`B=3.2\\times10^{19}\\sqrt{P\\dot P}` gives the threshold spin-down rate
    :math:`\\dot P_\\mathrm{death}=(b_\\mathrm{over\\,p2}/3.2\\times10^{19})^2\\,P^3`; pulsars with
    :math:`\\dot P` above this line are radio-loud (Bhattacharya & van den Heuvel 1991).
    """
    p = np.asarray(period_s, float)
    return (b_over_p2 / B_COEFF) ** 2 * p**3


def classify(period_s: np.ndarray, pdot: np.ndarray) -> np.ndarray:
    """Label each pulsar ``magnetar`` / ``msp`` / ``normal`` by surface field and period.

    Magnetars/high-B objects have :math:`B>` :data:`MAGNETAR_B_GAUSS`; millisecond pulsars are the
    short-period recycled population (``pulsarspec.is_millisecond``); everything else is ``normal``.
    """
    from .pulsarspec import is_millisecond

    p = np.asarray(period_s, float)
    b = magnetic_field(p, pdot)
    msp = np.asarray(is_millisecond(p))
    out = np.where(b > MAGNETAR_B_GAUSS, "magnetar", np.where(msp, "msp", "normal"))
    return out


def population_stats(period_s: np.ndarray, pdot: np.ndarray) -> dict:
    """Per-class counts and median surface field, plus the fraction above the death line."""
    p = np.asarray(period_s, float)
    pd = np.asarray(pdot, float)
    good = np.isfinite(p) & np.isfinite(pd) & (p > 0) & (pd > 0)
    p, pd = p[good], pd[good]
    b = magnetic_field(p, pd)
    cls = classify(p, pd)
    alive = pd > death_line(p)
    out: dict = {"n": int(p.size), "frac_above_death": float(np.mean(alive))}
    for name in ("normal", "msp", "magnetar"):
        m = cls == name
        out[name] = {
            "n": int(m.sum()),
            "median_log_b": float(np.median(np.log10(b[m]))) if m.any() else float("nan"),
        }
    return out


def synthetic_population(n_each: int = 400, *, seed: int = 0) -> dict:
    """Three injected pulsar populations (normal / millisecond / magnetar) with known truth labels.

    Each is a log-normal cloud in :math:`(\\log P, \\log\\dot P)` placed at the textbook location of its
    class, so the analysis recovers the three groups and their characteristic surface fields. Returns
    ``period_s``, ``pdot``, and the ``truth`` class labels.
    """
    rng = np.random.default_rng(seed)
    specs = {
        # class: (mean logP, sd logP, mean logPdot, sd logPdot)
        "normal": (-0.3, 0.3, -14.8, 0.6),
        "msp": (-2.5, 0.25, -19.6, 0.5),
        "magnetar": (0.6, 0.2, -11.0, 0.5),
    }
    periods, pdots, truth = [], [], []
    for name, (mlp, slp, mlpd, slpd) in specs.items():
        lp = rng.normal(mlp, slp, n_each)
        lpd = rng.normal(mlpd, slpd, n_each)
        periods.append(10.0**lp)
        pdots.append(10.0**lpd)
        truth.append(np.full(n_each, name))
    return {
        "period_s": np.concatenate(periods),
        "pdot": np.concatenate(pdots),
        "truth": np.concatenate(truth),
    }


def fetch_atnf_ppdot() -> dict:  # pragma: no cover - network
    """Fetch ATNF P0 (s) and P1 (Pdot) from VizieR ``B/psr``; return period and positive spin-down."""
    import numpy as _np
    from astroquery.vizier import Vizier

    v = Vizier(columns=["PSRJ", "P0", "P1"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("B/psr/psr")[0]
    p = _np.asarray(t["P0"], float)
    pd = _np.asarray(t["P1"], float)
    return {"period_s": p, "pdot": pd}


def run(out: str = ".", *, offline: bool = True) -> dict:
    """Full slice: derive B/age, classify the P--Pdot population, and write outputs."""
    import json
    from pathlib import Path

    if offline:
        pop = synthetic_population()
        source = "synthetic"
        truth: np.ndarray | None = pop["truth"]
    else:  # pragma: no cover - network
        pop = fetch_atnf_ppdot()
        source = "ATNF Pulsar Catalogue (B/psr)"
        truth = None

    stats = population_stats(pop["period_s"], pop["pdot"])
    metrics: dict = {
        "source": source,
        "n_pulsars": stats["n"],
        "frac_above_death": round(stats["frac_above_death"], 3),
        "n_normal": stats["normal"]["n"],
        "n_msp": stats["msp"]["n"],
        "n_magnetar": stats["magnetar"]["n"],
        "median_log_b_normal": round(stats["normal"]["median_log_b"], 2),
        "median_log_b_msp": round(stats["msp"]["median_log_b"], 2),
        "median_log_b_magnetar": round(stats["magnetar"]["median_log_b"], 2),
    }
    if truth is not None:
        good = (pop["period_s"] > 0) & (pop["pdot"] > 0)
        pred = classify(pop["period_s"][good], pop["pdot"][good])
        metrics["classify_accuracy"] = round(float(np.mean(pred == truth[good])), 3)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "ppdot_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _figure(pop["period_s"], pop["pdot"], op / "papers" / "ppdot" / "figures")
    _write_macros(metrics, op / "papers" / "ppdot" / "generated" / "macros.tex")
    return metrics


def _figure(period_s, pdot, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    p = np.asarray(period_s, float)
    pd = np.asarray(pdot, float)
    good = np.isfinite(p) & np.isfinite(pd) & (p > 0) & (pd > 0)
    p, pd = p[good], pd[good]
    cls = classify(p, pd)
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    colors = {"normal": "0.5", "msp": "C0", "magnetar": "C3"}
    for name, c in colors.items():
        m = cls == name
        ax.scatter(p[m], pd[m], s=6, c=c, label=name)
    pg = np.logspace(np.log10(p.min()), np.log10(p.max()), 100)
    ax.plot(pg, death_line(pg), "k--", lw=0.9, label="death line")
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="period $P$ (s)",
        ylabel=r"$\dot P$ (s s$^{-1}$)",
        title="The pulsar P--Pdot diagram",
    )
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "ppdot.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.ppdot._write_macros -- do not edit by hand.",
        rf"\newcommand{{\ppSource}}{{{m['source']}}}",
        rf"\newcommand{{\ppN}}{{{m['n_pulsars']}}}",
        rf"\newcommand{{\ppFracAlive}}{{{m['frac_above_death']}}}",
        rf"\newcommand{{\ppNnormal}}{{{m['n_normal']}}}",
        rf"\newcommand{{\ppNmsp}}{{{m['n_msp']}}}",
        rf"\newcommand{{\ppNmagnetar}}{{{m['n_magnetar']}}}",
        rf"\newcommand{{\ppLogBnormal}}{{{m['median_log_b_normal']}}}",
        rf"\newcommand{{\ppLogBmsp}}{{{m['median_log_b_msp']}}}",
        rf"\newcommand{{\ppLogBmagnetar}}{{{m['median_log_b_magnetar']}}}",
        rf"\newcommand{{\ppAccuracy}}{{{_fmt('classify_accuracy')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(description="The pulsar P--Pdot diagram (ATNF).")
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args(argv)
    print(json.dumps(run(args.out, offline=args.offline), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
