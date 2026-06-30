"""Synthesis: a solar type III electron beam from the corona to 0.4 AU, geometrically validated.

Four slices already track type III bursts with the same **drift-to-distance** idea --- the emission
sits near the local plasma frequency, so the frequency drift, inverted through a density model, gives
the beam's heliocentric distance: ``solarbursts`` (e-Callisto, corona, Newkirk model), ``windwaves``
(Wind/WAVES, to the Alfven surface, Leblanc model), ``swaves`` (STEREO/WAVES HFR, to 0.4 AU, Leblanc),
and ``triangulate`` (STEREO-A+B direction-finding, an **independent geometric** distance). This module
orchestrates the four into one figure and one macro set: a unified distance ladder from ~100 MHz
(corona) to 0.125 MHz (0.4 AU) over the Newkirk and Leblanc density models, plus the key cross-check
--- ``swaves`` and ``triangulate`` analyse the **same 2013-05-15 event**, so the density-model distance
the whole ladder rests on is independently confirmed by geometry (corr ~0.99).

No new physics: it calls each slice's tested ``run`` and reuses ``triangulate``'s track for the
cross-check. Pure-NumPy/matplotlib; offline it composes the four synthetic fixtures (so CI builds with
no network), and ``reproduce`` runs the four recover-a-known events on public data.
"""

from __future__ import annotations

import numpy as np

from . import solarbursts, swaves, triangulate, windwaves

__all__ = ["collect_metrics", "crosscheck_track", "run"]

# the canonical recover-a-known events each slice is reproduced on (see the Makefile reproduce target)
WINDWAVES_DATE = "20031028"
SWAVES_DATE = "20130515"
TRIANGULATE_DATE = "20130515"  # the SAME event as swaves -- the centrepiece cross-check


def collect_metrics(out: str, *, offline: bool) -> dict:
    """Run the four type III slices (offline-synthetic or on their real events) and collect metrics."""
    if offline:
        return {
            "solarbursts": solarbursts.run(out, offline=True),
            "windwaves": windwaves.run(out, offline=True),
            "swaves": swaves.run(out, offline=True),
            "triangulate": triangulate.run(out, offline=True),
        }
    # solarbursts.RECOVER_EVENT spelled out as literals (a mixed-type dict won't unpack into typed kwargs)
    return {  # pragma: no cover - network
        "solarbursts": solarbursts.run(
            out, offline=False, station="BIR", date="20110914", hhmm="1150", harmonic=2, fold=1.0
        ),
        "windwaves": windwaves.run(out, offline=False, date=WINDWAVES_DATE, receiver="rad2"),
        "swaves": swaves.run(out, offline=False, date=SWAVES_DATE, spacecraft="a"),
        "triangulate": triangulate.run(out, offline=False, date=TRIANGULATE_DATE),
    }


def crosscheck_track(*, offline: bool, harmonic: int = 2) -> dict:
    """The 2013-05-15 geometric-vs-plasma cross-check: per-frequency ``r_geom`` and ``r_plasma``.

    Offline this uses ``triangulate``'s synthetic event; on real data it triangulates the STEREO-A+B
    direction-finding for 2013-05-15. Returns the kept-channel arrays from
    :func:`triangulate.triangulate_track` (``freq_mhz``, ``r_geom``, ``r_plasma``).
    """
    if offline:
        ev = triangulate.synthetic_event(harmonic=harmonic)
        spec_a, spec_b = ev["spec_a"], ev["spec_b"]
    else:  # pragma: no cover - network
        spec_a = triangulate.fetch_stereo_df(TRIANGULATE_DATE, spacecraft="a")
        spec_b = triangulate.fetch_stereo_df(TRIANGULATE_DATE, spacecraft="b")
    return triangulate.triangulate_track(spec_a, spec_b, harmonic=harmonic)


def _model_curves(harmonic: int = 2) -> dict:
    """Newkirk (corona) and Leblanc (heliosphere) heliocentric radius vs emission frequency (MHz)."""
    from jansky import solar

    f_corona = np.logspace(np.log10(15.0), np.log10(300.0), 200)  # MHz
    r_corona = solar.newkirk_radius(solar.density_from_plasma_frequency(f_corona / harmonic))
    f_helio = np.logspace(np.log10(0.1), np.log10(20.0), 200)
    r_helio = windwaves.emission_radius(f_helio, harmonic=harmonic)
    return {"f_corona": f_corona, "r_corona": r_corona, "f_helio": f_helio, "r_helio": r_helio}


def run(out: str = ".", *, offline: bool = True, harmonic: int = 2) -> dict:
    """Full synthesis: orchestrate the four slices, build the ladder + cross-check, emit macros."""
    import json
    from pathlib import Path

    m = collect_metrics(out, offline=offline)
    track = crosscheck_track(offline=offline, harmonic=harmonic)

    rg = np.asarray(track.get("r_geom", []), float)
    rp = np.asarray(track.get("r_plasma", []), float)
    corr = (
        float(np.corrcoef(rg, rp)[0, 1])
        if rg.size >= 3 and np.ptp(rg) > 0 and np.ptp(rp) > 0
        else None
    )

    def _g(slice_name: str, key: str):
        return m[slice_name].get(key)

    r_au = windwaves.R_AU_RSUN
    # overall reach: smallest corona radius to the largest interplanetary radius
    r_lo = _g("solarbursts", "r_lo_rsun")
    r_hi = max(v for v in (_g("swaves", "r_hi_rsun"), _g("triangulate", "r_hi_rsun")) if v)
    metrics: dict = {
        "n_instruments": 4,
        "crosscheck_event": "2013-05-15 (STEREO/WAVES + STEREO-A+B triangulation)",
        "f_hi_mhz": _g("solarbursts", "f_hi_mhz"),  # corona, highest frequency
        "f_lo_mhz": _g("swaves", "f_lo_mhz"),  # interplanetary, lowest frequency
        "corona_r_lo": r_lo,
        "corona_r_hi": _g("solarbursts", "r_hi_rsun"),
        "corona_speed_c": _g("solarbursts", "speed_c"),
        "helio_r_hi": _g("windwaves", "r_hi_rsun"),
        "helio_speed_c": _g("windwaves", "speed_c"),
        "ip_r_hi_rsun": _g("swaves", "r_hi_rsun"),
        "ip_r_hi_au": _g("swaves", "r_hi_au"),
        "ip_speed_c": _g("swaves", "speed_c"),
        "geom_r_hi_rsun": _g("triangulate", "r_hi_rsun"),
        "geom_r_hi_au": _g("triangulate", "r_hi_au"),
        "geom_corr": round(corr, 3) if corr is not None else _g("triangulate", "corr_geom_plasma"),
        "geom_ratio": _g("triangulate", "ratio_geom_plasma"),
        "overall_r_hi_au": round(r_hi / r_au, 3) if r_hi else None,
    }

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    (op / "results" / "type3synthesis_metrics.json").write_text(
        json.dumps(metrics, indent=2) + "\n"
    )
    _figure(m, track, harmonic, op / "papers" / "type3synthesis" / "figures")
    _write_macros(metrics, op / "papers" / "type3synthesis" / "generated" / "macros.tex")
    return metrics


def _figure(m: dict, track: dict, harmonic: int, out_dir) -> None:
    from pathlib import Path

    from .report import _agg

    plt = _agg()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    curves = _model_curves(harmonic)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 4.0))

    # Left: the unified distance ladder (heliocentric radius vs emission frequency)
    ax1.plot(
        curves["f_corona"], curves["r_corona"], "-", color="0.6", lw=1, label="Newkirk (corona)"
    )
    ax1.plot(curves["f_helio"], curves["r_helio"], "--", color="0.4", lw=1, label="Leblanc (helio)")
    seg = [
        ("solarbursts", "e-Callisto", "C1"),
        ("windwaves", "Wind/WAVES", "C0"),
        ("swaves", "STEREO/WAVES", "C2"),
        ("triangulate", "STEREO A+B (geom)", "C3"),
    ]
    for name, label, color in seg:
        d = m[name]
        flo, fhi = d.get("f_lo_mhz"), d.get("f_hi_mhz")
        rlo, rhi = d.get("r_lo_rsun"), d.get("r_hi_rsun")
        if None in (flo, fhi, rlo, rhi):
            continue
        ax1.plot([fhi, flo], [rlo, rhi], "o-", color=color, ms=4, lw=2, label=label)
    ax1.set(
        xscale="log",
        yscale="log",
        xlabel="emission frequency (MHz)",
        ylabel=r"heliocentric distance ($R_\odot$)",
        title="Type III beam: corona to 0.4 AU",
    )
    ax1.axhline(windwaves.R_AU_RSUN, color="k", ls=":", lw=0.6)
    ax1.text(ax1.get_xlim()[0], windwaves.R_AU_RSUN * 1.05, "1 AU", fontsize=7)
    ax1.legend(fontsize=7, loc="lower left")

    # Right: the 2013-05-15 geometric-vs-plasma cross-check
    rg = np.asarray(track.get("r_geom", []), float)
    rp = np.asarray(track.get("r_plasma", []), float)
    if rg.size:
        ax2.plot(rp, rg, "o", color="C3", ms=4)
        lim = [min(rp.min(), rg.min()), max(rp.max(), rg.max())]
        ax2.plot(lim, lim, "k--", lw=0.8, label="1:1")
    ax2.set(
        xscale="log",
        yscale="log",
        xlabel=r"plasma-frequency distance ($R_\odot$)",
        ylabel=r"geometric (triangulated) distance ($R_\odot$)",
        title="2013-05-15 cross-check",
    )
    ax2.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "type3synthesis.pdf")
    plt.close(fig)


def _write_macros(m: dict, path) -> None:
    from pathlib import Path

    def _fmt(key: str) -> str:
        val = m.get(key)
        return "--" if val is None else str(val)

    lines = [
        "% Auto-generated by jansky_research.type3synthesis._write_macros -- do not edit by hand.",
        rf"\newcommand{{\synFhi}}{{{_fmt('f_hi_mhz')}}}",
        rf"\newcommand{{\synFlo}}{{{_fmt('f_lo_mhz')}}}",
        rf"\newcommand{{\synCoronaRlo}}{{{_fmt('corona_r_lo')}}}",
        rf"\newcommand{{\synCoronaRhi}}{{{_fmt('corona_r_hi')}}}",
        rf"\newcommand{{\synCoronaSpeed}}{{{_fmt('corona_speed_c')}}}",
        rf"\newcommand{{\synHelioRhi}}{{{_fmt('helio_r_hi')}}}",
        rf"\newcommand{{\synHelioSpeed}}{{{_fmt('helio_speed_c')}}}",
        rf"\newcommand{{\synIpRhi}}{{{_fmt('ip_r_hi_rsun')}}}",
        rf"\newcommand{{\synIpRhiAU}}{{{_fmt('ip_r_hi_au')}}}",
        rf"\newcommand{{\synIpSpeed}}{{{_fmt('ip_speed_c')}}}",
        rf"\newcommand{{\synGeomRhi}}{{{_fmt('geom_r_hi_rsun')}}}",
        rf"\newcommand{{\synGeomRhiAU}}{{{_fmt('geom_r_hi_au')}}}",
        rf"\newcommand{{\synGeomCorr}}{{{_fmt('geom_corr')}}}",
        rf"\newcommand{{\synGeomRatio}}{{{_fmt('geom_ratio')}}}",
        rf"\newcommand{{\synOverallRhiAU}}{{{_fmt('overall_r_hi_au')}}}",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="Synthesis: a type III beam from the corona to 0.4 AU, geometrically validated."
    )
    p.add_argument("--out", default=".")
    p.add_argument("--offline", action="store_true")
    p.add_argument("--harmonic", type=int, default=2)
    args = p.parse_args(argv)
    metrics = run(args.out, offline=args.offline, harmonic=args.harmonic)
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
