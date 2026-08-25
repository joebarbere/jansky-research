"""Synthesis: a reproducible type III drift-to-distance framework, corona to 0.4 AU.

Four slices already track type III bursts with the same **drift-to-distance** idea --- the emission
sits near the local plasma frequency, so the frequency drift, inverted through a density model, gives
the beam's heliocentric distance: ``solarbursts`` (e-Callisto, corona, Newkirk model), ``windwaves``
(Wind/WAVES, to the Alfven surface, Leblanc model), ``swaves`` (STEREO/WAVES HFR, to 0.4 AU, Leblanc),
and ``triangulate`` (STEREO-A+B direction-finding, an **independent geometric** distance). This module
orchestrates the four into one figure and one macro set: a unified distance ladder from ~100 MHz
(corona) to 0.125 MHz (0.4 AU) over the Newkirk and Leblanc density models, plus the key cross-check
--- ``swaves`` and ``triangulate`` analyse the **same 2013-05-15 event**, so the density-model
distance the ladder rests on is geometrically bounded: a constant ~13 R_sun additive offset (a
few-degree direction-finding bias), after which the Leblanc level is reproduced to ~12%.

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

    def _bracket(slice_name: str, key: str) -> tuple:
        grid = m[slice_name].get("speed_grid") or []
        vals = [g.get(key) for g in grid if g.get(key) is not None]
        return (min(vals), max(vals)) if vals else (None, None)

    # The ladder's reach is the PLASMA leg's band edge through the model (0.384 AU); the
    # geometric point (106 R_sun) carries the constant direction-finding offset and is
    # reported separately, never as the headline reach -- the paper cannot use a number as
    # its reach and disown its calibration three paragraphs later.
    metrics: dict = {
        # Provenance, so both merge guards can tell the two run modes apart. Without it this
        # was the only slice a forced offline rebuild could still overwrite.
        "source": (
            "synthetic type-III ladder (four synthesised offline slices)"
            if offline
            else "e-Callisto + STEREO/WAVES + Wind/WAVES + STEREO-A+B triangulation"
        ),
        "n_instruments": 4,
        "crosscheck_event": "2013-05-15 (STEREO/WAVES + STEREO-A+B triangulation)",
        "f_hi_mhz": _g("solarbursts", "f_hi_mhz"),  # corona, highest frequency
        "f_lo_mhz": _g("swaves", "f_lo_mhz"),  # interplanetary, lowest frequency
        "corona_r_lo": _g("solarbursts", "r_lo_rsun"),
        "corona_r_hi": _g("solarbursts", "r_hi_rsun"),
        "corona_speed_c": _g("solarbursts", "speed_c"),
        "corona_speed_c_min": _g("solarbursts", "speed_c_min"),
        "corona_speed_c_max": _g("solarbursts", "speed_c_max"),
        # the corona radii come from the FITTED band, so the figure/prose must use it too
        # (the "figure draws the fit it captions" fix, propagated here)
        "corona_fit_f_lo_mhz": _g("solarbursts", "fit_f_lo_mhz") or _g("solarbursts", "f_lo_mhz"),
        "corona_fit_f_hi_mhz": _g("solarbursts", "fit_f_hi_mhz") or _g("solarbursts", "f_hi_mhz"),
        "helio_r_hi": _g("windwaves", "r_hi_rsun"),
        "helio_speed_c": _g("windwaves", "speed_c"),
        "helio_speed_c_se": _g("windwaves", "speed_c_se"),
        "ip_r_hi_rsun": _g("swaves", "r_hi_rsun"),
        "ip_r_hi_au": _g("swaves", "r_hi_au"),
        "ip_speed_c": _g("swaves", "speed_c"),
        "ip_speed_c_se": _g("swaves", "speed_c_se"),
        "geom_r_hi_rsun": _g("triangulate", "r_hi_rsun"),
        "geom_r_hi_au": _g("triangulate", "r_hi_au"),
        "geom_corr": round(corr, 3) if corr is not None else _g("triangulate", "corr_geom_plasma"),
        "geom_ratio": _g("triangulate", "ratio_geom_plasma"),
        "overall_r_hi_au": _g("swaves", "r_hi_au"),
    }
    # the (harmonic x density) systematic brackets the siblings now propagate: the synthesis
    # must not quote bare per-leg numbers its own sources refuse to quote bare
    for leg, sl in (("helio", "windwaves"), ("ip", "swaves")):
        lo, hi = _bracket(sl, "speed_c")
        metrics[f"{leg}_grid_speed_lo"] = lo
        metrics[f"{leg}_grid_speed_hi"] = hi
        rlo, rhi = _bracket(sl, "r_hi_au")
        metrics[f"{leg}_grid_reach_au_lo"] = rlo
        metrics[f"{leg}_grid_reach_au_hi"] = rhi
    lo, hi = _bracket("solarbursts", "speed_c")
    metrics["corona_grid_speed_lo"] = lo
    metrics["corona_grid_speed_hi"] = hi
    # the geometric comparison in the additive framing (see papers/triangulate): a constant
    # offset at the scale of the ray miss, not a scale factor
    if rg.size >= 6:
        add = triangulate.additive_vs_multiplicative(rg, rp)
        metrics["geom_diff_med_rsun"] = add.get("diff_med_rsun")
        metrics["geom_diff_std_rsun"] = add.get("diff_std_rsun")
        metrics["geom_ols_slope"] = add.get("ols_slope")
        metrics["geom_rms_additive_rsun"] = add.get("rms_additive_rsun")
        lf, lg2, lp2 = np.log10(track["freq_mhz"]), np.log10(rg), np.log10(rp)
        _ = lf
        metrics["geom_loglog_slope"] = round(float(np.polyfit(lp2, lg2, 1)[0]), 3)
    # the Newkirk/Leblanc handoff discontinuity, computed at the two legs' ACTUAL overlap
    # band rather than hand-typed as "~50% at 15-20 MHz"
    from jansky import solar

    f_corona_lo = _g("solarbursts", "f_lo_mhz")
    f_helio_hi = _g("windwaves", "f_hi_mhz")
    if f_corona_lo and f_helio_hi:
        lo_f, hi_f = sorted((float(f_corona_lo), float(f_helio_hi)))
        ratios = []
        for fq in (lo_f, hi_f):
            r_new = float(solar.newkirk_radius(solar.density_from_plasma_frequency(fq / 2.0)))
            r_leb = float(windwaves.emission_radius(np.array([fq]), harmonic=2)[0])
            ratios.append(r_new / r_leb)
        metrics["handoff_f_lo_mhz"] = round(lo_f, 2)
        metrics["handoff_f_hi_mhz"] = round(hi_f, 2)
        metrics["handoff_pct_lo"] = round(100.0 * (min(ratios) - 1.0), 0)
        metrics["handoff_pct_hi"] = round(100.0 * (max(ratios) - 1.0), 0)

    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    from .report import write_results

    write_results(metrics, op / "results" / "type3synthesis_metrics.json")
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
        # the radii come from the FITTED band where one exists (solarbursts sigma-clips its
        # ridge), so the segment must span the fitted frequencies, not the full detection band
        # -- otherwise the corona segment sits a factor ~2 off the Newkirk curve it overplots
        flo = d.get("fit_f_lo_mhz") or d.get("f_lo_mhz")
        fhi = d.get("fit_f_hi_mhz") or d.get("f_hi_mhz")
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
        if rg.size >= 6:
            off = float(np.median(rg - rp))
            xs = np.linspace(lim[0], lim[1], 50)
            ax2.plot(xs, xs + off, "-", color="C0", lw=0.9, label=f"1:1 + {off:.0f} $R_\\odot$")
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

    # This slice synthesises four others, and its offline leg produces a different ladder from
    # the real one (\synCoronaSpeed 0.1347 -> 0.3002, \synHelioRhi 10.25 -> 48.37). It carried
    # NO provenance at all -- no `source` in its metrics and no `*Source` macro -- so it was the
    # one slice neither `preserve_live_results` nor `preserve_live_macros` could protect, and
    # the only results file a forced offline rebuild still changed. Emitting the marker is what
    # makes both guards able to see it.
    # Every ladder value is mode-dependent: the offline leg synthesises all four sibling
    # slices and produces a different ladder (corona speed 0.1347 real vs 0.3002 offline).
    # Namespaced so a synthetic rebuild cannot write its ladder under the names the paper
    # cites; the shared-name version was a live clobber this slice could not even detect,
    # because it emitted no provenance macro at all.
    real = not str(m.get("source", "")).lower().startswith("synthetic")
    ns, other = ("synReal", "synSyn") if real else ("synSyn", "synReal")
    lines = [
        "% Auto-generated by jansky_research.type3synthesis._write_macros -- do not edit by hand.",
        "% Ladder values are mode-dependent and namespaced; only the active run fills its side.",
        rf"\newcommand{{\synSource}}{{{m.get('source', '--')}}}",
        rf"\newcommand{{\{ns}Fhi}}{{{_fmt('f_hi_mhz')}}}",
        rf"\newcommand{{\{other}Fhi}}{{--}}",
        rf"\newcommand{{\{ns}Flo}}{{{_fmt('f_lo_mhz')}}}",
        rf"\newcommand{{\{other}Flo}}{{--}}",
        rf"\newcommand{{\{ns}CoronaRlo}}{{{_fmt('corona_r_lo')}}}",
        rf"\newcommand{{\{other}CoronaRlo}}{{--}}",
        rf"\newcommand{{\{ns}CoronaRhi}}{{{_fmt('corona_r_hi')}}}",
        rf"\newcommand{{\{other}CoronaRhi}}{{--}}",
        rf"\newcommand{{\{ns}CoronaSpeed}}{{{_fmt('corona_speed_c')}}}",
        rf"\newcommand{{\{other}CoronaSpeed}}{{--}}",
        rf"\newcommand{{\{ns}HelioRhi}}{{{_fmt('helio_r_hi')}}}",
        rf"\newcommand{{\{other}HelioRhi}}{{--}}",
        rf"\newcommand{{\{ns}HelioSpeed}}{{{_fmt('helio_speed_c')}}}",
        rf"\newcommand{{\{other}HelioSpeed}}{{--}}",
        rf"\newcommand{{\{ns}IpRhi}}{{{_fmt('ip_r_hi_rsun')}}}",
        rf"\newcommand{{\{other}IpRhi}}{{--}}",
        rf"\newcommand{{\{ns}IpRhiAU}}{{{_fmt('ip_r_hi_au')}}}",
        rf"\newcommand{{\{other}IpRhiAU}}{{--}}",
        rf"\newcommand{{\{ns}IpSpeed}}{{{_fmt('ip_speed_c')}}}",
        rf"\newcommand{{\{other}IpSpeed}}{{--}}",
        rf"\newcommand{{\{ns}GeomRhi}}{{{_fmt('geom_r_hi_rsun')}}}",
        rf"\newcommand{{\{other}GeomRhi}}{{--}}",
        rf"\newcommand{{\{ns}GeomRhiAU}}{{{_fmt('geom_r_hi_au')}}}",
        rf"\newcommand{{\{other}GeomRhiAU}}{{--}}",
        rf"\newcommand{{\{ns}GeomCorr}}{{{_fmt('geom_corr')}}}",
        rf"\newcommand{{\{other}GeomCorr}}{{--}}",
        rf"\newcommand{{\{ns}GeomRatio}}{{{_fmt('geom_ratio')}}}",
        rf"\newcommand{{\{other}GeomRatio}}{{--}}",
        rf"\newcommand{{\{ns}OverallRhiAU}}{{{_fmt('overall_r_hi_au')}}}",
        rf"\newcommand{{\{other}OverallRhiAU}}{{--}}",
    ]
    # round-8 additions: per-leg errors, the (harmonic x density) brackets the siblings
    # propagate, the additive geometric comparison, and the computed model handoff
    extra = (
        ("CoronaSpeedMin", "corona_speed_c_min"),
        ("CoronaSpeedMax", "corona_speed_c_max"),
        ("CoronaFitFlo", "corona_fit_f_lo_mhz"),
        ("CoronaFitFhi", "corona_fit_f_hi_mhz"),
        ("CoronaGridSpeedLo", "corona_grid_speed_lo"),
        ("CoronaGridSpeedHi", "corona_grid_speed_hi"),
        ("HelioSpeedErr", "helio_speed_c_se"),
        ("HelioGridSpeedLo", "helio_grid_speed_lo"),
        ("HelioGridSpeedHi", "helio_grid_speed_hi"),
        ("HelioGridReachLo", "helio_grid_reach_au_lo"),
        ("HelioGridReachHi", "helio_grid_reach_au_hi"),
        ("IpSpeedErr", "ip_speed_c_se"),
        ("IpGridSpeedLo", "ip_grid_speed_lo"),
        ("IpGridSpeedHi", "ip_grid_speed_hi"),
        ("IpGridReachLo", "ip_grid_reach_au_lo"),
        ("IpGridReachHi", "ip_grid_reach_au_hi"),
        ("GeomDiff", "geom_diff_med_rsun"),
        ("GeomDiffStd", "geom_diff_std_rsun"),
        ("GeomSlope", "geom_ols_slope"),
        ("GeomRmsAdd", "geom_rms_additive_rsun"),
        ("GeomLogSlope", "geom_loglog_slope"),
        ("HandoffFlo", "handoff_f_lo_mhz"),
        ("HandoffFhi", "handoff_f_hi_mhz"),
        ("HandoffPctLo", "handoff_pct_lo"),
        ("HandoffPctHi", "handoff_pct_hi"),
    )
    for macro, key in extra:
        lines.append(rf"\newcommand{{\{ns}{macro}}}{{{_fmt(key)}}}")
        lines.append(rf"\newcommand{{\{other}{macro}}}{{--}}")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Merge rather than overwrite: this run knows only its own mode's metrics and
    # would otherwise blank the other mode's macros with '--'. `make figures`
    # runs every slice offline in the repo root, so without this an offline
    # rebuild silently empties this paper. See report.preserve_live_macros.
    from .report import preserve_live_macros

    p.write_text(preserve_live_macros("\n".join(lines) + "\n", p))


def _main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin CLI
    import argparse
    import json

    p = argparse.ArgumentParser(
        description="Synthesis: a type III drift-to-distance framework, corona to 0.4 AU."
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
