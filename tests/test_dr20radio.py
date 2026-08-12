"""Tests for jansky_research.dr20radio — DR20 BHM radio census. No network."""

from __future__ import annotations

import numpy as np
import pytest

from jansky_research import dr20radio


def test_select_quasars_filters_and_flags_cartons():
    cls = np.array(["QSO", "GALAXY", "QSO ", "STAR", "QSO"])
    zw = np.array([0, 0, 0, 0, 5])
    carton = np.array(
        ["bhm_spiders_agn", "x", "openfibertargets_bhm_racsradio_boss", "y", "bhm_csc"]
    )
    quasar, radio = dr20radio.select_quasars(cls, zw, carton)
    assert quasar.tolist() == [True, False, True, False, False]  # ZWARNING!=0 excluded
    assert radio.tolist() == [False, False, True, False, False]


def test_crossmatch_known_geometry():
    ra_q = np.array([180.0, 181.0])
    dec_q = np.array([10.0, 10.0])
    # radio source 1" from quasar 0; nothing near quasar 1
    ra_r = np.array([180.0 + 1.0 / 3600.0 / np.cos(np.deg2rad(10.0)), 200.0])
    dec_r = np.array([10.0, -5.0])
    m, sep, idx = dr20radio.crossmatch(ra_q, dec_q, ra_r, dec_r, radius_arcsec=2.5)
    assert m.tolist() == [True, False]
    assert sep[0] == pytest.approx(1.0, abs=0.05)
    assert idx[0] == 0  # nearest-neighbour index enables counterpart flux lookup


def test_wilson_interval_properties():
    p, lo, hi = dr20radio.wilson_interval(10, 100)
    assert p == pytest.approx(0.1)
    assert lo < p < hi
    assert dr20radio.wilson_interval(0, 0)[0] != dr20radio.wilson_interval(0, 0)[0]  # nan


def test_synthetic_round_trip_recovers_fraction_and_carton_circularity():
    s = dr20radio.synthetic_survey(seed=1)
    quasar, radio_carton = dr20radio.select_quasars(s["cls"], s["zwarning"], s["firstcarton"])
    assert quasar.all()
    census = ~radio_carton
    m, _, _ = dr20radio.crossmatch(
        s["ra_q"][census],
        s["dec_q"][census],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
    )
    fm = dr20radio.false_match_rate(
        s["ra_q"][census],
        s["dec_q"][census],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
        n_trials=5,
        seed=2,
    )
    corrected = float(np.mean(m)) - fm["rate"]
    assert corrected == pytest.approx(s["true_fraction"], abs=0.03)
    # the radio-carton subset matches ~always (counterpart by construction) — the
    # circularity that must stay OUT of the census fraction
    mc, _, _ = dr20radio.crossmatch(
        s["ra_q"][radio_carton],
        s["dec_q"][radio_carton],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
    )
    assert float(np.mean(mc)) > 0.88  # 1" Rayleigh jitter puts ~4% beyond 2.5"
    assert float(np.mean(mc)) - float(np.mean(m)) > 0.5  # excluding them matters


def test_false_match_rate_scales_with_density():
    s = dr20radio.synthetic_survey(seed=3, n_radio=2000, counterpart_fraction=0.0)
    fm_lo = dr20radio.false_match_rate(
        s["ra_q"], s["dec_q"], s["ra_r"], s["dec_r"], n_trials=4, seed=4, radius_arcsec=30.0
    )
    s2 = dr20radio.synthetic_survey(seed=3, n_radio=20000, counterpart_fraction=0.0)
    fm_hi = dr20radio.false_match_rate(
        s2["ra_q"], s2["dec_q"], s2["ra_r"], s2["dec_r"], n_trials=4, seed=4, radius_arcsec=30.0
    )  # 30" radius gives enough chance matches for the scaling comparison to have power
    assert fm_hi["rate"] > fm_lo["rate"]  # denser radio catalog -> more chance matches
    assert fm_hi["rate"] == pytest.approx(10 * fm_lo["rate"], rel=0.6)


def test_detection_fraction_binning():
    z = np.array([0.5, 0.6, 1.5, 1.6, 3.0])
    matched = np.array([True, False, True, True, False])
    out = dr20radio.detection_fraction(z, matched, bins=np.array([0.0, 1.0, 2.0, 4.0]))
    assert out["n"] == [2, 2, 1]
    assert out["k"] == [1, 2, 0]
    assert out["frac"][1] == pytest.approx(1.0)


def test_read_spall_quasars_on_synthetic_fits(tmp_path):
    from astropy.io import fits

    n = 6
    cols = [
        fits.Column(name="RACAT", format="D", array=np.linspace(10, 15, n)),
        fits.Column(name="DECCAT", format="D", array=np.linspace(-50, 50, n)),
        fits.Column(name="Z", format="D", array=np.linspace(0.5, 3.0, n)),
        fits.Column(name="ZWARNING", format="J", array=np.array([0, 0, 0, 4, 0, 0])),
        fits.Column(
            name="CLASS",
            format="10A",
            array=np.array(["QSO", "QSO", "GALAXY", "QSO", "QSO", "QSO"]),
        ),
        fits.Column(
            name="FIRSTCARTON",
            format="40A",
            array=np.array(
                ["bhm_spiders", "openfibertargets_bhm_racsradio_boss", "x", "y", "z", "w"]
            ),
        ),
        fits.Column(
            name="OBS", format="3A", array=np.array(["APO", "LCO", "APO", "APO", "LCO", "APO"])
        ),
    ]
    path = tmp_path / "spall.fits"
    fits.BinTableHDU.from_columns(cols).writeto(path)
    q = dr20radio.read_spall_quasars(str(path))
    assert q["n_total_rows"] == 6
    assert q["ra"].size == 4  # 6 - GALAXY - ZWARNING!=0
    assert q["radio_carton"].sum() == 1
    assert set(q["obs"]) <= {"APO", "LCO"}


def test_two_survey_synthetic_fixes_the_carton_blind_spot():
    s = dr20radio.synthetic_two_surveys(seed=5, fade_fraction=0.35)
    _, radio_carton = dr20radio.select_quasars(s["cls"], s["zwarning"], s["firstcarton"])
    # vs the SELECTING survey: ~100% (counterpart by construction)
    m_sel, _, _ = dr20radio.crossmatch(
        s["ra_q"][radio_carton],
        s["dec_q"][radio_carton],
        s["ra_r"],
        s["dec_r"],
        radius_arcsec=s["radius_arcsec"],
    )
    assert float(np.mean(m_sel)) > 0.88
    # vs the other-frequency survey: ~fade_fraction (the increment-1 blind spot, now modeled)
    m_oth, _, _ = dr20radio.crossmatch(
        s["ra_q"][radio_carton],
        s["dec_q"][radio_carton],
        s["ra_r2"],
        s["dec_r2"],
        radius_arcsec=s["radius_arcsec"],
    )
    assert float(np.mean(m_oth)) == pytest.approx(s["fade_fraction"], abs=0.12)
    assert float(np.mean(m_sel)) - float(np.mean(m_oth)) > 0.3


def test_parse_racs_csv():
    text = "ra,dec,peak_flux\n3.14,-41.5,2.5\nbad,row\n10.0,-50.0,1.1\n"
    ra, dec, flux = dr20radio.parse_racs_csv(text)
    assert ra.tolist() == [3.14, 10.0]
    assert dec.tolist() == [-41.5, -50.0]
    assert flux.tolist() == [2.5, 1.1]


def test_log_luminosity_known_value():
    # z=1, 1 mJy at 1.4 GHz, alpha=-0.7: d_L(Planck18)=6791.3 Mpc = 2.096e26 m ->
    # 4 pi d_L^2 = 5.52e53 m^2, x 1e-29 W/m^2/Hz x 2^-0.3 -> 4.48e24 -> log10 = 24.65
    lg = dr20radio.log_luminosity_whz(np.array([1.0]), np.array([1.0]), freq_ghz=1.4, alpha=-0.7)
    assert lg[0] == pytest.approx(24.65, abs=0.03)
    # brighter flux -> proportionally higher luminosity (pure scaling)
    lg10 = dr20radio.log_luminosity_whz(np.array([1.0]), np.array([10.0]), freq_ghz=1.4, alpha=-0.7)
    assert lg10[0] - lg[0] == pytest.approx(1.0, abs=1e-9)


def test_luminosity_matched_fractions_applies_common_limit():
    rng = np.random.default_rng(7)
    n = 2000
    z = rng.uniform(0.5, 2.5, n)
    matched = rng.random(n) < 0.5
    # counterpart fluxes just above this survey's limit: 1-3x s_lim
    s_lim = 1.0
    s = s_lim * rng.uniform(1.0, 3.0, n)
    bins = np.array([0.0, 1.0, 2.0, 3.0])
    # other survey shallower by 10x at same freq -> common limit cuts everything below
    # 10*s_lim: only counterparts with s >= 10 survive -> fraction ~ 0
    out_deep_other = dr20radio.luminosity_matched_fractions(
        z,
        matched,
        s,
        freq_ghz=1.4,
        s_lim_this_mjy=s_lim,
        s_lim_other_mjy=10.0,
        other_freq_ghz=1.4,
        bins=bins,
    )
    assert sum(out_deep_other["k"]) == 0
    # other survey deeper -> common limit is our own -> all matched count
    out_shallow_other = dr20radio.luminosity_matched_fractions(
        z,
        matched,
        s,
        freq_ghz=1.4,
        s_lim_this_mjy=s_lim,
        s_lim_other_mjy=0.1,
        other_freq_ghz=1.4,
        bins=bins,
    )
    assert sum(out_shallow_other["k"]) == int(matched.sum())


def test_k_correction_sign_and_alpha_monotonicity():
    """Pin the sign of the K-correction and the direction of the alpha sweep.

    The paper's north/south contrast is dominated by `alpha`, so a sign error here would
    silently invert the published conclusion: flat spectra would become the gap-maximising
    case instead of the gap-minimising one, and every number in Section 4.3 would still
    render. Nothing in the suite covered a non-fiducial alpha before 2026-08-12.
    """
    # A steep-spectrum source is BRIGHTER at low frequency, so extrapolating an 888 MHz
    # limit UP to 1.4 GHz must LOWER it, while extrapolating a 3 GHz limit DOWN must RAISE
    # it. At alpha = 0 both are unchanged: that is the assumption-free case the paper quotes.
    for alpha, racs_expected, vlass_expected in [
        (0.0, 3.0, 1.0),
        (-0.7, 3.0 * (1400 / 888.0) ** -0.7, 1.0 * (1400 / 3000.0) ** -0.7),
    ]:
        racs = 3.0 * (1400 / 888.0) ** alpha
        vlass = 1.0 * (1400 / 3000.0) ** alpha
        assert racs == pytest.approx(racs_expected)
        assert vlass == pytest.approx(vlass_expected)
    # The published crossover: RACS binds the common limit until alpha ~ -0.9, then VLASS does.
    assert 3.0 * (1400 / 888.0) ** -0.7 > 1.0 * (1400 / 3000.0) ** -0.7
    assert 3.0 * (1400 / 888.0) ** -1.0 < 1.0 * (1400 / 3000.0) ** -1.0

    # Same sample, same fluxes, sweeping alpha: the count above the common limit must move
    # monotonically, and must NOT be flat (a dropped alpha would make every variant equal).
    rng = np.random.default_rng(11)
    n = 4000
    z = rng.uniform(0.5, 2.5, n)
    matched = rng.random(n) < 0.5
    s = rng.uniform(1.0, 12.0, n)
    bins = np.array([0.0, 1.0, 2.0, 3.0])
    counts = [
        sum(
            dr20radio.luminosity_matched_fractions(
                z,
                matched,
                s,
                freq_ghz=3.0,
                s_lim_this_mjy=1.0,
                s_lim_other_mjy=3.0,
                other_freq_ghz=0.888,
                bins=bins,
                alpha=alpha,
            )["k"]
        )
        for alpha in dr20radio.ALPHA_SWEEP
    ]
    assert counts == sorted(counts), f"steepening alpha must not reduce the count: {counts}"
    assert counts[0] < counts[-1], f"alpha had no effect at all: {counts}"


def test_spectral_index_round_trip_and_sign():
    """alpha is defined S ~ nu^alpha, so a source fainter at high frequency is NEGATIVE."""
    s_lo = np.array([10.0, 10.0, 10.0])
    # flat, steep (half the flux per e-fold up), and inverted
    for alpha in (0.0, -0.7, -1.5, 0.5):
        s_hi = s_lo * (3.0 / 0.8875) ** alpha
        got = dr20radio.spectral_index(s_hi, s_lo, freq_hi_ghz=3.0, freq_lo_ghz=0.8875)
        assert got == pytest.approx(alpha)
    # a steep spectrum really is fainter at the higher frequency
    assert (s_lo * (3.0 / 0.8875) ** -0.7 < s_lo).all()
    # non-positive fluxes drop out rather than raising or returning +-inf
    bad = dr20radio.spectral_index(
        np.array([1.0, 0.0, -1.0]), np.array([1.0, 1.0, 1.0]), freq_hi_ghz=3.0, freq_lo_ghz=0.9
    )
    assert np.isfinite(bad[0]) and np.isnan(bad[1]) and np.isnan(bad[2])


def test_alpha_complete_limit_removes_the_truncation():
    """Above the returned flux, no source with alpha >= alpha_min can be missing from the
    joint-detection sample. This is the whole basis of the unbiased sub-sample."""
    s_lim_hi, f_hi, f_lo, a_min = 1.0, 3.0, 0.8875, -1.5
    lim = dr20radio.alpha_complete_limit_mjy(
        s_lim_hi_mjy=s_lim_hi, freq_hi_ghz=f_hi, freq_lo_ghz=f_lo, alpha_min=a_min
    )
    # a source AT the limit with the steepest allowed index still clears the high-freq limit
    assert lim * (f_hi / f_lo) ** a_min == pytest.approx(s_lim_hi)
    # ...and anything steeper than alpha_min would not, which is why the bound is stated
    assert lim * (f_hi / f_lo) ** (a_min - 0.3) < s_lim_hi
    # below the limit, a steep source is lost -> that is the flat bias the paper reports
    assert 0.5 * lim * (f_hi / f_lo) ** a_min < s_lim_hi


def test_measured_alpha_constants_match_committed_measurement():
    """`ALPHA_MEASURED` is hard-coded so the census runs can sweep it, which means it can
    drift away from the measurement it came from. The paper quotes both, so a drift would
    print a measured index next to a contrast evaluated at a different one."""
    import json
    from pathlib import Path

    results = Path("results/dr20radio_alpha.json")
    # NOT a skip. This file is tracked committed evidence, and every alpha macro in the paper
    # reads from it; a skip here would disable the guard precisely in a fresh CI checkout,
    # which is where the hard-coded ALPHA_MEASURED most needs checking.
    assert results.exists(), "results/dr20radio_alpha.json is committed evidence and must exist"
    fc = json.loads(results.read_text())["flux_complete"]
    assert dr20radio.ALPHA_MEASURED == pytest.approx(fc["median"], abs=5e-5)
    assert dr20radio.ALPHA_MEASURED_SE == pytest.approx(fc["median_boot_se"], abs=5e-5)
    # the sweep the census runs must actually contain the three measured points
    for a in dr20radio.ALPHA_MEASURED_SWEEP:
        key = f"{a:g}"
        north = json.loads(Path("results/dr20radio_north.json").read_text())
        assert key in north["luminosity_matched_alpha"], f"census never evaluated alpha={key}"


def test_joint_detection_median_is_biased_flat():
    """The committed measurement must show the truncation bias in the direction the method
    section claims -- if it ever came out the other way, the flux-complete cut is wrong."""
    import json
    from pathlib import Path

    results = Path("results/dr20radio_alpha.json")
    # NOT a skip. This file is tracked committed evidence, and every alpha macro in the paper
    # reads from it; a skip here would disable the guard precisely in a fresh CI checkout,
    # which is where the hard-coded ALPHA_MEASURED most needs checking.
    assert results.exists(), "results/dr20radio_alpha.json is committed evidence and must exist"
    m = json.loads(results.read_text())
    assert m["joint_detection"]["median"] > m["flux_complete"]["median"], (
        "the joint-detection sample should be biased FLAT relative to the flux-complete one"
    )
    assert m["flux_complete"]["n"] < m["joint_detection"]["n"]


def test_per_source_alpha_is_a_bias_not_a_variance():
    """Per-source indices must change the ANSWER, and the realization spread must be small.

    If the spread across realizations were the headline, any breadth of index distribution
    would look harmless, because a fraction over N sources averages the draw away as
    1/sqrt(N). That is the shape of the `rmstructure` error -- an uncertainty estimated by a
    procedure that cannot see the effect it is supposed to bound.
    """
    rng = np.random.default_rng(3)
    n = 20000
    z = rng.uniform(0.5, 2.5, n)
    matched = rng.random(n) < 0.5
    s = rng.uniform(1.0, 12.0, n)
    bins = np.array([0.0, 1.0, 2.0, 3.0])
    kw = dict(
        freq_ghz=3.0,
        s_lim_this_mjy=1.0,
        s_lim_other_mjy=3.0,
        other_freq_ghz=0.888,
        bins=bins,
    )
    broad = rng.normal(-0.7, 0.9, 4000)
    out = dr20radio.luminosity_matched_per_source_alpha(
        z, matched, s, alpha_samples=broad, n_real=8, seed=1, **kw
    )
    single = dr20radio.luminosity_matched_fractions(z, matched, s, alpha=-0.7, **kw)
    single_frac = sum(single["k"]) / sum(single["n"])

    # the realization spread is small -- this is the number that would have been misleading
    assert out["fraction_realization_sd"] < 0.01
    # ...while the scatter genuinely shifts the answer away from the single-index result
    assert abs(out["fraction_mean"] - single_frac) > out["fraction_realization_sd"]
    # a degenerate "distribution" of one repeated value must reproduce the single-index case
    degenerate = dr20radio.luminosity_matched_per_source_alpha(
        z, matched, s, alpha_samples=np.full(500, -0.7), n_real=2, seed=1, **kw
    )
    assert degenerate["fraction_mean"] == pytest.approx(single_frac, abs=1e-12)
    assert degenerate["fraction_realization_sd"] == pytest.approx(0.0, abs=1e-12)


def test_alpha_systematics_are_recorded_and_the_floor_test_can_fail():
    """The three checks that bound the measured index must all be committed, and the
    completeness-floor one must be reported in a form that could have contradicted the paper.

    A completeness cut at alpha >= -1.5 is an assumption; the sample's own p16 is steeper. If
    the median steepens as the floor is lowered, the headline value is an upper bound on
    flatness, not an unbiased estimate -- which is what the committed run shows. A test that
    only asserted the cut's algebra would never have surfaced that.
    """
    import json
    from pathlib import Path

    m = json.loads(Path("results/dr20radio_alpha.json").read_text())
    floors = m["completeness_floor_sensitivity"]
    assert set(floors) == {"-1.5", "-2", "-2.5"}
    # a stricter floor must demand a brighter cut and retain fewer sources
    cuts = [floors[k]["flux_cut_mjy"] for k in ("-1.5", "-2", "-2.5")]
    ns = [floors[k]["n"] for k in ("-1.5", "-2", "-2.5")]
    assert cuts == sorted(cuts) and ns == sorted(ns, reverse=True)
    # the paper states the median steepens monotonically -- if that ever reverses, the
    # "upper bound on flatness" framing is wrong and must be rewritten
    meds = [floors[k]["median"] for k in ("-1.5", "-2", "-2.5")]
    assert meds == sorted(meds, reverse=True), f"floor sensitivity is no longer monotonic: {meds}"

    # both epochs and all flux bins present, with the max-of-epochs value bracketed
    assert set(m["per_epoch"]) == {"E2", "E3"}
    ep = [m["per_epoch"][k]["median"] for k in ("E2", "E3")]
    assert min(ep) <= m["flux_complete"]["median"] <= max(ep)
    assert len(m["flux_bins"]) == 3
    assert all(b["n"] > 0 and b["median"] is not None for b in m["flux_bins"])


def test_vlass_conservative_variant_can_move_the_ratio():
    """The RACS-side conservative variant provably cannot move the north/south ratio; the
    VLASS-side one must be able to, or it is the same vacuous check on a new axis."""
    import json
    from pathlib import Path

    n = json.loads(Path("results/dr20radio_north.json").read_text())
    s = json.loads(Path("results/dr20radio_south.json").read_text())["deep_south"]
    tot = lambda b: sum(b["k"]) / sum(b["n"])  # noqa: E731
    key = f"{dr20radio.ALPHA_MEASURED:g}"
    gap = tot(n["luminosity_matched_alpha"][key]) - tot(s["luminosity_matched_alpha"][key])
    for v in dr20radio.VLASS_S_LIM_CONSERVATIVE_MJY:
        vk = f"{v:g}"
        assert vk in n["luminosity_matched_vlass_conservative"]
        assert vk in s["luminosity_matched_vlass_conservative"]
    worst = f"{max(dr20radio.VLASS_S_LIM_CONSERVATIVE_MJY):g}"
    gap_c = tot(n["luminosity_matched_vlass_conservative"][worst]) - tot(
        s["luminosity_matched_vlass_conservative"][worst]
    )
    # it must actually change the answer -- a check that cannot fail is not a check
    assert abs(gap_c - gap) > 0.1 * gap, "the VLASS-limit variant barely moves the gap"
