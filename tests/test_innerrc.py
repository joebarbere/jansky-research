"""Tests for jansky_research.innerrc — Sofue & Kohno inner-RC replication. No network."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from jansky_research import innerrc

TABLES = Path(__file__).parent / "data" / "sofue2025"


def test_parse_paper_tables_vendored():
    t = innerrc.parse_paper_tables(str(TABLES))
    assert t["inner"]["R_pc"].size > 100
    assert t["unified"]["R_pc"].size > 50
    # sorted, physical ranges
    assert np.all(np.diff(t["inner"]["R_pc"]) >= 0)
    assert t["inner"]["R_pc"][0] == pytest.approx(50.0)
    assert 50 < np.nanmedian(t["inner"]["V_kms"]) < 350
    assert np.all(t["unified"]["dV_kms"] >= 0)


def test_gaussian_tvm_recovers_known_terminal_velocity():
    vel = np.arange(-50.0, 160.0, 0.5)
    for seed, vterm in [(0, 95.0), (1, 120.0), (2, 70.0)]:
        spec = innerrc.synthetic_spectrum(vel, vterm, seed=seed)
        est = innerrc.gaussian_tvm(vel, spec)
        assert est == pytest.approx(vterm, abs=6.0)


def test_threshold_tvm_reads_high_of_gaussian_tvm():
    # The documented hi-slice bias: the envelope crossing sits above the last component centre.
    vel = np.arange(-50.0, 160.0, 0.5)
    diffs = []
    for seed in range(5):
        spec = innerrc.synthetic_spectrum(vel, 100.0, seed=seed)
        diffs.append(innerrc.threshold_tvm(vel, spec) - innerrc.gaussian_tvm(vel, spec))
    assert np.mean(diffs) > 3.0  # threshold reads several km/s high on the same spectra
    assert innerrc.threshold_tvm(vel, np.zeros_like(vel)) != innerrc.threshold_tvm(vel, vel * 0 + 3)


def test_rc_from_terminal_and_dispersion_correction():
    r, v = innerrc.rc_from_terminal(np.array([30.0]), np.array([100.0]), sigma_v_kms=15.0)
    assert r[0] == pytest.approx(innerrc.R0_PC * 0.5)
    assert v[0] == pytest.approx((100.0 - 15.0) + innerrc.V0_KMS * 0.5)
    # fourth quadrant: blueshifted terminal velocity, same |R|, same V
    r4, v4 = innerrc.rc_from_terminal(np.array([-30.0]), np.array([-100.0]), sigma_v_kms=15.0)
    assert r4[0] == pytest.approx(r[0])
    assert v4[0] == pytest.approx(v[0])


def test_rotation_curve_weighted_recovers_flat_curve():
    rng = np.random.default_rng(0)
    r = rng.uniform(500, 8000, 500)
    v = 230.0 + rng.normal(0, 5.0, r.size)
    grid, vm, dv = innerrc.rotation_curve_weighted(r, v, dr_pc=250.0, half_width_pc=100.0)
    ok = np.isfinite(vm)
    assert np.nanmean(np.abs(vm[ok] - 230.0)) < 3.0
    assert np.nanmedian(dv[ok]) < 8.0


def test_decompose_recovers_injected_model():
    r = np.linspace(50, 15000, 300)
    truth = dict(
        v_bulge=400.0, a_bulge=330.0, v_disc=320.0, a_disc=5600.0, v_halo=150.0, h_halo=22000.0
    )
    v = innerrc.rc_model(r, *truth.values())
    fit = innerrc.decompose_rc(r, v)
    for k, val in truth.items():
        assert fit[k] == pytest.approx(val, rel=0.15), k
    assert fit["rms_kms"] < 1.0
    assert fit["rho_dm_gev"] > 0


def test_rho_dm_local_matches_hand_calculation():
    # rho0 = V^2/(4 pi G h^2); x = R0/h; NFW rho = rho0/(x(1+x)^2)
    v_h, h_pc = 150.0, 22000.0
    h_kpc = h_pc / 1000
    rho0 = v_h**2 / (4 * np.pi * innerrc.G_KPC * h_kpc**2)
    x = innerrc.R0_PC / h_pc
    expect = rho0 / (x * (1 + x) ** 2) / 1e9 * innerrc.GEV_PER_MSUN_PC3
    assert innerrc.rho_dm_local_gev(v_h, h_pc) == pytest.approx(expect)
    # sanity: in the ballpark of published values (0.05-1 GeV/cm^3)
    assert 0.01 < innerrc.rho_dm_local_gev(v_h, h_pc) < 2.0


def test_ew_asymmetry_fit_recovers_injected_sinusoid():
    r = np.linspace(500, 8000, 200)
    truth = (45.0, 3500.0, 4000.0, 4400.0)
    dv = truth[0] * np.exp(-r / truth[1]) * np.sin(2 * np.pi * (r - truth[2]) / truth[3])
    rng = np.random.default_rng(1)
    fit = innerrc.ew_asymmetry_fit(r, dv + rng.normal(0, 1.0, r.size))
    assert fit["amp_kms"] == pytest.approx(45.0, rel=0.25)
    assert fit["period_pc"] == pytest.approx(4400.0, rel=0.15)
    assert fit["rms_kms"] < 2.0


def test_paper_table1_rho_dm_reproduces_through_convention():
    # Their eq. 24-25 convention: rho0 = V_h^2/(G h^2); converting their V_h into ours by
    # sqrt(4 pi) must reproduce their published 0.107 GeV/cm^3 through our formula.
    v_ours = 64.4 * np.sqrt(4 * np.pi)
    assert innerrc.rho_dm_local_gev(v_ours, 22379.1) == pytest.approx(0.107, abs=0.005)


def test_lv_from_cube_and_tvm_spectrum_round_trip():
    # tiny synthetic (v, b, l) cube: one drifting terminal edge per longitude column
    vel = np.arange(-200.0, 200.0, 2.0)
    lat = np.arange(-5.0, 5.0, 0.5)
    lon = np.arange(20.0, 24.0, 1.0)
    cube = np.zeros((vel.size, lat.size, lon.size))
    vterms = [80.0, 90.0, 100.0, 110.0]
    for j, vt in enumerate(vterms):
        spec = innerrc.synthetic_spectrum(vel, vt, seed=j, peak_k=20.0)
        cube[:, :, j] = spec[:, None]  # same spectrum at every latitude
    glon, v, t_lv = innerrc.lv_from_cube(cube, 20.0, 1.0, -5.0, 0.5, -200.0, 2.0, b_max_deg=3.0)
    assert glon[0] == pytest.approx(20.0)
    assert t_lv.shape == (vel.size, lon.size)
    for j, vt in enumerate(vterms):
        est = innerrc.tvm_spectrum(v, t_lv[:, j], sign=+1, method="gaussian")
        assert est == pytest.approx(vt, abs=8.0)
        thr = innerrc.tvm_spectrum(v, t_lv[:, j], sign=+1, method="threshold")
        assert thr >= est - 2.0  # threshold sits at/above the outermost component centre
    # fourth quadrant: mirror the spectrum to negative velocities
    neg = innerrc.synthetic_spectrum(vel, 95.0, seed=9, peak_k=20.0)[::-1]
    est4 = innerrc.tvm_spectrum(v, neg, sign=-1, method="gaussian")
    assert est4 == pytest.approx(-95.0, abs=8.0)
    # empty spectrum -> nan
    assert np.isnan(innerrc.tvm_spectrum(v, np.zeros_like(v), sign=+1))


def test_run_anchor_offline_reproduces_paper_scale_rho_dm(tmp_path):
    m = innerrc.run_anchor(str(tmp_path), table_dir=str(TABLES))
    assert m["n_unified_rows"] > 50
    fit = m["anchor_fit"]
    # the decomposition must land in the paper's regime: a low (<0.3 consensus) halo-only DMD
    assert 0.02 < fit["rho_dm_gev"] < 0.5
    # n_converged counts variants with every parameter INTERIOR, not merely variants where
    # curve_fit returned. Most of this scan rails (excising R < 2 kpc leaves the bulge
    # unconstrained), so the bar is deliberately low — the point is that the quoted range is
    # built only from fits the data actually determined.
    sens = m["sensitivity"]
    assert sens["n_converged"] >= 2
    assert sens["n_converged"] <= sens["n_fitted"]
    assert sens["rho_dm_max_gev"] >= sens["rho_dm_min_gev"]
    for key, railed in sens["railed_variants"].items():
        assert railed, f"{key} listed as railed with no parameter named"
    assert fit["railed_params"] == []
    assert (tmp_path / "results" / "innerrc_anchor.json").exists()


def test_bound_contact_flags_railed_parameters():
    """A `curve_fit` that does not raise is not the same as a converged fit.

    Before 2026-08-12 the sensitivity scan reported `n_converged: 8` and a rho_DM range whose
    maximum came from a variant with v_bulge at exactly 800.0 km/s — the upper bound. SciPy
    reports success for a solution glued to a wall, so nothing caught it; the number then
    carried the paper's "fully compatible with the consensus density" claim.
    """
    from jansky_research.innerrc import FIT_LOWER, FIT_PARAM_NAMES, FIT_UPPER, bound_contact

    interior = dict(
        zip(
            FIT_PARAM_NAMES,
            [(lo + hi) / 2 for lo, hi in zip(FIT_LOWER, FIT_UPPER, strict=True)],
            strict=True,
        )
    )
    assert bound_contact(interior) == []

    for i, name in enumerate(FIT_PARAM_NAMES):
        at_upper = dict(interior, **{name: FIT_UPPER[i]})
        at_lower = dict(interior, **{name: FIT_LOWER[i]})
        assert bound_contact(at_upper) == [name]
        assert bound_contact(at_lower) == [name]
        # just inside the tolerance band still counts as railed
        span = FIT_UPPER[i] - FIT_LOWER[i]
        assert bound_contact(dict(interior, **{name: FIT_UPPER[i] - 0.005 * span})) == [name]
        # clearly interior does not
        assert bound_contact(dict(interior, **{name: FIT_UPPER[i] - 0.2 * span})) == []


def test_committed_anchor_fit_is_interior_and_scan_excludes_railed():
    """The published refit must not rest on a bound, and the quoted range must exclude ones
    that do. This is the assertion whose absence let a boundary artifact reach an abstract."""
    import json
    from pathlib import Path

    results = Path("results/innerrc_anchor.json")
    if not results.exists():  # pragma: no cover - real-results file absent in a bare checkout
        pytest.skip("committed anchor results not present")
    a = json.loads(results.read_text())

    assert a["anchor_fit"]["railed_params"] == [], (
        "the primary refit has a parameter on a bound; its rho_DM reports the bound, not the data"
    )
    quoted = [
        v["rho_dm_gev"]
        for v in a["sensitivity_variants"].values()
        if "rho_dm_gev" in v and not v["railed_params"]
    ]
    assert quoted, "no interior variants remain"
    assert a["sensitivity"]["rho_dm_min_gev"] == pytest.approx(min(quoted))
    assert a["sensitivity"]["rho_dm_max_gev"] == pytest.approx(max(quoted))
    assert a["sensitivity"]["n_converged"] == len(quoted)
    # every variant must carry the halo parameters its rho_DM is computed from
    for key, v in a["sensitivity_variants"].items():
        if "rho_dm_gev" in v:
            assert "v_halo" in v and "h_halo" in v, f"{key} drops the parameters behind its rho"
