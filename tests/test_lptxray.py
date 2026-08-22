"""Tests for the LPT / white-dwarf-candidate X-ray cross-match (plan 93)."""

from __future__ import annotations

import math

import pytest

from jansky_research import lptxray as lx


def _row(name: str, ra: float, dec: float, **extra: object) -> dict:
    return {"name": name, "ra": ra, "dec": dec, **extra}


class TestSeparation:
    def test_zero_separation(self) -> None:
        assert lx.angular_sep_arcsec(10.0, -20.0, 10.0, -20.0) == pytest.approx(0.0, abs=1e-6)

    def test_one_arcsec_in_declination(self) -> None:
        sep = lx.angular_sep_arcsec(10.0, -20.0, 10.0, -20.0 + 1 / 3600)
        assert sep == pytest.approx(1.0, abs=1e-3)

    def test_ra_separation_shrinks_with_declination(self) -> None:
        """One degree of RA is cos(dec) degrees on the sky -- the cos term must be present."""
        at_equator = lx.angular_sep_arcsec(10.0, 0.0, 11.0, 0.0)
        at_sixty = lx.angular_sep_arcsec(10.0, 60.0, 11.0, 60.0)
        assert at_sixty == pytest.approx(at_equator * math.cos(math.radians(60.0)), rel=1e-4)


class TestAssociate:
    def test_matches_inside_radius(self) -> None:
        rows = [_row("src", 10.0 + 5 / 3600, -20.0)]
        got = lx.associate(10.0, -20.0, rows, 15.0)
        assert got["matched"] is True
        assert got["n_within"] == 1
        assert got["sep_arcsec"] == pytest.approx(4.7, abs=0.5)

    def test_rejects_outside_radius(self) -> None:
        rows = [_row("src", 10.0 + 40 / 3600, -20.0)]
        got = lx.associate(10.0, -20.0, rows, 15.0)
        assert got["matched"] is False
        assert got["n_within"] == 0
        assert got["nearest"] == "src"

    def test_empty_cone(self) -> None:
        got = lx.associate(10.0, -20.0, [], 15.0)
        assert got == {"matched": False, "n_within": 0, "sep_arcsec": None, "nearest": None}


class TestChanceRate:
    def test_density_excludes_the_inner_region(self) -> None:
        """A source on top of the target must not inflate the *field* density."""
        rows = [_row("on-target", 10.0, -20.0)]
        d = lx.local_density_per_sqarcmin(10.0, -20.0, rows, cone_arcmin=10.0, inner_arcmin=2.0)
        assert d == 0.0

    def test_density_counts_the_annulus(self) -> None:
        rows = [_row("far", 10.0, -20.0 + 5 / 60)]
        d = lx.local_density_per_sqarcmin(10.0, -20.0, rows, cone_arcmin=10.0, inner_arcmin=2.0)
        assert d == pytest.approx(1.0 / (math.pi * (100 - 4)), rel=1e-6)

    def test_chance_expected_scales_as_radius_squared(self) -> None:
        a = lx.chance_expected(0.01, 15.0)
        b = lx.chance_expected(0.01, 30.0)
        assert b == pytest.approx(4.0 * a, rel=1e-9)

    def test_shift_trials_find_nothing_in_an_empty_field(self) -> None:
        got = lx.shift_trial_rate(10.0, -20.0, [], radius_arcsec=45.0, cone_arcmin=10.0)
        assert got["n_trials"] == 24
        assert got["n_hits"] == 0
        assert got["rate"] == 0.0

    def test_shift_trials_fire_in_a_crowded_field(self) -> None:
        """A field packed with sources must produce a non-zero measured chance rate."""
        rows = [
            _row(f"s{i}{j}", 10.0 + i / 60.0 / math.cos(math.radians(-20.0)), -20.0 + j / 60.0)
            for i in range(-9, 10)
            for j in range(-9, 10)
        ]
        got = lx.shift_trial_rate(10.0, -20.0, rows, radius_arcsec=45.0, cone_arcmin=10.0)
        assert got["n_hits"] > 0
        assert 0.0 < got["rate"] <= 1.0

    def test_shift_discs_stay_inside_the_cone(self) -> None:
        """A trial that ran off the edge of the cached cone would see false completeness."""
        got = lx.shift_trial_rate(
            10.0, -20.0, [], radius_arcsec=45.0, cone_arcmin=4.0, shifts_arcmin=(3.0, 5.0, 7.0)
        )
        assert got["shifts_arcmin"] == [3.0]
        assert got["n_trials"] == 8


class TestRecallGuard:
    def test_recall_of_a_working_crossmatch(self) -> None:
        recs = [{"name": f"c{i}", "lit_xray_detected": True, "matched": True} for i in range(16)]
        got = lx.catalogue_recall(recs)
        assert got["n_known_detections"] == 16
        assert got["recall"] == 1.0
        assert got["usable"] is True
        assert got["missed"] == []

    def test_low_recall_is_marked_unusable(self) -> None:
        """The LPT leg: 1 of 3 known detections recovered. A null here means nothing."""
        recs = [
            {"name": "J1832-0911", "lit_xray_detected": True, "matched": False},
            {"name": "J1448-6856", "lit_xray_detected": True, "matched": False},
            {"name": "J1745-5051", "lit_xray_detected": True, "matched": True},
            {"name": "quiet", "lit_xray_detected": False, "matched": False},
        ]
        got = lx.catalogue_recall(recs)
        assert got["n_known_detections"] == 3
        assert got["recall"] == pytest.approx(1 / 3, abs=1e-4)
        assert got["usable"] is False
        assert sorted(got["missed"]) == ["J1448-6856", "J1832-0911"]

    def test_no_known_detections_is_not_usable(self) -> None:
        got = lx.catalogue_recall([{"name": "x", "lit_xray_detected": False, "matched": False}])
        assert got["n_known_detections"] == 0
        assert got["usable"] is False
        assert math.isnan(got["recall"])


class TestFractionComparison:
    def test_reports_both_difference_and_ratio(self) -> None:
        got = lx.fraction_comparison(15, 21, 2, 35)
        assert got["frac_a"] == pytest.approx(15 / 21, abs=1e-5)
        assert got["frac_b"] == pytest.approx(2 / 35, abs=1e-5)
        assert got["difference_pp"] == pytest.approx(100 * (15 / 21 - 2 / 35), abs=1e-3)
        assert got["ratio"] == pytest.approx((15 / 21) / (2 / 35), abs=1e-4)

    def test_a_common_cut_moves_the_difference_but_not_the_ratio(self) -> None:
        """The `dr20radio` lesson, pinned: halving both arms leaves the ratio alone."""
        shallow = lx.fraction_comparison(20, 40, 10, 40)
        deep = lx.fraction_comparison(10, 40, 5, 40)
        assert deep["ratio"] == pytest.approx(shallow["ratio"], rel=1e-9)
        assert abs(deep["difference_pp"]) < abs(shallow["difference_pp"])

    def test_zero_denominator_ratio_is_none_not_infinite(self) -> None:
        got = lx.fraction_comparison(5, 10, 0, 10)
        assert got["ratio"] is None
        assert got["ci_b"][1] > 0.0

    def test_empty_sample_raises(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            lx.fraction_comparison(1, 0, 1, 5)


class TestSummarizePosition:
    def _cats(self, rows: list[dict]) -> dict:
        return {c: {"n": len(rows), "rows": list(rows)} for c in lx.CATALOGS}

    def test_detected_position_carries_flux_and_provenance(self) -> None:
        rows = [
            _row(
                "2RXS J1",
                10.0 + 3 / 3600,
                -20.0,
                count_rate=0.05,
                ep_flux=1e-13,
                b1_flux=2e-13,
                b_flux_ap=3e-13,
            )
        ]
        got = lx.summarize_position(
            {
                "name": "cand",
                "sample": "wdcand",
                "type": "polar",
                "ra": 10.0,
                "dec": -20.0,
                "lit_xray": "ROSAT",
            },
            self._cats(rows),
            cone_arcmin=10.0,
        )
        assert got["matched"] is True
        assert got["lit_xray_detected"] is True
        assert got["catalogs"]["rass2rxs"]["flux"] == 0.05
        assert got["catalogs"]["rass2rxs"]["flux_units"] == "ct/s"
        assert got["catalogs"]["xmmssc"]["flux_band"] == "0.2-12 keV"

    def test_blank_literature_flag_is_not_a_detection(self) -> None:
        got = lx.summarize_position(
            {
                "name": "q",
                "sample": "lpt",
                "type": "LPT",
                "ra": 10.0,
                "dec": -20.0,
                "lit_xray": "no",
            },
            self._cats([]),
            cone_arcmin=10.0,
        )
        assert got["matched"] is False
        assert got["lit_xray_detected"] is False

    def test_failed_query_is_recorded_not_silently_zero(self) -> None:
        """A catalogue that errored must not be counted as a non-detection."""
        got = lx.summarize_position(
            {"name": "q", "sample": "lpt", "type": "LPT", "ra": 10.0, "dec": -20.0, "lit_xray": ""},
            {"xmmssc": {"error": "service down"}},
            cone_arcmin=10.0,
        )
        assert got["catalogs"]["xmmssc"]["error"] == "service down"
        assert "matched" not in got["catalogs"]["xmmssc"]
        assert got["catalogs"]["rass2rxs"]["error"] == "not queried"


class TestClassifyCoverage:
    def test_on_axis_pointing_is_targeted(self) -> None:
        obs = [{"obsid": "1", "name": "ASKAP J1745", "ra": 266.287, "dec": -50.864}]
        got = lx.classify_coverage(266.287, -50.864, obs)
        assert got["n_targeted"] == 1
        assert got["n_serendipitous"] == 0
        assert got["targeted"][0]["offset_arcmin"] == pytest.approx(0.0, abs=1e-3)

    def test_nearby_pointing_at_another_source_is_serendipitous(self) -> None:
        """The ASKAP J1935+2148 / SGR 1935+2154 case: same RA digits, different source."""
        obs = [{"obsid": "9", "name": "SGR 1935+2154", "ra": 293.732, "dec": 21.8967}]
        got = lx.classify_coverage(293.7325, 21.8033, obs)
        assert got["n_targeted"] == 0
        assert got["n_serendipitous"] == 1
        assert got["serendipitous"][0]["offset_arcmin"] > lx.ON_AXIS_ARCMIN

    def test_rows_without_coordinates_are_skipped_not_counted(self) -> None:
        got = lx.classify_coverage(10.0, -20.0, [{"obsid": "x", "name": "no position"}])
        assert got["n_targeted"] == 0
        assert got["n_serendipitous"] == 0

    def test_unparseable_coordinates_do_not_raise(self) -> None:
        got = lx.classify_coverage(10.0, -20.0, [{"obsid": "x", "ra": "", "dec": "n/a"}])
        assert got["n_targeted"] == 0 and got["n_serendipitous"] == 0
