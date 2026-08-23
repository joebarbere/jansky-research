# --------------------------------------------------------------------------------------
# preserve_live_macros — the fix for the cross-run clobber that blanked four abstracts
# --------------------------------------------------------------------------------------


def _macros(path, **pairs):
    path.write_text("\n".join(rf"\newcommand{{\{k}}}{{{v}}}" for k, v in pairs.items()) + "\n")
    return path


def test_a_placeholder_never_overwrites_a_real_value(tmp_path):
    """The bug: an offline run blanked the real run's macros because its metrics dict had no
    such keys, and vice versa. A run may only ADD information."""
    from jansky_research.report import preserve_live_macros

    existing = _macros(tmp_path / "m.tex", tiiSynPurity="1.0", tiiRealDetRate="--")
    new = "\n".join([r"\newcommand{\tiiSynPurity}{--}", r"\newcommand{\tiiRealDetRate}{0.432}"])
    out = preserve_live_macros(new + "\n", existing)
    assert r"\newcommand{\tiiSynPurity}{1.0}" in out  # preserved
    assert r"\newcommand{\tiiRealDetRate}{0.432}" in out  # written


def test_a_real_value_still_overwrites_a_real_value(tmp_path):
    """Only blanking is prevented — a genuine recomputation must still take effect."""
    from jansky_research.report import preserve_live_macros

    existing = _macros(tmp_path / "m.tex", tiiSynPurity="0.90")
    out = preserve_live_macros("\n".join([r"\newcommand{\tiiSynPurity}{0.97}"]) + "\n", existing)
    assert r"\newcommand{\tiiSynPurity}{0.97}" in out


def test_no_existing_file_is_not_an_error(tmp_path):
    from jansky_research.report import preserve_live_macros

    new = "\n".join([r"\newcommand{\a}{1}"]) + "\n"
    assert preserve_live_macros(new, tmp_path / "absent.tex") == new


def test_merging_cannot_rescue_a_name_collision(tmp_path):
    """Why namespacing was also required, and the reason this helper is not sufficient alone.

    When one macro name means different things in the two run modes, BOTH runs write a real
    value for it and there is nothing left to arbitrate on. `\\tiiNEvents` meant 768 real
    observing days in the paper's prose and 48 synthetic events in the offline run; an offline
    rebuild silently turned "768 days, zero failures" into "48 days, zero failures". The
    helper below happily lets that through, which is correct behaviour for a *recomputation*
    and catastrophic for a *collision* — so every mode-dependent macro is now namespaced.
    """
    from jansky_research.report import preserve_live_macros

    existing = _macros(tmp_path / "m.tex", tiiNEvents="768")
    out = preserve_live_macros("\n".join([r"\newcommand{\tiiNEvents}{48}"]) + "\n", existing)
    assert r"\newcommand{\tiiNEvents}{48}" in out  # not rescued — hence the namespacing


# --------------------------------------------------------------------------------------
# preserve_live_results — the same rule for the results JSON, which had the same hole
# --------------------------------------------------------------------------------------


def _results(path, payload):
    import json

    path.write_text(json.dumps(payload) + "\n")
    return path


def test_results_no_existing_file_writes_the_new_payload(tmp_path):
    from jansky_research.report import preserve_live_results

    assert preserve_live_results({"is_real": False, "a": 1}, tmp_path / "m.json") == {
        "is_real": False,
        "a": 1,
    }


def test_results_real_replaces_synthetic(tmp_path):
    """A real run must be able to overwrite a synthetic placeholder."""
    from jansky_research.report import preserve_live_results

    p = _results(tmp_path / "m.json", {"is_real": False, "n": 48})
    got = preserve_live_results({"is_real": True, "n": 768}, p)
    assert got["n"] == 768 and got["is_real"] is True


def test_results_synthetic_does_not_replace_real(tmp_path):
    """The typeii clobber: run('.', offline=True) deleted a real census and flipped is_real."""
    from jansky_research.report import preserve_live_results

    p = _results(tmp_path / "m.json", {"is_real": True, "n_events": 768, "rows": [1, 2, 3]})
    got = preserve_live_results({"is_real": False, "n_events": 48}, p)
    assert got["n_events"] == 768
    assert got["is_real"] is True
    assert got["rows"] == [1, 2, 3]


def test_results_synthetic_over_synthetic_is_a_plain_rebuild(tmp_path):
    from jansky_research.report import preserve_live_results

    p = _results(tmp_path / "m.json", {"is_real": False, "n": 1})
    assert preserve_live_results({"is_real": False, "n": 2}, p)["n"] == 2


def test_results_partial_real_rerun_retains_the_other_leg(tmp_path):
    """The torchfdmt case: a CPU-only re-run must not drop a real GPU benchmark column."""
    from jansky_research.report import preserve_live_results

    p = _results(tmp_path / "m.json", {"is_real": True, "cpu_s": 44.12, "gpu_s": 1.5})
    got = preserve_live_results({"is_real": True, "cpu_s": 43.0, "device": "cpu"}, p)
    assert got["cpu_s"] == 43.0
    assert got["gpu_s"] == 1.5


def test_results_a_retained_key_is_recorded_as_cross_run(tmp_path):
    """Merging must not let a spliced row look like the product of one invocation."""
    from jansky_research.report import RESULTS_MERGE_KEY, preserve_live_results

    p = _results(tmp_path / "m.json", {"is_real": True, "cpu_s": 44.12, "gpu_s": 1.5})
    got = preserve_live_results({"is_real": True, "cpu_s": 43.0}, p)
    assert got[RESULTS_MERGE_KEY]["retained_from_previous_run"] == ["gpu_s"]


def test_results_no_merge_block_when_nothing_was_retained(tmp_path):
    from jansky_research.report import RESULTS_MERGE_KEY, preserve_live_results

    p = _results(tmp_path / "m.json", {"is_real": True, "a": 1})
    assert RESULTS_MERGE_KEY not in preserve_live_results({"is_real": True, "a": 2}, p)


def test_results_unreadable_existing_file_does_not_block_a_write(tmp_path):
    from jansky_research.report import preserve_live_results

    p = tmp_path / "m.json"
    p.write_text("{not json")
    assert preserve_live_results({"is_real": True, "a": 1}, p)["a"] == 1


def test_write_results_refuses_to_downgrade_on_disk(tmp_path):
    import json

    from jansky_research.report import write_results

    p = tmp_path / "sub" / "m.json"
    write_results({"is_real": True, "a": 1}, p)
    write_results({"is_real": False, "a": 999}, p)
    assert json.loads(p.read_text())["a"] == 1, "a synthetic rebuild overwrote real evidence"


def test_results_a_mixed_source_counts_as_real(tmp_path):
    """stokesv_discovery's source names both legs; the real census is in there, and treating
    any mention of "synthetic" as synthetic left it overwritable by its own offline rebuild."""
    from jansky_research.report import preserve_live_results

    p = _results(
        tmp_path / "m.json",
        {"source": "synthetic recover-a-known + real RACS-mid epoch pair", "n": 54},
    )
    got = preserve_live_results({"source": "synthetic epoch pair", "n": 400}, p)
    assert got["n"] == 54


def test_results_a_purely_synthetic_source_is_overwritable(tmp_path):
    """ecallisto_census is synthetic by design; its offline rebuild must still refresh it."""
    from jansky_research.report import preserve_live_results

    p = _results(tmp_path / "m.json", {"source": "synthetic", "n": 1})
    assert preserve_live_results({"source": "synthetic", "n": 2}, p)["n"] == 2


def test_results_an_unmarked_file_is_not_protected(tmp_path):
    """No source and no is_real means no provenance claim, so there is nothing to defend."""
    from jansky_research.report import preserve_live_results

    p = _results(tmp_path / "m.json", {"n": 1})
    assert preserve_live_results({"n": 2}, p)["n"] == 2
