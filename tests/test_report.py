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
