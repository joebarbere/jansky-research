"""Plan 99 step 3, E1: zdm p(z|DM) per burst under the published HoffmannEmin25 state.

Runs OUTSIDE the repo environment (zdm pulls ~114 packages incl. lalsuite and breaks on
pandas 3's copy-on-write), standalone -- it does not import jansky_research::

    uv run --python 3.12 \
      --with "git+https://github.com/FRBs/zdm@e0f985bb55ad03d8b3435bc6a334a58d3ed674ba" \
      --with "pandas<3" python -I scripts/dmzcal_e1_zdm.py

Reads ``results/dmzcal_provenance.json``; writes ``results/dmzcal_e1_zdm.json`` (per-burst
PIT, posterior median, central 68/95% intervals, survey model). Verdicts are computed
in-repo by ``scripts/dmzcal_real.py`` against the frozen rule.

Survey model per discovery instrument (fixed before any posterior, plan 99 step 3):
  DSA -> ``DSA``; ASKAP -> ``CRAFT_average_ICS`` (the 2024 certification bursts are
  Shannon2024 ICS detections, not CRACO); MeerKAT -> ``MeerTRAPcoherent`` (beam mode not
  recorded); CHIME -> sum of ``CHIME_decbin_{0..5}_of_6`` grids WITHOUT repeater modelling
  -- what a user gets from zdm out of the box; HoffmannEmin25 did not fit CHIME (its
  Sec. 2.1), so CHIME results are reported per telescope.
DM_EG = DM - DM_ISM(NE2001, from the FRBs/FRB JSON) - 50 (HoffmannEmin25's DMhalo).
"""

from __future__ import annotations

import json
import sys
import tempfile
from importlib import resources
from pathlib import Path

import numpy as np
from zdm import cosmology as cos
from zdm import misc_functions, states, survey

ROOT = Path(__file__).resolve().parents[1]
STATE = "HoffmannEmin25"
DM_HALO = 50.0
N_PLANT = 500  # planted-truth replications through this script's own PIT code (GATE-2 #7)
SURVEY = {"DSA": "DSA", "ASKAP": "CRAFT_average_ICS", "MeerKAT": "MeerTRAPcoherent"}


def versions() -> dict[str, str]:
    from importlib.metadata import version  # noqa: PLC0415

    return {k: version(k) for k in ("zdm", "numpy", "scipy", "pandas", "astropy")}


def _pruned_dir(prune: set[str]) -> str:
    """Copy every survey file used into a temp dir with the given bursts' rows removed
    (GATE-2 round 2, N3: certification bursts listed in zdm's own survey files can nudge
    the survey efficiency through a median DM_G)."""
    from astropy.table import Table  # noqa: PLC0415

    src = resources.files("zdm").joinpath("data/Surveys/")
    tmp = Path(tempfile.mkdtemp())
    (tmp / "CHIME").mkdir()
    names = list(SURVEY.values()) + ["MeerTRAPincoherent"]
    names += [f"CHIME/CHIME_decbin_{i}_of_6" for i in range(6)]
    for name in names:
        t = Table.read(str(src.joinpath(name + ".ecsv")), format="ascii.ecsv")
        keep = [str(x).removeprefix("FRB") not in prune for x in t["TNS"]]
        t[keep].write(str(tmp / (name + ".ecsv")), format="ascii.ecsv", overwrite=True)
    return str(tmp)


def grids(prune: set[str] | None = None) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    st = states.load_state(STATE)
    cos.set_cosmology(st)
    cos.init_dist_measures()
    zdm_grid, zvals, dmvals = misc_functions.get_zdm_grid(
        st, new=True, plot=False, method="analytic"
    )
    sdir = _pruned_dir(prune) if prune else None
    rates: dict[str, np.ndarray] = {}
    for tel, name in list(SURVEY.items()) + [("MeerTRAPincoherent", "MeerTRAPincoherent")]:
        s = survey.load_survey(name, st, dmvals, sdir=sdir)
        g = misc_functions.initialise_grids([s], zdm_grid, zvals, dmvals, st, wdist=True)[0]
        rates[tel] = np.array(g.rates)
    chime_dir = (
        str(Path(sdir) / "CHIME") + "/"
        if sdir
        else str(resources.files("zdm").joinpath("data/Surveys/CHIME/"))
    )
    tot = None
    for i in range(6):
        s = survey.load_survey(f"CHIME_decbin_{i}_of_6", st, dmvals, sdir=chime_dir)
        g = misc_functions.initialise_grids([s], zdm_grid, zvals, dmvals, st, wdist=True)[0]
        tot = np.array(g.rates) if tot is None else tot + np.array(g.rates)
    assert tot is not None
    rates["CHIME"] = tot
    return rates, np.asarray(zvals), np.asarray(dmvals)


def summarise(pz: np.ndarray, zvals: np.ndarray, z_true: float) -> dict:
    dz = float(zvals[1] - zvals[0])
    lo = zvals - 0.5 * dz
    frac = np.clip((z_true - lo) / dz, 0.0, 1.0)
    edges = np.concatenate([[lo[0]], zvals + 0.5 * dz])
    cdf = np.concatenate([[0.0], np.cumsum(pz)])

    def q(p: float) -> float:
        return float(np.interp(p, cdf, edges))

    return {
        "pit": float(np.sum(pz * frac)),
        "z_median": q(0.5),
        "ci68": [q(0.16), q(0.84)],
        "ci95": [q(0.025), q(0.975)],
    }


def dm_pit_row(row: np.ndarray, dmvals: np.ndarray, dm_eg: float) -> float:
    """P(DM_EG' < DM_EG | z) from one z-row of the survey rates (DM bins of width ddm
    centred on dmvals, each treated as uniform). Selection-free in z: tests the DM model
    at the burst's known redshift (GATE-2 #2)."""
    ddm = float(dmvals[1] - dmvals[0])
    p = row / row.sum()
    lo = dmvals - 0.5 * ddm
    frac = np.clip((dm_eg - lo) / ddm, 0.0, 1.0)
    return float(np.sum(p * frac))


def survey_file_frbs() -> dict[str, list[str]]:
    """Burst names listed in the zdm survey files used (GATE-2 #9: they enter the survey
    efficiency through a median DM_G; recorded as a known provenance path)."""
    import pandas as pd  # noqa: PLC0415
    from astropy.table import Table  # noqa: PLC0415

    sdir = resources.files("zdm").joinpath("data/Surveys/")
    out = {}
    for name in list(SURVEY.values()) + [f"CHIME/CHIME_decbin_{i}_of_6" for i in range(6)]:
        t = Table.read(str(sdir.joinpath(name + ".ecsv")), format="ascii.ecsv")
        col = "TNS" if "TNS" in t.colnames else t.colnames[0]
        out[name] = sorted(str(x) for x in pd.Series(t[col]).astype(str))
    return out


def main() -> None:
    prune_mode = "--prune" in sys.argv[1:]
    prov = json.loads((ROOT / "results" / "dmzcal_provenance.json").read_text())
    alts = prov["corrections"]["alternatives"]
    cert = {r["name"].removeprefix("FRB") for r in prov["bursts"] if r["side"] == "certification"}
    rates, zvals, dmvals = grids(prune=cert if prune_mode else None)
    bursts = [r for r in prov["bursts"] if r["telescope"] in rates]
    rak = prov["recover_a_known"]
    bursts.append(rak | {"side": "recover_a_known", "secure_host": True})
    rng = np.random.default_rng(0)
    out = []
    for b in bursts:
        dm_eg = float(b["DM"]) - float(b["DMISM"]) - DM_HALO
        key = b.get("survey_override") or b["telescope"]
        z_true = float(b["z"])  # provenance already applies primary-source corrections
        rec = {
            "name": b["name"],
            "telescope": b["telescope"],
            "survey_model": b.get("survey_override")
            or SURVEY.get(b["telescope"], "CHIME_decbin_*_of_6 (summed)"),
            "side": b["side"],
            "spec_z": b["spec_z"],
            "secure_host": b.get("secure_host", True),
            "repeater": b.get("repeater", False),
            "z_true": z_true,
            "dm_eg": dm_eg,
        }
        if dm_eg <= 0:
            rec["excluded"] = "DM_EG <= 0"
        else:
            col = rates[key][:, int(np.argmin(np.abs(dmvals - dm_eg)))]
            if col.sum() <= 0:
                rec["excluded"] = "zero survey rate at this DM_EG"
            else:
                pz = col / col.sum()
                rec |= summarise(pz, zvals, z_true)
                iz = int(np.argmin(np.abs(zvals - z_true)))
                rec["dm_pit"] = dm_pit_row(rates[key][iz], dmvals, dm_eg)
                if b["side"] == "certification" and b["spec_z"] and b.get("secure_host"):
                    # post-hoc halo sensitivity (GATE-2 #13): same z, DM_halo 25 / 75
                    rec["pit_halo"] = {}
                    for h in (25.0, 75.0):
                        dmh = float(b["DM"]) - float(b["DMISM"]) - h
                        ch = rates[key][:, int(np.argmin(np.abs(dmvals - dmh)))]
                        if dmh > 0 and ch.sum() > 0:
                            rec["pit_halo"][f"{h:g}"] = summarise(ch / ch.sum(), zvals, z_true)[
                                "pit"
                            ]
                    if b["name"] in alts:  # unresolved source conflict (GATE-2 r2 N4)
                        a = alts[b["name"]]
                        z_a = float(a["z"]["value"]) if "z" in a else z_true
                        dm_a = (
                            (float(a["DM"]["value"]) if "DM" in a else float(b["DM"]))
                            - float(b["DMISM"])
                            - DM_HALO
                        )
                        ca = rates[key][:, int(np.argmin(np.abs(dmvals - dm_a)))]
                        rec["pit_alt"] = summarise(ca / ca.sum(), zvals, z_a)["pit"]
                    # planted truth: z drawn from this burst's own column (cell-uniform)
                    dz = float(zvals[1] - zvals[0])
                    idx = rng.choice(zvals.size, size=N_PLANT, p=pz)
                    zp = zvals[idx] + (rng.random(N_PLANT) - 0.5) * dz
                    rec["planted_pits"] = [summarise(pz, zvals, float(z))["pit"] for z in zp]
        out.append(rec)
    res = {
        "source": "real: zdm p(z|DM) (E1, state HoffmannEmin25) on FRBs/FRB localized hosts",
        "plan": "plans/99-dmzcal-coverage.md step 3",
        "zdm_commit": prov["zdm_commit"],
        "state": STATE,
        "dm_halo": DM_HALO,
        "grid": {
            "nz": int(zvals.size),
            "zmax": float(zvals[-1]),
            "ndm": int(dmvals.size),
            "dmmax": float(dmvals[-1]),
        },
        "survey_models": SURVEY | {"CHIME": "CHIME_decbin_0..5_of_6 summed, no repeaters"},
        "versions": versions(),
        "planted_truth": {"n_rep": N_PLANT, "seed": 0},
        "survey_files_pruned_of_certification_bursts": prune_mode,
        "survey_file_frbs": survey_file_frbs(),
        "bursts": out,
    }
    path = ROOT / "results" / ("dmzcal_e1_zdm_pruned.json" if prune_mode else "dmzcal_e1_zdm.json")
    path.write_text(json.dumps(res, indent=1) + "\n")
    print("wrote", path, len(out), "bursts")


if __name__ == "__main__":
    main()
