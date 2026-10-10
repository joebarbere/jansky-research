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
from importlib import resources
from pathlib import Path

import numpy as np
from zdm import cosmology as cos
from zdm import misc_functions, states, survey

ROOT = Path(__file__).resolve().parents[1]
STATE = "HoffmannEmin25"
DM_HALO = 50.0
Z_OVERRIDE = {"FRB20231201A": 0.1119}  # Leung+2025 primary table (step-0 finding 7)
SURVEY = {"DSA": "DSA", "ASKAP": "CRAFT_average_ICS", "MeerKAT": "MeerTRAPcoherent"}


def grids() -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    st = states.load_state(STATE)
    cos.set_cosmology(st)
    cos.init_dist_measures()
    zdm_grid, zvals, dmvals = misc_functions.get_zdm_grid(
        st, new=True, plot=False, method="analytic"
    )
    rates: dict[str, np.ndarray] = {}
    for tel, name in SURVEY.items():
        s = survey.load_survey(name, st, dmvals)
        g = misc_functions.initialise_grids([s], zdm_grid, zvals, dmvals, st, wdist=True)[0]
        rates[tel] = np.array(g.rates)
    chime_dir = str(resources.files("zdm").joinpath("data/Surveys/CHIME/"))
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


def main() -> None:
    prov = json.loads((ROOT / "results" / "dmzcal_provenance.json").read_text())
    rates, zvals, dmvals = grids()
    bursts = [r for r in prov["bursts"] if r["telescope"] in rates]
    rak = prov["recover_a_known"]
    bursts.append(rak | {"side": "recover_a_known", "secure_host": True})
    out = []
    for b in bursts:
        dm_eg = float(b["DM"]) - float(b["DMISM"]) - DM_HALO
        z_true = Z_OVERRIDE.get(b["name"], float(b["z"]))
        rec = {
            "name": b["name"],
            "telescope": b["telescope"],
            "survey_model": SURVEY.get(b["telescope"], "CHIME_decbin_*_of_6 (summed)"),
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
            col = rates[b["telescope"]][:, int(np.argmin(np.abs(dmvals - dm_eg)))]
            if col.sum() <= 0:
                rec["excluded"] = "zero survey rate at this DM_EG"
            else:
                rec |= summarise(col / col.sum(), zvals, z_true)
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
        "z_overrides": Z_OVERRIDE,
        "bursts": out,
    }
    path = ROOT / "results" / "dmzcal_e1_zdm.json"
    path.write_text(json.dumps(res, indent=1) + "\n")
    print("wrote", path, len(out), "bursts")


if __name__ == "__main__":
    main()
