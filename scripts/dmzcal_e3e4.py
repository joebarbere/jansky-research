"""Plan 99 step 3, E3 (fruitbat 2.0.1) and E4 (FRBs/FRB z_from_DM) per burst.

Both run OUTSIDE the repo environment and do not import jansky_research.

E3 -- fruitbat no longer runs on a current stack (astropy>=6 removed
``astropy.cosmology.core``; it needs ``np.long``; it imports ``pygedm`` at load, which does
not build here). Working recipe (plan 99 "Estimators under audit")::

    uv venv -p 3.9 /tmp/venv-fb
    uv pip install -p /tmp/venv-fb/bin/python "numpy==1.22.4" "astropy==4.3.1" \
        "scipy<1.10" "h5py<3.8" "matplotlib<3.7" e13tools pandas "setuptools<70"
    uv pip install -p /tmp/venv-fb/bin/python --no-deps fruitbat==2.0.1
    /tmp/venv-fb/bin/python -I scripts/dmzcal_e3e4.py e3

A ``pygedm`` stub that RAISES if called is written to a temp dir and put on sys.path, so it
cannot silently supply a Milky Way DM: ``dm_galaxy`` is always passed explicitly.
Package defaults: method ``Batten2021``, ``prior="uniform"``, ``subtract_host=False``, no halo.

E4 (point estimate only)::

    uv run --python 3.12 --with "git+https://github.com/FRBs/FRB@996fcda9b0b22431e3171208b4c8e1cf2798e823" \
        python -I scripts/dmzcal_e3e4.py e4

``z_from_DM(DM - DM_ISM, corr_nuisance=True)`` -- its default -- no interval, so E4 enters
only the point-accuracy table.
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
Z_OVERRIDE = {"FRB20231201A": 0.1119}


def _bursts() -> list[dict]:
    prov = json.loads((ROOT / "results" / "dmzcal_provenance.json").read_text())
    rows = list(prov["bursts"])
    rows.append(prov["recover_a_known"] | {"side": "recover_a_known", "secure_host": True})
    return rows


def _base(b: dict) -> dict:
    return {
        "name": b["name"],
        "telescope": b["telescope"],
        "side": b["side"],
        "spec_z": b["spec_z"],
        "secure_host": b.get("secure_host", True),
        "repeater": b.get("repeater", False),
        "z_true": Z_OVERRIDE.get(b["name"], float(b["z"])),
        "dm": float(b["DM"]),
        "dm_ism": float(b["DMISM"]),
    }


def run_e3() -> dict:
    stub = Path(tempfile.mkdtemp())
    (stub / "pygedm.py").write_text(
        "def __getattr__(name):\n"
        "    raise RuntimeError(f'pygedm stub: {name} called; pass dm_galaxy explicitly')\n"
    )
    sys.path.insert(0, str(stub))
    import fruitbat  # noqa: PLC0415  (only importable in the py3.9 venv)

    out = []
    for b in _bursts():
        rec = _base(b)
        try:
            f = fruitbat.Frb(rec["dm"], dm_galaxy=rec["dm_ism"])
            zb, zpdf, dz = f.calc_redshift_pdf(method="Batten2021")
            zb, zpdf, dz = np.asarray(zb), np.asarray(zpdf), np.asarray(dz)
            mass = zpdf * dz
            if not np.isfinite(mass.sum()) or mass.sum() <= 0:
                raise ValueError("empty posterior (fruitbat sets dm_excess < 0 to 0)")
            mass = mass / mass.sum()
            lo = zb - dz  # z_bins are the upper bin edges
            frac = np.clip((rec["z_true"] - lo) / dz, 0.0, 1.0)
            edges = np.concatenate([[lo[0]], zb])
            cdf = np.concatenate([[0.0], np.cumsum(mass)])
            rec |= {
                "pit": float(np.sum(mass * frac)),
                "z_median": float(np.interp(0.5, cdf, edges)),
                "ci68": [float(np.interp(p, cdf, edges)) for p in (0.16, 0.84)],
                "ci95": [float(np.interp(p, cdf, edges)) for p in (0.025, 0.975)],
            }
        except Exception as e:  # noqa: BLE001 -- record, never drop silently
            rec["excluded"] = f"{type(e).__name__}: {e}"
        out.append(rec)
    return {
        "source": "real: fruitbat 2.0.1 Batten2021 p(z|DM) (E3) on FRBs/FRB localized hosts",
        "estimator": "E3",
        "settings": {
            "method": "Batten2021",
            "prior": "uniform",
            "subtract_host": False,
            "halo": "not subtracted (package default)",
            "dm_galaxy": "FRBs/FRB DMISM",
        },
        "env": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "fruitbat": getattr(fruitbat, "__version__", "2.0.1"),
        },
        "bursts": out,
    }


def run_e4() -> dict:
    import astropy.units as u  # noqa: PLC0415
    from frb.dm import igm  # noqa: PLC0415

    out = []
    for b in _bursts():
        rec = _base(b)
        try:
            z = igm.z_from_DM((rec["dm"] - rec["dm_ism"]) * u.pc / u.cm**3, corr_nuisance=True)
            rec["z_point"] = float(z)
        except Exception as e:  # noqa: BLE001
            rec["excluded"] = f"{type(e).__name__}: {e}"
        out.append(rec)
    return {
        "source": "real: FRBs/FRB z_from_DM point estimate (E4) on FRBs/FRB localized hosts",
        "estimator": "E4",
        "settings": {"corr_nuisance": True, "dm_input": "DM - DMISM(NE2001)"},
        "bursts": out,
    }


def main() -> None:
    which = sys.argv[1]
    res = {"e3": run_e3, "e4": run_e4}[which]()
    res["plan"] = "plans/99-dmzcal-coverage.md step 3"
    path = ROOT / "results" / f"dmzcal_{which}.json"
    path.write_text(json.dumps(res, indent=1) + "\n")
    n_ex = sum("excluded" in r for r in res["bursts"])
    print("wrote", path, len(res["bursts"]), "bursts,", n_ex, "excluded")


if __name__ == "__main__":
    main()
