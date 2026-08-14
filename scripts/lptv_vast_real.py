#!/usr/bin/env python
"""LPT Stokes-V forced photometry, VAST extension --- real leg driver (lptv increment).

Extends the plan-44 RACS census (`lptv_real.py`) to the VAST survey's 12-minute epochs:
same forced I+V photometry at every v3 LPT position, in every public VAST I+V observation
covering it (`obs_collection = 'VAST'`; identical filename conventions). VAST monitors its
fields roughly fortnightly, so covered sources gain tens of epochs over the ~10 RACS ones ---
enough for a per-source epoch count that supports duty-cycle statements.

Two deliberate differences from the RACS driver, both lessons from the round-2 referee:
- `epoch_mjd` is written at full precision (`t_min:.5f`, ~1 s), not `:.1f` (+/-72 min), so
  phase arithmetic against published ephemerides is possible from the committed CSV.
- `duration_s` (obscore t_exptime) is recorded per row, so phase-coverage arithmetic does not
  need an assumed integration time.

Run:  uv run --extra vlass python scripts/lptv_vast_real.py [--limit N]

Resumable: (name, obs_id) rows already measured are skipped; failed rows retry on the next run.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from jansky_research.lptv import lpt_positions  # noqa: E402
from jansky_research.stokesv import _casda_session  # noqa: E402
from wdpulsar_real import (  # noqa: E402 - reuse the plan-41 CASDA sweep helpers
    complete_iv_groups,
    measure_group,
    obscore_products,
)

CASDA_TAP = "https://casda.csiro.au/casda_vo_tools/tap"
CSV_PATH = REPO / "results" / "lptv_vast_epochs.csv"
CSV_FIELDS = [
    "name",
    "epoch",
    "obs_id",
    "band",
    "epoch_mjd",
    "duration_s",
    "i_mjy",
    "e_i",
    "v_mjy",
    "e_v",
    "offset_arcsec",
    "note",
]


def load_done(path: Path) -> set[tuple[str, str]]:
    """Rows already measured. Failed rows do NOT count: they retry on the next run."""
    if not path.exists():
        return set()
    with path.open() as fh:
        return {
            (r["name"], r["obs_id"])
            for r in csv.DictReader(fh)
            if not r.get("note", "").startswith("failed")
        }


def main() -> int:
    ap = argparse.ArgumentParser(description="LPT Stokes-V forced photometry, VAST extension")
    ap.add_argument("--limit", type=int, default=0, help="limit targets (0 = all 16)")
    ap.add_argument("--csv", default=str(CSV_PATH))
    args = ap.parse_args()

    import pyvo

    username = os.environ.get("CASDA_USERNAME", "joe.barbere@gmail.com")
    pw_path = "~/.casda_pw"
    out = Path(args.csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    done = load_done(out)
    print(f"resuming: {len(done)} (name, obs_id) rows already in {out}", flush=True)

    targets = lpt_positions()
    if args.limit:
        targets = targets[: args.limit]
    print(f"{len(targets)} LPT positions (v3 catalogue), collection=VAST", flush=True)

    tap = pyvo.dal.TAPService(CASDA_TAP)
    casda = None
    new_file = not out.exists() or out.stat().st_size == 0
    fh = out.open("a", newline="")
    writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
    if new_file:
        writer.writeheader()
        fh.flush()

    def emit(rowdict):
        writer.writerow(rowdict)
        fh.flush()
        os.fsync(fh.fileno())

    t0 = time.time()
    for it, tgt in enumerate(targets):
        print(
            f"[{it + 1}/{len(targets)}] {tgt['name']} (P={tgt['period_s'] / 60:.0f} min)",
            flush=True,
        )
        try:
            table = obscore_products(tap, tgt["ra_deg"], tgt["dec_deg"], collection="VAST")
            groups = complete_iv_groups(table)
        except Exception as exc:  # noqa: BLE001
            print(f"    obscore query failed: {exc!r} (will retry next run)", flush=True)
            continue
        pending = [g for g in groups if (tgt["name"], g["obs_id"]) not in done]
        print(f"    {len(groups)} complete I+V groups, {len(pending)} pending", flush=True)
        for g in pending:
            base = {
                "name": tgt["name"],
                "epoch": g["band"],
                "obs_id": g["obs_id"],
                "band": g["band"],
                "epoch_mjd": f"{g['t_min']:.5f}",
                "duration_s": f"{g.get('t_exptime', float('nan')):.1f}",
                "i_mjy": "nan",
                "e_i": "nan",
                "v_mjy": "nan",
                "e_v": "nan",
                "offset_arcsec": "nan",
                "note": "",
            }
            try:
                meas, casda = measure_group(
                    casda, username, pw_path, table, g, tgt["ra_deg"], tgt["dec_deg"]
                )
                base.update(
                    {
                        "i_mjy": f"{meas['i_peak']:.4f}",
                        "e_i": f"{meas['i_rms']:.4f}",
                        "v_mjy": f"{meas['v_peak']:.4f}",
                        "e_v": f"{meas['v_rms']:.4f}",
                        "offset_arcsec": f"{meas['offset_arcsec']:.2f}",
                    }
                )
            except Exception as exc:  # noqa: BLE001
                base["note"] = f"failed: {type(exc).__name__}"
                print(f"    {g['obs_id']}: FAILED {exc!r}", flush=True)
            emit(base)
        print(f"    elapsed {(time.time() - t0) / 60:.1f} min", flush=True)
    fh.close()
    print("sweep complete", flush=True)
    return 0


if __name__ == "__main__":
    _ = _casda_session  # kept importable for parity with the wdpulsar driver
    raise SystemExit(main())
