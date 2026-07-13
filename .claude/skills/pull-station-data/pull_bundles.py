#!/usr/bin/env python3
"""Pull codified observation bundles from a jansky-observe station into this repo.

The station (`../jansky-observe`) exports one documented bundle per observation —
schema ``jansky-observe.observation-bundle/1``: a ``bundle.json`` manifest (station
UUID, pointing, LST at each capture's start, timestamps, SDR settings/gain,
cal-epoch reference, classifier verdicts) plus one self-describing averaged-
spectrum ``capture-<id>.npz`` per npz capture. That is exactly the "averaged-
spectra format from the station's capture service" plan 78's ``hline.read_capture``
consumes.

This script is the **bulk, non-interactive** half of the ``pull-station-data``
skill: it hits the station's plain JSON+zip HTTP API (the same data the MCP tools
``list_observations`` / ``get_observation_bundle`` proxy) and unpacks each bundle
into ``<out>/<station-uuid>/observation-<id>/``. Stdlib only — no dependency, so it
runs anywhere the station is reachable on the LAN.

It only reads: it GETs identity, the observation list, and per-observation bundle
zips. It never mutates the station (the station's MCP surface carries no write or
delete verbs anyway, by design).

Exit codes:
  0  success: every requested bundle fetched and unpacked
  2  bad usage
  3  station unreachable / identity fetch failed
  4  the observation list could not be fetched
  5  one or more requested observations failed to fetch or unpack
  6  a bundle downloaded but was not a valid zip (station/network problem)
"""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

DEFAULT_STATION_URL = os.environ.get("JANSKY_STATION_URL", "http://raspberrypi.local:8000")


def _log(msg: str) -> None:
    print(f"[pull-station-data] {msg}", file=sys.stderr)


def _get(url: str, timeout: float) -> bytes:
    """GET a URL and return the raw body (raises urllib errors)."""
    request = urllib.request.Request(url, headers={"Accept": "*/*"})
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310 (LAN station)
        return response.read()


def _get_json(url: str, timeout: float) -> Any:
    return json.loads(_get(url, timeout).decode("utf-8"))


def fetch_identity(station_url: str, timeout: float) -> dict[str, Any]:
    """Fetch the station identity ({uuid, name, ...}); raise on failure."""
    return _get_json(f"{station_url}/api/station", timeout)


def select_observations(
    station_url: str,
    timeout: float,
    *,
    ids: list[int] | None,
    status: str,
    limit: int,
) -> list[dict[str, Any]]:
    """Return the observations to pull — explicit ids, or the list filtered by status."""
    rows: list[dict[str, Any]] = _get_json(f"{station_url}/api/observations", timeout)
    if ids:
        wanted = set(ids)
        return [r for r in rows if r.get("id") in wanted]
    if status != "any":
        rows = [r for r in rows if r.get("status") == status]
    return rows[:limit]


def pull_bundle(
    station_url: str,
    observation_id: int,
    out_dir: Path,
    timeout: float,
) -> dict[str, Any]:
    """Fetch + unpack one observation bundle; return its manifest.

    Raises
    ------
    zipfile.BadZipFile
        If the download is not a valid zip.
    urllib.error.URLError / HTTPError
        On network/HTTP failure.
    """
    url = f"{station_url}/api/observations/{observation_id}/bundle"
    payload = _get(url, timeout)
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        manifest = json.loads(archive.read("bundle.json"))
        out_dir.mkdir(parents=True, exist_ok=True)
        for name in archive.namelist():
            # Guard against path traversal from a hostile archive (defensive).
            target = (out_dir / name).resolve()
            if not str(target).startswith(str(out_dir.resolve())):
                raise ValueError(f"refusing to extract outside out dir: {name!r}")
            archive.extract(name, out_dir)
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Pull jansky-observe observation bundles into this repo.",
    )
    parser.add_argument(
        "--station-url",
        default=DEFAULT_STATION_URL,
        help=f"station base URL (default {DEFAULT_STATION_URL}; env JANSKY_STATION_URL)",
    )
    parser.add_argument(
        "--out",
        default="data/station",
        help="output root (default data/station); bundles land in <out>/<uuid>/observation-<id>/",
    )
    parser.add_argument(
        "--observation-id",
        type=int,
        action="append",
        dest="observation_ids",
        help="pull only this observation id (repeatable); default: pull the list",
    )
    parser.add_argument(
        "--status",
        default="done",
        help="filter the list by status when no ids are given (default 'done'; 'any' = no filter)",
    )
    parser.add_argument("--limit", type=int, default=50, help="max observations from the list")
    parser.add_argument(
        "--list",
        action="store_true",
        help="only list what would be pulled; download nothing",
    )
    parser.add_argument("--timeout", type=float, default=30.0, help="per-request timeout (s)")
    args = parser.parse_args(argv)

    station_url = args.station_url.rstrip("/")

    try:
        identity = fetch_identity(station_url, args.timeout)
    except (urllib.error.URLError, OSError, ValueError) as exc:
        _log(f"station unreachable at {station_url}: {exc}")
        return 3
    station_uuid = str(identity.get("uuid", "unknown"))
    _log(f"station: {identity.get('name')!r} uuid={station_uuid}")

    try:
        observations = select_observations(
            station_url,
            args.timeout,
            ids=args.observation_ids,
            status=args.status,
            limit=args.limit,
        )
    except (urllib.error.URLError, OSError, ValueError) as exc:
        _log(f"could not fetch the observation list: {exc}")
        return 4

    if not observations:
        _log("no matching observations to pull.")
        return 0

    if args.list:
        for row in observations:
            _log(f"  #{row.get('id')}  {row.get('status'):>8}  {row.get('name')}")
        _log(f"{len(observations)} observation(s) would be pulled (nothing downloaded).")
        return 0

    out_root = Path(args.out) / station_uuid
    failures = 0
    pulled = 0
    for row in observations:
        obs_id = row.get("id")
        if obs_id is None:
            continue
        dest = out_root / f"observation-{obs_id}"
        try:
            manifest = pull_bundle(station_url, int(obs_id), dest, args.timeout)
        except zipfile.BadZipFile as exc:
            _log(f"  #{obs_id}: not a valid bundle zip: {exc}")
            failures = max(failures, 1)
            return 6
        except (urllib.error.URLError, OSError, ValueError) as exc:
            _log(f"  #{obs_id}: fetch/unpack failed: {exc}")
            failures += 1
            continue
        n_spectra = sum(1 for c in manifest.get("captures", []) if c.get("spectrum_file"))
        _log(
            f"  #{obs_id}: {len(manifest.get('captures', []))} capture(s), {n_spectra} npz -> {dest}"
        )
        pulled += 1

    _log(f"done: {pulled} bundle(s) into {out_root} ({failures} failed).")
    return 5 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
