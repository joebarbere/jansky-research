"""PTE-II heavy-tail census -- real leg (plan 47, F10). Network + multi-GB SQLite; not run in CI.

Downloads the open PTE-II SQLite database (~1.5 GB zip, GitHub LFS, no auth), extracts it, loads the
per-pulsar single-pulse S/N arrays, runs the per-source giant-pulse test across all 363 pulsars, and
cross-matches the ATNF spin-down luminosity (Edot) to rank the giant-pulse excess against Edot. The
metric, its synthetic recover-a-known, and all unit tests run offline in core CI without any of this.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from .pte2 import (
    MIN_PULSES,
    PTE2_DB_NAME,
    PTE2_ZIP_URL,
    census,
    count_confound,
    load_pulse_sn,
    tail_vs_edot,
)

# Expected sha256 of Pulsar_fits_database_v1.zip (GATE-0 2026-07-10); verified after download.
PTE2_ZIP_SHA256_PREFIX = "775dffa0"


def fetch_pte2_db(cache_dir: str | Path) -> Path:  # pragma: no cover - large network download
    """Download + extract the PTE-II SQLite DB into ``cache_dir``; return the .db path (cached)."""
    import urllib.request
    import zipfile

    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    db = cache / PTE2_DB_NAME
    if db.exists():
        return db
    zip_path = cache / "Pulsar_fits_database_v1.zip"
    if not zip_path.exists():
        urllib.request.urlretrieve(PTE2_ZIP_URL, zip_path)  # noqa: S310 (trusted GitHub host)
    h = hashlib.sha256(zip_path.read_bytes()).hexdigest()
    if not h.startswith(PTE2_ZIP_SHA256_PREFIX):
        raise ValueError(f"PTE-II zip sha256 {h[:8]} != expected {PTE2_ZIP_SHA256_PREFIX}")
    with zipfile.ZipFile(zip_path) as z:
        name = next((n for n in z.namelist() if n.endswith(".db")), None)
        if name is None:
            raise ValueError(f"no .db in {zip_path.name}: {z.namelist()}")
        z.extract(name, cache)
        extracted = cache / name
        if extracted != db:
            extracted.rename(db)
    return db


def fetch_atnf_edot() -> dict[str, float]:  # pragma: no cover - network (VizieR)
    """ATNF spin-down luminosity by J-name: ``{jname: Edot [erg/s]}`` from VizieR ``B/psr``."""
    from astroquery.vizier import Vizier

    from .ppdot import spindown_luminosity

    v = Vizier(columns=["PSRJ", "P0", "P1"])
    v.ROW_LIMIT = -1
    t = v.get_catalogs("B/psr/psr")[0]
    out: dict[str, float] = {}
    for name, p0, p1 in zip(t["PSRJ"], t["P0"], t["P1"], strict=False):
        try:
            p, pd = float(p0), float(p1)
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(p) and np.isfinite(pd) and p > 0 and pd > 0):
            continue
        out[str(name).strip()] = float(spindown_luminosity(np.array([p]), np.array([pd]))[0])
    return out


def run_real_census(
    out: str, *, cache_dir: str | Path | None = None, min_pulses: int = MIN_PULSES
) -> dict:  # pragma: no cover - network + multi-GB SQLite
    """Full real census: fetch DB -> per-pulsar S/N -> giant-pulse test -> ATNF Edot cross-match."""
    import tempfile

    cache = Path(cache_dir) if cache_dir else Path(tempfile.gettempdir()) / "pte2_cache"
    db = fetch_pte2_db(cache)
    per = load_pulse_sn(db)
    n_total = len(per)
    rows = census(per, min_pulses=min_pulses)
    edot = fetch_atnf_edot()
    for r in rows:
        r["edot"] = edot.get(r["jname"])
    tv = tail_vs_edot(rows, edot)
    conf = count_confound(rows)
    n_fit = len(rows)
    n_heavy = int(sum(r["heavy_tailed"] for r in rows))
    # keep the top heavy-tailed sources in the artifact (drop the bulky per-pulse detail)
    top = [
        {
            k: r.get(k)
            for k in ("jname", "n", "n_giant", "excess", "gamma", "p_excess", "heavy_tailed")
        }
        for r in rows[:40]
    ]
    metrics = {
        "source": "PTE-II SQLite (Yang+2025); per-source giant-pulse test + ATNF Edot cross-match",
        "is_real": True,
        "n_pulsars_total": int(n_total),
        "n_fit": int(n_fit),
        "n_heavy": n_heavy,
        "heavy_fraction": round(n_heavy / n_fit, 4) if n_fit else float("nan"),
        "min_pulses": int(min_pulses),
        **tv,
        **conf,
        "top_heavy": top,
    }
    op = Path(out)
    (op / "results").mkdir(parents=True, exist_ok=True)
    import json

    (op / "results" / "pte2_census.json").write_text(
        json.dumps({**metrics, "all": rows}, indent=2) + "\n"
    )
    return metrics
