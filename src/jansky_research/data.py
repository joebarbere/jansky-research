"""Research-dataset registry + offline synthetic fallback.

Mirrors ``jansky.data`` in spirit: a small registry of *lightweight, openly
downloadable* public datasets, cached locally so analysis is reproducible, plus
a synthetic fallback so tests and CI never depend on the network.

The chosen-domain dataset is added to :data:`DATASETS` after GATE 1 (the gap
analysis picks ONE openly-downloadable dataset). Until then the registry seeds a
couple of representative public products and the synthetic generators used by the
offline tests.

Usage
-----
From Python::

    from jansky_research import data
    path = data.fetch("chime-frb-catalog")          # cached download
    arr = data.synthetic_dynamic_spectrum()          # offline fallback

From the command line::

    python -m jansky_research.data --list
    python -m jansky_research.data --fetch chime-frb-catalog
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "DATASETS",
    "Dataset",
    "data_dir",
    "fetch",
    "list_datasets",
    "synthetic_dynamic_spectrum",
]


@dataclass(frozen=True)
class Dataset:
    """Metadata for a downloadable research dataset.

    ``category`` is ``"small"`` for the lightweight products on the default path
    or ``"large"`` for opt-in bulk data kept off the default/offline path.
    """

    name: str
    url: str
    filename: str
    description: str
    size_hint: str = "unknown"
    category: str = "small"


# Seed registry. GATE 1 adds the chosen dataset; these are representative,
# openly-downloadable starting points. jansky already vendors several sub-MB real
# files (PSRFITS/filterbank, the NANOGrav NGC 6440E .par/.tim, a LAB HI slice) via
# ``jansky.data`` — reuse those directly rather than duplicating them here.
DATASETS: dict[str, Dataset] = {
    "chime-frb-catalog": Dataset(
        name="chime-frb-catalog",
        # Official source is https://www.chime-frb.ca/catalog (often 503); this is the
        # byte-identical chimefrbcat1.csv vendored by the maintained frbpoppy package.
        url="https://raw.githubusercontent.com/TRASAL/frbpoppy/master/data/frbcat/chimefrbcat1.csv",
        filename="chimefrbcat1.csv",
        description="CHIME/FRB Catalog 1 — public CSV of 600 fast radio bursts (DM, fluence, "
        "width, repeater name) from 18 repeaters + non-repeaters. The dataset for the FRB "
        "burst-statistics analysis (CHIME/FRB Collaboration 2021, ApJS 257, 59).",
        size_hint="~216 KB",
    ),
    "vlass-cirada-ql": Dataset(
        name="vlass-cirada-ql",
        # The full CIRADA VLASS Quick-Look component catalogues (~3.4M components/epoch) are bulk
        # FITS products; in practice the slice queries a sky region per epoch via CADC/CIRADA TAP
        # rather than downloading the whole table. Listed here for provenance.
        url="https://cirada.ca/vlasscatalogueql0",
        filename="CIRADA_VLASS_QL_catalogue.fits",
        description="VLA Sky Survey (VLASS) Quick-Look component catalogues (CIRADA; Gordon et al. "
        '2021), 2-4 GHz, 2.5", three epochs (2017-2024). Cross-matched across epochs for the '
        "multi-epoch radio-variability slice (jansky_research.vlass).",
        size_hint="~GB per epoch",
        category="large",
    ),
    "astrogeo-vlbi": Dataset(
        name="astrogeo-vlbi",
        # The Astrogeo VLBI archive is per-source pages, not one bulk table; the slice fetches each
        # source's S/X flux-density history over HTTP (no auth). Listed here for provenance.
        url="https://astrogeo.smce.nasa.gov/vlbi_images",
        filename="astrogeo_vlbi_index.html",
        description="Astrogeo VLBI image database (Petrov) — ~139k dual-band S/X (2.3/8.4 GHz) "
        "images of ~21k compact sources from decades of geodetic/astrometric VLBI. Per-source "
        "multi-decade flux histories drive the VLBI variability slice (jansky_research.vlbi).",
        size_hint="per-source pages",
        category="large",
    ),
    "ecallisto": Dataset(
        name="ecallisto",
        # The e-Callisto archive is per-station, per-15-minute gzipped FITS, not one table; the slice
        # fetches a file covering an event over HTTP (no auth). Listed here for provenance.
        url="http://soleil.i4ds.ch/solarradio/data/2002-20yy_Callisto",
        filename="ecallisto_index.html",
        description="e-Callisto network (i4ds/FHNW; Benz et al. 2009) — ground-based solar radio "
        "dynamic spectra (~20-900 MHz, 15-min FITS) from a global station network. Drives the "
        "type III burst drift / exciter-speed slice (jansky_research.solarbursts).",
        size_hint="per-station 15-min FITS",
        category="large",
    ),
    "wind-waves": Dataset(
        name="wind-waves",
        # Wind/WAVES Level-2 radio CDFs on SPDF, one file per day per receiver; the slice fetches by
        # date and parses with cdflib. Listed here for provenance.
        url="https://spdf.gsfc.nasa.gov/pub/data/wind/waves/rad2_l2",
        filename="wi_l2_wav_rad2.cdf",
        description="Wind/WAVES Level-2 radio dynamic spectra (Bougeret et al. 1995) — RAD1 "
        "(0.02-1.04 MHz) and RAD2 (1.075-13.825 MHz), space-based, public over SPDF (no auth). "
        "Drives the interplanetary type III slice (jansky_research.windwaves; needs the windwaves "
        "extra, cdflib).",
        size_hint="~1 MB per day per receiver",
        category="large",
    ),
    "atnf-psrcat": Dataset(
        name="atnf-psrcat",
        # The ATNF Pulsar Catalogue is on VizieR (B/psr); the pulsar slices query it directly.
        url="https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=B/psr",
        filename="atnf_psrcat.vot",
        description="ATNF Pulsar Catalogue (Manchester et al. 2005) via VizieR B/psr — periods, "
        "period derivatives, and flux densities for the known pulsars. Drives the pulsar spectral "
        "(jansky_research.pulsarspec) and P-Pdot (jansky_research.ppdot) slices.",
        size_hint="~3,500 rows",
    ),
    "taylor2009-rm": Dataset(
        name="taylor2009-rm",
        # The Taylor+2009 NVSS RM catalogue is on VizieR; the slice queries it directly rather than
        # caching the file. Listed here for provenance.
        url="https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/ApJ/702/1230",
        filename="taylor2009_rm.vot",
        description="Taylor, Stil & Sunstrum (2009) NVSS rotation-measure catalogue — RM for 37,543 "
        "polarised extragalactic sources (VizieR J/ApJ/702/1230). Drives the Galactic Faraday "
        "rotation-sky slice (jansky_research.rmsky).",
        size_hint="37,543 rows",
    ),
    "voyager1-h5": Dataset(
        name="voyager1-h5",
        url="http://blpd0.ssl.berkeley.edu/Voyager_data/Voyager1.single_coarse.fine_res.h5",
        filename="Voyager1.single_coarse.fine_res.h5",
        description="Breakthrough Listen GBT X-band fine-resolution observation of Voyager 1 "
        "(2015-12-30) — a real narrowband spacecraft carrier near 8420 MHz; the canonical SETI "
        "drift-search validation file. Bitshuffle-compressed HDF5 (needs the 'voyager' extra).",
        size_hint="~50 MB",
        category="large",
    ),
}


def data_dir() -> Path:
    """Return the directory where datasets are cached.

    Defaults to ``<repo>/data``; override with ``JANSKY_RESEARCH_DATA_DIR``.
    Created on first use.
    """
    env = os.environ.get("JANSKY_RESEARCH_DATA_DIR")
    if env:
        path = Path(env).expanduser()
    else:
        # src/jansky_research/data.py -> repo root is three parents up.
        path = Path(__file__).resolve().parents[2] / "data"
    path.mkdir(parents=True, exist_ok=True)
    return path


def list_datasets() -> list[str]:
    """Return the names of all registered datasets."""
    return sorted(DATASETS)


def fetch(name: str, *, force: bool = False) -> Path:
    """Download (and cache) a registered dataset, returning its local path.

    Raises
    ------
    KeyError
        If ``name`` is not registered.
    RuntimeError
        If the download fails (e.g. offline); the message points at the synthetic
        fallback for offline work.
    """
    if name not in DATASETS:
        raise KeyError(f"Unknown dataset {name!r}. Known datasets: {', '.join(list_datasets())}")
    spec = DATASETS[name]
    target = data_dir() / spec.filename
    if target.exists() and not force:
        return target
    try:
        _download(spec.url, target)
    except Exception as exc:  # noqa: BLE001 - re-raised with guidance
        raise RuntimeError(
            f"Failed to download {name!r} from {spec.url}: {exc}. "
            "For offline work use jansky_research.data.synthetic_dynamic_spectrum() "
            "(or a jansky-vendored sample via jansky.data)."
        ) from exc
    return target


def _download(url: str, target: Path) -> None:  # pragma: no cover - network
    """Stream ``url`` to ``target`` atomically, with a progress bar."""
    import requests
    from tqdm import tqdm

    tmp = target.with_suffix(target.suffix + ".part")
    with requests.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        with (
            open(tmp, "wb") as fh,
            tqdm(total=total, unit="B", unit_scale=True, desc=target.name) as bar,
        ):
            for chunk in resp.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
                bar.update(len(chunk))
    tmp.replace(target)


def synthetic_dynamic_spectrum(
    n_time: int = 256,
    n_chan: int = 128,
    seed: int | None = 0,
) -> np.ndarray:
    """Generate a small synthetic time-frequency array for offline tests/demos.

    Gaussian noise with a faint broadband transient near the middle time sample,
    enough to exercise downstream analysis (de-dispersion, SNR, statistics)
    without any network access.

    Returns
    -------
    numpy.ndarray
        A ``(n_time, n_chan)`` intensity array.
    """
    rng = np.random.default_rng(seed)
    arr = rng.normal(0.0, 1.0, size=(n_time, n_chan))
    arr[n_time // 2] += 5.0  # a faint broadband pulse
    return arr


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Fetch jansky-research datasets.")
    parser.add_argument("--list", action="store_true", help="list known datasets")
    parser.add_argument("--fetch", metavar="NAME", help="download a dataset by name")
    parser.add_argument("--force", action="store_true", help="re-download if cached")
    args = parser.parse_args(argv)

    if args.list or not args.fetch:
        print(f"Cache directory: {data_dir()}")
        for name in list_datasets():
            spec = DATASETS[name]
            print(f"  {name:22s} {spec.size_hint:>10s}  {spec.description}")
        if not args.fetch:
            return 0

    path = fetch(args.fetch, force=args.force)
    print(f"Ready: {path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
