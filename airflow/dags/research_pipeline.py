"""Airflow DAG: the FRB burst-statistics research pipeline.

The automation layer over ``jansky_research``. Every task is an **in-process PythonOperator** that
imports and calls ``jansky_research`` — *not* a DockerOperator/PodmanOperator spawning containers
(rootless-podman socket + nested-container volume relabeling is the most painful configuration,
and the package is already installed in this image). This is the **same code path** that
``make pipeline`` and the notebooks use, so the DAG and the Makefile produce identical artifacts.

Tasks:  fetch_dataset -> run_analysis -> make_figures -> assemble_paper_inputs

Outputs land under ``$JR_OUTPUT`` (a bind-mounted repo checkout): ``results/metrics.json``,
``paper/figures/*.pdf`` and ``paper/generated/macros.tex``. Manual trigger (``schedule=None``);
a research pipeline is run on demand, not on a cron.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pendulum
from airflow.decorators import dag, task

OUT = Path(os.environ.get("JR_OUTPUT", "/opt/airflow/project"))
OFFLINE = os.environ.get("JR_OFFLINE", "0") == "1"
CATALOG_NPZ = OUT / "results" / "catalog.npz"


@dag(
    dag_id="research_pipeline",
    schedule=None,
    catchup=False,
    start_date=pendulum.datetime(2026, 1, 1, tz="UTC"),
    tags=["jansky-research", "frb"],
    doc_md=__doc__,
)
def research_pipeline():
    @task()
    def fetch_dataset() -> str:
        """Build the catalogue (real CHIME CSV, or synthetic fixture offline) and cache it."""
        from jansky_research import pipeline

        catalog, source = pipeline.build_catalog(offline=OFFLINE)
        CATALOG_NPZ.parent.mkdir(parents=True, exist_ok=True)
        np.savez(CATALOG_NPZ, **catalog)
        return source

    @task()
    def run_analysis(source: str) -> None:
        """Run the burst-statistics analysis and write results/metrics.json."""
        from jansky_research import pipeline

        data = np.load(CATALOG_NPZ, allow_pickle=False)
        catalog = {k: data[k] for k in data.files}
        metrics = pipeline.analyze(catalog, source=source)
        (OUT / "results" / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    @task()
    def make_figures() -> None:
        """Emit the result figures into paper/figures/."""
        from jansky_research import frbstats, report

        data = np.load(CATALOG_NPZ, allow_pickle=False)
        catalog = {k: data[k] for k in data.files}
        report.make_figures(catalog, frbstats.summarise(catalog), OUT / "paper" / "figures")

    @task()
    def assemble_paper_inputs() -> None:
        """Write paper/generated/macros.tex (every headline number, for the LaTeX)."""
        from jansky_research import report

        metrics = json.loads((OUT / "results" / "metrics.json").read_text())
        report.write_macros(metrics, OUT / "paper" / "generated" / "macros.tex")

    source = fetch_dataset()
    analysis = run_analysis(source)
    analysis >> make_figures() >> assemble_paper_inputs()


research_pipeline()
