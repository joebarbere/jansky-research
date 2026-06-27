# Survey — open-source tooling & automation landscape

Synthesised from the cross-cutting ecosystem sweep + five domain briefings. All packages
URL-verified (June 2026); activity judged by latest release.

## Coverage map

**Well-covered (mature, maintained, tested):** core I/O & formats (Astropy, pyuvdata,
spectral-cube/radio-beam); pulsar timing (PINT, enterprise); CASA calibration/imaging; CPU pulsar
search (PRESTO, riptide, sigpyproc3); continuum source finding (PyBDSF, Aegean, SoFiA-2).

**Underserved / problematic:**
- **GPU-only, no CPU path:** AstroAccelerate, Heimdall (also quasi-abandoned, no CI), hyperseti.
- **Abandoned / stale:** blimpy (2022), turbo_seti (2022), fruitbat, burstfit, frbpa, pvextractor,
  PsrPopPy, ClaRAN, RASCIL (discontinued).
- **Heavy installs (no pip / compiled):** CASA, WSClean, Caesar, PRESTO, TEMPO2/libstempo, PyBDSF,
  SoFiA-2, BLISS.
- **Not offline by default:** frbpoppy, CASA (casadata), ClaRAN (fetch at runtime).
- **ML in radio astronomy:** thin, fragmented, GPU-assumed, mostly no-CI research code.

## The automation / orchestration gap (key finding)

Every serious radio pipeline assumes observatory HPC (Slurm/GPU + proprietary data) or a
domain-built engine. General orchestrators appear *only* at professional facilities:
- **Airflow** — AGLOW (LOFAR, arXiv:1808.10735); **CHIME/FRB baseband uses Airflow + Docker +
  Docker Swarm** (arXiv:2010.06748). HPC/Docker only.
- **Prefect** — Simons Observatory (arXiv:2406.10905). **Snakemake** — showyourwork (reproducible
  *papers*, laptop-scale — the one independent-researcher tool, but not data reduction).
  **Nextflow/CWL** — nf-core/meerpipe, LINC (HPC + Singularity). **DALiuGE** — SKA SDP.

**Airflow + Podman for radio astronomy: zero evidence exists.** Airflow is proven for radio
(CHIME/FRB) but always with Docker/Swarm; Podman appears in astronomy only via the unrelated
"Astronomer.io" vendor. A **rootless, daemonless, offline, CPU-only Airflow-on-Podman template for
public radio data on a laptop is a genuinely unfilled, citable niche** — and is exactly this
project's automation layer, so the automation contribution is novel regardless of the science
domain chosen.

**No maintained "awesome-radio-astronomy" list exists** (itself a gap). De-facto indexes: ASCL
(ascl.net), Astropy affiliated, radio-astro-tools, RadioAstronomySoftwareGroup.

## Per-domain tooling gaps (feeding the gap analysis)

- **FRB stats:** no pip-installable burst-statistics library — Weibull wait-times, power-law energy
  fits, activity-epoch detection live as one-off scripts inside papers; FRBSTATS is web-only.
- **Pulsar timing:** no lightweight offline `par,tim → residual + red-noise PSD` one-liner (PINT is
  full-API; la_forge is downstream of enterprise; enterprise/discovery assume HPC/GPU); no
  pip-installable small-glitch scanner.
- **HI/continuum:** no turnkey tested tangent-point rotation-curve extractor; `jansky.sourcecounts`
  has OLS slope but no Poisson/ML-error log N–log S fit; no small NVSS×FIRST spectral-index matcher.
- **RFI/ML:** no reproducible SK-vs-SumThreshold comparison harness; no CPU-only pip-installable
  FRB/RFI candidate classifier (FETCH is TF/GPU); no canonical <100 MB labelled RFI benchmark.
- **SETI:** no CPU-only pip drift-search for small files (turboSETI heavy, hyperseti/BLISS GPU/C++);
  setigen injection is not paired with a lightweight detector or a shared injection-recovery score.
