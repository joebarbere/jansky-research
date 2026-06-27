# Airflow image for the jansky-research pipeline.
# Build context is the repo root (see airflow/compose.yaml): the package is COPYed in and
# installed alongside Airflow so the DAG's in-process PythonOperators can import it.
#
# Fully-qualified image name: Fedora podman enforces short-name resolution and cannot prompt
# for a registry in a non-interactive (CI / background) build.
FROM docker.io/apache/airflow:2.9.3-python3.12

# The slim Airflow image has no usable git, which pip needs for the git+https jansky install.
USER root
RUN apt-get update && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*
USER airflow

# jansky (the library) from its pinned git tag; it brings numpy/scipy/matplotlib/astropy/etc.
# Then jansky_research from the build context with --no-deps (its deps are already satisfied).
COPY --chown=airflow:root . /opt/jansky-research
RUN pip install --no-cache-dir "jansky @ git+https://github.com/joebarbere/jansky@v0.1.0" \
    && pip install --no-cache-dir --no-deps /opt/jansky-research
