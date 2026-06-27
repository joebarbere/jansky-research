# Airflow image for the jansky-research pipeline.
# Build context is the repo root (see airflow/compose.yaml): the package is COPYed in and
# installed alongside Airflow so the DAG's in-process PythonOperators can import it.
FROM apache/airflow:2.9.3-python3.12

# jansky (the library) from its pinned git tag; it brings numpy/scipy/matplotlib/astropy/etc.
# Then jansky_research from the build context with --no-deps (its deps are already satisfied).
COPY --chown=airflow:root . /opt/jansky-research
RUN pip install --no-cache-dir "jansky @ git+https://github.com/joebarbere/jansky@v0.1.0" \
    && pip install --no-cache-dir --no-deps /opt/jansky-research
