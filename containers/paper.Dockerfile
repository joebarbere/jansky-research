# Reproducible LaTeX build for the paper, using tectonic (a single static binary that
# fetches the exact packages it needs from CTAN and caches them — far lighter than a full
# TeX Live install and hermetic once the cache is warm).
FROM docker.io/library/debian:bookworm-slim

ARG TECTONIC_VERSION=0.15.0
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl ca-certificates \
    && curl -fsSL \
       "https://github.com/tectonic-typesetting/tectonic/releases/download/tectonic%40${TECTONIC_VERSION}/tectonic-${TECTONIC_VERSION}-x86_64-unknown-linux-musl.tar.gz" \
       | tar -xz -C /usr/local/bin tectonic \
    && apt-get purge -y curl && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

WORKDIR /paper
# `make paper` mounts paper/ here; tectonic resolves \input{generated/macros} and the figures.
CMD ["tectonic", "main.tex"]
