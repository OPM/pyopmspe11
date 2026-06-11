FROM ghcr.io/astral-sh/uv:0.11.15 AS uv
FROM openporousmedia/opmreleases:2026.04_amd64

USER root

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends software-properties-common libhdf5-dev git python3 libpython3-dev \
  && rm -rf /var/lib/apt/lists/*

COPY --from=uv /uv /usr/local/bin/uv

ENV UV_PYTHON_INSTALL_DIR=/opt/uv/python
ENV UV_CACHE_DIR=/opt/uv/cache


RUN uv python install 3.14 && \
  uv venv /opt/venv --python 3.14 && \
  uv pip install --python /opt/venv/bin/python git+https://github.com/OPM/pyopmspe11.git

ENV PATH="/opt/venv/bin:${PATH}"

RUN chown -R opm:opm /opt/uv /opt/venv

USER opm
WORKDIR /shared_host
