# lrgsglib: Full environment with C backends and optional graph-tool
#
# Build:   docker build -t lrgsglib .
# Run:     docker run -it lrgsglib python
# Jupyter: docker run -p 8888:8888 lrgsglib jupyter notebook --ip=0.0.0.0 --allow-root

FROM condaforge/mambaforge:latest AS base

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY lrgsgenv.yml /tmp/lrgsgenv.yml
RUN mamba env create -f /tmp/lrgsgenv.yml && mamba clean -afy

SHELL ["conda", "run", "-n", "lrgsgnb", "/bin/bash", "-c"]

COPY . /workspace/lrgsglib
WORKDIR /workspace/lrgsglib

RUN pip install -e . --no-build-isolation

RUN python -c "import lrgsglib; print('lrgsglib OK')"
RUN python -c "import graph_tool; print('graph-tool OK')" || true

CMD ["conda", "run", "-n", "lrgsgnb", "python"]
