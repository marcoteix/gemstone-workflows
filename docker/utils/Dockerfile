FROM mambaorg/micromamba:1.5.1

LABEL version="1.0.0"
LABEL software="GEMSTONE utility scripts"
LABEL software.version="v1.0.0"
LABEL description="Image containing Python, Micromamba, and utility scripts for GEMSTONE."
LABEL maintainer="Marco Teixeira"
LABEL maintainer.email="mcarvalh@broadinstitute.org"
LABEL maintainer.institution="Broad Institute of MIT and Harvard"

COPY --chown=$MAMBA_USER:$MAMBA_USER requirements.yaml /tmp/requirements.yaml

RUN micromamba config append channels bioconda
RUN micromamba config append channels conda-forge

RUN micromamba install -y -n base -f /tmp/requirements.yaml && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER scripts /tmp

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV PATH="/opt/conda/bin:$PATH"