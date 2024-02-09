FROM mambaorg/micromamba:1.5.1

LABEL version="1.0.0"
LABEL software="checkM2"
LABEL software.version="checkM2 v1.0.1"
LABEL description="Image containing checkM2 v1.0.1 (bioconda installation) and micromamba. Includes the checkM2 database."
LABEL maintainer="Marco Teixeira"
LABEL maintainer.email="mcarvalh@broadinstitute.org"
LABEL maintainer.institution="Broad Institute of MIT and Harvard"

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes;

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN mkdir /tmp/checkm2-db && \
    checkm2 database --download --path /tmp/checkm2-db/ && \
    export CHECKM2DB="/tmp/checkm2-db/"

ENV PATH="/opt/conda/bin:$PATH"
#ENTRYPOINT ["micromamba", "run", "-n", "base"]

