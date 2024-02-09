FROM condaforge/mambaforge:23.3.1-1

LABEL version="1.0.0"
LABEL software="StrainGE"
LABEL software.version="StrainGE v1.3.7"
LABEL description="Image containing StrainGE v1.3.7 (bioconda installation)."
LABEL maintainer="Marco Teixeira"
LABEL maintainer.email="mcarvalh@broadinstitute.org"
LABEL maintainer.institution="Broad Institute of MIT and Harvard"

RUN conda config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge; \
    conda create -n strainge python=3.8 strainge pandas=1.5.3

SHELL ["conda", "run", "-n", "strainge", "/bin/bash", "-c"]

