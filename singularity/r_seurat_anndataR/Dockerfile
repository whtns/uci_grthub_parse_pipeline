FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL description="R environment with Seurat and anndataR for h5ad to RDS conversion"

RUN apt-get update -qq \
    && apt-get install -y --no-install-recommends \
        libhdf5-dev \
        libglpk-dev \
        libgsl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); BiocManager::install(c('Seurat', 'anndataR'), ask = FALSE, update = FALSE)"

ENV LC_ALL=C
ENV PATH=/usr/local/bin:${PATH}

ENTRYPOINT ["Rscript"]
