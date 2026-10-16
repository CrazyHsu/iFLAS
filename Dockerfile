# iFLAS core image (Python 3.11 + aligners + iFLAS).
# The cupcake, PacBio and R/rMATS environments (environment-cupcake.yml, environment-pacbio.yml,
# environment-diffsplice.yml) are not included; create them in the image with micromamba if needed.
#
#   docker build -t iflas .
#   docker run --rm -v "$PWD":/work iflas init-ref -f genome.fa --gtf genes.gtf -o ref/
FROM mambaorg/micromamba:1.5-jammy

# 1) core conda environment (installed into base for a simple entrypoint)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y

# 2) iFLAS itself (brings pymupdf etc. via pip)
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/iflas
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir -e /opt/iflas

WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "iflas"]
CMD ["--help"]
