FROM mambaorg/micromamba:0.17.0

# This dockerfile uses cached mounts, so to build use e.g.
# $ DOCKER_BUILDKIT=1 docker build .

# Don't need all of the dependencies of singlem, because only pipe is going to be run.
ADD singlem.yml /singlem.yml
RUN --mount=type=cache,target=/opt/conda/pkgs micromamba install -c conda-forge -c bioconda -c defaults -n base git -f /singlem.yml && \
    micromamba clean --all --yes

######### UPDATE ME !!!!!!!!!! #########
# NOTE: The following 2 hashes should be changed in sync.
RUN rm -rf singlem && git clone https://github.com/wwood/singlem && cd singlem && git checkout v1.0.0beta1
RUN echo '__version__ = "1.0.0beta1"' >singlem/singlem/version.py

# Remove bundled singlem packages
RUN rm -rfv singlem/singlem/data singlem/.git singlem/test singlem/appraise_plot.png

# Download default metapackage
RUN mkdir -p singlem/db && cd singlem/db && ls / && /tmp/singlem/bin/singlem data --output-directory .

######### UPDATE ME !!!!!!!!!! #########
# RUN ls /tmp/singlem/db/S3.0.5.metapackage20220806.smpkg.zb && fail
ENV SINGLEM_METAPACKAGE_PATH=/tmp/singlem/db/S3.0.5.metapackage20220806.smpkg.zb

WORKDIR /data
ENTRYPOINT ["/tmp/singlem/bin/singlem"]

# CMD /bin/bash
