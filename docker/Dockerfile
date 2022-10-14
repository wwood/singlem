FROM mambaorg/micromamba:0.17.0

# This dockerfile uses cached mounts, so to build use e.g.
# $ DOCKER_BUILDKIT=1 docker build .

# Don't need all of the dependencies of singlem, because only pipe is going to be run.
ADD singlem.yml /singlem.yml
RUN --mount=type=cache,target=/opt/conda/pkgs micromamba install -c conda-forge -c bioconda -c defaults -n base -f /singlem.yml && \
    micromamba clean --all --yes

######### UPDATE ME !!!!!!!!!! #########
# NOTE: The following 2 hashes should be changed in sync.
RUN rm -rf singlem && git clone https://github.com/wwood/singlem && cd singlem && git checkout 5d9042af
RUN echo '__version__ = "1.0.0beta1.5d9042af"' >singlem/singlem/version.py

# Remove bundled singlem packages
RUN rm -rfv singlem/singlem/data singlem/.git singlem/test singlem/appraise_plot.png

######### UPDATE ME !!!!!!!!!! #########
ADD S3.0.5.metapackage20220806.smpkg.zb /mpkg
ENV SINGLEM_METAPACKAGE_PATH=/mpkg

WORKDIR /data
ENTRYPOINT ["singlem"]

# CMD /bin/bash
