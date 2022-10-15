FROM continuumio/miniconda3:4.12.0

# This dockerfile uses cached mounts, so to build use e.g.
# $ DOCKER_BUILDKIT=1 docker build .

WORKDIR /tmp

# We cannot use micromamba because man-db doesn't appear to be available in any conda repo
ADD singlem.yml /singlem.yml
RUN --mount=type=cache,target=/opt/conda/pkgs conda install -c conda-forge mamba
RUN --mount=type=cache,target=/opt/conda/pkgs mamba env update -p /opt/conda git -f /singlem.yml
RUN --mount=type=cache,target=/opt/conda/pkgs mamba clean --all --yes
# Install man-db
RUN apt update && apt install -y man-db

######### UPDATE ME !!!!!!!!!! #########
# NOTE: The following 2 hashes should be changed in sync.
RUN rm -rf singlem && git clone https://github.com/wwood/singlem && cd singlem && git checkout v1.0.0beta1
# RUN echo '__version__ = "1.0.0beta1"' >singlem/singlem/version.py

# Remove bundled singlem packages
RUN rm -rfv singlem/singlem/data singlem/.git singlem/appraise_plot.png

# Download default metapackage
RUN mkdir -p singlem/db && cd singlem/db && /tmp/singlem/bin/singlem data --output-directory .

######### UPDATE ME !!!!!!!!!! #########
# RUN ls /tmp/singlem/db/S3.0.5.metapackage20220806.smpkg.zb && fail
ENV SINGLEM_METAPACKAGE_PATH=/tmp/singlem/db/S3.0.5.metapackage20220806.smpkg.zb

# Run some tests
RUN /tmp/singlem/bin/singlem --version
RUN /tmp/singlem/bin/singlem pipe --full-help 
# RUN cd singlem && python test/test_makedb_and_query.py

WORKDIR /data
ENTRYPOINT ["/tmp/singlem/bin/singlem"]

# CMD /bin/bash
