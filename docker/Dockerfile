FROM continuumio/miniconda3:latest

# SingleM only works with python-3.6 at the moment, because of dependency issues (conda fails to resolve)
RUN conda install -c bioconda -c conda-forge -c defaults python=3.6 singlem=0.13.2

WORKDIR /data
ENTRYPOINT ["singlem"]
