# Installing SingleM

There are several ways to install SingleM. 

## Installation via conda / mamba
SingleM can be installed through
[Bioconda](https://anaconda.org/bioconda/singlem). You can use the default conda
if you prefer but we recommend [mamba](https://mamba.readthedocs.io/) for faster
installation and better error messages.

```
mamba create -c bioconda -c conda-forge --name singlem singlem'>='0.18.0
``` 

Test if it works by running
```
conda activate singlem
singlem -h
```
After this, you'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data).


## Installation via DockerHub
A docker image generated from the conda package is [available](https://hub.docker.com/r/wwood/singlem) on DockerHub. After installing Docker, run the following:
```
docker pull wwood/singlem:0.18.0
```

Test if it works by running
```
docker run wwood/singlem:0.18.0 -h
```

If the sequence data to be analyzed is in the current working directory, SingleM `pipe` can be used like so:
```
docker run -v `pwd`:`pwd` wwood/singlem:0.18.0 pipe --sequences \
    `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```
Two things to note:

1. The default SingleM reference data is included in the docker image, so running [singlem data](/tools/data) is not necessary for this installation method - you can jump straight to using `pipe`.
2. You should not specify `singlem` in the command line, as this is automatic. Simply use e.g. `docker run wwood/singlem:0.18.0 pipe -h`.


## Installation via Singularity / Apptainer
SingleM can be installed via [Singularity](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org). After installing Singularity or Apptainer, run the following:
```
singularity pull docker://wwood/singlem:0.18.0
```

Test if it works by running
```
singularity run singlem_0.18.0.sif -h
```

If the sequence data to be analyzed is in the current working directory, SingleM `pipe` can be used like so:
```
singularity run -B `pwd`:`pwd` singlem_0.18.0.sif pipe --sequences \
    `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```
Two things to note:

1. The default SingleM reference data is included in the docker image, so running [singlem data](/tools/data) is not necessary for this installation method - you can jump straight to using `pipe`.
2. You should not specify `singlem` in the command line, as this is automatic. Simply use e.g. `singularity run singlem_0.18.0.sif pipe -h`.


## Installation via PyPI
To install the Python libraries required:
```
pip install singlem
```
You may need super-user privileges.

SingleM also has several non-Python dependencies, which are documented in the `singlem.yml` file in the base directory of the repository. You'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data).


## Installation via Github, with conda environment dependencies
SingleM can be installed from source together with its conda dependencies as follows.

```
git clone https://github.com/wwood/singlem
cd singlem
conda env create -n singlem -f singlem.yml
conda activate singlem
cd bin
export PATH=$PWD:$PATH
singlem -h
```

After this, you'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data).

# Containerised SingleM installation examples

To ensure that the instructions here work, they have been tested in containerised environments. Logs of this procedure are available at [https://github.com/wwood/singlem-installation](https://github.com/wwood/singlem-installation).
