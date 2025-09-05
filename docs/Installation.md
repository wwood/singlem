# Installing SingleM

There are several ways to install SingleM. On the assumption that a standard internet connection speed is available, each of these methods should take substantially less than 1 hour, hopefully less than 15 minutes.

SingleM can be installed through
[Bioconda](https://anaconda.org/bioconda/singlem). Lyrebird is installed by installing SingleM, and then using the `lyrebird` command instead of `singlem` e.g. use `lyrebird pipe ...` instead of `singlem pipe ...`.

```
conda create -c conda-forge -c bioconda --override-channels --name singlem singlem'>='0.20.0
```

Test if it works by running
```
conda activate singlem
singlem -h
lyrebird -h
```
After this, you'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data) or [lyrebird data](/tools/lyrebird_data).


## Installation via DockerHub
A docker image generated from the conda package is [available](https://hub.docker.com/r/wwood/singlem) on DockerHub. After installing Docker, run the following:
```
docker pull wwood/singlem:0.20.0
```

Test if it works by running
```
docker run wwood/singlem:0.20.0 -h
```

If the sequence data to be analyzed is in the current working directory, SingleM `pipe` can be used like so:
```
docker run -v `pwd`:`pwd` wwood/singlem:0.20.0 pipe --sequences \
    `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```
Two things to note:

1. The default SingleM reference data is included in the docker image, so running [singlem data](/tools/data) is not necessary for this installation method - you can jump straight to using `pipe`.
2. You should not specify `singlem` in the command line, as this is automatic. Simply use e.g. `docker run wwood/singlem:0.20.0 pipe -h`.

A similar procedure is true for Lyrebird, except that the docker image is [different](https://hub.docker.com/r/wwood/lyrebird), so you need to run:
```
docker pull wwood/lyrebird:0.20.0
docker run wwood/lyrebird:0.20.0 -h
```

## Installation via Singularity / Apptainer
SingleM can be installed via [Singularity](https://sylabs.io/singularity/) or [Apptainer](https://apptainer.org). After installing Singularity or Apptainer, run the following:
```
singularity pull docker://wwood/singlem:0.20.0
```

Test if it works by running
```
singularity run singlem_0.20.0.sif -h
```

If the sequence data to be analyzed is in the current working directory, SingleM `pipe` can be used like so:
```
singularity run -B `pwd`:`pwd` singlem_0.20.0.sif pipe --sequences \
    `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```
Two things to note:

1. The default SingleM reference data is included in the docker image, so running [singlem data](/tools/data) is not necessary for this installation method - you can jump straight to using `pipe`.
2. You should not specify `singlem` in the command line, as this is automatic. Simply use e.g. `singularity run singlem_0.20.0.sif pipe -h`.

A similar procedure is true for Lyrebird, except that the docker image is [different](https://hub.docker.com/r/wwood/lyrebird), so you need to run:
```
singularity pull docker://wwood/lyrebird:0.20.0
singularity run -B `pwd`:`pwd` lyrebird_0.20.0.sif -h
```

## Installation via PyPI
To install the Python libraries required for SingleM / Lyrebird:
```
pip install singlem
```
You may need super-user privileges.

SingleM also has several non-Python dependencies, which are documented in the `pixi.toml` file in the base directory of the repository. You'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data).


## Installation via Github, with pixi dependency management
SingleM can be installed from source together with its conda dependencies as follows. You will need to have [pixi](https://pixi.sh) installed first.

```
git clone https://github.com/wwood/singlem
cd singlem
pixi shell

singlem -h
lyrebird -h
```

After this, you'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data).

This procedure also installs lyrebird, but the reference metapackage is different. See [lyrebird data](/tools/lyrebird_data) for more information.

# Example data

To test the main subcommand of SingleM, [pipe](/tools/pipe) works, download a minimal dataset and generate a taxonomic profile like so:
```
wget 'https://github.com/wwood/singlem/raw/44e1f81404c12931742259088999290edbb271b3/test/data/methanobacteria/genomes/GCA_000309865.1_genomic.fna'
singlem pipe -1 GCA_000309865.1_genomic.fna -p /dev/stdout
```

This should output a profile similar to the below. When tested, the `pipe` took a little less than 2 minutes:
```
sample  coverage        taxonomy
GCA_000309865.1_genomic 0.39    Root; d__Archaea; p__Methanobacteriota; c__Methanobacteria; o__Methanobacteriales; f__Methanobacteriaceae; g__Methanobacterium
08/06/2024 04:50:25 PM INFO: Finished condense
```

# Containerised SingleM installation examples

To ensure that the instructions here work, they have been tested in containerised environments. Logs of this procedure are available at [https://github.com/wwood/singlem-installation](https://github.com/wwood/singlem-installation).
