# Installing SingleM

There are several ways to install SingleM. 

<!-- ### Installation via conda
SingleM can be installed through [Bioconda](https://anaconda.org/bioconda/singlem):

```
conda create -c bioconda --name singlem singlem
``` 
After this, you'll also need to procure the reference data (the "metapackage"). See [singlem data](usage/data).
-->

## Installation via DockerHub
A docker image generated from the conda package is available on DockerHub. After installing Docker, run the following, replacing `[RELEASE_TAG]` with a tag from [https://hub.docker.com/r/wwood/singlem/tags](https://hub.docker.com/r/wwood/singlem/tags):
```
docker pull wwood/singlem:[RELEASE_TAG]
```
If the sequence data to be analyzed is in the current working directory, SingleM can be used like so:
```
docker run -v `pwd`:`pwd` wwood/singlem:[RELEASE_TAG] pipe --sequences `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 14
```
The default SingleM reference data is included in the docker image, so [singlem data](/tools/data) is not necessary for this installation method.

<!-- ### Installation via PyPI
To install the Python libraries required:
```
pip install singlem
```
You may need super-user privileges.

SingleM also has several non-Python dependencies, which are documented in the `singlem.yml` file in 

* [OrfM](https://github.com/wwood/OrfM) >= 0.2.0 
* [HMMER](http://hmmer.janelia.org/) >= 3.1b1 
* [mfqe](https://github.com/wwood/mfqe) >= 0.5.0
* [KronaTools](http://sourceforge.net/p/krona/home/krona/) >= 2.4
* [diamond](https://github.com/bbuchfink/diamond) > 2.0.11
* sra-tools
* sqlite
* cd-hit -->

## Installation via Github, with conda environment dependencies
SingleM can be installed from source together with its conda dependencies as follows.

```
git clone https://github.com/wwood/singlem # Or download some specific release
cd singlem
conda env create -n singlem -f singlem.yml
conda activate singlem
cd bin
export PATH=$PWD:$PATH
singlem -h
```

After this, you'll also need to procure the reference data (the "metapackage"). See [singlem data](/tools/data).