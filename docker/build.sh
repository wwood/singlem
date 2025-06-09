
#!/bin/bash -eo pipefail

# To execute this script, ensure that the databases are up to date in db/, and the tag has been pushed to GitHub. Then:
# pixi run bash ./build.sh

export SINGLEM_VERSION=`singlem --version`
export SINGLEM_DOCKER_VERSION=wwood/singlem:$SINGLEM_VERSION
export SINGLEM_DB_BASENAME=`basename $SINGLEM_METAPACKAGE_PATH`

# Check SINGLEM_DB_BASENAME is not an empty string, else croak
if [ -z "$SINGLEM_DB_BASENAME" ]; then
    echo "SINGLEM_METAPACKAGE_PATH must be set to the path of the metapackage to use"
    exit 1
fi
##### Lyrebird #####
export LYREBIRD_DOCKER_VERSION=wwood/lyrebird:$SINGLEM_VERSION
export LYREBIRD_DB_BASENAME=`basename $LYREBIRD_METAPACKAGE_PATH`

# Check LYREBIRD_DB_BASENAME is not an empty string, else croak
if [ -z "$LYREBIRD_DB_BASENAME" ]; then
    echo "LYREBIRD_METAPACKAGE_PATH must be set to the path of the metapackage to use"
    exit 1
fi


sed 's/SINGLEM_VERSION/'$SINGLEM_VERSION'/g; s/SINGLEM_DB_BASENAME/'$SINGLEM_DB_BASENAME'/g; s/SINGLEM_COMMAND/singlem/g; s/SINGLEM_ENV_VAR/SINGLEM_METAPACKAGE_PATH/g' Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $SINGLEM_DOCKER_VERSION . && \
docker run -v `pwd`:`pwd` $SINGLEM_DOCKER_VERSION pipe --sequences `pwd`/test.fna --otu-table /dev/stdout && \




sed 's/SINGLEM_VERSION/'$SINGLEM_VERSION'/g; s/SINGLEM_DB_BASENAME/'$LYREBIRD_DB_BASENAME'/g; s/SINGLEM_COMMAND/lyrebird/g; s/SINGLEM_ENV_VAR/LYREBIRD_METAPACKAGE_PATH/g' Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $LYREBIRD_DOCKER_VERSION . && \
docker run -v `pwd`:`pwd` $LYREBIRD_DOCKER_VERSION pipe --genome-fasta-file `pwd`/lambda_phage.fna --otu-table /dev/stdout --output-extras && \
\
\
echo "Seems good - now you just need to 'docker push $LYREBIRD_DOCKER_VERSION && docker push $SINGLEM_DOCKER_VERSION' to upload the images to Docker Hub"
