
#!/bin/bash -eo pipefail

export SINGLEM_VERSION=`../bin/singlem --version`
export SINGLEM_DOCKER_VERSION=wwood/singlem:$SINGLEM_VERSION

cp ../singlem.yml . && \
sed 's/SINGLEM_VERSION/'$SINGLEM_VERSION'/g' Dockerfile.in > Dockerfile && \
DOCKER_BUILDKIT=1 docker build -t $SINGLEM_DOCKER_VERSION . && \
docker run -v `pwd`:`pwd` $SINGLEM_DOCKER_VERSION pipe --sequences `pwd`/test.fna --otu-table /dev/stdout && \
echo "Seems good - now you just need to 'docker push $SINGLEM_DOCKER_VERSION'"
