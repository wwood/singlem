
#!/bin/bash -eo pipefail

export SINGLEM_VERSION=wwood/singlem:1.0.0beta1

# Copy env
cp ../singlem.yml .

DOCKER_BUILDKIT=1 docker build -t $SINGLEM_VERSION .

docker run -v `pwd`:`pwd` $SINGLEM_VERSION pipe --sequences `pwd`/test.fna --otu-table /dev/stdout

echo "Remember to change version numbers and metapackage path in Dockerfile .. then 'docker push $SINGLEM_VERSION'"