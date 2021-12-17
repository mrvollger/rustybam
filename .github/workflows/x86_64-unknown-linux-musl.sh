#!/bin/sh
#apt-get install --yes libgsl0-dev
#apt-get install --yes libgsl-dev
#apk add --no-cache pkg-config
apk add --no-cache gsl-dev openssl
wget -qO- https:/micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba create -n rbt bcftools starcode -c conda-forge -c bioconda
yes | ./bin/micromamba shell init -s bash -p /home/gitpod/micromamba
