#!/bin/sh
sudo apt-get install --yes libgsl0-dev
sudo apt-get install libgsl-dev

#wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
#./bin/micromamba create -n rbt bcftools starcode -c conda-forge -c bioconda
#yes | ./bin/micromamba shell init -s bash -p /home/gitpod/micromamba
