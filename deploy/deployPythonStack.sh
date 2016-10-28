#!/bin/bash

# Module system
. /mnt/software/Modules/current/init/bash


module purge
module load gcc/4.8.4 python/2.7.9 virtualenv zlib/1.2.5 openblas/0.2.14

set -euo pipefail
set -x

#BAUHAUS_VE=/mnt/secondary/Share/MLG/Bauhaus/VE-Bauhaus
BAUHAUS_VE=/tmp/VE-Bauhaus


rm -rf $BAUHAUS_VE

set +u
virtualenv $BAUHAUS_VE
source $BAUHAUS_VE/bin/activate
set -u

pip install -r requirements.txt
python setup.py install
