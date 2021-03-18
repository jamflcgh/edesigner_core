#!/bin/bash
source /etc/profile.d/z00_lmod.sh
module load python
module load rdkit
export PYTHONPATH=$PYTHONPATH:../
python -m unittest

