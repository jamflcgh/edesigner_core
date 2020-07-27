#!/bin/bash
source /etc/profile.d/z00_lmod.sh
module load python
module load rdkit
[ ! -d ./resources ] && mkdir ./resources
[ ! -d ./comps ] && mkdir ./comps
[ ! -d ./data ] && mkdir ./data
[ ! -d ./logs ] && mkdir ./logs
[ ! -d ./results ] && mkdir ./results
python bbt_analyzer.py "$@"
