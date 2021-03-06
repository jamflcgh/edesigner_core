#!/bin/bash
source /etc/profile.d/z00_lmod.sh
module load python
[ ! -d ./resources ] && mkdir ./resources
[ ! -d ./comps ] && mkdir ./comps
[ ! -d ./data ] && mkdir ./data
[ ! -d ./logs ] && mkdir ./logs
[ ! -d ./results ] && mkdir ./results
python lib_designer.py "$@"
