#!/bin/bash
source /etc/cluster-setup.sh  # look if this is required
current=$(dirname "$0")
current=$(realpath ${current})
if [ -d ${current}/../eDESIGNER_venv ]
then
  source ${current}/../eDESIGNER_venv/bin/activate
else
  echo "WARNING: venv not installed (run install.sh at the first level of the repo to install the environment)."
  echo "Using the current active environment"
fi

export EDESIGNER_FOLDER=${current}
export EDESIGNER_PARFOLDER=${current}/resources
export EDESIGNER_TEST_FOLDER=${current}/test
export EDESIGNER_PREPS=${current}/preparations
export EDESIGNER_QUERIES=${current}/queries
export DEPROTECTION_FOLDER=${current}/deprotections
export PYTHONPATH=${PYTHONPATH}:${current}
export PYTHONPATH=${PYTHONPATH}:${current}/classes
export LILLYMOL_EXECUTABLES=$(head -1 ${current}/../lillymol.path)/bin/Linux


if [ -d ${LILLYMOL_EXECUTABLES} ]
then
  echo "using ${LILLYMOL_EXECUTABLES} as a source of LillyMol executables"
else
  echo "using GC3TK as a source of LillyMol executables"
fi

python ${current}/e_bbt_creator.py "$@"
