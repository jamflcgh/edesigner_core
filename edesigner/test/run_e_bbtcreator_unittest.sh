#!/bin/bash

current=$(dirname "$0")
current=$(realpath ${current})
if [ -d ../../eDESIGNER_venv ]
then
  source ${current}/../../eDESIGNER_venv/bin/activate
else
  echo "WARNING: venv not installed (run install.sh at the first level of the repo to install the environment)."
  echo "Using the current active environment"
fi

export EDESIGNER_FOLDER=${current}/..
export EDESIGNER_PARFOLDER=${current}/../resources
export EDESIGNER_TEST_FOLDER=${current}
export EDESIGNER_PREPS=${current}/../preparations
export EDESIGNER_QUERIES=${current}/../queries
export DEPROTECTION_FOLDER=${current}/../deprotections
export PYTHONPATH=${PYTHONPATH}:${current}/..
export PYTHONPATH=${PYTHONPATH}:${current}/../classes
export PYTHONPATH=${PYTHONPATH}:${current}


python -m unittest -v test_e_bbtcreator

