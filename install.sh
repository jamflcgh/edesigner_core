#!/bin/bash
current=$(dirname "$0")
wdir=$(pwd)
if [ $# -gt 2 ]
then
  echo "wrong number of arguments. Only one optional argument (LillyMol path) is allowed"
  exit
fi
if [ $# -eq 2 ]
then
  echo "$(realpath ${1})" > ${current}/lillymol.path
else
  if [ -d ${current}/../LillyMolPrivate ]
  then
    echo "$(realpath ${current}/../LillyMolPrivate)" > ${current}/lillymol.path
  fi
  if [ -d ${current}/../LillyMol ]
  then
    echo "$(realpath ${current}/../LillyMol)" > ${current}/lillymol.path
  fi
fi
LILLYMOL_REPO=$(head -1 ${current}/lillymol.path)

if [ ! -d ${current}/eDESIGNER_venv ]
then
  python3.8 -m venv ${current}/eDESIGNER_venv
  source ${current}/eDESIGNER_venv/bin/activate
  pip install -r ${current}/requirements.txt
else
  source ${current}/eDESIGNER_venv/bin/activate || exit
fi
cd ${LILLYMOL_REPO}
module load gcc12 || echo "WARNING::: unable to load gcc12, please ensure gcc12 is available"
module load bazelisk || echo "WARNING::: unable to load bazelisk, please ensure bazelisk is available"
make
deactivate
cd ${wdir}

