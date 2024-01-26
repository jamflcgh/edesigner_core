#!/bin/bash

current=$(dirname "$0")
wfolder=$(pwd)
${current}/../run_e_bbt_creator.sh -wF ${wfolder} -smi ${current}/dbs/enamine_test.smi -ots
${current}/../run_e_designer.sh -wF ${wfolder} -run R000000 -ots
${current}/../run_get_multireactions.sh -wF ${wfolder} -run R000000 -erun ED000000
for item in $(cat ${wfolder}/R000000/ED000000/results/enumeration_test_idxs.txt) 
do 
  ${current}/../run_enumerate_library.sh -wF ${wfolder} -run R000000 -erun ED000000 -lidx ${item} -n 100
done
${current}/run_check_enumerations.sh -eF ${wfolder}/R000000/ED000000/enumerations -oF ${wfolder}

