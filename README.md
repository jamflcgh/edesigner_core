
# LRL_edesigner_core

## Overview

eDESIGNER is an algorithm used to enumerate comprehensively DNA encoded library designs and select the ones giving a maximum number of library members that fulfill at the same time a pre-defined heavy atom distribution profile.

This is a new version of the original eDESIGNER. A detailed explanation of the original eDESIGNER can be found in Martín, A., Nicolaou, C.A. & Toledo, M.A. Navigating the DNA encoded libraries chemical space. Commun Chem 3, 127 (2020). https://doi.org/10.1038/s42004-020-00374-1. For explanation of differences and new features see the docs [what_is_new](./docs/what_is_new.md).

eDESIGNER is a python package which run under Linux with the following dependencies: python 3.8, rdkit 2022.09.1, tqdm 4.65.0, pandas 1.1.5, numpy 1.18.5 and LillyMol v7.0 (https://github.com/EliLillyCo/LillyMol). 
LillyMol has in turn the following dependencies: absl-py 2.0.0, pybind11 2.11.1 and protobuf 4.24.4.
LillyMol compilation requires gcc12 and bazelisk.
User is encouraged to read the installation guidelines of LillMol at https://github.com/EliLillyCo/LillyMol

## Installation

To install eDESIGNER use the following intructions:

1. cd to a folder that we will name in the future {Basefolder}
1. git clone {this_repo}
1. git clone https://github.com/EliLillyCo/LillyMol.git

If you do not use module load as method to load modules then load gcc12 and bazelisk manually, otherwise the install script will load modules for you.
Load python 3.8 (I have used 3.8.5 in my tests).
If you have downloaded LillyMol in a folder different than {Basefolder} then run step 4 else run step 5.

4. bash {Basefolder}/eDESIGNER_core/install.sh path_to_your_lillymol_repo
4. bash {Basefolder}/eDESIGNER_core/install.sh

This creates an environment with all python requirements and compiles Lillymol # modify the eDESIGNER_core with the name of the public repo. Process takes 1-2 hours.

In order to test that everything works you can run the following script.

6. bash {Basefolder}/eDESIGNER_core/edesigner/test/run_all_tests.sh

The testing will execute an eDESIGNER run with a set of BBs and standard parameters included in the repo and generate a set of enumerations of libraries. The process takes about 2 hours.

If you have run the script using {Basefolder} as your working directory, the testing creates two files: {Basefolder}/enum_failed and {Basefolder}/enum_incomplete and a folder named {Basefolder}/R000000. If the script is run from another folder replace {Basefolder} for your working directory. If any of the files or folder are not created or the file {Basefolder}/enum_failed is not empty, this is indicative of a failure.

The scripts were tested in an HPC environment running Red Hat Enterprise Linux Server v.7.9 .

## Usage overview

There are two main python scritps in the package: e_bbt_creator.py and e_designer.py. These scripts are run through bash wrapper scripts named: run_e_bbt_creator.sh and run_e_designer.sh. Typical line command examples would be:

bash {Basefolder}/LRL_eDESIGNER_core/edesigner/run_e_bbt_creator.sh

This will create a folder in your working directory of the form R{timestamp} which is the run ID. The simplest run for eDESIGNER would then be:

bash {Basefolder}/LRL_eDESIGNER_core/edesigner/run_e_designer.sh -run R{timestamp}

This will create another folder within the run ID folder of the form ED{timestamp} which is the eDESIGNER run ID. This key is used then to operate with the eDESIGNER output by other scripts.

It is possible to run different eDESIGNER runs in the same e_bbt_creator run (for example 2 or 3 cycle libraries). For each, run a different parameters file must be input. Each run generates a differernt eDESIGNER ID.

Once an eDESIGNER run exists different operations can be performed on it. The most typical are running an enumeration of a specific library with the script run_enumerate_libraries.sh or finding libraries using specific reactions with run_find_libraries.sh.

The output an the eDESIGNER run depends on the parameters and building block collecions used. The repo contains a set of parameters which can be used "as is" or edited. The most important file to edit is the file edesigner/resources/par.par which contains parameters of number of cycles and number of heavy atom limits.

Full instructions can be found in the docs. You can start with file [what_is_new](./docs/what_is_new.md)
