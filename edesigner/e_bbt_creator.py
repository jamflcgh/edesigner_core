# -*- coding: utf-8 -*-
# e_bbt_creator.v.12.0.0
# Jose Alfredo Martin

__version__ = 'e_bbt_creator.v.12.0.0'
__author__ = 'Alfredo Martin'

# Python modules
import os
import shutil
import sys
import _pickle as pic
import time
import argparse
# External modules
from classes.bbt import BBT
from classes.logger import Logger
from classes.parameter_reader import Parameters
from classes.bb_reader import BBReader

# Functions definitions

def parse_args():
    parser = argparse.ArgumentParser(description="""e_bbt_creator:
        This script searches the files of annotated compounds and classifies all compounds
        into the building block type (BBT) they belong to, or discard them if they do not belong
        to a specific BBT.
        The script starts generating by comprehension all building block types as a combination of
        up to three compatible functional groups (FGs).
        Once the list of BBTs are set the compounds from different collections are classified and
        deduplicated, compounds which smiles string cannot be parsed with rdkit or compounds not belonging to
        any BBT are eliminated.
        The effective number of atoms (number of atoms minus the number of atoms that will be lost upon
        reaction) of each building block is calculated, and the building blocks corresponding to
        each BBT stored in a file.
        Two files per BBT are created, one containing all building blocks corresponding to that BBT and the
        other containing building blocks in internal collections only.
        The characteristics of each BBT are stored in an object. This includes number of BBs containing each
        possible number of atoms, functional groups corresponding to this BBT etc, and a list of all BBTs
        is stored in a file for further use. The name of the list is BBTs and is used throughout all scripts
        in eDESIGNER.
        """)
    parser.add_argument('-pF', '--parameters_folder',
                        help="""folder where the parameters reside, if not given then the parameters from the 
                        repo (test folder) will be used""",
                        type=str,
                        default=None)
    parser.add_argument('-wF', '--wfolder',
                        help="""working folder. Folder where the data is going to be stored. Inside this folder a 
                        timestamp folder will be created which will be unique for this run. If not passed the 
                        current working directory is used as the working folder""",
                        type=str,
                        default=None)
    parser.add_argument('-smi', '--smiles_file',
                        help="""It overrides the use of db.par. Path to a single smiles file containing all 
                        building blocks.""",
                        type=str,
                        default=None)
    parser.add_argument('-ots', '--override_time_stamp',
                        help="""When invoked, the run folder will be set to R000000. This is used only for testing""",
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        help="""When invoked additional information is printed in the console.""",
                        action='store_true')
    args = parser.parse_args()
    if args.wfolder is None:
        args.wfolder = os.path.abspath(os.getcwd())
    assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
    return args


def compatible(A, fg):
    """Function compatible has as argument a list with three indexes corresponding to FGs for one BBT
    and returns a bool indicating whether these FGs are compatible in the same BBT
    A : BBT (list of three int corresponding to a BBT)
    fg : instance of the Parameters class (functional group parameters)
    returns : resultado (bool)"""
    resultado = True
    for i in range(3):
        for j in range(i+1, 3):
            if A[j] != 0 and A[i] != 0 and A[j] in fg.par[A[i]]['self_incompatibility']:
                resultado=False
    return resultado

def intialization(wfolder):
    """The initiallization function establishes the folder structure for this run. Creates the run name which
    will be used in all the other scritps. Reads the required parameters and creates the log file
    wfolder: str: path to the folder where the run will be created
    returns: tuple"""
    # Initialization
    current, _ = os.path.split(sys.argv[0])
    error_found = False
    tic = time.time()
    timestamp = time.strftime('%Y%m%d%H%M', time.gmtime())
    if args.override_time_stamp:
        RUNNAME = 'R000000'
    else:
        RUNNAME = 'R' + timestamp
    RUNFOLDER = os.path.join(wfolder, RUNNAME)
    LOGFOLDER = os.path.join(RUNFOLDER, 'logs')
    COMPSFOLDER = os.path.join(RUNFOLDER, 'comps')
    RESULTSFOLDER = os.path.join(RUNFOLDER, 'results')
    RESOURCESFOLDER = os.path.join(RUNFOLDER, 'resources')
    os.mkdir(RUNFOLDER)
    os.mkdir(LOGFOLDER)
    os.mkdir(COMPSFOLDER)
    os.mkdir(RESULTSFOLDER)
    os.mkdir(RESOURCESFOLDER)
    if args.parameters_folder is not None:
        for file in os.listdir(args.parameters_folder):
            if file.endswith('par'):
                shutil.copy(os.paht.join(args.parameters_folder, file), os.path.join(RESOURCESFOLDER, file))
    else:  # if parameters are not provided we use the repo parameters
        for file in os.listdir(os.path.join(os.environ['EDESIGNER_FOLDER'], 'resources')):
            if file.endswith('par'):
                shutil.copy(os.path.join(os.environ['EDESIGNER_FOLDER'], 'resources', file),
                            os.path.join(RESOURCESFOLDER, file))
        PARFOLDER = RESOURCESFOLDER  # in this case the function is not imported
    log = Logger(os.path.join(LOGFOLDER, 'e_bbt_creator.log'), prepend_timestamp=True)
    dbpar = Parameters(os.path.join(PARFOLDER, 'db.par'), fsource='list', how='to_list', multiple=True)
    for i in range(len(dbpar.par)):
        dbpar.par[i]['filename'] = os.path.expandvars(dbpar.par[i]['filename'])
    bblim = Parameters(os.path.join(PARFOLDER, 'bblim.par'), fsource='dict', how='to_dict', multiple=False)
    par = Parameters(os.path.join(PARFOLDER, 'par.par'), fsource='dict', how='to_dict', multiple=False)
    fg = Parameters(os.path.join(PARFOLDER, 'fg.par'), fsource='list', how='to_list', multiple=True)
    calcfg = Parameters(os.path.join(PARFOLDER, 'calcfg.par'), fsource='list', how='to_list', multiple=True)
    antifg = Parameters(os.path.join(PARFOLDER, 'antifg.par'), fsource='list', how='to_list', multiple=True)
    headpieces = Parameters(os.path.join(PARFOLDER, 'headpieces.par'), fsource='list', how='to_list', multiple=True)
    for item in [dbpar, par, fg, calcfg, antifg, headpieces]:
        if item.success is None:
            error_found = True
            log.update('    Error(s) found while reading parameters: ' + item.path)
            for error in item.errors:
                log.update('        ' + error)
    return tic, log, dbpar, bblim, par, fg, calcfg, antifg, headpieces, RUNFOLDER, RUNNAME, RESULTSFOLDER

def generate_bbts(fg, headpieces, bblim, log):
    """fg: instance of Parameters class
    headpieces: instance of Parameters class
    bblim: instance of Parameters class
    log: instance of Logger class
    teturn: BBTs: list of instances of BBT class"""
    # Generates the list of all BBTs by comprehension ignoring incompatible BBTs but including [0, 0, 0]
    log.update('Generating BB types...')
    BBT_list = [[i, j, k] for i in range(len(fg.par)) for j in range(i, len(fg.par)) for k in range(j, len(fg.par)) if
                compatible([i, j, k], fg)]
    BBTs = [BBT(BBT=BBT_list[i], fg=fg, headpieces=headpieces, index=i, maxna=bblim.par['max_bb_na']) for i in range(len(BBT_list))]
    return BBTs

if __name__ == '__main__':

    #Main body of the script
    args = parse_args()
    tic, log, dbpar, bblim, par, fg, calcfg, antifg, headpieces, RUNFOLDER, RUNNAME, RESULTSFOLDER = intialization(args.wfolder)
    log.update(__version__)
    log.update(__author__)
    BBTs = generate_bbts(fg, headpieces, bblim, log)

    # read compound sets and get valid compounds within valid BBTs
    reader = BBReader(RESULTSFOLDER, RUNFOLDER, BBTs, dbpar, fg, antifg, calcfg, bblim, log,
                      smi=args.smiles_file, verbose=args.verbose, debug=False)
    reader.run()

    # time and end the program
    tac = time.time()
    log.update(f'Run time = {round((tac - tic) / 60.0, 1)} min')
    log.update(f'the run ID for this run is {RUNNAME}')
    log.update('OK')
