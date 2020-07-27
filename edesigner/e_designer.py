# -*- coding: utf-8 -*-
# eDESIGNER
# Jose Alfredo Martin

version = 'e_designer.v.9.0.0'

# Python modules
import os
import pickle as pic
import time
import argparse
# Local modules
from classes.parameter_reader import Parameters
from classes.bbt import BBT
from classes.design import Design
from classes.design import create_designs
from classes.design import expand_designs
from classes.logger import Logger


if __name__ == '__main__':

    # Arg parser
    parser = argparse.ArgumentParser(description = """e_designer:
    This script uses a set of importer parameters related to functional groups and reactions, 
    and a list of possible BBTs generated with the e_bbt_creator script and generates a 
    list of eDESIGNS according with certain specifications.
    An argument is required indicating which is the working folder. A folder structure within this 
    folder is required for the scritp to work properly, please see the documentation to 
    build the approtrapte folder structre.
    The path parameters contain information about the names of folders where the current
    run of eDESIGNER is run and files are storeed so the files red are the corect ones and the
    files written have the correct name and are stored in the correct folder.
    The par are parameters that define how the current run is done and will be reported at the begining
    of the log file for each run. If the run name or other parameters in path are not changed the files are
    overwritten.
    The reaction and deprotection parameters contain information about primary reaction (used to connect BBTs in
    designs. These parameters contain the functional group(s) on-DNA and off-DNA that participate in the reactions
    and the functional groups that are created as a result. Also contain the name of the reaction and a list of
    fucntional groups that are incompatible either on-DNA or off-DNA.
    The enum_raction and enum_deprotection are names of reactions that usually comprises a set of primary reactions
    and can be performed phisically in the lab toguether. These names are connected univocally to specific enumeration
    instructions.
    The headpieces parameters contain the BBTs that can be used as headpieces in a library.
    """)
    parser.add_argument('-wf', '--wfolder', help='Working Folder', type=str, default='./')
    args = parser.parse_args()


    #initialization
    tic = time.time()
    error_found = False
    path = Parameters(os.path.join(args.wfolder, 'resources', 'path.par'), fsource='dict', how='to_dict', multiple=False)
    par = Parameters(os.path.join(args.wfolder, 'resources', 'par.par'), fsource='dict', how='to_dict', multiple=False)
    fg = Parameters(os.path.join(args.wfolder, 'resources', 'fg.par'), fsource='list', how='to_list', multiple=True)
    reaction = Parameters(os.path.join(args.wfolder, 'resources', 'reaction.par'), fsource='list', how='to_list',
                          multiple=True)
    enum_reaction = Parameters(os.path.join(args.wfolder, 'resources', 'enum_reaction.par'), fsource='list',
                               how='to_list', multiple=True)
    deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'deprotection.par'), fsource='list',
                              how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'enum_deprotection.par'), fsource='list',
                                   how='to_list', multiple=True)
    headpieces = Parameters(os.path.join(args.wfolder, 'resources', 'headpieces.par'), fsource='list', how='to_list',
                            multiple=True)
    token = '_'.join([path.par['Database_Run'], path.par['run']])
    log = Logger(os.path.join(args.wfolder, 'logs', token + '_e_designer.log'))
    log.update(version)
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'path.par'), 'path parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'par.par'), 'par parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'fg.par'), 'fg parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'reaction.par'), 'reaction parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'enum_reaction.par'), 'enum_reaction parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'deprotection.par'), 'deprotection parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'enum_deprotection.par'),
                           'enum_deprotection parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'headpieces.par'), 'headpieces parameters')
    for item in [path, par, fg, reaction, enum_reaction, deprotection, enum_deprotection, headpieces]:
        if item.success is None:
            error_found = True
            log.update('    Error(s) found while reading parameters: ' + item.path)
            for error in item.errors:
                log.update('        ' + error)

    if not error_found:
        # load BBTs object (list of BBT class instances)
        log.update('Loading BBTs...')
        with open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_BBTs.pic'), 'rb') as f:
            BBTs = pic.load(f)

    # Create designs
    create_designs(par, BBTs, reaction, deprotection, fg, log, args, path)

    tac = time.time()
    log.update(f'Running time: {round((tac-tic)/60, 1)} min.')
    log.update('OK')

