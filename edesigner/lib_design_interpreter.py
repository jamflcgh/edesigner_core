# -*- coding: utf-8 -*-
# lib_design_interpreter
# Jose Alfredo Martin

version = 'lib_design_interpreter.v.9.0.0'

# Python modules
import os
import pickle as pic
import copy
import time
import argparse
# Local modules
from classes.bbt import BBT
from classes.design import Design
from classes.libdesign import LibDesign
from classes.logger import Logger
from classes.parameter_reader import Parameters


if __name__ == '__main__':

    # Args parser
    parser = argparse.ArgumentParser(description="""lib_design_interpreter:
    This script reads a list of libraries that were created by lib_designer 
    and are contained in a file and translate them into enumeration instructions that will be
    written to a configuration text file which will be then interpreted for compound enumeration.
    The configuration file, at high level, records the id of the library, instructions on how to
    create files containing the building blocks that will be used in each cycle, the sequence of
    reactions and deprotections used to build the library, and a final optional deprotections of
    containing remaining BOC, Fmoc or ester protecting groups.
    It requires a parameter (wfolder) passed on command line which indicates the working folder
    and also a set of parameters imported using the Parameters class.
    The path parameters contain information about the names of folders where the current
    run of eDESIGNER is run and files are stored so the files red are the corect ones and the
    files written have the correct name and are stored in the correct folder.
    The par are parameters that define how the current run is done and will be reported at the begining
    of the log file for each run. If the run name or other parameters in path are not changed the files are
    overwritten.
    The reaction and deprotection parameters contain information about primary reaction (used to connect BBTs in
    designs. These parameters contain the functional group(s) on-DNA and off-DNA that participate in the reactions
    and the functional groups that are created as a result. Also contain the name of the reaction and a list of
    functional groups that are incompatible either on-DNA or off-DNA.
    The enum_raction and enum_deprotection are names of reactions that usually comprises a set of primary reactions
    and can be performed physically in the lab toguether. These names are connected univocally to specific enumeration
    instructions.
    The headpieces parameters contain the BBTs that can be used as headpieces in a library.""")
    parser.add_argument('-wf', '--wfolder', help='Working Folder', type=str, default='./')
    args = parser.parse_args()

    # Initialization
    error_found = False
    tic = time.time()
    path = Parameters(os.path.join(args.wfolder, 'resources', 'path.par'), fsource='dict', how='to_dict',
                      multiple=False)
    par = Parameters(os.path.join(args.wfolder, 'resources', 'par.par'), fsource='dict', how='to_dict', multiple=False)
    fg = Parameters(os.path.join(args.wfolder, 'resources', 'fg.par'), fsource='list', how='to_list', multiple=True)
    reaction = Parameters(os.path.join(args.wfolder, 'resources', 'reaction.par'), fsource='list', how='to_list',
                          multiple=True)
    deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'deprotection.par'), fsource='list',
                              how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(args.wfolder, 'resources', 'enum_reaction.par'), fsource='list',
                               how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'enum_deprotection.par'), fsource='list',
                                   how='to_list', multiple=True)
    headpieces = Parameters(os.path.join(args.wfolder, 'resources', 'headpieces.par'), fsource='list', how='to_list',
                            multiple=True)
    token = '_'.join([path.par['Database_Run'], path.par['run']])
    log = Logger(os.path.join(args.wfolder, 'logs', token + '_lib_design_interpreter.log'))
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

    # Body of the script
    if not error_found:
        log.update('Loading BBTs...')
        with open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_BBTs.pic'), 'rb') as f:
            BBTs = pic.load(f)
        # now we create a dictionary with the smiles of headpieces and the key is the bbt index corresponding to that
        # headpiece
        headpieces_dict = {}
        for headpiece in headpieces.par:
            headpieces_dict[[BBT.BBT for BBT in BBTs].index(headpiece['bbt'])] = headpiece['smiles']
        log.update('Loading pickled designs to object...')
        with open(os.path.join(args.wfolder, 'results', token + '_libDESIGNS.pic'), 'rb') as f:
            lib_list = []
            while True:
                try:
                    lib_list.append(pic.load(f))
                except:  # reached the end of the file
                    break
        log.update('Creating config.txt file...')
        config_path = os.path.join(args.wfolder, 'results', token + '_config.txt')
        try:  # If the file already exists it has to be removed because we will use it in append mode
            os.remove(config_path)
        except:
            pass
        end_deprotection_enumeration_indexes = [item['enum_index'] for item in deprotection.par if item['end_deprotect']]
        end_deprotection_enumeration_indexes = list(set(end_deprotection_enumeration_indexes))
        for library in lib_list:
            library.update_translation_file(config_path, par, headpieces_dict, enum_reaction, enum_deprotection,
                                            deprotection, end_deprotection_enumeration_indexes)
        log.update(f'Translated {len(lib_list)} library designs')
        # Time and end the sript
        tac = time.time()
        log.update(f'Running time: {round((tac-tic)/60, 1)} min.')
        log.update('OK')    
