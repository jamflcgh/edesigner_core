# -*- coding: utf-8 -*-
# lib_designer
# Jose Alfredo Martin

# Python modules
import os
import pickle as pic
import time
import argparse
# Local modules
from classes.bbt import BBT
from classes.design import Design
from classes.libdesign import LibDesign
from classes.libdesign import get_all_indexes
from classes.logger import Logger
from classes.parameter_reader import Parameters

version = 'lib_designer.v.9.0.0'

if __name__ == '__main__':
    # Arg parser
    parser = argparse.ArgumentParser(description="""lib:designer: 
    This script reads a list of e_designs that were created by e_designer 
    and are contained in a file, and then it transforms them in lib_desings (library designs) 
    It requires a parameter (wfolder) passed on command line which indicates the working folder 
    and also a set of parameters imported using the Parameters: 
    The path parameters contain information about the names of folders where the current 
    run of eDESIGNER is run and files are storeed so the files red are the corect ones and the 
    files written have the correct name and are stored in the correct folder. 
    The par parameters are parameters that define how the current run is done and will be reported at the
    begining of the log file for each run. If the run name or other parameters in path are not changed the
    files are overwritten.
    The script creates a object named global_indexes which is central to select the appropriate
    building blocks that goes into each library so the number of atoms profile matches the one
    set in the parameters. The object is stored as a global variable.
    Designs created with e_designer are red from a file and classified into libraries that contain
    the same library id. While doing this BBTs from different designs being used in each cycle are stored
    within the library but repeated ones are eliminated. Once the libraries are set the number of building
    blocks from each BBT are selected in such a way that the final distribution of number of atoms in the
    library meets the requirements and at the same time the library produces the maximum number of compounds.
    Once this distribution is set the libraries that do not contain enough compounds are eliminated.
    """)
    parser.add_argument('-wf', '--wfolder', help='Working Folder', type=str, default='./')
    args = parser.parse_args()

    #initialization
    error_found = False
    tic = time.time()
    path = Parameters(os.path.join(args.wfolder, 'resources', 'path.par'), fsource='dict', how='to_dict',
                      multiple=False)
    par = Parameters(os.path.join(args.wfolder, 'resources', 'par.par'), fsource='dict', how='to_dict', multiple=False)
    reaction = Parameters(os.path.join(args.wfolder, 'resources', 'reaction.par'), fsource='list', how='to_list',
                          multiple=True)
    deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'deprotection.par'), fsource='list', how='to_list',
                              multiple=True)
    token = '_'.join([path.par['Database_Run'], path.par['run']])
    log = Logger(os.path.join(args.wfolder, 'logs', token + '_lib_designer.log'))
    log.update(version)
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'path.par'), 'path parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'par.par'), 'par parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'reaction.par'), 'reaction parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'deprotection.par'), 'deprotection parameters')
    for item in [path, par, reaction, deprotection]:
        if item.success is None:
            error_found = True
            log.update('    Error(s) found while reading parameters: ' + item.path)
            for error in item.errors:
                log.update('        ' + error)

    # Main body of the script
    if not error_found:
        n_cycles = len(par.par['max_cycle_na'])
        global_indexes = get_all_indexes(par, n_cycles)
        # load BBTs object (list of BBT class instances)
        log.update('Loading BBTs...')
        with open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_BBTs.pic'), 'rb') as f:
            BBTs = pic.load(f)
        log.update('Loading dictionary of libraries...')
        with open(os.path.join(args.wfolder, 'data', token +  '_lib_id_list.pic'), 'rb') as f:
            lib_id_list = pic.load(f)
        log.update(f'Creating libraries for {len(lib_id_list)} library ids')
        lib_dict = {item : LibDesign(par) for item in lib_id_list} # this line creates the working dictionary of libraries
        log.update('Reading designs...')
        with open(os.path.join(args.wfolder, 'data', token + '_' + str(n_cycles) + '.pic'), 'rb') as f:
            continuar = True
            while continuar:
                try:
                    design = pic.load(f)
                except:
                    continuar = False # we have reached the end of the file
                else:
                    try:
                        lib_dict[design.lib_id].update_lib(design, reaction, deprotection)
                    except KeyError as error:
                        lib_dict[design.lib_id].update_lib(design, reaction, deprotection)
                        log.update(f'{design.id} with lib_id: {design.lib_id} is not in lib_dict')
                        error_found = True
                        continuar = False
    if not error_found:
        log.update('Processing libraries...')
        lib_list = [lib_dict[key] for key in lib_dict.keys() if lib_dict[key].lib_id is not None] # transform the lib dictionary in a lib list
        log.update(f'Processing {len(lib_list)} libraries...')
        for i, library in enumerate(lib_list): # mark the libraries not fulfilling the criteria as eliminable
            lib_list[i].validate_lib(BBTs, par, deprotection, global_indexes)
        lib_list = [library for library in lib_list if not library.eliminate] # eliminate libs marked as eliminable
        for i, library in enumerate(lib_list): # create lib id
            lib_list[i].id = i
        log.update(f'{len(lib_list)} libraries remaining after removing those that do not fulfill the criteria')
        log.update('Saving libraries to disk...')
        with open(os.path.join(args.wfolder, 'results', token + '_libDESIGNS.pic'), 'wb') as f:
            for item in lib_list:
                pic.dump(item, f)
        tac = time.time()
        log.update(f'Running time: {round((tac - tic) /60, 1)} min.')
        log.update('OK')
