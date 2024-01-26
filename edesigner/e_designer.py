# -*- coding: utf-8 -*-
# eDESIGNER
# Jose Alfredo Martin

__version__ = 'e_designer.v.12.0.0'
__author__ = 'Alfredo Martin'

# Python modules
import os
import shutil
import sys
import _pickle as pic
import time
import argparse
import numpy as np
# Local modules
from classes.parameter_reader import Parameters
from classes.bbt import BBT
from classes.logger import Logger
from classes.design import Design
from classes.libdesign import LibDesign
from tqdm import tqdm



def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description = """e_designer:
    This script uses a set of imported parameters related to functional groups and reactions, 
    and a list of possible BBTs generated with the e_bbt_creator script and generates a 
    list of eDESIGNS according with certain specifications.
    An argument is required indicating which is the working folder where all eDESEGNER runs are stored. 
    A rin ID is also required as an argument. This is an id created and provided by eDESIGNER.
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
    parser.add_argument('-wF', '--wfolder',
                        help="""Working Folder""",
                        type=str,
                        required=True)
    parser.add_argument('-run', '--run_id',
                        help="""name of the run id""",
                        type=str,
                        required=True)
    parser.add_argument('-ots', '--override_time_stamp',
                        help="""When invoked, eDESIGNER run folder will be set to R000000. 
                        This is used only for testing""",
                        action='store_true')
    parser.add_argument('-pf', '--par_file',
                        help="""path file containing the par.par parameters file. If not provided it will be exracted 
                        from the paramenters folder corresponding to this run""",
                        type=str,
                        default=None)

    args = parser.parse_args()
    assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
    args.wfolder = os.path.abspath(args.wfolder)
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id)), f'{os.path.join(args.wfolder, args.run_id)} does not exist'
    if args.par_file is not None:
        assert os.path.isfile(args.par_file), f'{args.par_file} does not exist'
    return args

def create_parallel_args(par, BBTs, reaction, deprotection, fg, log):
    """This function creates a tuple with a set of arguments that will be the input of the expansion cycles
        par : instance of Parameters class (par)
    BBTs : list of instances of BBT class
    reaction : instance of Paramters class (reaction)
    deprotection : instance of Parameters class (deprotection)
    fg : instance of Parameters class (fg)
    log : instance of Logger class
    returns : tuple""
    """
    # Initiallization of indexes
    log.update('Creating reaction indexes...')
    # the following variables contain list of tuples containign the pairs or incoming or outcoming FGs to a reaction for all reactions
    reaction_input = [(reaction.par[i]['fg_input_on_off'][0], reaction.par[i]['fg_input_on_off'][1]) for i in range(len(reaction.par))]
    reaction_output = [(reaction.par[i]['fg_output_on_off'][0], reaction.par[i]['fg_output_on_off'][1]) for i in range(len(reaction.par))]
    deprotection_input = [(deprotection.par[i]['fg_input_on_off'][0], deprotection.par[i]['fg_input_on_off'][1]) for i in range(len(deprotection.par))]
    deprotection_output = [(deprotection.par[i]['fg_output_on_off'][0], deprotection.par[i]['fg_output_on_off'][1]) for i in range(len(deprotection.par))]
    # the following calculates the list of indexes for available reactions
    reaction_indexes = [i for i in range(1, len(reaction.par)) if (par.par['include_designs'].upper() == 'BOTH' or reaction.par[i]['production'])]
    available_reaction_input = [reaction_input[i] for i in reaction_indexes]
    deprotection_indexes = [i for i in range(1, len(deprotection.par)) if (par.par['include_designs'].upper() == 'BOTH' or deprotection.par[i]['production'])]
    available_deprotection_input = [deprotection_input[i] for i in deprotection_indexes]
    # the following calculate lists of indexes related with BBTs
    hp_indexes = [i for i in range(len(BBTs)) if BBTs[i].headpiece is not None] # this extracts the indexes of the BBTs that are headpieces
    # indexes = [i for i in range(len(BBTs)) if BBTs[i].n_compounds[-1] > 0]  # todo eliminate after testing
    indexes = [i for i in range(len(BBTs)) if BBTs[i].n_compounds.sum() > 0]  # this extracts the indexes of the BBTs that contain at least one compound
    # the next line packages all the parameters required to expand a cycle into a tuple
    paralel_args = (BBTs, indexes, reaction, reaction_indexes, hp_indexes, available_reaction_input, reaction_input,
                    reaction_output, deprotection, deprotection_indexes, available_deprotection_input,
                    deprotection_input, deprotection_output, par, fg)
    return paralel_args

def create_sge_script(path):
    """creates a sge script in the working folder
    path: str
    cores: int (number of cores per node)
    return: str: path to the sge_script"""
    sge_script = os.path.join(path, 'sge_script.sh')
    # todo include in this script that we want a number of cores per node
    with open(sge_script, 'w') as f:
        f.write('#$ -S /bin/bash\n')
        f.write('#$ -o arrayjob_script.out\n')
        f.write('#$ -e arrayjob_script.err\n')
        f.write('#$ -j y\n')
        f.write('#$ -cwd\n')
        f.write('eval $( head -${SGE_TASK_ID} ${1} | tail -1)\n')
    command = f'chmod u+x {sge_script}'
    os.system(command)
    return sge_script

def create_qsub_script(numjobs, sge_script, path, cores=4):
    """creates the qsub_script in the working folder
    numjobs: int: number of jobs to be performed
    sge_script: str: path to the sge_script file
    cores: int (number of cores per node)
    returns: str: path to the qsub script"""
    qsub_script = os.path.abspath(os.path.join(path, 'qsub_script.sh'))
    with open(qsub_script, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('current=$(dirname "$0")\n')
        f.write(f'NUMJOBS={numjobs}\n')
        f.write(f'qsub -V -t 1-{numjobs}:1 -pe smp {cores} -sync y {sge_script} $1\n')
        # f.write(f'qsub -V -t 1-{numjobs}:1 -sync y {sge_script} $1\n')
    command = f'chmod u+x {qsub_script}'
    os.system(command)
    return qsub_script

def create_designs(par, BBTs, reaction, deprotection, fg, log, RUNFOLDER, current):
    """This is the master function to create an eDESIGN set using multiple cpu parallelization within a single node.
    It starts generating list of indexes from BBTs and reactions objects that will speed up the loops because they
    have less members than the original objects. Then it creates the different cycles storing intermediate eDESIGNS in
    disk and using parallelization to expand each design into the next cycle.
    par : instance of Parameters class (par)
    BBTs : list of instances of BBT class
    reaction : instance of Paramters class (reaction)
    deprotection : instance of Parameters class (deprotection)
    fg : instance of Parameters class (fg)
    RUNFOLDER: str: path to the run folder for this e_designer run
    log : instance of Logger class
    current : folder of this script
    returns : None"""
    paralel_args = create_parallel_args(par, BBTs, reaction, deprotection, fg, log)
    with open(os.path.join(RUNFOLDER, "results", "args.pic"), 'wb') as f:
        pic.dump(paralel_args, f)
    BBTs, indexes, reaction, reaction_indexes, hp_indexes, available_reaction_input, reaction_input, reaction_output, \
    deprotection, deprotection_indexes, available_deprotection_input, deprotection_input, deprotection_output, \
    par, fg = paralel_args
    n_cycles = len(par.par['max_cycle_na'])

    # create designs
    log.update('Creating designs...')
    # Create the first set of designs containing only the available headpieces
    log.update('Creating headpieces...')
    designs = [Design(par, BBTs, index, n_cycles) for index in hp_indexes]  # creates the first level designs containing just the headpieces
    log.update('        Saving designs to file...')
    with open(os.path.join(RUNFOLDER, 'results', 'eDESIGNs_0_0_0.pic'), 'wb') as f:
        for design in designs:
            pic.dump(design, f)
    for cycle in range(n_cycles):
        log.update(f'cycle {cycle + 1}: Creating config file for cycle {cycle +1}...')
        expand_files = os.listdir(os.path.join(RUNFOLDER, 'results'))
        expand_files = [file for file in expand_files if file.startswith('eDESIGNs_') and file.endswith('.pic')]
        config_file = os.path.join(RUNFOLDER, 'results', 'qsub.config')
        with open(config_file, 'w') as f:
            for i, file in enumerate(expand_files):
                command = f'python {os.path.join(current, "expander.py")} -wF {os.path.join(RUNFOLDER, "results")}'
                command += f' -a {i} -c {cycle + 1} -m {par.par["designs_in_memory"]}'
                command += f' -df {file} -af args.pic -nc {n_cycles}\n'
                f.write(command)
        if par.par["hpc"] and cycle > 0 and len(expand_files) > 1:
            log.update(f'cycle {cycle + 1}: Expanding {len(expand_files)} files through hpc...')
            sge_script = create_sge_script(os.path.join(RUNFOLDER, "results"))
            qsub_script = create_qsub_script(len(expand_files),
                                             sge_script,
                                             os.path.abspath(os.path.join(RUNFOLDER, "results")),
                                             cores=par.par['hpc_cores'])
            os.system(f'{qsub_script} {config_file}')
        else:
            log.update(f'cycle {cycle + 1}: Expanding {len(expand_files)} files in a single node...')
            with open(config_file, 'r') as f:
                commands = f.readlines()
            commands = [command.strip() for command in commands]
            for i, command in enumerate(commands):
                log.update(f'    cycle {cycle +1}: expandig file {i+1} out of {len(commands)}')
                os.system(command)

def get_all_indexes(par, bblim, ndim):
    """ this function gets a numpy array with the number of atoms in a library with every posibility of nuber of atoms
    coming from each cycle. The final number of atoms is increased by the number of atoms accounted by the headpiece.
    par: instance of Parameters class (par)
    ndim: number of cycles in the library
    return numpy array"""
    na_array = np.array(list(range(bblim.par['max_bb_na'] + 1)))
    if ndim == 2:
        na_dist = na_array[:, None] + na_array[None, :]
        na_dist += par.par['headpiece_na']
    elif ndim == 3:
        na_dist = na_array[:, None, None] + na_array[None, :, None] + na_array[None, None, :]
        na_dist += par.par['headpiece_na']
    else:
        return None
    return na_dist

def create_libdesigns(args, na_dist, par, deprotection, BBTs, reaction, log, RUNFOLDER):
    """creates lib_designs"""
    log.update(f'Processing eDESIGNs into libDESIGNs...')
    files = os.listdir(os.path.join(RUNFOLDER, 'results'))
    files = [file for file in files if file.startswith('eDESIGNs_') and file.endswith('.pic')]
    lib_dict = dict()
    n_designs = 0
    for file in tqdm(files):
        designs = []
        with open(os.path.join(RUNFOLDER, 'results', file), 'rb') as f:
            while True:
                try:
                    designs.append(pic.load(f))
                    n_designs += 1
                except:
                    break
        os.remove(os.path.join(RUNFOLDER, 'results', file))
        for design in designs:
            if design.lib_id not in lib_dict.keys():
                lib_dict[design.lib_id] = LibDesign(par)
            lib_dict[design.lib_id].update_lib(design, reaction, deprotection, args.run_id, RUNNAME)
    log.update(f'{n_designs} eDESIGNs were generated')
    log.update('**** Curating libDESIGNs ****')
    lib_list = []
    count = -1
    for key in tqdm(lib_dict.keys()):
        if lib_dict[key] is not None:
            lib_dict[key].validate_lib(BBTs, par, deprotection, na_dist)
            if not lib_dict[key].eliminate:
                count += 1
                lib_dict[key].id = count
                lib_list.append(lib_dict[key])
    log.update(f'{len(lib_list)} libDESIGNs remaining after removing those that do not fulfill the criteria')
    log.update('Saving libraries to disk...')
    with open(os.path.join(RESULTSFOLDER, 'libDESIGNs.pic'), 'wb') as f:
        for item in lib_list:
            pic.dump(item, f)


if __name__ == '__main__':
    args = parse_args()
    current, _ = os.path.split(sys.argv[0])
    current = os.path.abspath(current)
    #initialization
    tic = time.time()
    error_found = False
    RESOURCESFOLDER = os.path.abspath(os.path.join(args.wfolder, args.run_id, 'resources'))
    timestamp = time.strftime('%Y%m%d%H%M', time.gmtime())
    if args.override_time_stamp:
        RUNNAME = 'ED000000'
    else:
        RUNNAME = 'ED' + timestamp
    log = Logger(os.path.join(args.wfolder, args.run_id, 'logs', RUNNAME + '.log'))
    log.update(__version__)
    log.update(__author__)
    log.update('**** CREATING FOLDER SYSTEM ****')
    RUNFOLDER = os.path.abspath(os.path.join(args.wfolder, args.run_id, RUNNAME))
    os.mkdir(RUNFOLDER)
    PARFOLDER = os.path.abspath(os.path.join(RUNFOLDER, 'resources'))
    os.mkdir(PARFOLDER)
    RESULTSFOLDER = os.path.abspath(os.path.join(RUNFOLDER, 'results'))
    os.mkdir(RESULTSFOLDER)
    if args.par_file is None:
        shutil.copy(os.path.join(RESOURCESFOLDER, 'par.par'), os.path.join(PARFOLDER, 'par.par'))
    else:
        shutil.copy(args.par_file, os.path.join(PARFOLDER, 'par.par'))
    log.update('**** READING PARAMETERS ****')
    bblim = Parameters(os.path.join(RESOURCESFOLDER, 'bblim.par'), fsource='dict', how='to_dict', multiple=False)
    par = Parameters(os.path.join(PARFOLDER, 'par.par'), fsource='dict', how='to_dict', multiple=False)
    fg = Parameters(os.path.join(RESOURCESFOLDER, 'fg.par'), fsource='list', how='to_list', multiple=True)
    reaction = Parameters(os.path.join(RESOURCESFOLDER, 'reaction.par'), fsource='list', how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(RESOURCESFOLDER, 'enum_reaction.par'), fsource='list', how='to_list', multiple=True)
    deprotection = Parameters(os.path.join(RESOURCESFOLDER, 'deprotection.par'), fsource='list', how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(RESOURCESFOLDER, 'enum_deprotection.par'), fsource='list', how='to_list', multiple=True)
    headpieces = Parameters(os.path.join(RESOURCESFOLDER, 'headpieces.par'), fsource='list', how='to_list', multiple=True)
    for item in [par, bblim, fg, reaction, enum_reaction, deprotection, enum_deprotection, headpieces]:
        if item.success is None:
            error_found = True
            log.update('    Error(s) found while reading parameters: ' + item.path)
            for error in item.errors:
                log.update('        ' + error)

    if not error_found:
        # load BBTs object (list of BBT class instances)
        log.update('**** LOADING BBTs ****')
        with open(os.path.join(args.wfolder, args.run_id, 'results', 'BBTs.pic'), 'rb') as f:
            BBTs = pic.load(f)

    # Create designs
    log.update('**** CREATING eDESIGNs ****')
    create_designs(par, BBTs, reaction, deprotection, fg, log, RUNFOLDER, current)
    # Create lib_designs
    log.update('**** CREATING libDESIGNs ****')
    n_cycles = len(par.par['max_cycle_na'])
    na_dist = get_all_indexes(par, bblim, n_cycles)
    create_libdesigns(args, na_dist, par, deprotection, BBTs, reaction, log, RUNFOLDER)
    if os.path.isfile(os.path.join(RUNFOLDER, 'results', 'qsub_script.sh')):
        os.remove(os.path.join(RUNFOLDER, 'results', 'qsub_script.sh'))
    if os.path.isfile(os.path.join(RUNFOLDER, 'results', 'qsub.config')):
        os.remove(os.path.join(RUNFOLDER, 'results', 'qsub.config'))
    if os.path.isfile(os.path.join(RUNFOLDER, 'results', 'args.pic')):
        os.remove(os.path.join(RUNFOLDER, 'results', 'args.pic'))
    if os.path.isfile(os.path.join(RUNFOLDER, 'results', 'sge_script.sh')):
        os.remove(os.path.join(RUNFOLDER, 'results', 'sge_script.sh'))
    tac = time.time()
    log.update(f'Running time: {round((tac - tic) /60, 1)} min.')
    log.update(f'the run ID for this run is {args.run_id}')
    log.update(f'the eDESIGNER run name is {RUNNAME}')
    log.update('OK')

