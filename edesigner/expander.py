# -*- coding: utf-8 -*-
# expand_designs 2022
# Jose Alfredo Martin

__author__ = 'Alfredo Martin'
__version__ = 'expander.v12.0.0'


import _pickle as pic
import os
import sys
import argparse
from classes.parallel import starmap_parallel
from classes.design import Design
from classes.bbt import BBT
from classes.parameter_reader import Parameters


def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description = """Expander is a script used to expand a single file containing
    eDESINGs into a series of files containing eDESINGs icorporating an additional cycle. It paralelizes the work 
    withing cores in the same node. It can be called from qsub so it paralelizes the work in a cluster, or from 
    a loop if it eDESIGNER is used in a single machine. This is not a script to be called directly by the end-user.
    """)
    parser.add_argument('-wF', '--wfolder',
                        help="""Working Folder where all files are stored""",
                        type=str,
                        required=True)
    parser.add_argument('-a', '--agent',
                        help="""index of this agent. It is used to report files with differnet names when used with 
                        qsub""",
                        type=int,
                        required=True)
    parser.add_argument('-c', '--cycle',
                        help="""index of this expansion cycle""",
                        type=int,
                        required=True)
    parser.add_argument('-m', '--max_designs',
                        help="""maximum number of designs to be included in one output file""",
                        type=int,
                        required=True)
    parser.add_argument('-df', '--designs_file',
                        help="""name of file containing the designs to be expanded. It is a plckled file with binary 
                        format with one design in each record. The file must be in --wfolder""",
                        required=True)
    parser.add_argument('-af', '--args_file',
                        help="""name of the file containing the arguments used to expand the designs. It is a plckled 
                        file  with binary format with a tuple with all the arguments in the first record. The file must 
                        be in --wfolder""",
                        type=str,
                        required=True)
    parser.add_argument('-nc', '--n_cycles',
                        help="""Total number of cycles """,
                        type=int,
                        required=True)

    args = parser.parse_args()
    assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
    assert os.path.isfile(os.path.join(args.wfolder, args.designs_file)), f'{os.path.join(args.wfolder, args.designs_file)} does not exist'
    assert os.path.isfile(os.path.join(args.wfolder, args.args_file)), f'{os.path.join(args.wfolder, args.args_file)} does not exist'
    return args


def expand_designs(previous_designs, BBTs, indexes, reaction, reaction_indexes, hp_indexes, available_reaction_input,
                   reaction_input, reaction_output, deprotection, deprotection_indexes, available_deprotection_input,
                   deprotection_input, deprotection_output, par, fg):
    """This function generates an expansion of a list of designs with additional arguments:
    previous_designs : list of instances of Design class
    BBTs : list of instances of BBT class
    indexes : list of int containing
    reaction : instance of Paramters class (reaction)
    reaction_indexes : list of indexes for all available reactions
    hp: indexes: list of indexes of BBTs that are headpieces
    available_reaction_input : list of tuples containing the input FGs for all available reactions
    reaction_input : list of tuples containing the input FGs for all reactions
    reaction_output : list of tuples containing the output FGs for all reactions
    deprotection : instance of Parameters class (deprotection)
    deprotection_indexes : list of indexes for all available deprotections
    available_deprotection_input : ist of tuples containing the output FGs for all available deprotections
    deprotection_input : list of tuples containing the output FGs for all deprotections
    deprotection_output : list of tuples containing the output FGs for all deprotections
    par : instance of Parameters class (par)
    fg : instance of Parameters class (fg)
    returns result : list of instances of Design class"""
    result = []
    for design in previous_designs:
        result += design.add_cycle(BBTs, indexes,  reaction, reaction_indexes, available_reaction_input, reaction_input,
                                   reaction_output, deprotection, deprotection_indexes, available_deprotection_input,
                                   deprotection_input, deprotection_output, par, fg)
    return result


def expander(wfolder, cycle, agent, args_file, designs_file, maxdesigns, n_cycles):
    """This function generates an expansion of a list of designs with additional arguments:
    All arguments below and wrapped in a tuple named args:
    wfolder: wtr: working folder
    cycle: int: current cycle
    agent: int: number of this agent
    args_file: str: file with picled tuple of args for expansion
    designs_file: str: file containing pickled previous designs (one per record)
    maxdesigns: int: maximum number of designs in each file
    returns None"""
    with open(os.path.join(wfolder, args_file), 'rb') as f:
        paralel_args = pic.load(f)
    BBTs, indexes, reaction, reaction_indexes, hp_indexes, available_reaction_input, reaction_input, reaction_output, \
    deprotection, deprotection_indexes, available_deprotection_input, deprotection_input, deprotection_output, \
    par, fg = paralel_args
    previous_designs = []
    with open(os.path.join(wfolder, designs_file), 'rb') as f:
        while True:
            try:
                previous_designs.append(pic.load(f))
            except:
                break  # we reached the end of the file
    designs = starmap_parallel(previous_designs, expand_designs, args=paralel_args, is_list=True)
    counter = -1  # counter for files
    nd = -1  # counter for designs
    for design in designs:
        if cycle == n_cycles:
            if design.min_natoms >= par.par['max_na_absolute']:
                continue
            design.add_lib_id(reaction, deprotection)
        nd += 1
        if nd % maxdesigns == 0:
            if nd > 0:
                f.close()
            counter += 1
            filename = os.path.join(wfolder, f'eDESIGNs_{cycle}_{agent}_{counter}.pic')
            f = open(filename, 'wb')
        pic.dump(design, f)
    f.close()
    os.remove(os.path.join(wfolder, designs_file))


if __name__ == '__main__':
    args = parse_args()
    expander(args.wfolder, args.cycle, args.agent, args.args_file, args.designs_file, args.max_designs, args.n_cycles)
