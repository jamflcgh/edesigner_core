# -*- coding: utf-8 -*-
# find_libraries
# Jose Alfredo Martin


import os
from classes.libdesign import LibDesign
import _pickle as pic
from classes.parameter_reader import Parameters
import argparse
import copy

__author__ = 'Alfredo Martin 2023'
__version__ = 'find_libraries.v.12.0.0'


def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description="""finds all libraries containing one or more reactions and / or 
    deprotections. Alternatively searches for a library with a specific design id.""")
    parser.add_argument('-wF', '--wfolder',
                        help="""Working Folder""",
                        type=str,
                        required=True)
    parser.add_argument('-run', '--run_id',
                        help="""name of the run id""",
                        type=str,
                        required=True)
    parser.add_argument('-erun', '--ed_run_id',
                        help="""name of the eDESIGNER run""",
                        type=str,
                        default=None)
    parser.add_argument('-r', '--reaction',
                        help="""index of a required reaction. You can repeat this argument. All reactions would be 
                        required. If the same reaction is entered more than once it will find desings with teh 
                        reaction repeated as many times as entered.""",
                        type=int,
                        action='append')
    parser.add_argument('-d', '--deprotection',
                        help="""index of a required deprotection. You can repeat this argument. All deprotections would 
                        be required. If the same deprotection is entered more than once it will find desings with teh 
                        deprotection repeated as many times as entered.""",
                        type=int,
                        action='append')
    parser.add_argument('-lid', '--lib_id',
                        help="""library_id (numbers separated by underscores). It overrides --reactions and 
                        --deprotections.""",
                        type=str,
                        default=None)
    parser.add_argument('-lidx', '--lib_idx',
                        help="""library index. This overrides --lib_id and --reaction and --deprotection""",
                        type=int,
                        default=None)
    parser.add_argument('-n', '--n_libraries',
                        help="""How many libraries you want in the report. 0 means all the ones found. 
                        Default: 0""",
                        type=int,
                        default=0)


    args = parser.parse_args()
    assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
    args.wfolder = os.path.abspath(args.wfolder)
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id)), f'{os.path.join(args.wfolder, args.run_id)} does not exist'
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id, args.ed_run_id)), f'{os.path.join(args.wfolder, args.run_id, args.ed_nun_id)} does not exist'
    if args.lib_id is not None:
        args.lib_id = tuple([int(item) for item in args.lib_id.split('_')])
    return args

def print_design(design, enum_reaction, enum_deprotection):
    print('*************************************************')
    print(design.__dict__)
    print()
    design.print_summary_file(enum_reaction, enum_deprotection)
    print('*************************************************')

def main():
    args = parse_args()
    PARFOLDER = os.path.join(args.wfolder, args.run_id, 'resources')
    bbts_file = os.path.join(args.wfolder, args.run_id, 'results', 'BBTs.pic')
    libdesigns_file = os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'libDESIGNs.pic')
    multireaction = Parameters(os.path.join(PARFOLDER, 'multireaction.par'), fsource='list', how='to_list', multiple=True)
    preparations = Parameters(os.path.join(PARFOLDER, 'preparations.par'), fsource='list', how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(PARFOLDER, 'enum_deprotection.par'), fsource='list', how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(PARFOLDER, 'enum_reaction.par'), fsource='list', how='to_list',multiple=True)
    headpieces = Parameters(os.path.join(PARFOLDER, 'headpieces.par'), fsource='list', how='to_list', multiple=True)
    with open(bbts_file, 'rb') as f:
        BBTs = pic.load(f)
    headpieces_dict = {}
    for headpiece in headpieces.par:
        headpieces_dict[[BBT.BBT for BBT in BBTs].index(headpiece['bbt'])] = headpiece['smiles']
    designs = []
    with open(libdesigns_file, 'rb') as f:
        while True:
            try:
                design = pic.load(f)
                designs.append(design)
            except:
                break
    # search for lib_idx
    if args.lib_idx is not None and args.lib_idx < len(designs) and designs[args.lib_idx].id == args.lib_idx:
        print_design(designs[args.lib_idx], enum_reaction, enum_deprotection)
        return
    elif args.lib_idx is not None and args.lib_idx >= len(designs):
        print(f"ERROR::: --lib_idx {args.lib_idx} is out of range")
        return
    elif args.lib_idx is not None and designs[args.lib_idx].id != args.lib_idx:
        print(f"ERROR::: --lib_idx {args.lib_idx} does not match library id")
        return
    #search for lib_id
    if args.lib_id is not None:
        for design in designs:
            if args.lib_id == design.lib_id:
                print_design(design, enum_reaction, enum_deprotection)
                break
        else:
            print(f"ERROR::: --lib_id {args.lib_id} not found")
        return
    # search for reactions and deprotections ocurrence
    counter = 0
    for design in designs:
        reactions = copy.deepcopy(design.reactions)
        deprotections = copy.deepcopy(design.deprotections)
        keep = True
        if args.reaction is not None:
            for reaction in args.reaction:
                if reaction in reactions:
                    reactions.pop(reactions.index(reaction))
                else:
                    keep = False
                    break
        if args.deprotection is not None and keep:
            for deprotection in args.deprotection:
                if deprotection in deprotections:
                    deprotections.pop(deprotections.index(deprotection))
                else:
                    keep = False
                    break
        if keep:
            counter += 1
        if args.n_libraries > 0 and counter > args.n_libraries:
            keep = False
        if keep:
            print_design(design, enum_reaction, enum_deprotection)

if __name__ == '__main__':
    print(__version__)
    print(__author__)
    main()
