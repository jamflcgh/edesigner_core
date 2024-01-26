# -*- coding: utf-8 -*-
# eDESIGNER
# Jose Alfredo Martin

__version__ = 'enumerate_library.v.12.0.0'
__author__ = 'Alfredo Martin'

# Python modules
import os
import shutil
import sys
import argparse
import time
# Local modules
from classes.parameter_reader import Parameters
from classes.libdesign import LibDesign
from classes.enumerator import Enumerator
from classes.bbt import BBT
import _pickle as pic


def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description = """enumerate_library reads a pickled libDESIGN and creates an 
    enumeration for it. Alternatie it costructs a libDESIGN from a lib_id string and a set of building blocks and 
    runs a enumeration for it
    """)
    # Defines whether library is coming from eDESIGNER or is external
    parser.add_argument('-u', '--user',
                        help="""Defines whether library is coming from eDESIGNER or is external. When not set 
                        --run_id, --ed_run_id and --wfolder must be set. Also either --lib_id or --lib_idx must be set.
                         If this parameter is not set then --efolder --lib_id and --bbs_file must be set""",
                        action='store_true')
    # Not from eDESIGNER
    parser.add_argument('-eF', '--efolder',
                        help="""working folder for enumerations. It is required when a eDESIGNER run is not passed. 
                        The folder must exist""",
                        type=str,
                        default=None)
    parser.add_argument('-lid', '--lib_id',
                        help="""library_id (numbers separated by underscores). It is required when a eDESIGNER run is 
                        not passed.""",
                        type=str,
                        default=None)
    parser.add_argument('-bbf', '--bbs_file',
                        help="""path of smiles file containing building blocks. If not provided, then the building 
                            blocks from eDESIGNER will be used. If provided, it must be repeated as many times as cycles 
                            in the library so each call sets up a file for a cycle (in order). A file for the Headpiece 
                            is not passed since the headpiece smiles is retrieved from the --lib_id argument.
                            provided""",
                        type=str,
                        action='append')
    parser.add_argument('-pF', '--parfolder',
                        help="""folder containing eDESIGNER parameters. If not passed, the default parameters in the 
                        repo are used.""",
                        type=str,
                        default=None)
    # from eDESIGNER

    parser.add_argument('-wF', '--wfolder',
                        help="""eDESIGNER working folder""",
                        type=str,
                        default=None)
    parser.add_argument('-run', '--run_id',
                        help="""run id of BBT creator""",
                        type=str,
                        default=None)
    parser.add_argument('-erun', '--ed_run_id',
                        help="""run id of eDESIGNER""",
                        type=str,
                        default=None)
    parser.add_argument('-lidx', '--lib_idx',
                        help="""library index. This overrides --lib_id""",
                        type=int,
                        default=None)

    # common
    parser.add_argument('-n', '--nmols',
                        help="""number of molecules to sample. If set to 0 then the whole library is enumerated""",
                        type=int,
                        default=0)
    parser.add_argument('-nc', '--nc_per_run',
                        help="""How many compounds can be enumerated in a single core.
                        Default: 400000""",
                        type=int,
                        default=400000)
    parser.add_argument('-wj', '--write_json',
                        help="""When invoked the script will create the json config file for enumeration and gather 
                        building blocks but it will not conduct the actual enumeration.""",
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        help="""When invoked the script will provide additional information in the standart output.""",
                        action='store_true')

    args = parser.parse_args()
    if args.user:
        assert args.efolder is not None, f'--efolder must be passed'
        assert os.path.isdir(args.efolder), f'{args.efolder} does not exist'
        assert args.lib_id is not None, 'lib_id or lib_idx must be provided'
        lib_id = tuple([int(item) for item in args.lib_id.split('_')])
        assert args.bbs_file is not None, '--bbs_file must be passed as many times as number of cycles'
        assert len(args.bbs_file) == lib_id[0], '--bbs_file must be passed as many times as number of cycles'
        if args.parfolder is not None:
            assert os.path.isdir(args.parfolder), f'{args.parfolder} does not exist'
    else:
        assert args.wfolder is not None, f'--wfolder must be passed'
        assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
        assert args.run_id is not None, f'--run_id must be passed'
        assert(os.path.isdir(os.path.join(args.wfolder, args.run_id))), f'{os.path.join(args.wfolder, args.run_id)} does not exist'
        assert args.ed_run_id is not None, f'--ed_run_id must be passed'
        assert (os.path.isdir(
            os.path.join(args.wfolder, args.run_id, args.ed_run_id))), f'{os.path.join(args.wfolder, args.run_id, args.ed_run_id)} does not exist'
        assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
        assert args.lib_idx is not None or args.lib_id is not None, 'either lib_id or lib_idx must be provided'
    return args

def main():
    """
    runs the enumerator
    returns: None:
    """
    args = parse_args()
    if not args.user:
        print('eDESIGNER enumeration')
        ld_file = os.path.join(args.wfolder, args.run_id, args.ed_run_id, "results", "libDESIGNs.pic")
        ld = []
        with open(ld_file, 'rb') as f:
            while True:
                try:
                    ld.append(pic.load(f))
                except:
                    break
        PARFOLDER = os.path.abspath(os.path.join(args.wfolder, args.run_id, "resources"))
        if args.lib_idx:
            for lib in ld:
                if lib.id == args.lib_idx:
                    lib_id = lib.lib_id
                    lib_idx = args.lib_idx
                    break
            else:
                print(f'Could not find lib_idx {args.lib_idx} in eDESIGNER.')
                sys.exit(1)
        else:
            lib_id = tuple([int(item) for item in args.lib_id.split('_')])
            for lib in ld:
                if lib.lib_id == lib_id:
                    lib_idx = lib.id
                    break
            else:
                print(f'Could not find lib_id {args.lib_id} in eDESIGNER.')
                sys.exit(1)
        bwfolder = os.path.join(args.wfolder, args.run_id, args.ed_run_id, "enumerations")
        if not os.path.isdir(bwfolder):
            os.mkdir(bwfolder)
        increment = 1
        enum_id = "EN" + str(lib_idx) + '_' + str(increment)
        while True:
            if os.path.isdir(os.path.join(bwfolder, enum_id)):
                increment += 1
                enum_id = "EN" + str(lib_idx) + '_' + str(increment)
            else:
                break
        wfolder = os.path.abspath(os.path.join(bwfolder, enum_id))
        os.mkdir(wfolder)
    else:
        print('User enumeration')
        lib_idx = None
        lib = None
        lib_id = tuple([int(item) for item in args.lib_id.split('_')])
        if args.parfolder is not None:
            PARFOLDER = os.path.abspath(args.parfolder)
        else:
            PARFOLDER = os.path.abspath(os.environ['EDESIGNER_PARFOLDER'])
            print(f'INFO::: PARFOLDER has been set to {PARFOLDER}')
        increment = 1
        enum_id = "ENUSER_" + str(increment)
        while True:
            if os.path.isdir(os.path.join(args.efolder, enum_id)):
                increment += 1
                enum_id = "ENUSER_" + str(increment)
            else:
                break
        wfolder = os.path.abspath(os.path.join(args.efolder, enum_id))
        os.mkdir(wfolder)
    # read parameters
    multireaction = Parameters(os.path.join(PARFOLDER, 'multireaction.par'), fsource='list', how='to_list', multiple=True)
    preparations = Parameters(os.path.join(PARFOLDER, 'preparations.par'), fsource='list', how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(PARFOLDER, 'enum_deprotection.par'), fsource='list', how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(PARFOLDER, 'enum_reaction.par'), fsource='list', how='to_list',multiple=True)
    headpieces = Parameters(os.path.join(PARFOLDER, 'headpieces.par'), fsource='list', how='to_list', multiple=True)
    for headpiece in headpieces.par:
        if headpiece['bbt'][-1] == lib_id[-1]:
            hp_smiles = headpiece['smiles']
            break
    else:
        print('ERROR::: could not detect headpiede from lib_id')
        sys.exit(1)


    enumerator = Enumerator(wfolder, lib_id, hp_smiles, args.user, bbs=args.bbs_file, lib=lib, base_folder=args.wfolder,
                            verbose=args.verbose)
    enumerator.print_summary_file(enum_reaction, enum_deprotection)
    enumerator.run_graph_enumeration(multireaction, preparations, enum_deprotection, enum_reaction,
                                     n=args.nmols, chunksize=args.nc_per_run, just_json=args.write_json)


if __name__ == '__main__':
    main()
