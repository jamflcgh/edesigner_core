# -*- coding: utf-8 -*-
# stats_analyzer
# Jose Alfredo Martin 2020

__version__ = 'stats_analyzer.v.12.0.0'
__author__ = 'Alfredo Martin 2023'

import os
from tqdm import tqdm
import _pickle as pic
from classes.parameter_reader import Parameters
import argparse
import pandas as pd


def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description="""writes statistics of usage for bbts, reactons and deprotections""")
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
    parser.add_argument('-v', '--verbose',
                        help="""verbose output""",
                        action='store_true')
    args = parser.parse_args()
    assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
    args.wfolder = os.path.abspath(args.wfolder)
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id)), f'{os.path.join(args.wfolder, args.run_id)} does not exist'
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id, args.ed_run_id)), f'{os.path.join(args.wfolder, args.run_id, args.ed_nun_id)} does not exist'
    return args

def main():
    args = parse_args()
    PARFOLDER = os.path.join(args.wfolder, args.run_id, 'resources')
    libdesigns_file = os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'libDESIGNs.pic')
    enum_deprotection = Parameters(os.path.join(PARFOLDER, 'enum_deprotection.par'),
                                   fsource='list', how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(PARFOLDER, 'enum_reaction.par'),
                               fsource='list', how='to_list',multiple=True)
    bbts_file = os.path.join(args.wfolder, args.run_id, 'results', 'BBTs.pic')
    with open(bbts_file, 'rb') as f:
        BBTs = pic.load(f)
    designs = []
    with open(libdesigns_file, 'rb') as f:
        while True:
            try:
                design = pic.load(f)
                designs.append(design)
            except:
                break
    reaction_usage=[{'name': rxn['enum_name'], 'usage': 0} for rxn in enum_reaction.par]
    deprotection_usage=[{'name': rxn['enum_name'], 'usage': 0} for rxn in enum_reaction.par]
    bbts_usage=[{'name': ':'.join(bbt.BBT_name), 'usage': 0} for bbt in BBTs]
    for design in tqdm(designs):
        for idx in design.reactions:
            reaction_usage[idx]['usage'] += 1
        for idx in design.deprotections:
            if idx != 0:
                deprotection_usage[idx]['usage'] += 1
        for item in design.bbts:
            for idx in item:
                bbts_usage[idx]['usage'] += 1
    df = pd.DataFrame(reaction_usage)
    df.to_csv(os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'rxn_usage_report.csv'), index=False)
    print(f"rxn usage report: {os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'rxn_usage_report.csv')}")
    df = pd.DataFrame(deprotection_usage)
    df.to_csv(os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'depr_usage_report.csv'), index=False)
    print(f"depr usage report: {os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'depr_usage_report.csv')}")
    df = pd.DataFrame(bbts_usage)
    df.to_csv(os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'bbt_usage_report.csv'), index=False)
    print(f"bbt usage report: {os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'bbt_usage_report.csv')}")


if __name__ == '__main__':
    print(__version__)
    print(__author__)
    main()
