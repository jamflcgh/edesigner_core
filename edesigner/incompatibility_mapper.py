# -*- coding: utf-8 -*-
# incompatibility_mapper v.12.00
# Jose Alfredo Martin, 2023

# imports

__version__ = 'incompatibility_mapper.v.12.0.0'
__author__ = 'Alfredo Martin 2023'


# python modules
import os
import copy
import time
import argparse

# Local modules
from classes.logger import Logger
from classes.parameter_reader import Parameters

def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description = """This script creates several reports to inspect visually the 
    incompatibilities of reactions and de-protections with functional groups. It runs the report based on the 
    eDESIGNER parameters files. If a folder containing parameters is passed then the report is written in that 
    folder otherwise the parameters in the repo are used and the report is written in the folder "data" in the 
    repo""")
    # Defines whether library is coming from eDESIGNER or is external
    parser.add_argument('-pF', '--par_folder',
                        help="""folder containing parameters. If not passed the parameters in the repo are used""",
                        type=str,
                        default=None)
    args = parser.parse_args()
    if args.par_folder is not None:
        assert os.path.isdir(args.par_folder), f'{args.par_folder} does not exist'
    return args

def main():
    args = parse_args()
    if args.par_folder is not None:
        out_folder = args.par_folder
        RESOURCESFOLDER = args.par_folder
    else:
        out_folder = os.path.join(os.environ['EDESIGNER_FOLDER'], 'data')
        RESOURCESFOLDER = os.environ['EDESIGNER_PARFOLDER']
    fg = Parameters(os.path.join(RESOURCESFOLDER, 'fg.par'), fsource='list', how='to_list', multiple=True)
    reaction = Parameters(os.path.join(RESOURCESFOLDER, 'reaction.par'), fsource='list', how='to_list', multiple=True)
    deprotection = Parameters(os.path.join(RESOURCESFOLDER, 'deprotection.par'), fsource='list', how='to_list', multiple=True)
    reaction_report = os.path.join(out_folder, 'reaction_incompatibility_report.txt')
    deprotection_report = os.path.join(out_folder, 'deprotection_incompatibility_report.txt')
    reaction_on_graph = os.path.join(out_folder, 'reaction_incompatibility_on_graph.txt')
    reaction_off_graph = os.path.join(out_folder, 'reaction_incompatibility_off_graph.txt')
    deprotection_graph = os.path.join(out_folder, 'deprotection_incompatibility_graph.txt')
    fg_report = os.path.join(out_folder, 'fg_incompatibility_report.txt')
    fg_graph = os.path.join(out_folder, 'fg_incompatibility_graph.txt')
    # fg_report
    with open(fg_report, 'w') as f:
        f.write('FG INCOMPATIBILIT REPORT\n')
        f.write('\n')
        for i, FG1 in enumerate(fg.par):
            if i == 0:
                continue
            f.write(f"NAME:{FG1['name']} ({i})\n")
            for j, FG2 in enumerate(fg.par):
                if j == 0:
                    continue
                if j in FG1['self_incompatibility']:
                    f.write(f"  {FG2['name']} ({j})\n")
            f.write('\n')
    print(f'writen fg_incompatibility_report: {fg_report}')
    # reaction_report
    with open(reaction_report, 'w') as f:
        f.write('REACTION INCOMPATIBILIT REPORT\n')
        f.write('\n')
        for i, r in enumerate(reaction.par):
            if i == 0:
                continue
            f.write(f"NAME: {r['name']}\n")
            if r['fg_input_on_off'][0] != 0:
                f.write(f"INPUT FG: {fg.par[r['fg_input_on_off'][0]]['name']} ({r['fg_input_on_off'][0]})\n")
            if r['fg_input_on_off'][1] != 0:
                f.write(f"INPUT FG: {fg.par[r['fg_input_on_off'][1]]['name']} ({r['fg_input_on_off'][1]})\n")
            if r['fg_output_on_off'][0] != 0:
                f.write(f"OUTPUT FG: {fg.par[r['fg_output_on_off'][0]]['name']} ({r['fg_output_on_off'][0]})\n")
            if r['fg_output_on_off'][1] != 0:
                f.write(f"OUTPUT FG: {fg.par[r['fg_output_on_off'][1]]['name']} ({r['fg_output_on_off'][1]})\n")
            f.write(f"INCOMPATIBLE FGS ON:\n")
            for item in r['excluded_on']:
                f.write(f"  {fg.par[item]['name']} ({item})\n")
            f.write(f"INCOMPATIBLE FGS OFF:\n")
            for item in r['excluded_off']:
                f.write(f"  {fg.par[item]['name']} ({item})\n")
            f.write(f"COMPATIBLE FGS ON:\n")
            for j, item in enumerate(fg.par):
                if j != 0 and j not in r['excluded_on']:
                    f.write(f"  {item['name']} ({j})\n")
            f.write(f"COMPATIBLE FGS OFF:\n")
            for j, item in enumerate(fg.par):
                if j != 0 and j not in r['excluded_off']:
                    f.write(f"  {item['name']} ({j})\n")
            f.write('\n')
    print(f'writen reaction_incompatibility_report: {reaction_report}')
    # deprotection report
    with open(deprotection_report, 'w') as f:
        f.write('DEPROTECTION INCOMPATIBILIT REPORT\n')
        f.write('\n')
        for i, r in enumerate(deprotection.par):
            if i == 0:
                continue
            f.write(f"NAME: {r['name']}\n")
            if r['fg_input_on_off'][0] != 0:
                f.write(f"INPUT FG: {fg.par[r['fg_input_on_off'][0]]['name']} ({r['fg_input_on_off'][0]})\n")
            if r['fg_input_on_off'][1] != 0:
                f.write(f"INPUT FG: {fg.par[r['fg_input_on_off'][1]]['name']} ({r['fg_input_on_off'][1]})\n")
            if r['fg_output_on_off'][0] != 0:
                f.write(f"OUTPUT FG: {fg.par[r['fg_output_on_off'][0]]['name']} ({r['fg_output_on_off'][0]})\n")
            if r['fg_output_on_off'][1] != 0:
                f.write(f"OUTPUT FG: {fg.par[r['fg_output_on_off'][1]]['name']} ({r['fg_output_on_off'][1]})\n")
            f.write(f"INCOMPATIBLE FGS ON:\n")
            for item in r['excluded_on']:
                f.write(f"  {fg.par[item]['name']} ({item})\n")
            f.write(f"COMPATIBLE FGS ON:\n")
            for j, item in enumerate(fg.par):
                if j != 0 and j not in r['excluded_on']:
                    f.write(f"  {item['name']} ({j})\n")
            f.write('\n')
    print(f'writen deprotection_incompatibility_report: {deprotection_report}')
    # reaction_graph_on
    with open(reaction_on_graph, 'w') as f:
        linea = 'reaction\tFG\tincompatible\n'
        f.write(linea)
        for FG in fg.par:
            for R in reaction.par:
                linea = R['name'] + '\t' + FG['name'] + '\t'
                if FG['index'] in R['excluded_on']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)
    print(f'writen reaction_incompatibility_on_graph: {reaction_on_graph}')
    # reaction_graph_off
    with open(reaction_off_graph, 'w') as f:
        linea = 'reaction\tFG\tincompatible\n'
        f.write(linea)
        for FG in fg.par:
            for R in reaction.par:
                linea = R['name'] + '\t' + FG['name'] + '\t'
                if FG['index'] in R['excluded_off']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)
    print(f'writen reaction_incompatibility_off_graph: {reaction_off_graph}')
    # deprotection_graph
    with open(deprotection_graph, 'w') as f:
        linea = 'deprotection\tFG\tincompatible\n'
        f.write(linea)
        for FG in fg.par:
            for R in deprotection.par:
                linea = R['name'] + '\t' + FG['name'] + '\t'
                if FG['index'] in R['excluded_on']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)
    print(f'writen deprotection_incompatibility_graph: {deprotection_graph}')
    # fg_graph
    with open(fg_graph, 'w') as f:
        linea = 'FG1\tFG2\tincompatible\n'
        f.write(linea)
        for FG1 in fg.par:
            for FG2 in fg.par:
                linea = FG1['name'] + '\t' + FG2['name'] + '\t'
                if FG1['index'] in FG2['self_incompatibility']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)
    print(f'writen fg_incompatibility_graph: {fg_graph}')





if __name__ == '__main__':
    print(__version__)
    print(__author__)
    main()

    
        
    
