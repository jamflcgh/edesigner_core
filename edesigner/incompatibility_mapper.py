# -*- coding: utf-8 -*-
# incompatibility_mapper
# Jose Alfredo Martin, 2019

# imports

version = 'incompatibility_mapper.v.1.0.0'
print(version)
print('type "incompatibility_mapper.py -h" for help')
print('Importing modules...')

# python modules
import os
import copy
import time
import argparse

# Local modules
from classes.logger import Logger
from classes.parameter_reader import Parameters


if __name__ == '__main__':
    # Arg parser
    parser = argparse.ArgumentParser(description="""incompatibility_mapper: 
    This script creates files to visualize incompatibilities in reactions and functional groups.
    """)
    parser.add_argument('-wf', '--wfolder', help='Working Folder', type=str, default='./')
    args = parser.parse_args()
    # Initiallization
    tic = time.time()
    log = Logger(os.path.join(args.wfolder, 'logs', 'incompatibility_mapper.log'))
    log.update(version)
    fg = Parameters(os.path.join(args.wfolder, 'resources/fg.par'), fsource = 'list', how = 'to_list', multiple = True)
    reaction = Parameters(os.path.join(args.wfolder, 'resources/reaction.par'), fsource = 'list', how = 'to_list', multiple = True)
    deprotection = Parameters(os.path.join(args.wfolder, 'resources/deprotection.par'), fsource = 'list', how = 'to_list', multiple = True)

    log.update('creating FG incompatibility fileincompatibility...')
    with open(os.path.join(args.wfolder, 'data/fg_incompatibility.txt'), 'w') as f:
        linea = 'FG1\tFG2\tincompatible\n'
        f.write(linea)
        for FG1 in fg.par:
            for FG2 in fg.par:
                linea = FG1['name'] + '\t' + FG2['name'] + '\t'
                if FG2['index'] in FG1['self_incompatibility']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)


    log.update('creating reaction incompatibility fileincompatibility...')
    with open(os.path.join(args.wfolder, 'data/reaction_incompatibility.txt'), 'w') as f:
        linea = 'Reaction\tFG\tonDNAincompatible\toffDNAincompatible\n'
        f.write(linea)
        for item in reaction.par:
            for FG in fg.par:
                linea = item['name'] + '\t' + FG['name'] + '\t'
                if FG['index'] in item['excluded_on']:
                    linea += '1\t'
                else:
                    linea += '0\t'
                if FG['index'] in item['excluded_off']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)

    log.update('creating deprotection incompatibility fileincompatibility...')
    with open(os.path.join(args.wfolder, 'data/deprotection_incompatibility.txt'), 'w') as f:
        linea = 'Deprotection\tFG\tonDNAincompatible\n'
        f.write(linea)
        for item in deprotection.par:
            for FG in fg.par:
                linea = item['name'] + '\t' + FG['name'] + '\t'
                if FG['index'] in item['excluded_on']:
                    linea += '1\n'
                else:
                    linea += '0\n'
                f.write(linea)

    # time it and end script
    tac = time.time()
    log.update(f'Running time: {round((tac-tic)/60, 1)} min.')
    log.update('OK')
    
        
    
