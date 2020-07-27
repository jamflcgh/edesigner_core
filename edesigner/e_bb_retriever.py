# -*- coding: utf-8 -*-
# e_bb_retriever
# Jose Alfredo Martin

version = 'e_bb_retriever.v.9.0.0'

# Python modules
import os
import pickle as pic
import argparse
# External modules
import pandas as pd
# Local modules
from classes.libdesign import LibDesign
from classes.logger import Logger


if __name__ == '__main__':
    # Arg parser
    parser = argparse.ArgumentParser(description="""e_bb_retreaver collects all the building blocks for a giver
    libDESIGN into csv files in the folder that is given through arguments. The input arguments are wfolder (the
    folder system of eDESIGNER, tk: token, the combination of the Db_Run time stamp and run number separated by _,
    dn: the design number or a list of design numbers and of: the output folder. There are two files created for each 
    cycle and named ddn_Cn.csv and dn_Cn_all.csv where n is the cycle number and ddn is the design number. 
    These files contain the building block smiles and ids that should be used to obtain the maximum size of the 
    library while maintaining the heavy atom distribution. The output folder will be created if it does not already
    exist""")
    parser.add_argument('-wf', '--wfolder', help='Working Folder', type=str, default='./')
    parser.add_argument('-tk', '--token', help='combination of Db_Run time stamp and run number separated by _',
                        type=str, default=None)
    parser.add_argument('-df', '--design_number',
                        help='number of the design of interest or a list of designs , separated', type=str,
                        default=None)
    parser.add_argument('-of', '--output_folder', help='folder where the files will be saved', type=str,  default=None)
    args = parser.parse_args()
    assert args.token is not None, 'A token must be provided'
    assert args.design_number is not None, 'A design number or a list of design numbers must be provided'
    assert args.output_folder is not None, 'An output folder must be provided'
    args.design_number = args.design_number.split(',')
    args.design_number = [item.strip(' ') for item in args.design_number]
    args.design_number = [int(item) for item in args.design_number]

    # Initialization
    log = Logger(os.path.join(args.wfolder, 'logs', args.token + '_lib_design_interpreter.log'))
    log.update(version)
    try:  # Try to create the output folder but do not delete it if it aready exists
        os.mkdir(args.output_folder)
    except:
        pass

    # Body of the script
    log.update('Loading pickled designs to object...')
    with open(os.path.join(args.wfolder, 'results', args.token + '_libDESIGNS.pic'), 'rb') as f:
        lib_list = []
        while True:
            try:
                lib_list.append(pic.load(f))
            except:  # reached the end of the file
                break
    log.update('Creating dataframes and dumping them into files...')
    for lib_id in args.design_number:
        design = lib_list[lib_id]
        if design.id != lib_id:
            log.update(f'Design {lib_id} in list does not match with design.id {design.id}')
        else:
            for cycle in range(design.n_cycles):
                all_dfs = []
                int_dfs = []
                for bbt, limit in zip(design.bbts[cycle], design.int_limits[cycle]):
                    if limit > 0:
                        # It might be that a bbt does not have internal compounds but it does have all
                        # compounds and therefore the file for internal compounds does not exists.
                        # That is why we need to check for the number of compounds to avoid an I/O exception
                        df_int = pd.read_csv(os.path.join(args.wfolder, 'comps', args.token.split('_')[0], str(bbt) + '.int.smi'),
                                             sep=' ', header=None, names=['smiles', 'id'], nrows=limit)
                        int_dfs.append(df_int.copy())
                for bbt, limit in zip(design.bbts[cycle], design.all_limits[cycle]):
                    if limit > 0:  # This should not be necessary, but just in case
                        df_all = pd.read_csv(os.path.join(args.wfolder, 'comps', args.token.split('_')[0], str(bbt) + '.smi'),
                                             sep=' ', header=None, names=['smiles', 'id'], nrows=limit)
                        all_dfs.append(df_all.copy())
                all_dfs = pd.concat(all_dfs)
                int_dfs = pd.concat(int_dfs)
                all_dfs.to_csv(os.path.join(args.output_folder, args.token + '_' + str(lib_id) + '_C' + str(cycle + 1) + '_all.csv'), index=False)
                int_dfs.to_csv( os.path.join(args.output_folder, args.token + '_' + str(lib_id) + '_C' + str(cycle + 1) + '_int.csv'), index=False)
    log.update('OK')
