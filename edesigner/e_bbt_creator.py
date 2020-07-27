# -*- coding: utf-8 -*-
# e_bbt_creator.v.8.1.0
# Jose Alfredo Martin

version = 'e_bbt_creator.v.9.0.0'

# Python modules
import os
import pickle as pic
import time
import argparse
# External modules
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors as desc
import pandas as pd
import numpy as np
# Local modules
from classes.bbt import BBT
from classes.logger import Logger
from classes.parameter_reader import Parameters
from classes.desalter import SmilesDesalter
from classes.mol_property_calculator import MolPropertyCalculator
from classes.parallel import map_parallel

# Functions definitions

def compatible(A, fg):
    """Function compatible has as argument a list with three indexes corresponding to FGs for one BBT
    and returns a bool indicating whether these FGs are compatible in the same BBT
    A : BBT (list of three int corresponding to a BBT)
    fg : instance of the Parameters class (functional group parameters)
    returns : resultado (bool)"""
    resultado = True
    for i in range(3):
        for j in range(i+1, 3):
            if A[j] != 0 and A[i] != 0 and A[j] in fg.par[A[i]]['self_incompatibility']:
                resultado=False
    return resultado

def calculate_fg(row, fgsum, fgsub):
    """This function calculates a new functional group to be added to the
    dataframe
    fgsum : list of str (functional groups to add)
    fgsub : list of str (functional groups to substract)
    returns : count (int)"""
    count = 0
    for fg in fgsum:
        count += row[fg]
    for fg in fgsub:
        count -= row[fg]
    return count
def calculate_fg_parallel(df):
    """This function takes a pandas data frame and applies the calculate_fg function to it
    df : pandas dataframe
    returns : df (pandas dataframe)"""
    for fg in calcfg.par:
        df[fg['name']] = df.apply(calculate_fg, args=(fg['rule_add'], fg['rule_substract']), axis=1)
    return df

def desalt_smiles(row, smilesfield, desalter):
    """This function creates desalted smiles for a pandas dataframe row
    row : row of the dataframe
    smilesfiel : str (name of the smiles field in the row)
    desalter : instance of Smiles_desalter class
    reaction_list : list of dictionaries (rdkit reaction named rxn and reaction name named reaction required)
    returns : result (str, smiles string or np.nan)"""
    smiles = desalter.desalt_smiles(smiles = row[smilesfield])
    if smiles is None:
        return np.nan
    return smiles

def desalt_smiles_parallel(df):
    """This function takes a pandas data frame and applies the desalt_smiles function to it
    df : pandas dataframe
    returns : df (pandas dataframe)"""
    df['SMILES'] = df.apply(desalt_smiles, args=('SMILES', desalter), axis = 1)
    return df

def calculate_natoms(row, smilesfield, calculator):
    """This function calculates a number of heavy atoms based on a smiles stored in a field in a dataframe row
    smilesfield : str (name of the field containing the smiles
    calculator : instance of the Mol_property_calculator class
    return : natoms or np.nan (int or np.nan)"""
    calculator.calculate_properties(smiles=row[smilesfield], props=['py_na'])
    try:
        natoms = calculator.data['py_na']
    except:
        return np.nan
    else:
        return natoms

def calculate_natoms_parallel(df):
    """This function takes a pandas data frame and applies the calculate_natoms function to it
    df : pandas dataframe
    returns : out_df (pandas dataframe)"""
    df['N_ATOMS'] = df.apply(calculate_natoms, args=('SMILES', calculator), axis=1)
    return df

def clasify(row, BBTs, fg):
    """This function takes a row of a pandas data frame and finds which is the BBT it
    belongs to. If it does not find any then returns np.nan
    BBTs : list of instances of BBT class
    fg : instance of the Parameters class (fg parameters)
    returns : int or np.nan"""
    counter = 0
    # The first loop counts how many DEL FGs are there in this molecule
    # and fills a vector with the number of instances for each FG.
    # Only if they are between 1 and 3 it is considered.
    # Then it adds the NO_FG instances up to 3 and
    # checks that the vector is equal to the BBT_long attribute of a BBT
    this_comp_BBT_long = []
    for index, item in enumerate(fg.par):
        if index > 0:
            this_comp_BBT_long.append(row[item['name']])
    n_fgs = sum(this_comp_BBT_long)
    if 0 < n_fgs < 4:
        this_comp_BBT_long = [3 - n_fgs] + this_comp_BBT_long
        for item in BBTs:
            if item.BBT_multi == n_fgs:
                if item.BBT_long == this_comp_BBT_long:
                    return item.index
        else:
            return np.nan
    else:
        return np.nan

def clasify_parallel(df):
    """This function takes a pandas data frame and applies the clasify function to it
    df : pandas dataframe
    returns : df (pandas dataframe)"""
    df['BBT'] = df.apply(clasify, args=(BBTs, fg), axis=1)
    return df

def modify_natoms(row, BBTs, fg):
    """This function takes a row of a pandas data frame and calculates the new number of atoms
    based on the atom difference indicated in itw functional groups
    BBTs : list of instances of BBT class
    fg : instance of the Parameters class (fg parameters)
    returns : n_atoms (int)"""
    n_atoms = row['N_ATOMS']
    for i in BBTs[row['BBT']].BBT:
        n_atoms += fg.par[i]['atom_dif']
    if n_atoms < 1:
        return np.nan
    return n_atoms

def modify_natoms_parallel(df):
    """This function takes a pandas data frame and applies the modify_natoms function to it
    df : pandas dataframe
    returns : df (pandas dataframe)"""
    df['N_ATOMS'] = df.apply(modify_natoms, args = (BBTs, fg), axis = 1)
    return df

def read_compounds(db=None, BBTs=None, calcfg=None, antifg=None):
    """This function reads all databases in the dbpar parameters objec and
    classifies each compound to its BBT. Smiles are desalted and number of
    heavy atoms calculated for each compound. The compounds not classified are
    removed and duplicates consolidated in a single compound.
    db : instance of the Parameters class (database parameters)
    BBTs list of instances of the BBT class
    calcfg : instance of the Parameters class (calcfg parameters)
    antifg : isntance of the Parameters class (antifg parameters)
    returns : comp (pandas dataframe)"""
    com = None
    for db_index, database in enumerate(db.par):
        # read file to dataframe and rename ccolumns
        log.update('Reading ' + database['db'] + ' collection ...')
        separator = {'tab' : '\t', 'space' : ' ', 'comma' : ','}[database['separator']]
        df_chunk = pd.read_csv(os.path.join(database['path'], database['filename']), sep = separator, chunksize = 1000000)
        this_com = None
        for i, df in enumerate(df_chunk):
            log.update('Processing chunk ' + str(i+1) + ' with ' + str(df.shape[0]) + ' records...')
            df.rename(inplace=True, columns={column: column.upper() for column in df.columns})
            df.rename(inplace=True, columns={database['smiles_name']: 'SMILES', database['id_name']: 'ID'})
            df['SOURCE'] = str(db_index + 1) + '_' + database['db']
            df['EXTERNAL'] = database['external_db']
            # filter by properties
            if database['quantity_name'] is not None and database['quantity_filter'] is not None:
                df = df[df[database['quantity_name']] > database['quantity_filter']].copy()
                log.update('    ' + str(len(df)) + ' records remained after filtering by quantity')
            # keep the quatity field for all compounds
            if database['quantity_name'] is not None:
                df.rename(inplace = True, columns = {database['quantity_name'] : 'QUANTITY'})
            else:
                df['QUANTITY'] = np.nan
            if database['rb_filter'] is not None:
                df = df[df[database['rb_name']] < database['rb_filter']].copy()
                log.update('    ' + str(len(df)) + ' records remained after filtering by rb')
            if database['mw_filter'] is not None:
                df = df[df[database['mw_name']] < database['mw_filter']].copy()
                log.update('    ' + str(len(df)) + ' records remained after filtering by mw')
            # calculate calculated fgs
            df = map_parallel(df, calculate_fg_parallel)
            # filter by antifg
            for afg in antifg.par:
                if df.shape[0] > 0:
                    df = df[df[afg['name']] == 0].copy()
            log.update('    ' + str(len(df)) + ' records remained after filtering by antidel FGs')
            # clasify compounds and filter by classification
            df = map_parallel(df, clasify_parallel)
            df.dropna(inplace=True, subset=['BBT'])
            df['BBT'] = df['BBT'].astype('int32')
            log.update('    ' + str(len(df)) + ' records remained after filtering by belonging to a BBT')
            # desalt smiles
            df = map_parallel(df, desalt_smiles_parallel)
            df.dropna(inplace=True, subset=['SMILES'])
            log.update('    ' + str(len(df)) + ' records remained after desalting smiles')
            # count n_atoms
            df = map_parallel(df, calculate_natoms_parallel)
            df.dropna(inplace=True, subset=['N_ATOMS'])
            log.update('    ' + str(len(df)) + ' records remained after removing molecules with incalculable n atoms')
            # reduce columns
            df = df[['ID', 'SMILES', 'SOURCE', 'EXTERNAL', 'BBT', 'N_ATOMS', 'QUANTITY']].copy()
            # modify n_atoms
            df = map_parallel(df, modify_natoms_parallel)
            df.dropna(inplace=True, subset=['N_ATOMS'])
            log.update('    ' + str(len(df)) + ' records remained after removing molecules with incalculable excess n atoms')
            df = df[df['N_ATOMS'] >= dbpar.par[db_index]['na_filter'][0]].copy() #filter BBs where the number of atoms is outside the established range
            df = df[df['N_ATOMS'] < dbpar.par[db_index]['na_filter'][1]].copy()
            log.update('    ' + str(len(df)) + ' records remained after removing compounds with excess natoms')
            # append dataframe to this_master dataframe
            if this_com is None:
                this_com = df.copy()
            else:
                this_com = this_com.append(df)  
            this_com.sort_values(inplace=True, by='QUANTITY', ascending=False, na_position='last')
            log.update('    ' + str(len(this_com)) + ' records remained after pooling chunks')
        # append dataframe to master dataframe
        if com is None:
            com = this_com.copy()
        else:
            com = com.append(this_com)
    del df
    # eliminate duplicates
    log.update('A total of ' + str(len(com)) + ' records before combination and elimination of duplicates')
    com.sort_values(by=['SOURCE', 'SMILES'], inplace=True)
    id_list_df = com.groupby('SMILES')['ID'].apply(lambda X: ';'.join([str(item) for item in X.tolist()])).reset_index()
    id_list_df.rename(inplace=True, columns={'ID' : 'IDS'})
    com = com.groupby('SMILES').nth(0).reset_index()
    com = pd.merge(com, id_list_df, on='SMILES', how='left')
    com.sort_values(by=['N_ATOMS'], inplace=True)
    log.update('A total of ' + str(len(com)) + ' records remained after combination and elimination of duplicates')
    return com

def report_compound_files(comp_path, comp, BBTs):
    """report_compound_files create a smi file for each BBT and put in the file all compounds in compound
    assigned to that BBT_index. Warning: this function does not clea files used in previous script runs
    comp_path : str (folder where compounds files will be saved)
    comp : pandas dataframe containing all the compounds
    returns : None"""
    for BBT in BBTs:
        df = comp[comp['BBT'] == BBT.index]
        if df.shape[0] > 0:
            df['RID'] = df.apply(lambda row: ':'.join([str(row['N_ATOMS']), str(row['SOURCE']), str(row['ID'])]), axis = 1)
            dfi = df[df['EXTERNAL'] == False]
            df = df[['SMILES', 'RID']]
            df.to_csv(os.path.join(comp_path, str(BBT.index) + '.smi'), sep = ' ', index = False, header = False)
            if dfi.shape[0] > 0:
                dfi = dfi[['SMILES', 'RID']]
                dfi.to_csv(os.path.join(comp_path, str(BBT.index)+ '.int.smi'), sep = ' ', index = False, header = False)
    return None


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = """e_bbt_creator:
    This script searches the files of annotated compounds and classifies all compounds
    into the building block type (BBT) they belong to, or discard them if they do not belong
    to a specific BBT.
    The script starts generating by comprehension all building block types as a combination of
    up to three compatible functional groups (FGs).
    Once the list of BBTs are set the compounds from different collections are classified and
    duplicates, compounds which smiles string cannot be parsed with rdkit or compounds not belonging to
    any BBT are eliminated.
    The effective number of atoms (number of atoms minus the number of atoms that will be lost upon
    reaction) of each building block is calculatede, and the building blocks corresponding to
    each BBT stored in a file.
    Two files per BBT are created, one containing all building blocks corresponding to that BBT and the
    other containing building blocks in internal collections only.
    The characteristics of each BBT are stored in an object. This includes number of BBs containing each
    possible number of atoms, functional groups corresponding to this BBT etc, and a list of all BBTs
    is stored in a file for further use. The name of the list is BBTs and is used throughout all scripts
    in eDESIGNER.
    """)
    parser.add_argument('-wf', '--wfolder', help='Working Folder', type=str, default='./')
    args = parser.parse_args()

    # Initialization
    error_found = False
    tic = time.time()
    path = Parameters(os.path.join(args.wfolder, 'PAR/path.par'), fsource = 'dict', how = 'to_dict', multiple = False)
    log = Logger(os.path.join(args.wfolder, 'logs', 'e_bbt_creator.log'), prepend_timestamp=True)
    dbpar = Parameters(os.path.join(args.wfolder, 'resources', 'db.par'), fsource='list', how='to_list', multiple=True)
    par = Parameters(os.path.join(args.wfolder, 'resources', 'par.par'), fsource='dict', how='to_dict', multiple=False)
    fg = Parameters(os.path.join(args.wfolder, 'resources', 'fg.par'), fsource='list', how='to_list', multiple=True)
    calcfg = Parameters(os.path.join(args.wfolder, 'resources', 'calcfg.par'), fsource='list', how='to_list', multiple=True)
    antifg = Parameters(os.path.join(args.wfolder, 'resources', 'antifg.par'), fsource='list', how='to_list', multiple=True)
    headpieces = Parameters(os.path.join(args.wfolder, 'resources', 'headpieces.par'), fsource='list', how='to_list', multiple=True)
    log.update(version)
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'db.par'), 'db parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'par.par'), 'par parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'fg.par'), 'fg parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'calcfg.par'), 'calcfg parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'antifg.par'), 'antifg parameters')
    log.insert_file_in_log(os.path.join(args.wfolder, 'resources', 'headpieces.par'), 'headpieces parameters')
    desalter = SmilesDesalter()
    calculator = MolPropertyCalculator()
    for item in [dbpar, par, fg, calcfg, antifg, headpieces]:
        if item.success is None:
            error_found = True
            log.update('    Error(s) found while reading parameters: ' + item.path)
            for error in item.errors:
                log.update('        ' + error)

    # Body of the script
    if not error_found:
        log.update('Creating db run folder...')
        comps_path = os.path.join(args.wfolder, 'comps', log.strtimestamp)
        command = f'[ ! -d {comps_path} ] && mkdir {comps_path}'
        os.system(command)
        files_list = [f for f in os.listdir(comps_path)]
        for file_item in files_list:
            os.remove(os.path.join(comps_path, file_item))

        # Generates the list of all BBTs by comprehension ignoring incompatible BBTs but including [0, 0, 0]
        log.update('Generating BB types...')
        BBT_list = [[i, j, k] for i in range(len(fg.par)) for j in range(i, len(fg.par)) for k in range(j, len(fg.par)) if compatible([i, j, k], fg)]
        BBTs = [BBT(BBT=BBT_list[i], fg=fg, headpieces=headpieces, index=i) for i in range(len(BBT_list))]

        # read compound sets and get valid compounds within valid BBTs
        log.update('Reading compound sets...')
        comp = read_compounds(db=dbpar, BBTs=BBTs, calcfg=calcfg, antifg=antifg)

        # update BBTs
        log.update('Updating BBTs...')
        for i in range(len(BBTs)):
            BBTs[i].update(comp)
        BBTs.sort(key = lambda X: X.n_compounds[-1], reverse = True)
        BBTs.sort(key = lambda X: X.BBT_multi)
        for i in range(len(BBTs)):
            BBTs[i].order = i
        BBTs.sort(key = lambda X: X.index)

        # reporting files
        log.update('Reporting compound files...')
        report_compound_files(comps_path, comp, BBTs)
        log.update('pickling BBTs object...')
        with open(os.path.join(args.wfolder, 'data', log.strtimestamp + '_BBTs.pic'), 'wb') as f:
            pic.dump(BBTs, f, -1)
        log.update('Reporting BBTs file...')
        comp['FGs'] = comp.apply(lambda row: ';'.join([BBTs[row['BBT']].BBT_name[i] for i in range(3)]), axis=1)
        comp.to_csv(os.path.join(args.wfolder, 'data', log.strtimestamp + '_BBT_report.csv'), index=False)
        # time and end the program
        tac = time.time()
        log.update(f'Run time = {round((tac - tic) / 60.0, 1)} min')
        log.update('OK')
