# -*- coding: utf-8 -*-
# lib_bb_analyzer
# Jose Alfredo Martin 2020

version = 'lib_bb_analyzer.v.9.0.0'

# Python modules
import os
import pickle as pic
import copy
import time
import argparse
# Local modules
from classes.bbt import BBT
from classes.design import Design
from classes.libdesign import LibDesign
from classes.logger import Logger
from classes.parameter_reader import Parameters
# other modules
import numpy as np
import pandas as pd


def process_bbts(BBTs):
    """extracts key information from BBTs in the form of dictionaries of np arrays and scalars and puts them in
    a dictionary using the BBT indexes as keys
    BBTs: list of BBT class instance
    returns pBBTs: dict of dicts"""
    pBBTs = {}
    for bbt in BBTs:
        if sum(bbt.n_compounds) > 0:
            this_dict = {}
            this_dict['index'] = bbt.index
            this_dict['multi'] = bbt.BBT_multi
            this_dict['order'] = bbt.order
            this_dict['min_atoms'] = bbt.min_atoms
            this_dict['name'] = ';'.join(bbt.BBT_name)
            this_dict['fg_1'] = bbt.BBT_name[0]
            this_dict['fg_2'] = bbt.BBT_name[1]
            this_dict['fg_3'] = bbt.BBT_name[2]
            this_dict['n_all'] = max(bbt.n_compounds)
            this_dict['n_int'] = max(bbt.n_internal)
            this_dict['n_ext'] = max(bbt.n_external)
            n_compounds = np.array(bbt.n_compounds)
            n_compounds_s = [0] + bbt.n_compounds
            n_compounds_s = n_compounds_s[0:-1]
            n_compounds_s = np.array(n_compounds_s)
            this_dict['hist_all'] = n_compounds - n_compounds_s
            n_compounds = np.array(bbt.n_internal)
            n_compounds_s = [0] + bbt.n_internal
            n_compounds_s = n_compounds_s[0:-1]
            n_compounds_s = np.array(n_compounds_s)
            this_dict['hist_int'] = n_compounds - n_compounds_s
            n_compounds = np.array(bbt.n_external)
            n_compounds_s = [0] + bbt.n_external
            n_compounds_s = n_compounds_s[0:-1]
            n_compounds_s = np.array(n_compounds_s)
            this_dict['hist_ext'] = n_compounds - n_compounds_s
            this_dict['max_atoms'] = 0
            for i in range(100):
                if this_dict['hist_all'][i] > 0:
                    this_dict['max_atoms'] = i +1
            pBBTs[bbt.index] = this_dict
    return pBBTs

def analyze_reactions(lib_list, enum_reaction, enum_deprotection):
    """analyze_reactions takes as arguments a list of libDESIGNS and parameters objects defining linking reactions and
    deprotection reactions and couts how many times each reaction has been run overall in all designs and returns
    the results in two pandas dataframes
    reactionsliblist: list of libDESIGNS
    enum_reactions: instance of Parameters class
    enum_deprotections: instance of Parameters class
    retuns: rdf, ddf: tuple of pandas dataframes"""
    all_reactions = []
    all_deprotections = []
    for library in lib_list:
        all_reactions += library.reactions
        all_deprotections += library.deprotections
    rdf = pd.DataFrame({'rindex': all_reactions})
    rdf = rdf[rdf['rindex'] > 0]
    rdf['count'] = 1
    rdf = rdf.groupby('rindex')['count'].sum().reset_index()
    rdf['reaction name'] = rdf.apply(lambda row: enum_reaction.par[row['rindex']]['enum_name'], axis=1)
    rdf.sort_values(inplace=True, by=['count'], ascending=False)
    ddf = pd.DataFrame({'rindex': all_deprotections})
    ddf = ddf[ddf['rindex'] > 0]
    ddf['count'] = 1
    ddf = ddf.groupby('rindex')['count'].sum().reset_index()
    ddf['reaction name'] = ddf.apply(lambda row: enum_deprotection.par[row['rindex']]['enum_name'], axis=1)
    ddf.sort_values(inplace=True, by=['count'], ascending=False)
    return rdf, ddf

def analyze_bb_usage(lib_list, pBBTs):
    """analyze_bb_usage reads the takes as arguments a list of libDESIGNS and a processed list of BBT objects and
    extracts information on how many building blocks have been used. The analysis is run in different dimensions such
    as the source, the number of atoms or the multiplicity of the bb and results are returned as a pandas dataframe
    liblist: list of libDESIGNS
    pBBTs: dictionary of dictionaries
    returs: df: pandas dataframe"""
    bbs_array = np.array([[0 for _ in range(12)] for _ in range(100)])
    ubbs_array = np.array([[0 for _ in range(12)] for _ in range(100)])
    usedbbts = []
    usedatoms = []
    maxnatom = 0
    for library in lib_list:
        for cycle in range(library.n_cycles):
            for bbt in library.bbts[cycle]:
                for natom in range(library.best_index[cycle] + 1):
                    bbs_array[natom][0] += pBBTs[bbt]['hist_all'][natom]
                    bbs_array[natom][1] += pBBTs[bbt]['hist_int'][natom]
                    bbs_array[natom][2] += pBBTs[bbt]['hist_ext'][natom]
                    if pBBTs[bbt]['multi'] == 1:
                        bbs_array[natom][3] += pBBTs[bbt]['hist_all'][natom]
                        bbs_array[natom][4] += pBBTs[bbt]['hist_int'][natom]
                        bbs_array[natom][5] += pBBTs[bbt]['hist_ext'][natom]
                    elif pBBTs[bbt]['multi'] == 2:
                        bbs_array[natom][6] += pBBTs[bbt]['hist_all'][natom]
                        bbs_array[natom][7] += pBBTs[bbt]['hist_int'][natom]
                        bbs_array[natom][8] += pBBTs[bbt]['hist_ext'][natom]
                    else:
                        bbs_array[natom][9] += pBBTs[bbt]['hist_all'][natom]
                        bbs_array[natom][10] += pBBTs[bbt]['hist_int'][natom]
                        bbs_array[natom][11] += pBBTs[bbt]['hist_ext'][natom]
                    maxnatom = max(maxnatom, natom)
                if not bbt in usedbbts:
                    for natom in range(library.best_index[cycle] + 1):
                        ubbs_array[natom][0] += pBBTs[bbt]['hist_all'][natom]
                        ubbs_array[natom][1] += pBBTs[bbt]['hist_int'][natom]
                        ubbs_array[natom][2] += pBBTs[bbt]['hist_ext'][natom]
                        if pBBTs[bbt]['multi'] == 1:
                            ubbs_array[natom][3] += pBBTs[bbt]['hist_all'][natom]
                            ubbs_array[natom][4] += pBBTs[bbt]['hist_int'][natom]
                            ubbs_array[natom][5] += pBBTs[bbt]['hist_ext'][natom]
                        elif pBBTs[bbt]['multi'] == 2:
                            ubbs_array[natom][6] += pBBTs[bbt]['hist_all'][natom]
                            ubbs_array[natom][7] += pBBTs[bbt]['hist_int'][natom]
                            ubbs_array[natom][8] += pBBTs[bbt]['hist_ext'][natom]
                        else:
                            ubbs_array[natom][9] += pBBTs[bbt]['hist_all'][natom]
                            ubbs_array[natom][10] += pBBTs[bbt]['hist_int'][natom]
                            ubbs_array[natom][11] += pBBTs[bbt]['hist_ext'][natom]
                else:
                    usedatomlist = [usedatoms[i] for i in range(len(usedatoms)) if usedbbts[i] == bbt]
                    maxusedatoms = max(usedatomlist)
                    for natom in range(maxusedatoms + 1,  library.best_index[cycle] + 1):
                        ubbs_array[natom][0] += pBBTs[bbt]['hist_all'][natom]
                        ubbs_array[natom][1] += pBBTs[bbt]['hist_int'][natom]
                        ubbs_array[natom][2] += pBBTs[bbt]['hist_ext'][natom]
                        if pBBTs[bbt]['multi'] == 1:
                            ubbs_array[natom][3] += pBBTs[bbt]['hist_all'][natom]
                            ubbs_array[natom][4] += pBBTs[bbt]['hist_int'][natom]
                            ubbs_array[natom][5] += pBBTs[bbt]['hist_ext'][natom]
                        elif pBBTs[bbt]['multi'] == 2:
                            ubbs_array[natom][6] += pBBTs[bbt]['hist_all'][natom]
                            ubbs_array[natom][7] += pBBTs[bbt]['hist_int'][natom]
                            ubbs_array[natom][8] += pBBTs[bbt]['hist_ext'][natom]
                        else:
                            ubbs_array[natom][9] += pBBTs[bbt]['hist_all'][natom]
                            ubbs_array[natom][10] += pBBTs[bbt]['hist_int'][natom]
                            ubbs_array[natom][11] += pBBTs[bbt]['hist_ext'][natom]
                usedbbts.append(bbt)
                usedatoms.append(library.best_index[cycle])
    for i in range(bbs_array.shape[0]):
        if bbs_array[i][0] > 0:
            minnatom = i
            break
    bbs_array = bbs_array[minnatom: maxnatom + 1, :]
    ubbs_array = ubbs_array[minnatom: maxnatom + 1, :]
    times_array = bbs_array / ubbs_array
    natoms_array = np.array(list(range(minnatom, maxnatom +1)))
    df = pd.DataFrame({'number atoms': 12 * list(natoms_array)})
    n_instances = 3 * natoms_array.shape[0]
    df['multiplicity'] = ['all'] * n_instances + ['1'] * n_instances + ['2'] * n_instances + ['3'] * n_instances
    df['source'] = (['ALL'] * natoms_array.shape[0] + ['INTERNAL'] * natoms_array.shape[0] + ['EXTERNAL'] * natoms_array.shape[0]) * 4
    print(df.shape[0])
    df['number of BBs'] = list(bbs_array.T.flatten())
    df['number of unique BBs'] = list(ubbs_array.T.flatten())
    df['unique BBs usage'] = list(times_array.T.flatten())
    return df


if __name__ == '__main__':

    # Args parser
    parser = argparse.ArgumentParser(description="""This script runs an analysis on the reactions and building blocks
    used in the selected libraries by libDESIGNER. The wfolder argument indicates where to look for the souce file
    (a file containing a pickled libDESIGN instance in each record) the file must be in a folder named results within 
    wfolder. The timestamp and tag arguments define the name of the soruce file and the name of the output files.
    Output files are written in a folder named data within wfolder""")
    parser.add_argument('-wf', '--wfolder', help='Working Folder (default where this python scritp is run)',
                        type=str, default='./')
    parser.add_argument('-ts', '--timestamp', help='Timestamp (date when bbt creator was run with format YYYYMMDD)',
                        type=str, default=None)
    parser.add_argument('-tg', '--tag', help='identificator of the experiment', type=str, default=None)
    args = parser.parse_args()

    # Initialization
    tic = time.time()
    token = '_'.join([args.timestamp, args.tag])
    log = Logger(os.path.join(args.wfolder, 'logs', token + '_lib_bb_analyzer.log'))
    log.update(version)
    reaction = Parameters(os.path.join(args.wfolder, 'resources', 'reaction.par'), fsource='list', how='to_list',
                          multiple=True)
    deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'deprotection.par'), fsource='list',
                              how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(args.wfolder, 'resources', 'enum_reaction.par'), fsource='list',
                               how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(args.wfolder, 'resources', 'enum_deprotection.par'), fsource='list',
                                   how='to_list', multiple=True)

    # Body of the script
    log.update('Reading bbts..')
    with open(os.path.join(args.wfolder, 'data', args.timestamp + '_BBTs.pic'), 'rb') as f:
        BBTs = pic.load(f)
    pBBTs =  process_bbts(BBTs)   
    
    log.update('Reading lib designs...')
    with open(os.path.join(args.wfolder, 'results', token + '_libDESIGNS.pic'), 'rb') as f:
        lib_list = []
        while True:
            try:
                lib_list.append(pic.load(f))
            except:  # reached the end of the file
                break

    log.update('analyzing lib designs...')
    pBBTs = process_bbts(BBTs)
    rdf, ddf = analyze_reactions(lib_list, enum_reaction, enum_deprotection)
    rdf.to_csv(os.path.join(args.wfolder, 'data', token + '_reaction_analysis.csv'), index=False)
    ddf.to_csv(os.path.join(args.wfolder, 'data', token + '_deprotection_analysis.csv'), index=False)
    log.update('analyzing bb usage...')
    df = analyze_bb_usage(lib_list, pBBTs)
    df.to_csv(os.path.join(args.wfolder, 'data', token + '_bb_analysis.csv'), index=False)

    # Time and end the sript
    tac = time.time()
    log.update(f'Running time: {round(tac - tic, 1)} sec.')
    log.update('OK')
