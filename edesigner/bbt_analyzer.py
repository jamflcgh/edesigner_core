# -*- coding: utf-8 -*-
# lib_design_interpreter
# Jose Alfredo Martin 2020

version = 'bbt_analyzer.v.9.0.0'

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
# other modules
import numpy as np
import pandas as pd

def process_bbts(BBTs):
    """extracts key information from BBTs in the form of np arrays
    BBTs: list of BBT class instance
    returns pBBTs: list of dicts"""
    pBBTs = []
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
            pBBTs.append(this_dict)
    return pBBTs


def get_bbs_by_natoms(pBBTs):
    """gets a series of numpy arrays with an analysis of number of builing blocks by multiplicity and source
    and cretes an unpivoted pandas dataframe.
    BBTs: list of dictionaries
    plot: bool (whether to plot the results)
    returns: df: pandas dataframe"""
    n_atoms = np.array(list(range(100)))
    nbbs_mono_int = np.array([0 for _ in range(100)])
    nbbs_mono_ext = np.array([0 for _ in range(100)])
    nbbs_bi_int = np.array([0 for _ in range(100)])
    nbbs_bi_ext = np.array([0 for _ in range(100)])
    nbbs_tri_int = np.array([0 for _ in range(100)])
    nbbs_tri_ext = np.array([0 for _ in range(100)])
    max_atoms = 0
    min_atoms = 100
    for bbt in pBBTs:
        max_atoms = max(max_atoms, bbt['max_atoms'])
        min_atoms = min(min_atoms, bbt['min_atoms'])
        if bbt['multi'] == 1:
            nbbs_mono_int += bbt['hist_int']
            nbbs_mono_ext += bbt['hist_ext']
        if bbt['multi'] == 2:
            nbbs_bi_int += bbt['hist_int']
            nbbs_bi_ext += bbt['hist_ext']
        if bbt['multi'] == 3:
            nbbs_tri_int += bbt['hist_int']
            nbbs_tri_ext += bbt['hist_ext']

    n_atoms = n_atoms[min_atoms:max_atoms]
    nbbs_mono_int = nbbs_mono_int[min_atoms: max_atoms]
    nbbs_mono_ext = nbbs_mono_ext[min_atoms: max_atoms]
    nbbs_mono_all = nbbs_mono_int + nbbs_mono_ext
    nbbs_bi_int = nbbs_bi_int[min_atoms: max_atoms]
    nbbs_bi_ext = nbbs_bi_ext[min_atoms: max_atoms]
    nbbs_bi_all = nbbs_bi_int + nbbs_bi_ext
    nbbs_tri_int = nbbs_tri_int[min_atoms: max_atoms]
    nbbs_tri_ext = nbbs_tri_ext[min_atoms: max_atoms]
    nbbs_tri_all = nbbs_tri_int + nbbs_tri_ext
    df = pd.DataFrame({'number atoms': list(n_atoms) * 9})
    df['multiplicity'] = [1] * 3 * n_atoms.shape[0] + [2] * 3 * n_atoms.shape[0] + [3] * 3 * n_atoms.shape[0]
    df['source'] = (['ALL'] * n_atoms.shape[0] + ['INTERNAL'] * n_atoms.shape[0] + ['EXTERNAL'] * n_atoms.shape[0]) * 3
    data = list(nbbs_mono_all) + list(nbbs_mono_int) + list(nbbs_mono_ext)
    data += list(nbbs_bi_all) + list(nbbs_bi_int) + list(nbbs_bi_ext)
    data += list(nbbs_tri_all) + list(nbbs_tri_int) + list(nbbs_tri_ext)
    df['number bbs'] = data
    return df

def get_bbts_by_bbtype(pBBTs, minbbs=0):
    """Does an analysis of number of bbs by bbt and returns a pandas dataframe with the analysis
    BBTs: list of dictionaries
    minbbs: remove bbts with less than minbbs from the analysis (default is 0)
    returns None"""
    bbt_source =[]
    bbt_multi = []
    bbt_nbbs = []
    bbt_index = []
    bbt_name = []
    for bbt in pBBTs:
        if bbt['n_all'] > minbbs:
            bbt_source += ['ALL', 'INTERNAL', 'EXTERNAL']
            bbt_multi += [bbt['multi'] for _ in range(3)]
            bbt_nbbs += [bbt['n_all'], bbt['n_int'], bbt['n_ext']]
            bbt_index += [bbt['index'] for _ in range(3)]
            bbt_name += [bbt['name'] for _ in range(3)]
    df = pd.DataFrame({'bbt index': bbt_index})
    df['bbt name'] = bbt_name
    df['source'] = bbt_source
    df['multiplicity'] = bbt_multi
    df['number bbs'] = bbt_nbbs
    df = df[df['number bbs'] > 0]
    df['log(number bbs)'] = np.log10(df['number bbs'])
    df.sort_values(inplace=True, by=['multiplicity', 'number bbs'], ascending=True)
    dfall = df[df['source'] == 'ALL'].copy()
    dfall['order'] = list(range(dfall.shape[0]))
    dfint = df[df['source'] == 'INTERNAL'].copy()
    dfint['order'] = list(range(dfint.shape[0]))
    dfext = df[df['source'] == 'EXTERNAL'].copy()
    dfext['order'] = list(range(dfext.shape[0]))
    df = pd.concat([dfall, dfint, dfext])
    return df

if __name__ == '__main__':

    # Args parser
    parser = argparse.ArgumentParser(description="""This script runs an analysis on the number of building blocks by
    number of atoms and multiplicity and by bbt index and multiplicity. The source is the list of BBT onjects generated
    by bbt_creator and stored in a piclked file. The file is accessed by the timestamp passed as an argument to the
    script. This file is stored in the folder data relative to wfolder, which is also passed as argument.
    The resul of the analysis is stored in two files in the same data folder. The name of the files is defined also
    by the argument timestamp.
    used in the selected libraries by libDESIGNER. The wfolder argument indicates where to look for the souce file
    (a file containing a pickled libDESIGN instance in each record) the file must be in a folder named results within 
    wfolder. The timestamp and tag arguments define the name of the soruce file and the name of the output files.
    Output files are written in a folder named data within wfolder""")
    parser.add_argument('-wf', '--wfolder', help='Working Folder (default folder from which this script is run)',
                        type=str, default='./')
    parser.add_argument('-ts', '--timestamp', help='date in which the bbt_creator was run (format is YYYYMMDD)'
                        , type=str, default=None)
    args = parser.parse_args()

    # Initialization
    tic = time.time()
    log = Logger(os.path.join(args.wfolder, 'logs', args.timestamp + '_bbt_analyzer.log'))
    log.update(version)
    log.update('reading BBTs...')
    with open(os.path.join(args.wfolder, 'data', args.timestamp + '_BBTs.pic'), 'rb') as f:
        BBTs = pic.load(f)
    pBBTs = process_bbts(BBTs)
    log.update('performing nbbs by natoms and multi analysis...')
    df = get_bbs_by_natoms(pBBTs)
    df.to_csv(os.path.join(args.wfolder, 'data', args.timestamp + '_nbbs_by_natoms_and_multi_analysis.csv'), index=False)
    log.update('performing nbbs by natoms and multi analysis...')
    df = get_bbts_by_bbtype(pBBTs, minbbs=0)
    df.to_csv(os.path.join(args.wfolder, 'data', args.timestamp + '_nbbs_by_bbt_and_multi_analysis.csv'), index=False)

    # Time and end the sript
    tac = time.time()
    log.update(f'Running time: {round(tac - tic, 1)} sec.')
    log.update('OK')
