# -*- coding: utf-8 -*-
# parallel
# Jose Alfredo Martin

# Python modules
from multiprocessing import Pool
from multiprocessing import cpu_count
# External modules
import numpy as np
import pandas as pd

version = 'parallel.v.9.0.0'

def map_parallel(data, func, cores=-1, partitions=-1, n_threshold=1000, is_list=False):
    """This function parallelizes the processing of data in a pandas dataframe with a specific funcion func that
    applies changes to the dataframe
    data : pandas dataframe
    func : function to process the whole dataframe
    cores : number of cores to process the data (-1 means all cores)
    partitions : number of partitions to process the data (at leas the number of cores)
    n_threshold : (int) number of items in data that triguer paralelization or not
    is_list : (bool) whather the structure of the data is a list or not, if not it is assumed a pandas dataframe
    returns : out_data (pandas dataframe)"""
    # the instruction pool = Pool(cores) takes about 1.2 seconds, so parallelize is useful only
    # if the whole process takes more than 1.5 to 2 seconds. Time your routine and use n_threshold to
    # tune the desired behaviour
    if len(data) > n_threshold:
        if cores == -1:
            cores = cpu_count()
        if partitions == -1:
            partitions = cores
        data_split = np.array_split(data, partitions)
        pool = Pool(cores)
        if is_list:
            out_data = list(np.concatenate(pool.map(func, data_split)))
        else:
            out_data = pd.concat(pool.map(func, data_split), sort=True)
        pool.close()
        pool.join()
    else:
        out_data = func(data)
    return out_data


def starmap_parallel(data, func, cores=-1, partitions=-1, is_list=False, args=None):
    """This function parallelizes the processing of data in a pandas dataframe with a specific funcion func that
    applies changes to the dataframe
    data : pandas dataframe
    func : function to process the whole dataframe
    cores : number of cores to process the data (-1 means all cores)
    partitions : number of partitions to process the data (at leas the number of cores)
    is_list : (bool) whether the structure of the data is a list or not, if not it is assumed a pandas dataframe
    returns : out_data (pandas dataframe)"""
    if cores == -1:
        cores = cpu_count()
    if partitions == -1:
        partitions = cores
    if partitions > len(data):
        partitions = cores = len(data)
    data_split = np.array_split(data, partitions)
    data_split = [[item] + list(args) for item in data_split]
    pool = Pool(cores)
    if is_list:
        out_data = list(np.concatenate(pool.starmap(func, data_split)))
    else:
        out_data = pd.concat(pool.starmap(func, data_split), sort=True)
    pool.close()
    pool.join()

    return out_data


if __name__ == '__main__':
    print(version)
