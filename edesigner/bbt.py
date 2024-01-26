# -*- coding: utf-8 -*-
# bbt
# Jose Alfredo Martin

__version__ = 'bbt.v.12.0.0'
__author__ = 'Alfredo Martin 2023'

import numpy as np

class BBT:
    """BBT class instances store attributes of a building block type and also methods for its creation
    and modification"""

    def __init__(self, BBT=None, fg=None, headpieces=None, index=None, maxna=None):
        """This method initiallizes the BBT instance.
        BBT : A list of three integers corresponding to the functional group indexes
        fg: an instance of the Parameters class corresponding to the fg parameters
        headpieces : an instance of the Parameters class corresponding to the headpiece parameters
        index : int (BBT identifier)
        maxna : int (max number of atoms in a given building block
        returns : None"""
        self.maxna = maxna  # max number of atoms a bb assigned to an instance of this class can have
        self.BBT = BBT
        self.BBT_long = [0 for _ in range(len(fg.par))]
        for i in self.BBT:
            self.BBT_long[i] += 1
        self.BBT_name = [fg.par[i]['name'] for i in self.BBT]
        self.BBT_multi = sum([0 if i == 0 else 1 for i in self.BBT])
        self.n_compounds = np.array([0 for _ in range(maxna + 1)])  # todo new for testing
        self.smiles_example = None
        self.headpiece = None
        self.min_atoms = 0  # bb belonging to this bbt with min number of effective atoms
        self.max_atoms = 0  # bb belonging to this bbt with max number of efective atoms
        for item in headpieces.par:
            if self.BBT == item['bbt']:
                self.headpiece = item['index']
        self.index = index
        self.order = 0


if __name__ == '__main__':
    print(__version__)
    print(__author__)
