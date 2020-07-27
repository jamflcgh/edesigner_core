# -*- coding: utf-8 -*-
# bbt
# Jose Alfredo Martin

version = 'bbt.v.9.0.0'


class BBT:
    """BBT class instances store attributes of a building block type and also methods for its creation
    and modification"""

    def __init__(self, BBT=None, fg=None, headpieces=None, index=None):
        """This method initiallizes the BBT instance.
        BBT : A list of three integers corresponding to the functional group indexes
        fg: an instance of the Parameters class corresponding to the fg parameters
        headpieces : an instance of the Parameters class corresponding to the headpiece parameters
        index : int (BBT identifier)
        returns : None"""
        self.BBT = BBT
        self.BBT_long=[0 for x in range(len(fg.par))]
        for i in self.BBT:
            self.BBT_long[i] += 1
        self.BBT_name = [fg.par[i]['name'] for i in self.BBT]
        self.BBT_multi = sum([0 if i == 0 else 1 for i in self.BBT])
        self.n_compounds = [0 for j in range(100)]
        self.n_internal = [0 for j in range(100)]
        self.n_external = [0 for j in range(100)]
        self.smiles_example = None
        self.headpiece = None
        for item in headpieces.par:
            if self.BBT == item['bbt']:
                self.headpiece = item['index']
        self.index = index
        self.order = 0

    def update(self, comp):
        """This method updates the remaining fields of the BBT using the data extracted in the
        compounds dataframe
        com : pandas dataframe containing classified compounds
        returns : None"""
        df = comp[comp['BBT'] == self.index]
        if df.shape[0] > 0:
            for index, row in df.iterrows():
                self.n_compounds[row['N_ATOMS']] += 1
                if row['EXTERNAL']:
                    self.n_external[row['N_ATOMS']] += 1
                else:
                    self.n_internal[row['N_ATOMS']] += 1
                if index == 0:
                    self.smiles_example = row['SMILES']
            for i in range(1, 100):
                self.n_external[i] += self.n_external[i - 1]
                self.n_internal[i] += self.n_internal[i - 1]
                self.n_compounds[i] += self.n_compounds[i - 1]
            self.min_atoms = 0
            for i in range(100):
                if self.n_compounds[i] > 0:
                    self.min_atoms = i
                    break
        return None


if __name__ == '__main__':
    print(version)
