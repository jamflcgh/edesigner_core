# -*- coding: utf-8 -*-
# sample_enumerator
# Jose Alfredo Martin 2020

version = 'sample_enumerator.v.9.0.0'

# Python modules
import random
import copy
# External modules
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors as desc


class SampleEnumerator:
    def __init__(self, n_cycles=None, reactions=None, deprotections=None):
        """Instance constructor
        n_cycles: int: number of cycles in the library in the range[1, 3]. reactions and deprotections len must match n_cycles
        reactions: list of str: list of smirks reaction deffinitions (order of reactants on-DNA, off-DNA)
        deprotections: list of str or None: list of smirks deprotection definitions"""
        assert type(n_cycles) == int and 0 < n_cycles < 4, 'n_cycles must be an int between 1 and 3'
        assert type(reactions) == list and not None in reactions, 'reactions must be a list of smirks, None is not allowed'
        assert type(deprotections) == list, 'deprotections must be a list of smirks, if no deprotection is to be performed in a cycle replace the smirks by None'
        assert len(reactions) == n_cycles, 'len of reactions must match n_cycles'
        assert len(deprotections) == n_cycles, 'len of deprotections must match n_cycles'
        self.n_cycles = n_cycles
        self.reactions = [AllChem.ReactionFromSmarts(reaction) for reaction in reactions]
        self.deprotections = [AllChem.ReactionFromSmarts(deprotection) if deprotection is not None else None for deprotection in deprotections]
        final_deprotections = ['[#6:1][C:2](=[O:3])[O:4][C;$(C),$(CC):5]>>[#6:1][C:2](=[O:3])[O:4]',
                               '[#7;!$([#7](C=O)C=O):1][C:2](=[O:3])[O:4][C:5]([C:6])([C:7])[C:8]>>[#7:1]',
                               '[#7;!$([#7](C=O)C=O);$([#7]C(=O)OCC1c2ccccc2c3ccccc13):1][C:2](=[O:3])[O:4]>>[#7:1]']
        self.final_deprotections = [AllChem.ReactionFromSmarts(final_deprotection) for final_deprotection in final_deprotections]
        self.success = True

    def enumerate_molecule(self, headpiece, bbs):
        """enumerates a single molecule given the list of the bb used in each cycle
        hedpiece: str: smiles string representing the headpiece
        bbs: list str: list of building blocks, one per cycle"""
        hmol = Chem.MolFromSmiles(headpiece)
        mols = [Chem.MolFromSmiles(bb) for bb in bbs]
        mol = copy.deepcopy(hmol)
        for i in range(self.n_cycles):
            if self.deprotections[i] is not None:
                this_mol = self.deprotections[i].RunReactants((mol,))
                if len(this_mol) == 0:
                    return None
                this_mol = this_mol[0][0]
                try:
                    Chem.SanitizeMol(this_mol)
                except:
                    return None
                mol = this_mol
            this_mol = self.reactions[i].RunReactants((mol, mols[i]))
            if len(this_mol) == 0:
                return None
            this_mol = this_mol[0][0]
            try:
                Chem.SanitizeMol(this_mol)
            except:
                return None
            mol = this_mol
        for deprotection in self.final_deprotections:
            this_mol = deprotection.RunReactants((mol,))
            if len(this_mol) == 0:
                pass
            else:
                try:
                    Chem.SanitizeMol(this_mol)
                except:
                    pass
                else:
                    mol = this_mol
        return mol

    def calculate_properties(self, mol):
        """this method calculates basic properties for the smiles
        returns : list of int or float (properties)"""
        properties = []
        properties.append(mol.GetNumAtoms())
        properties.append(desc.CalcCrippenDescriptors(mol)[0])
        properties.append(desc.CalcTPSA(mol))
        properties.append(desc.CalcNumRotatableBonds(mol))
        properties.append(desc.CalcFractionCSP3(mol))
        return properties

    def enumerate_set(self, headpiece, bbs, n_samples=1, seed=1, props=True):
        """Enumerates a set of molecules
        headpiece: str: smiles string used to represent the headpiece
        bbs: list of lists of str: one list per cycle containing a list of bb smiles
        n_samples: int: number of samples to be enumerated (default 1)
        seed: int: seed used for the random generator (default 1)
        probs: bool: whether to add calculated properties (na, clogp, tpsa, nrb, fcsp3) (default True)
        returns: tuple of list of smiles strings, list of lists of float opr integer (molecules, properties)"""
        random.seed(seed)
        indexes = []
        properties = []
        for _ in range(n_samples):
            indexes.append([random.randrange(0, len(bbs[j])) for j in range(self.n_cycles)])
        molecules = []
        for i in range(n_samples):
            mol = self.enumerate_molecule(headpiece, [bbs[j][indexes[i][j]] for j in range(self.n_cycles)])
            if mol is not None:
                molecules.append(Chem.MolToSmiles(mol, isomericSmiles=True))
                if props:
                    properties.append(self.calculate_properties(mol))
                else:
                    properties.append(None)
        return molecules, properties

    def enumerate_sets(self, headpiece, bbs, n_sets=1, n_samples=1, props=True):
        """enumerates one or mores sets corresponding to a different random sequence each and return all of them in a
        single pandas dataframe
        headpiece: str: smiles string used to represent the headpiece
        bbs: list of lists of str: one list per cycle containing a list of bb smiles
        nsets: int: number of set to be enumerated
        n_samples: int: number of samples to be enumerated (default 1)
        probs: bool: whether to add calculated properties (na, clogp, tpsa, nrb, fcsp3) (default True)
        returns: pandas dataframe
        """
        df = []
        for i in range(n_sets):
            molecules, properties = self.enumerate_set(headpiece, bbs, n_samples=n_samples, seed=i, props=props)
            na = [prop[0] for prop in properties]
            clogp = [prop[1] for prop in properties]
            tpsa = [prop[2] for prop in properties]
            nrb = [prop[3] for prop in properties]
            fcsp3 = [prop[4] for prop in properties]
            this_df = pd.DataFrame({'smiles': molecules,
                                    'na': na,
                                    'clogp': clogp,
                                    'tpsa': tpsa,
                                    'nrb': nrb,
                                    'fcsp3': fcsp3,
                                    'set': [i + 1 for _ in range(len(molecules))]})
            df.append(this_df)
        df = pd.concat(df)
        return df


if __name__ == '__main__':
    print(version)
