# -*- coding: utf-8 -*-
# mol_property_calculator
# Jose Alfredo Martin

version = 'mol_property_calculator.v.9.0.0'

# External modules
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as desc


class MolPropertyCalculator:
    def __init__(self):
        self.data = {}
        self.success = True

    def calculate_properties(self, smiles=None, mol=None, props=[]):
        """this method calculates basic properties for the mol
        returns : error (bool)"""
        if len(props) == 0:
            return True
        if mol is None:
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return True
        if 'py_formula' in props:
            self.data['py_formula'] = desc.CalcMolFormula(mol)
        if 'py_em' in props:
            self.data['py_em'] = round(desc.CalcExactMolWt(mol), 5)
        if 'py_n_Cl_Br' in props:
            all_atoms = []
            for atom in mol.GetAtoms():
                all_atoms.append(atom.GetSymbol())
            n_Cl = all_atoms.count('Cl')
            n_Br = all_atoms.count('Br')
            self.data['py_n_Cl_Br'] = n_Cl + n_Br
        if 'py_na' in props:
            self.data['py_na'] = mol.GetNumAtoms()
        if 'py_mw' in props:
            self.data['py_mw'] = desc._CalcMolWt(mol)
        if 'py_fsp3' in props:
            self.data['py_fsp3'] = desc.CalcFractionCSP3(mol)
        if 'py_rb' in props:
            self.data['py_rb'] = desc.CalcNumRotatableBonds(mol)
        if 'py_tpsa' in props:
            self.data['py_tpsa'] = desc.CalcTPSA(mol)
        if 'py_clogp' in props:
            self.data['py_clogp'] = desc.CalcCrippenDescriptors(mol)[0]
        if 'py_nar' in props:
            self.data['py_nar'] = desc.CalcNumAromaticRings(mol)
        if 'py_nhba' in props:
            self.data['py_nhba'] = desc.CalcNumHBA(mol)
        if 'py_nhbd' in props:
            self.data['py_nhbd'] = desc.CalcNumHBD(mol)
        return False


if __name__ == '__main__':
    print(version)
