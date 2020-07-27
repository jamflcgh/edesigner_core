# -*- coding: utf-8 -*-
# desalter
# Jose Alfredo Martin

# External modules
from rdkit import Chem

version = 'desalter.v.9.0.0'


class SmilesDesalter:

    def __init__(self, salt_list=None):
        if salt_list is None:
            self.salt_list = ['I', 'Cl', 'Br', '[Cl-]', '[Br-]', '[I-]',
                              'O=S(=O)(O)O', 'Cc1ccc(S(=O)(=O)O)cc1', 'O=S(O)O',
                              'O=S(=O)(O)c1cc(Cl)ccc1C', 'CS(=O)(=O)O', 'O=C(O)O', 'O=P(O)(O)O', 'O',
                              '[NH4+]', '[Li+]', '[Na+]', '[K+]', '[Ca+2]', '[Zn+2]', 'O=[N+]([O-])O',
                              'O=[N+]([O-])c1cc([N+](=O)[O-])c(O)c([N+](=O)[O-])c1', 'O=C(O)C(=O)O', 'O=C(O)[C@H](O)c1ccccc1',
                              'O=C(O)C(O)C(O)C(=O)O', 'O=C(O)[C@H](O)[C@H](O)C(=O)O', 'O=C(O)[C@@H](O)[C@@H](O)C(=O)O',
                              'O=C(O)C(F)(F)F', 'O=CO', 'O=C(O)CO', 'O=C(O)/C=C/C(=O)O', 'O=C(O)/C=C\C(=O)O',
                              'O=C(O)C=CC(=O)O', 'CC(=O)O', 'O=C(O)CCC(=O)O', 'O=C(O)O', 'O=C(O)[C@H](O)[C@@H](O)C(=O)O']
        else:
            self.salt_list = salt_list

    def desalt_smiles(self, smiles=None, mol=None):
        """this method removes all the salts in a smiles or a mol that are contained in a salt dictionary
        if all smiles components are removed or more that two smiles remain it returns None
        smiles: str (smiles string)
        mol : rdkit mol object
        returns : new smiles (str) or None"""
        if mol is None:
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        new_smiles = smiles.split('.')
        if len(new_smiles) > 1:
            if self.salt_list is None:
                new_smiles = None
            else:
                for i in range(len(self.salt_list)):
                    while self.salt_list[i] in new_smiles:
                        new_smiles.remove(self.salt_list[i])
                if len(new_smiles) == 1:
                    new_smiles = new_smiles[0]
                else:
                    new_smiles = None
        else:
            new_smiles = new_smiles[0]
        return new_smiles


if __name__ == '__main__':
    print(version)
