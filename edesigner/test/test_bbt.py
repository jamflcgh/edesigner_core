# -*- coding: utf-8 -*-
# test_ebbtcreator
# Jesus Blas 2021

import os
import unittest

from classes.bbt import BBT
from classes.parameter_reader import Parameters
from e_bbt_creator import read_compounds
from classes.logger import Logger
from classes.desalter import SmilesDesalter
from classes.mol_property_calculator import MolPropertyCalculator

fg = Parameters(os.path.join('resources', 'fg.par'), fsource='list', how='to_list', multiple=True)
headpieces = Parameters(os.path.join('resources', 'headpieces.par'), fsource='list', how='to_list',
                        multiple=True)
db = Parameters(os.path.join('resources', 'db.par'), fsource='list', how='to_list', multiple=True)
calcfg = Parameters(os.path.join('resources', 'calcfg.par'), fsource='list', how='to_list', multiple=True)
antifg = Parameters(os.path.join('resources', 'antifg.par'), fsource='list', how='to_list', multiple=True)
log = Logger(os.path.join("logs", "log.log"))
desalter = SmilesDesalter()
calculator = MolPropertyCalculator()

class TestBbt(unittest.TestCase):

    def test_bbt_init(self):
        """Test of BBT instantiator. The instance of BBT is generated for a given BBtype. The script must ensure that
         the attribute BBT has been created, that it has a defined lengh matchng with the fg.par and that a new index
         is generated for the attribute"""
        bbt = BBT(BBT=[0, 0, 1], fg=fg, headpieces=headpieces, index=0)
        self.assertTrue(bbt.BBT is not None, "BBT attribute was not initialized")
        self.assertEqual(len(bbt.BBT_long), len(fg.par), "BBT_long attribute was not initialized")
        self.assertEqual(bbt.index, 0, "index attribute was not initialized")

    def test_bbt_update(self):
        """Test of BBT updater. The instance of BBT is generated for a given BBtype and the df comp is generated from the
        database. This will generate None but will change the attributes of BBT. The script must find some bbt,
        generate an example smiles and ensure that the molecules has at list 1 atom"""
        bbt = BBT(BBT=[0, 0, 1], fg=fg, headpieces=headpieces, index=0)
        comp = read_compounds(db=db, BBTs=[bbt], fg=fg, desalter=desalter, calculator=calculator, calcfg=calcfg,
                              antifg=antifg, log=log)
        bbt.update(comp)
        self.assertTrue(sum(bbt.n_compounds) > 0, "n_compounds has not initialized in bbt.update")
        self.assertTrue(bbt.smiles_example is not None, "smiles_example has not initialized in bbt.update")
        self.assertTrue(bbt.min_atoms > 0, "there in no molecules with n_atoms > 0")
if __name__ == "__main__":
    unittest.main()

