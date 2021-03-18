import unittest
import sys
import os
import e_bbt_creator


wfolder = os.path.abspath(os.path.dirname(sys.argv[0]))
"""It gives you back a list with the arguments of this script: path of the script"""
subfolders = os.listdir(os.path.join(wfolder, 'comps'))
"""This will delete the folder comps each time the test is run """
for folder in subfolders:
    for file in os.listdir(os.path.join(wfolder, 'comps', folder)):
        os.remove(os.path.join(wfolder, 'comps', folder, file))
    os.rmdir(os.path.join(wfolder, 'comps', folder))


class MyTestCase(unittest.TestCase):
    def test_initialization(self):
        """This function test the initialization of the working folder in the test environment. This initialization will
            create new directories"""

        tic, path, log, dbpar, par, fg, calcfg, antifg, headpieces, desalter, calculator, comps_path = e_bbt_creator.intialization(wfolder)
        self.assertTrue(dbpar is not None, "the script has not run")
        self.assertTrue(comps_path is not None, "the script has been initialized")

    def test_generate_bbts(self):
        """"""

        tic, path, log, dbpar, par, fg, calcfg, antifg, headpieces, desalter, calculator, comps_path = e_bbt_creator.intialization(wfolder)
        BBTs = e_bbt_creator.generate_bbts(fg, headpieces, log)
        self.assertEqual(len(BBTs),6815, "the number of BBTs in not correct")

    def test_read_compounds(self):
        """"""

        tic, path, log, dbpar, par, fg, calcfg, antifg, headpieces, desalter, calculator, comps_path = e_bbt_creator.intialization(wfolder)
        BBTs = e_bbt_creator.generate_bbts(fg, headpieces, log)
        comp = e_bbt_creator.read_compounds(db=dbpar, BBTs=BBTs, fg=fg, desalter=desalter, calculator=calculator, calcfg=calcfg,
                              antifg=antifg, log=log)
        self.assertEqual(comp.shape, (865,8), "the shape of comps df is not correct")


if __name__ == '__main__':
    unittest.main()
