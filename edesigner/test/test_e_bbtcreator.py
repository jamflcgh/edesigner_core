import unittest
import sys
import os
import shutil
from classes.bbt import BBT
from classes.logger import Logger
from classes.parameter_reader import Parameters
from classes.bb_reader import BBReader
from e_bbt_creator import generate_bbts

LOGFOLDER = os.path.join(os.environ['EDESIGNER_TEST_FOLDER'], 'logs')
RUNFOLDER = os.path.join(os.environ['EDESIGNER_TEST_FOLDER'], 'results', 'R000000000')
print(RUNFOLDER)
print(os.environ['EDESIGNER_TEST_FOLDER'])
COMPSFOLDER = os.path.join(RUNFOLDER, 'comps')
RESULTSFOLDER = os.path.join(RUNFOLDER, 'results')
RESOURCESFOLDER = os.path.join(RUNFOLDER, 'resources')
os.mkdir(RUNFOLDER)
os.mkdir(COMPSFOLDER)
os.mkdir(RESULTSFOLDER)
os.mkdir(RESOURCESFOLDER)
PARFOLDER = os.path.join(os.environ['EDESIGNER_TEST_FOLDER'], '..', 'resources')

log = Logger(os.path.join(LOGFOLDER, 'e_bbt_creator.log'), prepend_timestamp=False)
dbpar = Parameters(os.path.join(PARFOLDER, 'db.par'), fsource='list', how='to_list', multiple=True)
for i in range(len(dbpar.par)):
    dbpar.par[i]['filename'] = os.path.expandvars(dbpar.par[i]['filename'])
bblim = Parameters(os.path.join(PARFOLDER, 'bblim.par'), fsource='dict', how='to_dict', multiple=False)
par = Parameters(os.path.join(PARFOLDER, 'par.par'), fsource='dict', how='to_dict', multiple=False)
fg = Parameters(os.path.join(PARFOLDER, 'fg.par'), fsource='list', how='to_list', multiple=True)
calcfg = Parameters(os.path.join(PARFOLDER, 'calcfg.par'), fsource='list', how='to_list', multiple=True)
antifg = Parameters(os.path.join(PARFOLDER, 'antifg.par'), fsource='list', how='to_list', multiple=True)
headpieces = Parameters(os.path.join(PARFOLDER, 'headpieces.par'), fsource='list', how='to_list', multiple=True)
for item in [dbpar, par, fg, calcfg, antifg, headpieces]:
    if item.success is None:
        error_found = True
        log.update('    Error(s) found while reading parameters: ' + item.path)
        for error in item.errors:
            log.update('        ' + error)

BBTs = generate_bbts(fg, headpieces, bblim, log)


class MyTestCase(unittest.TestCase):
    def test_run_bbt_creator(self):
        """test_run_bbt_creator"""
        reader = BBReader(RESULTSFOLDER, RUNFOLDER, BBTs, dbpar, fg, antifg, calcfg, bblim, log,
                          smi=None, verbose=False, debug=False)
        reader.run()
        success = os.path.isfile(os.path.join(RESULTSFOLDER, 'bbt_report.csv'))
        shutil.rmtree(RUNFOLDER)
        self.assertTrue(success, "e_bbt_creator has error(s)")

if __name__ == '__main__':
    unittest.main()
