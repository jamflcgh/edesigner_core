# -*- coding: utf-8 -*-
# libdesign
# Jose Alfredo Martin

# Python modules
import sys
# Local Modules
from classes.bbt import BBT
# external modules
import numpy as np


version = 'libdesign.v.12.0.0'

# CLASS DEFINITION

class LibDesign:
    """LibDesign instances store attributes of a library. A library is a combination of eDESIGNS
    that are compatible (this means that the BBTs it contains can be mixed and react with the same
    enumeration reaction), and also that they share topology (BBTs are attached with the same topology
    and deptortections / scaffold additions are conducted to the topology BBTs). The class also contains
    methods to build and filter out libraries and to generate enumeration instructions and append them
    into the configuration file.
    The topology is embeded in the id vector that has the structure:
    n_cycles: int (number of cycles)
    d0, d1, d2...: (deprotection enumeration indexes at each cycle from 0 to n_cycles -1)
    r1, r2, r3...: (reaction enumeration indexes at the incorporation of each cycle from 1 to n_cycles)
    sd0, sd1, sd2...: (source each deprotection from d0 to d[n_cycles-1]
    sc1, sc2, sc2...: (source of each cycle or reaction, which is the same, from dr1 to r[n_cycles]
    headpiece: int: index of the bbt coding the headpiece
    The posible values for the sources, which are the same for deprotections and cycles (reactions) are:
    0: c0
    1: d0
    2: r1
    3: c1
    4: d1
    5: r2
    6: c2
    7: d2
    8: r3
    9: c3
    10: d3
    ...
    """

    def __init__(self, par):
        """ This method creates an instance of the LibDesign class.
        it requires the par instance of Parameters to create some of the
        attributes of the instance. Initially the lib_id attribute remains as
        None so the instance can be eliminated if there is no design with this lib_id.
        par : instance of the Parameters class (par)
        returns : None"""
        self.id = None
        self.lib_id = None
        self.run_id = None
        self.ed_run_id = None
        self.n_cycles = len(par.par['max_cycle_na'])
        self.bbts = [[] for i in range(self.n_cycles)]
        self.headpiece = None
        self.eliminate = False
        self.n_all = None
        self.best_all_index = None  # list of maximum number of atoms in a building block per cycle for all bbs
        self.all_limits = None  # list of maximum number of bbs selected per cycle from a file that mix int bbts bbs and sorts by number of atoms
        self.all_bbt_limits = None  # list of lists decoupling the number of bbs in self.all_limits by bbt in this design. This facilitates selection of bbs
        self.reactions = None  # list enumeration reacion indexes
        self.deprotections = None  # list of enumeration deprotection reaction indexes
        self.scaffold_reactions = None  # list of deprotections indexes with increased mass (not that this is not from enumeration deprotections) since it it used to calculat atom diff

    def update_lib(self, design, reaction, deprotection, run_id, ed_run_id):
        """This method updates the LibDesign instance with data coming from a design that
        matches its lib_id with the lib_id of this instance (this is must be enforced
        upon the call of this mehtod.
        design: instance of Design class
        deprotection: instance of Parameters class (deprotection)
        reaction: instance of Parameters class (reaction)
        run_id: str: id for the current bbt_creator run
        ed_run_id: str: id for the current eDESIGNER run
        returns None:"""
        if self.lib_id is None:  # this happens the first time the library is updated with a design
            self.lib_id = design.lib_id
            self.headpiece = design.bbts[0]
            self.reactions = [reaction.par[i]['enum_index'] for i in design.reactions]
            self.deprotections = [deprotection.par[i]['enum_index'] for i in design.deprotections]
            self.scaffold_reactions = [i for i in design.deprotections if deprotection.par[i]['atom_dif'] > 0]  # track of deprotections that incorporate scaffolds
            self.bbts = [[design.bbts[i]] for i in range(1, len(design.bbts))]  # BBTs of identical designs will be piled in lists in macrodesigns
            self.run_id = run_id
            self.ed_run_id = ed_run_id
        else:
            if design.lib_id != self.lib_id:
                print('ERROR::: attempted to update a libDESIGN with an dDESINGN where lib_id does not match')
                print(f'current lib_id {self.lib_id}')
                print(f'design lib_id {design.lib_id}')
                sys.exit(1)
            for i in range(1, len(design.bbts)):
                if design.bbts[i] not in self.bbts[i - 1]:
                    self.bbts[i - 1].append(design.bbts[i])
        return None

    def validate_lib(self, BBTs, par, deprotection, na_dist):
        """This method validates the library and marks it for elimination if it does not meet
        the appropriate criteria
        BBTs : list of instances of BBT class
        par: instance of Parameters class (par)
        deprotection: instance of Parameters class (deprotection)
        na_dist: numpy array with as many dimensions of cycles in the library (contains number of atoms in library in
           each position of the array for the combination of number of atoms coming fro each cycle)
        returns: None"""
        if self.n_cycles not in [2, 3]:
            self.eliminate = True
            return None
        if len(self.scaffold_reactions) > 0:
            scaffolds_natoms_list = [deprotection.par[item]['atom_dif'] for item in self.scaffold_reactions]
            scaffolds_natoms = sum(scaffolds_natoms_list)
        else:
            scaffolds_natoms = 0
        if scaffolds_natoms > par.par['max_scaffolds_na']:
            self.eliminate = True
            return None
        na_mdist = na_dist + scaffolds_natoms  # we create a new array because we dont want to change na_dist wich will be used in other lib_designs

        # Next creates a distribution of number of molecules based on atom count for each cycle, both for internal and
        # all BBs
        # first we sum all bbs coming from different BBTs in one single array by number of atoms for each cycle
        all_ncomps = [np.array([BBTs[bbt].n_compounds for bbt in self.bbts[i]]) for i in range(self.n_cycles)]
        all_ncomps = [item.sum(axis=0) for item in all_ncomps]
        all_cum_ncomps = [item.cumsum() for item in all_ncomps]

        # next calculates the matrix distribution of number of compounds
        if self.n_cycles == 2:
            nc_all_dist = all_cum_ncomps[0][:, None] * all_cum_ncomps[1][None, :]
        else:
            nc_all_dist = all_cum_ncomps[0][:, None, None] * all_cum_ncomps[1][None, :, None] * all_cum_ncomps[2][None, None, :]
        best_na_dist = np.where(na_mdist <= par.par['max_na_percentile'])
        if any([len(item) == 0 for item in best_na_dist]):
            self.eliminate = True  # there are not any indexes that fulfil the median
            return None
        perc_all_ncomp = nc_all_dist[best_na_dist].max()
        max_all_ncomp = int(perc_all_ncomp / par.par['percentile'])
        if max_all_ncomp < par.par['min_count']:
            self.eliminate = True  # if lib design does not meet ncompound criteria it is eliminated before computing the best index
            return None
        self.n_all = nc_all_dist[np.where(nc_all_dist <= max_all_ncomp)].max()
        best_all_indexes = np.where(nc_all_dist == self.n_all)
        # From all possible best indexes, we select the one that ensures the most homogenous number of compounds per cycle
        # and on tying, the one that gives more diversity in the las cycles
        if len(best_all_indexes[0]) > 0:
            self.best_all_index = [item[0] for item in best_all_indexes]  # Initialize the best index for all
            for i in range(len(best_all_indexes[0])):  # this is clunky because it serves for both two and three cycles
                this_index_nc = np.array([all_cum_ncomps[j][best_all_indexes[j][i]] for j in range(len(best_all_indexes))])
                best_index_nc = np.array([all_cum_ncomps[j][self.best_all_index[j]] for j in range(len(best_all_indexes))])
                if this_index_nc.std() < best_index_nc.std():
                    self.best_all_index = [item[i] for item in best_all_indexes]
            self.all_limits = [all_cum_ncomps[j][self.best_all_index[j]] for j in range(len(best_all_indexes))]
            self.all_bbt_limits = [[BBTs[bbt_index].n_compounds[self.best_all_index[i]] for bbt_index in self.bbts[i]] for i in range(self.n_cycles)]

        return None

    def print_summary_file(self, enum_reaction, enum_deprotection):
        """prints a summary of the library desing to be understood by a human
        enum_reaction: instance of Par class coding the enum_reaction parameters
        enum_deprotection: instance of Par class coding the enum_deprotection parameters
        returns: None"""
        deprotections = [self.lib_id[i] for i in range(1, self.n_cycles + 1)]
        reactions = [self.lib_id[i] for i in range(1 + self.n_cycles, 2 * self.n_cycles + 1)]
        print()
        print('LIBRARY DESCRIPTION:')
        print(f'Library_id: {self.lib_id}')
        print(f'Run id: {self.run_id}')
        print(f'eDESIGNER run id: {self.ed_run_id}')
        print(f'Headpiece_smiles: {self.headpiece}')
        for i, (deprotection, reaction) in enumerate(zip(deprotections, reactions)):
            if deprotection > 0:
                print(f'Add Deprotection at cycle {i}: {enum_deprotection.par[deprotection]["enum_name"]}')
            print(f'Add cycle {i+1} through reaction: {enum_reaction.par[reaction]["enum_name"]}')
        print('')


if __name__ == '__main__':
    print(version)
