# -*- coding: utf-8 -*-
# libdesign
# Jose Alfredo Martin

# Python modules
import copy
import pickle as pic
import os
# Local Modules
from classes.parallel import starmap_parallel

version = 'libdesign.v.9.0.0'

# CLASS DEFINITION

class LibDesign:
    """LibDesign instances store attributes of a library. A library is a combination of eDESIGNS
    that are compatible (this means that the BBTs it contains can be mixed and react with the same
    enumeration reaction), and also that they share topology (BBTs are attached with the same topology
    and deptortections / scaffold additions are conducted to the topology BBTs). The class also contains
    methods to build and filter out libraries and to generate enumeration instructions and append them
    into the configuration file.
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
        self.design_id = []
        self.n_cycles = len(par.par['max_cycle_na'])
        self.bbts = [[] for i in range(self.n_cycles)]
        self.headpiece = None
        self.eliminate = False
        self.n_int = None
        self.best_index = None
        self.int_limits = None
        self.all_limits = None

    def update_lib(self, design, reaction, deprotection):
        """This method updates the LibDesign instance with data coming from a design that
        matches its lib_id with the lib_id of this instance (this is must be enforced
        upon the call of this mehtod.
        design: instance of Design class
        deprotection: instance of Parameters class (deprotection)
        reaction: instance of Parameters class (reaction)
        returns None:"""
        if self.lib_id is None:  # this happens the first time the library is updated with a design
            self.lib_id = design.lib_id
            self.headpiece = design.bbts[0]
            self.reactions = [reaction.par[design.reactions[i]]['enum_index'] for i in range(1, len(design.reactions))]
            self.deprotections = [deprotection.par[design.deprotections[i]]['enum_index'] for i in range(1, len(design.deprotections))]
            self.scaffold_reactions = [i for i in design.deprotections if deprotection.par[i][
                'atom_dif'] > 0]  # track of deprotections that incorporate scaffolds
            self.bbts = [[design.bbts[i]] for i in range(1, len(design.bbts))]  # BBTs of identical designs will be piled in lists in macrodesigns
        else:
            for i in range(1, len(design.bbts)):
                if design.bbts[i] not in self.bbts[i - 1]:
                    self.bbts[i - 1].append(design.bbts[i])
        self.design_id.append(design.id)
        return None

    def get_best_index_from_indexes(self, distribution, indexes):
        """ This method returns a list of integers and an integer (number of compounds in a library)
        in such a way that the list of integers in the list represent the
        maximum number of atoms in each cycle building block so the number of compounds is maximized
        maintaining the addition of atoms constant
        distribution : list of lists of int
        indexes : list of lists of int
        ndim : int
        neturns : tuple of n_comp, best index (integer, list of integers)"""

        if self.n_cycles == 2:
            n_comp_list = [distribution[0][index[0]] * distribution[1][index[1]] for index in indexes]
            # this is the number of molecules for each number of atoms total
        elif self.n_cycles == 3:
            n_comp_list = [distribution[0][index[0]] * distribution[1][index[1]] * distribution[2][index[2]] for index in indexes]
        n_comp = max(
            n_comp_list)
        # this is the number of compounds of the largest design for this speficic maximum number of atoms
        best_index = indexes[n_comp_list.index(n_comp)]  # this is the index set that gives the largest design
        return n_comp, best_index

    def get_best_index_from_all_indexes(self, distribution, all_indexes, percentile):
        """ This method returns a list of integers and an integer (number of compounds in a library)
        in such a way that the list of integers in the list represent the
        maximum number of atoms in each cycle building block so the number of compounds augmented
        with a percentile, is maximized maintaining the addition of atoms constant
        distribution : list of lists of int
        all_indexes : list of lists of lists of int
        percentile : float
        returns : tuple of n_comp, best index (integer, list of integers)"""
        # calculates the best indexes that maximizes the library size for a number of atoms equal or less than the
        # required at the parameter percentile
        indexes = all_indexes[0]
        best_n_comp, best_index = self.get_best_index_from_indexes(distribution, indexes)
        # Now, since this was just a fraction of the maximum size of the library is calculated
        max_n_comp = int(best_n_comp / float(percentile))
        for indexes in all_indexes[1:]:
            # we increase the max number of atoms one by one by iterating on the rest of all_indexes list and for each
            # # number of atoms calculate the best index and the number of compounds at that best index
            new_n_comp, new_index = self.get_best_index_from_indexes(distribution, indexes)
            # if the number of compounds surpases the max_comp value then the fraction of library below the initial
            # percentile would not be valid and therefore the loop is terminated
            if new_n_comp > max_n_comp:
                break
            else:
                best_n_comp = new_n_comp
                best_index = copy.deepcopy(new_index)
        return best_n_comp, best_index

    def validate_lib(self, BBTs, par, deprotection, global_indexes):
        """This method validates the library and marks it for elimination if it does not meet
        the appropriate criteria
        BBTs : list of instances of BBT class
        par: instance of Parameters class (par)
        deprotection: instance of Parameters class (deprotection)
        returns: None"""
        if len(self.scaffold_reactions) > 0:
            scaffolds_natoms_list = [deprotection.par[item]['atom_dif'] for item in self.scaffold_reactions]
            scaffolds_natoms = sum(scaffolds_natoms_list)
        else:
            scaffolds_natoms = 0
        if scaffolds_natoms < par.par['max_scaffolds_na'] + 1:
            all_indexes = copy.deepcopy(global_indexes[ scaffolds_natoms])  # copy the indexes relevant for the
            # scaffolds of this macrodesign based on the number of atoms for these scaffolds
        else:
            self.eliminate = True
            return None
        # Next creates a distribution of number of molecules based on atom count for each cycle, both for internal and
        # all BBs
        internal_distribution = [[sum([BBTs[bbt].n_internal[m] for bbt in self.bbts[i]]) for m in range(100)] for i in range(self.n_cycles)]
        all_distribution = [[sum([BBTs[bbt].n_compounds[m] for bbt in self.bbts[i]]) for m in range(100)] for i in range(self.n_cycles)]
        # Next finds the best indexes fulfilling the number of compounds
        self.n_int, self.best_index = self.get_best_index_from_all_indexes(internal_distribution, all_indexes, par.par['percentile'])
        if self.n_int < par.par['min_count']:
            # macrodesigns containing a number of compounds, that fulfill the required distribution, less than the
            # threshold are eliminated
            self.eliminate = True
            return None
        # Next we derive the number of BBs that each BBT in a cycle contributes to the design (list of lists, cycle is
        # outer, BBT is inner)
        self.int_limits = [[BBTs[bbt_index].n_internal[self.best_index[i]] for bbt_index in self.bbts[i]] for i in range(self.n_cycles)]
        self.all_limits = [[BBTs[bbt_index].n_compounds[self.best_index[i]] for bbt_index in self.bbts[i]] for i in range(self.n_cycles)]
        # Now we calculate the number of compounds in the macrodesign using internal and external BBs
        self.n_all = 1
        for i in range(self.n_cycles):
            self.n_all *= all_distribution[i][self.best_index[i]]
            # Now we remove the bbts that dont contribute with any compound to the all group (if contribute to all bun
            # not o int we keep it)
            self.bbts = [[self.bbts[i][j] for j in range(len(self.all_limits[i])) if self.all_limits[i][j] > 0] for i in range(self.n_cycles)]
            self.int_limits = [[self.int_limits[i][j] for j in range(len(self.int_limits[i])) if self.all_limits[i][j] > 0] for i in range(self.n_cycles)]
            self.all_limits = [[self.all_limits[i][j] for j in range(len(self.all_limits[i])) if self.all_limits[i][j] > 0] for i in range(self.n_cycles)]
        return None

    def update_translation_file(self, config_path, par, headpieces_dict, enum_reaction, enum_deprotection, deprotection,
                                end_deprotection_enumeration_indexes):
        """This method updates the tranlation to enumeration isntructions file with a new design
        config_path : str (complete path and filename of the file holding the enumeration instructions)
        par: instance of Parameters class (par)
        headpieces_dict: dictionary (of smiles)
        BBTs: list of instaces of BBT class
        enum_reaction : instance of Parameters class (enum_reaction)
        enum_deprotection : instance of Parameters class (enum_deprotection)
        deprotection : instance of Parameters class (deprotection)
        end_deprotection_enumeration_indexes : list of int (enumeration indexes for the end deprotections)
        returns :  None"""
        with open(config_path, 'a') as f:
            for scope in ['INTERNAL', 'ALL']:
                linea = '# Start enumeration instructions\n'
                f.write(linea)
                linea = '# Design number ' + str(self.id) + '\n'
                f.write(linea)
                linea = '# Design fingerprint ' + str(self.lib_id) + '\n'
                f.write(linea)
                linea = '# Design scope ' + str(self.id) + '.' + scope + '\n'
                f.write(linea)
                if scope == 'INTERNAL':
                    linea = '# Design size ' + str(self.n_int) + '\n'
                else:
                    linea = '# Design size ' + str(self.n_all) + '\n'
                f.write(linea)
                linea = '# Design number of cycles ' + str(self.n_cycles) + '\n'
                f.write(linea)
                for i in range(self.n_cycles):
                    if scope == 'INTERNAL':
                        linea = '# MAKE C' + str(i + 1) + '.int.smi WITH {' + ','.join(
                            ['\'' + str(self.bbts[i][j]) + '.int.smi\':' + str(self.int_limits[i][j]) for j in range(len(self.bbts[i]))]) + '}\n'
                    else:
                        linea = '# MAKE C' + str(i + 1) + '.smi WITH {' + ','.join(
                            ['\'' + str(self.bbts[i][j]) + '.smi\':' + str(self.all_limits[i][j]) for j in range(len(self.bbts[i]))]) + '}\n'
                    f.write(linea)
                linea = 'START: ' + headpieces_dict[self.headpiece] + ' core\n'
                f.write(linea)
                for i in range(len(self.reactions)):
                    if self.deprotections[i] != 0:
                        linea = 'AND:\n'
                        f.write(linea)
                        linea = par.par['final_reactions_folder']
                        linea += str(enum_deprotection.par[self.deprotections[i]]['enum_name']) + '\n'
                        f.write(linea)
                        linea = '|\n'
                        f.write(linea)
                    linea = 'AND:\n'
                    f.write(linea)
                    linea = par.par['final_reactions_folder']
                    linea += str(enum_reaction.par[self.reactions[i]]['enum_name'])
                    linea += '||file=' + par.par['final_compounds_folder']
                    if scope == 'INTERNAL':
                        linea += 'C' + str(i + 1) + '.int.smi\n'
                    else:
                        linea += 'C' + str(i + 1) + '.smi\n'
                    f.write(linea)
                    linea = '|\n'
                    f.write(linea)
                # now we include the optional deprotections
                for i, item in enumerate(end_deprotection_enumeration_indexes):
                    # first instance of the deprotection
                    linea = 'AND:\n'
                    f.write(linea)
                    linea = par.par['final_reactions_folder']
                    linea += str(enum_deprotection.par[item]['enum_name']) + '\n'
                    f.write(linea)
                    linea = 'NOOP\n'
                    f.write(linea)
                    linea = '|\n'
                    f.write(linea)
                    # second instance of the deprotection
                    linea = 'AND:\n'
                    f.write(linea)
                    linea = par.par['final_reactions_folder']
                    linea += str(enum_deprotection.par[item]['enum_name']) + '\n'
                    f.write(linea)
                    linea = 'NOOP\n'
                    f.write(linea)
                    if i != len(end_deprotection_enumeration_indexes) - 1:
                        linea = '|\n'
                        f.write(linea)
                linea = '# End enumeration instructions\n'
                f.write(linea)
        return None


#  FUNCTION DEFINITIONS

def get_all_indexes(par, ndim):
    """ this function gets a list of lists of lists with all indexes from a minimum number of atoms to a maximum number of atoms such as they add up to that number of atoms
    the first(outer) level in the list covers all possible atoms in the scaffolds, the second level of the list represents the number of atoms from the threshold
    at the percentile specified in par to the maximum number of atoms in the library. the third level of the list represents all possible combinations of the maximum
    number of atoms in each cycle so that when added and added to the scaffolds atoms and added to the headpiece atoms the total number of atoms is equal to that scpeficified
    at the second level. Finally the fourht level (inner) is a list of the number of atoms for each cycle.
    par: instance of Parameters class (par)
    ndim: number of cycles in the library
    return all_indexes (list of lists of lists of int)"""
    if ndim == 2:
        all_indexes = [[[[i, j] for i in range(1, n) for j in range(1, n) if i + j + ba + par.par['headpiece_na'] == n] for n in range(par.par['max_na_percentile'], par.par['max_na_absolute'])] for ba in range(0, par.par['max_scaffolds_na'] + 1)]
    elif ndim == 3:
        all_indexes = [[[[i, j, k] for i in range(1, n) for j in range(1, n) for k in range(1, n) if i + j + k + ba + par.par['headpiece_na'] == n] for n in range(par.par['max_na_percentile'], par.par['max_na_absolute'])] for ba in range(0, par.par['max_scaffolds_na'] + 1)]
    else:
        return None
    return all_indexes


if __name__ == '__main__':
    print(version)
