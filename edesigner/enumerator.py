# -*- coding: utf-8 -*-
# enumerator
# Jose Alfredo Martin

# Python modules
import copy
import pandas as pd
import os
import shutil
import sys
import json
# Local Modules
from classes.synth_graph import SynthGraph



version = 'libdesign.v.12.0.0'

# CLASS DEFINITION


class Enumerator:
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

    def __init__(self, wfolder, lib_id, hp_smiles, user, bbs=None, lib=None, base_folder=None, verbose=False):
        """
        Initiallizes the instance and writes the building blocks files in the working folder
        :param wfolder: str: path to enumeration folder
        :param lib_id: tuple: library id
        :param hp_smiles: str: smiles of the headpiece
        :param user: str: bool: whether librabry comes from user or from eDESIGNER
        :param bbs: list of str: path to the files containing building blocks to enumerate or None
        :param lib: instance of libDEDIGN or None
        :param base_folder: str: base_folder for the run in which the lib was generated or None
        :param verbose: str: bool: whether to provide additional information on the console
        """
        self.wfolder = wfolder
        self.lib_id = lib_id
        self.hp_smiles = hp_smiles
        self.n_cycles = lib_id[0]
        self.user = user
        self.verbose = verbose
        self.bbs = bbs
        self.lib = lib
        self.base_folder = base_folder

    def write_bbs_files(self):
        """
        Writes the building block files into the working folder. Picks the bbs from bbs if passed and if not from lib,
        but one of them is required
        :return: None
        """
        if self.bbs is not None:
            if len(self.bbs) != self.n_cycles:
                print(f'ERROR::: number of bbs files passed must match the number of cycles in the library excluding headpiece')
                with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                    f.write(f'ERROR::: number of bbs files passed must match the number of cycles in the library excluding headpiece\n')
                sys.exit(1)
        if self.bbs is None:
            if self.base_folder is None or self.lib is None:
                print(f'ERROR::: if not bbs are passed then a libDESIGN and a base folder must be passed')
                with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                    f.write(f'ERROR::: if not bbs are passed then a libDESIGN and a base folder must be passed\n')
                sys.exit(1)
        with open(os.path.join(self.wfolder, 'C0.smi'), 'w') as f:
            f.write(self.hp_smiles + ' headpiece\n')
        for i in range(self.n_cycles):
            cycle = i + 1
            out_bbs_file = f'C{cycle}.smi'
            if self.bbs is not None:
                shutil.copy(self.bbs[i], os.path.join(self.wfolder, out_bbs_file))
            else:
                for j in range(len(self.lib.bbts[i])):
                    bbs_file = str(self.lib.bbts[i][j]) + ".smi"
                    command = f'head -{self.lib.all_bbt_limits[i][j]} '
                    command += f'{os.path.join(self.base_folder, self.lib.run_id, "comps", bbs_file)} '
                    command += f'>> {os.path.join(self.wfolder, out_bbs_file)}'
                    os.system(command)

    def get_valences(self, mrd, preparations):
        """
        get_valences extracts the number of connecting edges and the number of cyclic edges to be used in each node
        involved in the sub-synthesis represented by mrd
        :param mrd: dictionary containing the multireaction entry for the current sun-synthesis
        :param preparations: instance of Parameters class containing preparations parameters
        :return: tuple of linking_valences, cyclic_valences
           linking_valences: list of number of linking edges for each node in the sequence in the order it incorporates
           cyclic_valences: list of number of linking edges for each node in the sequence in the order it incorporates
        """
        valences = []
        for i, prep_idx in enumerate(mrd['preparations']):
            if i == mrd['out_isotopes'] - 1:
                out_isotope = True
            else:
                out_isotope = False
            valences.append({'idx': i,
                             'out_isotope': out_isotope,
                             'n_edges': preparations.par[prep_idx]['n_edges']})
        l_valences = [item['n_edges'] for item in valences]
        c_valences = [0 for _ in range(len(valences))]
        lv = [l_valences[i] + 1 if valences[i]['out_isotope'] else l_valences[i] for i in range(len(l_valences))]
        ilv = [item for item in lv]
        bonds = []
        counter = 0
        while True:
            counter += 1
            for i in range(len(lv) - 1):
                for j in range(i + 1, len(lv)):
                    if lv[i] > 0 and lv[j] > 0:
                        lv[i] -= 1
                        lv[j] -= 1
                        if (i, j) in bonds:
                            c_valences[i] += 1
                            c_valences[j] += 1
                        bonds.append((i, j))
                        break
            if sum(lv) == 0:
                break
            if counter > 20 or sum(lv) < 0:
                print("ERROR::: too many cycles to find bonds")
                print(mrd)
                print('bonds:', bonds)
                print('remaining valences', lv)
                print('initial valences', ilv)
                print('l_valences', l_valences)
                print('c_valences', c_valences)
                with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                    f.write("ERROR::: too many cycles to find bonds\n")
                    f.write(f'{mrd}\n')
                    f.write(f'bonds: {bonds}\n')
                    f.write(f'remanining valences: {lv}\n')
                    f.write(f'initial valences: {ilv}\n')
                    f.write(f'l_valences: {l_valences}\n')
                    f.write(f'c_valences: {c_valences}\n')
                sys.exit(1)
        if sum(c_valences) > 2:
            print("ERROR::: Only one cycle edge is allowed per sub-synthesis")  # todo is this really required?
            with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                f.write("ERROR::: Only one cycle edge is allowed per sub-synthesis")
            sys.exit(1)
        l_valences = [item[0] - item[1] for item in zip(l_valences, c_valences)]
        return l_valences, c_valences

    def get_sub_synthesis(self, s_reactions, nodes, last_cycle_edge, multireaction, preparations, enum_deprotection):
        """sub_synthesis is an object coded as a list of dictionaries. It comprises all the information for a simple
        (single reaction) or multiple (multi reaction) node connecting systems. The list of dictionaries
        (one per node in the sub-synthesis) is to be used to tranlate dicectly into a json config file that defines
        an isotopic graph enumeration system to enumerate a library.
        The keys of the nodes objects comprising the sub-synthesis list has the same keys of the final config file:
        node, bbs_files, edges, out_edges, bond_orders, cycle_edges, cycle_bond_orders and preparations.
        Preparations is a list of dictionaries, each one having the following keys: edges, keep unprepared, reaction.
        Parameters:
        s_reactions: tuple of str. Defines the sub-synthesis
        nodes: list of int, nodes involved in this sub-synthesis, in increasing order
        last_cycle_edge: int, index of the last cycle edge used, if no cycle edges used then it would be the last cycle
        multireactions: instance of Parameters class containing multireactions parameters
        preparations: instance of Parameters class containing preparations parameters
        enum_deprotection: instance of Parameters class containing enum_deprotection parameters
        """
        used_cycle_edges = 0  # to keep track of how many cycle edges are found
        sub_synthesis = list()
        # esr is the reaction sequence eliminating the initial deprotection (if any)
        d0 = []
        esr = copy.deepcopy(s_reactions)
        while True:
            if esr[0].startswith('d'):  # remove the first de-protection for matching if it exists
                d0.append(esr[0])
                esr = [item for item in esr[1:]]
                if len(esr) == 0:
                    print(f"ERROR esr had no reactions")
            else:
                break
        esr = tuple(esr)
        for mrd in multireaction.par:  # mrd is the multi-reaction dictionary that matches the esr
            if tuple(mrd['s_reactions']) == esr:
                break
        else:
            with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                f.write(f'ERROR::: No match for this sub-synthesis {s_reactions} in multireactions parameters\n')
            print(f'ERROR::: No match for this sub-synthesis {s_reactions} in multireactions parameters')
            sys.exit(1)

        linking_valences, cyclic_valences = self.get_valences(mrd, preparations)

        n_edges = len([r for r in esr if r.startswith('r')])  # How many reactions (edges)
        deprotections = [d for d in esr if d.startswith('d')]
        n_deprotections = len(deprotections)  # How many deprotections
        for i, node in enumerate(nodes):
            # compute the edges that will be present at this node for this sub-synthesis
            cyclic_edges = [last_cycle_edge + 1 for item in range(cyclic_valences[i])]
            # we use the original last_cycle edge because only one cycle is allowed per sub-synthesis
            if i == 0:
                edges = [item for item in nodes[i + 1: i + 1 + linking_valences[i]]]
            else:
                edges = [item for item in nodes[i: i + linking_valences[i]]]
            # create an object for this node
            this_dict = dict()
            this_dict['node'] = node
            this_dict['bbs_file'] = None
            this_dict['preparations'] = []
            if len(d0) > 0 and i == 0:  # We are in a reaction where the first node needs deprotection before reaction
                for di0 in d0:
                    depr_index = int(di0.replace('d', ''))
                    this_preparation = dict()
                    this_preparation['edges'] = []
                    this_preparation['keep_unprepared'] = False
                    this_preparation['reaction'] = enum_deprotection.par[depr_index]['enum_rxn']
                    this_dict['preparations'].append(copy.deepcopy(this_preparation))
            # Complete the information for the node object
            this_dict['edges'] = copy.deepcopy(edges)
            this_dict['cycle_edges'] = copy.deepcopy(cyclic_edges)
            this_dict['bond_orders'] = [mrd['bond_orders'][j] for j in range(len(edges))]
            this_dict['cycle_bond_orders'] = [mrd['cycle_bond_order'] for _ in cyclic_edges]
            if mrd['out_isotopes'] - 1 != i:  # this node does not give an out isotope
                this_dict['out_edges'] = [0 for _ in this_dict['edges']]
            else:
                this_dict['out_edges'] = [nodes[2]]  # out isotopes are allowed only in preparations that prepare a single edge
            # Now we work in the reaction itself
            this_preparation = dict()
            this_preparation['keep_unprepared'] = False
            this_preparation['reaction'] = preparations.par[mrd['preparations'][i]]['preparation']
            this_preparation['edges'] = copy.deepcopy(edges)
            this_preparation['edges'] += this_dict['cycle_edges']  # we append the cycle edges to the edges for the preparation
            this_dict['preparations'].append(copy.deepcopy(this_preparation))
            sub_synthesis.append(copy.deepcopy(this_dict))
        if sum([cyclic_valences[i] for i, node in enumerate(nodes)]) > 0:
            last_cycle_edge += 1  # there was a cycle in the sub_synthesis so we increase the last_cycle_edge
        if used_cycle_edges == 1:
            print(f'ERROR::: found one and only one node with a cycle edge in sub-synthesis {s_reactions} ')
            with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                f.write(f'ERROR::: found one and only one node with a cycle edge in sub-synthesis {s_reactions}\n')
            sys.exit(1)
        if used_cycle_edges == 2:
            last_cycle_edge += 1
        return sub_synthesis, last_cycle_edge

    def merge_nodes(self, node1, node2):
        """merge one node into another to build the config file for enumerations
        node1: dict
        node2: dict
        return: merged_node. Dict"""
        new_node = dict()
        new_node['node'] = node1['node']
        new_node['bbs_file'] = node1['bbs_file']
        new_node['edges'] = node1['edges']
        new_node['out_edges'] = node1['out_edges']
        new_node['bond_orders'] = node1['bond_orders']
        new_node['cycle_edges'] = node1['cycle_edges']
        new_node['cycle_bond_orders'] = node1['cycle_bond_orders']
        new_node['preparations'] = node1['preparations']
        for edge, out_edge, bond_order in zip(node2['edges'], node2['out_edges'], node2['bond_orders']):
            if edge not in new_node['edges']:
                new_node['edges'].append(edge)
                new_node['out_edges'].append(out_edge)
                new_node['bond_orders'].append(bond_order)
        for cycle_edge, cycle_bond_order in zip(node2['cycle_edges'], node2['cycle_bond_orders']):
            if cycle_edge not in new_node['cycle_edges']:
                new_node['cycle_edges'].append(cycle_edge)
                new_node['cycle_bond_orders'].append(cycle_bond_order)
        # the next part of the code removes deprotections that do not keep unprepared to avoid using the same
        # required deptorection in the same node more than once because it was the start of two different sub_synthesis.
        # this could happend in cyanuric acid addition to an amine or glycine ester reductive amination with an aldehyde.
        excluded_preps = []
        for i, preparation2 in enumerate(node2['preparations']):
            excluded_prep = False
            for j, preparation1 in enumerate(new_node['preparations']):
                deprs = len(preparation2['edges']) == 0 and len(preparation1['edges']) == 0
                req_deprs = not preparation2['keep_unprepared'] and not preparation1['keep_unprepared']
                same_preps = preparation2['reaction'] == preparation1['reaction']
                if deprs and req_deprs and same_preps:
                    excluded_prep = True
                    break
            if not excluded_prep:
                new_node['preparations'].append(preparation2)
        return new_node

    def get_graph_config(self, sub_synthesis_list):
        """creates the config file to enumerate the library from the list of sub_synthesis
        sub_synthesis_list
        sub_synthesis_list: list of dictionaries containing raw config file
        returns: config: list of dictionaries containing the config file"""
        config = [None for _ in range(self.n_cycles + 1)]
        for node in sub_synthesis_list:
            if config[node['node']] is None:  # this is the first time that this node appears
                config[node['node']] = copy.deepcopy(node)
            else:  # merge is needed
                if min(node['edges']) < min(config[node['node']]['edges']):  # todo figure out if < or <=
                    # merge existing node into incoming node
                    config[node['node']] = copy.deepcopy(self.merge_nodes(node, config[node['node']]))
                else:
                    # merge incoming node into existing node
                    config[node['node']] = copy.deepcopy(self.merge_nodes(config[node['node']], node))
        return config

    def find_s_reactions_and_nodes(self):
        """This code finds the s_reactions list and nodes for a design given its lib_id.
        returns tuple of c_nodes and c_s_reactions
            c_nodes: list of lists of int (contains the nodes in each sub synthesis)
            c_s_reactions: list of lists of strings (contains the reaction sequence for each sub-synthesis)"""
        def find_source(n):
            """returns source type and cycle 0 based counting headpiece"""
            if self.lib_id[n] % 3 == 0:
                st = 'cycle'
                si = self.lib_id[n] // 3
            elif (self.lib_id[n] - 1) % 3 == 0:
                st = 'deprotection'
                si = (self.lib_id[n] - 1) // 3
            elif (self.lib_id[n] - 2) % 3 == 0:
                st = 'reaction'
                si = (self.lib_id[n] - 2) // 3
            return st, si

        def add_source(cycle, st, si, j, c_nodes, c_s_reactions, lib_id, iteration=0, verbose=False):
            iteration += 1
            if iteration > 20:
                print('ERROR too many iterations in add_source function')
                return c_nodes, c_s_reactions, iteration
            if st == 'cycle':
                c_nodes[j] = [si] + c_nodes[j]
                if (verbose):
                    print(f'cycle {cycle}: found direct final source in cyle {si}')
                return c_nodes, c_s_reactions, iteration
            elif st == 'reaction':
                c_s_reactions[j] = [f'r{lib_id[r_pos[si]]}'] + c_s_reactions[j]
                c_nodes[j] = [si + 1] + c_nodes[j]  # add upstream node for this reaction
                if verbose:
                    print(f'cycle {cycle}: found source in reaction {si + 1}')
                st, si = find_source(rt_pos[si])
                c_nodes, c_s_reactions, iteration = add_source(cycle, st, si, j, c_nodes, c_s_reactions, lib_id,
                                                               iteration=iteration, verbose=verbose)
                return c_nodes, c_s_reactions, iteration
            elif st == 'deprotection':
                c_s_reactions[j] = [f'd{lib_id[d_pos[si]]}'] + c_s_reactions[j]
                # dont add node here because it will be added downstream
                if verbose:
                    print(f'cycle {cycle}: found source in deprotection {si + 1}')
                st, si = find_source(dt_pos[si])
                c_nodes, c_s_reactions, iteration = add_source(cycle, st, si, j, c_nodes, c_s_reactions, lib_id,
                                                               iteration=iteration, verbose=verbose)
                return c_nodes, c_s_reactions, iteration
            else:
                print(f"ERROR::: no clear source found for cycle {cycle} in source index {si}")
                with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                    f.write(f"ERROR::: no clear source found for cycle {cycle} in source index {si}\n")
                return None, None, iteration

        rt_pos = [1 + 3 * self.n_cycles + i for i in range(self.n_cycles)]
        dt_pos = [1 + 2 * self.n_cycles + i for i in range(self.n_cycles)]
        r_pos = [1 + self.n_cycles + i for i in range(self.n_cycles)]
        d_pos = [1 + i for i in range(self.n_cycles)]
        c_nodes = [[] for _ in range(self.n_cycles)]
        c_s_reactions = [[] for _ in range(self.n_cycles)]
        for i in range(self.n_cycles):
            cycle = self.n_cycles - i  # go backwards
            j = self.n_cycles - i - 1  # backward index
            if i > 0 and any([cycle in item[1:] for item in c_nodes[j + 1:]]):
                continue  # this allows not to track the same sub-route twice
            c_nodes[j] = [cycle] + c_nodes[j]  # add the cycle to the cycles
            c_s_reactions[j] = [f'r{self.lib_id[r_pos[j]]}'] + c_s_reactions[j]  # add the cycle to the cycles
            st, si = find_source(rt_pos[j])
            c_nodes, c_s_reactions, iteration = add_source(cycle, st, si, j, c_nodes, c_s_reactions, self.lib_id,
                                                           verbose=self.verbose)
            if self.verbose:
                print(f'cycle {cycle}: found all sources in {iteration} iterations')

        c_nodes = [item for item in c_nodes if len(item) > 0]
        c_s_reactions = [item for item in c_s_reactions if len(item) > 0]
        return c_nodes, c_s_reactions

    def get_sub_synthesis_list(self, multireaction, preparations, enum_deprotection):
        """gets a list of sub_synthesis for enumeration purposes with graph_enumerator.
        multireaction: instance of Par class coding the multireaction parameters
        preparations: instance of Par class coding the preparations parameters
        enum_deprotection: instacne of Par class coding the enumeration deprotections
        returns: sub_synthesis_list, list of dictionaries. One dictionary per node per reaction
        """
        sub_synthesis_list = []  # list of sub_synthesis objects, one per cycle
        last_cycle_edge = self.lib_id[0]
        c_nodes, c_s_reactions = self.find_s_reactions_and_nodes()
        with open(os.path.join(self.wfolder, 's_reactions.txt'), 'w') as f:
            for item in c_s_reactions:
                f.write(' '.join(item) + '\n')
        if c_nodes is None or c_s_reactions is None:
            sys.exit(1)
        for i, (nodes, s_reactions) in enumerate(zip(c_nodes, c_s_reactions)):
            if len(s_reactions) == 0:
                print(f'ERROR::: could not get s_reaction for node {self.n_cycles - i}')
                with open(os.path.join(self.wfolder, 'error.txt'), 'w') as f:
                    f.write(f'ERROR::: could not get s_reaction for node {self.n_cycles - i}\n')
                sys.exit(1)
            if self.verbose:
                print('s_reactions, nodes', s_reactions, nodes)
            sub_synthesis, last_cycle_edge = self.get_sub_synthesis(s_reactions, nodes, last_cycle_edge,
                                                                    multireaction, preparations, enum_deprotection)
            sub_synthesis_list += copy.deepcopy(sub_synthesis)
        return sub_synthesis_list

    def write_graph_enumeration_json(self, multireaction, preparations, enum_deprotection):
        """Creates the json config file for the graph enumeration for this specific lib design
        multireaction: instance of Par class coding the multireaction parameters
        preparations: instance of Par class coding the preparations parameters
        enum_deprotection: instacne of Par class coding the enumeration deprotections
        headpieces_dict: dictionary of smiles strings for the headpieces
        returns: str: path to the enumerated library"""
        # Create the draft config file
        sub_synthesis_list = self.get_sub_synthesis_list(multireaction, preparations, enum_deprotection)
        config = self.get_graph_config(sub_synthesis_list)
        # add standard deprotections to each cycle
        std_boc_deprotection = {"edges": [],
                                "keep_unprepared": True,
                                "reaction": "${DEPROTECTION_FOLDER}/6.1.1_N-Boc_deprotection_FROM_boc_AND_Null.rxn"}
        std_ester_hydrolysis = {"edges": [],
                                "keep_unprepared": True,
                                "reaction": "${DEPROTECTION_FOLDER}/9.7.61_Ester_hydrolysis_FROM_esters-methyl-ethyl_AND_Null.rxn"}
        std_fmoc_deprotection = {"edges": [],
                                 "keep_unprepared": True,
                                "reaction": "${DEPROTECTION_FOLDER}/6.1.6_N-Fmoc_deprotection_FROM_fmoc_AND_Null.rxn"}
        for i, node in enumerate(config):
            config[i]["preparations"].append(copy.deepcopy(std_boc_deprotection))
            #config[i]["preparations"].append(copy.deepcopy(std_boc_deprotection))
            config[i]["preparations"].append(copy.deepcopy(std_ester_hydrolysis))
            #config[i]["preparations"].append(copy.deepcopy(std_ester_hydrolysis))
            config[i]["preparations"].append(copy.deepcopy(std_fmoc_deprotection))
            #config[i]["preparations"].append(copy.deepcopy(std_fmoc_deprotection))
        # Add the appropriate bbs files to the config file
        for i, node in enumerate(config):
            out_bbs_file = f'C{i}.smi'  # todo assumes that the headpiece is stored in C0.smi and len(config) = 1 + len(self.n_cycles)
            config[i]["bbs_file"] = os.path.abspath(os.path.join(self.wfolder, out_bbs_file))

        # write config to file
        with open(os.path.join(self.wfolder, "config.json"), 'w') as f:
            json.dump(config, f, indent=2)

    def write_summary_file(self, enum_reaction, enum_deprotection):
        """writes a file with a summary of the library desing to be understood by a human
        enum_reaction: instance of Par class coding the enum_reaction parameters
        enum_deprotection: instance of Par class coding the enum_deprotection parameters
        returns: None"""
        deprotections = [self.lib_id[i] for i in range(1, self.n_cycles + 1)]
        reactions = [self.lib_id[i] for i in range(1 + self.n_cycles, 2 * self.n_cycles + 1)]
        with open(os.path.join(self.wfolder, 'libDESIGN_summary.txt'), 'w') as f:
            f.write(f'LIBRARY SUMMARY:\n')
            f.write(f'Library_id: {self.lib_id}\n')
            f.write(f'Working_folder: {self.wfolder}\n')
            f.write(f'Headpiece_smiles: {self.hp_smiles}\n')
            for i, (deprotection, reaction) in enumerate(zip(deprotections, reactions)):
                if deprotection > 0:
                    f.write(f'Add Deprotection at cycle {i}: {enum_deprotection.par[deprotection]["enum_name"]}\n')
                f.write(f'Add cycle {i+1} through reaction: {enum_reaction.par[reaction]["enum_name"]}\n')
            f.write('\n')

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
        print(f'Working_folder: {self.wfolder}')
        print(f'Headpiece_smiles: {self.hp_smiles}')
        for i, (deprotection, reaction) in enumerate(zip(deprotections, reactions)):
            if deprotection > 0:
                print(f'Add Deprotection at cycle {i}: {enum_deprotection.par[deprotection]["enum_name"]}')
            print(f'Add cycle {i+1} through reaction: {enum_reaction.par[reaction]["enum_name"]}')
        print('')

    def run_bb_analysis(self):
        """runs an analysis of the building blocks and writes a set of files with the compounds in the original
        bbs files that were discarded during preparation
        returns: None"""
        for i in range(self.n_cycles):
            cycle = i + 1
            idf = pd.read_csv(os.path.join(self.wfolder, f'C{cycle}.smi'), sep=' ', header=None, names=['smiles', 'id'])
            odf = pd.read_csv(os.path.join(self.wfolder, f'R0{cycle}.smi'), sep=' ', header=None, names=['smiles', 'id'])
            odf = pd.merge(idf, odf, on='id', how='inner')
            idf = idf[~ idf['id'].isin(odf['id'].unique().tolist())].copy()
            if idf.shape[0] > 0:
                idf[['smiles', 'id']].to_csv(os.path.join(self.wfolder, f'C{cycle}_badBBs.smi'),
                                             header=False, index=False, sep=' ')


    def run_graph_enumeration(self, multireaction, preparations, enum_deprotection, enum_reaction,
                              n=0, chunksize=400000, just_json=False):
        """Runs an enumeration for the instance of the class through graph_enumerator
        base_foldr: str: path to the base folder
        multireaction: instance of Par class coding the multireaction parameters
        preparations: instance of Par class coding the preparations parameters
        enum_deprotection: instacne of Par class coding the enumeration deprotections
        n: str: number of compounds to enumerate, 0 means the whole library
        chunksize: int: maximum number of molecules to enumerate in each core
        returns: str: path to the enumerated library"""
        self.write_graph_enumeration_json(multireaction, preparations, enum_deprotection)
        self.write_summary_file(enum_reaction, enum_deprotection)
        if not just_json:
            # write the source building block files
            self.write_bbs_files()
            # instantiate the graph enumerator
            gen = SynthGraph(self.wfolder, os.path.join(self.wfolder, "config.json"), n=n, chunksize=chunksize)
            gen.run_graph()
            if not gen.success:
                pass
            print(f'{os.path.join(self.wfolder, "enumeration.smi")} has been created')
            command = f'wc -l {self.wfolder}/*.smi'
            os.system(command)
            self.run_bb_analysis()


if __name__ == '__main__':
    print(version)
