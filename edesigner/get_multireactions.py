# -*- coding: utf-8 -*-
# get_multireactions
# Jose Alfredo Martin


import os
import sys
from classes.libdesign import LibDesign
from classes.enumerator import Enumerator
from tqdm import tqdm
import _pickle as pic
from classes.parameter_reader import Parameters
import argparse
import copy

__author__ = 'Alfredo Martin 2023'
__version__ = 'get_multireactions.v.12.0.0'



class MultiReactionFinder:

    def __init__(self, design, verbose=False):
        self.design = design
        self.verbose = verbose

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
        linking_valences = [item['n_edges'] for item in valences]
        cyclic_valences = [0 for _ in range(len(valences))]
        active_valences = 0
        hist_active_valences = []
        for i, item in enumerate(valences):
            active_valences = abs(active_valences - item['n_edges'])
            if item['out_isotope']:
                active_valences += 1
            hist_active_valences.append(active_valences)
        for i, active_valence in enumerate(hist_active_valences):
            if i < len(hist_active_valences) - 1 and active_valence == 0:
                print("ERROR::: No valence available to continue")
                sys.exit(1)
            if i == len(hist_active_valences) - 1 and active_valence > 0:
                print("ERROR::: Subsynthesis ends with an open valence")
                sys.exit(1)
            if i < len(hist_active_valences) - 1:
                if abs(active_valence - hist_active_valences[i + 1]) > 1:
                    cyclic_valences[i] += abs(active_valence - hist_active_valences[i + 1]) - 1
                    cyclic_valences[i + 1] += abs(active_valence - hist_active_valences[i + 1]) - 1
        if sum(cyclic_valences) > 2:
            print("ERROR::: Only one cycle edge is allowed per sub-synthesis")  # todo is this really required?
            sys.exit(1)
        linking_valences = [item[0] - item[1] for item in zip(linking_valences, cyclic_valences)]
        return linking_valences, cyclic_valences

    def find_s_reactions_and_nodes(self):
        """This code finds the s_reactions list and nodes for a design given its lib_id.
        returns tuple of c_nodes and c_s_reactions
            c_nodes: list of lists of int (contains the nodes in each sub synthesis)
            c_s_reactions: list of lists of strings (contains the reaction sequence for each sub-synthesis)"""
        def find_source(n):
            """returns source type and cycle 0 based counting headpiece"""
            if self.design.lib_id[n] % 3 == 0:
                st = 'cycle'
                si = self.design.lib_id[n] // 3
            elif (self.design.lib_id[n] - 1) % 3 == 0:
                st = 'deprotection'
                si = (self.design.lib_id[n] - 1) // 3
            elif (self.design.lib_id[n] - 2) % 3 == 0:
                st = 'reaction'
                si = (self.design.lib_id[n] - 2) // 3
            return st, si

        def add_source(cycle, st, si, j, c_nodes, c_s_reactions, lib_id, iteration=0, verbose=False):
            iteration += 1
            if iteration > 20:
                print('ERROR too many iterations in add_source function')
                return c_nodes, c_s_reactions, iteration
            if st == 'cycle':
                c_nodes[j] = [si] + c_nodes[j]
                return c_nodes, c_s_reactions, iteration
            elif st == 'reaction':
                c_s_reactions[j] = [f'r{lib_id[r_pos[si]]}'] + c_s_reactions[j]
                c_nodes[j] = [si + 1] + c_nodes[j]  # add upstream node for this reaction
                st, si = find_source(rt_pos[si])
                c_nodes, c_s_reactions, iteration = add_source(cycle, st, si, j, c_nodes, c_s_reactions, lib_id,
                                                               iteration=iteration, verbose=verbose)
                return c_nodes, c_s_reactions, iteration
            elif st == 'deprotection':
                c_s_reactions[j] = [f'd{lib_id[d_pos[si]]}'] + c_s_reactions[j]
                # dont add node here because it will be added downstream
                st, si = find_source(dt_pos[si])
                c_nodes, c_s_reactions, iteration = add_source(cycle, st, si, j, c_nodes, c_s_reactions, lib_id,
                                                               iteration=iteration, verbose=verbose)
                return c_nodes, c_s_reactions, iteration
            else:
                print(f"ERROR::: no clear source found for cycle {cycle} in source index {si}")
                return None, None, iteration

        rt_pos = [1 + 3 * self.design.n_cycles + i for i in range(self.design.n_cycles)]
        dt_pos = [1 + 2 * self.design.n_cycles + i for i in range(self.design.n_cycles)]
        r_pos = [1 + self.design.n_cycles + i for i in range(self.design.n_cycles)]
        d_pos = [1 + i for i in range(self.design.n_cycles)]
        c_nodes = [[] for _ in range(self.design.n_cycles)]
        c_s_reactions = [[] for _ in range(self.design.n_cycles)]
        for i in range(self.design.n_cycles):
            cycle = self.design.n_cycles - i  # go backwards
            j = self.design.n_cycles - i - 1  # backward index
            if i > 0 and any([cycle in item[1:] for item in c_nodes[j + 1:]]):
                continue  # this allows not to track the same sub-route twice
            c_nodes[j] = [cycle] + c_nodes[j]  # add the cycle to the cycles
            c_s_reactions[j] = [f'r{self.design.lib_id[r_pos[j]]}'] + c_s_reactions[j]  # add the cycle to the cycles
            st, si = find_source(rt_pos[j])
            c_nodes, c_s_reactions, iteration = add_source(cycle, st, si, j, c_nodes, c_s_reactions, self.design.lib_id,
                                                           verbose=self.verbose)
        c_nodes = [item for item in c_nodes if len(item) > 0]
        c_s_reactions = [item for item in c_s_reactions if len(item) > 0]
        return c_nodes, c_s_reactions


def parse_args():
    # Arg parser
    parser = argparse.ArgumentParser(description="""finds all the multireacion sequences generated by an eDESIGNER run 
    and reports the ones not coded in the multireaction parameters. If a multireaction is generated but not coded in 
    the multireaction parameters, the library containing that multireaction cannot be enumerated, therefore these 
    multireactions should be added to the parameters file (see documentation for guidance).""")
    parser.add_argument('-wF', '--wfolder',
                        help="""Working Folder""",
                        type=str,
                        required=True)
    parser.add_argument('-run', '--run_id',
                        help="""name of the run id""",
                        type=str,
                        required=True)
    parser.add_argument('-erun', '--ed_run_id',
                        help="""name of the eDESIGNER run""",
                        type=str,
                        default=None)
    parser.add_argument('-v', '--verbose',
                        help="""verbose output""",
                        action='store_true')

    args = parser.parse_args()
    assert os.path.isdir(args.wfolder), f'{args.wfolder} does not exist'
    args.wfolder = os.path.abspath(args.wfolder)
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id)), f'{os.path.join(args.wfolder, args.run_id)} does not exist'
    assert os.path.isdir(os.path.join(args.wfolder, args.run_id, args.ed_run_id)), f'{os.path.join(args.wfolder, args.run_id, args.ed_nun_id)} does not exist'
    return args

def main():
    args = parse_args()
    PARFOLDER = os.path.join(args.wfolder, args.run_id, 'resources')
    bbts_file = os.path.join(args.wfolder, args.run_id, 'results', 'BBTs.pic')
    libdesigns_file = os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'libDESIGNs.pic')
    multireaction = Parameters(os.path.join(PARFOLDER, 'multireaction.par'),
                               fsource='list', how='to_list', multiple=True)
    preparations = Parameters(os.path.join(PARFOLDER, 'preparations.par'),
                              fsource='list', how='to_list', multiple=True)
    enum_deprotection = Parameters(os.path.join(PARFOLDER, 'enum_deprotection.par'),
                                   fsource='list', how='to_list', multiple=True)
    enum_reaction = Parameters(os.path.join(PARFOLDER, 'enum_reaction.par'),
                               fsource='list', how='to_list',multiple=True)
    headpieces = Parameters(os.path.join(PARFOLDER, 'headpieces.par'),
                            fsource='list', how='to_list', multiple=True)
    with open(bbts_file, 'rb') as f:
        BBTs = pic.load(f)
    headpieces_dict = {}
    for headpiece in headpieces.par:
        headpieces_dict[[BBT.BBT for BBT in BBTs].index(headpiece['bbt'])] = headpiece['smiles']
    designs = []
    with open(libdesigns_file, 'rb') as f:
        while True:
            try:
                design = pic.load(f)
                designs.append(design)
            except:
                break
    sub_synthesis_list = []
    design_idx_list = []
    for design in tqdm(designs):
        finder = MultiReactionFinder(design, verbose=args.verbose)
        _, this_sub_synthesis_list = finder.find_s_reactions_and_nodes()

        for item in this_sub_synthesis_list:
            while True:
                if len(item) == 0:
                    print(f'ERROR::: found a multireaction containing only deprotections {this_sub_synthesis_list}')
                    break
                if item[0].startswith('d'):
                    item.pop(0)
                else:
                    if tuple(item) not in sub_synthesis_list:
                        sub_synthesis_list.append(tuple(item))
                        design_idx_list.append(design.id)
                    break
    if args.verbose:
        print()
        print('Not_found_multireacton : design_idx (Code these before enumerating the associaed designs)')
    not_found = []
    mrd_list = []
    for mrd in multireaction.par:  # mrd is the multi-reaction dictionary that matches the esr
        mrd_list.append(tuple(mrd['s_reactions']))
    for item, design_idx in zip(sub_synthesis_list, design_idx_list):
        if item not in mrd_list:
            not_found.append(item)
            if args.verbose:
                print(';'.join(list(item)), ':', design_idx)
    if args.verbose and len(not_found) == 0:
        print('  All multireaction in this eDESIGNER run were found in the parameters. You are OK to enumerate')

    out_file = os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'full_multiereaction_report.csv')
    with open(out_file, 'w') as f:
        f.write('multireaction,design_idx,found')
        for item, design_idx in zip(sub_synthesis_list, design_idx_list):
            if item in not_found:
                f.write(';'.join(list(item)) +',' + str(design_idx) + ',' + 'FALSE\n')
            else:
                f.write(';'.join(list(item)) + ',' + str(design_idx) + ',' + 'TRUE\n')
    print(f' multireaction_report: {out_file}')

    if args.verbose:
        print()
        print('Recommended test set of desings indexes for enumeration:')
    test_idxs = list(set(design_idx_list))
    for item in test_idxs:
        if args.verbose:
            print(item)

    out_file = os.path.join(args.wfolder, args.run_id, args.ed_run_id, 'results', 'enumeration_test_idxs.txt')
    with open(out_file, 'w') as f:
        for item in test_idxs:
            f.write(str(item) + '\n')
    print(f' enumeration test idxs: {out_file}')


if __name__ == '__main__':
    print(__version__)
    print(__author__)
    main()
