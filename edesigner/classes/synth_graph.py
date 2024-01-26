# -*- coding: utf-8 -*-
# graph_enumerator
# Jose Alfredo Martin 2023

__version__ = 'graph_enumerator_for_eDESIGNER.v.1.0.0'
__author__ = 'Alfredo Martin 2023'

# Python modules
import copy
import os
import sys
import json
import shutil

class SynthNode:
    """class that holds a node in the graph containing a set of builing blocks"""

    def __init__(self, wfolder, node):
        """Constructor of the instance
        wfolder: str: path to the working folder to store working files
        node: dict: contains the following keys:
            node: int: index of the node
            bbs_file: str: path to the bbs file
            edges: list of int: edges for this node:
            cycle_edges: list of int: cyclic edges for this node:
            preparations: list of dictionaries: (empty list is ok). Each dict contain the following keys:
                edges: list of int: edges arising from this preparation
                keep_unprepared: bool: whether to keep the bbs that fail preparation or deprotection
                reaction: str: path to the reaction file
        """
        self.node = node['node']
        self.edges = [item for item in node['edges']]
        self.bond_orders = [item for item in node['bond_orders']]
        self.cycle_edges = [item for item in node['cycle_edges']]
        self.cycle_bond_orders = [item for item in node['cycle_bond_orders']]
        self.preps = copy.deepcopy(node['preparations'])
        self.wfolder = wfolder
        filename = os.path.join(self.wfolder, "R" + str(self.node).rjust(2, "0")) + '.smi'
        command = f'cp {node["bbs_file"]} {filename}'  # makes a copy of the input file to avoid loss by corruption
        os.system(command)
        command = f"sed -i '/^$/d' {filename}"
        os.system(command)
        self.bbsfile = filename
        if os.path.isdir(os.environ['LILLYMOL_EXECUTABLES']):
            self.trxn = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'trxn')
            self.make_these_molecules = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'make_these_molecules')
        else:
            self.trxn = 'trxn.sh'
            self.make_these_molecules = 'make_these_molecules.sh'
        self.success = True

    def countlines(self, file):  # todo this is new
        """counts the number of lines in a file
        file: str: path to a file
        returns: int"""
        count = 0
        with open(file, 'r') as f:
            for line in f:
                if line.strip():
                    count += 1
        return count

    def add_preparations(self):
        """ Adds all the deprotections and preparations for this node
        """
        for i, prep in enumerate(self.preps):
            if prep['reaction'].endswith('.rxn'):
                preparation = os.path.join(self.wfolder, 'preparation.rxn')
                preptype = '-r'
            elif prep['reaction'].endswith('.prxn'):
                preparation = os.path.join(self.wfolder, 'preparation.proto')
                preptype = '-P'
            else:
                preptype = None
                print(f"ERROR::: {prep['reaction']} is not of a valid type to be used as a preparation")
                self.success = False
                sys.exit(1)
            command = f'cp {prep["reaction"]} {preparation}'
            os.system(command)
            command = "sed -i s%'${EDESIGNER_TEST_FOLDER}'%" + os.environ['EDESIGNER_TEST_FOLDER'] + "%g " + preparation
            os.system(command)
            for j, edge in enumerate(prep['edges']):
                command = f"sed -i s/'$isotope_{j+1}'/{edge}/ {preparation}"
                os.system(command)
            stem = os.path.join(self.wfolder, "preparation")
            outfile = stem + '.smi'
            if os.path.isfile(outfile):
                os.remove(outfile)
            if prep['keep_unprepared']:
                tag = 'w'
            else:
                tag = 'i'
            if len(prep['edges']) == 0:
                mtag = '-m do=0 -M do=0'
            else:
                mtag = '-m RMX -M RMX'

            ini_lines = self.countlines(self.bbsfile)
            if not prep['keep_unprepared']:
                print("****")
                print('preparation', prep["reaction"])
                print('bbsfile', self.bbsfile)
            command = f'{self.trxn}'
            command += f' {preptype} {preparation} -z {tag} {mtag} -W "" -i smi -o smi'
            command += f' -S {stem}'
            command += f' {self.bbsfile}'
            os.system(command)
            end_lines = self.countlines(outfile)
            if end_lines == 0:
                print('ERROR::: The output file is empty')
                self.success = False
                print('command:', command)
                sys.exit(1)
            if not prep['keep_unprepared']:
                print(f'Keeping unprepared: {prep["keep_unprepared"]}. Processed {ini_lines} compounds. Obtained {end_lines} compounds.')
                if ini_lines != end_lines:
                    print(f'WARNING::: There were {ini_lines - end_lines} compounds lost in this reaction')
                    print('command:', command)
                print("*****")
            command = f'mv {outfile} {self.bbsfile}'
            os.system(command)
            if self.node == 0:
                shutil.copy(self.bbsfile, os.path.join(self.wfolder, "RP00.smi"))


class SynthGraph:
    """class that specifies the graph based enumeration"""

    def __init__(self, wfolder, config_file, n=0, chunksize=400000):
        """constructor of the instance
        wfolder: str: path to a folder where files are stored
        config_file: str: path to the json config file
        n: int: number of compounds to enumerate from the graph. A negative number or 0 indicates enumerate
            all of them
        chunksize: int: number of compounds enumerated in each core using hpc"""
        self.nodes = []  # list of node instances
        self.edges = []  # list of edge dictionaries (edge_id, orig_node, dest_node, orig_isotope, dest_isotope and
        # bond_order are the keys)
        # the goal of the self.edges object is to allow a single atom to form more than one different bond to
        # different nodes, for example reductive amination of aldehyde with primary amine followed by acylation.
        # this only affects to the main edges and not the cyclization edges. Since these are the last ones to
        # compute the isotopes should be set. Double macro-cyclizations on the same atom are not supported.
        self.products = None
        self.wfolder = wfolder
        self.n = n
        self.chunksize = chunksize
        with open(config_file, 'r') as f:
            self.par = json.load(f)
        self.par_quality_control()
        if os.path.isdir(os.environ['LILLYMOL_EXECUTABLES']):
            self.trxn = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'trxn')
            self.make_these_molecules = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'make_these_molecules')
        else:
            self.trxn = 'trxn.sh'
            self.make_these_molecules = 'make_these_molecules.sh'
        self.success = True

    def add_node(self, node):
        """adds a node to the graph. Nodes must be added in a sequential fashion with first node being 0. edges must
        contain the index given to this node except for node 0, if node contains cycle_edges, all their indexes must be
        greater than the index of the last node. cycle_edge indexes are independent of nodes and must be unique per
        formed bond
        node: dict: contains the following keys:
            node: int: index of the node
            bbs_file: str: path to the bbs file
            edges: list of int: edges for this node
            bond_orders: list of int: order of the bonds formed
            out_edges: list of int: resulting isotope after attaching through each edge
            cycle_edges: list of int: cyclic edges for this node:
            cyclic_bond_orders: list of int: bond orders of the cyclic bonds formed
            preparations: list of dictionaries: (empty list is ok). Each dict contain the followitg keys:
                edges: list of int: edges arising from this preparation
                keep_unprepared: bool: whether to keep the bbs that fail preparation or deprotection
                reaction: str: path to the reaction file

        """
        self.nodes.append(SynthNode(self.wfolder, node))
        if len(self.nodes) == 1:
            self.products = copy.deepcopy(self.nodes[0])
        # update the self.edges object
        # values: orig_node: original node. it is the node with highest index for this edge
        #         dest_node: destination node. it is the node with lowest index for this edge
        #         bond_order: is the bond order for the edge
        #         orig_isotope: is the out isotope for the atom of dest_node for the orig_node
        #         dest_isotope: is the out isotope for the atom of the dest_node
        for bond_order, edge, out_edge in zip(node['bond_orders'], node['edges'], node['out_edges']):
            if edge not in [item['edge_id'] for item in self.edges]:
                # This is the first time we work on this edge but we know that
                # the origin is not the current node (otherwise it would be not new).
                # So therefore teh destination is this node.
                self.edges.append({'edge_id': edge,
                                   'values': {'orig_node': edge,
                                              'dest_node': node['node'],
                                              'bond_order': bond_order,
                                              'orig_isotope': 0,  # placeholder
                                              'dest_isotope': out_edge}
                                   })
            else:
                # This is an edge that is seen before and corresponds to this specific node
                # so we have to update the orig_isotope
                for i, item in enumerate(self.edges):
                    if item['edge_id'] == edge:
                        self.edges[i]['values']['orig_isotope'] = out_edge
                        break
            self.edges.sort(key=lambda x: x['edge_id'])  # sorts the edges so its index is its id-1

    def create_reaction(self, node):
        """creates a trxn reaction for this speficic edge
        node: int: node index
        returns: str: path to the reaction file"""
        reaction_file = 'join_by_isotope_' + str(node).rjust(2, "0") + '.rxn'
        edge = self.edges[node-1]  # node 0 is not the origin or any edges
        with open(os.path.join(self.wfolder, reaction_file), 'w') as f:
            linea = '(0 Reaction\n'
            f.write(linea)
            linea = f'  (A C Comment "join_by_isotope_{str(node).rjust(2, "0")}")\n'
            f.write(linea)
            linea = '  (0 Scaffold\n'
            f.write(linea)
            linea = f'    (A  C smarts "[{node}*]")\n'
            f.write(linea)
            linea = f'    (A I isotope (0 {edge["values"]["dest_isotope"]}))\n'  # the scaffold is the destination because edge direction is outer to inner
            f.write(linea)
            linea = '  )\n'
            f.write(linea)
            linea = '  (1 Sidechain\n'
            f.write(linea)
            linea = f'    (A C smarts "[{node}*]")\n'
            f.write(linea)
            linea = f'    (A I isotope (0 {edge["values"]["orig_isotope"]}))\n'  # the side chain is the origin because edge direction is outer to inner
            f.write(linea)
            linea = f'    (A I join (0 0 {edge["values"]["bond_order"]}))\n'
            f.write(linea)
            linea = '  )\n'
            f.write(linea)
            linea = ')\n'
            f.write(linea)
        return os.path.join(self.wfolder, reaction_file)

    def create_cyclization(self, edge, cycle_bond_order):
        """creates a trxn reaction for cyclization at this speficic edge
        edge: int: edge number
        cycle_bond_order: int: order of the bond to be created
        returns: str: path to the reaction file"""
        reaction_file = 'cyclize_by_isotope_' + str(edge).rjust(2, "0") + '.rxn'
        with open(os.path.join(self.wfolder, reaction_file), 'w') as f:
            linea = '(0 Reaction\n'
            f.write(linea)
            linea = f'  (A C Comment "cyclize_by_isotopes_{str(edge).rjust(2, "0")}")\n'
            f.write(linea)
            linea = '  (0 Scaffold\n'
            f.write(linea)
            linea = f'    (A  C smarts "[{edge}*]...[{edge}*]")\n'
            f.write(linea)
            linea = '    (A I one_embedding_per_start_atom 1)\n'
            f.write(linea)
            linea = '    (A I isotope (0 0))\n'
            f.write(linea)
            linea = '    (A I isotope (1 0))\n'
            f.write(linea)
            if cycle_bond_order == 1:
                linea = '    (A I single_bond (0 1))\n'
            elif cycle_bond_order == 2:
                linea = '    (A I double_bond (0 1))\n'
            elif cycle_bond_order == 3:
                linea = '    (A I triple_bond (0 1))\n'
            else:
                print(f'ERROR::: in edge {edge}, bond order {cycle_bond_order} is not supported')
                self.success = False
                sys.exit(1)
            f.write(linea)
            linea = '  )\n'
            f.write(linea)
            linea = ')\n'
            f.write(linea)
        return os.path.join(self.wfolder, reaction_file)

    def create_remove_isotope(self, edge):
        """creates a trxn reaction that removes the isotope for this speficic edge
        edge: int: edge number
        returns: str: path to the reaction file"""
        reaction_file = 'remove_isotope_' + str(edge).rjust(2, "0") + '.rxn'
        with open(os.path.join(self.wfolder, reaction_file), 'w') as f:
            linea = '(0 Reaction\n'
            f.write(linea)
            linea = f'  (A C Comment "remove_isotope_{str(edge).rjust(2, "0")}")\n'
            f.write(linea)
            linea = '  (0 Scaffold\n'
            f.write(linea)
            linea = f'    (A  C smarts "[{edge}*]")\n'
            f.write(linea)
            linea = '    (A I isotope (0 0))\n'
            f.write(linea)
            linea = '  )\n'
            f.write(linea)
            linea = ')\n'
            f.write(linea)
        return os.path.join(self.wfolder, reaction_file)

    def count_compounds(self, file1, file2=None, file3=None):
        """counts the number of building blocks in one to three.
        file1: str: path for file 1
        file2: str or None: path for file 2
        returns: int or tuple of two ints"""
        wc_file = os.path.join(self.wfolder, 'wc.txt')
        command = f'wc {file1} '
        if file2 is not None:
            command += f'{file2} '
        if file3 is not None:
            command += f'{file3} '
        command += f'> {wc_file}'
        os.system(command)
        with open(wc_file, 'r') as f:
            lineas = f.readlines()
        os.remove(wc_file)
        wc1 = str(lineas[0].strip().split(' ')[0])
        if file3 is not None and file2 is not None:
            wc2 = str(lineas[1].strip().split(' ')[0])
            wc3 = str(lineas[2].strip().split(' ')[0])
            return int(wc1), int(wc2), int(wc3)
        elif file2 is not None:
            wc2 = str(lineas[1].strip().split(' ')[0])
            return int(wc1), int(wc2)
        else:
            return int(wc1)

    def create_sge_script(self):
        """creates a sge script in the working folder
        return: str: path to the sge_script"""
        sge_script = os.path.join(self.wfolder, 'sge_script.sh')
        with open(sge_script, 'w') as f:
            f.write('#$ -S /bin/bash\n')
            f.write('#$ -o arrayjob_script.out\n')
            f.write('#$ -e arrayjob_script.err\n')
            f.write('#$ -j y\n')
            f.write('#$ -cwd\n')
            f.write('source /usr/share/Modules/init/bash\n')  # is this required?
            f.write('module load gc3tk\n')  # is this required?
            f.write('eval $( head -${SGE_TASK_ID} ${1} | tail -1)\n')
        command = f'chmod u+x {sge_script}'
        os.system(command)
        return sge_script

    def create_qsub_script(self, numjobs, sge_script):
        """creates the qsub_script in the working folder
        numjobs: int: number of jobs to be performed
        sge_script: str: path to the sge_script file
        returns: str: path to the qsub script"""
        qsub_script = os.path.join(self.wfolder, 'qsub_script.sh')
        with open(qsub_script, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('source /etc/profile.d/z00_lmod.sh\n')  # is this required?
            f.write('current=$(dirname "$0")\n')
            f.write(f'NUMJOBS={numjobs}\n')
            f.write(f'qsub -V -t 1-{numjobs}:1 -sync y {sge_script} $1\n')
        command = f'chmod u+x {qsub_script}'
        os.system(command)
        return qsub_script

    def enumerate_target_file(self, target_file):
        """enumerates a predefined set of compounds given in a target file. The target file is a text file
        with columns separated by spaces and no header containing ids of the building blocks to add. The columns are
        not headed.
        target_file: str: path to a target file"""
        # todo if target file has more lines than self.chunksize then code the enumeration through qsub
        with open(target_file, 'r') as f:
            linea = f.readline()
        if len(linea.strip().split(' ')) != len(self.nodes):
            print('The number of building blocks in the target_file does not match the number of nodes. Exiting.')
            sys.exit(1)
        reactions = [self.create_reaction(i+1) for i, node in enumerate(self.nodes[1:])]
        bb_files = " ".join([os.path.join(self.wfolder, "R" + str(i).rjust(2, "0") + ".smi") for i in range(len(self.nodes))])
        command = f'{self.make_these_molecules}'
        command = f' -R {" -R ".join([reaction for reaction in reactions])}'
        command += f' -M {target_file}'
        command += f' -S {os.path.join(self.wfolder, "products")}'
        command += f' -z f -z i -W "+" -l -i smi -g all -A D'
        command += f' {bb_files}'
        print(command)
        os.system(command)

        command = f'mv {os.path.join(self.wfolder, "products.smi")} {os.path.join(self.wfolder, "R00.smi")}'
        os.system(command)

    def join_node(self, node):
        """conducts enumeration using the edge joining a node and products
        node: int: index for the node to join
        """
        print('')
        print(f'INFO::: Joining node {node}')
        # set file names
        scaffolds_file = os.path.join(self.wfolder, 'R00.smi')
        reactive_file = os.path.join(self.wfolder, 'R00reactive.smi')
        inert_file = os.path.join(self.wfolder, 'R00inert.smi')
        intermediate_file = os.path.join(self.wfolder, 'R00int.smi')
        products_file = os.path.join(self.wfolder, 'PRD.smi')
        stem = os.path.join(self.wfolder, 'PRD')
        # create reactive and inert files and also create reaction file
        if node < 10:
            pattern = f'"\[{node}[a-z,A-Z]"'
        else:
            pattern = f'"\[{node}"'
        command = f'grep {pattern} {scaffolds_file} > {reactive_file}'
        os.system(command)
        command = f'grep -v {pattern} {scaffolds_file} > {inert_file}'
        os.system(command)
        reaction = self.create_reaction(node)
        # compute number of compounds to enumerate
        n1, n2, n3 = self.count_compounds(inert_file, reactive_file, self.nodes[node].bbsfile)
        n2t = n2
        # reduce input scaffold file if required
        reduce_output = False
        if self.n > 0:
            n2t = min(n2, int(self.n / n3))  # how many molecules should I keep in R00reactive.smi
            if n2t < n2:
                command = f'shuf -n {n2t} {reactive_file} > {intermediate_file} '
                command += f'&& mv {intermediate_file} {reactive_file}'
                os.system(command)
            if n1 + n2t * n3 > self.n:
                reduce_output = True
        # determine if the reaction must be performed in parallel through qsub and eventually run it
        if n2t * n3 > self.chunksize:
            # run reaction in multiple nodes through qsub
            print('INFO::: Running job in the HPC')
            print('')
            if n3 >= self.chunksize:  # assumes that n3 must be always less than self.chunksize
                print(f'there are more building blocks in node {node} than the required number of compounds to ')
                print(f'enumerate {self.chunksize}. Please increase the required number of compounds to enumerate')
                sys.exit(1)
            n2chunk = self.chunksize // n3
            command = f'split -l {n2chunk} -a 4 {reactive_file} {reactive_file}.'
            os.system(command)
            # rename the split products and create the lists files
            oldfiles = os.listdir(self.wfolder)
            oldfiles = [file for file in oldfiles if file.startswith('R00reactive.smi.')]
            products_files = [file[-4:] + 'PRD.smi' for file in oldfiles]
            reactive_files = [file[-4:] + 'R00reactive.smi' for file in oldfiles]
            stems = [file.replace('.smi', '') for file in products_files]
            oldfiles = [os.path.join(self.wfolder, file) for file in oldfiles]
            products_files = [os.path.join(self.wfolder, file) for file in products_files]
            reactive_files = [os.path.join(self.wfolder, file) for file in reactive_files]
            stems = [os.path.join(self.wfolder, file) for file in stems]
            for oldfile, file in zip(oldfiles, reactive_files):
                os.rename(oldfile, file)
            # create the qsub config file
            config_file = os.path.join(self.wfolder, 'qsub.config')
            with open(config_file, 'w') as f:
                for qstem, qreactive_file in zip(stems, reactive_files):
                    command = f'{self.trxn}'
                    command += f' -r {reaction} -z i -m RMX -M RMX -W "+" -L -i smi -o smi'
                    command += f' -S {qstem} {qreactive_file} {self.nodes[node].bbsfile}\n'
                    f.write(command)
            # create the sge script
            sge_script = self.create_sge_script()
            # create the qsub script
            numjobs = len(stems)
            qsub_script = self.create_qsub_script(numjobs, sge_script)
            # run the qsub
            command = f'{qsub_script} {config_file}'
            os.system(command)
            # combine all the files in a single file and remove the files
            command = f'cat {" ".join([file for file in products_files])} > {products_file}'
            os.system(command)
            command = f'rm {" ".join([file for file in products_files])}'
            os.system(command)
            command = f'rm {" ".join([file for file in reactive_files])}'
            os.system(command)
        else:
            # run reaction in a single core
            print('')
            command = f'{self.trxn}'
            command += f' -r {reaction} -z i -m RMX -M RMX -W "+" -L -i smi -o smi'
            command += f' -S {stem} {reactive_file} {self.nodes[node].bbsfile}'
            os.system(command)
        # combine products_file and inert file and remove unnecesary files
        command = f'cat {products_file} {inert_file} > {scaffolds_file} && '
        command += f'rm {reactive_file} {inert_file} {products_file}'
        os.system(command)
        # reduce the size of the produts_file if required
        if reduce_output:
            command = f'shuf -n {self.n} {scaffolds_file} > {intermediate_file} '
            command += f'&& mv {intermediate_file} {scaffolds_file}'
            os.system(command)

    def cycle_node(self, edge, cycle_bond_order):
        """conducts a cyclization reaction for the specified edge
        edge: int: index for this edge
        cycle_bond_order: int, bond order of the bond to be created"""
        print('')
        print(f'INFO::: cycling edge {edge}')
        # set file names
        scaffolds_file = os.path.join(self.wfolder, 'R00.smi')
        reactive_file = os.path.join(self.wfolder, 'R00reactive.smi')
        inert_file = os.path.join(self.wfolder, 'R00inert.smi')
        products_file = os.path.join(self.wfolder, 'PRD.smi')
        stem = os.path.join(self.wfolder, 'PRD')
        # keep only the lines that contain either 0 or 2 instances of the cycle_edge tag
        if edge < 10:
            pattern_reactive = f'"\[{edge}[a-z,A-Z].*\[{edge}[a-z,A-Z]"'
            pattern_inert = f'"\[{edge}[a-z,A-Z]"'
        else:
            pattern_reactive = f'"\[{edge}.*\[{edge}"'
            pattern_inert = f'"\[{edge}"'
        command = f'grep {pattern_reactive} {scaffolds_file} > {reactive_file}'
        os.system(command)
        command = f'grep -v {pattern_inert} {scaffolds_file} > {inert_file}'
        os.system(command)
        reaction = self.create_cyclization(edge, cycle_bond_order)
        # determine if the reaction must be performed in parallel through qsub and eventually run it
        n1 = self.count_compounds(reactive_file)
        if n1 > self.chunksize:
            # run reaction in multiple nodes through qsub
            print('INFO::: Running job in the HPC')
            print('')
            n1chunk = self.chunksize // n1
            command = f'split -l {n1chunk} -a 4 {reactive_file} {reactive_file}.'
            os.system(command)
            # rename the split products and create the lists files
            oldfiles = os.listdir(self.wfolder)
            oldfiles = [file for file in oldfiles if file.startswith('R00.smi.')]
            reactive_files = [file[-4:] + 'R00reactive.smi' for file in oldfiles]
            products_files = [file[-4:] + 'PRD.smi' for file in oldfiles]
            stems = [file.replace('.smi', '') for file in products_files]
            oldfiles = [os.path.join(self.wfolder, file) for file in oldfiles]
            reactive_files = [os.path.join(self.wfolder, file) for file in reactive_files]
            products_files = [os.path.join(self.wfolder, file) for file in products_files]
            stems = [os.path.join(self.wfolder, file) for file in stems]
            for oldfile, file in zip(oldfiles, reactive_files):
                os.rename(oldfile, file)
            # create the qsub config file
            config_file = os.path.join(self.wfolder, 'qsub.config')
            with open(config_file, 'w') as f:
                for qstem, qreactive_file in zip(stems, reactive_files):
                    command = f'{self.trxn}'
                    command += f' -r {reaction} -z i -m RMX -M RMX -W "+" -i smi -o smi'
                    command += f' -S {qstem} {qreactive_file}'
                    f.write(command)
            # create the sge script
            sge_script = self.create_sge_script()
            # create the qsub script
            numjobs = len(stems)
            qsub_script = self.create_qsub_script(numjobs, sge_script)
            # run the qsub
            command = f'{qsub_script} {config_file}'
            os.system(command)
            # combine all the files in a single file and remove the files
            command = f'cat {" ".join([file for file in products_files])} > {products_file}'
            os.system(command)
            command = f'rm {" ".join([file for file in products_files])}'
            os.system(command)
            command = f'rm {" ".join([file for file in reactive_files])}'
            os.system(command)
        else:
            print('')
            command = f'{self.trxn}'
            command += f' -r {reaction} -z i -m RMX -M RMX -W "+" -i smi -o smi'
            command += f' -S {stem} {reactive_file}'
            os.system(command)
        command = f'cat {products_file} {inert_file} > {scaffolds_file} && '
        command += f'rm {reactive_file} {inert_file} {products_file}'
        os.system(command)

    def par_quality_control(self):
        """runs quality control of par"""
        edges = list(range(len(self.par)))
        for i, node in enumerate(self.par):
            if 'node' not in node.keys():
                print('config file quality control:')
                print(f'node {i}: missing "node" key. Exiting.')
                sys.exit(1)
            if type(node['node']) != int:
                print('config file quality control:')
                print(f'node {i}: value of "node" must be an int. Exiting.')
                sys.exit(1)
            if i != node['node']:
                print('config file quality control:')
                print(f'node index of node {i} is {node["node"]}. It must be {i}. Exiting.')
                sys.exit(1)
            if 'edges' not in node.keys():
                print('config file quality control:')
                print(f'node {i}: missing "edges" key. Exiting.')
                sys.exit(1)
            if type(node['edges']) != list:
                print('config file quality control:')
                print(f'node {i}: value of "edges" must be a list. Exiting.')
                sys.exit(1)
            if len(set(edges).intersection(set(node['edges']))) != len(node['edges']):
                print('config file quality control:')
                print(f'node {i} contains edges out of scope. Exiting.')
                sys.exit(1)

            if 'out_edges' not in node.keys():
                print('config file quality control:')
                print(f'node {i}: missing "out_edges" key. Exiting.')
                sys.exit(1)
            if type(node['out_edges']) != list:
                print('config file quality control:')
                print(f'node {i}: value of "out_edges" must be a list. Exiting.')
                sys.exit(1)
            if len(node['edges']) != len(node['out_edges']):
                print('config file quality control:')
                print(f'number of members in edges list in node {i} must match number of members in out_edges list {len(node["edges"])} is not equal to {len(node["out_edges"])}. Exiting.')
                sys.exit(1)
            for j, k in zip(node['edges'], node['out_edges']):
                if k <= j and k != 0:
                    print('config file quality control:')
                    print(f'node {i}: values of out_edges must be 0 or a value that is higher to the corresponding value in the edges list. Exiting.')
                    sys.exit(1)
            if node['node'] != 0 and not node['node'] in node['edges']:
                print('config file quality control:')
                print(f'node {i}: {i} must be included in the node edges ({node["edges"]}).')
                sys.exit(1)
            if 'cycle_edges' not in node.keys():
                print('config file quality control:')
                print(f'node {i}: missing "cycle_edges" key. Exiting.')
                sys.exit(1)
            if type(node['cycle_edges']) != list:
                print('config file quality control:')
                print(f'node {i}: value of "cycle_edges" must be a list. Exiting.')
                sys.exit(1)
            if len(set(edges).intersection(set(node['cycle_edges']))) > 0:
                print('config file quality control:')
                print(f'node {i} contains cycle edges in scope. Exiting.')
                sys.exit(1)
            if 'bbs_file' not in node.keys():
                print('config file quality control:')
                print(f'node {i}: missing "nodes_file" key. Exiting.')
                sys.exit(1)
            if type(node['bbs_file']) != str:
                print('config file quality control:')
                print(f'node {i}: value of "bbs_file" must be a list. Exiting.')
                sys.exit(1)
            if not os.path.isfile(os.path.expandvars(node['bbs_file'])):
                print(os.path.expandvars(node['bbs_file']))  # todo remove after testing
                print('config file quality control:')
                print(f'node {i}: bbs file {node["bbs_file"]} does not exist. Exiting.')
                sys.exit(1)
            if 'preparations' not in node.keys():
                print('config file quality control:')
                print(f'node {i}: missing "preparations" key. Exiting.')
                sys.exit(1)
            if type(node['preparations']) != list:
                print('config file quality control:')
                print(f'node {i}: value of "preparations" key must be a list. Exiting.')
                sys.exit(1)
            for j, preparation in enumerate(node['preparations']):
                if type(preparation) != dict:
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: all elements of "preparations" key must be dictionaries. Exiting.')
                    sys.exit(1)
                if 'edges' not in preparation.keys():
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: missing "edges" key. Exiting.')
                    sys.exit(1)
                if type(preparation['edges']) != list:
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: value of "edges" must be list. Exiting.')
                    sys.exit(1)
                if len(set(preparation['edges']).intersection(set(node['edges'] + node['cycle_edges']))) != len(preparation['edges']):
                    print('config file quality control:')
                    print(f'node {i}, preparation {j} contains edges out of scope. Exiting.')
                    sys.exit(1)
                if 'keep_unprepared' not in preparation.keys():
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: missing "keep_unprepared" key. Exiting.')
                    sys.exit(1)
                if type(preparation['keep_unprepared']) != bool:
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: value of "keep_unprepared" must be bool. Exiting.')
                    sys.exit(1)
                if node['node'] in preparation['edges'] and preparation['keep_unprepared']:
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: when the preparation generates the node isotope, keep_unprepared must be set to false. Exiting.')
                    sys.exit(1)
                if 'reaction' not in preparation.keys():
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: missing "reaction" key. Exiting.')
                    sys.exit(1)
                if type(preparation['reaction']) != str:
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: value of "reaction" must be str. Exiting.')
                    sys.exit(1)
                if not os.path.isfile(os.path.expandvars(preparation['reaction'])):
                    print(os.path.expandvars(preparation['reaction']))  # todo remoce after testing
                    print('config file quality control:')
                    print(f'node {i}, preparation {j}: reaction file {preparation["reaction"]} does not exist. Exiting.')
                    sys.exit(1)

    def run_graph(self, target_file=None):
        """Runs all the reactions specified in the graph
        target_file: str: path to a target file, The target file is a text file
            with columns separated by spaces and no header containing ids of the building blocks to add.
            The columns are not headed. If not passed then the enumeration goes node by node
        """
        # generate the node objects and append them in the node list, then prepare bbs if required
        for node in self.par:
            self.add_node(node)
            self.nodes[-1].add_preparations()
            self.success = self.success and self.nodes[-1].success
        # conduct enumeration of the main graph
        if target_file is not None:
            self.enumerate_target_file(target_file)
            return
        for node_index in range(len(self.nodes)):
            if node_index != 0:
                self.join_node(node_index)
        # conduct cyclizations
        cycles = []
        for node_index, node in enumerate(self.nodes):
            if len(node.cycle_edges) > 0:
                for cycle_edge, cycle_bond_order in zip(node.cycle_edges, node.cycle_bond_orders):
                    if cycle_edge not in cycles:
                        cycles.append(cycle_edge)
                        self.cycle_node(cycle_edge, cycle_bond_order)
        os.rename(os.path.join(self.wfolder, "R00.smi"), os.path.join(self.wfolder, "enumeration.smi"))
        os.rename(os.path.join(self.wfolder, "RP00.smi"), os.path.join(self.wfolder, "R00.smi"))


if __name__ == '__main__':
    print(__version__)
    print(__author__)
