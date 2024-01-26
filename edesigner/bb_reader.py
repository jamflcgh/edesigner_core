# -*- coding: utf-8 -*-
# bb_reader
# Jose Alfredo Martin

__version__ = 'bb_reader.v.12.0.0'
__author__ = 'Alfredo Martin 2023'

import os
import pandas as pd
import numpy as np
import copy
import _pickle as pic
from tqdm import tqdm


class BBReader:
    """class that handles the reading and processing of building blocks and update s the already created
    BBTs list. It also creates BB files required to enumeration of designs."""

    def __init__(self, wf, runfolder, BBTs, db, fg, antifg, calcfg, bblim, log, smi=None, verbose=False, debug=False):
        """
        Initiallizes the class
        :param wf: str: working folder
        :param runfolder: str: folder tor the current run
        :param BBTs: list of BBT objects
        :param dbpar: instance of the Parameters class corresponding to the db parameters
        :param fg: instance of the Parameters class corresponding to the fg parameters
        :param antifg: instance of the Parameters class corresponding to the antifg parameters
        :param calcfg: instance of the Parameters class corresponding to the calcfg parameters
        :param bblim: instance of the Parameters class corresponding to the bblim parameters
        :param log: instance of logger class
        :param smi: str: path to a smiles file containing all building blocks
        :param verbose: report Lyllymol detailed info
        :param debug: report file lengths and headers for debugging purposes
        :return None
        """
        self.wf = wf
        self.success = True
        self.reason = ""
        self.BBTs = BBTs
        self.db = db
        self.smi = smi
        self.fg = fg
        self.antifg = antifg
        self.calcfg = calcfg
        self.bblim = bblim
        self.log = log
        self.verbose = verbose
        self.debug = debug
        self.runfolder = runfolder
        if os.path.isdir(os.environ['LILLYMOL_EXECUTABLES']):
            self.tsubstructure = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'tsubstructure')
            self.fileconv = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'fileconv')
            self.unique_molecules = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'unique_molecules')
            self.iwdescr = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'iwdescr')
            self.iwcut = os.path.join(os.environ['LILLYMOL_EXECUTABLES'], 'iwcut')
        else:
            self.tsubstructure = 'tsubstructure.sh'
            self.fileconv = 'fileconv.sh'
            self.unique_molecules = 'unique_molecules.sh'
            self.iwdescr = 'iwdescr.sh'
            self.iwcut = 'iwcut.sh'

    def _annotate_bbs(self, bb_file, how='fg'):
        """
        Annotates the file of bbs file with the vector of fgs for anti-fgs (how='antifg') or fgs (how='fg')
        :param bb_file: str: path to smiles file containing bbs (two columns space separated with smiles an id fields)
        :param how: str: what type of fgs to annotate
        :return: str: path to the annotated file
        """
        qrys = []
        protos = []
        if how == 'antifg':
            for fg in self.antifg.par:
                qrys += [os.path.expandvars(file) for file in fg['base_queries'] if file.endswith('qry')]
                protos += [os.path.expandvars(file) for file in fg['base_queries'] if file.endswith('proto')]
        else:
            for fg in self.fg.par:
                qrys += [os.path.expandvars(file) for file in fg['base_queries'] if file is not None and file.endswith('qry')]
                protos += [os.path.expandvars(file) for file in fg['base_queries'] if file is not None and file.endswith('proto')]
        with open(os.path.join(self.wf, "qrys.txt"), 'w') as f:
            for file in qrys:
                f.write(file + '\n')
        with open(os.path.join(self.wf, "protos.txt"), 'w') as f:
            for file in protos:
                f.write(file + '\n')
        if len(qrys) > 0:
            qry_stamp = f'-q F:{os.path.join(self.wf, "qrys.txt")}'
        else:
            qry_stamp = ''
        if len(protos) > 0:
            proto_stamp = f'-q PROTOFILE:{os.path.join(self.wf, "protos.txt")}'
        else:
            proto_stamp = ''
        if len(qrys) + len(protos) == 0:
            self.log.update("ERROR::: no fg qry or proto files found")
        command = f'{self.tsubstructure}'
        command += f' {qry_stamp} {proto_stamp}'
        if self.verbose:
            command += f' -v'
        command += f' -A D -u -r -a -l {bb_file} > {os.path.join(self.wf, "annotated.txt")}'
        os.system(command)
        return os.path.join(self.wf, "annotated.txt")

    def _read_bbs(self):
        """
        reads, canonicalize, deduplicate, desalt and compute properties for all bbs
        :return: tuple of str: file1, file2
           file1: path to file smiles and id of returned compounds
           file2: path to file with header containing nha and nhbd for the same compounds in the smiles file
        """
        # 1. Read compounds from a single smiles
        if self.smi is not None:
            self.log.update(f'INFO::: Working in bbs file {self.smi}')
            file = os.path.expandvars(self.smi)
            _, name = os.path.split(self.smi)
            name = name.replace('.smi', '')
            command = f'{self.fileconv}'
            command += f' -O B -i smi -o usmi -V -E autocreate -f lod -S'
            command += f' {os.path.join(self.wf, "c_desalted")} {file}'
            os.system(command)
            df = pd.read_csv(os.path.join(self.wf, "c_desalted.smi"), sep=' ', header=None, names=['smiles', 'id'])
            df['id'] = df['id'].apply(lambda x: name + ':' + x)
            df[['smiles', 'id']].to_csv(os.path.join(self.wf, "c_desalted.smi"), sep=' ', header=False, index=False)
            self.log.update(f'INFO::: Adding {df.shape[0]} compounds to the desalted compounds pool')
            command = f'cat {os.path.join(self.wf, "c_desalted.smi")} >> {os.path.join(self.wf, "desalted.smi")}'
            os.system(command)
        # 1. Read compounds from all the sets, desalt them and tag them with its source. Creates a single file
        else:
            for db in self.db.par:
                self.log.update(f'INFO::: Working in bbs file {db["db"]}')
                file = os.path.expandvars(db['filename'])
                command = f'{self.fileconv}'
                command += f' -O B -i smi -o usmi -V -E autocreate -f lod -S'
                command += f' {os.path.join(self.wf, "c_desalted")} {file}'
                os.system(command)
                df = pd.read_csv(os.path.join(self.wf, "c_desalted.smi"), sep=' ', header=None, names=['smiles', 'id'])
                df['id'] = df['id'].apply(lambda x: db['db'] + ':' + x)
                df[['smiles', 'id']].to_csv(os.path.join(self.wf, "c_desalted.smi"), sep=' ', header=False, index=False)
                self.log.update(f'INFO::: Adding {df.shape[0]} compounds to the desalted compounds pool')
                command = f'cat {os.path.join(self.wf, "c_desalted.smi")} >> {os.path.join(self.wf, "desalted.smi")}'
                os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "desalted.smi")}')
            os.system(f'head -n 1 {os.path.join(self.wf, "desalted.smi")}')
        # 2. Eliminate repeated smiles
        self.log.update('INFO::: eliminating duplicated smiles')
        command = f'{self.unique_molecules}'
        command += f' -l -a -I -A D -g all -S {os.path.join(self.wf, "unique")}'
        command += f' {os.path.join(self.wf, "desalted.smi")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "unique.smi")}')
            os.system(f'head -n 1 {os.path.join(self.wf, "unique.smi")}')
        command = f'rm {os.path.join(self.wf, "desalted.smi")}'
        os.system(command)
        # 3. filter among number of heavy atoms and rotatable bonds
        self.log.update('INFO::: level 1 filetering by heavy atoms and rotatable bonds')
        max_atoms = self.bblim.par['raw_na_filter']
        max_rb = self.bblim.par['raw_rb_filter']
        if max_atoms is None:
            max_atoms = 10000
        if max_rb is None:
            max_rb = 10000
        command = f'{self.iwdescr}'
        command += f' -F "w_natoms<{max_atoms}" -F "w_rotbond<{max_rb}" {os.path.join(self.wf, "unique.smi")}'
        command += f' > {os.path.join(self.wf, "f_unique.smi")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "f_unique.smi")}')
            os.system(f'head -n 1  {os.path.join(self.wf, "f_unique.smi")}')
        # 4. annotate compounds with all properties
        self.log.update('INFO::: Annotating compounds with properties')
        command = f'{self.iwdescr}'
        command += f' {os.path.join(self.wf, "f_unique.smi")} > {os.path.join(self.wf, "f_unique.txt")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "f_unique.txt")}')
            os.system(f'head -n 1 {os.path.join(self.wf, "f_unique.txt")}')
        # 4. eliminate all properties except for number of heavy atoms and number of rotatable bon
        command = f'{self.iwcut}'
        command += f' -d w_natoms,w_rotbond {os.path.join(self.wf, "f_unique.txt")}'
        command += f' > {os.path.join(self.wf, "s_unique.txt")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "s_unique.txt")}')
            os.system(f'head -n 1 {os.path.join(self.wf, "s_unique.txt")}')
        # 5. separate the properties file so it only contains properties
        command = f'cut {os.path.join(self.wf, "s_unique.txt")} -d " " -f 2,3 > {os.path.join(self.wf, "iwdata.txt")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "iwdata.txt")}')
            os.system(f'head -n 1 {os.path.join(self.wf, "iwdata.txt")}')
        return os.path.join(self.wf, "f_unique.smi"), os.path.join(self.wf, "iwdata.txt")

    def _add_anitfg(self, in_smiles, in_data):
        """
        creates annotated files for bbs
        in_smiles: str: path to file smiles and id of returned compounds
        in_data: str: path to file with header containing nha and nhbd for the same compounds in the smiles file
        :return: str: file with annotated antifg
        """
        self.log.update('INFO::: Annotating anti FG')
        # 1. compute the vector for the antitargets
        annotated = self._annotate_bbs(in_smiles, how='antifg')
        if self.debug:
            os.system(f'wc -l {annotated}')
            os.system(f'head -n 1 {annotated}')
        # 2. add smiles header
        os.system(f'echo "smiles id" > {os.path.join(self.wf, "header")}')
        os.system(f'cat {os.path.join(self.wf, "header")} {in_smiles} > {os.path.join(self.wf, "h_smiles")}')
        os.system(f'cut -d " " -f 1 {os.path.join(self.wf, "h_smiles")} > {os.path.join(self.wf, "h_smiles.txt")}')
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "h_smiles.txt")}')
            os.system(f'head -n 1 {os.path.join(self.wf, "h_smiles.txt")}')
        # 3. add properties to file
        command = f'paste -d " " {os.path.join(self.wf, "h_smiles.txt")} {annotated} {in_data} '
        command += f'> {os.path.join(self.wf, "processed.txt")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "processed.txt")}')
            os.system(f'head -1 {os.path.join(self.wf, "processed.txt")}')
        return os.path.join(self.wf, "processed.txt")

    def _add_calc_fgs(self, file, how='fg'):
        """
        calculates the calc_fgs
        :param file: file with annotations
        :param how: str: which ar the base fgs ['antifg' | 'fg']
        :return: df: pandas dataframe
        """
        self.log.update(f'INFO::: Adding calc_FG to {how}')
        df = pd.read_csv(file, sep=' ')
        self.log.update(f'INFO::: Working on {df.shape[0]} compounds')
        fields = set(df.columns.tolist())
        if how == 'fg':
            domain = copy.deepcopy(self.fg.par)
        else:
            domain = copy.deepcopy(self.antifg.par)
        for fg in domain:
            if fg['name'].startswith('calc_'):
                for calcfg in self.calcfg.par:
                    if fg['name'] == calcfg['name']:
                        rule_add = calcfg['rule_add']
                        rule_substract = calcfg['rule_substract']
                        break
                else:
                    self.log.update(f'WARNING::: could not find calculated fg {fg["name"]}')
                    continue
                c1 = len(set(rule_add).intersection(fields)) == len(rule_add)
                c2 = len(set(rule_substract).intersection(fields)) == len(rule_substract)
                if c1 and c2:
                    if self.verbose:
                        self.log.update(f'INFO::: adding fg {fg["name"]}')
                    df[fg['name']] = 0
                    for field in rule_add:
                        df[fg['name']] += df[field]
                    for field in rule_substract:
                        df[fg['name']] -= df[field]
        return df

    def _filter_by_fgs(self, df, how='fg'):
        """
        filters dataframe by the number on instances of fgs for every given record. If filtering by antifgs remove all
        fg counts columns if not keeps only the ones with 1 to 3 FG counts of different FGs.
        :param df: pandas dataframe
        :param how: str: which ar the base fgs ['antifg' | 'fg']
        :return: pandas dataframe
        """
        self.log.update(f'INFO::: Filtering componds by {how}')
        self.log.update(f'INFO::: {df.shape[0]} compounds to process')
        if how == 'fg':
            domain = copy.deepcopy(self.fg.par)
        else:
            domain = copy.deepcopy(self.antifg.par)
        base = ['smiles', 'Name', 'w_natoms', 'w_rotbond']
        expected_fgs = [fg['name'] for fg in domain if fg['name'] if fg['name'].upper() != 'NO_FG']
        keep_fields = base + expected_fgs
        remove = [item for item in df.columns.tolist() if item not in keep_fields]
        df.drop(inplace=True, columns=remove)
        fgs = [item for item in df.columns.tolist() if item not in base]

        if len(fgs) != len(expected_fgs):
            self.log.update(f"WARNING::: The number of fields in the dataframe ({len(fgs)}) does not match the number of fgs in {how} ({len(domain)})")
        if how == 'fg':
            keep1 = np.greater(np.array(df[fgs]).sum(axis=1), 0)  # at least onf FG
            keep2 = np.less_equal(np.array(df[fgs]).sum(axis=1), 3)  # at most three FG
            keep3 = np.less(np.array(df[fgs]).max(axis=1), 2)  # all fgs are different
            keep = np.logical_and(keep1, keep2, keep3)
        else:
            keep = np.equal(np.array(df[fgs]).sum(axis=1), 0)

        if self.verbose:
            self.log.update('INFO::: list of excluded compounds')
            for index, row in df.iterrows():
                if not keep[index]:
                    linea = row['smiles']
                    for column in fgs:
                        if row[column] > 0:
                            linea += ' ' + column
                    self.log.update(linea, np.array(df[fgs])[index])  # check compounds excluded
        df = df[keep].copy()
        self.log.update(f'INFO::: {df.shape[0]} compounds remaining')
        return df

    def _add_fg(self, df):
        """
        adds the fgs to the dataframe
        :param df: pandas dataframe
        :return: str: filename
        """
        df[['smiles', 'Name']].to_csv(os.path.join(self.wf, "f_unique.smi"), sep=' ', header=False, index=False)
        df[['smiles']].to_csv(os.path.join(self.wf, "smiles.txt"), sep=' ', index=False)
        df[['w_natoms', 'w_rotbond']].to_csv(os.path.join(self.wf, "iwdata.txt"), sep=' ', index=False)
        file = self._annotate_bbs(os.path.join(self.wf, "f_unique.smi"), how='fg')
        command = f'paste -d " " {os.path.join(self.wf, "smiles.txt")} {file} '
        command += f'{os.path.join(self.wf, "iwdata.txt")} > {os.path.join(self.wf, "processed.txt")}'
        os.system(command)
        if self.debug:
            os.system(f'wc -l {os.path.join(self.wf, "processed.txt")}')
            os.system(f'head -1 {os.path.join(self.wf, "processed.txt")}')
        return os.path.join(self.wf, "processed.txt")

    def _annotate_with_bbts(self, df):
        """
        annotates each bb with the index of the bbt it belongs to
        :param df: pandas dataframe
        :return: pandas dataframe
        """
        self.log.update('INFO::: Annotatingg BBTs')
        self.log.update(f'INFO::: annotating {df.shape[0]} bbs')
        # Annotate with index
        cum_found = np.zeros(df.shape[0]).astype('int32')
        fgs = [item['name'] for item in self.fg.par if item['name'].upper() != 'NO_FG']
        arr = df[fgs].to_numpy()  # numpy array with all items in the dataframe
        for i in tqdm(range(len(self.BBTs))):
            bbt_query = np.array(self.BBTs[i].BBT_long[1:])
            found = np.equal(np.equal(arr, bbt_query[None, :]).sum(axis=1), bbt_query.shape[0])
            cum_found += found * self.BBTs[i].index
        df['bbt_index'] = cum_found
        df = df[df['bbt_index'] > 0].copy()
        self.log.update(f'INFO::: {df.shape[0]} bbs remain after annotation')
        self.log.update(f'INFO::: {len(df.bbt_index.unique().tolist())} unique BBTs found in compounds')
        arr = df[fgs].to_numpy()
        df.drop(inplace=True, columns=fgs)
        # compute the effective number of heavy atoms and rotatable bonds
        self.log.update(f'INFO::: filtering by effective number heavy atoms and rotatable bonds')
        excess_rb = np.array([item['excess_rb'] for item in self.fg.par][1:])
        atom_dif = np.array([item['atom_dif'] for item in self.fg.par][1:])
        delta_rb = (arr * excess_rb[None, :]).sum(axis=1)
        delta_ha = (arr * atom_dif[None, :]).sum(axis=1)
        df['delta_rb'] = delta_rb
        df['delta_ha'] = delta_ha
        df['nrb'] = np.maximum(0, np.array(df['w_rotbond']) - delta_rb)
        df['nha'] = np.maximum(0, np.array(df['w_natoms'] + delta_ha))
        df = df[df['nrb'] <= self.bblim.par['rb_filter']].copy()
        df = df[df['nha'] >= self.bblim.par['min_bb_na']].copy()
        df = df[df['nha'] <= self.bblim.par['max_bb_na']].copy()
        self.log.update(f'INFO::: {df.shape[0]} bbs remain after filtering')
        self.log.update(f'INFO::: {len(df.bbt_index.unique().tolist())} unique BBTs remain after filtering')
        df.drop(inplace=True, columns=['w_natoms', 'w_rotbond', 'nrb'])
        df.sort_values(by=['nha'], inplace=True)
        return df

    def _report_compound_files(self, df):
        """report_compound_files create a smi file for each BBT and put in the file all compounds in compound
        assigned to that BBT_index.
        :param df: pandas dataframe
        :param df: pandas dataframe
        :return: None
        returns : None"""
        self.log.update('INFO::: Reporting compounds...')
        counter = 0
        for bbt in tqdm(self.BBTs):
            cdf = df[df['bbt_index'] == bbt.index].copy()
            if cdf.shape[0] > 0:
                counter += 1
                cdf['id'] = cdf.apply(lambda row: ':'.join([str(row['nha']), str(row['Name'])]), axis=1)
                cdf = cdf[['smiles', 'id']].copy()
                cdf.to_csv(os.path.join(self.runfolder, 'comps', str(bbt.index) + '.smi'), sep=' ',
                           index=False, header=False)
        self.log.update(f'INFO::: Reported {counter} BBT files...')

    def _update_bbts(self, df):
        """
        updates the BBTs object
        :param df:
        :return: list of BBT objects
        """
        self.log.update('INFO::: Updating BBTs...')
        df['count'] = 0
        dfc = df.groupby(['bbt_index', 'nha'])['count'].count().reset_index()
        n_comps = np.zeros((len(self.BBTs), self.bblim.par['max_bb_na'] + 1)).astype('int32')
        for index, row in dfc.iterrows():
            n_comps[row['bbt_index']][row['nha']] = row['count']
        for i in tqdm(range(len(self.BBTs))):
            self.BBTs[i].maxna = self.bblim.par['max_bb_na']
            self.BBTs[i].n_compounds = n_comps[self.BBTs[i].index, :].flatten()
            if self.BBTs[i].n_compounds.sum() > 0:
                self.BBTs[i].smiles_example = df[df['bbt_index'] == self.BBTs[i].index]['smiles'].tolist()[0]
                for j, nat in enumerate(self.BBTs[i].n_compounds):
                    if nat > 0:
                        if self.BBTs[i].min_atoms == 0:
                            self.BBTs[i].min_atoms = j
                        self.BBTs[i].max_atoms = j
        self.BBTs.sort(key=lambda x: x.n_compounds.sum(), reverse=True)
        self.BBTs.sort(key=lambda x: x.BBT_multi)
        for i in range(len(self.BBTs)):
            self.BBTs[i].order = i
        self.BBTs.sort(key=lambda x: x.index)
        self.log.update('INFO::: BBTs updated...')
        counter = 0

    def _report_bbt_info(self):
        """creates a csv file to inspect that bbs have been assigned correctly. It includes FGs required, smiles and
        other information
        returns: None"""
        report = []
        for bbt in self.BBTs:
            if bbt.n_compounds.sum() > 0:
                this_dict = dict()
                this_dict['FGs'] = ':'.join([item for item in bbt.BBT_name if item.upper() != 'NO_FG'])
                this_dict['n_compounds'] = bbt.n_compounds.sum()
                this_dict['smiles'] = bbt.smiles_example
                this_dict['index'] = bbt.index
                this_dict['order'] = bbt.order
                this_dict['multi'] = bbt.BBT_multi
                this_dict['max_atoms'] = bbt.max_atoms
                this_dict['min_atoms'] = bbt.min_atoms
                report.append(this_dict)
        df = pd.DataFrame(report)
        df.to_csv(os.path.join(self.wf, 'bbt_report.csv'), index=False)

    def run(self):
        """
        runs the complete workflow
        :return: str: path to the bbs file
        """
        file1, file2 = self._read_bbs()
        file = self._add_anitfg(file1, file2)
        df = self._add_calc_fgs(file, how='antifg')
        df = self._filter_by_fgs(df, how='antifg')
        file = self._add_fg(df)
        df = self._add_calc_fgs(file, how='fg')
        df = self._filter_by_fgs(df, how='fg')
        df = self._annotate_with_bbts(df)
        self._report_compound_files(df)
        self._update_bbts(df)
        self._report_bbt_info()
        with open(os.path.join(self.wf, 'BBTs.pic'), 'wb') as f:
            pic.dump(self.BBTs, f)


if __name__ == '__main__':
    print(__version__)
    print(__author__)
