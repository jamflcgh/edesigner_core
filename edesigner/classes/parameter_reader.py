# -*- coding: utf-8 -*-
# parameter_reader
# Jose Alfredo Martin

version = 'parameter_reader.v.9.0.0'

#Python modules
import copy

############ PARAMETERS CLASS #############

class Parameters:
    def __init__(self, path, fsource='dict', how='to_list', multiple=False):
        """this method initiallizes the paramters instance
        the instance will contain a list of dictionaries, a dictionary of lists or a dictionary
        path : str (complete path and filename to the file containing the parameters
        fsource : str (format of the parameters table that is read) options are:
            dict ; Format of the file is one parameter per row and one or more values per column.
            list ; Format of the file is on parameter per column and one or more values per row
        how : str (how to read the parameters file) options are:
            to_list ; The instance will contain a list of dictionaries if multiple = True, or a dictionary if multiple = False)
            to_dict : The instance will contain a dictionary of lists if multiple = True or a dictionary if multiple = False)
        multiple : bool
            True : Each parameter is a list of what is indicated in the input file taking each set of parameters to each list item
            False : Each parameter is a value of what is inidcated in the input file, taking the first set of parameters to that value
        input file format:
            Input file must be a text tab delimited file
            first line is ignored if fsource is dict, first column is ignored if fsource is list
            second column/row contains the dictionary keys for fsource dict/list
            third column/row contains the data type (str, int, float, bool or dict (dicts contain str keys and float values separated by :) for fsource dict/list
            fourth column/row contains the list separator mark for lists for fsource dict/list
            fifth column/row contains an explanation of the data for fsource dict/list
            following columns/rows contain data (if empty should contain None or Null, cannot be blank) fsource dict/list
        returns : None"""
        self.path = path
        self.errors = []
        self.success = True
        self.par = {} # initial format will always be a dictionary of lists and will be changed afterwards depending on the values of how and multiple arguments
        self.type = {}
        self.separator = {}
        self.description = {}
        try:
            f = open(path, 'r')
        except:
            self.errors.append('    file ' + path + ' not found:')
            self.success = False
        else:
            lineas = f.readlines()
            f.close()
            row_lineas = copy.deepcopy([linea.rstrip('\r\n').split('\t') for linea in lineas])
            lineas = []
            if fsource == 'dict':
                lineas = copy.deepcopy(row_lineas)
            elif fsource == 'list':
                for j in range(len(row_lineas[0])):
                    lineas.append([])
                    for i in range(len(row_lineas)):
                        lineas[-1].append(copy.deepcopy(row_lineas[i][j]))
            else:
                self.success = False
                self.errors.append('    fsource not stated properly')
            if self.success:
                lineas = copy.deepcopy(lineas[1:]) # eliminates the first row
                for linea in lineas:
                    self.type[linea[0]] = linea[1]
                    self.separator[linea[0]] = linea[2]
                    self.description[linea[0]] = linea[3]
                    self.par[linea[0]] = linea[4:]
                    # retrieve and split the data
                    if self.separator[linea[0]] == '': # applied only to inmutable vaules
                       for i in range(len(self.par[linea[0]])): 
                           self.par[linea[0]][i] = [self.par[linea[0]][i]] # inmutable values are converted to lists to assign the data types at the same time than lists
                    elif len(self.separator[linea[0]]) == 1 or (len(self.separator[linea[0]]) > 1 and ';' not in self.separator[linea[0]]): # applied to lists
                        for i in range(len(self.par[linea[0]])):
                            self.par[linea[0]][i] = self.par[linea[0]][i].split(self.separator[linea[0]])
                    elif len(self.separator[linea[0]]) > 1 and ';' in self.separator[linea[0]]: # applied to lists of lists
                        this_separator = self.separator[linea[0]].split(';')
                        for i in range(len(self.par[linea[0]])):
                            self.par[linea[0]][i] = self.par[linea[0]][i].split(this_separator[0])
                            for j in range(len(self.par[linea[0]][i])):
                                self.par[linea[0]][i][j] = self.par[linea[0]][i][j].split(this_separator[1])
                    # convert data to the appropriate data types
                    if len(self.separator[linea[0]]) < 2 or (len(self.separator[linea[0]]) > 1 and ';' not in self.separator[linea[0]]): # for inmutables and lists
                        for i in range(len(self.par[linea[0]])):
                            for j in range(len(self.par[linea[0]][i])):
                                if self.par[linea[0]][i][j].upper() in ['', 'NONE', 'NULL']: # everything that is empty or marked as none will be set to None
                                    self.par[linea[0]][i][j] = None
                                elif self.type[linea[0]].upper() == 'STR':
                                    pass
                                elif self.type[linea[0]].upper() == 'INT':
                                    self.par[linea[0]][i][j] = int(self.par[linea[0]][i][j])
                                elif self.type[linea[0]].upper() == 'FLOAT':
                                    self.par[linea[0]][i][j] = float(self.par[linea[0]][i][j])
                                elif self.type[linea[0]].upper() == 'BOOL':
                                    self.par[linea[0]][i][j] = self.par[linea[0]][i][j].upper() == 'Y' or self.par[linea[0]][i][j].upper() == 'TRUE'
                                elif self.type[linea[0]].upper() == 'DICT':
                                    if j == 0:
                                        this_dict = {}
                                    this_dict[self.par[linea[0]][i][j].split(':')[0]] = float(self.par[linea[0]][i][j].split(':')[1])
                                    if j == len(self.par[linea[0]][i]) -1:
                                        self.par[linea[0]][i] = copy.deepcopy(this_dict)
                                else:
                                    self.errors.append('    unable to parse type ' + self.type[linea[0]] + ' with data ' + str(self.par[linea[0]][i][j]) + ' for key ' + linea[0])
                                    self.success = False
                            if self.separator[linea[0]] == '':
                                self.par[linea[0]][i] = self.par[linea[0]][i][0] # this turns back mutable values to inmutable values
                    if len(self.separator[linea[0]]) > 1 and ';' in self.separator[linea[0]]: # for lists of lists
                        for i in range(len(self.par[linea[0]])):
                            for j in range(len(self.par[linea[0]][i])):
                                for k in range(len(self.par[linea[0]][i][j])):
                                    if self.par[linea[0]][i][j][k].upper() in ['', 'NONE', 'NULL']: # everything that is empty or marked as none will be set to None
                                        self.par[linea[0]][i][j][k] = None
                                    elif self.type[linea[0]].upper() == 'STR':
                                        pass
                                    elif self.type[linea[0]].upper() == 'INT':
                                        self.par[linea[0]][i][j][k] = int(self.par[linea[0]][i][j][k])
                                    elif self.type[linea[0]].upper() == 'FLOAT':
                                        self.par[linea[0]][i][j][k] = float(self.par[linea[0]][i][j][k])
                                    elif self.type[linea[0]].upper() == 'BOOL':
                                        self.par[linea[0]][i][j][k] = self.par[linea[0]][i][j][k].upper() == 'Y' or self.par[linea[0]][i][j][k].upper() == 'TRUE'
                                    elif self.type[linea[0]].upper() == 'DICT':
                                        if k == 0:
                                           this_dict = {} 
                                        this_dict[self.par[linea[0]][i][j][k].split(':')[0]] = float(self.par[linea[0]][i][j][k].split(':')[1])
                                        if k == len(self.par[linea[0]][i]) - 1:
                                            self.par[linea[0]][i][j] = copy.deepcopy(this_dict)
                                    else:
                                        self.errors.append('    unable to parse type ' + self.type[linea[0]] + ' with data ' + str(self.par[linea[0]][i][j][k]) + ' for key ' + linea[0])
                                        self.success = False
                if how == 'to_list': # the dictionary of lists is converted into a list of dictionaries
                    this_list = []
                    for key in list(self.par.keys()):
                        len_key_0 = len(self.par[key])
                        break
                    for i in range(len_key_0):
                        this_dict = {}
                        for key in list(self.par.keys()):
                            this_dict[key] = copy.deepcopy(self.par[key][i])
                        this_list.append(this_dict)
                    self.par = copy.deepcopy(this_list)
                    if not multiple: # if multiple is false we keep just one dict
                        self.par = self.par[0]
                if how == 'to_dict': # the format is already correct, a dictionary of lists
                    if not multiple: # if multiple is false we keep just one dict
                        for key in list(self.par.keys()):
                            self.par[key] = copy.deepcopy(self.par[key][0])


if __name__ == '__main__':
    print(version)
