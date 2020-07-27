# -*- coding: utf-8 -*-
# design
# Jose Alfredo Martin

# Python modules
import copy
import pickle as pic
import os
import time
# Local Modules
from classes.parallel import starmap_parallel

version = 'design.v.9.0.0'

# CLASS DEFINITION


class Design:
    """Design instances store attributes of an eDESIGN and methods and attributes for its creation
    and modification. The most important methods are the 'add_BBT_to_designs' that creates and returns
    a list of Design instances based on the addition of a BBT instance to the current Design instance,
    and the 'add_deprotection_to_design' method that creates a list of Design instances upon the application
    of a deptortection / scaffold addition reaction.
    """
    def __init__(self, par, BBTs, hp_index):
        """This method initiallizes the Design class
        par : instance of Parameters class (parameters)
        BBTs : list of instances of the BBT class
        hp_index : int (index of the BBT in BBTs used as headpiece)"""
        self.id = None
        self.lib_id = None
        self.n_cycles = 0
        self.bbts = [hp_index]
        self.btopology = [-1]
        self.dtopology = [-1]
        self.reactions = [0]
        self.deprotections = [0]
        self.fgs = [BBTs[hp_index].BBT[i] for i in range(3) if BBTs[hp_index].BBT[i] > 0]
        self.fg_sources = [0 for item in self.fgs]
        self.min_natoms = par.par['headpiece_na']

    def clone(self):
        """This method returns a clone of the current instance of the Design class"""
        return copy.deepcopy(self)

    def add_BBT_to_design(self, BBT, reaction_indexes, available_reaction_input, reaction_input, reaction_output, reaction, par, fg):
        """This methods will try to add a scaffold to a design. If it succeed returns a new design
        BBT : instance of BBT class (incoming BBT)
        reaction_indexes : list of int representing indexes of all the available reactions
        available_reaction_input : list of tuples with the incoming FGs for only the available reactions
        reaction_input : list of tuples with the incoming FGs for all the reactions
        reaction_output : list of tuples with the outcoming FGs for all the reactions
        reaction : instance of Parameters class (reaction)
        par : instance of Parameters class (parameters)
        fg : instance of Parameters class (fg)
        returns : resultato (list of instances of Design class)"""
        resultado = []
        if self.fgs != []: # if design is closed no further reactions can be made
            for fg_off in [bbt_fg for bbt_fg in BBT.BBT if bbt_fg != 0]: # We loop into all the non Null FGs of the incoming BBT
                for fg_on in self.fgs: # We loop into all the FGs exposed from the initial design
                    if (fg_on, fg_off) in available_reaction_input: # the FG tuple is in the list of the tuples for incoming FGs for all the available reactions
                        indices = [i for i in reaction_indexes if (reaction_input[i] == (fg_on, fg_off))] # indicates all the reactions that have been found within the available ones
                        for i in indices: # We loop in all the reactions that have been found
                            new_design = self.clone()
                            fg_on_index = new_design.fgs.index(fg_on)
                            new_design.btopology.append(new_design.fg_sources[fg_on_index]) # we update the topology of the new design
                            del new_design.fgs[fg_on_index] # In the new design the reacting FG is removed
                            del new_design.fg_sources[fg_on_index] #In the new design the source corresponding to the removed fg is removed
                            new_BBT_used = copy.deepcopy(BBT.BBT) # We make a copy of the incoming BBT FG list
                            new_BBT_used.remove(fg_off) # and then we remove the FG that is reacting
                            excluded = any([on_ex in new_design.fgs for on_ex in reaction.par[i]['excluded_on']]) # becames True if any of the FGs becoming exposed is in the excluded_on list for this reaction
                            if not excluded: # excluded if any of the fgs in reactions excluded_on for this reaction is in the remaining FGs of the incoming BBT
                                excluded = any([off_ex in new_BBT_used for off_ex in reaction.par[i]['excluded_off']])
                            if not excluded: # excluded if any of the resulting fgs for this reaction is in the self incompatibility list of any of the fgs that were already on dna except the one that reacted
                                excluded = any([any([new_fg in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs]) for new_fg in reaction_output[i]])
                            if not excluded: # excluded if any of the fgs in the incoming BBT, except the one that reacts is in the self incompatibility list of any of the fgs that were already on dna except the one that reacted
                                excluded = any([any([new_fg in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs]) for new_fg in new_BBT_used])
                            if not excluded:
                                new_design.n_cycles += 1
                                new_design.bbts.append(BBT.index) # The index of the incoming BBT is incorporated to the design
                                new_design.reactions.append(i) # The reaction used to incorporate the incoming BBT is incorporated to the design
                                new_design.fgs += [new_BBT_used[j] for j in range(len(new_BBT_used)) if  new_BBT_used[j] > 0] # New exposed fgs are incorporated to the design
                                new_design.fg_sources += [3 * (new_design.n_cycles + 1) for j in range(len(new_BBT_used)) if  new_BBT_used[j] > 0] # the sources of the incorporated FGs coming from the BB are inserted
                                if reaction_output[i][0] != 0: # append reaction output first FG to design and corresponding FG source if it is not the null FG
                                    new_design.fgs.append(reaction_output[i][0])
                                    new_design.fg_sources.append(3 * (new_design.n_cycles + 1) - 1)
                                if reaction_output[i][1] != 0: # append reaction output second FG and FG corresponding source to design if it is not the null FG
                                    new_design.fgs.append(reaction_output[i][1])
                                    new_design.fg_sources.append(3 * (new_design.n_cycles + 1) - 1)
                                if len(new_design.fgs) == 0 and len(par.par['max_cycle_na']) > new_design.n_cycles:
                                    excluded = True # this eliminates any design that is closed but only if the cycle is not the last one
                            if not excluded: # this checks for the size of the molecules in this design and excludes it if the smalest molecule is too large already
                                new_design.min_natoms  += BBT.min_atoms
                                excluded = new_design.min_natoms > par.par['max_cycle_na'][new_design.n_cycles - 1] # need to substract 1 because the list ends in n_cycles-1 index
                            if (not excluded) and len(par.par['max_cycle_na']) == new_design.n_cycles: #however when comparing len there is no problem
                                excluded = not all([fg.par[efg]['allowed_end_exposed'] for efg in new_design.fgs]) # if there are any reactive FGs exposed in the last cycle the design is excluded
                            if not excluded:
                                resultado.append(new_design) # the new design is appended to the result list
        return resultado

    def add_deprotections_to_design(self, deprotection_input, deprotection_indexes, available_deprotection_input, deprotection_output, deprotection, par, fg):
        """add_deprotections_to_design will try to deprotect this instance of the Design class by every possible deprotection reaction
        it returns a list of the deprotected Design instances toghether with the original design
        deprotection_input : list of tuples with the incoming FGs for all the reactions
        deprotection_indexes : list of int representing indexes of all the available reactions
        available_deprotection_input : list of tuples with the incoming FGs for only the available reactions
        deprotection_output : list of tuples with the outcoming FGs for all the reactions
        deprotection : instance of Parameters class (deprotection)
        par : instance of Parameters class (par)
        fg : instance of Parameters class (fg)
        returns : resultado (list of Design class instances)"""
        resultado = [self.clone()] # adding the original design to the results list
        resultado[0].deprotections.append(0) # adding the no deprotection reaction to the original design
        resultado[0].dtopology.append(-1) # we update the topology of the design where no deprotection is included
        if self.fgs != []: # if design is closed no further reactions can be made
            for fg_on in self.fgs: # We loop into all the FGs exposed from the initial design
                if (fg_on, 0) in available_deprotection_input: # the FG tuple is in the list of the tuples for incoming FGs for all the available deprotections
                    indices = [i for i in deprotection_indexes if deprotection_input[i][0] == fg_on] # indicates all the possible deprotections that have been found
                    for i in indices: # We loop in all the reactions that have been found
                        new_design = self.clone()
                        fg_on_index = new_design.fgs.index(fg_on) # Here we get the index of the FG that is going to be transformed
                        new_design.dtopology.append(new_design.fg_sources[fg_on_index]) # we update the topology of the new design
                        del new_design.fgs[fg_on_index] # eliminate the FG that is going to be transformed
                        del new_design.fg_sources[fg_on_index] # eliminate the source of the FG being eliminated
                        excluded = any([ex_on in new_design.fgs for ex_on in deprotection.par[i]['excluded_on']]) # becames True if any of the FGs becoming exposed is in the excluded_on list for this reaction
                        if not excluded: # excluded if the transformed GF at position 1 is incompatible with any of the FGs that were already in the design
                            excluded = any([deprotection_output[i][0] in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs])
                        if not excluded: # excluded if the transformed GF at position 2 is incompatible with any of the FGs that were already in the design
                            excluded = any([deprotection_output[i][1] in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs])
                        if not excluded:
                            new_design.deprotections.append(i)
                            if deprotection_output[i][0] != 0:
                                new_design.fgs.append(deprotection_output[i][0]) # Now we add the incoming FG and incorporate also its source
                                new_design.fg_sources.append(3 * (new_design.n_cycles + 1) - 2)
                            if deprotection_output[i][1] != 0: # in most cases this will be 0 but in cases in which deprotection is actually an insertion of a scaffold this could not be the case
                                new_design.fgs.append(deprotection_output[i][1])
                                new_design.fg_sources.append(3 * (new_design.n_cycles + 1) - 2)
                            if deprotection.par[i]['atom_dif'] > 0: # new atoms are added to the design only if a deprotection adds an scaffold, otherwise deprotectionatoms would have been accounted for within the BBTs creation
                                new_design.min_natoms += deprotection.par[i]['atom_dif']
                                 # the design is discarded if is too large for this cycle
                                 # compared with the equivalent line in the addition of the BB, we use new_design.n_cycles instead of new_design.n_cycles -1
                                 # this is because the deprotection (and scaffold incorporation) happens before the BB addition, where new_design.n_cycles is updated
                                excluded = new_design.min_natoms > par.par['max_cycle_na'][new_design.n_cycles]
                            if not excluded:
                                resultado.append(new_design)
        return resultado

    def add_cycle(self, BBTs, indexes,
                  reaction, reaction_indexes, available_reaction_input, reaction_input, reaction_output,
                  deprotection, deprotection_indexes, available_deprotection_input, deprotection_input, deprotection_output,
                  par, fg):
        """add_cycle is a method that will run all possible deprotections (and scaffolds) followed
        by incorporation of compatible BBTs to the design. Both steps are run as a tree, in such
        a way that all designs generated from this one after adding the deprotections will be
        subjected to the incorporation of every compatible building block.
        It returns a new list of designs (instances of Design class)
        BBTs : list of instance of BBT class
        indexes : list of indexes of BBTs that are available (at least one compound has been found in the databases)
        reaction : instance of Parameters class (reactions)
        reaction_indexes : list of int representing indexes of all the available reactions
        available_reaction_input : list of tuples with the incoming FGs for only the available reactions
        reaction_input : list of tuples with the incoming FGs for all the reactions
        reaction_output : list of tuples with the outcoming FGs for all the reactions
        deprotection : instance of Parameters class (deptrotections)
        deprotection_indexes: list of int representing indexes of all the available deprotections
        acailable_deprotection_input : list of tuples with the incoming FGs for only the available deprotections
        deprotection_input : list of tuples with the incoming FGs for all the deprotections
        deptotection_output : list of tuples with the outcoming FGs for all the deprotection
        par : instance of Parameters class (par)
        fg : instance of Parameters class (fg)
        returns : resultado (list of instances of Design class)"""
        resultado = []
        pre_resultado = self.add_deprotections_to_design(deprotection_input, deprotection_indexes, available_deprotection_input, deprotection_output, deprotection, par, fg)
        for design in pre_resultado:
            for i in indexes:
                resultado += design.add_BBT_to_design(BBTs[i], reaction_indexes, available_reaction_input, reaction_input, reaction_output, reaction, par, fg)
        return resultado

    def add_lib_id(self, reaction, deprotection):
        """adds macrodesign ids to this instance of Design by taking all the headpieces reactions and
        deprotections (except final ones) and converting them in a tuple
        reaction : instance of Parameters class (reaction)
        deprotection : instance of Parameters class (deprotection)
        returns : None"""
        lista = [self.n_cycles]
        lista += [self.bbts[0]]
        lista += [reaction.par[self.reactions[i]]['enum_index'] for i in range(1, len(self.reactions))]
        lista += [deprotection.par[self.deprotections[i]]['enum_index'] for i in range(1, len(self.deprotections))]
        lista += self.btopology
        lista += self.dtopology
        self.lib_id = tuple(lista)


# FUNCTION DEFINITIONS

def create_designs(par, BBTs, reaction, deprotection, fg, log, args, path):
    """This is the master function to create an eDESIGN set. It starts generating list of indexes from BBTs and
    reactions objects that will speed up the loops because they have less members than the original objects. Then
    it creates the different cycles storing intermediate eDESIGNS in disk and using parallelization to expand
    each design into the next cycle.
    par : instance of Parameters class (par)
    BBTs : list of instances of BBT class
    reaction : instance of Paramters class (reaction)
    deprotection : instance of Parameters class (deprotection)
    fg : instance of Parameters class (fg)
    args : arguments passed from the main script
    log : instance of Logger class
    path : instance of Parameters class (path)
    returns : None"""

    # Initiallization of indexes
    n_cycles = len(par.par['max_cycle_na'])
    log.update('Creating reaction indexes...')
    # the following variables contain list of tuples containign the pairs or incoming or outcoming FGs to a reaction for all reactions
    reaction_input = [(reaction.par[i]['fg_input_on_off'][0], reaction.par[i]['fg_input_on_off'][1]) for i in range(len(reaction.par))]
    reaction_output = [(reaction.par[i]['fg_output_on_off'][0], reaction.par[i]['fg_output_on_off'][1]) for i in range(len(reaction.par))]
    deprotection_input = [(deprotection.par[i]['fg_input_on_off'][0], deprotection.par[i]['fg_input_on_off'][1]) for i in range(len(deprotection.par))]
    deprotection_output = [(deprotection.par[i]['fg_output_on_off'][0], deprotection.par[i]['fg_output_on_off'][1]) for i in range(len(deprotection.par))]
    # the following calculates the list of indexes for available reactions
    reaction_indexes = [i for i in range(1, len(reaction.par)) if (par.par['include_designs'].upper() == 'BOTH' or reaction.par[i]['production'])]
    available_reaction_input = [reaction_input[i] for i in reaction_indexes]
    deprotection_indexes = [i for i in range(1, len(deprotection.par)) if (par.par['include_designs'].upper() == 'BOTH' or deprotection.par[i]['production'])]
    available_deprotection_input = [deprotection_input[i] for i in deprotection_indexes]
    # the following calculate lists of indexes related with BBTs
    hp_indexes = [i for i in range(len(BBTs)) if BBTs[i].headpiece is not None] # this extracts the indexes of the BBTs that are headpieces
    indexes = [i for i in range(len(BBTs)) if BBTs[i].n_compounds[-1] > 0] # this extracts the indexes of the BBTs that contain at least one compound

    # the next line packages all the parameters required to expand a cycle into a tuple
    parallel_args = (BBTs, indexes, reaction, reaction_indexes, available_reaction_input, reaction_input, reaction_output,
                     deprotection, deprotection_indexes, available_deprotection_input, deprotection_input, deprotection_output,
                     par, fg)

    # create designs
    log.update('Creating designs...')
    # Create the first set of designs containing only the available headpieces
    log.update('Creating headpieces...')
    designs = [Design(par, BBTs, index) for index in hp_indexes] # creates the first level designs containing just the headpieces
    log.update(f'        this round number designs = {len(designs)}')
    log.update(f'        cycle 0 number designs = {len(designs)}')
    log.update('        Saving designs to file...')
    with open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_' + path.par['run'] + '_0.pic'), 'wb') as f:
        for design in designs:
            pic.dump(design, f)

    lib_id_list = []  # this list will hold a unique copy of the lib_design_ids
    # Create designs
    for cycle in range(n_cycles):
        log.update(f'Creating cycle {cycle + 1}...')
        f = open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_' + path.par['run'] + '_' + str(cycle) + '.pic'), 'rb')
        g = open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_' + path.par['run'] + '_' + str(cycle + 1) + '.pic'), 'wb')
        continuar = True
        counter_read = 0
        counter = 0
        last_counter = 0
        while continuar:
            previous_designs = []
            log.update(f'            Reading from disk from design number {counter_read}...')
            for _ in range(par.par['designs_in_memory']):
                try:
                    previous_designs.append(pic.load(f))
                except:  # We have reached the end of the file
                    continuar = False
                    break
                else:
                    counter_read +=1
            else:  # We have not reached the end of the file
                continuar = True
            log.update(f'            Processing {len(previous_designs)} designs...')
            if len(previous_designs) == 0:
                designs = []
            else:
                tic = time.time()
                designs = starmap_parallel(previous_designs, expand_designs, args=parallel_args, is_list=True)
                tac = time.time()
                log.update(f'            Time to expand designs: {round(tac - tic, 0)} seconds')

            counter += len(designs)
            if cycle == n_cycles - 1:
                log.update('            Removing desings with not enough size...')
                designs = [design for design in designs if design.min_natoms < par.par['max_na_absolute']]
                log.update('            Adding lib_id to designs...')
                for i in range(len(designs)):
                    designs[i].add_lib_id(reaction, deprotection)
                log.update('            Maintaining lib_id dictionary...')
                lib_id_list += [design.lib_id for design in designs]
                lib_id_list = list(set(lib_id_list))
                log.update('            Adding desing_id to desings...')
                for i in range(len(designs)):
                    designs[i].id = i + last_counter
                last_counter += len(designs)
            log.update('            Writting to disk ' + str(len(designs)) + ' designs...')
            for design in designs:
                pic.dump(design, g)
        f.close()
        g.close()
        log.update(f'        cycle {cycle + 1} number designs = {counter}')
    log.update(f'Numebr of designs after removing low density desings: {last_counter}')
    log.update(f'Number of library ids generated: {len(lib_id_list)}')
    log.update('Saving lib_id_list to file...')
    with open(os.path.join(args.wfolder, 'data', path.par['Database_Run'] + '_' + path.par['run'] + '_lib_id_list.pic'), 'wb') as f:
        pic.dump(lib_id_list, f)
    return None


def expand_designs(previous_designs, BBTs, indexes, reaction, reaction_indexes, available_reaction_input,
                   reaction_input, reaction_output, deprotection, deprotection_indexes, available_deprotection_input,
                   deprotection_input, deprotection_output, par, fg):
    """This function generates an expansion of a list of designs with additional arguments:
    previous_designs : list of instances of Design class
    BBTs : list of instances of BBT class
    indexes : list of int containing
    reaction : instance of Paramters class (reaction)
    reaction_indexes : list of indexes for all available reactions
    available_reaction_input : list of tuples containing the input FGs for all available reactions
    reaction_input : list of tuples containing the input FGs for all reactions
    reaction_output : list of tuples containing the output FGs for all reactions
    deprotection : instance of Parameters class (deprotection)
    deprotection_indexes : list of indexes for all available deprotections
    available_deprotection_input : ist of tuples containing the output FGs for all available deprotections
    deprotection_input : list of tuples containing the output FGs for all deprotections
    deprotection_output : list of tuples containing the output FGs for all deprotections
    par : instance of Parameters class (par)
    fg : instance of Parameters class (fg)
    returns result : list of instances of Design class"""
    result = []
    for design in previous_designs:
        result += design.add_cycle(BBTs, indexes,  reaction, reaction_indexes, available_reaction_input, reaction_input,
                                   reaction_output, deprotection, deprotection_indexes, available_deprotection_input,
                                   deprotection_input, deprotection_output, par, fg)
    return result


if __name__ == '__main__':
    print(version)
