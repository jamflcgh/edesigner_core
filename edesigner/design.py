# -*- coding: utf-8 -*-
# design
# Jose Alfredo Martin

# Python modules
import copy


__version__ = 'design.v.9.0.0'
__author__ = 'Alfredo Martin'

# CLASS DEFINITION


class Design:
    """Design instances store attributes of an eDESIGN and methods and attributes for its creation
    and modification. The most important methods are the 'add_BBT_to_designs' that creates and returns
    a list of Design instances based on the addition of a BBT instance to the current Design instance,
    and the 'add_deprotection_to_design' method that creates a list of Design instances upon the application
    of a deprotection / scaffold addition reaction.
    """
    def __init__(self, par, BBTs, hp_index, total_cycles):
        """This method initializes the Design class
        par : instance of Parameters class (parameters)
        BBTs : list of instances of the BBT class
        hp_index : int (index of the BBT in BBTs used as headpiece)
        total_cycles: how many cycles will have this design"""
        self.id = None
        self.lib_id = None
        self.total_cycles = total_cycles
        self.n_cycles = 0
        self.bbts = [hp_index]
        self.btopology = []
        self.dtopology = []
        self.reactions = []
        self.deprotections = []
        self.n_deprotections = 0
        self.n_unpr_deprotections = 0
        self.fgs = [BBTs[hp_index].BBT[i] for i in range(3) if BBTs[hp_index].BBT[i] > 0]
        self.fg_sources = [0 for _ in self.fgs]
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
        if len(self.fgs) > 0:  # if design is closed no further reactions can be made
            # We loop into all the non Null FGs of the incoming BBT
            for fg_off in [bbt_fg for bbt_fg in BBT.BBT if bbt_fg != 0]:
                # We loop into all the FGs exposed from the initial design
                for fg_on in self.fgs:
                    # the FG tuple is in the list of the tuples for incoming FGs for all the available reactions
                    if (fg_on, fg_off) in available_reaction_input:
                        # indicates all the reactions that have been found within the available ones
                        indices = [i for i in reaction_indexes if (reaction_input[i] == (fg_on, fg_off))]
                        # We loop in all the reactions that have been found
                        for i in indices:
                            new_design = self.clone()
                            new_design.n_cycles += 1
                            fg_on_index = new_design.fgs.index(fg_on)
                            # Check whether the new reaction would act in a FG that comes from de-protection
                            if (new_design.fg_sources[fg_on_index] - 1) % 3 == 0:  # The source is a de-protection
                                new_design.n_unpr_deprotections -= 1  # we deactivate an unproductive de-protection
                            # now we exclude the design if it will end with at least one unproductive de-protection
                            excluded = new_design.total_cycles - new_design.n_cycles - new_design.n_unpr_deprotections < 0
                            if excluded:
                                continue
                            # we update the topology of the new design
                            new_design.btopology.append(new_design.fg_sources[fg_on_index])
                            # In the new design the reacting FG is removed
                            del new_design.fgs[fg_on_index]
                            # In the new design the source corresponding to the removed fg is removed
                            del new_design.fg_sources[fg_on_index]
                            # We make a copy of the incoming BBT FG list
                            new_BBT_used = copy.deepcopy(BBT.BBT)
                            # and then we remove the FG that is reacting
                            new_BBT_used.remove(fg_off)
                            # excluded if any of the FGs becoming exposed is in the excluded_on list for this reaction
                            excluded = any([on_ex in new_design.fgs for on_ex in reaction.par[i]['excluded_on']])
                            if excluded:
                                continue
                            # excluded if any of the fgs in reactions excluded_on for this reaction is in the
                            # remaining FGs of the incoming BBT
                            excluded = any([off_ex in new_BBT_used for off_ex in reaction.par[i]['excluded_off']])
                            if excluded:
                                continue
                            # excluded if any of the resulting fgs for this reaction is in the self incompatibility
                            # list of any of the fgs that were already on dna except the one that reacted
                            excluded = any([any([new_fg in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs]) for new_fg in reaction_output[i]])
                            if excluded:
                                continue
                            # excluded if any of the fgs in the incoming BBT, except the one that reacts is in the self
                            # incompatibility list of any of the fgs that were already on dna
                            # except the one that reacted
                            excluded = any([any([new_fg in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs]) for new_fg in new_BBT_used])
                            if excluded:
                                continue
                            # The index of the incoming BBT is incorporated to the design
                            new_design.bbts.append(BBT.index)
                            # The reaction used to incorporate the incoming BBT is incorporated to the design
                            new_design.reactions.append(i)
                            # New exposed fgs are incorporated to the design
                            new_design.fgs += [item for item in new_BBT_used if item > 0]
                            # the sources of the incorporated FGs coming from the BB are inserted
                            # fg_sources must be new_design.n_cycles * 3 because we have increased already the
                            # n_cycles of the new design
                            new_design.fg_sources += [3 * new_design.n_cycles for item in new_BBT_used if item > 0]
                            # Now we append the FGs that are an output of the reaction if any. Also we set their sources
                            # Currently the both output FGs have the source of the same reaction, which is the incoming
                            # cycle index
                            if reaction_output[i][0] != 0:
                                new_design.fgs.append(reaction_output[i][0])
                                new_design.fg_sources.append(3 * new_design.n_cycles - 1)
                            if reaction_output[i][1] != 0:
                                new_design.fgs.append(reaction_output[i][1])
                                new_design.fg_sources.append(3 * new_design.n_cycles - 1)
                            # Now we exclude the design if it is closed and in the last cycle
                            excluded = len(new_design.fgs) == 0 and new_design.total_cycles - new_design.n_cycles > 0
                            if excluded:
                                continue
                            # this checks for the size of the molecules in this design and excludes it if the smalest
                            # molecule is too large already
                            new_design.min_natoms += BBT.min_atoms
                            # need to substract 1 because the list ends in n_cycles-1 index
                            excluded = new_design.min_natoms > par.par['max_cycle_na'][new_design.n_cycles - 1]
                            if excluded:
                                continue

                            # Exclude if there are any reactive FGs exposed in the last cycle the design is excluded
                            excluded = new_design.total_cycles == new_design.n_cycles and not all([fg.par[efg]['allowed_end_exposed'] for efg in new_design.fgs])
                            if excluded:
                                continue
                            resultado.append(new_design)
        return resultado

    def add_deprotections_to_design(self, deprotection_input, deprotection_indexes, available_deprotection_input, deprotection_output, deprotection, par, fg):
        """add_deprotections_to_design will try to deprotect this instance of the Design class by every possible deprotection reaction
        it returns a list of the deprotected Design instances together with the original design
        deprotection_input : list of tuples with the incoming FGs for all the reactions
        deprotection_indexes : list of int representing indexes of all the available reactions
        available_deprotection_input : list of tuples with the incoming FGs for only the available reactions
        deprotection_output : list of tuples with the outcoming FGs for all the reactions
        deprotection : instance of Parameters class (deprotection)
        par : instance of Parameters class (par)
        fg : instance of Parameters class (fg)
        returns : resultado (list of Design class instances)"""
        resultado = [self.clone()]  # adding the original design to the results list
        resultado[0].deprotections.append(0)  # adding the no de-protection reaction to the original design
        resultado[0].dtopology.append(0)  # we update the topology of the design where no de-protection is included
        if self.fgs != []:  # if design is closed no further reactions can be made
            for fg_on in self.fgs:  # We loop into all the FGs exposed from the initial design
                if (fg_on, 0) in available_deprotection_input:  # the FG tuple is in the list of the tuples for incoming FGs for all the available deprotections
                    indices = [i for i in deprotection_indexes if deprotection_input[i][0] == fg_on]  # indicates all the possible deprotections that have been found
                    for i in indices:  # We loop in all the de-protections that have been found
                        new_design = self.clone()
                        new_design.n_deprotections += 1
                        new_design.n_unpr_deprotections += 1
                        # now we exclude the design if it will end with at least one unproductive de-protection
                        excluded = new_design.total_cycles - new_design.n_cycles - new_design.n_unpr_deprotections < 0
                        if excluded:
                            continue

                        # Next we get the index of the FG that is going to be transformed
                        fg_on_index = new_design.fgs.index(fg_on)
                        # Next we update the d topology of the new design with the source of the FG it operates on
                        new_design.dtopology.append(new_design.fg_sources[fg_on_index])
                        del new_design.fgs[fg_on_index]  # eliminate the FG that is going to be transformed
                        del new_design.fg_sources[fg_on_index]  # eliminate the source of the FG being eliminated
                        # Eliminated if any of the FGs becoming exposed is in the excluded_on list for this de-protection
                        excluded = any([ex_on in new_design.fgs for ex_on in deprotection.par[i]['excluded_on']])
                        if excluded:
                            continue
                        # excluded if the transformed GF at position 1 is incompatible with any of the FGs that were already in the design
                        excluded = any([deprotection_output[i][0] in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs])
                        if excluded:
                            continue
                        # excluded if the transformed GF at position 2 is incompatible with any of the FGs that were already in the design
                        excluded = any([deprotection_output[i][1] in fg.par[old_fg]['self_incompatibility'] for old_fg in new_design.fgs])
                        if excluded:
                            continue
                        new_design.deprotections.append(i)
                        # Now we add the incoming FG(s) and incorporate also their source, which is this deprotection
                        if deprotection_output[i][0] != 0:
                            new_design.fgs.append(deprotection_output[i][0])
                            new_design.fg_sources.append(1 + 3 * new_design.n_cycles)
                        if deprotection_output[i][1] != 0:
                            new_design.fgs.append(deprotection_output[i][1])
                            new_design.fg_sources.append(1 + 3 * new_design.n_cycles)
                        # new atoms are added to the design only if a de-protection adds an scaffold, otherwise
                        # de-protection atoms would have been accounted for within the BBTs creation
                        if deprotection.par[i]['atom_dif'] > 0:
                            # the design is discarded if is too large for this cycle
                            # compared with the equivalent line in the addition of the BB, we use
                            # new_design.n_cycles instead of new_design.n_cycles -1
                            # this is because the deprotection (and scaffold incorporation) happens before the
                            # BB addition, where new_design.n_cycles is updated
                            new_design.min_natoms += deprotection.par[i]['atom_dif']
                            excluded = new_design.min_natoms > par.par['max_cycle_na'][new_design.n_cycles]
                        if excluded:
                            continue
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
        pre_resultado = self.add_deprotections_to_design(deprotection_input, deprotection_indexes,
                                                         available_deprotection_input, deprotection_output,
                                                         deprotection, par, fg)
        for design in pre_resultado:
            for i in indexes:
                resultado += design.add_BBT_to_design(BBTs[i], reaction_indexes, available_reaction_input,
                                                      reaction_input, reaction_output, reaction, par, fg)
        return resultado

    def add_lib_id(self, reaction, deprotection):
        """adds macrodesign ids to this instance of Design by taking all the headpieces reactions and
        deprotections (except final ones) and converting them in a tuple
        reaction : instance of Parameters class (reaction)
        deprotection : instance of Parameters class (deprotection)
        returns : None"""
        lista = [self.total_cycles]
        lista += [deprotection.par[i]['enum_index'] for i in self.deprotections]
        lista += [reaction.par[i]['enum_index'] for i in self.reactions]
        lista += self.dtopology
        lista += self.btopology
        lista += [self.bbts[0]]
        self.lib_id = tuple(lista)


if __name__ == '__main__':
    print(__version__)
    print(__author__)
