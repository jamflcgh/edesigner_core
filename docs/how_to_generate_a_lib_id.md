# How to generate a lib_id

lib_id is the name used throughout the different scripts in eDESINGER to refer to a topology token that uniquely identifies an eDESIGNER library.
As eDESIGNER generates libraries it computes their corresponding lib_id but a lib_id can be created manually for any custom library that uses reactions coded in eDESIGNER parameters.
This token can be used to enumerate any custom library using the eDESIGNER enumeration engine. 
This file explais step by step how to create a lib_id for any custom library.

In order to create a lib_id you have to use the information in the parameters enum_deprotection.par, enum_reaction.par and headpieces.par which are located in the folder resouces under edesigner folder

The lib_id is a tuple of integer numbers than when passed to scripts, is passed in a form of a string with the numbers separated by underscores, for example 1_0_0_0_1_2_3_0_0_0_0_1_2_3

In the following, cycles are numbered as 0 (headpiece) and then 1, 2, 3 ...

The topology is embeded in the lib_id vector that has the structure:

    n_cycles: int (number of cycles)
    
    d0, d1, d2...: (deprotection enumeration indexes at each cycle from 0 to n_cycles -1)
    
    r1, r2, r3...: (reaction enumeration indexes at the incorporation of each cycle from 1 to n_cycles)
    
    sd0, sd1, sd2...: (source each deprotection from d0 to d[n_cycles-1]). Sources are defined by the origin of the FG in which deprotection is applied. It can be a cycle or a reaction.
    
    sc1, sc2, sc2...: (source of each cycle or reaction, which is the same, from r1 to r[n_cycles]. Sources are defined by the origin of the FG in which deprotection is applied. It can be a cycle, a deprotection or a reaction.
    
    headpiece: int: index of the bbt coding the headpiece. It is extracted from the hp parameters as the last number (the index of the last FG in the three number tuple describing the headpiece)

    In summary, the structure of the lib_id for a 2 cycle libary would be:

    2_d0_d1_r1_r2_sd0_sd1_sc1_sc2_hp

    and for a 3 cycle library would be:

    3_d0_d1_d2_r1_r2_r3_sd0_sd1_sd2_sc1_sc2_sc3_hp

    From this, the least straightforward to derive are the sources (both for deprotections and reactions). 
    As a rule of thumb you have to look to the functional group that was used in a deprotection or reaction.
    If the functional group used by a reaction (the first of the two used) or a deprotection came from a building block then the cycle corresponding to the building block incorporation is the source.
    If the functional group comes from a deprotection then that deprotection is the source.
    If a functional group was generated in a reaction then that reaction is the source.
    
    The posible values for the sources, which are the same for deprotections and cycles (reactions) are:
    
    0: c0 (Headpiece)
    
    1: d0 (deprotection conduced to the headpiece before the first reaction, 0 if no deprotection conducted)
    
    2: r1 (reaction that joined headpeice and cycle 1 building block)
    
    3: c1 (functional group was incorporated with building block in cycle 1)
    
    4: d1 (deprotection conducted before incorporation of cycle 2 but after reaction that incorporated cycle 1, 0 if no deprotection conducted)
    
    5: r2 (reaction that joined the cycle 2 with the intemrediate at that stage)
    
    6: c2 (functional group incorporated with building block in cycle 2)
    
    7: d2 (deprotection conducted before incorporation of cycle 3 but after reaction that incorporated cycle 2, 0 if no deprotection conducted)
    
    8: r3 (reaction that joined the cycle 3 with the intemrediate at that stage)
    
    9: c3 (functional group incorporated with building block in cycle 2)
    
    ...
