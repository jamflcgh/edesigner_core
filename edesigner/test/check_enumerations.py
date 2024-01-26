# -*- coding: utf-8 -*-
# check_enumerations.v.12.0.0
# Jose Alfredo Martin 2023

__version__ = 'check_enumerations.v.12.0.0'
__author__ = 'Alfredo Martin 20203'

# Python modules
import os
import argparse
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="""Checks for quality of enumerations by ensuring that all the 
    compounds have been enumerated
        """)
    parser.add_argument('-eF', '--enumeration_folder',
                        help="""folder where the enumerations are stored""",
                        type=str,
                        required=True)
    parser.add_argument('-oF', '--out_folder',
                        help="""folder where the output files will be stored""",
                        type=str,
                        required=True)
    parser.add_argument('-v', '--verbose',
                        help="""When invoked additional information is printed in the console.""",
                        action='store_true')
    args = parser.parse_args()
    assert os.path.isdir(args.enumeration_folder), f'{args.enumeration_folder} does not exist'
    assert os.path.isdir(args.out_folder), f'{args.out_folder} does not exist'
    return args


def main():
    args = parse_args()
    dirlist = os.listdir(args.enumeration_folder)
    failed = []
    incomplete = []
    for folder in tqdm(dirlist):
        file_list = os.listdir(os.path.join(args.enumeration_folder, folder))
        if not 'enumeration.smi' in file_list:
            failed.append(folder)
            continue
        file_list = [os.path.join(args.enumeration_folder, folder, file) for file in file_list if file.startswith('R0')]
        ncomps = 1
        for file in file_list:
            with open(file, 'r') as f:
                ncomps *= len(f.readlines())
        with open(os.path.join(args.enumeration_folder, folder, 'enumeration.smi'), 'r') as f:
            enum_comps = len(f.readlines())
        if ncomps != enum_comps:
            incomplete.append(folder)
    with open(os.path.join(args.out_folder, 'enum_failed'), 'w') as f:
        for line in failed:
            f.write(line + '\n')
    with open(os.path.join(args.out_folder, 'enum_incomplete'), 'w') as f:
        for line in incomplete:
            f.write(line + '\n')


if __name__ == '__main__':
    print(__version__)
    print(__author__)
    main()
