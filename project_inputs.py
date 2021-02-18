

# SBI / PYT PROJECT EXERCISE

# 1. INPUT / OUTPUT INTERFACE

# SUGGESTIONS
# Check correct filename with REGEX?

import argparse
import os
import sys

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def createdir(dirpath):
    if options.forcing:
        try:
            return dirpath
        except NotADirectoryError:
            os.mkdir(dirpath)



parser = argparse.ArgumentParser(description="This program does something",
                                 prog="SBI/PYT project",
                                 usage="Calculating PPIs")

parser.add_argument('-i', '--input-directory',
                    dest="infile",
                    required=True,
                    action="store",
                    default=None,
                    type=dir_path,
                    help="Input set of PDB structure files.") # The files must follow a specific structure

parser.add_argument('-s', '--stechiometry',
                    dest="stech",
                    action="store",
                    type=dir_path,
                    default=None,
                    help="Chain stechiometry file.")

parser.add_argument('-f', '--force',
                    dest="forcing",
                    action="store",
                    type=dir_path,
                    default=False,
                    help="Chain stechiometry file.")

parser.add_argument('-o', '--output-directory',
                    dest="outdir",
                    required=True,
                    action="store",
                    type=dir_path,
                    default=None,
                    help="Output file.")

parser.add_argument('-v', '--verbose',
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Print log in stderr.")

try:
    options = parser.parse_args()
except NotADirectoryError:
    createdir(sys.argv[1])

