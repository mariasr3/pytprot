
##############################
# SBI / PYT PROJECT EXERCISE #
##############################

import argparse
import os
import sys
import re

# 1. INPUT / OUTPUT INTERFACE
##############################


# Functions to check input and output directories, and filenames
# Lo podríamos meter en una clase

def dir_path(string):
    """
    This function takes as an input a string indicating the PDB files source directory. It checks
    that it is, in fact, a directory.
    """
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

def create_outdir(string):
    """
    This functions checks if the output directory indicated in the command-line exists. If It does not,
    it creates one. If it does, it is reported to the user. Returns the ouptut path as a string.
    """
    if not os.path.exists(string):
        os.mkdir(string)
    else:
        print("This output directory already exists...")
    return string

def check_infile_names(string):
    """
    This function checks that the input directory .pdb files are named correctly. It returns
    a list with the correct filenames.
    """
    pattern = r"^\w+\_\w\_\w\.pdb+(\.gz)?"
    matches = [str(string)+"/"+f for f in os.listdir(string) if re.match(pattern, f) is not None]
    return matches


# Argparse setup
#
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
                    dest="stechiometry",
                    action="store",
                    type=str,
                    default=None,
                    help="Chain stechiometry file.")

parser.add_argument('-f', '--force',
                    dest="forcing",
                    action="store_true",
                    default=False,
                    help="Forces creation of output directory.")

parser.add_argument('-o', '--output-directory',
                    dest="outdir",
                    required=True,
                    action="store",
                    type=create_outdir,
                    default=None,
                    help="Output file.")

parser.add_argument('-v', '--verbose',
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Print log in stderr.")

# Store parser arguments in options object
#
options = parser.parse_args()

# Forcing setting
#
if options.forcing:
    create_outdir(options.outdir)
    print("-f flag used. The output directory has been overwritten.")


# An example
infiles = check_infile_names(options.infile)
for file in infiles:
    print(file)

