
##############################
# SBI / PYT PROJECT EXERCISE #
##############################

import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages') # Add anaconda path

import argparse
import os
import re
import numpy

# 1. INPUT / OUTPUT INTERFACE
##############################


# Functions to check input and output directories, and filenames
# Lo podríamos meter en una clase, todo esto

def dir_path(string):
    """
    This function takes as an input a string indicating the PDB files source directory. It checks
    that it is, in fact, a directory.
    """
    if os.path.isdir(string):
        sys.stderr.write("\tInput path is correctly set.\n")
        sys.stderr.flush()
        return string
    else:
        raise NotADirectoryError(string)

def create_outdir(string, forcing = False, exiting = False):
    """
    This functions checks if the output directory indicated in the command-line exists. If It does not,
    it creates one. If it does, it is reported to the user. Returns the ouptut path as a string.
    """
    path = string
    if forcing == False:
        if not os.path.exists(path):
            os.mkdir(string)
            sys.stderr.write("\tCreating output directory...\n")
            sys.stderr.flush()
        else:
            sys.stderr.write("\tOutput directory already exists\n")
            sys.stderr.flush()
            if exiting == True:
                sys.stderr.write("\tExiting the program...\n")
                sys.stderr.flush()
                return exit()
        return path

    if forcing == True:
        if os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            sys.stderr.write("\t-f flag used: The output directory will be overwritten.\n")
            sys.stderr.flush()

def check_infile_names(string):
    """
    This function checks that the input directory .pdb files are named correctly. It returns
    a list with the correct filenames.
    """
    pattern = r"^\w+\_\w\_\w\.pdb+(\.gz)?"
    matches = [str(string)+"/"+f for f in os.listdir(string) if re.match(pattern, f) is not None]
    sys.stderr.write("\tPDB input files named correctly!\n\n")
    sys.stderr.flush()
    return matches

# Argparse setup
#
def output_argparser():
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
    args = parser.parse_args()

    # Forcing setting
    #
    if args.forcing:
        create_outdir(args.outdir, forcing=True)
    else:
        create_outdir(args.outdir, exiting=True) # If -f flag not used and exiting is stated

    return args # returns the argparse object


options = output_argparser() # Object that contains all argument parsing functions
a = check_infile_names(options.infile) # Checking correct file-naming


