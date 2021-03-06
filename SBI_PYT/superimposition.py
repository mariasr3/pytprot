import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages')
sys.path.append('/Users/Maria/Desktop/Uni/sbi/SBI-Project/SBI_PYT/PPI_main')
import os
import numpy
import scipy
import I_O_args as args
from Bio.PDB import *
from Bio.Seq import Seq
from Bio import pairwise2 as pw2
import dash_test
import test_project



# Argument parsing functions
options = args.output_argparser() # Object that contains all argument parsing functions
a = args.check_infile_names(options.infile) # Checking correct file-naming

# Using 1gzx as an input and always -f flag active


# List of files
pdblist = []
for pdbfile in os.listdir(options.infile):
    pdblist.append(str(options.infile+pdbfile)) # List of PDBs with their path


# Dictionary with files and their PDB models
i = 0
dictmodels = {}
for pairs in pdblist:
    #final_name = pairs.split("/")[-1] # Better to not split it, actually
    final_name = pairs
    dictmodels[final_name]= []
    parser = PDBParser(PERMISSIVE=True, QUIET=False)
    structure = parser.get_structure(f"TEST_{i}", pairs)
    model = structure[0]
    dictmodels[final_name] = model # Path+file as KEY, model object as VALUE
    i += 1

#print(test_project.alignment(dictmodels))


def get_atoms_list(chain):
    """Creates a list of the atoms only taking CA or P for protein and acid nucleics, respectively.
    This list of atoms will be lately used in the superimposition process."""
    atom_id = "CA"
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        atoms_list.append(atom)
    return atoms_list



pairs_list=[]
for id, pairs in dictmodels.items():
    p=[]
    for chain in pairs:
        p.append(chain)
    pairs_list.append(p)


for single_pair in pairs_list:
    fixed=single_pair[0]
    moving=single_pair[1]
    fixed_atoms = get_atoms_list(fixed)
    moving_atoms = get_atoms_list(moving)

    super_imposer = Superimposer()
    super_imposer.set_atoms(fixed_atoms, moving_atoms)
    super_imposer.apply(model.get_atoms())
    print (super_imposer.rms)


"""
#prova amb cadenes = length
fixed=pairs_list[1][0]
moving=pairs_list[1][1]
fixed_atoms = get_atoms_list(fixed)
moving_atoms = get_atoms_list(moving)

super_imposer = Superimposer()
super_imposer.set_atoms(fixed_atoms, moving_atoms)
super_imposer.apply(model.get_atoms())
print (super_imposer.rms)
"""
