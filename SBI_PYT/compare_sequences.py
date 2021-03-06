# Adding module location to sys.path
import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages')
sys.path.append('/Users/Maria/Desktop/Uni/sbi/SBI_PYT/PPI_main') # Path to I_O_args script
sys.path.append('/home/oth/BHS/PYT/1.project/dash_test.py')
import os
import numpy
import scipy
import I_O_args as args
from Bio.PDB import *
from Bio.Seq import Seq
from Bio import pairwise2 as pw2
import dash_test



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
    structure = parser.get_structure(f"TYPE_{i}", pairs)
    model = structure[0]
    dictmodels[final_name] = model # Path+file as KEY, model object as VALUE
    i += 1


####################
#   FUNCTIONS      #
####################
similar_chains = {}
pairs_list=[]
for id, pairs in dictmodels.items():
    p=[]
    for chain in pairs:
        p.append(chain)
    pairs_list.append(p)


#Compare within the pair
for single_pair in pairs_list:
    ppb = PPBuilder()
    for pp1 in ppb.build_peptides(sigle_pair[0]):
        seq1=pp1.get_sequence()
    for pp2 in ppb.build_peptides(sigle_pair[1]):
        seq2=pp2.get_sequence()
        alignments = pw2.align.globalxx(seq, seq2, score_only=True)
                            #print(pw2.format_alignment(*alignments[0]))
        identity_perc = round(alignments/len(seq),2)
        if identity_perc > 0.95:
            fixed=single_pair[0] #first seq of the pair will be the fixed
            moving=single_pair[1] #second seq of the pair will be the moving
            fixed_atoms = get_atoms_list(fixed)
            moving_atoms = get_atoms_list(moving)
            fixed_atoms = get_atoms_list(fixed)
            moving_atoms = get_atoms_list(moving)

        #Superimposition
            super_imposer = Superimposer()
            super_imposer.set_atoms(fixed_atoms, moving_atoms)
            super_imposer.apply(model.get_atoms())
            print (super_imposer.rms)
