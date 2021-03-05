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

superimpose=[]
superimp={}
for id, pairs in dictmodels.items():
    superimp={}
    for chain in pairs:
        id=chain.id
        ppb=PPBuilder().build_peptides(chain)
        seq=ppb[0].get_sequence()
        superimp[id]=seq
    superimpose.append(superimp)


fixed_atoms=[]
moving_atoms=[]
for pair in superimpose:
    for id1, fixed in pair.items():
        print(id1, fixed)

        #for ref_res, alt_res, amino, allow in zip(fixed, moving, use) :
        #    assert ref_res.resname == alt_res.resname
        #    assert ref_res.id      == alt_res.id
        #    if allow :
                #CA = alpha carbon
        #        ref_atoms.append(ref_res['CA'])
        #        alt_atoms.append(alt_res['CA'])
