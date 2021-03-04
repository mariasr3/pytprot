
# Adding module location to sys.path
import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages')
sys.path.append('/home/oth/BHS/PYT/1.project/SBI_PYT/PPI_main') # Path to I_O_args script
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
    structure = parser.get_structure(f"TEST_{i}", pairs)
    model = structure[0]
    dictmodels[final_name] = model # Path+file as KEY, model object as VALUE
    i += 1


####################
#   FUNCTIONS      #
####################
def alignment(input_pdb_dict):
    """
    Given an dictionary of PDB files as Keys, and their Structure objects as values, this function
    makes pairwise alignments between the chains, keepin only those with a 95% or higher
    similarity.
    ¿¿¿¿¿ Vale la pena añadir los parámetros del alignment?????

    PBuilder need as an input structure object, we use
    dictmodels where the key is the filename and the value is
    the model object obtained from PDBParser.
    The output is the sequence in 1 letter code.
    Pairwsie2 given the "parameter" globalsxx it finds the
    best global aligntment between two seqquences. Identical characters are
    given 1 point. No points are deductes for mismatches or gaps"
    """
    similar_chains = {}
    for id, model in input_pdb_dict.items():
        ppb = PPBuilder()
        for chain in model:
            for pp in ppb.build_peptides(chain):
                seq = pp.get_sequence()
                for id2, model2 in input_pdb_dict.items():
                    for chain2 in model2:
                        for pp in ppb.build_peptides(chain2):
                            seq2 = pp.get_sequence()
                            alignments = pw2.align.globalxx(seq, seq2, score_only=True)
                            #print(pw2.format_alignment(*alignments[0]))
                            identity_perc = round(alignments/len(seq),2)
                            if identity_perc > 0.95:
                                #print(f"This is the alignment between {chain.id}, {chain2.id}: ")
                                #print(identity_perc)
                                #similar_chains.setdefault(chain.id, []).append(chain2.id)
                                try:
                                    similar_chains[chain.id] = chain2.id
                                except KeyError:
                                    similar_chains[chain.id].append(chain2.id)
    return similar_chains

# Alignment funciton test:
for x,y in alignment(dictmodels).items():
    print(x, y)



# IMPROVEMENTS: Crear un archivo con los alineamientos. O a lo mejor ocupa mucha memoria.


def check_type(input_pdb_file):
    """
    Given a PDB file, it checks if it pertains to a DNA or protein structure.
    :param input_pdb_file:
    :return: structure object
    """
    i = 0
    final_name = input_pdb_file.split("/")[-1] # Only keep name, not full path
    structure = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(f"model_{i}", input_pdb_file)
    atoms = structure[0].get_atoms()
    for x in atoms:
        if x.get_id() == "CA":
            print(f"{final_name} is a PROTEIN file") # Print optional
            return structure
        if x.get_id() == "P":
            print(f"{final_name} is a DNA file") # Print optional
            return structure

    # IMPROVEMENT: Maybe compute the sequences within this function?? Porque de momento
    # sólo comprueba si es prot y dna

def change_chain_id(input_pdb_files):
    """
    Given a list of .pdb files, the Chain IDs are converted from a letter format to a numeric format
    :param input_pdb_file: 
    :return: A list of Structures
    """
    struct_list = []
    i = 1
    j = 0
    while i < len(input_pdb_files)*2: # As we will always have 2 chains per file
        for file in input_pdb_files:
            structure = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(f"model_{j}", file)
            chains = structure[0].get_chains()
            j += 1
            for ch in chains:
                ch.id = i # replace chain ID with number
                i += 1 # i will keep increasing depending on the number of models
            struct_list.append(structure)
            break

    # IMPROVEMENT: Cambiar la letra por su posición en el abecedario? De esa manera
    # si tenemos varias "A" todas se cambiarían por 1, pero no sé si esto es
    # realmente necesario o no

    return struct_list



# Chain ID change test
#
"""for x in change_chain_id(pdblist):
    print(x.get_id())
    for ch in x.get_chains():
        print(ch.get_id())"""

"""for x,y in dictmodels.items():
    print(x)
    for ch in y.get_chains():
        print(ch)"""

# FOR THE DASH APP TEST #

# Now it automatically opens the browser

#if __name__ == '__main__':
#    dash_test.Timer(1, dash_test.open_browser).start();
#    dash_test.app.run_server(debug=True, port=dash_test.port) # Hot reloading







