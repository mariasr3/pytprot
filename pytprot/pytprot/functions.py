
# Adding module location to sys.path
import sys
sys.path.append('/home/oth/anaconda3/lib/python3.8/site-packages')
sys.path.append('/Users/Maria/Desktop/Uni/sbi/SBI-Project/SBI_PYT/PPI_main') # Path to I_O_args script
sys.path.append('/home/oth/BHS/PYT/1.project/SBI_PYT/PPI_main') # Path to I_O_args script
sys.path.append('/home/oth/BHS/PYT/1.project/dash_test.py')
import os
import numpy
import scipy
import argparse as args
from Bio.PDB import *
from Bio.Seq import Seq
from Bio.PDB.PDBIO import PDBIO, Select
from Bio import pairwise2 as pw2
#import dash_test


####################
#   FUNCTIONS      #
####################

def check_type(input_chain):
    """
    Given a chain, it checks if it pertains to a DNA or protein structure.
    :param input_pdb_file:
    :return: structure object
    """

    atoms = set(Selection.unfold_entities(input_chain, "A"))
    for at in atoms:
        if at.id == "CA":
            return "Protein"
        elif at.id == "P":
            return "DNA"


class SelectNonHet(Select):
    """
    Select sub-class that saves a pdb file without heteroatoms.
    Not used for now.
    """
    def accepted_residues(self, residue):
        return 1 if residue.id[0] != " " else 0


def mutate_dna_chain(input_chain):
    """
    Given an input DNA chain, its C1' is transformed into a CA in order
    to make it a suitable input for the Superimposer and clashes functions.
    :param input_chain:
    :return:
    """

    if check_type(input_chain) == "Protein":
        pass
        return input_chain
    elif check_type(input_chain) == "DNA":
        print("Transforming DNA chains.....")
        new_chain = Chain.Chain(input_chain.id)
        for res in input_chain: # Should add copy of chain??? I don't really see why
            new_chain.add(res.copy()) # Shallow copy to actually copy the object attributes

        for res in new_chain:
            for atom in res:
                if atom.id == "C1'":
                    atom.id = "CA"

        return new_chain
    else:
        print("The chain introduce was not PROT nor ACNUC.")


def chain_align(chain1, chain2):
    """Pairwise alignment between two chains.
    NOT USED FOR NOW.
    """
    ppb = PPBuilder()
    for pp in ppb.build_peptides(chain1):
        seq1 = pp.get_sequence()
        print(f"From {chain1} we get {seq1}")
    for pp in ppb.build_peptides(chain2):
        seq2 = pp.get_sequence()
        print(f"From {chain2} we get {seq2}")
    alignments = pw2.align.globalxx(seq1, seq2, score_only=True) # Only scores to save time
    min_length = max(len(seq1), len(seq2)) # Why min and not max?
    identity_perc = round(alignments / min_length, 2)
    print(f"Between {chain1} and {chain2} there is a {identity_perc}% of identity")
    if identity_perc > 0.95: # Should be modificable
        return chain1, chain2

def acnucseq_from_pdb(input_chain):
    acnucseq = ""
    for res in input_chain:
        if len(res.get_resname()) > 1:
            acnucseq += str(res.get_resname().rstrip()[2:])
        else:
            acnucseq += str(res.get_resname()).rstrip()

    return acnucseq


def similar_chains(input_model_dict, type):
    """
    Given an dictionary of PDB files as Keys, and their Structure objects as values, this function
    makes pairwise alignments between the chains, keeping only those with a 95% or higher
    similarity.
    """
    print("Looking for similar chains...")
    similar_chains = {}
    ppb = PPBuilder()
    model_list = list(input_model_dict.values())
    seq1 = ""
    seq2 = ""
    for index, model in enumerate(model_list):
        for chain1 in model:
            if check_type(chain1) == type:
                #print(chain1, type)
                for model2 in model_list[index+1:]:
                    for chain2 in model2:
                        if check_type(chain2) == type:
                            #print(chain2)

                            if type == "DNA":
                                seq1 = acnucseq_from_pdb(chain1)
                                seq2 = acnucseq_from_pdb(chain2)
                                alignments = pw2.align.globalxx(seq1, seq2, score_only=True)
                                min_length = max(len(seq1), len(seq2))
                                identity_perc = round(alignments / min_length, 2)
                                #print(identity_perc)
                                if identity_perc > 0.95:
                                    similar_chains.setdefault(chain2, chain1)


                            elif type == "Protein":
                                for pp1 in ppb.build_peptides(chain1):
                                    seq1 = pp1.get_sequence()
                                for pp2 in ppb.build_peptides(chain2):
                                    seq2 = pp2.get_sequence()
                                alignments = pw2.align.globalxx(seq1, seq2, score_only=True)
                                min_length = max(len(seq1), len(seq2))
                                identity_perc = round(alignments / min_length, 2)
                                #print(identity_perc)
                                if identity_perc > 0.95:
                                    similar_chains.setdefault(chain2, chain1)
    return similar_chains


def common_chain_res(chain1, chain2):
    """
    Given a pair of chains, obtains the common residues and their respective CA atoms.
    Returns a tuple with two lists: The first one contains the list of atoms corresponding
    to the first chain. The second list contains the list of atoms of the second chain.
    :param chain1:
    :param chain2:
    :return:
    """
    print("Getting the common chains...")
    #res1 = [res for res in chain1 if res["CA"] or res["P"]
    #res2 = [res for res in chain2 if res["CA"] pr res["P"]]

    res1 = []
    res2 = []

    for res in chain1:
        for atom in res:
            if atom.id == "CA":
                res1.append(res)

    for res in chain2:
        for atom in res:
            if atom.id == "CA":
                res2.append(res)

    common_res1 = [res1 for res1, res2 in zip(res1, res2) if res1.get_resname() == res2.get_resname()]
    common_res2 = [res2 for res1, res2 in zip(res1, res2) if res1.get_resname() == res2.get_resname()]

    #chain1_atoms_list = [res["CA"] for res in common_res1]
    #chain2_atoms_list = [res["CA"] for res in common_res2]

    chain1_atoms_list = []
    chain2_atoms_list = []

    for res in common_res1:
        for atom in res:
            if atom.id == "CA":
                chain1_atoms_list.append(atom)

    for res in common_res2:
        for atom in res:
            if atom.id == "CA":
                chain2_atoms_list.append(atom)

    common_atoms = (chain1_atoms_list, chain2_atoms_list)

    return common_atoms


def get_atoms_list(chain):
    """Creates a list of the atoms only taking CA or P for protein and acid nucleics, respectively.
    This list of atoms will be lately used in the superimposition process."""
    type_chain = check_type(chain)
    if type_chain == "Protein":
        atom_id = "CA"
    elif type_chain == "DNA":
        atom_id = "P"
    atoms = chain.get_atoms()
    atoms_list = []
    for atom in atoms:
        if atom.id == atom_id:
            atoms_list.append(atom)
    return atoms_list


def Superimpose_modified(chain1, chain2):
    """
    Given a pair of chains, these are superimposed. The first input is the
    fixed chain, and the second one the moving chain. The rotran matrix is computed, applied to the
    second chain, and then the RMSD is computed.
    :param :
    :return: RMSD score between structures
    """

    # Set fixed, moving models
    fixed_chain = chain1
    moving_chain = chain2

    # Create the LIST OF ATOMS to be aligned
    fixed_atoms = get_atoms_list(fixed_chain)
    moving_atoms = get_atoms_list(moving_chain)

    # When superimposing chains are not equally sized
    if len(fixed_atoms) != len(moving_atoms):
        if len(fixed_atoms) > len(moving_atoms):
            fixed_atoms=fixed_atoms[:len(moving_atoms)]
        else:
            moving_atoms=moving_atoms[:len(fixed_atoms)]

    imposer = Superimposer()
    imposer.set_atoms(fixed_atoms, moving_atoms)
        #imposer.apply(moving_chain.get_atoms())
    return(imposer.rms)


# FOR THE SUPERIMPOSE: If we're mutating the DNA chains, we don't really need to
# separate them, right?

def Superimpose(chain1, chain2):
    """
    Given a pair of chains, these are superimposed. The first input is the
    fixed chain, and the second one the moving chain. The rotran matrix is computed, applied to the
    second chain, and then the RMSD is computed.
    :param :
    :return: RMSD score between structures
    """

    print("Superimposing the structures...")

    # Set fixed, moving models
    fixed_chain = chain1
    moving_chain = chain2

    # Create the LIST OF ATOMS to be aligned
    fixed_atoms = []
    moving_atoms = []
    if check_type(chain1) == "Protein" and check_type(chain2) == "Protein":
        for res in fixed_chain:
            for atom in res:
                if atom.id == "CA":
                    fixed_atoms.append(atom)

        for res in moving_chain:
            for atom in res:
                if atom.id == "CA":
                    moving_atoms.append(atom)

        #fixed_atoms = [res["CA"] or res["P"] for res in chain1 if res.id[0] == " "]
        #moving_atoms = [res["CA"] or res["P"] for res in chain2 if res.id[0] == " "]

    elif check_type(chain1) == "DNA" and check_type(chain2) == "DNA":
        for res in fixed_chain:
            for atom in res:
                if atom.id == "P":
                    fixed_atoms.append(atom)

        for res in moving_chain:
            for atom in res:
                if atom.id == "P":
                    moving_atoms.append(atom)

        # fixed_atoms = [res["CA"] or res["P"] for res in chain1 if res.id[0] == " "]
        # moving_atoms = [res["CA"] or res["P"] for res in chain2 if res.id[0] == " "]

    # When superimposing chains are not equally sized
    if len(fixed_atoms) != len(moving_atoms):
        common_atoms = common_chain_res(fixed_chain, moving_chain)

        imposer = Superimposer()
        imposer.set_atoms(common_atoms[0], common_atoms[1])
        #imposer.apply(moving_chain.get_atoms())
        return(imposer.rms)

    # If they are the same size
    else:
        imposer = Superimposer()
        imposer.set_atoms(fixed_atoms, moving_atoms)
        #imposer.apply(moving_chain.get_atoms())
        return(imposer.rms)



def clashes(chain1, chain2, dist=1.9):
    """
    Computes the steric clashes between two chains.
    Change the clash distance cutoff?
    :param chain1:
    :param chain2:
    :return:
    """

    print("Obtaining the clashes....")

    atoms1 = [atom for atom in chain1.get_atoms()]
    atoms2 = [atom for atom in chain2.get_atoms()]

    neighborsearch = NeighborSearch(atoms2)

    neighbors = set()

    for atom in atoms1:
        center = atom.get_coord()
        neighbs = neighborsearch.search(center, dist, level='C')
        for x in neighbs:
            if x != chain1:
                neighbors.add(x)

    if len(neighbors) != 0:
        return True
    else:
        return False

def change_chain_id(file, i, j):
    """
    Given a PDB MODEL object, the Chain IDs are converted from a letter format to a numeric format.
    It also checks for heteroatoms, and detaches them from the file. ¿Pero esto no funciona?
    :param:
    :return:
    """
    #struct_list = []
    #structure = PDBParser(PERMISSIVE=True, QUIET=True).get_structure(f"{j}", file)
    #chains = structure[0].get_chains()

    current_model = file
    current_model.id = j
    j += 1
    for chain in current_model:
        heteroatoms = list(filter(lambda x: x.id[0] != " ", chain.get_residues()))
        for heteroatom in heteroatoms:
            chain.detach_child(heteroatom.id)
        chain.id = i # replace chain ID with number
        i += 1 # i will keep increasing depending on the number of models
    return current_model