import os, copy
import sys
from Bio.PDB import *
import chainfunctions, modelfunctions


################
# Functions employed for processing the different input files
################


##### Input file parsing functions
## These are different functions that retrieve the uploaded files, be it multiple interacting chain-pair PDBs,
## a multicomplex input or the stoichiometry.txt file

def macrocomplex_parser(pdbfile, interact_dist=8, interact_contacts=8, redundant_dist=2.5, verbose = False):
    """
    Given a macrocomplex PDB file, it computes the interacting pairs, deletes the redundant ones
    and converts each pair of interacting chains into a specific Bio.Structure.

    Input parameters:
    - pdbfile: The macrocomplex PDB file
    - interact_dist: indicates the distance in Angstroms at which two atoms of two chains are considered to be
    interacting.
    - interact_contacts: Number of contacts ad the distance indicated by interact_dist above which two chains
    are considered to be interacting. In conjunction with interact_dist, these parameters dictate the interacting
    chains found in the macrocomplex.
    - redundant_dist: Distance at which the pair-chain of two interacting-chain-pairs that share a similar chains
    and have been superimposed, in order to consider it to be a redundant interaction.
    - verbose: Prints information to terminal.

    Returns a dictionary with a placeholder name (instead of the path, as obtained in pdb_parser) as
    key and a Bio.PDB.Structure object with pairs of interacting chains as a value.

    """

    # PDB parsing

    pdbname = [file for file in os.listdir(pdbfile) if file.endswith(".pdb")][0] # turn to string
    structure_name = (pdbfile.split("/"))[-1]  # Get name from file
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure(structure_name, pdbfile+"/"+pdbname)
    structure.id = 0  # placeholder ID



    # Obtaining the interacting pairs

    if verbose == True:
        print("Looking for interacting pairs...")

    chains_list = [ch for ch in structure.get_chains()]
    pairs = []
    pairs2 = set()  # Same as pairs but a plain list. Useful for debugging purposes.

    for idx, ch in enumerate(chains_list):
        for ch2 in chains_list[idx+1:]:
            if ch.id != ch2.id: # Avoid comparison between the same IDs (same chains)

                interaction = chainfunctions.interacting_pairs(ch, ch2, dist=interact_dist, num_contacts=interact_contacts)

                if interaction == True:
                    pairs.append([ch, ch2]) # Pairs that do interact

    if verbose == True:
        print(f"We have {len(pairs)} interacting pairs\n")

    # Eliminating redundant pairs

    pairs_without_redundant = chainfunctions.redundant_pairs(pairs, dist=redundant_dist, verbose=verbose)

    # Building the "fake" pdb_dict

    if verbose == True:
        print("Creating the pdb_dict file...")

    pairs_dict = {}  # Will have the same format as the str_dict from the "pdb_parser" function
    pdbname = pdbname.split(".")[0]

    i = 0
    for pair in pairs_without_redundant: # Loop to "fill in" the pairs_dict dictionary
        str_pair = Structure.Structure(f"{i}") # Need to create a placeholder Bio.PDB Structure and Model
        model_pair = Model.Model(f"{i}")
        model_pair.add(pair[0])
        model_pair.add(pair[1])
        str_pair.add(model_pair)
        pairs_dict[f"Structure_{pdbname}_{i}"] = str_pair # Placeholder name as key
        i += 1

    if verbose == True:
        print("FILE-MODEL dictionary successfully created.")

    return pairs_dict


def pdb_parser(list_pathpdbs, verbose=False):
    """
    Given a list of paths pointing to a specific PDB file, it is read, converted to a specific Bio.Structure.
    "verbose" is a parameter that, when set to True, prints the information to the Terminal.
    Returns a dictionary with the path as key and a Bio.PDB.Structure that contains a pair of interacting chains as a value.
    """

    # Get dictionary of structures
    final_name = (list_pathpdbs[0].split("/")[-1]).split("_")[0]
    str_dict = {}

    i = 0
    for str in list_pathpdbs:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure(f"Structure_{final_name}_{i}", str)
        str_dict[str] = structure  # Path+file as KEY, structure object as VALUE
        i += 1

    if verbose == True:
        print("FILE-MODEL dictionary successfully created.")

    return str_dict


def stoichiometry_parser(stoichiometry_file):
    """
    Given a .txt file parsed from the argparse with a specific , it returns a dictionary with the chain ID as a key and the number of
    subunits of said chains as values.
    """

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz" # needed to change chain IDs

    i = 0
    stechio_dict = {}
    for line in open(stoichiometry_file):
        letter = str(line.split(":")[0].rstrip())
        num = int(line.split(":")[1].rstrip())
        new_id = alphabet.find(letter)  # Use the index in the alphabet string of the chain ID as the new chain ID
        stechio_dict[new_id] = num
        i += 1

    return stechio_dict

################
# Chain processing
# These functions produce the necessary modifications to the str_dict returned by
# the macrocomplex_parser and pdb_parser functions
################

def chain_processing(str_dict, verbose=True):
    """
    Given a dictionary with pdb structures of interacting pairs (a "str_dict" type of dictionary), this function:
    1. Creates the equivalence dictionary between alphabetical and numerical chain IDs.
    2. Prints, if indicated through the "verbose" parameter, the number of DNA and Protein chains.
    3. Changes the chain IDs to numbers according to an alphabet.

    :param pdbdict:
    :return: quivalence_chains, str_dict_final
    """

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    # alphabet includes both upper and lowercase letters in order to be able to include more chains
    # in the model


    # Creating the dictionary for equivalence between chain IDs

    equivalence_chains = {}

    eq_str_dict = copy.deepcopy(str_dict)  # deepcopy in order to have different id() for the object and all references

    #for struct in eq_str_dict.values():
    #    for chains in struct.get_chains():
    #        chains2 = chains.copy()  # Work with the copies to not modify the original chain
    #        new_chain = chains.copy()
    #        if not isinstance(chains.id, int):  # Filtering any possible numerical ID
    #            chains.id = alphabet.find(str(new_chain.id))
    #            str_dict_final[name] = struct_copy

    i = 0
    for struct in eq_str_dict.values():
        for chains in struct.get_chains():
            i += 1
            chains2 = chains.copy()  # Work with the copies to not modify the original chain
            new_chain = chains.copy()
            new_chain.id = alphabet.find(str(chains.id))
            if new_chain.id == -1:
                new_chain.id = 100 + i

            equivalence_chains.setdefault(new_chain, chains2)


    # Counting DNA and nucleotide chains
    # This is only useful if the verbose flag (-v) is active, as these objects are not used later-on in the script
    dnachainlist = set()
    protchainlist = set()

    for structure in str_dict.values():
        for model in structure:
            for chain in model:
                if check_type(chain) == "Protein":
                    if type(chain.id) == str:
                        protchainlist.add(chain)
                elif check_type(chain) == "DNA":
                    if type(chain.id) == str:
                        dnachainlist.add(chain)

    if verbose == True:
        print(f"There are {len(dnachainlist)} DNA chains and {len(protchainlist)} protein chains")

    del dnachainlist, protchainlist  # To free memory space

    # Changing the chain IDs

    str_dict_final = {}
    new_str_dict = copy.copy(str_dict)  # Avoid modifying the original str_dict

    i = 0
    for name, structure in new_str_dict.items():
        struct_copy = copy.deepcopy(structure)
        for chain in struct_copy.get_chains():
            i += 1
            new_chain = chain.copy()
            chain.id = alphabet.find(str(new_chain.id))
            if chain.id == -1:
                chain.id = 100 + i

            str_dict_final[name] = struct_copy


    if verbose == True:
        print("Chain IDs have been converted from letters to numbers.\n")

    return [equivalence_chains, str_dict_final]



################
# Nucleotide chain processing
# These functions are specifically designed to treat the input nucleotide Chains
################

def check_type(input_chain):
    """
    Given a Bio.PDB.Chain object, it checks if it pertains to a DNA or protein structure.
    Returns a string that indicates the structure type.
    """

    atoms = set(Selection.unfold_entities(input_chain, "A"))  # We do not need a list with all the atoms, a set is memory savvy
    atoms_id = [at.id for at in atoms]

    if "CA" in atoms_id:
        return "Protein"
    if "P" in atoms_id:
        return "DNA"


def mutate_dna_chain(input_chain):
    """
    Given an input DNA chain, its C1' is transformed into a CA in order
    to make it a suitable input for the Superimposer and clashes functions (see modelfunctions.py).
    Returns a copy of the chain with the mutated C1'.
    """

    if check_type(input_chain) == "Protein":
        return input_chain

    elif check_type(input_chain) == "DNA":
        new_chain = Chain.Chain(input_chain.id)
        for res in input_chain:
            new_chain.add(res.copy())  # shallow copy is enough to copy the residue information

        for res in new_chain:
            for atom in res:
                if atom.id == "C1'":
                    atom.id = "CA"

        return new_chain
    else:
        print("The chain introduced was not a Protein chain nor an Nucleotide chain.")


def acnucseq_from_pdb(input_chain):
    """ Given an input DNA chain, returns a string with its nucleotide sequence. """

    acnucseq = ""
    for res in input_chain:
        if len(res.get_resname()) > 1:
            acnucseq += str(res.get_resname().rstrip()[2:])
        else:
            acnucseq += str(res.get_resname()).rstrip()

    return acnucseq
