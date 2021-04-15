# Adding module location to sys.path
from Bio.PDB import *
from Bio import pairwise2 as pw2
import inputfunctions, modelfunctions


################
# Functions employed for processing and re-arranging the different chains in
# order to obtain a suitable input for the model building
################


# Chain-processing functions for the case of a macrocomplex input


def interacting_pairs(chain1, chain2, dist, num_contacts, verbose=False):
    """
    Given a pair of PDB Chain objects as an input, it indicates if they are interacting, based on two thresholds:
    num_contacts and dist, which represent the number of "contacts" between these two chains and the
    distance at which these contacts are being analysed.

    Returns a Boolean object indicating if they are interacting or not.

    """

    # Check type and mutate DNA chain whenever necessary
    if inputfunctions.check_type(chain1) == "DNA":
        chain1 = inputfunctions.mutate_dna_chain(chain1)
    elif inputfunctions.check_type(chain2) == "DNA":
        chain2 = inputfunctions.mutate_dna_chain(chain2)

    # Set atom lists
    chain1_alpha = []
    chain11 = [[chain1_alpha.append(at) for at in res if at.id == "CA"] for res in chain1.get_residues()]
    chain2_alpha = []
    chain22 = [[chain2_alpha.append(at) for at in res if at.id == "CA" or at.id == "C1'"] for res in chain2.get_residues()]

    # Small bug, sometimes the chain2 is not properly mutated (C1' to CA), so this works for now.

    del chain11, chain22  # To free memory space, as these will not be used afterwards


    total_contacts = 0  # Contact counter

    if len(chain1_alpha) != 0 and len(chain2_alpha) != 0: # Filter for the cases in which the input chains are empty
        neighborsearch = NeighborSearch(chain2_alpha)

        for atom in chain1_alpha:
            center = atom.get_coord()
            neighbs = neighborsearch.search(center, dist, level='R')  # Look at every CA of the residues at the provided distance

            if len(neighbs) != 0: # If there is any residue at that distance, we can consider it a contact
                total_contacts += 1

    else:
        if verbose == True:
            print("Input chains are empty")

    if total_contacts >= num_contacts:
        return True
    else:
        return False


def redundant_pairs(interacting_pairs, dist=2.5, verbose=False):
    """
    Given a list of lists of lists that contains interacting chain-pairs, it deletes the redundant ones.
    Returns another list, similar to the input one, but with the redundant pairs eliminated.

    In order to classify two interacting chain-pairs as redundant, we will need to superimpose those that
    share a common chain, which we describe as fixed and moving pairs. Then, after superimposition, we
    compute the clashes between the chain of the fixed pair which is not the one with the common
    similar chain in the moving, and the fixed pair at the distance indicated. If there is a clash, we
    consider both of these interacting pairs to be redundant, and thus one of them is eliminated.
    """

    distance = dist
    chains_to_delete = []  # List of redundant interacting pairs that can be left out


    if verbose == True:
        print("Looking for redundant pairs...")


    for index, pair1 in enumerate(interacting_pairs):
        for pair2 in interacting_pairs[index+1:]:
            ident1 = chain_align(pair1[0], pair2[0])
            if ident1 is not None and ident1 > 0.95:

                moving_chain = modelfunctions.superimpose_pair(pair1, pair2, 0, 0)[1]  # Superimpose both pairs

                clashes = modelfunctions.clashes(pair1[1], moving_chain, dist=distance)
                # Check clashes of other two chains

                if clashes == True:                     # A clash indicates a pair of chains with similar sequence
                                                        # and a pair of chains with the same 3D disposition
                    chains_to_delete.append(pair2)
                break  # If the redundant pair was found for that first pair, exit the loop

            # Comparing first chain of first pairs, second chain of second pair
            # The pipeline is the same as for the first case
            ident2 = chain_align(pair1[0], pair2[1])

            if ident2 is not None and ident2 > 0.95:
                moving_chain = modelfunctions.superimpose_pair(pair1, pair2, 0, 1)[1]
                clashes = modelfunctions.clashes(pair1[1], moving_chain, dist=distance)
                if clashes == True:
                    chains_to_delete.append(pair2)
                break

            # Comparing second chain of first pair, first chain of second pair
            ident3 = chain_align(pair1[1], pair2[0])

            if ident3 is not None and ident3 > 0.95:
                moving_chain = modelfunctions.superimpose_pair(pair1, pair2, 1, 0)[1]
                clashes = modelfunctions.clashes(pair1[0], moving_chain, dist=distance)
                if clashes == True:
                    chains_to_delete.append(pair2)

                break

            # Comparing second chains of both pairs
            #
            ident4 = chain_align(pair1[1], pair2[1])

            if ident4 is not None and ident4 > 0.95:
                #print(f"{pair1[1]} and {pair2[1]} are similar, these will be superimposed.\n")
                moving_chain = modelfunctions.superimpose_pair(pair1, pair2, 1, 1)[1]
                clashes = modelfunctions.clashes(pair1[0], moving_chain, dist=distance)
                if clashes == True:
                    chains_to_delete.append(pair2)

                break

    if verbose == True:
        print(f"{len(chains_to_delete)} redundant pairs deleted.")

    return list(pair for pair in interacting_pairs if pair not in chains_to_delete)




# Chain-processing functions common to both input types


def chain_align(chain1, chain2):
    """
    Given a pair of Bio.PDB.Chain objects as an input, it makes a global pairwise alignment between two chains.
    Returns the identity percentage between both chains, normalizing by the length of the longest chain.
    If one of the input chains is empty, returns 0.
    """

    seq1 = ''
    seq2 = ''
    ppb = PPBuilder()

    if inputfunctions.check_type(chain1) == "Protein" and inputfunctions.check_type(chain2) == "Protein":
        for pp in ppb.build_peptides(chain1):
            seq1 = pp.get_sequence()
        for pp in ppb.build_peptides(chain2):
            seq2 = pp.get_sequence()
        alignments = pw2.align.globalxx(seq1, seq2, score_only=True)  # Only scores to save time
        max_length = max(len(seq1), len(seq2))
        if type(alignments) != list:  # If the input chains are empty, the alignment returns an empty list. This is a workaround.
            identity_perc = round(alignments / max_length, 2)
            return identity_perc
        else:
            return 0


def similar_chains(input_model_dict, type, verbose=False):
    """
    Given an dictionary of PDB files as Keys, and their Structure objects as values (str_dict-like object), this function
    makes pairwise alignments between the chains, keeping only those with a 95% or higher sequence
    similarity.

    Returns a dictionary that contains, as key-pairs, highly sequence-similar chains.

    """

    if verbose:
        print("Looking for similar chains...")

    similar_chains = {}

    ppb = PPBuilder()  # set ppb builder variables
    seq1 = ""
    seq2 = ""

    structure_list = list(input_model_dict.values())  # Fetch pdb_dict Structure objects

    model_list = []
    if structure_list[0].level == "S":  # Reassure that it is in fact a Structure object
        try:
            model_list = list(str[0] for str in structure_list)

        except KeyError:  # Managing the KeyError sometimes raised due to the macrocomplex input
            model_list = []
            model_list_tmp = [[model_list.append(model) for model in str] for str in structure_list]
            del model_list_tmp

    elif structure_list[0].level == "M": # If a Model is provided, instead
        model_list = list(str for str in structure_list)



    # Loop to look for similar chains

    for index, model in enumerate(model_list):
        for chain1 in model:
            if inputfunctions.check_type(chain1) == type: # Filter only protein or DNA chains
                for model2 in model_list[index+1:]: # Avoid repeated comparisons
                    for chain2 in model2:
                        if inputfunctions.check_type(chain2) == type:

                            # For nucleotide chains, specific function provides sequence

                            if type == "DNA":
                                seq1 = inputfunctions.acnucseq_from_pdb(chain1)
                                seq2 = inputfunctions.acnucseq_from_pdb(chain2)
                                alignments = pw2.align.globalxx(seq1, seq2, score_only=True) # Global alignment
                                max_length = min(len(seq1), len(seq2)) # Normalize on longest sequence
                                identity_perc = round(alignments / max_length, 2)
                                if identity_perc > 0.90:
                                    similar_chains.setdefault(chain2, chain1)


                            # For protein chains, PPBuilder provides sequences

                            else:
                                for pp1 in ppb.build_peptides(chain1):
                                    seq1 = pp1.get_sequence()
                                for pp2 in ppb.build_peptides(chain2):
                                    seq2 = pp2.get_sequence()
                                alignments = pw2.align.globalxx(seq1, seq2, score_only=True) # Global alignment
                                max_length = max(len(seq1), len(seq2)) # Normalize on longest sequence
                                if not isinstance(alignments, list):
                                    identity_perc = round(alignments / max_length, 2)
                                    if identity_perc > 0.95:
                                        similar_chains.setdefault(chain2, chain1)

    return similar_chains


def unique_common_chains(similar_chains, verbose=False):
    """
    Given a dictionary of pairwise sequence-similar chains (output of the similar_chains function),
    returns a "unique-common" chains-type dictionary, in which the keys are those chains to which most of the
    chains in the model are common-to, and as values a list of all the similar chains to said key-chain.
    """

    unicommon = {}

    if verbose == True:
        print("Looking for unique-common chains...")

    unique_chains = set(list(similar_chains.values())) # Use values as keys

    for chain in unique_chains:
        unicommon[chain] = []
        for keys, val in similar_chains.items():
            if val == chain and (keys != val):  # Fetch keys that have that specific value
                unicommon[chain].append(keys)

    return unicommon

def unicommon_completer(unicommon, stoichiometry_input, diff_chains):
    """ Small function that 'completes' the eunicommons dictionary with the stoichiometry provided """


    unicommon_keys = set([ch.id for ch in unicommon.keys()])
    print(f"These are the different chain IDs in unicommons: {unicommon_keys}")

    print(f"These are the different chain IDs: {diff_chains}")



    for stoich_input_key in list(stoichiometry_input.keys()):
        print(f"This is the input key:{stoich_input_key}")
        if stoich_input_key not in unicommon_keys:
            print(f"{stoich_input_key} is missing from unicommons")
            new_chain = [ch for ch in diff_chains if ch.id == stoich_input_key]
            print(new_chain)
            if len(new_chain) != 0:
                unicommon[new_chain[0]] = []

    return unicommon


##### STECHIOMETRY CHECK

def stoichiometry_check_simple(unicommon, stoichiometry_dict):
    """
    Given the parsed stoichiometry .txt file, this function checks if it is actually possible or not,
    based on the "unicommon" dictionary type. It checks if the sum of the required number of subunits
    can be achieved with the amount of unique and similar chains the program has found.

    Returns a Boolean that indicates if the stoichiometry is possible or not.
    """

    stoich_subunits = sum(int(x) for x in stoichiometry_dict.values())
    present_subunits = len(unicommon) + (sum(len(x) for y, x in unicommon.items()))

    if present_subunits >= stoich_subunits:
        return True
    else:
        return False
