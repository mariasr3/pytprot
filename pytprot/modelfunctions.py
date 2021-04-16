import time, random
from Bio.PDB import *
from Bio.PDB.PDBIO import PDBIO, Select
import inputfunctions, chainfunctions

################
# Functions employed for building the macrocomplex through interacting chains
################


######### Functions employed for the superimposition

def get_atoms_list(chain):
    """ Given a Bio.PDB.Chain object, returns a list of its CA or P atoms for protein and nucleic acids, respectively. """

    type_chain = inputfunctions.check_type(chain)
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


def common_chain_res(model1, model2):
    """
    Given a pair of Bio.PDB.Chain objects, returns the common residues and their respective CA atoms.
    Returns a tuple with two lists: The first one contains the list of atoms corresponding
    to the first chain. The second list contains the list of atoms of the second chain, both with the
    same length.
    """


    res1 = []
    res2 = []
    chain1_atoms_list = []
    chain2_atoms_list = []

    for res in model1:
        res1.append(res)

    for res in model2:
        res2.append(res)

    # Looking for the common atoms
    common_res1 = [res1 for res1, res2 in zip(res1, res2) if res1.get_resname() == res2.get_resname()]
    common_res2 = [res2 for res1, res2 in zip(res1, res2) if res1.get_resname() == res2.get_resname()]


    for res in common_res1:
        for atom in res:
            if atom.id == "CA" or atom.id == "C1'":
                chain1_atoms_list.append(atom)

    for res in common_res2:
        for atom in res:
            if atom.id == "CA" or atom.id == "C1'":
                chain2_atoms_list.append(atom)

    common_atoms = (chain1_atoms_list, chain2_atoms_list)  # Tuple with same-length lists

    return common_atoms

def extract_key(chain, unicommon):
    """
    Given a specific chain, it finds the "unique" chain to which it is similar to, as it looks for it
    within the "unicommon" type dictionary.

    If the provided chain is already a key, just returns the chain,
    If not, returns the actual "unique" chain (its corresponding key within the unicommon dictionary).
    """

    list_keys = [ch.id for ch in unicommon.keys()]

    if chain.id in list_keys:
        pair_key = chain
        return pair_key
    else:
        for key, value in unicommon.items():
            list_values = [ch.id for ch in value]
            if chain.id in list_values:
                pair_key = key
                return pair_key

def superimpose_pair(pair1, pair2, i, j):
    """
    Given a pair of chains, these are superimposed. The first input is the
    fixed chain, and the second one the moving chain. The rotran matrix is computed, applied to the
    second chain, and then the RMSD is computed.
    :param :
    :return: RMSD score between structures
    """

    #print("Superimposing the structures...")

    equal_chains_moved=[]
    # Set fixed, moving models
    fixed_chain = pair1[i]
    moving_chain = pair2[j]
    pair_fixed_chain = pair1[abs(i-1)]
    pair_moving_chain = pair2[abs(j-1)]

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
        imposer.apply(pair_moving_chain.get_atoms())
        return (imposer.rms, pair_moving_chain)

    # If they are the same size
    else:
        imposer = Superimposer()
        imposer.set_atoms(fixed_atoms, moving_atoms)
        imposer.apply(pair_moving_chain.get_atoms())
        return (imposer.rms, pair_moving_chain)

def superimpose(chain1, chain2):
    """
    Given a pair of Bio.PDB.Chain objects, it superimposes their CA backbone common atoms.
    The first input is the fixed chain, and the second one the moving chain.
    Returns a Bio.PDB,Superimposer object that encodes the rotran matrix.
    If the input chains are empty, it returns None.
    """

    #print("Superimposing the structures...")


    # Set fixed, moving models
    fixed_model = chain1
    moving_model = chain2

    # Create the LIST OF ATOMS to be aligned
    fixed_atoms = []
    moving_atoms = []

    for res in fixed_model:
        for atom in res:
            if atom.id == "CA" or atom.id == "C1'":
                fixed_atoms.append(atom)

    for res in moving_model:
        for atom in res:
            if atom.id == "CA" or atom.id == "C1'":
                moving_atoms.append(atom)

    # When superimposing chains are not equally sized
    if len(fixed_atoms) != len(moving_atoms):
        common_atoms = common_chain_res(fixed_model, moving_model)
        imposer = Superimposer()
        if len(common_atoms[0]) != 0 or len(common_atoms[1]) != 0:
            imposer.set_atoms(common_atoms[0], common_atoms[1])
            return imposer
        else:
            return False

    # If they are the same size
    else:
        imposer = Superimposer()
        if len(fixed_atoms) != 0 or len(moving_atoms) != 0:
            imposer.set_atoms(fixed_atoms, moving_atoms)
            return imposer
        else:
            return False

def clashes(chain1, chain2, dist=2.5):
    """
    Given two Bio.PDB.Chain objects, this function will look for any possible steric clashes between them.
    - dist: Refers to the distance in Å at which we consider the two chains to be clashing.
    Returns a Boolean that indicates if there are clashes (True) or not (False).
    """

    atoms1 = [atom for atom in chain1.get_atoms()]
    atoms2 = [atom for atom in chain2.get_atoms()]

    neighborsearch = NeighborSearch(atoms2)

    neighbors = set()
    for atom in atoms1:
        center = atom.get_coord()
        # Look for Chain objects at the specified distance from each atom of the chain2
        neighbs = neighborsearch.search(center, dist, level='C')
        for x in neighbs:
            if x != atoms1:
                neighbors.add(x)

    if len(neighbors) != 0:  # If there is just ONE contact between both chains, we consider them to be clashing
        return True
    else:
        return False


def clash_list (chain1, list_chains, dist=2.5):
    """
    Given a Bio.PDB.Chain object, and a list of Bio.PDB.Chains, computes the clashes between said
    chain and all the chains located in the model, at the distance indicated.

    Returns an integer representing the number of structures with which this input chain finds clashes.
    """

    clash_count = 0

    atoms1 = [atom for atom in chain1.get_atoms()]

    for chain2 in list_chains:
        atoms2 = [atom for atom in chain2.get_atoms()]
        neighborsearch = NeighborSearch(atoms2)

        neighbors = set()

        for atom in atoms1:
            center = atom.get_coord()
            neighbs = neighborsearch.search(center, dist, level='A')
            for x in neighbs:
                if x != atoms1:
                    neighbors.add(x)

        if len(neighbors) != 0:
            #print(f"The chain with which there are clashes is {chain2}")
            clash_count += 1

    return clash_count




def convert_key (unicommon, current_key):
    """
    Given the "unicommon" dictionary and a given Bio.PDB.Chain object, this function looks for the
    corresponding key within this dictionary of the provided chain.

    Returns the "unique" key.
    """

    i = 0
    for key in unicommon.keys():
        if inputfunctions.check_type(key) == "Protein":
            if key.id == current_key.id:
                new_key = i
                return new_key
            i += 1

def stoich_count(current_stoich, stoich_dict, unicommon, key_to_add):
    """
    Given the provided stoichiometry, the current model stoichiometry while being built,
    and the "unicommon" dictionary, this function checks, before adding a new chain (key_to_add),
    if the stoichiometry provided has been fulfilled or not.
    """

    current_stoich_id = {}
    i = 0
    for key, values in current_stoich.items():
        if key_to_add is not None and key == key_to_add.id:  # Quick fix: key_to_add is None sometimes
            new_key = convert_key(unicommon, key_to_add)
            current_stoich_id[new_key] = values
            try:
                if current_stoich_id[new_key] >= stoich_dict[new_key]:
                    return True
                else:
                    return False
            except KeyError: # If the key is not found in the stoich_dict
                continue


def model_construction(unicommon, pdb_dict, equivalent_chains, stoichiometry_input=None, verbose=False, forcing=False):
    """
    Central pytprot multicomplex protein assembly function. With the result of all the chain pre-processing, it will
    build the final model. The input files are the following:

        unicommon: Dictionary with unique chain keys and similar chains to said "key-chains" as values.

        pdb_dict: Dictionary with Path to files as keys (or placeholder names for when we have a multicomplex)
        and Bio.PDB.Structure objects as values that contain pairs of interacting chains.

        equivalent_chains: Dictionary with the chains with letters as IDs and their equivalent numerical Chain ID.

        stoichiometry_input: By default, a stoichiometry file is not considered. On the contrary, if added
        when running the script, it takes the dictionary with the letters-transformed-to-numbers of the
        provided stoichiometry as keys, and the number of subunits of said number the model needs to have.

        verbose: If true, prints to log information regarding the building process.

        forcing: If true, the model will try to force any Chain that has been indicated on the stoichiometry file,
        but, due to clashes, could not be introduced in the model while building it.


    Returns a Bio.Pdb.Model.Model object with the built model.

    """

    new_model = Model.Model('X')

    stoichiometry = False
    if stoichiometry_input is not None:
        stoichiometry = True

    if stoichiometry:

        # Check, based on the Unicommons dictionary, if the provided stoichiometry is possible
        stch_check = chainfunctions.stoichiometry_check_simple(unicommon, stoichiometry_input)

        if stch_check == True:
            if verbose:
                print("Correct stoichiometry.\n")
            pass
        #If the stoichiometry is incorrect the model will be build as if no stoichiometry has been provided
        else:
            if verbose:
                print("Stoichiometry is incorrect, the model will add as many chains as possible\n")

    if verbose:
        print(f"\n############### MODEL CONSTRUCTION ############# \n")

    # Iterating over the models
    model_pairs = []
    listmodels = [[model_pairs.append(model) for model in str] for str in pdb_dict.values()]
    del listmodels

    #List of list of pairs of interacting chains that will be superimpose to construct the model
    listpairs = []
    for model in model_pairs:
        pairs = []
        for ch in model:
            pairs.append(ch)
        listpairs.append(pairs)


    if verbose:
        print(f"We have {len(listpairs)} total pairs")

    # Equivalent chains: dctionary where the keys are the chain id in numbers and the values the chain with id in letters
    equivalent_chains_id = {}
    added_chains = []
    for x, y in equivalent_chains.items():
        equivalent_chains_id[x.id] = y

    ## SEED FOR MODEL: Random.

    current_stoich = {}

    for chain in listpairs[int(random.randrange(0, len(listpairs)))]:
        added_chains.append(chain) # A list of chains is created to look for clashes when a new chain is added
        chain_to_add = equivalent_chains_id[chain.id] # Change the chain for its equivalent with a letter ID.
        new_model.add(chain_to_add) # The first pair will be added to the model and used as a starting pair
        key = extract_key(chain, unicommon)
        if inputfunctions.check_type(chain) == "Protein" and key != None:
            if key.id not in current_stoich:
                current_stoich[key.id] = 1
            else:
                current_stoich[key.id] += 1

    for pair2 in listpairs:
        model_list = list(ch for ch in added_chains)

        for idx, chain1 in enumerate(pair2):
            key1 = extract_key(chain1, unicommon)  # To know if two chains are equal and so they can be superimposed
                                                   # their keys of unicommon dictionray will be compared
            for chain2 in model_list:
                key2 = extract_key(chain2, unicommon)

                # Only if the keys are equal we can superimpose two chains
                if key1 == key2:
                    imposer = superimpose(chain1, chain2) # Superimpose a chain of each pair that share a key
                    moved_chain = pair2[abs(idx - 1)].copy()
                    if imposer is not False:
                        imposer.apply(moved_chain.get_atoms()) # Apply the rotran matrix to the other chain of pair2

                        # The chain will just be added to the model if there are no clashes with other chains of the model
                        clash_count = clash_list(moved_chain, added_chains,
                                                 dist=0.5)  # As seen on the Chimera documentation
                        if clash_count == 0:
                            chain_to_add = equivalent_chains_id[moved_chain.id]

                            #If stoichiometry is provided and correct it will check it before adding the chain
                            if stoichiometry:
                                if stch_check == True:
                                    # Only checks the stoichiometry if it is a protein chain
                                    if inputfunctions.check_type(moved_chain) == "Protein":
                                        key_to_add = extract_key(moved_chain, unicommon)
                                        try:
                                            # If the stoichiometry for this chain is already fullfilled the program will stop and continue with the next chain
                                            if stoich_count(current_stoich, stoichiometry_input, unicommon, key_to_add) == True:
                                                break
                                        except:
                                            pass
                            if chain_to_add not in list([ch for ch in new_model.get_chains()]): # If the chain is not already in the model
                                added_chains.append(moved_chain) # Add it to the list of chains
                                new_model.add(chain_to_add) # Add the chain with letter ID to the model
                                key_to_add = extract_key(moved_chain, unicommon)

                                # Update the current stoichiometry dictionary with the new chain added
                                if inputfunctions.check_type(moved_chain) == "Protein":
                                    if key_to_add is not None and key_to_add.id not in current_stoich:
                                        current_stoich[key_to_add.id] = 1
                                        if verbose:
                                            print(f"chain {chain_to_add} NEW KEY\n")
                                    else:
                                        if key_to_add is not None:
                                            current_stoich[key_to_add.id] += 1



    new_model_chains = [ch for ch in new_model.get_chains()]

    # If stoichiometry is provided and correct it will "force" chains of unicommons to be added to the model until the stoichiometry is completed
    if forcing and stoichiometry:
        if stch_check == True and verbose:
            print("Force the model")
            chains_to_model = []
            for key, values in unicommon.items(): #iterate through keys
                if key != None:
                    chain_to_add = equivalent_chains_id[key.id]
                    if chain_to_add not in new_model_chains: # check if the chain is not in the model

                        # If the stoichiometry for this chain is already fullfilled the program will stop and continue with the next chain
                        if inputfunctions.check_type(chain_to_add) == "Protein":
                            try:
                                if stoich_count(current_stoich, stoichiometry_input, unicommon, key) == True:
                                    break
                            except:
                                pass

                        # If not the chain will be added to the model
                        new_model_chains.append(chain_to_add)
                        new_model.add(chain_to_add)
                        if inputfunctions.check_type(chain_to_add) == "Protein": # Update the stoichiometry
                            current_stoich[key.id] = 1

                for val in values: # Iterate through values
                    chain_to_add = equivalent_chains_id[val.id]
                    if chain_to_add not in new_model_chains:
                        key_val = extract_key(val, unicommon)
                        if key_val is not None:

                            # If the stoichiometry for this chain is already fullfilled the program will stop and continue with the next chain
                            if inputfunctions.check_type(chain_to_add) == "Protein":
                                try:
                                    if stoich_count(current_stoich, stoichiometry_input, unicommon, key_val) == True:
                                        break
                                except:
                                    pass

                        # If not the chain will be added to the model
                        new_model_chains.append(chain_to_add)
                        new_model.add(chain_to_add)
                        #print(f"{chain_to_add} has been added to the model")
                        if inputfunctions.check_type(chain_to_add) == "Protein": # Update the stoichiometry
                            if key_val.id not in current_stoich:
                                current_stoich[key_val.id] = 1
                            else:
                                if key_val is not None:
                                    current_stoich[key_val.id] += 1

    if verbose:
        print("\n")
        print(f"NEW MODEL WITH:\n{len(new_model_chains)} CHAINS\nWHICH ARE: {new_model_chains}")


    print("The model has been correctly built.")

    return new_model


def save_model(final_model, pdblist, outdir, verbose=False, macrocomplex=False):
    """
    Given the built, final model, this function saves it in its correct format in the provided output directory.
        final_model: Refers to the returned model of the model_construction() function
        pdblist: List of path input pdb files. Employed to fetch the model name.
        outdir: To obtain the relevant output directory in which the model will be saved.
        verbose: If True, prints on the terminal the progress of the function.
    """

    new_model_chains = [ch for ch in final_model.get_chains()]

    #### SAVING THE STRUCTURE
    parser=PDBParser()
    #model_structure=parser.get_structure(final_chains)
    io = PDBIO()
    io.set_structure(final_model)
    if macrocomplex:
        model_name = (str(pdblist[0])).split("/")[-2]
    else:
        model_name = (str(pdblist[0]).split("/")[-1]).split("_")[0]
    # timestamp
    model_out_path = outdir
    t = time.localtime()
    timestamp = time.strftime('%d%m%Y_%H:%M', t)

    # Saving the structure
    if len(new_model_chains) != 0:
        io.save(f"{model_out_path}/{model_name}_complex_{len(new_model_chains)}_chains_{timestamp}.pdb")

        if verbose:
            print(f"\nSaving the structure...\n")

    print("The model has been correctly saved.")
