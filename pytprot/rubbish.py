
def clashes(model1, model2, type, dist=1.9):
    """
    Computes the steric clashes between two chains.
    Change the clash distance cutoff?
    :param chain1:
    :param chain2:
    :return:
    """

    #print("Obtaining the clashes....")

    if type == "atom":
        atoms1 = [atom for atom in model1.get_atoms()]
        atoms2 = [atom for atom in model2.get_atoms()]

        neighborsearch = NeighborSearch(atoms2)

        neighbors = set()

        for atom in atoms1:
            center = atom.get_coord()
            neighbs = neighborsearch.search(center, dist, level='C')
            for x in neighbs:
                if x != atoms1:
                    neighbors.add(x)

        if len(neighbors) != 0:
            return True
        else:
            return False

    if type == "residue": # Looking at the clashes between two chains # ESTO NO LO UTILIZAMOS CREO, ASÍ QUE LO PODEMOS QUITAR
        chain1 = [atom for atom in model1.get_atoms()]
        chain2 = [atom for atom in model2.get_atoms()]

        neighborsearch = NeighborSearch(chain2) # Is it relevant which one we choose?

        neighbors = set()

        for atom in chain1:
            center = atom.get_coord()
            neighbs = neighborsearch.search(center, dist, level='R')
            #print(neighbs)
            for x in neighbs:
                if x != chain1:
                    neighbors.add(x)

        if len(neighbors) != 0:
            return True
        else:
            return False

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

# FOR THE SUPERIMPOSE: If we're mutating the DNA chains, we don't really need to
# separate them, right?

def stoichiometry_check(unicommon_chains, stoichiometry_dict):
    """
    It returns a boolean argument that indicates if the provided and actual stoichiometry are
    similar (True) or not (False).
    :param unicommon_chains:
    :return:
    """
    print("Checking stoichiometry....")

    # Comparing both dictionaries based on number of keys and length of each value

    new_unicommons = unicommon_chains.copy()

    for x in unicommon_chains.keys():
        if inputfunctions.check_type(x) == "DNA":
            del new_unicommons[x]


    c = 0
    if len(new_unicommons) == len(stoichiometry_dict):
        c += 1

    unicommon_values = list(unicommon_chains.values())
    stechio_values = list(stoichiometry_dict.values())

    for num, values in zip(stechio_values, unicommon_values):
        if num <= len(values)+1:
            c += 1

    if c >= len(unicommon_chains) + 1:
        return True
    else:
        return False

def protein_interactions(dictmodels):
    """Gets all possible interactions given a list of pdb files. It computes the distance between chains and pairs the interacting
    ones. In order to interct they need to fulfill the following condition: at least 8 residues are found at a distanc of less than ess than 8A"""
    set_models = set()
    main_dict = dict()

    for model in dictmodels.values():
        for chain in model:
            set_models.add(chain)

    for ind, chain1 in enumerate(list(set_models)):
        main_dict[chain1]={}
        for chain2 in list(set_models)[ind+1:]:
            if chain1.get_id() != chain2.get_id():
                if chain2 not in main_dict:
                    main_dict[chain1][chain2]= modelfunctions.get_distance(chain1, chain2)
    #print(main_dict)

    inter_list = []
    for k1, v1 in main_dict.items():
        for k2,v2 in v1.items():
            if v2 >= 8:
                if [k1, k2] not in inter_list:
                    inter_list.append([k1, k2])
    return(inter_list)

# PREVIOUS check_type
"""
for at in atoms:
    if at.id == "CA":
        return "Protein"
    elif at.id == "P":
        return "DNA"
"""

"""

i = 0
value_list = list(unicommon_chains.values()) # Iterate over similar chains
#print(value_list)
unicommon_chains2 = copy.deepcopy(unicommon_chains)

for key in unicommon_chains.keys():
    current_chains = value_list[i]

    for chain in current_chains:
        rmsd = modelfunctions.Superimpose(key, chain, chains=True)
        if rmsd == False:
            break
        else:
            rmsd = rmsd.rms
            clash = modelfunctions.clashes(key, chain, type="atom")
            if rmsd < 0.001 and clash == True:
                unicommon_chains2[key].remove(chain)
                #print("CCH", current_chains)
        #print(f"RESULTING CHAINS {current_chains}\n")

    for index, ch in enumerate(current_chains):
        for ch2 in current_chains:
            #print(ch2)
            if ch.id != ch2.id:
                rmsd = modelfunctions.Superimpose(ch, ch2, chains=True)
                rmsd = rmsd.rms
                clash = modelfunctions.clashes(ch, ch2, type="atom")
                #print(f"Between {ch} and {ch2} the RMSD is {rmsd} and the CLASH status is {clash}")
                if rmsd < 0.001 and clash == True:
                    if ch2 in unicommon_chains2[key]:
                        unicommon_chains2[key].remove(ch2)
    i+=1
"""

## AN ALTERNATIVE TO LOOKING FOR CLASHES: COMPUTING THE DISTANCE

def get_dist(list1, list2):
    """Given two lists of atom coordinates, it calculates the distance."""
    import math
    distance_result = (list1[0] - list2[0]) ** 2 + (list1[1] - list2[1]) ** 2 + (list1[2] - list2[2]) ** 2
    return math.sqrt(abs(distance_result))

def get_distance(chain1, chain2):
    distances=[]
    chain1_atoms=chain1.get_atoms()
    atoms1=list()
    for atom in chain1_atoms:
        if atom.id == "CA":
            atoms1.append(atom)
    chain2_atoms=chain2.get_atoms()
    atoms2=list()
    for atom in chain2_atoms:
        if atom.id == "CA":
            atoms2.append(atom)

    for combination in itertools.product(atoms1, atoms2):
        coords1=combination[0].get_coord()
        coords2=combination[1].get_coord()
        dist=get_dist(list(coords1.tolist()), list(coords2.tolist()))
        distances.append(dist)

     #print(distances)
    num_min = 0
    if len(distances) != 0:
        dist_arr = numpy.array(distances)
        for dist in dist_arr:
            if dist < 8:
                num_min += 1
        return num_min



def macrocomplex_interacting_chains(pdbfile, verbose = False):
    """
    Given a macrocomplex PDB file, it computes the interacting pairs, deletes the redundant ones
    and builds a "fake" pdb_dict (similar to the object built with the interacting chain-pair
    PDBs. This dictionary has a placeholder name (instead of the path) as key and a MODEL PDB
    object with pairs of interacting chains as a value.

    This dictionary is also employed to: Check for heteroatoms, check for DNA, Protein chains and
    changes the chain IDs according to an alphabet.

    :param pdbfile: The input PDB macrocomplex file
    :return: input_pairs_dict
    """

    # PDB parsing
    pdbname = os.listdir(pdbfile)[0] # turn to string
    structure_name = (pdbfile.split("/"))[-1]
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure(structure_name, pdbfile+pdbname)
    chains_list = [ch for ch in structure.get_chains()]


    ### Interacting pairs

    if verbose == True:
        print("Looking for interacting pairs...")

    pairs = []
    pairs2 = []  # Same as pairs but a plain list. Useful for debugging purposes.

    for idx, ch in enumerate(chains_list):
        for ch2 in chains_list[idx+1:]:
            if ch.id != ch2.id: # Avoid comparison between the same IDs (same chains)

                interaction = modelfunctions.interacting_pairs(ch, ch2, dist=8, num_contacts=8)

                print(f"{ch} - {ch2} - {interaction}")

                if interaction == True:
                    pairs.append([ch, ch2]) # Pairs that do interact

                    pairs2.append(ch) # debugging, not necessary
                    pairs2.append(ch2)

    if verbose == True:
        print(f"We have {len(pairs)} interacting pairs")


    ## Non-interacting pairs
    # Is this even possible?
    no_pair = []
    no_pair2 = [(no_pair.append(ch) for ch in chain) for chain in chains_list]

    if verbose == True:
        print(f"We have {len(no_pair)} non-interacting pairs")


    ## Eliminating redundant pairs

    pairs_without_redundant = chainfunctions.redundant_pairs(pairs, dist=1) # Careful with the distance here



    # Creating "fake" pdb_dict

    print("Creating fake pdb_dict...")

    input_pairs_dict = {}
    pdbname = pdbname.split(".")[0]

    i = 0
    for pair in pairs_without_redundant:
        str_pair = Structure.Structure(f"{i}")
        model_pair = Model.Model(f"0")
        model_pair.add(pair[0])
        model_pair.add(pair[1])
        str_pair.add(model_pair)
        input_pairs_dict[f"Structure_{pdbname}_{i}"] = str_pair
        i += 1


    # Check Prot and ACNUC chains

    dna_chain_count = 0
    prot_chain_count = 0
    #dnachainlist = []
    #protchainlist = []
    for structure in input_pairs_dict.values():
        for chain in structure.get_chains():
            if check_type(chain) == "Protein":
                prot_chain_count += 1
                #protchainlist.append(chain)
            elif check_type(chain) == "DNA":
                dna_chain_count += 1
                #dnachainlist.append(chain)

    if verbose == True:
        print(f"There are {dna_chain_count} NT chains and {prot_chain_count} Protein chains")


    ### CHECKING FOR HETEROATOMS

    hetatm_count = 0
    for str in input_pairs_dict.values():
        for chain in str.get_chains():
            for res in chain:
                if res.id[0] != " ":
                    hetatm_count += 1

    if verbose == True:
        print(f"Number of heteroatoms: {hetatm_count}")


    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    fake_pdbdict = {}

    for name, structure in input_pairs_dict.items():
        for chain in structure.get_chains():
            if type(chain.id) != int:
                new_id = alphabet.find(chain.id)
                chain.id = new_id
                fake_pdbdict[name] = structure


    return fake_pdbdict



def pdb_parser(list_pathpdbs, verbose=False):
    """
    Given a list of paths pointing to a specific PDB file, it is read, converted to a specific Bio.Structure
    object, checks the protein and DNA chains, checks for heteroatoms and changes chain ID's
    :param list_pathpdbs:
    :return:
    """

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    # Get dictionary of structures
    #
    final_name = (list_pathpdbs[0].split("/")[-1]).split("_")[0]
    str_dict = {}
    i = 0
    for str in list_pathpdbs:
        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        structure = parser.get_structure(f"Structure_{final_name}_{i}", str) # Is this a good name???
        str_dict[str] = structure  # Path+file as KEY, structure object as VALUE

        str_dict[str] = structure
        i += 1

    if verbose == True:
        print("FILE-MODEL dictionary successfully created.")

    ## create equivalence dict

    equivalence_chains = {}

    eq_str_dict = str_dict.copy()

    for str in eq_str_dict.values():
        for chains in str.get_chains():
            chains2 = chains.copy()
            new_chain = chains.copy()
            new_chain.id = alphabet.find(chains.id)
            equivalence_chains.setdefault(new_chain, chains2)


    # Check Prot and ACNUC chains

    dna_chain_count = 0
    prot_chain_count = 0
    dnachainlist = []
    protchainlist = []
    for structure in str_dict.values():
        for chain in structure[0]:
            if check_type(chain) == "Protein":
                prot_chain_count += 1
                protchainlist.append(chain)
            elif check_type(chain) == "DNA":
                dna_chain_count += 1
                dnachainlist.append(chain)

    if verbose == True:
        print(f"There are {dna_chain_count} DNA chains and {prot_chain_count} protein chains")


    ### CHECKING FOR HETEROATOMS

    hetatm_count = 0
    for str in str_dict.values():
        for chain in str[0]:
            for res in chain:
                if res.id[0] != " ":
                    hetatm_count += 1

    if verbose == True:
        print(f"Number of heteroatoms: {hetatm_count}")


    ### Changing chain IDs

    str_dict_final = {}



    """
    i = 0
    j = 0
    for file, model in str_dict.items():
        str_dict_final[file] = change_chain_id(model, i, j)
        i += 2
        j += 1
    """

    new_str_dict = str_dict.copy()

    for name, structure in new_str_dict.items():
        for chain in structure.get_chains():
            #new_chain = chain.copy()
            new_id = alphabet.find(chain.id)
            chain.id = new_id
            str_dict_final[name] = structure

    if verbose == True:
        print("Chain IDs have been converted from letters to numbers.\n")

    return [str_dict_final, dna_chain_count, prot_chain_count, equivalence_chains]


"""
for index, pair1 in enumerate(listpairs):
    #print(f"PAIR 1: {pair1}")

    for pair2 in listpairs[index+1:]:

        #print(f"PAIR 2: {pair2}")
        model_list = list(ch for ch in new_model.get_chains())
        print(f"\nCURRENT MODEL WITH CHAINS: {model_list}")


        for chain1 in model_list:
            key1 = modelfunctions.extract_key(chain1, unicommon)

            #print(f"Chain 1: {chain1}\tKey 1: {key1}")

            for idx, chain2 in enumerate(pair2):
                print(f"Comparing between {chain1} and {chain2}")

                key2 = modelfunctions.extract_key(chain2, unicommon)
                #print(f"Chain 2: {chain2}\tKey 2: {key2}")

                if key1 == key2:
                    print("\nCHAIN1 and CHAIN2", chain1, chain2)
                    if inputfunctions.check_type(chain2) == "Protein":
                        imposer = modelfunctions.Superimpose(chain1, chain2, chains=True)
                        imposer.apply(pair2[abs(idx + 1)])
                        moved_chain = pair2[abs(idx + 1)].copy()
                        clash_count = modelfunctions.clash_list(moved_chain, model_list)
                        print("CLASH COUNT", clash_count)
                        if clash_count == 0:
                            new_model.add(moved_chain)
                            print(f"chain {moved_chain} has been added")
    break
"""

"""

for idx, pair1 in enumerate(listpairs):
    print("\n")
    print("Seed pair is", pair1)
    pair1_key = ""
    pair2_key = ""

    for pair2 in listpairs[i + 1:]:
        print("\nCURRENT MODEL CHAINS", [chain for chain in new_model.get_chains()])
        print("Checking second pair:",  pair2,)

        if pair1 == pair2:
            break
        chainlist_model = (chain for chain in new_model.get_chains())
        i = 0
        j = 0

        clash_count = 0

        pair1_key = modelfunctions.extract_key(pair1[i], unicommon_prot)
        pair2_key = modelfunctions.extract_key(pair2[j], unicommon_prot)

        # If they do have a key in common (first chain and second chain

        if pair1_key == pair2_key:
            print(f"Between both pairs there is a common key, which is {pair2_key}")

            moved_chain = modelfunctions.superimpose_pair(new_model, pair2, i, j)[1]

            for chain in chainlist_model:
                clash = modelfunctions.clashes(moved_chain, chain, "atom")
                if clash == True:
                    clash_count += 1

            if clash_count == 0:
                if moved_chain not in list(new_model.get_chains()):
                    print(f"No clashes found. {moved_chain} added to the model")
                    new_model.add(moved_chain)

            break

        # If not
        else:
            print("Looking at the other chain for a common key...")
            j = 1
            pair1_key = modelfunctions.extract_key(pair1[i], unicommon_prot)
            pair2_key = modelfunctions.extract_key(pair2[j], unicommon_prot)
            if pair1_key == pair2_key:
                moved_chain = modelfunctions.superimpose_pair(pair1, pair2, i, j)[1]
                for chain in chainlist_model:
                    clash = modelfunctions.clashes(chain, moved_chain, "atom")
                    if clash == True:
                        clash_count += 1
                if clash_count == 0:
                    if pair2[j-1] not in list(new_model.get_chains()):
                        new_model.add(moved_chain)
            else:
                j = 0
                i = 1
                moved_chain = modelfunctions.superimpose_pair(pair1, pair2, i, j)[1]
                for chain in chainlist_model:
                    clash = modelfunctions.clashes(chain, moved_chain, "atom")
                    if clash == True:
                        clash_count += 1
                if clash_count == 0:
                    if pair2[j+1] not in list(new_model.get_chains()):
                        new_model.add(moved_chain)
                else:
                    j = 1
                    moved_chain = modelfunctions.superimpose_pair(pair1, pair2, i, j)[1]
                    for chain in chainlist_model:
                        clash = modelfunctions.clashes(chain, moved_chain, "atom")
                        if clash == True:
                            clash_count += 1

                    if clash_count == 0:
                        if moved_chain not in list(new_model.get_chains()):
                            new_model.add(moved_chain)
                    else:
                        continue
    break

#for chain in new_model.get_chains():
#    print("MODEL", chain)

"""

# etc

#### DIFFERENT CHAINS
# Chains that are NOT present in similar chains
# These can NOT be superimposed. I don't see how could we include them


"""
ch1 = [ch for ch in listmodels[0].get_chains() if ch not in [ch for ch in new_model.get_chains()]]

for chain in ch1:
    if ch1[0] in keys_list or ch1[1] in keys_list:  # Only introduce the ones that contain unique chains first
        new_model.add(chain)
        print(f"{chain} has been added to the model")

for idx, model in enumerate(listmodels):

    #print("Chain 1", ch1)

    for model2 in listmodels[idx+1:]:
        ch2 = [ch for ch in model2.get_chains()] # List of chains of the model that is going to be incorporated
        #print("Chain 2", ch2)
        #print(newdict_prot[ch1[0]])

        # Looking for the common chain, if there is any
        print(f"Looking for the common chains between {ch1} and {ch2}")
        try:
            if ch2[0].id in newdict_prot[ch1[0]] or ch2[1].id in newdict_prot[ch1[0].id]:
                print("Common chain found")
                common_chain1 = ch2[0]
                common_chain2 = ch1[0]
                print("Common chains", common_chain1, common_chain2)
                exit()
                rmsd = modelfunctions.Superimpose(model, model2, apply=True)
                clash = modelfunctions.clashes(model, model2, type="atom")
                # print(f"Between {key} and {chain} the RMSD is {rmsd} and the CLASH status is {clash}")
                print(rmsd, clash)
                if rmsd > 5:
                    try:
                        for chain in model2:
                            new_model.add(chain)
                        new_model_stech[key] += 1
                        print(f"###### This {chain} has been ADDED")
                    except:
                        pass
        except:
            pass




# Then use the unique-common dict to build the model
for key in unicommon_prot.keys():
    #print(f"\nFor unique chain {key}")
    try:
        new_model.add(key)
    except:
        pass
    new_model_stech[key] = 1
    current_chains = value_list[i]
    j = 0
    for chain in current_chains:
        key = inputfunctions.mutate_dna_chain(key)
        chain = inputfunctions.mutate_dna_chain(chain)
        rmsd = modelfunctions.Superimpose(key, chain, apply=False, chains=True)
        clash = modelfunctions.clashes(key, chain, type="atom", dist=0.01)
        print(f"Between {key} and {chain} the RMSD is {rmsd} and the CLASH status is {clash}")
        if rmsd > 0.001 and clash == False:
            if chain.id not in [ch.id for ch in new_model.get_chains()]:
                try:
                    new_model.add(chain)
                    new_model_stech[key] += 1
                    print(f"###### This {chain} has been ADDED")
                except:
                    pass

    for index, ch in enumerate(current_chains):
        for ch2 in current_chains[index+1:]:
            print(ch, ch2)
            rmsd = modelfunctions.Superimpose(ch, ch2, apply=False, chains=True)
            clash = modelfunctions.clashes(ch, ch2, type="atom", dist=0.01)
            print(f"Between {ch} and {ch2} the RMSD is {rmsd} and the CLASH status is {clash}")
            if rmsd > 0.001 and clash == False:
                if ch.id not in [ch.id for ch in new_model.get_chains()]: # We don't need this?
                    try: # Added the try-except block again because it wasn't working with the last model
                        new_model.add(ch)
                        new_model_stech[key] += 1
                        print(f"###### This {ch} has been ADDED")
                    except:
                        pass

    i += 1

"""

"""
if options.verbose:
    print(f"Looking for the different chains")

different_chains = []
for model in total_models.values():
    for chain in model:
        if chain not in similar_chains_prot.values() and chain not in similar_chains_prot.keys():
            different_chains.append(chain)

if options.verbose:
    print("These are the DIFFERENT CHAINS:")
    print(different_chains, "\n")

"""

"""

# First add the different chains
for chain in different_chains:
    new_model.add(chain) # No clash check??


unique_keys=set()
for index, key1 in enumerate(keys_list):
    for key2 in keys_list[index+1:]:
        clash = functions.clashes(key1, key2)
        if clash == False:
            unique_keys.add(key1)
print(unique_keys)

unique_dict={}
for key in unique_keys:
    unique_dict[key]=newdict_prot[key]
print(unique_dict.items())

value_list_prot = list(unique_dict.values())
"""