#!/usr/bin/python3

import parser as args # argparser script
import inputfunctions, chainfunctions, modelfunctions


# Import the argparser object with the necessary arguments
options = args.output_argparser()


################
# INPUT PDBs
# Processing the input PDB files
################

verb = False
if options.verbose: # If the verbose flag is activated
    verb = True
    print("\tReading PDBs\n")

pdblist = args.check_infile_names(options.infile)  # List of PATHS to the pdb files correctly named



### MACROCOMPLEX INPUT
# If only 1 PDB, containing a full model is provided
if options.macrocomplex and len(pdblist) == 1:  # -m flag activated, only 1 input file
    pdb_macro = pdblist[0]  # Convert to string


    if options.verbose:
        print("\n\nSplitting and re-combining input chains into a pdb_dict...")

    #Parsing the macrocomplex into Structure with interacting pairs and return a file list to pdblist type dictionary
    strdict = inputfunctions.macrocomplex_parser(pdb_macro, interact_dist=options.contact_distance, interact_contacts=options.contact_num, redundant_dist=1.9, verbose=verb)

### INTERACTING PAIRS INPUT
# A list of PDB chain pairs is provided
else:

    # If -m indicated
    if options.macrocomplex:
        if options.verbose:
            sys.stderr.write("-m flag wrongly indicated. The input is a list of PDB files")
            sys.stderr.flush()


    # Convert input file list to pdblist type dictionary
    strdict = inputfunctions.pdb_parser(pdblist, verbose=True)
    #print("STR DICT", strdict)



################
# PROCESSING CHAINS
# Re-format list of interacting chains to properly construct the model
################


##### EQUIVALENT CHIAN ID
### Obtain the equivalence between chain id in numbers and letters
##### PDB DICT
### Obtain the same dict as strdict but with the chain id changed to numbers
equivalent_chains, pdb_dict = inputfunctions.chain_processing(strdict)

##### DIFFERENT CHAINS
### Check the number of "different" chains of a macrocomplex (debugging purposes)
if options.macrocomplex:
    diff_chains = set()
    for str in pdb_dict.values():
        for chain in str.get_chains():
            diff_chains.add(chain)

    print(f"We have {len(diff_chains)} different chains.")



##### SIMILAR SEQUENCE CHAINS
### Obtain high-sequence (>95%) similarity chains from the interacting pairs
similar_chains_prot = chainfunctions.similar_chains(pdb_dict, "Protein", verbose=verb) # Protein
similar_chains_dna = chainfunctions.similar_chains(pdb_dict, "DNA", verbose=verb) # Nucleic acids

# Merge both dictionaries
similar_chains = similar_chains_prot
similar_chains.update(similar_chains_dna)

if options.verbose:
    print(f"Similar chains: ")
    for x, y in similar_chains.items():
        print(f"{x}\t{y}")
    print("\n")

##### UNIQUE-COMMON CHAINS
### Build specific type of dictionary, necessary for model building

unicommon = chainfunctions.unique_common_chains(similar_chains, verbose=verb)


if options.verbose:
    print("UNIQUE-COMMON Chains:")
    for x, y in sorted(unicommon.items(), key=lambda x:x[0]):
        print(f"UNIQUE CHAIN: {x}\tSIMILAR CHAINS: {y}")
    print("\n")

################
# STOICHIOMETRY CHECK
# Read stoichiometry file
################

if options.stoichiometry:
    if options.verbose:
        print(f"Reading stechiometry...")

    # Fetch stoichiometry from input

    stoichiometry_input = inputfunctions.stoichiometry_parser(options.stoichiometry)

    print(f"The provided stechiometry is:")
    for id, count in stoichiometry_input.items():
        print(f"{id}\t{count}")
    print("\n")


    # Complete unicommons

    diff_chains = set()
    for str in pdb_dict.values():
        for chain in str.get_chains():
            diff_chains.add(chain)

    unicommon = chainfunctions.unicommon_completer(unicommon, stoichiometry_input, diff_chains)


    print("MODIFIED UNICOMMONS")

    for x, y in sorted(unicommon.items(), key=lambda x:x[0]):
        print(x, y)

else:
    print("No stoichiometry provided... The model will add as many input chains as possible.")


## model construction

if options.stoichiometry:
    final_model = modelfunctions.model_construction(unicommon, pdb_dict, equivalent_chains, stoichiometry_input=stoichiometry_input, verbose=verb, forcing=True)
else:
    final_model = modelfunctions.model_construction(unicommon, pdb_dict, equivalent_chains, verbose=verb, forcing=True)


modelfunctions.save_model(final_model, pdblist, outdir=options.outdir, verbose=verb, macrocomplex=options.macrocomplex)
