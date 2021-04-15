
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/mariasr3/pytprot/blob/main/LICENSE.txt)

# pyprot: PPI-based multicomplex builder

**pytprot** is a web-app and command-line Python3 package to build protein multicomplexes from binary chain interactions.

*Joana Llauradó, Maria Sopena, Othmane Hayoun*



- [pyprot: PPI-based multicomplex builder](#pyprot-ppi-based-multicomplex-builder)
  - [0. Pre-requisites](#0-pre-requisites)
  - [1. Installation](#1-installation)
      - [1.3. Through PyPi (Recommended)](#13-through-pypi-recommended)
      - [1.2. Through github](#12-through-github)
  - [2. Module structure](#2-module-structure)
    - [2.1. Script detail: The different pytprot submodules](#21-script-detail-the-different-pytprot-submodules)
  - [3. The pytprot pipeline](#3-the-pytprot-pipeline)
  - [4. How to use *pytprot*: A quick tutorial](#4-how-to-use-pytprot-a-quick-tutorial)
    - [4.1. Running pytprot from command line](#41-running-pytprot-from-command-line)
      - [4.1.1. Example 1.a: 1GZX, interacting chains input](#411-example-1a-1gzx-interacting-chains-input)
      - [4.1.2. Example 2.a: 6GMH, multicomplex input](#412-example-2a-6gmh-multicomplex-input)
    - [4.2. Running pytprot from Dash](#42-running-pytprot-from-dash)

## 0. Pre-requisites

In order to run **pytprot**, it is compulsory to have a Python3 interpreter, and strongly recommended to have any biomolecule visualization software, such as Chimera or ICM. This module relies heavily on the **biopython** Python module, so it is also strongly advised to download said module.

## 1. Installation

#### 1.3. Through PyPi (Recommended)

We have added our module into PyPi to make it even easier to install. In this case, you only have to execute the following command on your terminal:

```
pip3 install pytprot

```

For Unix users, _pytprot_ can run without specifically invoking the Python3 interpreter with:

```
chmod +x ./path/to/pytprot/pytprot.py

```

Once this is done, we can just run _pytprot_ from any directory within the system with:

```
pytprot.py -i input/path -o output/path 

```

With, this, the package and its dependencies will be automatically installed within your PYTHONPATH, where all the pytprot folder, with the examples and the corresponding scripts, will be installed.

#### 1.2. Through github

First, you have to copy the GitHub repository to a given folder, move into it, and execute the `setup.py` file.

```
git clone https://github.com/mariasr3/pytprot.git
cd pytprot
sudo python3 setup.py install

```
Using the `chmod` command, as explained in Section 1.1., will allow for running _pytprot_ from any directory.



## 2. Module structure

The `pytprot` module consists of the following folders:

* dash_example_files: Contains example inputs for when *pytprot* is run from Dash.
* pyptrot: Contains the different Python scripts that allow the module to run properly. These will be explained later-on.
* pytprot.egg-info: Information folder formed after the compilation of the module.
* example_files: Contains example inputs for when *pytprot* is run from command-line.

And the following files:
* LICENSE.txt: Contains the MIT license for the distribution of the *pytprot* module source code.
* README.md: MarkDown document that contains the Tutorial, and some biological information background.
* setup.py: Python script that contains the information needed in order to install the package.

### 2.1. Script detail: The different pytprot submodules

Within the `pytprot` folder, we find the different scripts that, altogether form the whole package. These are:

- **__init__.py**: The *init* script that indicates which scripts from this same folder have to be loaded.
- **pytprot.py**: The main module of the script, where the actual scripting pipeline is layed-down. This includes the checking for correct inputs, the processing of the input interacting chains or macrocomplex, the building of the macrocomplex and saving the structure. All the functions employed in this script are imported from the other scripts within the folder.
- **inputfunctions.py**: The module with which the inputs will be processed in order to end up having a common Python dictionary object.
- **chainfunctions.py**: The module with which the similar chains are found.
- **modelfunctions.py**: The module with which the actual model is built.
- **parser.py**: The *argparser* module script. Allocates the different command-line arguments.
- **pytprot_dash.py**: This script allows for opening the Dash pytprot app. The pipeline is similar to that of the main script, and it also references the different function scripts.


Nevertheless, a more detailed explanation of all the submodels, and their corresponding functions, is detailed on the Theory.pdf.



## 3. The pytprot pipeline

The *pytprot* module follows a specific pipeline in order to build the models.

1. **Input recognition**: In order for it to correctly work, *pytprot* only allows a specific type of input files. It has to be a folder with one or more PDB files. If there is only one PDB file, it is understood as a macrocomplex, independently of its filename. If there is a list of PDB files, these are understood to be interacting chain pairs, which have to be named as 1gzx_A_B.pbd or 1gzx_A_ED.pdb, referring to the structure name (1gzx) and the pair of chains (A and B), or the pair of interacting protein-double-DNA chains (A and ED). The stoichiometry input file is also read at this point. This **.txt** must also have a specific format.
 
2. **Chain processing**: Given a macrocomplex, it is broken into single chains in order to obtain the interacting chains¹, based on the distance at which there is a contact between chains and the number of contacts. At this point, the program also checks and deletes redundant interactions². If a set of interacting chains is already provided, this step is not needed. After this, the program searches for similar chains, that is, chains with a sequence identity higher than 95%, as these will be the ones that will allow for superimposition.

3. **Model building**: If a stoichiometry file is provided, it checks if it is actually correct or not. If it is not, the program will proceed with the model building without the stoichiometry, and a warning message is raised. Then, it takes into account the information from the similar chains, all the provided chains, the stoichiometry (if provided) in order to build the model. This is an iterative process that adds pairs to the model if one chain of the pair finds a similar chain within the model. Then, it superimposes both pairs, obtains a rotran matrix, and applies it to the interacting chain pair to the one superimposed. If there are no clashes at 1.9 Å with any chain of the model, it is added. This is repeated until all the provided interacting chain pairs are checked.

The model is then saved with a specific filename that includes the model name, the number of chains and the timestamp.




(1) Here, **interacting chains** refers to a distance of 12 Å and 8 number of contacts to consider two chains to be interacting, by default, although this can be changed.
(2) **Redundant chains** are pairs of interacting chains that have one common similar chain, and the other pair of chains produce clashes between them, at a 1.9 Å distance.



## 4. How to use *pytprot*: A quick tutorial

### 4.1. Running pytprot from command line

To run the module, you need to execute the main script (**pytprot.py**), located in the *pytprot* folder, found after installation. In order to view the different command-line arguments, we can execute:

```
python3 /path/to/pytprot.py -h (or --help)

```

Which will display the following arguments:

- **-i INFILE or --input-directory INFILE**: Required. The user must provide a directory which contains the input files. **NOTE**: Said input files need to follow the correct naming convention. If not, these will not be accepted by the program.
- **-s or --stoichiometry**: Optional. The user can handle an input **.txt file** with the model stoichiometry. It also has to follow a specific structure, as shown in the sample stoichiometry file in the package directory.
- **-f or --force**: Optional. It forces the creation of an output directory. If it already exists, the contents of said directory are overwritten.
- **-o OUTDIR or --output-directory OUTDIR**: Required. The user must provide an ouput directory where the output files will be stored. If the output directory already exists, and the "force" flag has not been activated, the program will terminate.
- **-v or --verbose**: Optional. It prints in the terminal some information about the process. 
- **-m or --macrocomplex**: Optional. If the input file is a macrocomples, the user should use this option to indicate it. 
- **-d or --contact_dist**: Optional. Only employed if the macrocomplex flag is active. By default the distance between the CA backbone of two chains to consider they are interacting is set to 12 Å. 
- **-cn or --contact_num**: Optitonal. Only employed if the macrocomplex flag is active. By default there needs to be a minimum of 8 residues to interact at a distance *d* between two chains to interact.

#### 4.1.1. Example 1: 1GZX, interacting chains input

In order to show the basic functioning of the program, we will use the [1gzx](https://www.rcsb.org/structure/1GZX) structure, which is the T-state haemoglobin, which is a hetero-4-mer. As an input, we have three .pdb files with interacting chains (1gzx_A_B.pdb, 1gzx_A_C.pdb, 1gzx_A_D.pdb) provided by our teachers, within a folder named "1gzx", and a stoichiometry file (1gzx_stoichiometry.txt). We will assume that pytprot is being executed from the installation folder named *pytprot*.

In order to run pytprot, we will need to open the terminal and write:

```
python3 ./pytprot/pytprot.py -i ./example_files/1gzx -o ./example_output_files -s ./example_files/1gzx_stoichiometry.txt -f -v 

```
Here, the PDB files allocated in the 1gzx folder will be read, and the final built model will be stored in the "example_output_files" folder. As the forcing flag (**-f**) is activated, the output will be stored in the output folder, even if it already exists. As the **-v** verbose flag is activated, the progress of the program will be printed on command-line. After the program runs, the program will exit, and the output file will be available for analysis. 

Note that, if the stoichiometry file is not available, we can just run it without the stoichiometry flag:

```
python3 ./pytprot/pytprot.py -i ./example_files/1gzx -o ./example_output_files -f -v 

```

This process is identical, even for larger proteins that also contain DNA.

**Note**: All these input and output files for 1gzx are included in the package in order to show the correct naming format for the PDB files and the specific format of the stoichiometry file.


#### 4.1.2. Example 2: 6GMH, multicomplex input

In order to also show how the program can work with macrocomplexes, we will use [**6gmh**](https://www.rcsb.org/structure/6gmh), which is the structure for the Activated Transcription Complex Pol II, and a hetero-20-mer, that we have retrieved from the PDB. Here, we will only use one file, "6gmh.pdb", in a folder named "6gmh", and a stoichiometry file, "6gmh_stoichiometry.txt". We will assume that pytprot is being executed from the installation folder named *pytprot*. 

To run it:

```
python3 ./pytprot/pytprot.py -i ./example_files/6gmh -o ./example_output_files -s ./example_files/6gmh_stoichiometry.txt -f -v -m -d 8 -nc 5

```

We have indicated the input, output file and the **-f** flag the same way as before. As this is a macrocomplex, we will need to indicate the **-m** flag. Within the macrocomplex processsing, we can indicate a specific contact distance (Å) with **-d** and a specific number of contacts with **-nc**. For this example, we will assume two chains to be interacting if there are 5 or more contacts at 8 Å between their CA backbone. As the **-v** flag is set, when running the program, we will see the processing information on the Terminal. After it finishes running, the output model will be saved with a timestamp and the number of chains that the model has, in the provided output folder.



### 4.2. Running pytprot from Dash

The other option to run *pytprot* is through a local Dash app, which provides a more user-friendly GUI. In order to launch the the app, we need to execute in the terminal:

```

python3 /path/to/pytprot/pytprot_dash.py

```

This will prompt, on the command line, a link that, when clicked, will open a tab in your default browser, that looks like this:

![pytprot Dash GUI](pytprot.png)

Using the app is quite simple, we will first need to upload the corresponding PDB and stoichiometry files (*It is advisable to upload them all at once*). It automatically detects, depending on the number of files, if the input is a set of interacting chain-pairs or a macrocomplex. If it is the latter case, we can also indicate the contact distance or the number of contacts in the **Multicomplex assembly parameters**. 

**NOTE:** The contact distance and number of contact parameters, if different from the default ones, must be indicated **before** uploading the files.

The **Upload information** window shows some basic information regarding the input files.

After processing the model, its filename will appear under the **Models built** section. In order to access it, we will need to retrieve it from the "app_uploaded_files/built_models" folder located in the current working directory, or through the link that Dash returns after building the model.

