
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
  - [3. The pytprot pipeline](#3-the-pytprot-pipeline)
  - [4. How to use *pytprot*: A quick tutorial](#4-how-to-use-pytprot-a-quick-tutorial)

## 0. Pre-requisites

In order to run **pytprot**, it is compulsory to have a Python3 interpreter, and strongly recommended to have any biomolecule visualization software, such as Chimera or ICM. This module relies heavily on the **biopython** Python module, so it is also strongly advised to download said module.

## 1. Installation

#### 1.3. Through PyPi (Recommended)

We have added our module into PyPi to make it even easier to install. In this case, you only have to execute the following command on your terminal:

```
pip3 install pytprot

```

With, this, the package and its dependencies will be automatically installed within your PYTHONPATH, where all the pytprot folder, with the examples and the corresponding scripts, will be installed.

#### 1.2. Through github

First, you have to copy the GitHub repository to a given folder, move into it, and execute the `setup.py` file.

```
git clone https://github.com/mariasr3/pytprot.git
cd pytprot
sudo python3 setup.py install
```

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


## 3. The pytprot pipeline



## 4. How to use *pytprot*: A quick tutorial


To know which arguments needs the program to run or which ones can be determined by the user, type in the terminal: 
python3 main.py -h (or --help)
The arguments are:
- **-i INFILE or --input-directory INFILE**: Required. The user must provide a directory which contains the input files. 
- **-s or --stoichiometry**: Optional. The user can handle an input file with the chain stoichiometry. 
- **-f or --force**: Optional. It forces the creation of an output directory. If it already exists, the contents are overwritten
- **-o OUTDIR or --output-directory OUTDIR**: Required. The user must provide an ouput directory where the output files will be stored.
- **-v or --verbose**: Optional. It prints in the terminal some information about the process. 
- **-m or --macrocomplex**: Optional. If the input file is a macrocomples, the user should use this option to indicate it. 
- **-opt or --optimization**: Optional. It refines the model through MODELLER
- **-d or --contact_dist**: Optional. By default the distance between two chains to consider they are interacting is set to 12. 
- **-cn or --contact_num**: Optitonal. By default there needs to be a minimum of 8 residues to interact at a distance d between two chains to interact.
