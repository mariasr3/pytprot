
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![GitHub license](https://img.shields.io/github/license/Naereen/StrapDown.js.svg)](https://github.com/mariasr3/pytprot/blob/main/LICENSE.txt)


# pyprot: PPI-based multicomplex builder

**pytprot** is a web-app and command-line Python3 package to build protein multicomplexes from binary chain interactions.

*Joana Llauradó, Maria Sopena, Othmane Hayoun*

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
