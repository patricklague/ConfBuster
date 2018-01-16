![N|Solid](http://132.203.89.236/ConfBuster/confbuster.png)
# ConfBuster - Open source tools for macrocycles conformational search and analysis

&nbsp;

# Dependencies
**CORE**
- Python 2.7 (https://www.python.org)
- NetworkX (tested with 1.11) (https://networkx.github.io)
- Pymol (≥ 1.8) (https://sourceforge.net/projects/pymol)
- Open Babel (≥ **2.4.1**) (http://openbabel.org/wiki/Main_Page)
The version of Open Babel should **imperatively be 2.4.1**.
As indicated in the instruction, Eigen version 2 is requiered to use properly Open Babel (https://openbabel.org/docs/dev/Installation/install.html).
Link to download: https://sourceforge.net/projects/openbabel/files/

For a clean and functional installation:
 ```{sh}
In openbabel-2.4.1 directory:
$ mkdir build
$ cd build
$ cmake ../ -DCMAKE_INSTALL_PREFIX=/path/to/where/you/want
$ make -j4
$ make install
(This is essential for Open Babel to work properly)
```
**OPTIONAL (VISUALISATION)**
- R (≥ 3.0.0) (https://cran.r-project.org/index.html)
- R package ComplexHeatmap (≥ 1.14) (from Bioconductor release 3.5 (bioconductor.org))

**R shell**
 ```{r}
source("http://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
```
- R package circlize (≥ 0.4.0)

**R shell**
 ```{r}
install.packages("circlize")
```
# Installation
The most simple way to install ConfBuster is by downloading the official release in the GitHub section "release". Alternatively, one may clone the ConfBuster repository using git :
 ```sh
$ git clone https://github.com/patricklague/ConfBuster.git
  ```
At your choice, put the scripts in a directory present in your $PATH OR add the directory where the scripts are to your $PATH:
   ```sh
 export PATH=$PATH:/full-path-here/ConfBuster/
   ```
 PyMOL and Open Babel should also be included in your $PATH. If you want to install them locally, you should change the headers of the ConfBuster scripts to add their correct absolute paths.
 
# ConfBuster tools
- **ConfBuster-Single-Molecule-Minimization**
Performs a simple minimization of the given molecule (recommended before).
 ```sh
 -i input filename [mandatory]
 -o output name prefix [default: replace input file]
 ```
- **ConfBuster-Rotamer-Search**
Identify rotational isomers of a molecule.
 ```sh
-i input filename in mol2 [mandatory] 
-g number of generations [default: 100] 
-e the energy cutoff used to discriminate conformations in units of kcal/mol [default: 50]
-d output directory name [default: use the prefix of the input filename] 
-f format of outputted molecules [xyz or default: mol2]
 ```
- **ConfBuster-Macrocycle-Linear-Sampling**
Performs a conformational search of a cyclic molecule.
 ```sh
-i input filename [mandatory]
-r rmsd cutoff in Angstrom [default: 0.5]
-n for each cleaving point, number of rotamer searches performed [default: 5] 
-N for each cleaving point, number of molecules extracted from each rotamer search [default: 5] 
-o output directory name [default: prefix of the input filename]
```

- **ConfBuster-Analysis**
Perform post-analyses to visualize a clustering based on RMSD values between the conformations.
 ```sh 
 -i directory name of the search results [mandatory] 
 -r rmsd cut off [default: none]
 -n number of conformations to include in the analysis [default: all]
 -e mid-point value of the energy color scale [default: 0]
  ```
# Tutorial

Several examples including all the command lines and required files to run macrocyle conformational searches are included with the distribution, in the **examples** folder and in the **examples/Instructions.pdf** file.The user can follow this tutorial by using the molecule 1w96.pdb available in the directory **examples/1w96**.
The initial structure is:

![minipic](http://132.203.89.236/ConfBuster/macro-1w96-1.png)

**Step 1**
Perform a single point optimisation on macro-1w96.pdb:
 ```sh
$ ConfBuster-Single-Molecule-Minimization.py -i macro-1w96.pdb
  ```
**Step 2**
Start Macrocyclic linear sampling:
 ```sh
 $ ConfBuster-Macrocycle-Linear-Sampling.py -i macro-1w96.mol2 -n 5 -N 5 -r 0.5
  ```
**Step 3 (Optional)**
View and analyse results:
 ```sh
$ cd macro-1w96
$ pymol Follow-macro-1w96.py
  ```
 or alternatively, directly in Pymol:
   ```sh
  run Follow-macro-1w96.py
  ```
 The result should looks like:

![minipic](http://132.203.89.236/ConfBuster/conformational-range-3.png)

Lowest energy conformation is in green (conf-16.mol2 in the matrix below). PDB structure is in cyan. RMSD between conf-16.mol2 and PDB structure is of 0.405 Angstrom over all atoms. All sampled conformations are in black thin lines to show the conformational range of the macrocycle.
 
Using the ConfBuster-Analysis.py module ($ ConfBuster-Analysis.py -i macro-1w96 -n 20) allows to compare the n best results with a reference or the initial structure. Results will be shown as a 2d RMSD matrix of n conformations. On the right, an energy bar (from green to purple) between the best and worst energy is displayed.

![minipic](http://132.203.89.236/ConfBuster/Heatmap_20.png)

**Other examples are provided in the directory "Examples"**

# How to cite ConfBuster
Barbeau, X., Vincent, A.T. & Lagüe, P., (2018). ConfBuster: Open-Source Tools for Macrocycle Conformational Search and Analysis. Journal of Open Research Software. 6(1), p.1. DOI: http://doi.org/10.5334/jors.189
# License
>ConfBuster Suite
>Copyright (C) 2017  Xavier Barbeau, Antony T. Vincent and Patrick Lagüe

>This program is free software: you can redistribute it and/or modify
>it under the terms of the GNU General Public License as published by
>the Free Software Foundation, either version 3 of the License, or
>(at your option) any later version.

>This program is distributed in the hope that it will be useful,
>but WITHOUT ANY WARRANTY; without even the implied warranty of
>MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
>GNU General Public License for more details.

>You should have received a copy of the GNU General Public License
>along with this program.  If not, see http://www.gnu.org/licenses/.


