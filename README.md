![N|Solid](http://132.203.89.236/ConfBuster/confbuster.png)
# ConfBuster - Open source tools for macrocycles conformational search and analysis

&nbsp;

# Dependencies
**CORE**
- Python 2.7 (https://www.python.org)
- R (≥ 3.0.0) (https://cran.r-project.org/index.html)
- NetworkX (tested with 1.11) (https://networkx.github.io)
- Pymol (≥ 1.8) (https://sourceforge.net/projects/pymol)
- OpenBabel (≥ **2.4.1**) (http://openbabel.org/wiki/Main_Page)
The version of OpenBabel should **imperatively be 2.4.1**.
Link to download it: https://sourceforge.net/projects/openbabel/files/
For a clean and functional installation:
 ```{sh}
In openbabel-2.4.1 directory:
$ mkdir build
$ cd build
$ cmake ../ -DCMAKE_INSTALL_PREFIX=/path/to/where/you/want
$ make -j4
$ make install
(This is essential for OpenBabel to work)
```
**OPTIONAL (VISUALIZATION)**
- R package ComplexHeatmap (≥ 3.5) (from bioconductor.org)

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
$ git clone https://github.com/lague_patrick/ConfBuster.git
  ```
At your choice, put the scripts in a directory present in your $PATH OR add the directory where the scripts are to your $PATH:
   ```sh
 export PATH=$PATH:/full-path-here/ConfBuster/
   ```
 PyMOL and OpenBabel should also be recognized in your $PATH. If you want to install them locally, you should change the headers of the ConfBuster scripts to add their correct respective paths.
 
# ConfBuster tools
- **ConfBuster-Single-Molecule-Minimization**
Perform a simple minimization of the given molecule (recommended before all analysis).
 ```sh
 -i inputfile [mandatory]
 -o output name prefix [auto: replaces initial file]
  ```

  - **ConfBuster-Macrocycle-Linear-Sampling**
 Find the conformation of a cyclic molecule.
  ```sh
-i inputfile [mandatory]
-r rms deviation cutoff [0.5]
-n for each cliving point, number of rotamer searches performed [5] 
-N for each cliving point, number of molecules extracted from each rotamer search [5] 
-o output directory name
```
- **ConfBuster-Rotamer-Search**
Identify rotational isomers of a molecule.
 ```sh
-i inputfile in mol2 [mandatory] 
-g # of generations [100] 
-e accepted delta E in kcal/mol [50] 
-d output directory [automatic]
-f format of outputted molecules [xyz,mol2(default)]
 ```
 - **ConfBuster-Analysis**
 Perform post-analyses to visualize a clustering based on RMSD values between the conformations.
 ```sh 
 -i inputfile in mol2 [mandatory] 
 -r rmsd cut off 
 -n use n conformation
 -e adjust e_midpoint variable
  ```
# Tutorial
The can follow this tutorial by using the molecule 1w96.pdb available in the directory "example".
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
   ```python
  run Follow-macro-1w96.py
  ```
 The result should looks like:

![minipic](http://132.203.89.236/ConfBuster/conformational-range-3.png)

Lowest energy conformation is in green (conf-16.mol2 in the matrix below). PDB structure is in cyan. RMSD between conf-16.mol2 and PDB structure is of 0.405 over all atoms. All sampled conformations are in black thin lines to show the conformational range of the macrocycle.
 
Using the ConfBuster-Analysis.py module allows to compare the n best results with a reference or the initial structure. Results will be shown as a 2d RMSD matrix of n conformations. On the right, an energy bar (from green to purple) between the best and worst energy is displayed.

![minipic](http://132.203.89.236/ConfBuster/Heatmap_20.png)

**Other examples are provided in the directory "Examples"**

# How to cite ConfBuster
A publication describing ConfBuster is submitted. This section will be updated as soon as possible. However, if using ConfBuster, we recommend to cite the present GitHub page.
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


