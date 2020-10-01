# RACER
Rapid Coarse-grained Epitope TCR Model

## An implementation of the RACER model for TCR recognition


## Installation
* Clone RACER repository
```
git clone https://github.com/XingchengLin/RACER.git
```

* Download and install python using conda: https://www.anaconda.com/products/individual
* Download and install Modeller using conda: https://anaconda.org/salilab/modeller
* Download and install Biopython using conda: https://biopython.org/wiki/Packages

## Example demo
* Optimize a energy model given available PDB of 3QIB.pdb : https://www.rcsb.org/structure/3QIB
* Calculate the binding energies of native sequence (strong binders) of 3QIB.pdb and testing sequences (weak binders) generated by Lanzarotti et. al. DOI: 10.1016/j.molimm.2017.12.019

* Step 1 Clean the PDB files
* Clean the PDB files so that they contain only one TCR-p-MHC pair with defined Chain IDs for TCR as well as p-MHC. An example is given in the folder data/. The native.pdb is cleaned from 3QIB.pdb. The testBinders were built using Modeller based on the template of native.pdb, with replaced peptide sequences.

* Step 2 Optimize the energy model
* This step includes 1. Processing the PDB files so that they follow the format of AWSEM model developed by Wolynes group DOI: 10.1021/jp212541y; 2. Generate 1000 decoy peptide sequences from the strong binder; 3. Optimize a force field by maximize the Z score between the binding energies of strong and decoy sequences.
```
bash cmd.preprocessing.sh 3qib C D 782 794
```
* Note: "C" and "D" are the chain ID of TCR alpha and beta chains. 782 and 794 are the starting and ending residue IDs of the presented peptide.






