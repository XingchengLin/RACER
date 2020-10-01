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
* Optimize an energy model given available PDB of 3QIB.pdb : https://www.rcsb.org/structure/3QIB
* Calculate the binding energies of native sequence (strong binders) of 3QIB.pdb and testing sequences (weak binders) generated by Lanzarotti et. al. DOI: 10.1016/j.molimm.2017.12.019

* Step 1 Clean the PDB files
* Clean the PDB files so that they contain only one TCR-p-MHC pair with defined Chain IDs for TCR as well as p-MHC. An example is given in the folder data/. The native.pdb is cleaned from 3QIB.pdb. The testBinders were built using Modeller based on the template of native.pdb, with replaced peptide sequences.

* Step 2 Optimize the energy model
* This step includes 1. Processing the PDB files so that they follow the format of AWSEM model developed by Wolynes group DOI: 10.1021/jp212541y; 2. Generate 1000 decoy peptide sequences from the strong binder; 3. Evaluate a "Phi" map between the TCR and the peptide, based on their contact probability; 4. Optimize a force field by maximizing the Z score between the binding energies of strong and decoy sequences (also c.f. DOI: 10.1002/ijch.201300145)
```
bash cmd.preprocessing.sh 3qib C D 782 794
bash cmd.optimize.sh
```
* Note: "C" and "D" are the chain IDs of TCR alpha and beta chains. 782 and 794 are the starting and ending residue IDs of the presented peptide.
* Note: For now, one needs to hard-set the cutoff for noise-filtering of the eigenvalues of the interaction matrix. Please make sure the file gammas/randomized_decoy/01022019/direct_contact/proteins_list_phi_pairwise_contact_well4.5_6.5_5.0_10_lamb_filtered has no zero terms.

* Step 3 Use the optimized energy model to evaluate the effective binding energies of strong/weak binders
* This step uses the optimized energy model (a 20 by 20 matrix for different amino acid types) to evaluate the binding energies of the strong (native) and weak (testBinder) binders. A lower binding energy corresponds to a stronger binding affinity
```
bash cmd.evaluate_bindingE.sh 3qib C D
```
* All the codes in this demo can be executed by one command in the folder molecular_demo (one check button was added to remind users that they need to check the chain IDs of their rebuilt testing structures):
```
bash cmd.sh
```
The final results (binding energies) are reported in the folder evaluated_binding_E/


## Reference:
* Rapid Assessment of T-Cell Receptor Specificity of the Immune Repertoire, Xingcheng Lin, Jason T. George, Nicholas P. Schafer, Kevin Ng Chau, Cecilia Clementi, José N. Onuchic, Herbert Levine (https://www.biorxiv.org/content/10.1101/2020.04.06.028415v2)
