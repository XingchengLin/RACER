# RACER
Rapid Coarse-grained Epitope TCR Model

## An implementation of the RACER model for TCR recognition


## Installation
* Clone RACER repository
```
git clone https://github.com/XingchengLin/RACER.git
```

## Molecular Demo

## System requirement:
Note: The code showing here is built for a Linux system. For a Mac user, please use the one line command:
```
cmd_4mac.sh
```


* Download and install python using conda: https://www.anaconda.com/products/individual
* Download and install Modeller using conda: https://anaconda.org/salilab/modeller
* Note: One needs to obtain a Modeller academic license in order to run Modeller: https://salilab.org/modeller/registration.html
* Download and install Biopython using conda: https://biopython.org/wiki/Packages

## Example
* Optimize an energy model given available PDB of 3QIB.pdb : https://www.rcsb.org/structure/3QIB
* In this example, we will generate 1000 decoy sequences for optimization, and calculate the binding energies of native sequence (strong binders) of 3QIB.pdb and testing sequences (weak binders) generated by Lanzarotti et. al. DOI: 10.1016/j.molimm.2017.12.019

* Step 1 Clean the PDB files
* Clean the PDB files so that they contain only one TCR-p-MHC pair with defined Chain IDs for TCR as well as p-MHC. An example is given in the folder data/. The native.pdb is cleaned from 3QIB.pdb. The testBinders were built using Modeller based on the template of native.pdb, with replaced peptide sequences.

* Step 2 Optimize the energy model
* This step includes 1. Processing the PDB files so that they follow the format of the AWSEM model (DOI: 10.1021/jp212541y) developed by Wolynes group; 2. Generate 1000 decoy peptide sequences from the strong binder; 3. Evaluate a "Phi" map between the TCR and the peptide, based on their contact probability; 4. Optimize a force field by maximizing the Z score between the binding energies of strong and decoy sequences 
```
bash cmd.preprocessing.sh 3qib C D 782 794
bash cmd.optimize.sh
```
* Note: "C" and "D" are the chain IDs of TCR alpha and beta chains. 782 and 794 are the starting and ending residue IDs of the presented peptide.
* Note: For now, one needs to hard-set the cutoff for noise-filtering of the eigenvalues of the interaction matrix. Please make sure the file: ./gammas/randomized_decoy/01022019/direct_contact/proteins_list_phi_pairwise_contact_well4.5_6.5_5.0_10_lamb_filtered has no zero terms. Here, the cutoff was set as 74 in the file cmd.optimize.sh

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

* Explanation for the output:
* epitopeE.txt -- Binding energies of the strong binder
* non-epitopeE.txt -- Binding energies of the weak binders



## Statistical Demo

* Note: The output from this demo is proper subset of (and so different from) the larger data analyzed in the referenced RACER manuscript. In particular, we focused on 1000 of the original 10^5 T cells and truncate the thymic selection to be performed on 100 of the original 10^4 self-peptides. 

## System requirement
The code uses MATLab script and is compiled on version R2017b. 

~ ~ ~ ~ ~ ~ ~ ~
* DESCRIPTION for the files:
- RACERMATLab.m: 		Script file which provides an example of thymic selection and T-cell recognition of foreign peptides and point-mutated self-peptides.

- PairwiseAffinity.mat:		MATLab data file containing pairwise binding energy values for 100 thymic self-peptides and 1000 T-cells (peptides delineated by column and T-cells by row)

- PairwiseAffinityMutant.mat:	MATLab data file containing pairwise binding energy values for 1000 point-mutated (non-self) peptides and 1000 T-cells (peptides delineated by column and T-cells by row). The binding energy came from the output of the molecular module of RACER.

- PairwiseAffinityRandom.mat:	MATLab data file containing pairwise binding energy values for 1000 randomly-generated foreign (non-self) peptides and 1000 T-cells (peptides delineated by column and T-cells by row). The binding energy came from the output of the molecular module of RACER.

- RACERMATLab.m:		Function file which generates the T-cell activation energy cutoff on the normalized affinity interval of [0,10] yielding 50% thymic negative selection (variable En50); outputs plots of the distribution of binding energies, thymic selection, and post-selection T-cell recognition profiles of point mutant and random peptides.


* NOTES:
- The T-cells and their corresponding indices are the same across all three input arrays. In other words, binding energies for the j^th T-cell to self-peptides, mutant peptide, and foreign peptide are located at the j^th column of PairwiseAffinity.mat, PairwiseAffinityMutant.mat, and PairwiseAffinityRandom.mat respectively.


* Explanation for the output figures: 
* Figure 1a. Empirical maximum binding energy distributions of T-cells with their self-peptides (maximum for each T-cell taken over all self-peptides)
* Figure 1b. Thymic selection curve (T-cell deletion probability as a function of thymic selection energy cutoff.
* Figure 2a. Post-selection individual T-cell recognition of foreign peptides as a function of T-cell survival probability
* Figure 2b. Post-selection T-cell repertoire recognition of foreign peptides as a function of T-cell survival probability
* Figure 2c. Post-selection individual T-cell recognition of mutant peptides as a function of T-cell survival probability
* Figure 2d. Post-selection T-cell repertoire recognition of mutant peptides as a function of T-cell survival probability



## Reference:
* Rapid Assessment of T-Cell Receptor Specificity of the Immune Repertoire, Xingcheng Lin, Jason T. George, Nicholas P. Schafer, Kevin Ng Chau, Cecilia Clementi, José N. Onuchic, Herbert Levine (https://www.biorxiv.org/content/10.1101/2020.04.06.028415v2)
* RACER borrows the idea of the Principle of Minimal Frustration in protein folding, illustrated here: Learning To Fold Proteins Using Energy Landscape Theory, Nicholas P. Schafer, Bobby L. Kim, Weihua Zheng, Peter G. Wolynes, DOI: 10.1002/ijch.201300145
