RACERMATLabOutput
~ ~ ~ ~ ~ ~ ~ ~
DESCRIPTION:
- RACERMATLab.m: 		Script file which provides an example of thymic selection and T-cell recognition of foreign peptides and point-mutated self-peptides.

- PairwiseAffinity.mat:		MATLab data file containing pairwise binding energy values for 100 thymic self-peptides and 1000 T-cells (peptides delineated by column and T-cells by row)

- PairwiseAffinityMutant.mat:	MATLab data file containing pairwise binding energy values for 1000 point-mutated (non-self) peptides and 1000 T-cells (peptides delineated by column and T-cells by row)

- PairwiseAffinityRandom.mat:	MATLab data file containing pairwise binding energy values for 1000 randomly-generated foreign (non-self) peptides and 1000 T-cells (peptides delineated by column and T-cells by row)

- RACERMATLab.m:		Function file which generates the T-cell activation energy cutoff on the normalized affinity interval of [0,10] yielding 50% thymic negative selection (variable En50); outputs plots of the distribution of binding energies, thymic selection, and post-selection T-cell recognition profiles of point mutant and random peptides.


NOTES:
- The T-cells and their corresponding indices are the same across all three input arrays. In other words, binding energies for the j^th T-cell to self-peptides, mutant peptide, and foreign peptide are located at the j^th column of PairwiseAffinity.mat, PairwiseAffinityMutant.mat, and PairwiseAffinityRandom.mat respectively.

- MATLab source script and functions compiled on version R2017b.