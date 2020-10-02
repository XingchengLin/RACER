%RACERMatlabScript

%% 1. Load necessary data
load('PairwiseAffinity.mat');
load('PairwiseAffinityMutant.mat');
load('PairwiseAffinityRandom.mat');

% Data files contains pairwise energies of (column, row):
%   PairwiseAffinity:      (self-peptide,T-cell)    10^2-by-10^3
%   PairwiseAffinityMutant:(mutant-peptide,T-cell)  10^3-by-10^3
%   PairwiseAffinityRandom:(random-peptide,T-cell)  10^3-by-10^3

%Run RACERMATLab function
RACERMATLab(PairwiseAffinity,PairwiseAffinityMutant,PairwiseAffinityRandom)

