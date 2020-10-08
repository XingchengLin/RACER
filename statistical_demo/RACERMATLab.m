function En50=RACERMATLab(SelfAffinity,MutantAffinity,RandomAffinity)

%DESCRIPTION: This code inputs pairwise affinity values for self-peptides,
%point-mutated self peptides, and random peptides, all to the same
%collection of T-cells, and then outputs plots of the distribution of
%binding energies, thymic selection, as well as post-selection T-cell 
%recognition profiles of point mutant and random peptides.

%INPUT:
%   SelfAffinity  : Pairwise energy array between self-peptides 
%                   (columns) and T-cells (rows)
%   MutantAffinity: Pairwise energy array between mutant peptides 
%                   (columns) and T-cells (rows)
%   RandomAffinity: Pairwise energy array between random peptides 
%                   (columns) and T-cells (rows)

%OUTPUT:
%   En50=

%   1a. Empirical maximum binding energy distributions of T-cells with their
%      self-peptides (maximum for each T-cell taken over all self-peptides)
%   1b. Thymic selection curve (T-cell deletion probability as a function of
%      thymic selection energy cutoff.

%   2a. Post-selection individual T-cell recognition of foreign peptides as
%      a function of T-cell survival probability
%   2b. Post-selection T-cell repertoire recognition of foreign peptides as
%      a function of T-cell survival probability
%   2c. Post-selection individual T-cell recognition of mutant peptides as
%      a function of T-cell survival probability
%   2d. Post-selection T-cell repertoire recognition of mutant peptides as
%      a function of T-cell survival probability
tic;
AffinityMD=SelfAffinity;
AffinityMDMutant=MutantAffinity;
AffinityMDRandom=RandomAffinity;


%% 2. Thymic selection of donor TCRs.
disp('Performing thymic selection');
[Nn Nt]=size(AffinityMD);
[NnMutant Nt]=size(AffinityMDMutant);
[NnRandom Nt]=size(AffinityMDRandom);

AffinityMDTemp=AffinityMD;


% Normalization: Energy --> affinity
%Affinity=-Affinity-globalminima
AffinityMD=AffinityMDTemp;

AffinityMD=-AffinityMD;
AffinityMDMutant=-AffinityMDMutant;
AffinityMDRandom=-AffinityMDRandom;

globalMinima=min(...
[min(min(AffinityMD)) min(min(AffinityMDRandom))...
min(min(AffinityMDMutant))]);

AffinityMD=AffinityMD-globalMinima;
AffinityMDRandom=AffinityMDRandom-globalMinima;
AffinityMDMutant=AffinityMDMutant-globalMinima;


%Normalized to compact domain [0,10]
globalMaxima=max(...
    [max(max(AffinityMD)) max(max(AffinityMDRandom))...
    max(max(AffinityMDMutant))]);

AffinityMD=AffinityMD/globalMaxima*10;
AffinityMDRandom=AffinityMDRandom/globalMaxima*10;
AffinityMDMutant=AffinityMDMutant/globalMaxima*10;

Delta=max(max(AffinityMD))-min(min(AffinityMD));
DeltaRandom=max(max(AffinityMDRandom))-min(min(AffinityMDRandom));
DeltaMutant=max(max(AffinityMDMutant))-min(min(AffinityMDMutant));

Partition=1000;

xMD=[0 : ...
    min([Delta DeltaRandom DeltaMutant])/Partition : ...
    max([max(max(AffinityMD)) max(max(AffinityMDRandom))...
    max(max(AffinityMDMutant))])];

% MD Selection Curve
RecognitionFraction=0.5; %Defines En cutoff
z=1;
for i=1:Nt %For each TCR
    MaxAffinityMD(i)=max(AffinityMD(:,i)); %<-Maximum affinity for TCR i
end


for i=1:length(xMD)
    RecognitionMD(i,z)=sum(MaxAffinityMD>xMD(i));
end
RecognitionMD(:,z)=RecognitionMD(:,z)/Nt;
RecognitionMean=mean(RecognitionMD');

IndexSelection=find(abs(RecognitionMD-RecognitionFraction)...
    ==min(abs(RecognitionMD-RecognitionFraction)));
En=xMD(IndexSelection); %Energy cutoff for 50% negative selection survival
En50=En;
disp('Thymic selection complete'); toc;
%Check robustness of deriving selection curve from less TCRs


%% Random Peptide Recognition
disp('Performing random peptide recognition');

CountMDRandom=zeros(length(xMD),1);
PercentMDRandom=zeros(length(xMD),1);
CountMDRepertoireRandom=zeros(length(xMD),1);
PercentMDRepertoireRandom=zeros(length(xMD),1);

toc; 

for w=1:length(xMD) %For each E_n threshold
    En=xMD(w); %En scan
    TCRMD=MaxAffinityMD<=En; %1 if TCRs survive negative selection
    TCRMDRandom=AffinityMDRandom<=En;  %1 if TCRs recognie random peptide
    
    for y=1:NnRandom %For each Random peptide
        CountMDRandom(w) = CountMDRandom(w) +...
            sum( TCRMD & ~TCRMDRandom(y,:) ); %Survives selection and 
%                                              recognizes random peptide
        if sum(TCRMD & ~TCRMDRandom(y,:) )>=1
            CountMDRepertoireRandom(w) = CountMDRepertoireRandom(w) + 1;
        end
    end
    if sum(TCRMD==1)==0 %No TCRs survive negative selection
        PercentMDRandom(w)=NaN;
        PercentMDRepertoireRandom(w)=NaN;
    else
        PercentMDRandom(w)=CountMDRandom(w)/(NnRandom*sum(TCRMD));
        PercentMDRepertoireRandom(w)=CountMDRepertoireRandom(w)/NnRandom;
    end                   
end



disp('Random peptide recognition complete'); toc;

%% MD Mutant Peptide Recognition
disp('Performing mutant peptide recognition');

CountMDMutant=zeros(length(xMD),1);
PercentMDMutant=zeros(length(xMD),1);
CountMDRepertoireMutant=zeros(length(xMD),1);
PercentMDRepertoireMutant=zeros(length(xMD),1);

for w=1:length(xMD) %For each E_n threshold
    En=xMD(w); %En scan
    TCRMD=MaxAffinityMD<=En; %1 if TCRs survive negative selection
    TCRMDMutant=AffinityMDMutant<=En;  %1 if TCRs recognize mutant peptide
    
    for y=1:NnMutant %For each mutant peptide
        CountMDMutant(w) = CountMDMutant(w) + ...
            sum(TCRMD & ~TCRMDMutant(y,:) ); %Single TCR Recognition
        if sum(TCRMD & ~TCRMDMutant(y,:) )>=1
            CountMDRepertoireMutant(w) = CountMDRepertoireMutant(w) + 1;
        end
    end
    if sum(TCRMD==1)==0 %No TCRs survive negative selection
        PercentMDMutant(w)=NaN;
        PercentMDRepertoireMutant(w)=NaN;
    else
        PercentMDMutant(w)=CountMDMutant(w)/(NnMutant*sum(TCRMD));
        PercentMDRepertoireMutant(w)=CountMDRepertoireMutant(w)/NnMutant;
    end   
end



disp('Mutant peptide recognition complete');

%% Figure 1:
figure; 
subplot(1,2,1); hold on; box on;
histogram(MaxAffinityMD,'BinLimits',...
    [min(min(AffinityMD)) max(max(AffinityMD))],...
    'normalization','probability');
xlabel('Maximum T-cell affinity'); ylabel('Probability');
title('Empirical distribution: max binding affinity over self-peptides');

subplot(1,2,2); hold on; box on;
plot(xMD,RecognitionMD,'--','LineWidth',2,'Color',[0.494 0.184 0.557]);
title('Thymic selection curve');
xlabel('Negative selection threshold');
ylabel('T-cell recognition probbaility')
ylim([0,1]);
set(gcf, 'Position',  [200, 200, 900, 500]);

savefig('MaxEnergiesSelectionCurve.fig');
saveas(gcf,'MaxEnergiesSelectionCurve.png');

%% Figure 2
%Individual TCR Recognition vs PSurvival (1-RecognitionMD)
figure; subplot(2,2,1); box on;
semilogy(1-RecognitionMD,PercentMDRandom,...
    '-*','Color',[0.494 0.184 0.557]);
hold on;
xlabel('Survival probability'); ylabel('Recognition probability');
title(sprintf(...
    'Single T-cell recognition rate of %.2g foreign peptides',...
    NnRandom));

%Repertoire TCR Recognition vs PSurvival (1-RecognitionMD)
subplot(2,2,2); box on;
semilogy(1-RecognitionMD,PercentMDRepertoireRandom,...
    '-*','Color',[0.494 0.184 0.557]);
hold on;
xlabel('Survival Probability'); ylabel('Recognition Probability');
title(sprintf(...
    'T-cell repertoire recognition rate of %.2g foreign peptides',...
    NnRandom));
ylim([0 1]);

%Individual TCR Recognition vs PSurvival (1-RecognitionMD)
subplot(2,2,3); box on;
semilogy(1-RecognitionMD,PercentMDMutant,...
    '*','Color',[0.6 0.4 1.0]);
hold on;
xlabel('Survival Probability'); ylabel('Recognition Probability');
title(sprintf(...
    'Single T-cell recognition rate of %.2g mutant peptides',...
    NnRandom));

%Repertoire TCR Recognition vs PSurvival (1-RecognitionMD)
subplot(2,2,4); box on;
semilogy(1-RecognitionMD,PercentMDRepertoireMutant,...
    '*','Color',[0.6 0.4 1.0]);
hold on;
xlabel('Survival Probability'); ylabel('Recognition Probability');
title(sprintf(...
    'T-cell repertoire recognition rate of %.2g mutant peptides',...
    NnRandom));
ylim([0 1]);

set(gcf, 'Position',  [200, 200, 900, 500]);


savefig('TCellRecognitionRates.fig');
saveas(gcf,'TCellRecognitionRates.png');


end