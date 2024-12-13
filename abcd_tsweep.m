%% Louisa Schilling - updated Dec 2024
% find optimal T for largest negative correlation between transition energy
% and transition probability
clear all; close all;

% set paths 
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(dataDir);
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy'; 

% set inputs
numClusters = 4; k = 4; 

sc_type = 'avg';
parc = 'fs86'; nparc = 86;
ts_type = 'bp_gsr_gmnorm_exclout';

% Load subjectInfo
load([dataDir '/subjectInfo_SUDcohort.mat']);
nsubjs = height(subjectInfo);

% Load concTS
load([dataDir '/concTS_' parc '_' ts_type '_' num2str(nsubjs) 'subj.mat']);

% Load avg SC
load('/Users/louisaschilling/Desktop/Datasets/ABCD/Data/SC/scavg_fs86_abcd_noself.mat');
SC = sc_avg_noSelf;

% Load parititon
load([resultsDir '/Partition_bp_',parc, '_k',num2str(numClusters),'_', ts_type, num2str(nsubjs) 'subj.mat']);

% Load subjcentroids
load([resultsDir '/subcentroids/subjcentroids_k',num2str(numClusters), '_',ts_type,'_' parc '_' num2str(nsubjs) '.mat']);

%% Calculate transition probability 
%load(fullfile(resultsDir,'TransProbs',['TransProbsData_', parc, '_k',num2str(numClusters),'_',ts_type,'.mat']));
[transitionProbability2D,transitionProbabilityMats] = GET_TRANS_PROBS(partition,subjInd);

%% Calculate TE for sweep of T 
c = 1; 

% for range of T, determine which gives greatest neg corr for TE x TP 
T_rng = [0.001:0.5:10]; nT = length(T_rng);

Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

x0 = centroids(:,Xo_ind);
xf = centroids(:,Xf_ind);

% now each column of x0 and xf represent state transitions
E_full_grpavg_T = NaN(nT,numClusters^2);
E_w_grpavg_T = NaN(nT,numClusters^2);

Anorm = NORMALIZE(SC,c); % normalize A by maximum eigenvalue - eye(N) to make marginally stable

for i=1:nT
    T=T_rng(i);
    WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
    E_full_grpavg_T(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each state transition
end

%% corrs Persist
for i=1:nT
    [R(i),P(i)] = corr(E_full_grpavg_T(i,:)',mean(transitionProbability2D)','type','Spearman');
end

%% plot R vs T (no abs)
figure; 
plot(T_rng,R); title('Transition Energy vs Transition Probability');

[m,i]=min(R); 
disp(['Largest negative correlation: T = ' num2str(T_rng(i)) ' r = ' num2str(R(i))]); 

%% save optimal T 
T = T_rng(i); 
save([saveDir '/optimalT_k',num2str(numClusters),'_c',num2str(c), '_' ts_type, '_sc_' sc_type '_noself_' num2str(nsubjs) '.mat'],'T'); 
