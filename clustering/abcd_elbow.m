%% Louisa Schilling - updated Dec 2024
% Makes elbow plot of total variance explained by clusters at different k values
% Run after abcd_repkmeans to get clustering across range of k 
% Dependencies: GET_CENTROIDS, VAREXPLAINED (from EJC: https://github.com/ejcorn/brain_states) 

clear all; close all;
parc = 'fs86';  ts_type = 'bp_gsr_gmnorm_exclout';
 
basedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(basedir);
savedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/repeatkmeans'; 
figDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures';

% Load subject info 
load('subjectInfo_SUDcohort.mat'); 
subjkeys = subjectInfo.subjectkey; 
nsubjs = height(subjectInfo);

% Load concTS 
load(['concTS_fs86_' ts_type '_' num2str(nsubjs) 'subj.mat']);
if strcmp(parc, 'fs68') % select only cortical regions 
    concTS = concTS(:,19:86);
end

% Find minimum number of frames from a single subject 
framesPerSubject = accumarray(subjInd, 1); % Count occurrences of each subjInd
minFrames = min(framesPerSubject);
maxk = floor(sqrt(minFrames)); % k*k must be less than min number of frames across subjects

k_rng = 2:maxk;% range of k

%% Calculate variance explained across range of k 
VarianceExplained = zeros(length(k_rng),1);
count = 0; 
for numClusters = k_rng
	count = count + 1;
    disp(['K = ',num2str(numClusters)])
    load(fullfile(dataDir,['kmeans_k_',num2str(numClusters),'_', parc, '_' ts_type '_' num2str(nsubjs) 'subj.mat']));
	kClusterCentroids = GET_CENTROIDS(concTS,partition,numClusters);
	[VarianceExplained(count), x, y] = VAREXPLAINED(concTS,partition,kClusterCentroids,numClusters);
end

mkdir(savedir); cd(savedir); 

save(fullfile(savedir,['VarianceExplained_' parc '_' ts_type '.mat']),'VarianceExplained','k_rng');

%% Determine k where R^2 gain drops below 1%
R2_gain = diff(VarianceExplained);
k_drop_threshold = find(R2_gain < 0.01, 1, 'first') + 2;
disp(['K where R^2 gain drops below 1%: ', num2str(k_drop_threshold)]);

%% Plot 
f=figure;
plot(k_rng,VarianceExplained,'.-r');
xlabel('\it{k}'); ylabel('R^2');
title('Variance Explained by Clustering');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(figDir,['ElbowPlot_',ts_type,'.pdf']));

f = figure;
plot(k_rng(2:end),diff(VarianceExplained),'.-b')
xlabel('\it{k}'); ylabel('R^2 Gain');
title('R^2 Gain by increasing \it{k}');
prettifyEJC;
f.PaperUnits = 'centimeters';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
saveas(f,fullfile(figDir,['GainInVarianceExplained_' ts_type '.pdf']));
