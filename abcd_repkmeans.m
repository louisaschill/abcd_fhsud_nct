%% Louisa Schilling - updated Dec 2024 
% Run k-means clustering across full range of possible k's to find optimal k 
% Dependencies: run_kmeans, GET_CENTROIDS, NAME_CLUSTERS_ANGLE

clear all; close all;

% Parameters
parc = 'fs86'; 
ts_type = 'bp_gsr_gmnorm_exclout';
distanceMethod = 'correlation'; % Distance metric for k-means
nreps = 5;                     % # of repetitions 
maxIter = 1000;                 % Maximum iterations

% Paths
basedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data/';
cd(basedir);

% Load Subject Information
load('subjectInfo_SUDcohort.mat');
nsubjs = height(subjectInfo);

% Load Concatenated Time Series
load(['concTS_' parc '_' ts_type '_' num2str(nsubjs) 'subj.mat'], 'concTS', 'subjInd', 'scanInd', 'subj_scanInd', 'subjkeys');

% Find minimum number of frames from a single subject 
framesPerSubject = accumarray(subjInd, 1); % Count occurrences of each subjInd
minFrames = min(framesPerSubject);
kMax = floor(sqrt(minFrames)); % k*k must be less than min number of frames across subjects

kRange = 2:kMax;% range of k

% K-means Clustering
for numClusters = kRange
    disp(['Running K-means for k = ' num2str(numClusters)]);
    
    % Perform K-means clustering
    [partition, D, overallClusters, subjectClusters] = run_kmeans(...
        concTS, numClusters, distanceMethod, nreps, maxIter, subjInd, []);
    
    % Get cluster centroids and assign names based on Yeo networks
    centroids = GET_CENTROIDS(concTS, partition, numClusters);
    clusterNames = NAME_CLUSTERS_ANGLE(centroids); % Align clusters with prior Yeo labels
    
    % Save Results
    savename = sprintf('%s_%s_%dsubj', parc, ts_type, nsubjs);
    saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/';
    save(fullfile(saveDir, sprintf('kmeans_k_%d_%s.mat', numClusters, savename)), ...
        'partition', 'D', 'subjectClusters', 'overallClusters', 'subjkeys', ...
        'clusterNames', 'numClusters');
end

disp('K-means clustering completed.');
