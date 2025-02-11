%% Louisa Schilling - updated Dec 2024
% Transition energy calculation
clear; close all;

% Directories and parameters
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';

sc_type = 'avg'; % SC type: 'avg' (group average) or 'indiv' (individual)
parc = 'fs86'; nparc = 86; % Parcellation type and number of parcels
numClusters = 4; k = 4; % Number of clusters/k
ts_type = 'bp_gsr_gmnorm_exclout'; % rsfmri timeseries type

% Load subject info
cd(dataDir);
load('subjectInfo_SUDcohort.mat');
nsubjs = height(subjectInfo);

% Load concatenated time series (concTS)
load(['concTS_fs86_', ts_type, '_', num2str(nsubjs), 'subj.mat']);
nscans = length(subj_scanInd);

% Ensure subjectInfo matches subjkeys
if ~all(strcmp(subjkeys, subjectInfo.subjectkey))
    error('Subj keys from imaging do not match subjectInfo order');
else
    [sortedKeys, sortIdx] = sort(subjkeys);
    subjectInfo = subjectInfo(sortIdx, :);
end

c = 1; % Constant parameter for energy calculation

% Load same optimal T caluclated using group avg SC
load([saveDir, '/optimalT_k', num2str(numClusters), '_c', num2str(c),...
    '_', ts_type, '_sc_avg_noself_', num2str(nsubjs), '.mat'], 'T');

%% Load SC matrix
if strcmp(sc_type, 'avg')
    disp('Using average SC');
    load('/Users/louisaschilling/Desktop/Datasets/ABCD/Data/SC/scavg_fs86_abcd_noself.mat');
    SC = sc_avg_noSelf;
elseif strcmp(sc_type, 'indiv')
    disp('Using individual SC');
    concTS = concTS(:,19:86);
    load(['/Users/louisaschilling/Desktop/Datasets/ABCD/Data/SC/scmats_' ...
        'sc2080_ifod2act_DK68_orig_volnorm_noself.mat']);

    i = ismember(subjects, subjectInfo.subjectkey);
    SC = SC(i); sc_subjects = subjects(i); clear subjects;

    [~,m] = ismember(sc_subjects,subjectInfo.subjectkey);
    sc_inds = ismember(subjInd,m);
    
    subjInd = subjInd(sc_inds); concTS = concTS(sc_inds,:); scanInd = scanInd(sc_inds);
    subj_scanInd = subj_scanInd(ismember(subj_scanInd,m));

    [subj_inc_inds,i] = ismember(subjectInfo.subjectkey, sc_subjects);
    subjectInfo = subjectInfo(subj_inc_inds,:);

    if ~all(strcmp(sc_subjects, subjectInfo.subjectkey))
        error('SC subjkeys do not match subjectInfo order')
    end
end

%% For k range, load subj centroids, calculate TE and save
for numClusters = k

    % Load partition 
    cd(resultsDir); 
    load(['Partition_bp_fs86', '_k', num2str(numClusters), '_', ts_type, ...
        num2str(nsubjs), 'subj.mat']);
    
    if strcmp(sc_type, 'indiv')
        partition = partition(sc_inds);
    end

    % Load subject centroids
    load([resultsDir, '/subcentroids/subjcentroids_k', num2str(numClusters), ...
        '_', ts_type, '_', parc, '_', num2str(nsubjs), '.mat']);

    if strcmp(sc_type, 'indiv')
        sub_centroids = sub_centroids(subj_inc_inds,:,:);
        subjects = subjects(subj_inc_inds);
    end

    % % Calculate subject-specific centroids from partition
    % [sub_centroids,subjects,clusterNames_sub,net9angle_Up, net9angle_Down] = abcd_generate_subcentroids(concTS,partition, ...
    %     ts_type,numClusters,subjInd,scanInd,subj_scanInd,length(subj_scanInd),[]);

    % Calculate transition energy
    [E_full, E_region, globalTE, regionalTE] = calculate_subjenergy(SC, sub_centroids, subj_scanInd,c,T,numClusters);

    nsubjs = length(unique(subj_scanInd));

    % Save results
    subjkeys = subjectInfo.subjectkey;
    savename = [ts_type, '_', sc_type, 'noselfSC_', parc];
    disp('Saving subject TEs');
    nsubjs = height(subjectInfo);
    save(fullfile(saveDir, ['subjenergy_k', num2str(numClusters), '_T', num2str(T), '_c', num2str(c), '_', savename, '_', num2str(nsubjs), '.mat']),...
        'E_full', 'globalTE', 'E_region', 'regionalTE', 'subjkeys', 'clusterNames', 'numClusters', '-v7.3');
end
