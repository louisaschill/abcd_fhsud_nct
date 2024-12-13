% Louisa Schilling - Dec 2024 
% Transition energy calculation 
clear all; close all; 

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy'; 

sc_type = 'avg'; % group average or individual 
parc = 'fs86'; nparc = 86; 
numClusters = 4; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout';

% Load subjectInfo
cd(dataDir);
load('subjectInfo_SUDcohort.mat');
nsubjs = height(subjectInfo);

% Load concTS 
load(['concTS_fs86_' ts_type '_' num2str(nsubjs) 'subj.mat']);
nscans = length(subj_scanInd);

% check subjectInfo in same order as subjkeys 
if ~all(strcmp(subjkeys, subjectInfo.subjectkey))
    error('Subj keys from imaging does not match subjectInfo order');
    
    disp('Putting in same order');
    i = ismember(subjectInfo.subjectkey, subjkeys);
    subjectInfo = subjectInfo(i,:);
    i = ismember(subjkeys,subjectInfo.subjectkey);  subjkeys = subjkeys(i);
    [i,m] = ismember(subjkeys,subjectInfo.subjectkey);
    subjectInfo = subjectInfo(m,:);
    nsubjs = height(subjectInfo);
end

% Load group average SC
if strcmp(sc_type,'avg')
    disp('Using average SC')
    load('/Users/louisaschilling/Desktop/Datasets/ABCD/Data/SC/scavg_fs86_abcd_noself.mat');
    SC = sc_avg_noSelf;

% Load individual SC 
elseif strcmp(sc_type,'indiv')
    disp('Using individual SC (FS68)')
    concTS = concTS(:, 19:86);
    load('/Users/louisaschilling/Desktop/Datasets/ABCD/Data/SC/scmats_sc2080_ifod2act_DK68_orig_volnorm_noself.mat');
    % Exclude subjects without SC and put in same order
    SC = SC(ismember(subjects, subjectInfo.subjectkey));
    sc_subjects = subjects(ismember(subjects, subjectInfo.subjectkey));
    
    [~,m] = ismember(subjects,subjectInfo.subjectkey);
    sc_inds = ismember(subjInd,m); 
    subj_scanInd = subj_scanInd(ismember(subj_scanInd,m));
    concTS = concTS(sc_inds,:);
    scanInd = scanInd(sc_inds); subjInd = subjInd(sc_inds);
    nscans = length(subj_scanInd);nsubjs = length(unique(subj_scanInd));
    [subj_inc_inds,~] = ismember(subjectInfo.subjectkey, subjects);
    subjectInfo = subjectInfo(subj_inc_inds,:);subjkeys = subjkeys(subj_inc_inds);
    [i,m] = ismember(sc_subjects,subjectInfo.subjectkey);
    sc_subjects= sc_subjects(m); SC = SC(m,:);
    % check in same order 
    if ~all(strcmp(sc_subjects, subjectInfo.subjectkey))
        error('SC subjkeys do not match subjectInfo order')
    end

end

% Load optimal T from abcd_tsweep.mat
c = 1;
load([saveDir '/optimalT_k',num2str(numClusters),'_c',num2str(c), '_' ts_type, '_sc_' sc_type '_noself_' num2str(nsubjs) '.mat'],'T'); 

%% Loop through k range
for numClusters = k
    
    % Load parition 
    cd(resultsDir);
    load(['Partition_bp_',parc, '_k',num2str(numClusters),'_', ts_type, num2str(nsubjs) 'subj.mat']);
    
    if strcmp(sc_type,'indiv')
        partition = partition(sc_inds);
        parts = parts(sc_inds,:);
    end 

    load([resultsDir '/subcentroids/subjcentroids_k',num2str(numClusters), '_', ts_type, '_', parc, '_', num2str(nsubjs) '.mat']);

    % if individual SC, get only subcentroids of subjs with SC 
    if strcmp(sc_type, 'indiv')
        sub_centroids = sub_centroids(subj_inc_inds,:,:);
        subjects = subjects(subj_inc_inds);
    end 

    % Calculate TE 
    [E_full,E_region,globalTE,regionalTE] = calculate_subjenergy(SC,sub_centroids,subj_scanInd,c, T,numClusters);
    
    % Save 
    subjkeys = subjectInfo.subjectkey; 
    savename = [ts_type '_' sc_type 'noselfSC_' parc]; 
    disp('Saving subject TEs');
    save(fullfile(saveDir, ['subjenergy_k',num2str(numClusters),'_T' num2str(T),'_c',num2str(c), '_' savename, '_' num2str(nsubjs) '.mat']),...
        'E_full','globalTE','E_region','regionalTE','subjkeys', 'numClusters','clusterNames', '-v7.3');
end
