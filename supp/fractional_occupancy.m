%% Fractional Occupancy 
% Generate subject specific centroids from partition for NCT analysis 
clear all; close all;

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(dataDir);
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures/Reviews'; 

parc = 'fs86'; nparc=86; k = 4; numClusters = k; 
ts_type = 'bp_gsr_gmnorm_exclout';

% Load subject info 
load([dataDir '/subjectInfo_SUDcohort.mat']); 
subjkeys = subjectInfo.subjectkey; 
nsubjs = height(subjectInfo);

% Load concTS 
load(['concTS_fs86_' ts_type '_' num2str(nsubjs) 'subj.mat']);
if strcmp(parc,'fs68'); concTS = concTS(:,19:86); end 

% Load FD
load(fullfile(dataDir, 'MeanFD_PerSubject.mat'), 'resultsTable');
% Check subjectInfo and TE results in same order 
if ~all(strcmp(resultsTable.subjectkey,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD mean');
end 
subjectInfo.FD_mean = resultsTable.FD_mean; 
subjectInfo.FD_mean_nanoutlier = resultsTable.FD_mean_nanoutlier; 

nscans = length(subj_scanInd);
savedir = fullfile('/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/subcentroids'); 
cd(savedir); 

% Load parition
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results'; cd(resultsDir);
load(fullfile(resultsDir,['Partition_bp_' parc, ...
    '_k',num2str(numClusters),'_', ts_type,  num2str(nsubjs), 'subj.mat']));

%% Calculate fractional occupancy
% FO by scan 
nscans = length(unique(scanInd)); 
fo1 = zeros(numClusters, nscans);
for i = 1:nscans
    scan_partition = partition(scanInd == i);
    scan_tr = length(scan_partition); 
    scan_time = scan_tr * 0.8;
    scan_count = hist(scan_partition,1:numClusters);
    scan_fo1 = scan_count/scan_tr; 
    fo1(:,i) = scan_fo1; 
end 

% Average FO per subject
subjects = unique(subj_scanInd,'stable'); 
nsubjs = length(subjects); 
fo = NaN(numClusters,nsubjs);
for i = 1:nsubjs
    subject = subjects(i);
    fo(:,i) = mean(fo1(:,subj_scanInd == subject),2);
end

fo = fo'; 

%% ANCOVA
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'); 
anovaVarNames = {'Sex','Age','FHSUD','FD','Model','Income',...
    'Parent Education','Race','Prenatal','Parent MH', 'Puberty'};
interaction_term_sets = {{'Sex', 'FHSUD'},{'FHSUD','Income'},{'Sex','Puberty'}};
contVar = [2,4];
anovaVars = {subjectInfo.sex(inds)=='F', subjectInfo.age(inds),...
    subjectInfo.FHSUD(inds) == 'FH+', subjectInfo.FD_mean(inds),...
    subjectInfo.model(inds), subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds), categorical(subjectInfo.race(inds)),...
    categorical(subjectInfo.subDuringPregAfter(inds)),...
    subjectInfo.parentMH(inds) == 1,...
    categorical(subjectInfo.pds_mod(inds))};

for i= 1:4
    [all_results,stats] = run_ANOVA(fo(inds,i),anovaVars,contVar,...
        anovaVarNames,interaction_term_sets,'off');
    [F_stats_sexsud,p_val_sexsud, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
    [F_stats_sud, p_val_sud, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');
    disp(all_results)
end

%% Pie charts 

f=figure; subplot(2,2,1); pie(mean(fo(subjectInfo.sexFHSUD=='F-FH+',:),1));
title('Female FH+');
legend(clusterNames,"Location","bestoutside");

subplot(2,2,2);pie(mean(fo(subjectInfo.sexFHSUD=='F-FH-',:),1));
title('Female FH-');

subplot(2,2,3);pie(mean(fo(subjectInfo.sexFHSUD=='M-FH+',:),1));
title('Male FH+');

subplot(2,2,4);pie(mean(fo(subjectInfo.sexFHSUD=='M-FH-',:),1));
title('Male FH-');

sgtitle('Fractional Occupancy')
saveas(f,fullfile(figureDir,'fractionalocc_pie.png'));
