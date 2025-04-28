%% CBCL correlation analysis 
clear all; close all; clc;

%% Set Directories
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures/Final';

%% Load Data
load(fullfile(dataDir, 'subjectInfo_SUDcohort.mat'));
nsubjs = height(subjectInfo);

parc = 'fs86';
numClusters = 4; k = 4; ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg';
c = 1;

% Load optimal T
load(fullfile(resultsDir, ['optimalT_k', num2str(numClusters), '_c', num2str(c), '_', ts_type, '_sc_', sc_type, '_noself_', num2str(nsubjs), '.mat']), 'T');

% Load Energy Data
load(fullfile(resultsDir, ['subjenergy_k', num2str(numClusters), '_T', num2str(T), '_c', num2str(c), '_', ts_type, '_', sc_type, 'noselfSC_', parc, '_', num2str(nsubjs), '.mat']));

% Verify subject key match
assert(all(strcmp(subjkeys, subjectInfo.subjectkey)), 'Subject keys mismatch.');

%% Calculate Network TE
networks = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
YeoLUT = table2array(readtable('fs86_to_yeo.csv'));
numNetworks = numel(networks);

E_network = nan(nsubjs, length(clusterNames)^2, numNetworks);
for i = 1:numNetworks
    E_network(:,:,i) = sum(E_region(:,:,YeoLUT == i), 3);
end
networkTE = squeeze(mean(E_network, 2));

%% Merge FD data
load(fullfile(dataDir, 'MeanFD_PerSubject.mat'), 'resultsTable');
assert(all(strcmp(resultsTable.subjectkey, subjectInfo.subjectkey)), 'FD subject keys mismatch.');
subjectInfo.FD_mean = resultsTable.FD_mean;
subjectInfo.FD_mean_nanoutlier = resultsTable.FD_mean_nanoutlier;

%% Remove outliers based on globalTE
outlier = isoutlier(globalTE, 'quartiles');
subjectInfo = subjectInfo(~outlier, :);
globalTE = globalTE(~outlier);
networkTE = networkTE(~outlier, :);
TE_data = [globalTE, networkTE]; % Predictor matrix

%% Behavioral Variables
varNames = {'Externalizing', 'Internalizing', 'WithDep', 'AnxDep', 'Attention', 'Somatic', 'Social', 'Thought', ...
    'RuleBreak', 'Aggression', 'Total Problems', '-[BIS]', 'BAS Reward','BAS Fun', 'BAS Drive', ...
    'UPPS Negative Urgency', 'UPPS Lack Planning', 'UPPS Sensation Seeking', 'UPPS Positive Urgency', 'UPPS Lack Perseverance', 'UPPS Sum'};

D = table2array(table(...
    subjectInfo.external, subjectInfo.internal, subjectInfo.withdep, subjectInfo.anxdep, subjectInfo.attention, ...
    subjectInfo.somatic, subjectInfo.social, subjectInfo.thought, subjectInfo.rulebreak, subjectInfo.aggressive, ...
    subjectInfo.totprob, -subjectInfo.bis_sum, subjectInfo.bas_reward_mod, subjectInfo.bas_fun, subjectInfo.bas_drive_mod, ...
    subjectInfo.upps_neg_urgency, subjectInfo.upps_lack_plan, subjectInfo.upps_sensation_seek, subjectInfo.upps_pos_urgency, ...
    subjectInfo.upps_lack_preservance, subjectInfo.upps_sum, ...
    'VariableNames', varNames));

%% Covariates
income_dummy = dummyvar(categorical(subjectInfo.income_cat)); income_dummy(:,1) = [];
parentEd_dummy = dummyvar(categorical(subjectInfo.parentEd_cat)); parentEd_dummy(:,1) = [];
race_dummy = dummyvar(categorical(subjectInfo.race)); race_dummy(:,1) = [];
model_dummy = dummyvar(categorical(subjectInfo.model)); model_dummy(:,1) = [];
puberty_dummy = dummyvar(categorical(subjectInfo.pds_mod)); puberty_dummy(:,1) = [];

covariates = [subjectInfo.age, subjectInfo.FD_mean_nanoutlier, ...
    model_dummy, income_dummy, parentEd_dummy, race_dummy, ...
    subjectInfo.subDuringPreg == 1, subjectInfo.parentMH == 1, puberty_dummy];

%% Settings
corrType = 'spearman'; 
figSize = [100 100 1400 1000]; 
map = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1); ones(128,1), linspace(1,0,128)', linspace(1,0,128)']; % blue-white-red

group = 'All';
partial = false; % Partial correlation toggle

if strcmp(group, 'All')
    indsm = subjectInfo.sex == 'M';
    indsf = subjectInfo.sex == 'F';
elseif strcmp(group, 'FH+')
    indsm = subjectInfo.sexFHSUD == 'M-FH+';
    indsf = subjectInfo.sexFHSUD == 'F-FH+';
end

networkLabels = [{'Global'}, networks];

%% Males

f = figure('Position', figSize);

if partial
    [corrMatrix, pValues] = partialcorr(TE_data(indsm,:), D(indsm,:), covariates(indsm,:), 'type', corrType);
else
    [corrMatrix, pValues] = corr(TE_data(indsm,:), D(indsm,:), 'type', corrType);
end

% FDR correction
pcorr = mafdr(pValues(:), 'BHFDR', true);
pcorr = reshape(pcorr, size(pValues));

% Plot
imagesc(corrMatrix);
colormap(map);
colorbar; ylabel(colorbar,'R');
lim = max(abs(corrMatrix(:))) + 0.01; 
caxis([-lim lim]);
set(gca, 'XTick', 1:numel(varNames), 'XTickLabel', varNames, 'XTickLabelRotation', 45, 'FontSize', 16);
set(gca, 'YTick', 1:numel(networkLabels), 'YTickLabel', networkLabels, 'FontSize', 16);
xlabel('Behavioral/Psychological Measure', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Global and Network TE', 'FontSize', 16, 'FontWeight', 'bold');
title([group ' Males: Correlation Heatmap']);

% Annotate
for i = 1:size(corrMatrix,1)
    for j = 1:size(corrMatrix,2)
        valStr = sprintf('%.2f', corrMatrix(i,j));
        if pValues(i,j) < 0.05
            if pcorr(i,j) < 0.05
                text(j, i, [valStr '**'], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            else
                text(j, i, [valStr '*'], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        else
            text(j, i, valStr, 'FontSize', 14, 'HorizontalAlignment', 'center');
        end
    end
end

saveas(f, fullfile(figureDir, 'all_males_corr_cbcl.png'));

%% Females

f = figure('Position', figSize);

if partial
    [corrMatrix, pValues] = partialcorr(TE_data(indsf,:), D(indsf,:), covariates(indsf,:), 'type', corrType);
else
    [corrMatrix, pValues] = corr(TE_data(indsf,:), D(indsf,:), 'type', corrType);
end

pcorr = mafdr(pValues(:), 'BHFDR', true);
pcorr = reshape(pcorr, size(pValues));

imagesc(corrMatrix);
colormap(map);
colorbar; ylabel(colorbar,'R');
lim = max(abs(corrMatrix(:))) + 0.01; 
caxis([-lim lim]);      
set(gca, 'XTick', 1:numel(varNames), 'XTickLabel', varNames, 'XTickLabelRotation', 45, 'FontSize', 16);
set(gca, 'YTick', 1:numel(networkLabels), 'YTickLabel', networkLabels, 'FontSize', 16);
xlabel('Behavioral/Psychological Measure', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Global and Network TE', 'FontSize', 16, 'FontWeight', 'bold');
title([group ' Females: Correlation Heatmap']);

% Annotate
for i = 1:size(corrMatrix,1)
    for j = 1:size(corrMatrix,2)
        valStr = sprintf('%.2f', corrMatrix(i,j));
        if pValues(i,j) < 0.05
            if pcorr(i,j) < 0.05
                text(j, i, [valStr '**'], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            else
                text(j, i, [valStr '*'], 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
            end
        else
            text(j, i, valStr, 'FontSize', 14, 'HorizontalAlignment', 'center');
        end
    end
end

saveas(f, fullfile(figureDir, 'all_females_corr_cbcl.png'));
