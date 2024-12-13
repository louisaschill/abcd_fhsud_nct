% Louisa Schilling - updated December 2024 
% Concatenate timeseries of all subjects for clustering 

clear all; close all; 
parc = 'fs86'; ts_type = 'bp_gsr';

gm_norm = true; exclout = true; 

% Load subject information and imaging data
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(dataDir);
load('subjectInfo_SUDcohort.mat'); 
subjkeys = subjectInfo.subjectkey; 
nsubjs = height(subjectInfo);

load([dataDir '/rsfMRI_' parc '_' ts_type '_' num2str(nsubjs) 'subj.mat']);
numScans = height(imagingInfo);

if gm_norm, ts_type = [ts_type '_gmnorm']; end
if exclout, ts_type = [ts_type '_exclout']; end

% Preallocate
maxEntries = numScans * 375; 
concTS = nan(maxEntries, 86);
scanInd = nan(maxEntries, 1);
subjInd = nan(maxEntries, 1);

length_count = 0;
scanCount = 0;
subj_scanInd = imagingInfo.subjInd; 

% Grab data from each subj and add to concTS 
subjects = unique(imagingInfo.subjectkey, 'stable');
for i = 1:nsubjs
    subj = subjects{i};
    scans = imagingInfo(strcmp(imagingInfo.subjectkey, subj), :);
    
    scansTS = [];
    for s = 1:height(scans)
        scanCount = scanCount + 1;

        % Select the appropriate time series
        if gm_norm && exclout
            ts = scans.ts_gmnorm_nanoutlier{s};
        elseif exclout && ~gm_norm
            ts = scans.ts_nanoutlier{s};
        elseif ~exclout && ~gm_norm
            ts = scans.ts{s};
        elseif ~exclout && gm_norm
            ts = scans.ts_gmnorm{s};
        end

        % Remove rows with NaNs
        ts = ts(~all(isnan(ts), 2), :);

        scansTS = [scansTS; ts]; % Concatenate scan time series
        rowsToUpdate = length_count + (1:size(ts, 1));
        scanInd(rowsToUpdate) = scanCount;
        subjInd(rowsToUpdate) = i;
        length_count = length_count + size(ts, 1);
    end

    rowsToUpdate = length_count - size(scansTS, 1) + 1:length_count;
    concTS(rowsToUpdate, :) = scansTS;
end

% Trim 
concTS = concTS(1:length_count, :);
scanInd = scanInd(1:length_count);
subjInd = subjInd(1:length_count);

% Remove rows with all NaNs
validRows = ~all(isnan(concTS), 2);
concTS = concTS(validRows, :);
scanInd = scanInd(validRows);
subjInd = subjInd(validRows);

% Save concatenated time series
savename = sprintf('concTS_%s_%s_%dsubj.mat', parc, ts_type, nsubjs);
save(fullfile(dataDir, savename), 'concTS', 'subjInd', 'scanInd','subj_scanInd', 'subjkeys', '-v7.3');
disp(['Saved: ' savename]);
