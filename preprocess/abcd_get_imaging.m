% Louisa Schilling - last updated 12/24 

% This script pre-processes and compiles imaging data from all subjects of
% the ABCD dataset for baseline visit with available data
% into a consolidated table for analysis. 

% Output: saves table imagingInfo as a .mat file

%% 
clear all; close all;
TR = 375; parc = 'fs86'; nparc = 86;

basedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(basedir);

% Load subjectInfo
load('subjectInfo_SUDcohort.mat'); 
subjkeys = subjectInfo.subjectkey;
nsubjs = height(subjectInfo);

% Load GM mean signal data
gm_means = readtable('abcd_wb_wm_csf_mean.txt', 'ReadVariableNames', false);
gm_means = gm_means(:,1:3);
gm_means.Properties.VariableNames = {'subject', 'scan', 'gm_mean'};
gm_means = gm_means(ismember(gm_means.subject, subjkeys), :);

dataDir = '/Users/louisaschilling/Desktop/Datasets/ABCD/Data/';
cd(dataDir);

% Load MRI information
mri = readtable([dataDir 'non_imaging/Release 5.1/core/imaging/mri_y_adm_info.csv']);
mri = extractInstrument(mri, {'src_subject_id', 'eventname', 'mri_info_manufacturer', ...
    'mri_info_manufacturersmn', 'mri_info_deviceserialnumber'}, ...
    {'subjectkey', 'eventname', 'manufacturer', 'model', 'serial_number'}, ...
    'baseline_year_1_arm_1', subjkeys);

% Prepare imagingInfo table
imagingInfo = table();
scanCount = 0;

% Pre-index available files for faster access
rsfmriFiles = dir(fullfile(dataDir, 'imaging', 'rsfMRI', '*.mat'));
rsfmriFileMap = containers.Map();
for f = 1:length(rsfmriFiles)
    rsfmriFileMap(rsfmriFiles(f).name) = fullfile(rsfmriFiles(f).folder, rsfmriFiles(f).name);
end

% Process rsfMRI data
dataDir = '/Users/louisaschilling/Desktop/Datasets/ABCD/Data/imaging/rsfMRI';
cd(dataDir);
n_excl = 0; 

for i = 1:nsubjs
    subject = subjkeys{i};
    subjData = gm_means(strcmp(gm_means.subject, subject), :);
    subjMri = mri(strcmp(mri.subjectkey, subject), :);

    subjScanCount = 0; % Track scans per subject
    for s = 1:4
        scanname = sprintf('bld00%d', s);
        ts_filename = sprintf('%s_%s_rest_mc_skip_residc_interp_FDRMS0.3_DVARS50_bp_0.009_0.08_%s_ts.mat', subject, scanname, parc);

        % Check if file exists
        if isKey(rsfmriFileMap, ts_filename)
            tsData = load(rsfmriFileMap(ts_filename), 'ts'); % Load timeseries data

            %if size(tsData.ts, 1) == TR % Check TR
                scanCount = scanCount + 1;
                subjScanCount = subjScanCount + 1;

                % Get GM mean
                gmData = subjData(strcmp(subjData.scan, scanname), :);
                if isempty(gmData)
                    continue;
                end
                gm_mean = gmData.gm_mean;

                % Normalize time series by GM mean
                ts_GMnorm = tsData.ts / gm_mean;

                % Load FD and Outlier mask
                FD = load(sprintf('%s_%s_FD.txt', subject, scanname));
                outlierMask = load(sprintf('%s_%s_FDRMS0.3_DVARS50_motion_outliers.txt', subject, scanname));
                outlierMask = (outlierMask == 0); % Mark outliers as 1

                % Populate imagingInfo
                imagingInfo = [imagingInfo; table(...
                    scanCount, ...
                    i, ...
                    string(subject), ... 
                    string(scanname), ... 
                    TR, ...
                    subjScanCount, ...
                    categorical(subjMri.manufacturer(1)), ...
                    categorical(subjMri.model(1)), ...
                    categorical(subjMri.serial_number(1)), ...
                    gm_mean, ...
                    {tsData.ts}, ... 
                    {ts_GMnorm}, ...
                    {FD}, ...
                    {outlierMask}, ...
                    'VariableNames', {'scanInd', 'subjInd', 'subjectkey', 'scanName', 'TR', ...
                    'subjScanNum', 'manufacturer', 'model', 'serial_number', ...
                    'meanGM', 'ts', 'ts_gmnorm', 'FD', 'outliermask'})];

        end
    end
end

% Replace Outliers with NaNs in Time Series
imagingInfo.ts_nanoutlier = imagingInfo.ts;
imagingInfo.ts_gmnorm_nanoutlier = imagingInfo.ts_gmnorm;
imagingInfo.FD_nanoutlier = imagingInfo.FD;

for row = 1:height(imagingInfo)
    mask = imagingInfo.outliermask{row};
    imagingInfo.ts_nanoutlier{row}(mask, :) = NaN;
    imagingInfo.ts_gmnorm_nanoutlier{row}(mask, :) = NaN;
    imagingInfo.FD_nanoutlier{row}(mask) = NaN;
end

% Compute Mean FD
imagingInfo.meanFD = cellfun(@mean, imagingInfo.FD);
imagingInfo.meanFD_nanoutlier = cellfun(@nanmean, imagingInfo.FD_nanoutlier);

% Save Processed Data
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
save(fullfile(saveDir, sprintf('rsfMRI_%s_bp_gsr_%dsubj.mat', parc, nsubjs)), 'imagingInfo', '-v7.3');
disp('imagingInfo saved.');
