% Load saved imagingInfo data
clear all; close all;
TR = 375; parc = 'fs86'; nparc = 86;

basedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(basedir);

% Load subjectInfo
load('subjectInfo_SUDcohort.mat'); 
subjkeys = subjectInfo.subjectkey;
nsubjs = height(subjectInfo);

load(['/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data/' sprintf('rsfMRI_%s_bp_gsr_%dsubj.mat', parc, nsubjs)]),

% Initialize variables to store results
subjectKeys = unique(subjectInfo.subjectkey); % Get unique subjects
numSubjects = numel(subjectKeys);
meanFD_perSubject = zeros(numSubjects, 1);
meanFD_nanoutlier_perSubject = zeros(numSubjects, 1);

% Calculate mean FD and mean FD nan outlier for each subject
for i = 1:numSubjects
    subjKey = subjectKeys{i}; % Current subject
    subjRows = imagingInfo(strcmp(imagingInfo.subjectkey, subjKey), :); % Rows corresponding to this subject

    % Mean FD across all scans for the subject
    meanFD_perSubject(i) = mean(subjRows.meanFD);

    % Mean FD with NaN outliers across all scans for the subject
    meanFD_nanoutlier_perSubject(i) = mean(subjRows.meanFD_nanoutlier, 'omitnan');
end

% Create a table for results
resultsTable = table(subjectKeys, meanFD_perSubject, meanFD_nanoutlier_perSubject, ...
    'VariableNames', {'subjectkey', 'FD_mean', 'FD_mean_nanoutlier'});

% Check subjectInfo and TE results in same order 
if ~all(strcmp(resultsTable.subjectkey,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD mean');
end 
 
% Save the results
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
save(fullfile(saveDir, 'MeanFD_PerSubject.mat'), 'resultsTable');


