%% Subcentroids: generate subject specific centroids to make subj-specific comparisons and an TE comparison matrix
function [sub_centroids,subjects, clusterNames, net9angle_Up, net9angle_Down] = generate_subcentroids(concTS,partition, ts_type, numClusters, subjInd,...
    scanInd,subj_scanInd,nscans,savedir)

% Set up
nparc = size(concTS,2);
clusterNames={};
net9angle_up = [];
net9angle_Down = []; 

% For each scan, compute subject-specifc centroids 
scanNums = unique(scanInd,'stable');
centroids = NaN(nscans,nparc,numClusters);
for scan = 1:nscans
    scanNum = scanNums(scan);
    disp(['Generating centroids for scan ',num2str(scan)]);
    centroids(scan,:,:) = GET_CENTROIDS(concTS(scanInd==scanNum,:),partition(scanInd==scanNum),numClusters);
    [clusterNames{scan},~,~,net9angle(scan,:,:)] = NAME_CLUSTERS_ANGLE(squeeze(centroids(scan,:,:)));
    [~,~,net9angle_Up(scan,:,:),net9angle_Down(scan,:,:)] = NAME_CLUSTERS_UP_DOWN(squeeze(centroids(scan,:,:)));
end

% Average across scans 
subjects = unique(subjInd,'stable'); nsubjs = length(subjects);
sub_centroids=NaN(length(subjects),nparc,numClusters);
for subj = 1:length(subjects)
    subjectInd = subjects(subj);
    sub_centroids(subj,:,:) = mean(centroids(subj_scanInd == subjectInd,:,:),1);
end

% save
save(fullfile(savedir,['subjcentroids_k',num2str(numClusters), '_',ts_type,'_fs', num2str(nparc),'_', num2str(nsubjs), '.mat']),'sub_centroids','subjects')
