%% Louisa Schilling - updated Dec 2024
% run after choosing k (or a range of k)
% generate N partitions and compare pair-wise for mutual information 
% compute centroids and plot 

clear all; close all;
parc = 'fs86';  ts_type = 'bp_gsr_gmnorm_exclout';
 
basedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(basedir);
savedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';

% Load subject info 
load('subjectInfo_SUDcohort.mat'); 
subjkeys = subjectInfo.subjectkey; 
nsubjs = height(subjectInfo);

% Load concTS 
load(['concTS_fs86_' ts_type '_' num2str(nsubjs) 'subj.mat']);
if strcmp(parc, 'fs68') % select only cortical regions 
    concTS = concTS(:,19:86);
end 

%% inputs
distanceMethod = 'correlation'; % distance metric for clustering, we used correlation
nreps = 10;	% how many times to repeat clustering. will choose lowest error solution
maxI = 10000; % how many times to let kmeans try to converge before aborting rep
N = 10; % # of partitions to make and compare
[T,nparc] = size(concTS);

%% generate N partitions to be compared pair-wise for mutual information
for numClusters = 4 % can run single k or range   
    
    parts = NaN(T,N); % partitions will be stored here
    D = NaN(N,T,numClusters); % distance matrices will be stored here
    for i=1:N

        disp(['Clusters: ',num2str(numClusters),'. Starting kmeans number: ',num2str(i)]);
        [parts(:,i),D(i,:,:),~,~] = run_kmeans(concTS,numClusters,distanceMethod,nreps,maxI,subjInd,[]);

    end
    %% calculate adjusted mutual information for every pair of partitions    
    ami_results = NaN(N,N);   
    for i=1:N
        for j=1:N
            ami_results(i,j) = ami(parts(:,i),parts(:,j));
        end
    end
    
    % assess
    [m,td_inds] = max(sum(ami_results,1)); %ind corresponds to the partition which has the highest mutual information with all other partitions
    partition = parts(:,td_inds); % take partition that has most agreement with all other for further analysis
    
    f = figure;
    
    imagesc(ami_results); title(['Adjusted Mutal Information between Partitions k=',num2str(numClusters)]); colorbar;
    axis square; set(gca,'FontSize',8);
    f.PaperUnits = 'inches';
    f.PaperSize = [4 2];
    f.PaperPosition = [0 0 4 2];
    saveas(f,fullfile(savedir,['AMI_bp_', parc, '_k',num2str(numClusters),'_', ts_type, '_', num2str(nsubjs), '.jpg']));   
   
    %% compute centroids and plot
    centroids = GET_CENTROIDS(concTS,partition,numClusters);
    
    % name clusters based on alignment with Yeo resting state networks
    clusterNames = NAME_CLUSTERS_ANGLE(centroids);  % need to add a prior Yeo partition labels for your parcellation
    [clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(centroids);  % need to add a prior Yeo partition labels for your parcellation
    
    [c,r]= corr(centroids);

    if numClusters == 4
        % sort order DMN+, DMN-, VIS+, VIS-
        DMN_pos_ind = find(strcmp(clusterNames,'DMN+'));
        DMN_neg_ind = find(strcmp(clusterNames,'DMN-'));
        VIS_pos_ind = find(strcmp(clusterNames,'VIS+'));
        VIS_neg_ind = find(strcmp(clusterNames,'VIS-'));
        td_inds = [DMN_pos_ind,DMN_neg_ind,VIS_pos_ind,VIS_neg_ind];
        sorted = c(td_inds,td_inds);
        states_ordered = clusterNames(td_inds);
    else 
        td_inds = 1:5; 
    end

    % Centroids, 
    f = figure('Position', [100, 100, 800, 600]);
    imagesc(centroids(:,td_inds)); title('Centroids'); 
    xticks(1:numClusters); xticklabels(states_ordered);
    colormap('plasma'); axis square; colorbar; set(gca,'FontSize',15); 
    COLOR_TICK_LABELS(true,false,numClusters);
    xtickangle(90);

    saveas(f,fullfile(savedir, ['Centroids_bp_',parc, '_k',...
        num2str(numClusters),'_', ts_type,  '_', num2str(nsubjs), '.png']));
    
   
    f = figure('Position', [100, 100, 800, 600]);
    imagesc(sorted); title('centroid similarity'); h=colorbar;
    ylabel(h,'correlation strength'); caxis([-1 1]);
    colormap(flipud(cbrewer2('RdYlBu')));
    axis square; set(gca,'FontSize',15,'Fontname','Calibri');
    xticks(1:numClusters); yticks(1:numClusters);
    xticklabels(states_ordered); yticklabels(states_ordered); xtickangle(90);
    COLOR_TICK_LABELS(true,true,numClusters);

    saveas(f,fullfile(savedir, ['Centroids_similarity_bp_',parc, '_k',...
        num2str(numClusters),'_', ts_type,  '_', num2str(nsubjs), '.png']));
    
   
    %% save partition 
    save(fullfile(savedir,['Partition_bp_',parc, '_k',num2str(numClusters),...
        '_', ts_type, num2str(nsubjs), 'subj.mat']), 'parts', 'ami_results', ...
        'partition', 'clusterNames','centroids','td_inds','-v7.3');
end
