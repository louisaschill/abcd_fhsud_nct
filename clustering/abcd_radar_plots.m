%% Radar Plots -- updated Louisa Schilling 12/24
% Calculate cosine similarity to Yeo networks (plus subcortex and cerebellum)
% Plot low and high amp activity in radar plots 
% Assign each cluster to a network based on max value 

% Dependencies: GET_CLUSTER_COLORS, NAME_CLUSTERS_UP_DOWN (modified from EJC: https://github.com/ejcorn/brain_states) 

clear all; close all;
parc = 'fs86';  ts_type = 'bp_gsr_gmnorm_exclout';
 
k = 4; % Set k to optimal value as found via abcd_repkmeans and abcd_elbow 
 
baseDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(baseDir);
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/'; 
figDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures';
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/radar'; mkdir(saveDir); 

% Load subject info 
load('subjectInfo_SUDcohort.mat'); 
subjkeys = subjectInfo.subjectkey; 
nsubjs = height(subjectInfo);

for numClusters = k
    load(fullfile(dataDir,['Partition_bp_',parc, '_k',num2str(numClusters),'_', ts_type, num2str(nsubjs) 'subj.mat']));

    overallNames = clusterNames;
    [nparc,numClusters] = size(centroids);

    YeoColors = [0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0;0 0 0];
    YeoColors = [YeoColors;YeoColors];

    [~,~,net9angle_Up,net9angle_Down] = NAME_CLUSTERS_UP_DOWN(centroids);
    YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN', 'SUB','CER'};
    numNets = numel(YeoNetNames);

    % Plot 
    clusterColors = GET_CLUSTER_COLORS(numClusters);

    clusterColors = hex2rgb(clusterColors);
    netAngle = linspace(0,2*pi,numNets+1);
    thetaNames = YeoNetNames; thetaNames{10} = '';

    if numClusters == 4 % hacky way to re-order clusters 
        [c,r]= corr(centroids);
        DMN_pos_ind = find(strcmp(clusterNames,'DMN+'));
        DMN_neg_ind = find(strcmp(clusterNames,'DMN-'));
        VIS_pos_ind = find(strcmp(clusterNames,'VIS+'));
        VIS_neg_ind = find(strcmp(clusterNames,'VIS-'));
        td_inds = [DMN_pos_ind,DMN_neg_ind,VIS_pos_ind,VIS_neg_ind];
    end

    f=figure('Position', [100 100 1200 400]);
    for k = 1:numClusters
        if numClusters == 4
            K = td_inds(k);
        else 
            K = k; 
        end 

        ax = subplot(1,numClusters,k,polaraxes); hold on

        polarplot(netAngle,[net9angle_Up(K,:) net9angle_Up(K,1)],'r','LineWidth', 2);
        polarplot(netAngle,[net9angle_Down(K,:) net9angle_Down(K,1)],'b','LineWidth', 2);

        thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
        rlim([0.0 0.7]);
        rticks([0.2 0.4 0.8]); rticklabels({'','0.4','0.8'});
        for L = 1:numNets
            ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
                YeoColors(L,:), ax.ThetaTickLabel{L});
        end
        set(ax,'FontSize',20);
        set(gca,'Fontname','calibri');
        title(overallNames{K},'Color',clusterColors(K,:),'FontSize',20);
    end
    saveas(f,fullfile(figDir,['Radar_plots_',parc,'_k',num2str(numClusters),'_', ts_type, '.png']));

end
