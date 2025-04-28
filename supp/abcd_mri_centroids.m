%% Comparing centroids across scanners 
clear; close all;

% Directories and parameters
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results';
saveDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures/Reviews';

sc_type = 'avg'; % SC type: 'avg' (group average) or 'indiv' (individual)
parc = 'fs86'; nparc = 86; % Parcellation type and number of parcels
numClusters = 4; k = 4; % Number of clusters/k
ts_type = 'bp_gsr_gmnorm_exclout'; % rsfmri timeseries type

% Load subject info
cd(dataDir);
load('subjectInfo_SUDcohort.mat');
nsubjs = height(subjectInfo);

c = 1; % Constant parameter for energy calculation

% Load subject centroids
load([resultsDir, '/subcentroids/subjcentroids_k', num2str(numClusters), ...
    '_', ts_type, '_', parc, '_', num2str(nsubjs), '.mat']);

% Load partition for cluster names
cd(resultsDir);
load(['Partition_bp_fs86', '_k', num2str(numClusters), '_', ts_type, ...
    num2str(nsubjs), 'subj.mat'],'clusterNames');

%%
YeoNetNames = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
numNets = numel(YeoNetNames);
clusterColors = GET_CLUSTER_COLORS(numClusters);
clusterColors = hex2rgb(clusterColors);
YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117; 80 80 255; 255 165 0] / 255;
YeoColors = [YeoColors;YeoColors];
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; 

if numClusters == 4 % hacky way to re-order clusters
    DMN_pos_ind = find(strcmp(clusterNames,'DMN+'));
    DMN_neg_ind = find(strcmp(clusterNames,'DMN-'));
    VIS_pos_ind = find(strcmp(clusterNames,'VIS+'));
    VIS_neg_ind = find(strcmp(clusterNames,'VIS-'));
    td_inds = [DMN_pos_ind,DMN_neg_ind,VIS_pos_ind,VIS_neg_ind];
end

% Siemens 
Siemens_centroids = mean(sub_centroids(subjectInfo.manufacturer == "SIEMENS",:,:),1);
[clusterNamesSiemens,~,~,net9angleSiemens] = NAME_CLUSTERS_ANGLE(squeeze(Siemens_centroids)); 
[~,~,net9angle_Up_Siemens,net9angle_Down_Siemens] = NAME_CLUSTERS_UP_DOWN(squeeze(Siemens_centroids)); 

overallNames = clusterNamesSiemens;
inet9angle_Up = squeeze(net9angle_Up_Siemens);
inet9angle_Down = squeeze(net9angle_Down_Siemens);

f=figure("Position", [100 100 1200 400]);
for k = 1:numClusters
    if numClusters == 4
        K = td_inds(k);
    else
        K = k;
    end
    ax = subplot(2,numClusters,k,polaraxes); hold on
    polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'r');
    polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'b');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
            YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
sgtitle('Siemens');

% GE  
GE_centroids = mean(sub_centroids(subjectInfo.manufacturer == "GE MEDICAL SYSTEMS",:,:),1);
[clusterNamesGE,~,~,net9angleGE]= NAME_CLUSTERS_ANGLE(squeeze(GE_centroids));
[~,~,net9angle_Up_GE,net9angle_Down_GE] = NAME_CLUSTERS_UP_DOWN(squeeze(GE_centroids));

overallNames = clusterNamesGE;
inet9angle_Up = squeeze(net9angle_Up_GE);
inet9angle_Down = squeeze(net9angle_Down_GE);

for k = 1:numClusters
    if numClusters == 4
        K = td_inds(k);
    else
        K = k;
    end
    ax = subplot(2,numClusters,k+4,polaraxes); hold on
    polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'r');
    polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'b');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
            YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
sgtitle('Siemens (top) and GE (bottom) Centroids');
%f.PaperUnits = 'inches';
%f.PaperSize = [8 1.5];
%f.PaperPosition = [0 0 8 1.5];

saveas(f,fullfile(figureDir,'radarplots_MRIscanners.png'));

[rcent,~]=corr(squeeze(GE_centroids),squeeze(Siemens_centroids));
f= figure; imagesc(rcent);
xticks(1:numClusters); yticks(1:numClusters); 
colormap(flipud(cbrewer2('RdYlBu')));
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Siemens Centroids'); xlabel('GE Centroids');
title('Correlation between MRI Model-Specific Centroids');
caxis_bound = 1;
h = colorbar; ylabel(h,'R'); clim([-caxis_bound caxis_bound]);
h.Ticks = [-caxis_bound 0 caxis_bound];
h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',12);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','calibri');

% Add text of r values
for i = 1:numClusters
    for j = 1:numClusters
        text(j, i, sprintf('%.4f', rcent(i, j)), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
    end
end

saveas(f,fullfile(figureDir,'centroid_similarity_MRIscanners.png'));
