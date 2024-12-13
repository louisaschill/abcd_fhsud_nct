%% Louisa Schilling - updated Dec 2024 
% Generate subject specific centroids from partition for NCT analysis 
clear all; close all;

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(dataDir);
parc = 'fs86'; nparc=86; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout';

% Load subject info 
load([dataDir '/subjectInfo_SUDcohort.mat']); 
subjkeys = subjectInfo.subjectkey; 
nsubjs = height(subjectInfo);

% Load concTS 
load(['concTS_fs86_' ts_type '_' num2str(nsubjs) 'subj.mat']);
if strcmp(parc,'fs68'); concTS = concTS(:,19:86); end 

nscans = length(subj_scanInd);
savedir = fullfile('/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/subcentroids'); 
cd(savedir); 

% Loop through k range
for i = 1:length(k)
    numClusters = k(i);
    disp(['Starting k = ' num2str(numClusters)]);

    % Load parition
    resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results'; cd(resultsDir);
    load(fullfile(resultsDir,['Partition_bp_' parc, ...
        '_k',num2str(numClusters),'_', ts_type,  num2str(nsubjs), 'subj.mat']));
        
    % Calculate subject-specific centroids from partition
    [sub_centroids,subjects,clusterNames_sub, net9angle_Up, net9angle_Down] = abcd_generate_subcentroids(concTS,partition, ...
        ts_type, numClusters, subjInd, scanInd,subj_scanInd,nscans,savedir);
end 

%% Compare FH+ and FH- centroids 
YeoNetNames = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
numNets = numel(YeoNetNames);
clusterColors = GET_CLUSTER_COLORS(numClusters);
clusterColors = hex2rgb(clusterColors);
YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117; 80 80 255; 255 165 0] / 255;
YeoColors = [YeoColors;YeoColors];
netAngle = linspace(0,2*pi,numNets+1);
thetaNames = YeoNetNames; 

FHpos_centroids = mean(sub_centroids(subjectInfo.FHSUD == "FH+",:,:),1);
FHneg_centroids = mean(sub_centroids(subjectInfo.FHSUD == "FH-",:,:),1);

[clusterNamesFHpos,~,~,net9angleFHpos] = NAME_CLUSTERS_ANGLE(squeeze(FHpos_centroids)); 
[clusterNamesFHneg,~,~,net9angleFHneg] = NAME_CLUSTERS_ANGLE(squeeze(FHneg_centroids));
[~,~,net9angle_Up_FHpos,net9angle_Down_FHpos] = NAME_CLUSTERS_UP_DOWN(squeeze(FHpos_centroids)); 
[~,~,net9angle_Up_FHneg,net9angle_Down_FHneg] = NAME_CLUSTERS_UP_DOWN(squeeze(FHneg_centroids));

overallNames = clusterNamesFHpos;
inet9angle_Up = squeeze(net9angle_Up_FHpos);
inet9angle_Down = squeeze(net9angle_Down_FHpos);

f=figure("Position", [100 100 1200 400]);
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'k');
    polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
            YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];

overallNames = clusterNamesFHneg;
inet9angle_Up = squeeze(net9angle_Up_FHneg);
inet9angle_Down = squeeze(net9angle_Down_FHneg);

f=figure("Position", [100 100 1200 400]);
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'k');
    polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
            YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];

[rcent,~]=corr(squeeze(FHpos_centroids),squeeze(FHneg_centroids));

figure; imagesc(rcent);
xticks(1:numClusters); yticks(1:numClusters);
colormap(flipud(cbrewer2('RdYlBu')));
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('FH+ Centroids'); xlabel('FH- Centroids');
title('Correlation between FH of SUD Group-Specific Centroids');
caxis_bound = 1;
h = colorbar; ylabel(h,'R'); clim([-caxis_bound caxis_bound]); 
h.Ticks = [-caxis_bound 0 caxis_bound]; 
h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','calibri');

% Add text of r values
for i = 1:numClusters
    for j = 1:numClusters
        text(j, i, sprintf('%.4f', rcent(i, j)), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
    end
end

%% Compare centroids from each scanner type  
% Siemens 
Siemens_centroids = mean(sub_centroids(subjectInfo.manufacturer == "SIEMENS",:,:),1);
[clusterNamesSiemens,~,~,net9angleSiemens] = NAME_CLUSTERS_ANGLE(squeeze(Siemens_centroids)); 
[~,~,net9angle_Up_Siemens,net9angle_Down_Siemens] = NAME_CLUSTERS_UP_DOWN(squeeze(Siemens_centroids)); 

overallNames = clusterNamesSiemens;
inet9angle_Up = squeeze(net9angle_Up_Siemens);
inet9angle_Down = squeeze(net9angle_Down_Siemens);

f=figure("Position", [100 100 1200 400]);
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'k');
    polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'r');
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
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];

% GE  
GE_centroids = mean(sub_centroids(subjectInfo.manufacturer == "GE MEDICAL SYSTEMS",:,:),1);
[clusterNamesGE,~,~,net9angleGE]= NAME_CLUSTERS_ANGLE(squeeze(GE_centroids));
[~,~,net9angle_Up_GE,net9angle_Down_GE] = NAME_CLUSTERS_UP_DOWN(squeeze(GE_centroids));

overallNames = clusterNamesGE;
inet9angle_Up = squeeze(net9angle_Up_GE);
inet9angle_Down = squeeze(net9angle_Down_GE);

f=figure("Position", [100 100 1200 400]);
for K = 1:numClusters
    ax = subplot(1,numClusters,K,polaraxes); hold on
    polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'k');
    polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'r');
    thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
    rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
    for L = 1:numNets
        ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
            YeoColors(L,:), ax.ThetaTickLabel{L});
    end
    set(ax,'FontSize',10);
    title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
end
sgtitle('GE');
f.PaperUnits = 'inches';
f.PaperSize = [8 1.5];
f.PaperPosition = [0 0 8 1.5];

[rcent,~]=corr(squeeze(GE_centroids),squeeze(Siemens_centroids));

figure; imagesc(rcent);
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
