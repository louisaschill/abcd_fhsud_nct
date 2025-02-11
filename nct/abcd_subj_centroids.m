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
% YeoNetNames = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
% numNets = numel(YeoNetNames);
% clusterColors = GET_CLUSTER_COLORS(numClusters);
% clusterColors = hex2rgb(clusterColors);
% YeoColors = [137 72 155; 128 162 199; 91 153 72; 202 114 251;250 218 94; 218 166 86; 199 109 117; 80 80 255; 255 165 0] / 255;
% YeoColors = [YeoColors;YeoColors];
% netAngle = linspace(0,2*pi,numNets+1);
% thetaNames = YeoNetNames; 
% 
% FHpos_centroids = mean(sub_centroids(subjectInfo.FHSUD == "FH+",:,:),1);
% FHneg_centroids = mean(sub_centroids(subjectInfo.FHSUD == "FH-",:,:),1);
% 
% [clusterNamesFHpos,~,~,net9angleFHpos] = NAME_CLUSTERS_ANGLE(squeeze(FHpos_centroids)); 
% [clusterNamesFHneg,~,~,net9angleFHneg] = NAME_CLUSTERS_ANGLE(squeeze(FHneg_centroids));
% [~,~,net9angle_Up_FHpos,net9angle_Down_FHpos] = NAME_CLUSTERS_UP_DOWN(squeeze(FHpos_centroids)); 
% [~,~,net9angle_Up_FHneg,net9angle_Down_FHneg] = NAME_CLUSTERS_UP_DOWN(squeeze(FHneg_centroids));
% 
% overallNames = clusterNamesFHpos;
% inet9angle_Up = squeeze(net9angle_Up_FHpos);
% inet9angle_Down = squeeze(net9angle_Down_FHpos);
% 
% if numClusters == 4 % hacky way to re-order clusters
%     DMN_pos_ind = find(strcmp(clusterNames,'DMN+'));
%     DMN_neg_ind = find(strcmp(clusterNames,'DMN-'));
%     VIS_pos_ind = find(strcmp(clusterNames,'VIS+'));
%     VIS_neg_ind = find(strcmp(clusterNames,'VIS-'));
%     td_inds = [DMN_pos_ind,DMN_neg_ind,VIS_pos_ind,VIS_neg_ind];
% end
% 
% f=figure("Position", [100 100 1200 400]);
% for k = 1:numClusters
%     if numClusters == 4
%         K = td_inds(k);
%     else
%         K = k;
%     end
%     ax = subplot(2,numClusters,k,polaraxes); hold on
%     polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'r');
%     polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'b');
%     thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
%     rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
%     for L = 1:numNets
%         ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
%             YeoColors(L,:), ax.ThetaTickLabel{L});
%     end
%     set(ax,'FontSize',10);
%     title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
% end
% 
% overallNames = clusterNamesFHneg;
% inet9angle_Up = squeeze(net9angle_Up_FHneg);
% inet9angle_Down = squeeze(net9angle_Down_FHneg);
% 
% for k = 1:numClusters
%     if numClusters == 4
%         K = td_inds(k);
%     else
%         K = k;
%     end
%     ax = subplot(2,numClusters,k+4,polaraxes); hold on
%     polarplot(netAngle,[inet9angle_Up(K,:) inet9angle_Up(K,1)],'r');
%     polarplot(netAngle,[inet9angle_Down(K,:) inet9angle_Down(K,1)],'b');
%     thetaticks(rad2deg(netAngle)); thetaticklabels(thetaNames);
%     rticks([0.4 0.8]); rticklabels({'0.4','0.8'});
%     for L = 1:numNets
%         ax.ThetaTickLabel{L} = sprintf('\\color[rgb]{%f,%f,%f}%s', ...
%             YeoColors(L,:), ax.ThetaTickLabel{L});
%     end
%     set(ax,'FontSize',10);
%     title(overallNames{K},'Color',clusterColors(K,:),'FontSize',12);
% end
% sgtitle('FH+ (top) and FH- (bottom) Centroids');
% f.PaperUnits = 'inches';
% f.PaperSize = [8 1.5];
% f.PaperPosition = [0 0 8 1.5];
% 
% [rcent,~]=corr(squeeze(FHpos_centroids),squeeze(FHneg_centroids));
% 
% figure; imagesc(rcent);
% xticks(1:numClusters); yticks(1:numClusters);
% colormap(flipud(cbrewer2('RdYlBu')));
% xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
% COLOR_TICK_LABELS(true,true,numClusters);
% ylabel('FH+ Centroids'); xlabel('FH- Centroids');
% title('Correlation between FH of SUD Group-Specific Centroids');
% caxis_bound = 1;
% h = colorbar; ylabel(h,'R'); clim([-caxis_bound caxis_bound]); 
% h.Ticks = [-caxis_bound 0 caxis_bound]; 
% h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
% COLOR_TICK_LABELS(true,true,numClusters);
% set(gca,'FontSize',8);
% set(gca,'TickLength',[0 0]);
% set(gca,'Fontname','calibri');
% 
% % Add text of r values
% for i = 1:numClusters
%     for j = 1:numClusters
%         text(j, i, sprintf('%.4f', rcent(i, j)), 'HorizontalAlignment', 'center', ...
%              'VerticalAlignment', 'middle', 'FontSize', 10, 'Color', 'k', 'FontWeight', 'bold');
%     end
% end
