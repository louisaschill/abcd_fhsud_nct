% Louisa Schilling - updated 
% Global TE analysis 
clear all; close all;

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures/Final'; mkdir(figureDir);

parc = 'fs86'; k = 4; numClusters = k; 
ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg'; 

% Load subject info
load([dataDir '/subjectInfo_SUDcohort.mat']);
nsubjs = height(subjectInfo);

c = 1;

% Load optimal T from abcd_tsweep.mat
load([resultsDir '/optimalT_k4_c',num2str(c), '_' ts_type, '_sc_' sc_type '_noself_' num2str(nsubjs) '.mat'],'T'); 

% Load subject energy
load(fullfile(resultsDir, ['subjenergy_k',num2str(numClusters), '_T' num2str(T),...
    '_c',num2str(c),'_', ts_type, '_', sc_type, 'noselfSC_' parc '_' num2str(nsubjs) '.mat']));

% Rearrange cluster order to be DMN+, DMN-, VIS+, VIS-
if numClusters == 4 
    DMN_pos_ind = find(strcmp(clusterNames,'DMN+'));
    DMN_neg_ind = find(strcmp(clusterNames,'DMN-'));
    VIS_pos_ind = find(strcmp(clusterNames,'VIS+'));
    VIS_neg_ind = find(strcmp(clusterNames,'VIS-'));
    clusterOrder = [DMN_pos_ind,DMN_neg_ind,VIS_pos_ind,VIS_neg_ind];
end

% Check subject order 
assert(all(strcmp(subjkeys, subjectInfo.subjectkey)), 'Subject keys do not match.');

% Load FD
load(fullfile(dataDir, 'MeanFD_PerSubject.mat'), 'resultsTable');
% Check subjectInfo and TE results in same order 

if ~all(strcmp(resultsTable.subjectkey,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD mean');
end 
subjectInfo.FD_mean = resultsTable.FD_mean; 
subjectInfo.FD_mean_nanoutlier = resultsTable.FD_mean_nanoutlier; 

outlier = isoutlier(globalTE,"quartiles") | isnan(globalTE);

[h,p,c,s] = ttest2(subjectInfo.FD_mean(outlier == 1),subjectInfo.FD_mean(outlier == 0));

disp(['N = ' num2str(length(outlier(outlier==1))) ' subjects with outlier globalTE removed']);
subjectInfo = subjectInfo(~outlier,:);
globalTE = globalTE(~outlier); 
E_full = E_full(~outlier,:); 

%% ANCOVA 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'); 
anovaVarNames = {'Sex','Age','FHSUD','FD','Model','Income',...
    'Parent Education','Race','Prenatal','Parent MH', 'Puberty'};
interaction_term_sets = {{'Sex', 'FHSUD'},{'FHSUD','Income'},{'Sex','Puberty'}};
contVar = [2,4];
anovaVars = {subjectInfo.sex(inds)=='F', subjectInfo.age(inds),...
    subjectInfo.FHSUD(inds) == 'FH+', subjectInfo.FD_mean(inds),...
    subjectInfo.model(inds), subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds), subjectInfo.race(inds),...
    categorical(subjectInfo.subDuringPregAfter(inds)),...
    subjectInfo.parentMH(inds) == 1,...
    categorical(subjectInfo.pds_mod(inds))};

[all_results, stats, ~, ~, ~]  = run_ANOVA(globalTE(inds),anovaVars,contVar,...
    anovaVarNames,interaction_term_sets,'off');
disp(all_results)
[F_stats_sexsud,p_val_sexsud, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
[F_stats_sud, p_val_sud, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');

%% T-tests with Effect Size (Cohen's d)
pvals = []; 

% Helper function for Cohen's d
computeCohensD = @(x1, x2) (mean(x1) - mean(x2)) / sqrt(((std(x1)^2) + (std(x2)^2)) / 2);

% F-FH+ vs F-FH- 
group1 = globalTE(subjectInfo.sexFHSUD == 'F-FH+');
group2 = globalTE(subjectInfo.sexFHSUD == 'F-FH-');
[~,p,~,s] = ttest2(group1, group2); 
t = s.tstat; 
d = computeCohensD(group1, group2);
disp(['Global TE FH+ vs FH- females t-test: t = ', num2str(t,3), ', p = ', num2str(p,3), ', Cohen''s d = ', num2str(d,3)]);
pvals = [pvals p];

% M-FH+ vs M-FH- 
group1 = globalTE(subjectInfo.sexFHSUD == 'M-FH+');
group2 = globalTE(subjectInfo.sexFHSUD == 'M-FH-');
[~,p,~,s] = ttest2(group1, group2); 
t = s.tstat; 
d = computeCohensD(group1, group2);
disp(['Global TE FH+ vs FH- males t-test: t = ', num2str(t,3), ', p = ', num2str(p,4), ', Cohen''s d = ', num2str(d,3)]);
pvals = [pvals p];

% F vs M
group1 = globalTE(subjectInfo.sex == 'F');
group2 = globalTE(subjectInfo.sex == 'M');
[~,p,~,s] = ttest2(group1, group2); 
t = s.tstat; 
d = computeCohensD(group1, group2);
disp(['Global TE females vs males t-test: t = ', num2str(t,3), ', p = ', num2str(p,3), ', Cohen''s d = ', num2str(d,3)]);
pvals = [pvals p];

% Parent MH vs no parent MH
group1 = globalTE(subjectInfo.parentMH == 1);
group2 = globalTE(subjectInfo.parentMH == 0);
[~,p,~,s] = ttest2(group1, group2); 
t = s.tstat; 
d = computeCohensD(group1, group2);
disp(['Global TE parent MH vs no parent MH t-test: t = ', num2str(t,3), ', p = ', num2str(p,4), ', Cohen''s d = ', num2str(d,3)]);
pvals = [pvals p];

% Multiple comparison correction
pcorr = mafdr(pvals, 'bhfdr', 'true'); 
disp(pcorr);

%% Violin plots and FHD correlation 
clear g; 
f = figure('Position', [0 0 1500 650]);
g(1,1) = gramm('x', subjectInfo.sexFHSUD, 'y', globalTE, 'color', subjectInfo.FHSUD,...
    'subset', (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'));
g(1,1).stat_boxplot('width',0.15);
g(1,1).stat_violin('normalization','width','dodge', 0,'fill','edge');
g(1,1).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
g(1,1).set_names('x', 'sex:family history of SUD', 'y', 'global transition energy', 'color', '');
g(1,1).set_text_options('base_size', 20,'Font','calibri');
g(1,1).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);

g(1,2) = gramm('x', subjectInfo.familydensitySUD, 'y', globalTE,...
    'color', subjectInfo.sex);
g(1,2).geom_point('dodge', 0.2);
g(1,2).stat_glm();
g(1,2).axe_property('YGrid', 'on','GridColor',[0.5,0.5,0.5]);
g(1,2).set_text_options('base_size',20,'Font','calibri');
g(1,2).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]); 
g(1,2).set_names('x', 'family history density', 'y', 'global transition energy','color', '');
g.set_title('  ');  % to fix bug where axis is cut off 
g.draw(); 
exportgraphics(f, [figureDir '/Fig2ab_violin_fhd_global_k' num2str(numClusters) '.png']);

%% Correlation FHD x globalTE by sex 

inds = subjectInfo.sex == 'F'; 
[r1,p1] = corr(subjectInfo.familydensitySUD(inds),globalTE(inds),'type','Spearman','rows','complete'); 
inds = subjectInfo.sex == 'M'; 
[r2,p2] = corr(subjectInfo.familydensitySUD(inds),globalTE(inds),'type','Spearman','rows','complete'); 
pcorr = mafdr([p1, p2],'bhfdr','true'); 

disp(['Females: r = ' num2str(r1,2) ', p = ' num2str(p1,2)]);
disp(['Males: r = ' num2str(r2,2) ', p = ' num2str(p2,2)]);


%% Sex:FHSUD matrix
if k==5; clusterOrder = [1, 3, 5,4,2];end
clusterNamesOrdered = clusterNames(clusterOrder); 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');

anovaVarNames = {'Sex','Age','FHSUD','FD','MRI Model','Income Category',...
    'ParentEd','Race','In Utero Substances','Parental History Mental Health'};
interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
contVar = [2,4]; 
anovaVars = {subjectInfo.sex(inds)=='F',subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
    subjectInfo.FD_mean(inds),categorical(cellstr(subjectInfo.model(inds))),subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds),subjectInfo.race(inds), ...
    subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds)};

[all_results,stats] = run_ANOVA(E_full(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[F_sexsud,p_val_sexsud, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
p_val_sexsud_corr = mafdr(p_val_sexsud,'bhfdr','true'); 

% Sex*FHSUD 
p_sexsud_corr = reshape(p_val_sexsud_corr,[numClusters numClusters])';
p_sexsud = reshape(p_val_sexsud,[numClusters numClusters])';
fsexsud = reshape(squeeze(F_sexsud),[numClusters numClusters])';
p_sexsud = p_sexsud(clusterOrder,clusterOrder); 
p_sexsud_corr = p_sexsud_corr(clusterOrder,clusterOrder); 
fsexsud = fsexsud(clusterOrder,clusterOrder); 

l_caxis_bound = 0;  u_caxis_bound = max(max(fsexsud));
f=figure('Position', [100 100 700 600]); 
imagesc(fsexsud); 
colormap(cbrewer2('Greens')); fontsize = 25;
xticks(1:numClusters); xticklabels(clusterNamesOrdered); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNamesOrdered); axis square
ylabel('initial state'); xlabel('final state'); title('sex:family history of SUD')
sig_thresh = 0.05;
[y,x] = find(p_sexsud < sig_thresh);
text(x-.25,y+.2,'*','Color',[1,1,1],'Fontsize', 60);
[y,x] = find(p_sexsud_corr < sig_thresh);
text(x-.25,y+.2,'**','Color',[1,1,1],'Fontsize',60);
h = colorbar('FontName','calibri'); ylabel(h,'f-statistic');
clim([l_caxis_bound u_caxis_bound]);
h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound];
h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,2,'significant') round(u_caxis_bound,1,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
set(gca,'FontSize',25); set(gca,'TickLength',[0 0]); 
set(gca,'Fontname','calibri');
exportgraphics(f, [figureDir '/Fig2c_sexfhsud_matrix_global_k' num2str(numClusters) '.png']);

%% Sex-specific t-stat matrics 
[~,pf,~,s] = ttest2(E_full(subjectInfo.sexFHSUD == 'F-FH+',:),...
    E_full(subjectInfo.sexFHSUD == 'F-FH-',:)); tf = s.tstat; 
[~,pm,~,s] = ttest2(E_full(subjectInfo.sexFHSUD == 'M-FH+',:),...
    E_full(subjectInfo.sexFHSUD == 'M-FH-',:)); tm = s.tstat; 

lLim =  -3; uLim = 3.0; 
colorLim = [lLim uLim];

saveName = []; titletext = []; bhfdr = 1; 

% Within Females 
inds = subjectInfo.sexFHSUD == 'F-FH+' | subjectInfo.sexFHSUD == 'F-FH-';

[f,grpAvg1, grpAvg2] = plotMatricesDiff(E_full(inds,:), subjectInfo.sexFHSUD(inds),...
    numClusters,clusterNames,saveName, ["F-FH+", "F-FH-"],...
    'females only: FH+ vs FH-',pf,bhfdr,colorLim,clusterOrder,1,'t-statistic');
exportgraphics(f, [figureDir '/Fig2d_females_pairwise_global_k' num2str(numClusters) '.png']);

% Within Males 
inds = subjectInfo.sexFHSUD == 'M-FH+' | subjectInfo.sexFHSUD == 'M-FH-'; 
[f,grpAvg1, grpAvg2] = plotMatricesDiff(E_full(inds,:), subjectInfo.sexFHSUD(inds),...
    numClusters,clusterNames,saveName, ["M-FH+", "M-FH-"],'males only: FH+ vs FH-',...
    pm,bhfdr,colorLim,clusterOrder,1,'t-statistic');
exportgraphics(f, [figureDir '/Fig2d_males_pairwise_global_k' num2str(numClusters) '.png']);
