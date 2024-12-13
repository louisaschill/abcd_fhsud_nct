%% Louisa Schilling - updated Dec 2024  
% Figure 3: Network level 
clear all; close all; 

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures'; 

parc = 'fs86'; numClusters = 4; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg'; 

clusterOrder = [1,2,4,3]; % put in order: DMN+, DMN-, VIS+, VIS-

% Load subject info
load([dataDir '/subjectInfo_SUDcohort.mat']);
nsubjs = height(subjectInfo);

% Load optimal T from abcd_tsweep.mat
c = 1;
load([resultsDir '/optimalT_k',num2str(numClusters),'_c',num2str(c), '_' ts_type, '_sc_' sc_type '_noself_' num2str(nsubjs) '.mat'],'T'); 

% Load subject energy
load(fullfile(resultsDir, ['subjenergy_k',num2str(numClusters), '_T' num2str(T),...
    '_c',num2str(c),'_', ts_type, '_', sc_type, 'noselfSC_' parc '_' num2str(nsubjs) '.mat']));

% Check subjectInfo and TE results in same order 
if ~all(strcmp(subjkeys,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD mean');
end 

% Network names and assignments 
networks = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
YeoLUT = readtable('fs86_to_yeo.csv'); YeoLUT = table2array(YeoLUT);
numNetworks= 9;

% Network TE calc 
E_network = nan(length(subjkeys), length(clusterNames)^2, numNetworks);
for i = 1:numNetworks
    network = E_region(:,:,YeoLUT == i); E_network(:, :,i) = sum(network,3);
end
networkTE = squeeze(mean(E_network,2));

% Load FD
load(fullfile(dataDir, 'MeanFD_PerSubject.mat'), 'resultsTable');

% Check subjectInfo and TE results in same order 
if ~all(strcmp(resultsTable.subjectkey,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD mean');
end 
subjectInfo.FD_mean = resultsTable.FD_mean; 
subjectInfo.FD_mean_nanoutlier = resultsTable.FD_mean_nanoutlier; 

outlier = isoutlier(globalTE,"median"); % Remove outliers
subjectInfo = subjectInfo(~outlier,:);
networkTE = networkTE(~outlier,:);  E_network = E_network(~outlier,:,:); 

%% ANCOVA 
subjectInfo.model = categorical(cellstr(subjectInfo.model));
% Variables for ANCOVA model 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'); %not FH+/-
anovaVarNames = {'Sex','Age','FHSUD','FD','MRI Model','Income Category',...
    'Parent Edu Category','Race','In Utero Substances',...
    'Parental History Mental Health'};
interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
contVar = [2,4]; 
anovaVars = {subjectInfo.sex(inds)=='F', subjectInfo.age(inds),...
    subjectInfo.FHSUD(inds) == 'FH+', subjectInfo.FD_mean(inds),...
    subjectInfo.model(inds), subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds), subjectInfo.race(inds),...
    subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds)};

[all_results] = run_ANOVA(networkTE(inds,:),anovaVars,contVar,...
    anovaVarNames,interaction_term_sets,'off');
[F_stats, p_val, ~, ~, ~, sig_ind, ~] = extract_anova_results(all_results,'Sex:FHSUD');
p_val_corr = mafdr(p_val, 'BHFDR',true);

disp('Network TE Sex:FHSUD sig :'); disp(networks(sig_ind))
disp('Network TE Sex:FHSUD sig post-FDR:'); disp(networks(p_val_corr <= 0.05));

disp(['VIS TE sex:FHSUD : F = ', num2str(F_stats(1),3), ' p = ', num2str(p_val(1),4), ', pFDR = ' num2str(p_val_corr(1),4)]);  
disp(['DAT TE sex:FHSUD : F = ', num2str(F_stats(3),3), ' p = ', num2str(p_val(3),4), ', pFDR = ' num2str(p_val_corr(3),4)]);  
disp(['VAT TE sex:FHSUD : F = ', num2str(F_stats(4),3), ' p = ', num2str(p_val(4),4), ', pFDR = ' num2str(p_val_corr(4),4)]);  
disp(['DMN TE sex:FHSUD : F = ', num2str(F_stats(7),3), ' p = ', num2str(p_val(7),4), ', pFDR = ' num2str(p_val_corr(7),4)]);  
disp(['CER TE sex:FHSUD : F = ', num2str(F_stats(9),3), ' p = ', num2str(p_val(9),4), ', pFDR = ' num2str(p_val_corr(9),4)]);  

[F_stats, p_val, ~, eta_squared, partial_eta_squared, sig_ind, ~] = extract_anova_results(all_results, 'FHSUD');
p_val_corr = mafdr(p_val, 'BHFDR',true);
disp('Network TE FHSUD sig :'); disp(networks(sig_ind))
disp('Network TE FHSUD sig post-FDR:');  disp(networks(p_val_corr <= 0.05));

%% T-test p-values and correlations 
pvals = []; 
[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'F-FH+',3), networkTE(subjectInfo.sexFHSUD == 'F-FH-',3)); t= s.tstat; 
disp(['DAT TE females ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'M-FH+',3), networkTE(subjectInfo.sexFHSUD == 'M-FH-',3)); t= s.tstat; 
disp(['DAT TE males ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'F-FH+',4), networkTE(subjectInfo.sexFHSUD == 'F-FH-',4)); t= s.tstat; 
disp(['VAT TE females ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'M-FH+',4), networkTE(subjectInfo.sexFHSUD == 'M-FH-',4)); t= s.tstat; 
disp(['VAT TE males ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'F-FH+',7), networkTE(subjectInfo.sexFHSUD == 'F-FH-',7)); t= s.tstat; 
disp(['DMN TE females ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'M-FH+',7), networkTE(subjectInfo.sexFHSUD == 'M-FH-',7)); t= s.tstat; 
disp(['DMN TE males ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);  pvals = [pvals, p]; 


[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'F-FH+',1), networkTE(subjectInfo.sexFHSUD == 'F-FH-',1)); t= s.tstat; 
disp(['VIS TE females ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'M-FH+',1), networkTE(subjectInfo.sexFHSUD == 'M-FH-',1)); t= s.tstat; 
disp(['VIS TE males ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);  pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'F-FH+',9), networkTE(subjectInfo.sexFHSUD == 'F-FH-',9)); t= s.tstat; 
disp(['CER TE females ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);   pvals = [pvals, p]; 

[~,p,~,s] = ttest2(networkTE(subjectInfo.sexFHSUD == 'M-FH+',9), networkTE(subjectInfo.sexFHSUD == 'M-FH-',9)); t= s.tstat; 
disp(['CER TE males ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);  pvals = [pvals, p]; 

pcorr = mafdr(pvals, 'BHFDR', 'true');
disp(['Corrected pvals: ', num2str(pcorr,3)]); 

pvals = [];
for n = [3,4,7] 
    inds = subjectInfo.sex == 'F'; 
    [r1,p1] = corr(subjectInfo.familydensitySUD(inds),networkTE(inds,n),'type','Spearman');
    inds = subjectInfo.sex == 'M'; 
    [r2,p2] = corr(subjectInfo.familydensitySUD(inds),networkTE(inds,n),'type','Spearman');
    pvals = [pvals, p1, p2]; 
    disp([networks{n} ' x FHD Females: r = ' num2str(r1) ',p = ' num2str(p1)])
    disp([networks{n} ' x FHD Males: r = ' num2str(r2) ',p = ' num2str(p2)])
end 
pcorr = mafdr(pvals, 'bhfdr', 'true');
disp(['Corrected pvals: ', num2str(pcorr,3)]); 

%% Bar plots: network F-stats (FHSUD, FHSUD:Sex)
networks = categorical({'vis','som','dat','vat','lim','fpn','dmn','sub','cer'});
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'); %not FH+/-

% Variables for ANCOVA model 
[all_results] = run_ANOVA(networkTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[F_stats, p_val, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');
p_val_corr = mafdr(p_val, 'bhfdr', true);
sig_sud= makeSigCategory(p_val, p_val_corr); 

[F_stats_sexsud, p_val, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
all_F = [F_stats, F_stats_sexsud];

% PLOT SUD F STATS
clear g; f= figure('Position', [100,100,1500,700]);
g(1,1) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sud);
g(1,1).set_order_options('x',{'cer', 'sub', 'dmn', 'fpn', 'lim', 'vat','dat','som','vis'}); 
g(1,1).coord_flip();
g(1,1).geom_bar();
g(1,1).set_text_options('base_size', 25,'font','calibri');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],'YLim', [-0.05, max(all_F(:))+0.2]);
g(1,1).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706]]);
g(1,1).set_names('x','','y', 'f-statistic: family history of SUD', 'color', 'Significance');

[F_stats, p_val, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
p_val_corr = mafdr(p_val, 'bhfdr', true);
sig_sexsud= makeSigCategory(p_val, p_val_corr); 

% PLOT SEX*SUD F STATS
g(1,2) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sexsud);
g(1,2).set_order_options('x',{'cer','sub','dmn','fpn','lim','vat','dat','som','vis'}); 
g(1,2).coord_flip(); 
g(1,2).geom_bar();
g(1,2).set_text_options('base_size', 25,'font','calibri');
g(1,2).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],'YLim', [-0.05, max(all_F(:))+0.2]);
g(1,2).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
g(1,2).set_names('x','','y', 'f-statistic: sex * family history of SUD', 'color', 'Significance');
g(1,2).no_legend(); 
g.draw();
exportgraphics(f, [figureDir '/Fig4ab_bar_network_k ' num2str(numClusters) '.png']);

%% ETA Squared 
% networks = categorical({'vis','som','dat','vat','lim','fpn','dmn','sub','cer'});
% [all_results, stats, overall_eta_squared,overall_r_squared] = run_ANOVA(networkTE(inds,:), anovaVars, ...
%     contVar, anovaVarNames, interaction_term_sets, 'off');
% 
% % Create a table with network names and overall eta squared
% overall_effect_table = table;
% overall_effect_table.Network = networks';
% overall_effect_table.Overall_Eta_Squared = overall_eta_squared;
% overall_effect_table.Overall_R_Squared = overall_r_squared;
% 
% % Display the table
% disp(overall_effect_table);
% 
% 
% [F_stats, p_val, ~, eta_sq, partial_eta_sq_sud] = extract_anova_results(all_results, 'FHSUD');
% p_val_corr = mafdr(p_val, 'bhfdr', true);
% sig_sud= makeSigCategory(p_val, p_val_corr); 
% 
% D = table; 
% D.network = networks'; 
% D.Fstats = F_stats; 
% D.p_values = p_val; 
% D.p_values_corr = p_val_corr; 
% D.eta_sq = eta_sq;
% D.partial_eta_sq_sud = partial_eta_sq_sud;
% disp(D);
% 
% [~, ~, ~, eta_sq_sexsud, partial_eta_sq_sexsud] = extract_anova_results(all_results, 'Sex:FHSUD');
% %all_eta = [partial_eta_sq_sexsud, partial_eta_sq_sud];
% 
% % PLOT SUD F STATS
% clear g; f= figure('Position', [100,100,1500,700]);
% g(1,1) = gramm('x', categorical(networks), 'y', partial_eta_sq_sud, 'color', sig_sud);
% g(1,1).set_order_options('x',{'cer', 'sub', 'dmn', 'fpn', 'lim', 'vat','dat','som','vis'}); 
% g(1,1).coord_flip();
% g(1,1).geom_bar();
% g(1,1).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
% g(1,1).set_text_options('base_size', 25,'font','calibri');
% g(1,1).set_names('x','','y', 'partial eta sq: family history of SUD', 'color', 'Significance');
% 
% [F_stats, p_val, ~, eta_sq_sexsud, partial_eta_sq_sexsud] = extract_anova_results(all_results, 'Sex:FHSUD');
% p_val_corr = mafdr(p_val, 'bhfdr', true);
% sig_sexsud= makeSigCategory(p_val, p_val_corr); 
% 
% D = table; 
% D.network = networks'; 
% D.Fstats = F_stats; 
% D.p_values = p_val; 
% D.p_values_corr = p_val_corr; 
% D.eta_sq = eta_sq_sexsud;
% D.partial_eta_sq_sud = partial_eta_sq_sexsud;
% disp(D);
% 
% % PLOT SEX*SUD F STATS
% g(1,2) = gramm('x', categorical(networks), 'y', partial_eta_sq_sexsud, 'color', sig_sexsud);
% g(1,2).set_order_options('x',{'cer','sub','dmn','fpn','lim','vat','dat','som','vis'}); 
% g(1,2).coord_flip(); 
% g(1,2).geom_bar();
% g(1,2).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
% g(1,2).set_text_options('base_size', 25,'font','calibri');
% g(1,2).set_names('x','','y', 'partial eta sq: sex * family history of SUD', 'color', 'Significance');
% g(1,2).no_legend(); 
% g.draw();

%% Violin Plots
FHSUD = subjectInfo.FHSUD; sexFHSUD = subjectInfo.sexFHSUD;

f = figure('Position', [0 0 1600 500]); clear g; 
g(1,1) = gramm('x', sexFHSUD, 'y', networkTE(:,7), 'color', FHSUD,...
    'subset', subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
g(1,1).stat_boxplot('width',0.15);
g(1,1).stat_violin('normalization','width','dodge', 0,'fill','edge');
g(1,1).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
g(1,1).set_names('x', '', 'y', 'dmn transition energy', 'color', '');
g(1,1).set_text_options('base_size', 20,'Font','calibri');
g(1,1).no_legend();
g(1,1).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);

g(1,2) = gramm('x', sexFHSUD, 'y', networkTE(:,3), 'color', FHSUD,...
    'subset', subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
g(1,2).stat_boxplot('width',0.15);
g(1,2).stat_violin('normalization','width','dodge', 0,'fill','edge');
g(1,2).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
g(1,2).set_names('x', '', 'y', 'dat transition energy', 'color', '');
g(1,2).set_text_options('base_size', 20,'Font','calibri');
g(1,2).no_legend();
g(1,2).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);

g(1,3) = gramm('x', sexFHSUD, 'y', networkTE(:,4), 'color', FHSUD,...
    'subset', subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
g(1,3).stat_boxplot('width',0.15);
g(1,3).stat_violin('normalization','width','dodge', 0,'fill','edge');
g(1,3).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
g(1,3).set_names('x', '', 'y', 'vat transition energy', 'color', '');
g(1,3).set_text_options('base_size', 20,'Font','calibri');
g(1,3).no_legend();
g(1,3).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);

g.set_title('  '); g.draw(); 
exportgraphics(f, [figureDir 'Fig4cde_networkviolin_k' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% FHD plots 
f = figure('Position', [0 0 1600 450]); clear g; 
g(1,1) = gramm('x', subjectInfo.familydensitySUD,'y', networkTE(:,7), 'color', subjectInfo.sex);
g(1,1).geom_point(); 
g(1,1).stat_glm();
g(1,1).set_names('x', 'family history density', 'y', 'dmn transition energy','color', 'Sex');
g(1,1).set_text_options('base_size', 20,'font','calibri');
g(1,1).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g(1,1).no_legend();

g(1,2) = gramm('x', subjectInfo.familydensitySUD,'y', networkTE(:,3), 'color', subjectInfo.sex);
g(1,2).geom_point();
g(1,2).stat_glm();
g(1,2).set_names('x', 'family history density', 'y', 'dat transition energy','color', 'Sex');
g(1,2).set_text_options('base_size', 20,'font','calibri');
g(1,2).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g(1,2).no_legend();

g(1,3) = gramm('x', subjectInfo.familydensitySUD,'y', networkTE(:,4), 'color', subjectInfo.sex);
g(1,3).geom_point();
g(1,3).stat_glm();
g(1,3).set_names('x', 'family history density', 'y', 'vat transition energy','color', 'Sex');
g(1,3).set_text_options('base_size', 20,'font','calibri');
g(1,3).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g(1,3).no_legend();

g.set_title('  '); g.draw(); 
exportgraphics(f, [figureDir '/Fig4fgh_networkFHD_k' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% Network pairwise matrices with tstat pvals 
lLim = -3; uLim = 3; colorLim = [lLim uLim]; 
saveName = []; titletext = []; bhfdr = 1; 
clusterOrder = [1,2,4,3];
 
% DAT males 
n = 3;
[~,p,~,~] = ttest2(squeeze(E_network(subjectInfo.sexFHSUD=='M-FH+',:,n)),...
    squeeze(E_network(subjectInfo.sexFHSUD=='M-FH-',:,n))); 
f = plotMatricesDiff(squeeze(E_network(subjectInfo.sex=='M',:,n)),...
    subjectInfo.sexFHSUD(subjectInfo.sex=='M'),...
    numClusters,clusterNames,saveName, ["M-FH+", "M-FH-"],...
    'males only: FH+ vs FH-',p,bhfdr,colorLim,clusterOrder,1,...
    'FH+ vs. FH- T-Statistic');
exportgraphics(f, [figureDir '/Fig4ia_malesDATpairwise_k' num2str(numClusters) '.png']);

% DAT females 
[~,p,~,~] = ttest2(squeeze(E_network(subjectInfo.sexFHSUD=='F-FH+',:,n)),...
    squeeze(E_network(subjectInfo.sexFHSUD=='F-FH-',:,n))); 
f = plotMatricesDiff(squeeze(E_network(subjectInfo.sex=='F',:,n)),...
    subjectInfo.sexFHSUD(subjectInfo.sex=='F'),...
    numClusters,clusterNames,saveName, ["F-FH+", "F-FH-"],...
    'females only: FH+ vs FH-',p,bhfdr,colorLim,clusterOrder,1,...
    'FH+ vs. FH- T-Statistic');
exportgraphics(f, [figureDir '/Fig4ib_femalesDATpairwise_k' num2str(numClusters) '.png']);

% VAT males 
n = 4;
[~,p,~,~] = ttest2(squeeze(E_network(subjectInfo.sexFHSUD=='M-FH+',:,n)),...
    squeeze(E_network(subjectInfo.sexFHSUD=='M-FH-' ,:,n))); 
f = plotMatricesDiff(squeeze(E_network(subjectInfo.sex=='M',:,n)), ...
    subjectInfo.sexFHSUD(subjectInfo.sex=='M'),...
    numClusters,clusterNames,saveName, ["M-FH+", "M-FH-"],...
    'males only: FH+ vs FH-',p,bhfdr,colorLim,clusterOrder,1,...
    'FH+ vs. FH- T-Statistic');
exportgraphics(f, [figureDir '/Fig4ja_malesVATpairwise_k' num2str(numClusters) '.png']);

% VAT females 
[~,p,~,~] = ttest2(squeeze(E_network(subjectInfo.sexFHSUD=='F-FH+',:,n)),...
    squeeze(E_network(subjectInfo.sexFHSUD=='F-FH-',:,n))); 
f = plotMatricesDiff(squeeze(E_network(subjectInfo.sex=='F',:,n)),...
    subjectInfo.sexFHSUD(subjectInfo.sex=='F'),...
    numClusters,clusterNames,saveName, ["F-FH+", "F-FH-"],...
    'females only: FH+ vs FH-',p,bhfdr,colorLim,clusterOrder,1,...
    'FH+ vs. FH- T-Statistic');
exportgraphics(f, [figureDir '/Fig4jb_femalesVATpairwise_k' num2str(numClusters) '.png']);

% DMN
n = 7;
[~,p,~,~] = ttest2(squeeze(E_network(subjectInfo.sexFHSUD=='M-FH+',:,n)),...
    squeeze(E_network(subjectInfo.sexFHSUD=='M-FH-' ,:,n))); 
f = plotMatricesDiff(squeeze(E_network(subjectInfo.sex=='M' ,:,n)),...
    subjectInfo.sexFHSUD(subjectInfo.sex=='M'),...
    numClusters,clusterNames,saveName, ["M-FH+", "M-FH-"],...
    'males only: FH+ vs FH-',p,bhfdr,colorLim, clusterOrder,1,...
    'FH+ vs. FH- T-Statistic');
exportgraphics(f, [figureDir '/Fig4ka_malesDMNpairwise_k' num2str(numClusters) '.png']);

[~,p,~,~] = ttest2(squeeze(E_network(subjectInfo.sexFHSUD=='F-FH+',:,n)),...
    squeeze(E_network(subjectInfo.sexFHSUD=='F-FH-',:,n))); 
f = plotMatricesDiff(squeeze(E_network(subjectInfo.sex=='F',:,n)),...
    subjectInfo.sexFHSUD(subjectInfo.sex=='F'),...
    numClusters,clusterNames,saveName, ["F-FH+", "F-FH-"],...
    'females only: FH+ vs FH-',p,bhfdr,colorLim, clusterOrder,1,...
    'FH+ vs. FH- T-Statistic');
exportgraphics(f, [figureDir '/Fig4kb_femalesDMNpairwise_k' num2str(numClusters) '.png']);
