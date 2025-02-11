% Supplemental figure: single site 
clear all; close all; 

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures'; 

parc = 'fs86'; numClusters = 4; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg'; 

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
globalTE = globalTE(~outlier,:); 

%%
networks = {'vis','som','dat','vat','lim','fpn','dmn','sub','cer'};

site_inds = subjectInfo.site == 'site16';
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & site_inds;
anovaVarNames = {'Sex','Age','FHSUD','FD','Income Category', ...
    'ParentEd','Race','In Utero Substances'};
interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
contVar = [2,4];
anovaVars = {subjectInfo.sex(inds),subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
    subjectInfo.FD_mean(inds),subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds),subjectInfo.race(inds), ...
    subjectInfo.subDuringPreg(inds)};
[all_results,stats] = run_ANOVA(globalTE(inds),anovaVars,contVar,...
    anovaVarNames,interaction_term_sets,'off');
disp(all_results)
[F_stats_sexsud, p_val_sexsud, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
[F_stats_sud, p_val_sud, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');

clear g; f= figure('Position', [100,100,1500,700]);
g(1,1) = gramm('x', subjectInfo.sexFHSUD, 'y', globalTE, 'color', subjectInfo.FHSUD,...
    'subset', (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & site_inds);
g(1,1).stat_boxplot('width',0.15);
g(1,1).stat_violin('normalization','width','dodge', 0,'fill','edge');
g(1,1).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
g(1,1).set_names('x', '', 'y', 'global transition energy', 'color', '');
g(1,1).set_text_options('base_size', 18,'Font','calibri','title_scaling', 1);
g(1,1).no_legend();
g(1,1).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);

% Variables for ANCOVA model
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & site_inds;
anovaVarNames = {'Sex','Age','FHSUD','FD','Income Category',...
    'ParentEd','Race','In Utero Substances','Parental History Mental Health'};
interaction_term_sets = {{'Sex', 'FHSUD'},{'FHSUD', 'Income Category'}};
contVar = [2,4];
anovaVars = {subjectInfo.sex(inds),subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
    subjectInfo.FD_mean(inds),subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds),subjectInfo.race(inds), ...
    subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds)};
[all_results] = run_ANOVA(networkTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[F_stats, p_val, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
p_val_corr = mafdr(p_val, 'bhfdr', true);
sig_sexsud = makeSigCategory(p_val, p_val_corr);

% PLOT SEX*SUD F STATS
g(1,2) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sexsud);
g(1,2).set_order_options('x',{'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
g(1,2).coord_flip();
g(1,2).geom_bar();
g(1,2).set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
g(1,2).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,2).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],'YLim', [-0.1, max(F_stats)+0.2]);
g(1,2).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
g(1,2).set_names('x','','y', 'f-statistic: sex * family history of SUD', 'color', 'Significance');
g(1,2).no_legend();

g(2,1) = gramm('x', subjectInfo.sexFHSUD,'subset', ...
    (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & ...
    subjectInfo.site == 'site16');
g(2,1).stat_bin();
g(2,1).set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
g(2,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(2,1).set_names('x', '','y', 'count', 'color', 'Significance','column', '');
g(2,1).set_color_options('map','brewer_dark');

g.set_title(' ');
g.draw();
exportgraphics(f,[figDir '/suppFig_bysite_globalnetworkTE_k' num2str(numClusters) '.png']);

