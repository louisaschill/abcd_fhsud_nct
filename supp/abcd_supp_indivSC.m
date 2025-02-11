% Supplemental materials: individual SC 
clear all; close all;
parc = 'fs68'; numClusters = 4;
c = 1; T = 1.501; 
ts_type = 'bp_gsr_gmnorm_exclout';
sc_type = 'indivnoselfSC'; 

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures'; 

% Load subject info
load([dataDir '/subjectInfo_SUDcohort.mat']);

% Load subject energy
nsubjs = 1682; % number of subjs with SC 
cd(resultsDir); 
load(fullfile(resultsDir, ['subjenergy_k', num2str(numClusters), '_T', num2str(T),...
    '_c',num2str(c),'_', ts_type, '_', sc_type, '_' parc , '_' num2str(nsubjs),'.mat' ]));

subjectInfo = subjectInfo(ismember(subjectInfo.subjectkey, subjkeys),:); 
[~, i] =  ismember(subjectInfo.subjectkey, subjkeys); 
subjectInfo = subjectInfo(i,:);

% Check subject order 
assert(all(strcmp(subjkeys, subjectInfo.subjectkey)), 'Subject keys do not match.');

% Network names and assignments 
networks = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN'};
YeoLUT = readtable('fs86_to_yeo.csv'); YeoLUT = table2array(YeoLUT);
numNetworks= 7; YeoLUT = YeoLUT(19:86); 

% Network TE calc 
E_network = nan(length(subjkeys), length(clusterNames)^2, numNetworks);
for i = 1:numNetworks
    network = E_region(:,:,YeoLUT == i); E_network(:, :,i) = sum(network,3);
end
networkTE = squeeze(mean(E_network,2));

% Load FD
load(fullfile(dataDir, 'MeanFD_PerSubject.mat'), 'resultsTable');
resultsTable = resultsTable(ismember(resultsTable.subjectkey,...
    subjectInfo.subjectkey),:); 
[~, i] = ismember(resultsTable.subjectkey, subjectInfo.subjectkey); 
resultsTable = resultsTable(i,:);

% Check subjectInfo and TE results in same order 
if ~all(strcmp(resultsTable.subjectkey,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD mean');
end 
subjectInfo.FD_mean = resultsTable.FD_mean; 
subjectInfo.FD_mean_nanoutlier = resultsTable.FD_mean_nanoutlier; 

% Exclude outliers 
outlier = isoutlier(globalTE,"median"); % Remove outliers
subjectInfo = subjectInfo(~outlier,:);
globalTE = globalTE(~outlier);
networkTE = networkTE(~outlier,:); E_network = E_network(~outlier,:,:); 

%% ANCOVA 
subjectInfo.model = categorical(cellstr(subjectInfo.model));
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
anovaVarNames = {'Sex','Age','FHSUD','FD','MRI Model','Income Category',...
    'ParentEd','Race','In Utero Substances','Parental History Mental Health'};
interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
contVar = [2,4]; 
anovaVars = {subjectInfo.sex(inds),subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
    subjectInfo.FD_mean(inds),subjectInfo.model(inds),subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds),subjectInfo.race(inds), ...
    subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds)};

[all_results,stats] = run_ANOVA(globalTE(inds),anovaVars,contVar,...
    anovaVarNames,interaction_term_sets,'off');
disp(all_results)
[F_stats_sexsud, p_val_sexsud, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
[F_stats_sud, p_val_sud, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');

inds = subjectInfo.sex == 'F'; 
[r1,p1] = corr(subjectInfo.familydensitySUD(inds),globalTE(inds),...
    'type','Spearman','rows','complete'); 
inds = subjectInfo.sex == 'M'; 
[r2,p2] = corr(subjectInfo.familydensitySUD(inds),globalTE(inds),...
    'type','Spearman','rows','complete'); 
pcorr = mafdr([p1, p2],'bhfdr','true'); 

disp(['Females: r = ' num2str(r1,2) ', p = ' num2str(p1,2)]);
disp(['Males: r = ' num2str(r2,2) ', p = ' num2str(p2,2)]);
%% Global and network figures
clear g; 
f = figure('Position', [0 0 1200 1200]);
g(1,1) = gramm('x', subjectInfo.sexFHSUD, 'y', globalTE, 'color', subjectInfo.FHSUD,...
    'subset', (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'));
g(1,1).stat_boxplot('width',0.15);
g(1,1).stat_violin('normalization','width','dodge', 0,'fill','edge');
g(1,1).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
g(1,1).set_names('x', 'sex * family history of SUD', 'y', 'global transition energy', 'color', '');
g(1,1).set_text_options('base_size', 20,'Font','calibri');
g(1,1).no_legend();
g(1,1).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);

g(1,2) = gramm('x', subjectInfo.familydensitySUD, 'y', globalTE,...
    'color', subjectInfo.sex);
g(1,2).geom_point();
g(1,2).stat_glm();
g(1,2).axe_property('YGrid', 'on','GridColor',[0.5,0.5,0.5]);
g(1,2).set_text_options('base_size',20,'Font','calibri');
g(1,2).no_legend();
g(1,2).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]); 
g(1,2).set_names('x', 'family history density', 'y', 'global transition energy','color', 'Sex');

networks = categorical({'vis','som','dat','vat','lim','fpn','dmn'});
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'); 
[all_results] = run_ANOVA(networkTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[F_stats, p_val, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');
p_val_corr = mafdr(p_val, 'bhfdr', true);
sig_sud= makeSigCategory(p_val, p_val_corr); 

% PLOT SUD F STATS
g(2,1) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sud);
g(2,1).set_order_options('x',{'cer', 'sub', 'dmn', 'fpn', 'lim', 'vat','dat','som','vis'}); 
g(2,1).coord_flip();
g(2,1).geom_bar();
g(2,1).set_text_options('base_size', 25,'font','calibri');
g(2,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(2,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5]);
g(2,1).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706]]);
g(2,1).set_names('x','','y', 'f-statistic: family history of SUD', 'color', 'Significance');

[F_stats, p_val, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
p_val_corr = mafdr(p_val, 'bhfdr', true);
sig_sexsud= makeSigCategory(p_val, p_val_corr); 

% PLOT SEX*SUD F STATS
g(2,2) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sexsud);
g(2,2).set_order_options('x',{'cer','sub','dmn','fpn','lim','vat','dat','som','vis'}); 
g(2,2).coord_flip(); 
g(2,2).geom_bar();
g(2,2).set_text_options('base_size', 25,'font','calibri');
g(2,2).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(2,2).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5]);
g(2,2).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
g(2,2).set_names('x','','y', 'f-statistic: sex * family history of SUD', 'color', 'Significance');
g(2,2).no_legend(); 

g.set_title('Individual SC: DK68')
g.draw();
exportgraphics(f, [figureDir '/suppFig_indivSC_globalnetwork_k' num2str(numClusters) '.png']);

%% by income boxplot
% incomes = {'1', '2', '3'}; 
% clear g; f= figure('Position', [0,0,1200,1200]);
% for i = 1:3
% 
%     income_inds = subjectInfo.income_cat == incomes{i}; 
% 
%     g(1,i) = gramm('x', subjectInfo.sexFHSUD, 'y', globalTE, 'color', subjectInfo.FHSUD,...
%         'subset', (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & income_inds);
%     g(1,i).stat_boxplot('width',0.15);
%     g(1,i).stat_violin('normalization','width','dodge', 0,'fill','edge');
%     g(1,i).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5]);
%     g(1,i).set_names('x', '', 'y', 'global transition energy', 'color', '');
%     g(1,i).set_text_options('base_size', 18,'Font','calibri','title_scaling', 1);
%     g(1,i).no_legend();
%     g(1,i).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);
%     g(1,i).set_title([' Income: ' incomes{i}]);
% 
%     inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & income_inds;
%     anovaVarNames = {'Sex','Age','FHSUD','FD','ParentEd','Race','In Utero Substances','Parental History Mental Health'};
%     interaction_term_sets = {{'Sex', 'FHSUD'}};
%     contVar = [2,4];
%     anovaVars = {subjectInfo.sex(inds),subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
%         subjectInfo.FD_mean(inds),...
%         subjectInfo.parentEd_cat(inds),subjectInfo.race(inds), ...
%         subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds)};
%     [all_results] = run_ANOVA(networkTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
%     [F_stats, p_val, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
%     p_val_corr = mafdr(p_val, 'bhfdr', true);
%     sig_sexsud = makeSigCategory(p_val, p_val_corr);
% 
%     g(2,i) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sexsud);
%     g(2,i).set_order_options('x',{'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
%     g(2,i).coord_flip();
%     g(2,i).geom_bar();
%     g(2,i).set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
%     g(2,i).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5]);
%     g(2,i).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
%     g(2,i).set_names('x','','y', 'f-statistic: sex * family history of SUD', 'color', 'Significance');
%     g(2,i).no_legend();
% 
%     g.draw();
% end
% g.set_title('Individual SC: DK68'); 
% g.draw();
% exportgraphics(f, [figureDir '/suppFig_indivSC_byincome_k' num2str(numClusters) '.png']);
