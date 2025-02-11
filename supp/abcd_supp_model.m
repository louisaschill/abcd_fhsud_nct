% Supplemental figure: by MRI model
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
models = categorical(cellstr({"Prisma_fit", "Prisma", "DISCOVERY MR750"}));
for i = 1:3
    model_inds = subjectInfo.model == models(i);

    %% ANCOVA
    inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & model_inds;
    anovaVarNames = {'Sex','Age','FHSUD','FD','Income Category','Site' ...
        'ParentEd','Race','In Utero Substances'};
    interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
    contVar = [2,4];
    anovaVars = {subjectInfo.sex(inds),subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
        subjectInfo.FD_mean(inds),subjectInfo.income_cat(inds),subjectInfo.site(inds),...
        subjectInfo.parentEd_cat(inds),subjectInfo.race(inds), ...
        subjectInfo.subDuringPreg(inds)};
    [all_results,stats] = run_ANOVA(globalTE(inds),anovaVars,contVar,...
        anovaVarNames,interaction_term_sets,'off');
    disp(models(i))
    disp(all_results)
    [F_stats_sexsud, p_val_sexsud, ~, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
    [F_stats_sud, p_val_sud, ~, ~, ~] = extract_anova_results(all_results, 'FHSUD');

    inds = subjectInfo.sex == 'F';
    [r1,p1] = corr(subjectInfo.familydensitySUD(inds),globalTE(inds),'type','Spearman');
    inds = subjectInfo.sex == 'M';
    [r2,p2] = corr(subjectInfo.familydensitySUD(inds),globalTE(inds),'type','Spearman');
    pcorr = mafdr([p1, p2],'bhfdr','true');
end

%% By MRI model: Global TE 
models = categorical(cellstr({"Prisma_fit", "Prisma", "DISCOVERY MR750"}));

clear g; f= figure('Position', [100,100,1200,500]);
for i = 1:3

    model_inds = subjectInfo.model == models(i);

    g(1,i) = gramm('x', subjectInfo.sexFHSUD, 'y', globalTE, 'color', subjectInfo.FHSUD,...
        'subset', (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') ...
        & model_inds);
    g(1,i).stat_boxplot('width',0.15);
    g(1,i).stat_violin('normalization','width','dodge', 0,'fill','edge');
    g(1,i).axe_property('YGrid', 'on', 'GridColor', [0.5,0.5,0.5],'YLim', [min(globalTE), mean(globalTE) + 0.00035]);
    g(1,i).set_names('x', '', 'y', 'global transition energy', 'color', '');
    g(1,i).set_text_options('base_size', 18,'Font','calibri','title_scaling', 1);
    g(1,i).no_legend();
    g(1,i).set_color_options('map',[[0.200,0.627,0.173];[0.122,0.471,0.706]]);
    g(1,i).set_title(models(i));
    g.draw(); 
end
exportgraphics(f, [figDir '/suppFig_bymodel_globalTE_k' num2str(numClusters) '.png']);

%% By MRI model: Network TE 
networks = categorical({'vis','som','dat','vat','lim','fpn','dmn','sub','cer'});
clear g; f= figure('Position', [100,100,1200,1200]);
for i = 1:3

    model_inds = subjectInfo.model == models(i);
    % Variables for ANCOVA model
    inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & model_inds;
    anovaVarNames = {'Sex','Age','FHSUD','FD','Income Category',...
        'ParentEd','Race','In Utero Substances','Parental History Mental Health'};
    interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
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
    g(1,i) = gramm('x', categorical(networks), 'y', F_stats, 'color', sig_sexsud);
    g(1,i).set_order_options('x',{'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
    g(1,i).coord_flip();
    g(1,i).geom_bar();
    g(1,i).set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
    g(1,i).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
    g(1,i).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],'YLim', [-0.1, 10.5]);
    g(1,i).set_color_options('map',[[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706];]);
    g(1,i).set_names('x','','y', 'f-statistic: sex * family history of SUD', 'color', 'Significance');
    g(1,i).no_legend();

    % PLOT SEX*SUD F STATS
    inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-') & model_inds;
    g(2,i) = gramm('x',subjectInfo.sexFHSUD,'subset',inds);
    g(2,i).stat_bin();
    g(2,i).set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
    g(2,i).set_names('x', '','y', '# of subject', 'color', '','column', '');
    g(2,i).set_names('x','','y', '', 'color', '');
    g(2,i).axe_property('YLim', [0, 250]);
    g(2,i).no_legend();
    g(2,i).set_color_options('map','brewer_dark');

end
g.draw();
exportgraphics(f, [figDir '/suppFig_bymodel_networkTE_k' num2str(numClusters) '.png']);

%% demographics
% figure; 
% g = gramm('x', subjectInfo.sexFHSUD,'subset', ...
%     (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'));
% g.stat_bin();
% g.facet_grid([], subjectInfo.model);
% g.set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
% g.geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
% g.set_names('x', '','y', 'count', 'color', 'Significance','column', '');
% g.set_color_options('map','brewer_dark');
% g.draw();
% 
% figure;
% g = gramm('x', subjectInfo.income_cat,'y', subjectInfo.age, 'subset', ...
%     (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'));
% g.stat_bin();
% g.facet_grid([], subjectInfo.model);
% g.set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
% g.geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
% g.set_names('x', '','y', 'income category (1 = low, 3 = high)', 'color', 'Significance','column', '');
% g.set_color_options('map','brewer_dark');
% g.draw();
% 
% figure;
% g = gramm('x', subjectInfo.age,'y', subjectInfo.age, 'subset', ...
%     (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'));
% g.stat_bin();
% g.facet_grid([], subjectInfo.model);
% g.set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
% g.geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
% g.set_names('x', '','y', 'age', 'color', 'Significance','column', '');
% g.set_color_options('map','brewer_dark');
% g.draw();
% 
% [h,p,c,s]= ttest2(subjectInfo.age(subjectInfo.model=='DISCOVERY MR750'),subjectInfo.age(subjectInfo.model=='Prisma'))
% [h,p,c,s]= ttest2(subjectInfo.age(subjectInfo.model=='DISCOVERY MR750'),subjectInfo.age(subjectInfo.model=='Prisma_fit'))
% [h,p,c,s]= ttest2(subjectInfo.age(subjectInfo.model=='Prisma'),subjectInfo.age(subjectInfo.model=='Prisma_fit'))
% 
% figure;
% g = gramm('x', subjectInfo.model,'y', subjectInfo.FD_mean, 'subset', ...
%     (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'));
% g.stat_boxplot('width',0.15);
% g.stat_violin('normalization','width','dodge', 0,'fill','edge');
% g.set_text_options('base_size', 18,'font','calibri','title_scaling', 1);
% g.set_names('x', '','y', 'Framewise Displacement', 'color', 'Significance','column', '');
% g.set_color_options('map','brewer_dark');
% g.draw();
% 
% [h,p,c,s]= ttest2(subjectInfo.FD_mean(subjectInfo.model=='DISCOVERY MR750'),subjectInfo.FD_mean(subjectInfo.model=='Prisma'))
% [h,p,c,s]= ttest2(subjectInfo.FD_mean(subjectInfo.model=='DISCOVERY MR750'),subjectInfo.FD_mean(subjectInfo.model=='Prisma_fit'))
% [h,p,c,s]= ttest2(subjectInfo.FD_mean(subjectInfo.model=='Prisma'),subjectInfo.FD_mean(subjectInfo.model=='Prisma_fit'))
% 
