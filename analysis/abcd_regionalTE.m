%% Louisa Schilling - updated Dec 2024  
% Figure 3: Regional level TE analysis 
clear all; close all; 

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures'; 

parc = 'fs86'; numClusters = 4; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg'; 

% Network names and assignments 
networks = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
YeoLUT = readtable('fs86_to_yeo.csv'); YeoLUT = table2array(YeoLUT);
numNetworks= 9;

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
regionalTE = regionalTE(~outlier,:);  E_region = E_region(~outlier,:,:); 

% Load region names
regionNames = readtable('fs86_yeo7_lobe.txt'); 
regionNames = regionNames.Var2;

% Names of transitions between clusters
transitions = {};c = 0;
for i = 1:numClusters
    for j = 1:numClusters 
        c = c+1; 
        transitions{c} = [clusterNames{j} ' to ' clusterNames{i}];
    end
end
transitions = categorical(transitions); 

%% ANOVA Variable set up
subjectInfo.model = categorical(cellstr(subjectInfo.model));
% Variables for ANCOVA model 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-'); %not FH+/-
anovaVarNames = {'Sex','Age','FHSUD','FD','MRI Model','Income Category',...
    'Parent Edu Category','Race','In Utero Substances',...
    'Parental History Mental Health','Puberty'};
interaction_term_sets = {{'Sex', 'FHSUD'}, {'FHSUD','Income Category'}};
contVar = [2,4]; 
anovaVars = {subjectInfo.sex(inds)=='F', subjectInfo.age(inds),...
    subjectInfo.FHSUD(inds) == 'FH+', subjectInfo.FD_mean(inds),...
    subjectInfo.model(inds), subjectInfo.income_cat(inds),...
    subjectInfo.parentEd_cat(inds), subjectInfo.race(inds),...
    subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds),...
    subjectInfo.puberty_combined(inds)};

[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsud, psud, ~,eta_sq_sud,partial_eta_sq_sud] = extract_anova_results(all_results, 'FHSUD');

[Fsexsud, psexsud, ~,eta_sq_sexsud,partial_eta_sq_sexsud, sigInd, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
psexsudcorr = mafdr(psexsud,'bhfdr',true); 

disp('Regional TE Sex:FHSUD sig :'); disp(regionNames(psexsud < 0.05));
disp('Regional TE Sex:FHSUD sig post-FDR:'); disp(regionNames(psexsudcorr < 0.05));
save('Fsexsud.mat', 'Fsexsud'); 
save('partial_eta_sq_sud.mat', 'partial_eta_sq_sud'); 
save('partial_eta_sq_sexsud', 'partial_eta_sq_sexsud'); 

% FHSUD
[Fsud, psud,~,~,~] = extract_anova_results(all_results, 'FHSUD');
psudcorr = mafdr(psud,'bhfdr',true);
disp('Regional TE FHSUD sig :'); disp(regionNames(psud < 0.05))
disp('Regional TE FHSUD sig post-FDR:'); disp(regionNames(psudcorr < 0.05)); 
save('Fsud.mat', 'Fsud'); 

%% SexSUD t-tests 
sig_ind_corr = find(psexsud <0.05);

pvals = []; r = []; 
for i = 1:length(sig_ind_corr)
    ind = sig_ind_corr(i);
    [all_results,stats] = run_ANOVA(regionalTE(inds,ind),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
    [F_stats, p_val, p_val_corr, ~, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
    disp([regionNames{ind} ' TE sex*FHSUD : F = ', num2str(F_stats,3), ' p = ', num2str(p_val,3), ', pFDR = ' num2str(p_val_corr,3)]);

    [~,p,~,s] = ttest2(regionalTE(subjectInfo.sexFHSUD == 'F-FH+',ind), ...
        regionalTE(subjectInfo.sexFHSUD == 'F-FH-',ind)); t= s.tstat;
    disp([regionNames{ind} ' TE females ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);
    pvals = [pvals p];
    [~,p,~,s] = ttest2(regionalTE(subjectInfo.sexFHSUD == 'M-FH+',ind), ...
        regionalTE(subjectInfo.sexFHSUD == 'M-FH-',ind)); t= s.tstat;
    disp([regionNames{ind} ' TE males ttest : t =  ', num2str(t,3), ' p = ', num2str(p,3)]);
    pvals = [pvals p];

end 
pcorr = mafdr(pvals, 'bhfdr', 'true');

count = 0; 
for i = 1:length(sig_ind_corr)
    ind = sig_ind_corr(i);
    disp(regionNames{ind})
    disp([pcorr(i+ count), pcorr(i+count+1)])
    count = count + 1; 
end 

%% Variables for ANCOVA model 
SUDregions = {'L HIPP', 'L BK', 'L ST', 'R cMF', 'r pOB', 'r ST'};
SUD_order = [1,5,4,3,6,2]; 


% FHSUD 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsud, psud, ~, eta_sq_sexsud, partial_eta_sq_sexsud, sigind,~] = extract_anova_results(all_results, 'FHSUD');
indSUD = find(psud < 0.05); 

[~,pSUD,~,s]=ttest2(regionalTE(subjectInfo.FHSUD == 'FH+' ,...
    indSUD),regionalTE(subjectInfo.FHSUD == 'FH-',indSUD));
tSUD = s.tstat; pSUD_corr = mafdr(pSUD, 'BHFDR',true);
sig_sud= makeSigCategory(pSUD, pSUD_corr);
regions = regionNames(indSUD);

% PLOT SUD t-stats 
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',SUDregions, 'y', tSUD); %,'lightness', sig_sud);
g(1,1).set_order_options('x',SUDregions(SUD_order)); 
g(1,1).coord_flip(); g(1,1).geom_bar('width',0.35);
g(1,1).set_text_options('base_size', 25,'font','calibri');
g(1,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5]);
g(1,1).set_color_options('map',[0.122,0.471,0.706]);
g(1,1).set_names('x','','y', 'FH+ vs FH- t-statistic');
g.draw(); 
exportgraphics(f, ['/Users/louisaschilling/Desktop/ABCD_FHSUD/Figures/Fig3_bar_sud_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% SEX:FHSUD 
SEXSUDregions = categorical({'r. cerebellum', 'l. isthmus cingulate', ...
     'l. pars orbitalis',...
    'l. pericalcarine', 'l. superior parietal', 'r. pars orbitalis',...
    'r. superior parietal', 'r. supramarginal'}); 
SEXSUD_order = [1,4,7,5,8,6,3,2];% flipped to be plotted in top down order

sexSUDregions = {'L IST', 'L lOFC', 'R lOFC', 'L MT',...
    'L pOR', 'R pOR', 'L PCL', 'L SP', 'R CER'};
sexSUD_order = [2, 3, 5, 6, 4, 8, 9, 7, 1]; % flipped to be plotted in top down order 


indSEXSUD = find(psexsud < 0.05); 
tSEXSUD= []; pSEXSUD= []; regions = {}; MF = {};
for i = 1:length(indSEXSUD)
    [~,pSEXSUD_F,~,s]=ttest2(regionalTE(subjectInfo.sexFHSUD == 'F-FH+',indSEXSUD(i)),...
        regionalTE(subjectInfo.sexFHSUD == 'F-FH-',indSEXSUD(i)));
    tSEXSUD_F = s.tstat; pSEXSUD = [pSEXSUD pSEXSUD_F]; tSEXSUD = [tSEXSUD tSEXSUD_F];
    [~,pSEXSUD_M,~,s]=ttest2(regionalTE(subjectInfo.sexFHSUD == 'M-FH+',indSEXSUD(i)),...
        regionalTE(subjectInfo.sexFHSUD == 'M-FH-',indSEXSUD(i)));
    tSEXSUD_M = s.tstat; pSEXSUD = [pSEXSUD pSEXSUD_M]; tSEXSUD = [tSEXSUD tSEXSUD_M];
    regions = [regions SEXSUDregions(i),SEXSUDregions(i)];
    MF = [MF, "F", "M"];
end 

pSEXSUD_corr = mafdr(pSEXSUD, 'BHFDR',true);
sig_sexsud= makeSigCategory(pSEXSUD, pSEXSUD_corr);

x = min(tSEXSUD)-0.6; 
% PLOT SEXSUD t-stats 
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',regions, 'y', tSEXSUD, 'color', categorical(MF))% ,'lightness',sig_sexsud);
g(1,1).set_order_options('x',SEXSUDregions(SEXSUD_order)); 
g(1,1).coord_flip(); g(1,1).geom_bar('dodge',0.5);
g(1,1).set_text_options('base_size', 25,'font','calibri');
g(1,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],...
    'YLim', [x, abs(x)]);
g(1,1).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g(1,1).no_legend(); 
g(1,1).set_names('x','','y', 'FH+ vs FH- t-statistic', 'color', '');
g.draw(); 
exportgraphics(f, ['/Users/louisaschilling/Desktop/ABCD_FHSUD/Figures/Fig3_bar_sexsud_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% FHD Correlations: FHSUD SIG REGIONS 
SUDregions = {'r. amygdala', 'l. banks of sts', 'l. paracentral', 'l. superior temporal',...
    'r. banks of sts', 'r. paracentral', 'r. superior temporal'}; 
SUD_order = [1,7,4,6,3,5,2]; % flipped to be plotted in top down order 

% FHSUD 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsud, psud, ~,sigind,~] = extract_anova_results(all_results, 'FHSUD');
indSUD = find(psud < 0.05); 

[rFHP,pFHP]=corr(regionalTE(:,indSUD),subjectInfo.familydensitySUD,'type','spearman');
pFHP_corr = mafdr(pFHP,'BHFDR',true);
sig_sud= makeSigCategory(pFHP, pFHP_corr);
regions = regionNames(indSUD);

% PLOT SUD t-stats 
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',SUDregions, 'y', rFHP); %,'lightness',sig_sud);
g(1,1).set_order_options('x',SUDregions(SUD_order)); 
g(1,1).coord_flip(); g(1,1).geom_bar('width',0.35);
g(1,1).set_text_options('base_size', 25,'font','calibri');
g(1,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5]);
g(1,1).set_color_options('map',[0.122,0.471,0.706]);
g(1,1).set_names('x','','y', 'r','column','','row','');
g(1,1).no_legend(); 
g.draw(); 
exportgraphics(f, ['/Users/louisaschilling/Desktop/ABCD_FHSUD/Figures/Fig3_bar_sud_FHDcorr_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% FHD correlations: Sex*SUD regions 
SEXSUDregions = categorical({'r. cerebellum', 'l. isthmus cingulate', ...
     'l. pars orbitalis',...
    'l. pericalcarine', 'l. superior Parietal', 'r. pars orbitalis',...
    'r. superior parietal', 'r. supramarginal'}); 
SEXSUD_order = [1,4,7,5,8,6,3,2];% flipped to be plotted in top down order

% SEX:FHSUD 
[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsexsud, psexsud, ~,sigInd,~] = extract_anova_results(all_results, 'Sex:FHSUD');
indSEXSUD = find(psexsud < 0.05); 

rSEXSUD= []; pSEXSUD= []; regions = {}; MF = {};
for i = 1:length(indSEXSUD)
    [rf,pf] =corr(regionalTE(subjectInfo.sex == 'F',indSEXSUD(i)),...
        subjectInfo.familydensitySUD(subjectInfo.sex=='F'),'type','Spearman');
    pSEXSUD = [pSEXSUD pf]; rSEXSUD = [rSEXSUD rf];  
    [rm, pm] =corr(regionalTE(subjectInfo.sex == 'M',indSEXSUD(i)),...
        subjectInfo.familydensitySUD(subjectInfo.sex=='M'),'type','Spearman');
    pSEXSUD = [pSEXSUD pm]; rSEXSUD = [rSEXSUD rm]; 
    regions = [regions SEXSUDregions(i),SEXSUDregions(i)];
    MF = [MF, "F", "M"];
end 

pSEXSUD_corr = mafdr(pSEXSUD, 'BHFDR',true);
sig_sexsud = makeSigCategory(pSEXSUD, pSEXSUD_corr);

x = min(rSEXSUD) - 0.025; 
% PLOT SEXSUD t-stats 
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',regions, 'y', rSEXSUD, 'color', categorical(MF)); %,'lightness',sig_sexsud);
g(1,1).set_order_options('x',SEXSUDregions(SEXSUD_order)); 
g(1,1).coord_flip(); g(1,1).geom_bar('dodge',0.5);
g(1,1).set_text_options('base_size',25,'font','calibri');
g(1,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],...
    'YLim', [x, abs(x)]);
g(1,1).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g(1,1).no_legend(); 
g(1,1).set_names('x','','y', 'r', 'color', '');
g.draw(); 
exportgraphics(f, ['/Users/louisaschilling/Desktop/ABCD_FHSUD/Figures/Fig3_bar_sexsud_FHDcorr_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');


%% FHSUD SIG REGIONS 

[~,psud,~,s] = ttest2(squeeze(E_region(subjectInfo.FHSUD=='FH+',:,indSUD)),...
    squeeze(E_region(subjectInfo.FHSUD=='FH-',:,indSUD))); 
t = squeeze(s.tstat); psud = squeeze(psud); 
colors = plasma(length(indSUD)); 

tSUD= []; pSUD= [];regions = {}; trans=[];  
for i = 1:length(indSUD)
    for j = 1:16
        ind = indSUD(i); 
        [~,p,~,s]=ttest2(squeeze(E_region(subjectInfo.FHSUD=='FH+',j,ind)),...
            squeeze(E_region(subjectInfo.FHSUD=='FH-',j,ind))); t = squeeze(s.tstat); 
        pSUD = [pSEXSUD pSEXSUD_F]; tSUD = [tSUD t];
        regions = [regions SUDregions(i)];
        trans = [trans j]; 
    end
end

% barplots
f = figure('Position', [100 100 1000 1600]); clear g;
g = gramm('x',categorical(regions),'y', tSUD);
g.facet_wrap(trans,'ncols', 4,'column_labels', false); 
g.coord_flip(); g.geom_bar('dodge',0.5);
g.set_text_options('font','calibri');
g.set_order_options('x', SUDregions(SUD_order)); 
g.set_text_options('base_size',18,'font','calibri');
g.geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g.set_names('column','','row','', 'x','','y','','color','')
g.set_color_options('map',[0.122,0.471,0.706]);
g.draw(); 

exportgraphics(f, ['/Users/louisaschilling/Desktop/ABCD_FHSUD/Figures/Fig3_bar_matrix_sud_k ' num2str(numClusters) '.png']);

col1 = sum(tSUD(ismember(trans, [1,5,9,13])));
col2 = sum(tSUD(ismember(trans, [2,6,10,14])));
col3 = sum(tSUD(ismember(trans, [3,7,11,15])));
col4 = sum(tSUD(ismember(trans, [4,8,12,16])));

td = sum(tSUD(ismember(trans, [3, 4, 7,8]))); 
bu = sum(tSUD(ismember(trans, [9,10, 13, 14])));

%% SEXSUD 

[~,pSEXSUD_F,~,s]=ttest2(squeeze(E_region(subjectInfo.sexFHSUD=='F-FH+',:,indSEXSUD(i))),...
squeeze(E_region(subjectInfo.sexFHSUD=='F-FH-',:,indSEXSUD(i)))); tSEXSUD_F = squeeze(s.tstat);

tSEXSUD= []; pSEXSUD= []; MF = {};regions = {}; trans=[];  
for i = 1:length(indSEXSUD)
    for j = 1:16
        ind = indSEXSUD(i); 
        [~,pSEXSUD_F,~,s]=ttest2(squeeze(E_region(subjectInfo.sexFHSUD=='F-FH+',j,ind)),...
            squeeze(E_region(subjectInfo.sexFHSUD=='F-FH-',j,ind))); tSEXSUD_F = squeeze(s.tstat); 
        pSEXSUD = [pSEXSUD pSEXSUD_F]; tSEXSUD = [tSEXSUD tSEXSUD_F];
        [~,pSEXSUD_M,~,s]=ttest2(squeeze(E_region(subjectInfo.sexFHSUD=='M-FH+',j,ind)),...
            squeeze(E_region(subjectInfo.sexFHSUD=='M-FH-',j,ind))); tSEXSUD_M = squeeze(s.tstat); 
        pSEXSUD = [pSEXSUD pSEXSUD_M]; tSEXSUD = [tSEXSUD tSEXSUD_M];
        MF = [MF, "F", "M"];
        regions = [regions SEXSUDregions(i) SEXSUDregions(i)];
        trans = [trans j, j]; 
    end
end

% barplots
f = figure('Position', [100 100 1000 1600]); clear g;
g = gramm('x',categorical(regions),'y', tSEXSUD,'color',categorical(MF));
g.facet_wrap(trans,'ncols', 4,'column_labels', false); 
g.coord_flip(); g.geom_bar('dodge',0.5);
g.set_text_options('font','calibri');
g.set_order_options('x', SEXSUDregions(SEXSUD_order));
%g.axe_property('YLim', [min(tSEXSUD) - 0.2, max(tSEXSUD) + 0.2]);
g.set_text_options('base_size',18,'font','calibri');
g.geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g.set_names('column','','row','', 'x','','y','','color','')
g.no_legend(); 
g.set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g.draw(); 

exportgraphics(f, ['/Users/louisaschilling/Desktop/ABCD_FHSUD/Figures/Fig3_bar_matrix_sexsud_k ' num2str(numClusters) '.png']);

col1f = sum(tSEXSUD(MF== 'F' & ismember(trans, [1,5,9,13])));
col2f = sum(tSEXSUD(MF== 'F' & ismember(trans, [2,6,10,14])));
col3f = sum(tSEXSUD(MF== 'F' & ismember(trans, [3,7,11,15])));
col4f = sum(tSEXSUD(MF== 'F' & ismember(trans, [4,8,12,16])));

col1m = sum(tSEXSUD(MF== 'M' & ismember(trans, [1,5,9,13])));
col2m = sum(tSEXSUD(MF== 'M' & ismember(trans, [2,6,10,14])));
col3m = sum(tSEXSUD(MF== 'M' & ismember(trans, [3,7,11,15])));
col4m = sum(tSEXSUD(MF== 'M' & ismember(trans, [4,8,12,16])));

td_f = sum(tSEXSUD(MF== 'F' & ismember(trans, [3,4,7,8])));
bu_f = sum(tSEXSUD(MF== 'F' & ismember(trans, [9,10,13,14])));

td_m = sum(tSEXSUD(MF== 'M' &ismember(trans, [3,4,7,8])));
bu_m = sum(tSEXSUD(MF== 'M' &ismember(trans, [9,10,13,14])));

%% dopamine 
[h,p,c,s]=ttest2(regionalTE(subjectInfo.sex=='F' ,:), ...
    regionalTE(subjectInfo.sex=='M',:)); 
t= s.tstat; 

[r1,p1]= corr(t',D1_norm,'type','spearman');
[r2,p2]= corr(t',D2_norm,'type','spearman');

figure; subplot(2,1,1); scatter(D1_norm,t); lsline; 
title(['r = ' num2str(r1), ', p = ' num2str(p1)]);
ylabel('t-stat F vs M'); 
subplot(2,1,2); scatter(D2_norm,t);lsline; 
title(['r = ' num2str(r2), ', p = ' num2str(p2)]);
ylabel('t-stat F vs M'); 

%% f-stats 
inds = (subjectInfo.sex == 'F' | subjectInfo.sex == 'M');
anovaVarNames = {'Sex', 'Age','FHSUD','FD','MRI Model','Income Category',...
    'Parent Education','Race','In Utero Substances','Parental History Mental Health'};
anovaVars = {subjectInfo.sex(inds), subjectInfo.age(inds),subjectInfo.FHSUD(inds)== 'FH+',...
    subjectInfo.FD_mean(inds),categorical(subjectInfo.model(inds)),subjectInfo.income_cat(inds), subjectInfo.parentEd_cat(inds),...
    subjectInfo.race(inds), subjectInfo.subDuringPreg(inds), subjectInfo.parentMH(inds)};
contVar = [2,4];
[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsex, psex, ~,~,~] = extract_anova_results(all_results, 'Sex');
psexcorr = mafdr(psex,'bhfdr',true); 

[r1,p1]= corr(Fsex,D1_norm,'type','spearman');
[r2,p2]= corr(Fsex,D2_norm,'type','spearman');

figure; subplot(2,1,1); scatter(D1_norm,Fsex); lsline; 
title(['r = ' num2str(r1), ', p = ' num2str(p1)]);
ylabel('Sex f-stat'); 
subplot(2,1,2); scatter(D2_norm,Fsex);lsline; 
title(['r = ' num2str(r2), ', p = ' num2str(p2)]);
ylabel('Sex f-stat'); 
