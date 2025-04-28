%% Louisa Schilling - updated Dec 2024  
% Figure 3: Regional level TE analysis 
clear all; close all; 

dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures/Final'; 

parc = 'fs86'; numClusters = 4; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg'; 

% Network names and assignments 
networks = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
YeoLUT = readtable('fs86_to_yeo.csv'); YeoLUT = table2array(YeoLUT);
numNetworks= 9;

% Load subject info
load([dataDir '/subjectInfo_SUDcohort.mat']);
nsubjs = height(subjectInfo);

% Load optimal T from abcd_tsweep.mat
c = 1;
load([resultsDir '/optimalT_k',num2str(numClusters),'_c',num2str(c), '_' ts_type, '_sc_' sc_type '_noself_' num2str(nsubjs) '.mat'],'T'); 

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

outlier = isoutlier(globalTE,"quartiles"); % Remove outliers
disp(['N = ' num2str(length(outlier(outlier==1))) ' subjects with outlier globalTE removed']);
subjectInfo = subjectInfo(~outlier,:);
regionalTE = regionalTE(~outlier,:); E_region = E_region(~outlier,:,:);  
globalTE = globalTE(~outlier);

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

% ANCOVA 
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

[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsud, psud, ~,eta_sq_sud,partial_eta_sq_sud] = extract_anova_results(all_results, 'FHSUD');

[Fsexsud, psexsud, ~,eta_sq_sexsud,partial_eta_sq_sexsud, sigInd, ~] = extract_anova_results(all_results, 'Sex:FHSUD');
psexsudcorr = mafdr(psexsud,'bhfdr',true); 

disp('Regional TE Sex:FHSUD sig :'); disp(regionNames(psexsud < 0.05));
disp('Regional TE Sex:FHSUD sig post-FDR:'); disp(regionNames(psexsudcorr < 0.05));
save('/Users/louisaschilling/Desktop/Fsexsud.mat', 'Fsexsud'); 

% FHSUD
[Fsud, psud, ~, ~, ~, indSUD, ~] = extract_anova_results(all_results, 'FHSUD');
psudcorr = mafdr(psud,'bhfdr',true);
disp('Regional TE FHSUD sig :'); disp(regionNames(psud < 0.05))
disp('Regional TE FHSUD sig post-FDR:'); disp(regionNames(psudcorr < 0.05)); 
save('/Users/louisaschilling/Desktop/Fsud.mat', 'Fsud'); 

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
SUDregions = {'rAMYG', 'lBK', 'lPAC', 'lST', 'rPAC', 'rST'};
SUD_Yeo = [8, 7, 2, 2, 2, 2];
SUD_order = [1, 3, 5, 4, 6, 2]; 

% FHSUD 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,...
    anovaVarNames,interaction_term_sets,'off');

[~,pSUD,~,s]=ttest2(regionalTE(subjectInfo.FHSUD == 'FH+' ,...
    indSUD),regionalTE(subjectInfo.FHSUD == 'FH-',indSUD));
tSUD = s.tstat; pSUD_corr = mafdr(pSUD, 'BHFDR',true);
sig_sud= makeSigCategory(pSUD, pSUD_corr);
regions = regionNames(indSUD);

star_labels = cell(size(tSUD));
for i = 1:length(pSUD)
    if pSUD_corr(i) < 0.05
        star_labels{i} = '**';
    elseif pSUD(i) < 0.05
        star_labels{i} = '*';
    else
        star_labels{i} = '';
    end
end

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

t = tSUD(SUD_order);
star = star_labels(SUD_order);
% Add stars to bar tips (corrected for coord_flip)
ax = g(1,1).facet_axes_handles;
for i = 1:length(tSUD)
    if ~isempty(star{i})
        text(ax, i, t(i) + sign(t(i))*0.0005, star{i},...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
            'FontSize', 30, 'FontWeight', 'bold', 'Color', 'k');
    end
end

exportgraphics(f, [figureDir '/Fig5_tstat_bar_sud_regionalTE_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% SEX:FHSUD 
SEXSUDregions = {'rCER', 'lIST', 'lpOR', 'lSP', 'lSMG', 'rpOR', 'rSP', 'rSMG'};
SEXSUD_YEO = [9, 7, 7, 3, 4, 7, 3, 4];
SEXSUD_order = [1, 4, 7, 5, 8, 2, 3, 6]; % flipped to be plotted in top down order 

SEXSUD_order_2 = [];
for i = 1:length(SEXSUD_order)
    n = SEXSUD_order(i) + (SEXSUD_order(i)-1); 
    SEXSUD_order_2 = [SEXSUD_order_2, n];
    SEXSUD_order_2 = [SEXSUD_order_2, (SEXSUD_order(i)*2)];
end 

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


% Identify significant regions with stars
star_labels = cell(size(pSEXSUD));
for i = 1:length(pSEXSUD)
    if pSEXSUD_corr(i) < 0.05
        star_labels{i} = '**';
    elseif pSEXSUD(i) < 0.05
        star_labels{i} = '*';
    else
        star_labels{i} = '';
    end
end

x = min(tSEXSUD)-0.6; 
% PLOT SEXSUD t-stats 
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',regions, 'y', tSEXSUD, 'color', categorical(MF));
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

% Correct star positions accounting for dodge and coord_flip
ax = g(1,1).facet_axes_handles;

t = tSEXSUD(SEXSUD_order_2);
star = star_labels(SEXSUD_order_2);

count = 0; 
for i = 1:8
    for s = 1:2
        count = count + 1; 
        if ~isempty(star{count})
            if mod(count,2) == 0
                text(ax, i + 0.1, t(count) - 0.2, star{count},...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
                    'FontSize', 30, 'FontWeight', 'bold', 'Color', 'k');
            else
                text(ax, i - 0.20, t(count) + 0.2, star{count},...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
                    'FontSize', 30, 'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end
end

exportgraphics(f, [figureDir '/Fig5_bar_sexsud_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% FHD Correlations: FHSUD SIG REGIONS 

% FHSUD 
inds = (subjectInfo.FHSUD == 'FH+' | subjectInfo.FHSUD == 'FH-');
[all_results] = run_ANOVA(regionalTE(inds,:),anovaVars,contVar,anovaVarNames,interaction_term_sets,'off');
[Fsud, psud, ~,sigind,~] = extract_anova_results(all_results, 'FHSUD');
indSUD = find(psud < 0.05); 

[rFHP,pFHP]=corr(regionalTE(:,indSUD),subjectInfo.familydensitySUD,'type','spearman');
pFHP_corr = mafdr(pFHP,'BHFDR',true);
sig_sud= makeSigCategory(pFHP, pFHP_corr);
regions = regionNames(indSUD);

star_labels = cell(size(rFHP));
for i = 1:length(pFHP)
    if pFHP_corr(i) < 0.05
        star_labels{i} = '**';
    elseif pFHP(i) < 0.05
        star_labels{i} = '*';
    else
        star_labels{i} = '';
    end
end

% Generate the bar plot with stars
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',SUDregions, 'y', rFHP,'lightness',sig_sud);
g(1,1).set_order_options('x',SUDregions(SUD_order)); 
g(1,1).coord_flip(); 
g(1,1).geom_bar('width',0.35);
g(1,1).set_text_options('base_size', 25,'font','calibri');
g(1,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5]);
g(1,1).set_color_options('map',[0.122,0.471,0.706]);
g(1,1).set_names('x','','y', 'r','column','','row','');
g(1,1).no_legend(); 
g.draw(); 

r = rFHP(SUD_order);
star = star_labels(SUD_order);
% Add stars to bar tips (corrected for coord_flip)
ax = g(1,1).facet_axes_handles;
for i = 1:length(rFHP)
    if ~isempty(star{i})
        disp(r(i))
        text(ax, i, r(i) + sign(r(i))*0.0005, star{i},...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
            'FontSize', 30, 'FontWeight', 'bold', 'Color', 'k');
    end
end

% Save figure
exportgraphics(f, [figureDir '/Fig5_bar_sud_FHDcorr_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% FHD correlations: Sex*SUD regions 

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

x = min(rSEXSUD) - 0.02; 

% Identify significant regions with stars
star_labels = cell(size(rSEXSUD));
for i = 1:length(pSEXSUD)
    if pSEXSUD_corr(i) < 0.05
        star_labels{i} = '**';
    elseif pSEXSUD(i) < 0.05
        star_labels{i} = '*';
    else
        star_labels{i} = '';
    end
end

% PLOT SEXSUD r-stats
clear g; f = figure('Position', [100,100,900,1000]);
g(1,1) = gramm('x',regions, 'y', rSEXSUD, 'color', categorical(MF));
g(1,1).set_order_options('x',SEXSUDregions(SEXSUD_order)); 
g(1,1).coord_flip(); 
g(1,1).geom_bar('dodge',0.5);
g(1,1).set_text_options('base_size',25,'font','calibri');
g(1,1).geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
g(1,1).axe_property('TickDir','out', 'Ygrid','on','GridColor',[0.5 0.5 0.5],...
    'YLim', [x, abs(x)]);
g(1,1).set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g(1,1).no_legend(); 
g(1,1).set_names('x','','y', 'r', 'color', '');
g.draw(); 

% Correct star positions accounting for dodge and coord_flip
ax = g(1,1).facet_axes_handles;

% Adjust positions for dodged bars (Male before Female)
unique_regions = SEXSUDregions(SEXSUD_order);
num_regions = length(unique_regions);
dodge_offset = [0.15, -0.15]; % Adjust dodge offset if necessary (Male before Female)

r = rSEXSUD(SEXSUD_order_2);
star = star_labels(SEXSUD_order_2);

count = 0; 
for i = 1:8
    for s = 1:2
        count = count + 1; 
        if ~isempty(star{count})
            if mod(count,2) == 0
                text(ax, i + 0.08, r(count) - 0.005, star{count},...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
                    'FontSize', 30, 'FontWeight', 'bold', 'Color', 'k');
            else
                text(ax, i - 0.2, r(count) + 0.005, star{count},...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
                    'FontSize', 30, 'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end
end
% Save figure
exportgraphics(f, [figureDir '/Fig5_bar_sexsud_FHDcorr_region_k ' num2str(numClusters) '.png'], 'ContentType', 'vector');

%% FHSUD pairwise regional TE 
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

% f = figure('Position', [100 100 1000 1600]); clear g;
% g = gramm('x',trans,'y', tSUD);
% g.facet_wrap(categorical(regions),'ncols', 4,'column_labels', true); 
% g.coord_flip(); g.geom_bar('dodge',0.5);
% g.set_text_options('font','calibri');
% g.set_text_options('base_size',15,'font','calibri');
% g.geom_abline('intercept', 0, 'slope', 0, 'style', 'k--');
% g.set_names('column','','row','', 'x','','y','','color','')
% g.set_color_options('map',[0.122,0.471,0.706]);
% g.draw(); 

exportgraphics(f, [figureDir '/Fig5_bar_matrix_sud_k ' num2str(numClusters) '.png']);

col1 = sum(tSUD(ismember(trans, [1,5,9,13])));
col2 = sum(tSUD(ismember(trans, [2,6,10,14])));
col3 = sum(tSUD(ismember(trans, [3,7,11,15])));
col4 = sum(tSUD(ismember(trans, [4,8,12,16])));

td = sum(tSUD(ismember(trans, [3, 4, 7,8]))); 
bu = sum(tSUD(ismember(trans, [9,10, 13, 14])));

%% SEX*FHSUD pairwise regional TE

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
g.set_names('row','', 'x','','y','','color','');
g.no_legend(); 
g.set_color_options('map',[[0.596, 0.306, 0.639];[0.105, 0.620, 0.467]]);
g.draw(); 

exportgraphics(f, [figureDir '/Fig5_bar_matrix_sexsud_k ' num2str(numClusters) '.png']);
 
