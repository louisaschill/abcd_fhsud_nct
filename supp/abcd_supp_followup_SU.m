%% FOLLOW UP SUBSTANCE USE
clear all; close all; clc;

%% Set-up 
dataDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
resultsDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Results/energy';
figureDir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Figures/Final'; 
basedir = '/Users/louisaschilling/Desktop/Datasets/ABCD/Data/non_imaging/Release 5.1/core/';

% Set parameters
parc = 'fs86'; numClusters = 4; k = 4; 
ts_type = 'bp_gsr_gmnorm_exclout'; sc_type = 'avg';

% Load subject info
load([dataDir '/subjectInfo_SUDcohort.mat']);
nsubjs = height(subjectInfo);

% Load optimal T value
c = 1;
load([resultsDir '/optimalT_k',num2str(numClusters),'_c',num2str(c), '_' ts_type, '_sc_' sc_type '_noself_' num2str(nsubjs) '.mat'],'T'); 

%% Load data 
% Subject energy
load(fullfile(resultsDir, ['subjenergy_k', num2str(numClusters), '_T', num2str(T),...
    '_c', num2str(c),'_', ts_type, '_', sc_type, 'noselfSC_', parc, '_', num2str(nsubjs) '.mat']));

% Confirm subject key alignment
if ~all(strcmp(subjkeys,subjectInfo.subjectkey))
    error('Subject keys not in same order for energy results');
end 

% Load network info
networks = {'VIS','SOM','DAT','VAT','LIM','FPN','DMN','SUB','CER'};
YeoLUT = readtable('fs86_to_yeo.csv'); YeoLUT = table2array(YeoLUT);
numNetworks = length(networks);

% Calculate network-level TE
E_network = nan(length(subjkeys), length(clusterNames)^2, numNetworks);
for i = 1:numNetworks
    network = E_region(:,:,YeoLUT == i); 
    E_network(:,:,i) = sum(network,3);
end
networkTE = squeeze(mean(E_network,2));

% Load framewise displacement (FD)
load(fullfile(dataDir, 'MeanFD_PerSubject.mat'), 'resultsTable');
if ~all(strcmp(resultsTable.subjectkey,subjectInfo.subjectkey))
    error('Subject keys not in same order for FD results');
end 
subjectInfo.FD_mean = resultsTable.FD_mean;
subjectInfo.FD_mean_nanoutlier = resultsTable.FD_mean_nanoutlier;

% Add energy metrics to subjectInfo
subjectInfo.globalTE = globalTE;
subjectInfo.networkTE = networkTE;
subjectInfo.regionalTE = regionalTE;

% Remove outliers
outlier = isoutlier(subjectInfo.globalTE,"quartiles");
subjectInfo = subjectInfo(~outlier,:);
disp(['N = ' num2str(sum(outlier)) ' outlier subjects removed.']);

%% Load Substance Use (SU) data
%% --- Get Follow-up Substance Use (Years 1-3 only) ---

basedir = '/Users/louisaschilling/Desktop/Datasets/ABCD/Data/non_imaging/Release 5.1/core/';
cd(basedir);

subjs = subjectInfo.subjectkey;
SU_all = readtable([basedir '/substance-use/su_y_sui.csv']);

% Define substance use variables
su_var = {'src_subject_id', 'isip_1b_yn_l', 'isip_2_l', 'tlfb_tob_puff_l', 'tlfb_chew_use_l', ...
    'tlfb_cigar_use_l', 'tlfb_hookah_use_l', 'tlfb_pipes_use_l', 'tlfb_nicotine_use_l', ...
    'tlfb_mj_puff_l', 'tlfb_blunt_use_l', 'tlfb_mj_conc_use_l', 'tlfb_mj_drink_use_l', ...
    'tlfb_tincture_use_l', 'tlfb_mj_synth_use_l', 'tlfb_coc_use_l', 'tlfb_bsalts_use_l', ...
    'tlfb_meth_use_l', 'tlfb_mdma_use_l', 'tlfb_ket_use_l', 'tlfb_ghb_use_l', 'tlfb_opi_use_l', ...
    'tlfb_lsd_use_l', 'tlfb_shrooms_use_l', 'tlfb_salvia_use_l', 'tlfb_steroids_use_l', ...
    'tlfb_bitta_use_l', 'tlfb_inhalant_use_l', 'tlfb_amp_use_l', 'tlfb_tranq_use_l', ...
    'tlfb_vicodin_use_l', 'tlfb_cough_use_l', 'tlfb_other_use_l'};
su_var_names = su_var;
su_var_names{1} = 'subjectkey';

% Initialize
events = {'1_year_follow_up_y_arm_1', '2_year_follow_up_y_arm_1', '3_year_follow_up_y_arm_1'};
years = {'1','2','3'};
su_years = cell(1,3);

% Loop through 1,2,3 year follow-ups
for i = 1:3
    event = events{i};
    year = years{i};
    SU = extractInstrument(SU_all, su_var, su_var_names, event, subjs);

    % Create year-specific table
    T = table;
    T.subjectkey = SU.subjectkey;
    T.(['alcohol_' year 'yr']) = double(SU.isip_1b_yn_l == 1 | SU.isip_2_l == 1);
    T.(['nicotine_' year 'yr']) = (SU.tlfb_tob_puff_l == 1) + (SU.tlfb_chew_use_l == 1) + ...
                                  (SU.tlfb_hookah_use_l == 1) + (SU.tlfb_nicotine_use_l == 1) + ...
                                  (SU.tlfb_cigar_use_l == 1) + (SU.tlfb_pipes_use_l == 1);
    T.(['cannabis_' year 'yr']) = (SU.tlfb_mj_puff_l == 1) + (SU.tlfb_blunt_use_l == 1) + ...
                                  (SU.tlfb_mj_conc_use_l == 1) + (SU.tlfb_mj_drink_use_l == 1) + ...
                                  (SU.tlfb_tincture_use_l == 1) + (SU.tlfb_mj_synth_use_l == 1);
    T.(['rec_drug_' year 'yr']) = (SU.tlfb_coc_use_l == 1) + (SU.tlfb_bsalts_use_l == 1) + ...
                                  (SU.tlfb_meth_use_l == 1) + (SU.tlfb_mdma_use_l == 1) + ...
                                  (SU.tlfb_ket_use_l == 1) + (SU.tlfb_ghb_use_l == 1) + (SU.tlfb_opi_use_l == 1);
    T.(['hallucinogen_' year 'yr']) = (SU.tlfb_lsd_use_l == 1) + (SU.tlfb_shrooms_use_l == 1) + ...
                                      (SU.tlfb_salvia_use_l == 1);
    T.(['prescription_' year 'yr']) = (SU.tlfb_amp_use_l == 1) + (SU.tlfb_tranq_use_l == 1) + ...
                                      (SU.tlfb_vicodin_use_l == 1) + (SU.tlfb_cough_use_l == 1) + ...
                                      (SU.tlfb_steroids_use_l == 1);
    T.(['other_drug_' year 'yr']) = (SU.tlfb_other_use_l == 1) + (SU.tlfb_bitta_use_l == 1) + ...
                                    (SU.tlfb_inhalant_use_l == 1);

    % Overall SU at year
    T.(['SU_' year 'year']) = nan(height(T),1);
    T.(['SU_' year 'year'])(T.(['alcohol_' year 'yr']) > 0 | T.(['nicotine_' year 'yr']) > 0 | ...
                            T.(['cannabis_' year 'yr']) > 0 | T.(['rec_drug_' year 'yr']) > 0 | ...
                            T.(['hallucinogen_' year 'yr']) > 0 | T.(['prescription_' year 'yr']) > 0 | ...
                            T.(['other_drug_' year 'yr']) > 0) = 1;
    T.(['SU_' year 'year'])(T.(['alcohol_' year 'yr']) == 0 & T.(['nicotine_' year 'yr']) == 0 & ...
                            T.(['cannabis_' year 'yr']) == 0 & T.(['rec_drug_' year 'yr']) == 0 & ...
                            T.(['hallucinogen_' year 'yr']) == 0 & T.(['prescription_' year 'yr']) == 0 & ...
                            T.(['other_drug_' year 'yr']) == 0) = 0;
    
    su_years{i} = T;
end

%% --- Merge across all years 1-3

tabs = su_years;
SU = mergeTables(tabs);

% Summarize use across all years
SU.alcohol = ((SU.alcohol_1yr == 1) + (SU.alcohol_2yr == 1) + (SU.alcohol_3yr == 1)) >= 1;
SU.nicotine = ((SU.nicotine_1yr == 1) + (SU.nicotine_2yr == 1) + (SU.nicotine_3yr == 1)) >= 1;
SU.cannabis = ((SU.cannabis_1yr == 1) + (SU.cannabis_2yr == 1) + (SU.cannabis_3yr == 1)) >= 1;
SU.rec_drug = ((SU.rec_drug_1yr == 1) + (SU.rec_drug_2yr == 1) + (SU.rec_drug_3yr == 1)) >= 1;
SU.hallucinogen = ((SU.hallucinogen_1yr == 1) + (SU.hallucinogen_2yr == 1) + (SU.hallucinogen_3yr == 1)) >= 1;
SU.prescription = ((SU.prescription_1yr == 1) + (SU.prescription_2yr == 1) + (SU.prescription_3yr == 1)) >= 1;
SU.other_drug = ((SU.other_drug_1yr == 1) + (SU.other_drug_2yr == 1) + (SU.other_drug_3yr == 1)) >= 1;

% Overall SU
SU.SU = SU.alcohol | SU.nicotine | SU.cannabis | SU.rec_drug | SU.hallucinogen | SU.prescription | SU.other_drug;

% SU_ANY across 1/2/3 year
SU.SU_ANY = nan(height(SU),1);
SU.SU_ANY(SU.SU_1year == 1 | SU.SU_2year == 1 | SU.SU_3year == 1) = 1;
SU.SU_ANY(SU.SU_1year == 0 & SU.SU_2year == 0 & SU.SU_3year == 0) = 0;

% Keep only subjects with SU_ANY info
SU = SU(~isnan(SU.SU_ANY),:);

% --- Merge into subjectInfo

subjectInfo = subjectInfo(ismember(subjectInfo.subjectkey, SU.subjectkey),:);
subjectInfo = mergeTables({subjectInfo, SU});
subjectInfo.sexSU = categorical(cellstr(strcat(char(subjectInfo.sex), '-', num2str(subjectInfo.SU))));

%% --- Plot 1: Histogram of future SU by sex and FHSUD ---

f = figure;
g = gramm('x', categorical(subjectInfo.sexSU), 'color', subjectInfo.FHSUD, 'subset', subjectInfo.FHSUD ~= 'FH-');
g.stat_bin();
g.set_names('x', 'Future Substance Use', 'color', 'FH of SUD');
g.set_text_options('base_size', 20, 'Font', 'calibri');
g.draw();
exportgraphics(f, fullfile(figureDir, 'supp_followup_SU_histogram.png'));

%% --- Make substance counts summary table ---

substance_counts = sum([subjectInfo.SU > 0, ...
    subjectInfo.alcohol > 0, ...
    subjectInfo.nicotine > 0, ...
    subjectInfo.cannabis > 0, ...
    subjectInfo.rec_drug > 0, ...
    subjectInfo.hallucinogen > 0, ...
    subjectInfo.prescription > 0, ...
    subjectInfo.other_drug > 0, ...
    subjectInfo.SU == 0], 1);

substance_names = {'Any Substances', 'Alcohol', 'Nicotine', 'Cannabis', ...
    'Recreational', 'Hallucinogen', 'Prescription', 'Other', 'None'};

total_subjects = height(subjectInfo);
substance_percentages = (substance_counts / total_subjects) * 100;

substance_table = table(substance_names', substance_counts', substance_percentages', ...
    'VariableNames', {'Substance', 'Number_of_Subjects', 'Percentage_of_Subjects'});
disp('Summary of substance use:');
disp(substance_table);

%% --- Group counts by sexFHSUD ---

group_names = unique(subjectInfo.sexFHSUD);
grouped_counts = zeros(length(group_names), length(substance_names));
group_sizes = zeros(length(group_names), 1);

for g = 1:length(group_names)
    idx = subjectInfo.sexFHSUD == group_names(g);
    group_sizes(g) = sum(idx);

    grouped_counts(g, :) = sum([subjectInfo.SU(idx) > 0, ...
        subjectInfo.alcohol(idx) > 0, ...
        subjectInfo.nicotine(idx) > 0, ...
        subjectInfo.cannabis(idx) > 0, ...
        subjectInfo.rec_drug(idx) > 0, ...
        subjectInfo.hallucinogen(idx) > 0, ...
        subjectInfo.prescription(idx) > 0, ...
        subjectInfo.other_drug(idx) > 0, ...
        subjectInfo.SU(idx) == 0], 1);
end

% Display grouped counts
disp('Grouped counts by sexFHSUD:');
grouped_table = array2table(grouped_counts, 'VariableNames', substance_names, 'RowNames', cellstr(group_names));
disp(grouped_table);

% Grouped percentages
grouped_percentages = (grouped_counts ./ group_sizes) * 100;
disp('Grouped percentages by sexFHSUD:');
grouped_table_percent = array2table(grouped_percentages, 'VariableNames', substance_names, 'RowNames', cellstr(group_names));
disp(grouped_table_percent);

%% --- Plot 2: Violin plots of Global TE by sex and future SU ---

f = figure('Position', [100, 100, 800, 1000]);
g = gramm('x', categorical(subjectInfo.sexSU), 'y', subjectInfo.globalTE, ...
    'color', subjectInfo.sex, 'subset', subjectInfo.FHSUD ~= 'FH-');
g.stat_boxplot('width', 0.15);
g.stat_violin('normalization', 'width', 'dodge', 0, 'fill', 'edge');
g.set_names('x', 'Follow-up SU', 'y', 'Global TE', 'color', 'Sex');
g.set_color_options('map', [[0.596, 0.306, 0.639]; [0.105, 0.620, 0.467]]);
g.set_text_options('base_size', 20, 'Font', 'calibri');
g.axe_property('YGrid', 'on', 'GridColor', [0.5, 0.5, 0.5]);
g.draw();
exportgraphics(f, fullfile(figureDir, 'all_followup_SU_globalTE.png'));

%% --- t-tests Global TE by group (Females, Males, All) ---

inds = subjectInfo.FHSUD ~= 'FH-'; % exclude FH- if needed

groups = {
    'Females', subjectInfo.sexSU == 'F-1' & inds, subjectInfo.sexSU == 'F-0' & inds;
    'Males',   subjectInfo.sexSU == 'M-1' & inds, subjectInfo.sexSU == 'M-0' & inds;
    'All',     subjectInfo.SU == 1 & inds,        subjectInfo.SU == 0 & inds
};

results = cell(size(groups,1), 5); % Group, t, p, d, p_fdr

raw_p = zeros(size(groups,1),1);

for i = 1:size(groups,1)
    label = groups{i,1};
    g1 = subjectInfo.globalTE(groups{i,2});
    g2 = subjectInfo.globalTE(groups{i,3});
    
    [~, p, ~, s] = ttest2(g1, g2);
    tstat = s.tstat;
    
    % Cohen's d
    n1 = length(g1); n2 = length(g2);
    spooled = sqrt(((n1-1)*var(g1) + (n2-1)*var(g2)) / (n1+n2-2));
    d = (mean(g1) - mean(g2)) / spooled;
    
    results{i,1} = label;
    results{i,2} = tstat;
    results{i,3} = p;
    results{i,4} = d;
    
    raw_p(i) = p;
end

% FDR correct
p_fdr = mafdr(raw_p, 'BHFDR', true);
for i = 1:length(p_fdr)
    results{i,5} = p_fdr(i);
end

% Display
summary = cell2table(results, 'VariableNames', {'Group', 't', 'p', 'Cohen_d', 'p_fdr'});
disp('T-tests Global TE results:');
disp(summary);

%% --- ANOVA on global TE and network TE ---

% Define indices for subjects (excluding FH-)
inds = (subjectInfo.FHSUD ~= 'FH-');

% Define ANOVA variables
anovaVarNames = {'Sex','Age','SUB','FD','Model','Income','Parent_Ed','Race','Prenatal','Parent_MH','Puberty'};
interaction_term_sets = {{'Sex','SUB'},{'Sex','Puberty'},{'SUB','Income'}};
contVar = [2,4]; % Age and FD

anovaVars = {subjectInfo.sex(inds) == 'F', ...
             subjectInfo.age(inds), ...
             subjectInfo.SU(inds), ...
             subjectInfo.FD_mean(inds), ...
             subjectInfo.model(inds), ...
             subjectInfo.income_cat(inds), ...
             subjectInfo.parentEd_cat(inds), ...
             subjectInfo.race(inds), ...
             categorical(subjectInfo.subDuringPregAfter(inds)), ...
             subjectInfo.parentMH(inds) == 1, ...
             categorical(subjectInfo.pds_mod(inds))};

% --- Global TE ANOVA
disp('ANOVA on Global TE:');
globalTE_results = run_ANOVA(subjectInfo.globalTE(inds), anovaVars, contVar, anovaVarNames, interaction_term_sets, 'off');
disp(globalTE_results);

% --- Network TE ANOVA
disp('ANOVA on Network TE:');
networkTE_results = run_ANOVA(subjectInfo.networkTE(inds,:), anovaVars, contVar, anovaVarNames, interaction_term_sets, 'off');

% Extract SU and Sex:SU results
[F_stats_SU, p_val_SU, ~,~,~, sig_su, ~] = extract_anova_results(networkTE_results, 'SUB');
[F_stats_sexSU, p_val_sexSU, ~,~,~, sig_sexsu, ~] = extract_anova_results(networkTE_results, 'Sex:SUB');

% Correct for multiple comparisons
p_val_corr_SU = mafdr(p_val_SU, 'BHFDR', true);
p_val_corr_sexSU = mafdr(p_val_sexSU, 'BHFDR', true);

% Define significance categories
sig_su_cat = makeSigCategory(p_val_SU, p_val_corr_SU);
sig_sexsu_cat = makeSigCategory(p_val_sexSU, p_val_corr_sexSU);

networks_cat = categorical({'vis','som','dat','vat','lim','fpn','dmn','sub','cer'});

%% --- Plot 3: Network TE F-stats (future SU and Sex×SU)

clear g;
f = figure('Position', [100,100,1500,700]);

% Future SU
g(1,1) = gramm('x', networks_cat, 'y', F_stats_SU, 'color', sig_su_cat);
g(1,1).set_order_options('x', {'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
g(1,1).coord_flip();
g(1,1).geom_bar();
g(1,1).set_text_options('base_size', 25, 'Font', 'calibri');
g(1,1).axe_property('TickDir','out','Ygrid','on','GridColor',[0.5,0.5,0.5], 'YLim', [-0.05, max([F_stats_SU;F_stats_sexSU])+0.2]);
g(1,1).set_color_options('map', [[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706]]);
g(1,1).set_names('x','','y','F-statistic (Future SU)', 'color', 'Significance');

% Sex×Future SU
g(1,2) = gramm('x', networks_cat, 'y', F_stats_sexSU, 'color', sig_sexsu_cat);
g(1,2).set_order_options('x', {'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
g(1,2).coord_flip();
g(1,2).geom_bar();
g(1,2).set_text_options('base_size', 25, 'Font', 'calibri');
g(1,2).axe_property('TickDir','out','Ygrid','on','GridColor',[0.5,0.5,0.5], 'YLim', [-0.05, max([F_stats_SU;F_stats_sexSU])+0.2]);
g(1,2).set_color_options('map', [[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706]]);
g(1,2).set_names('x','','y','F-statistic (Sex × Future SU)', 'color', 'Significance');
g(1,2).no_legend();

g.draw();
exportgraphics(f, fullfile(figureDir, ['supp_network_followup_k' num2str(numClusters) '.png']));

%% --- t-tests Network TE separately for females and males

% Initialize
t_stats_females = zeros(1,length(networks_cat));
t_stats_males = zeros(1,length(networks_cat));
p_values_females = zeros(1,length(networks_cat));
p_values_males = zeros(1,length(networks_cat));

% Loop over networks
for i = 1:length(networks_cat)
    % Female
    [~, pf, ~, s_f] = ttest2(subjectInfo.networkTE(subjectInfo.sexSU == 'F-1' & inds, i), ...
                             subjectInfo.networkTE(subjectInfo.sexSU == 'F-0' & inds, i));
    t_stats_females(i) = s_f.tstat;
    p_values_females(i) = pf;
    
    % Male
    [~, pm, ~, s_m] = ttest2(subjectInfo.networkTE(subjectInfo.sexSU == 'M-1' & inds, i), ...
                             subjectInfo.networkTE(subjectInfo.sexSU == 'M-0' & inds, i));
    t_stats_males(i) = s_m.tstat;
    p_values_males(i) = pm;
end

% Correct p-values
p_corr_females = mafdr(p_values_females, 'BHFDR', true);
p_corr_males = mafdr(p_values_males, 'BHFDR', true);

sig_females = makeSigCategory(p_values_females, p_corr_females);
sig_males = makeSigCategory(p_values_males, p_corr_males);

%% --- Plot 4: t-statistics of Network TE (Females vs Males)

clear g;
f = figure('Position', [100,100,1500,700]);

% Females
g(1,1) = gramm('x', networks_cat, 'y', t_stats_females, 'color', sig_females);
g(1,1).set_order_options('x', {'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
g(1,1).coord_flip();
g(1,1).geom_bar();
g(1,1).set_text_options('base_size', 25, 'Font', 'calibri');
g(1,1).axe_property('TickDir','out', 'Ygrid','on', 'GridColor',[0.5,0.5,0.5], ...
    'YLim', [-max(abs([t_stats_females t_stats_males]))-0.05, max(abs([t_stats_females t_stats_males]))+0.05]);
g(1,1).set_color_options('map', [[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706]]);
g(1,1).set_names('x','','y','t-statistic (Females: Future SU 1 vs 0)', 'color', 'Significance');

% Males
g(1,2) = gramm('x', networks_cat, 'y', t_stats_males, 'color', sig_males);
g(1,2).set_order_options('x', {'cer','sub','dmn','fpn','lim','vat','dat','som','vis'});
g(1,2).coord_flip();
g(1,2).geom_bar();
g(1,2).set_text_options('base_size', 25, 'Font', 'calibri');
g(1,2).axe_property('TickDir','out', 'Ygrid','on', 'GridColor',[0.5,0.5,0.5], ...
    'YLim', [-max(abs([t_stats_females t_stats_males]))-0.05, max(abs([t_stats_females t_stats_males]))+0.05]);
g(1,2).set_color_options('map', [[1,1,1];[0.651,0.808,0.890];[0.122,0.471,0.706]]);
g(1,2).set_names('x','','y','t-statistic (Males: Future SU 1 vs 0)', 'color', 'Significance');

g.draw();
exportgraphics(f, fullfile(figureDir, 'Fig_tstats_networks.png'));
