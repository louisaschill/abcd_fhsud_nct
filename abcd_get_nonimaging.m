% abcd_extract_nonimaging.m

% Louisa Schilling - last updated 12/24 
% This script processes and compiles non-imaging (behavioral, demographics, 
% family history) data from all subjects of the ABCD dataset for baseline visit 
% into a consolidated table for analysis. It extracts variables related to 
% demographics, behavior, and family history.

% Loads CSV files from the ABCD nonimaing data 
% Output: saves a table (`subjectInfo`) saved as a .mat file

% Dependencies: `extractInstrument`, `mergeTables`

clear all; close all;

subjs = 'all'; event = 'baseline_year_1_arm_1'; 

basedir = '/Users/louisaschilling/Desktop/Datasets/ABCD/Data/non_imaging/Release 5.1/core/'; 
cd(basedir); 

%% MRI (model, manufacturer and serial number) 
mri_var = {'src_subject_id','mri_info_manufacturer','mri_info_manufacturersmn','mri_info_deviceserialnumber'};
mri_var_names = {'subjectkey','manufacturer','model','serial_number'};
mri = readtable('/Users/louisaschilling/Desktop/Datasets/ABCD/Data/imaging/mri_y_adm_info.csv');
mri = extractInstrument(mri, mri_var, mri_var_names, event, subjs);
mri.manufacturer = categorical(mri.manufacturer);
mri.model = categorical(mri.model);


%% Site, family ID and age 

site_var = {'src_subject_id','site_id_l','rel_family_id','interview_age'};
site_var_names = {'subjectkey','site','rel_family_id','age'};
site = readtable([basedir '/abcd-general/abcd_y_lt.csv']);
site = extractInstrument(site, site_var, site_var_names, event, subjs);
site.site = categorical(site.site);

%% Sex and gender 
sex_var = {'src_subject_id','demo_sex_v2','demo_gender_id_v2'};
sex_var_names = {'subjectkey','sex','gender'};
sex = readtable([basedir '/gender-identity-sexual-health/gish_p_gi.csv']);
sex = extractInstrument(sex, sex_var, sex_var_names, event, subjs);
sex.sex= categorical(sex.sex); sex.gender = categorical(sex.gender);  

% Recode 1 = M, 2 = F, 3 = IM (intersex male), 4 = IF (intersex female)
sex.sex(sex.sex== '1') = categorical("M"); 
sex.sex(sex.sex=='2') = categorical("F");
sex.sex(sex.sex== '3') = categorical("IM"); 
sex.sex(sex.sex== '4') = categorical("IF"); 

% Recode 1 = M, 2 = F, 3 = TM (trans-male), 4 = TF (trans-female), 5 = GQ (gender queer) 
sex.gender(sex.gender == '1') = categorical("M"); 
sex.gender(sex.gender == '2') = categorical("F");
sex.gender(sex.gender == '3') = categorical("TM"); 
sex.gender(sex.gender == '4') = categorical("TF"); 
sex.gender(sex.gender == '5') = categorical("GQ"); 

%% Demographics: race, hispanic, income, parent education, adopted 
% Race:  1 = White; 2 = Black; 3 = Hispanic; 4 = Asian; 5 = Other

dem_var = {'src_subject_id','eventname','race_ethnicity','demo_comb_income_v2',...
    'demo_prnt_income_v2','demo_prtnr_income_v2',...
    'demo_prnt_ed_v2', 'demo_prtnr_ed_v2','demo_prim'}; 
dem_var_names = {'subjectkey','eventname','race','comb_income',...
    'parent_income', 'partner_income', 'parentEd', 'partnerEd', 'prim_id'};
dem = readtable([basedir '/abcd-general/abcd_p_demo.csv']);
dem = extractInstrument(dem, dem_var, dem_var_names, event, subjs);

% Household income (use combined income or if unavailable use one of parent)
dem.income = dem.comb_income;
missing_comb_income = isnan(dem.income); % Logical array for missing combined income
dem.income(missing_comb_income) = dem.parent_income(missing_comb_income);
missing_comb_income = isnan(dem.income);
dem.income(missing_comb_income) = dem.partner_income(missing_comb_income);

% Parent Education level: Recode parent education: 21 levels --> 5 categories:
% 0 = Never attended/Kindergarten only, 1 = 1st grade; 2 = 2nd grd; 
% 3 = 3rd grd; 4 = 4th grade; 5 = 5th grade
% 6 = 6th grade; 7 = 7th grade; 8 = 8th grade; 9 = 9th grade; 
% 10 = 10th grade; 11 = 11th grade; 12 = 12th grade;
% 13 = High school graduate; 14 = GED/equiv. diploma; 15 = Some college; 
% 16 = Associate degree: occupational; 
% 17 = Associate degree: academic; 18 = Bachelor's degree; 
% 19 = Master's degree; 20 = Professional School degree

% Recode to: %1: < High School, 2: High School/GED, 3: Some College,
%4: Associate/Bachelor, 5: Postgrad

parentEd_cat = nan(height(dem),1);
for i = 1:height(dem)
    parentEd = max(double(dem.parentEd(i)), double(dem.partnerEd(i))); % take higher of two 
    if parentEd < 13
        parentEd_cat(i) = 1;
    elseif parentEd == 13 || parentEd == 14
        parentEd_cat(i) = 2;
    elseif parentEd == 15
        parentEd_cat(i) = 3;
    elseif parentEd >= 16 && parentEd <= 18
        parentEd_cat(i) = 4;
    elseif parentEd >= 19
        parentEd_cat(i) = 5; 
    end 
end 
dem.parentEd_cat = categorical(parentEd_cat); 

% Household income: Recode household income recoded 10 --> 3 levels 
% 1= <$5,000, 2=$5,000-$11,999, 3=$12,000-$15,999, 
% 4=$16,000-$24,999, 5=$25,000- $34,999, 6=$35,000- $49,999,
% 7=$50,000-$74,999; 8= $75,000 - $99,999, 9=$100,000 - $199,999, 10=>$200,000. 
% Recode to: <50 K, 50-<100 K, ≥100 K

income_cat = nan(height(dem),1);
for i = 1:height(dem)
    income = max(dem.parent_income(i), dem.partner_income(i)); 
    income = max(income, dem.income(i));
   
    if income <= 6
        income_cat(i) = 1;
    elseif income >= 7 && income <= 8
        income_cat(i) = 2;
    elseif income == 9 || income == 10
        income_cat(i) = 3;
    end
end
dem.income_cat = categorical(income_cat); 

%% Family history of SUD and other mental illness 
% FH+ if 1+ parents and/or 2+ grandparents have illness (even if missing some info)
% FH- if not missing any info and 0 parents/grandparents have illness
% FH+/- if not missing any info and 1 grandparent only has illness 
% NaN if missing any info and not qualify as FH+ 

fhxssp_var = {'src_subject_id', 'famhx_ss_fath_prob_alc_p','famhx_ss_fath_prob_dg_p','famhx_ss_moth_prob_alc_p',...
    'famhx_ss_moth_prob_dg_p','famhx_ss_patgf_prob_alc_p','famhx_ss_patgf_prob_dg_p', 'famhx_ss_patgm_prob_alc_p',...
    'famhx_ss_patgm_prob_dg_p', 'famhx_ss_matgf_prob_alc_p','famhx_ss_matgf_prob_dg_p','famhx_ss_matgm_prob_alc_p',...
    'famhx_ss_matgm_prob_dg_p','famhx_ss_momdad_scd_p','famhx_ss_momdad_dprs_p','famhx_ss_momdad_ma_p', 'famhx_ss_momdad_prf_p',...
    'famhx_ss_momdad_nrv_p', 'famhx_ss_momdad_hspd_p', 'famhx_ss_momdad_vs_p', 'famhx_ss_momdad_trb_p'}; 
fhxssp_var_names = {'subjectkey','fatherAlc','fatherDrug','motherAlc','motherDrug','patGFAlc',...
    'patGFDrug','patGMAlc','patGMDrug', 'matGFAlc', 'matGFDrug', 'matGMAlc', 'matGMDrug',...
    'parentSCD', 'parentDEP','parentMA', 'parentPRF', 'parentNRV', 'parentHSPD', 'parentVS', 'parentTRB'}; 
fhxssp = readtable([basedir '/mental-health/mh_p_fhx.csv']);
fhxssp = extractInstrument(fhxssp, fhxssp_var, fhxssp_var_names,event,subjs);

fhxssp.missingSUD = isnan(fhxssp.fatherAlc)| isnan(fhxssp.fatherDrug)|...
    isnan(fhxssp.motherAlc)|isnan(fhxssp.motherDrug)| ...
    isnan(fhxssp.patGFAlc)| isnan(fhxssp.patGFDrug)|...
    isnan(fhxssp.patGMAlc)| isnan(fhxssp.patGMDrug)|...
    isnan(fhxssp.matGFAlc)| isnan(fhxssp.matGFDrug)|...
    isnan(fhxssp.matGMAlc)|isnan(fhxssp.matGMDrug);
fhxssp.missingAUD = isnan(fhxssp.fatherAlc)|isnan(fhxssp.motherAlc)| isnan(fhxssp.patGFAlc)| ...
    isnan(fhxssp.patGMAlc)| isnan(fhxssp.matGFAlc)| isnan(fhxssp.matGMAlc);
fhxssp.missingDUD = isnan(fhxssp.fatherDrug)|isnan(fhxssp.motherDrug)| isnan(fhxssp.patGFDrug)| ...
    isnan(fhxssp.patGMDrug)| isnan(fhxssp.matGFDrug)| isnan(fhxssp.matGMDrug);

fhxssp.momAlc = fhxssp.motherAlc > 0; fhxssp.momDrug = fhxssp.motherDrug > 0;
fhxssp.dadAlc = fhxssp.fatherAlc > 0; fhxssp.dadDrug = fhxssp.fatherDrug > 0;

fhxssp.parentAlc = fhxssp.motherAlc > 0 | fhxssp.fatherAlc > 0;
fhxssp.parentDrug = fhxssp.motherDrug > 0 | fhxssp.fatherDrug > 0;

fhxssp.grandparentAlc = (fhxssp.patGFAlc > 0) + (fhxssp.matGFAlc > 0) + ...
    (fhxssp.patGMAlc > 0) + (fhxssp.matGMAlc > 0);

fhxssp.grandparentDrug = (fhxssp.patGFDrug > 0) + (fhxssp.matGFDrug > 0) + ...
    (fhxssp.patGMDrug > 0) + (fhxssp.matGMDrug > 0); 

fhxssp.dadSUD = (fhxssp.fatherAlc) > 0 | (fhxssp.fatherDrug) > 0;
fhxssp.momSUD = (fhxssp.motherAlc) > 0 | (fhxssp.motherDrug) > 0;
fhxssp.parentSUD = fhxssp.momSUD + fhxssp.dadSUD; 

fhxssp.patGpSUD = (fhxssp.patGFAlc) > 0 | (fhxssp.patGFDrug) > 0;
fhxssp.patGmSUD= (fhxssp.patGMAlc) > 0 | (fhxssp.patGMDrug) > 0;
fhxssp.matGpSUD = (fhxssp.matGFAlc) > 0 | (fhxssp.matGFDrug) > 0;
fhxssp.matGmSUD = (fhxssp.matGMAlc) > 0 | (fhxssp.matGMDrug) > 0;
fhxssp.grandparentsSUD = fhxssp.patGpSUD + fhxssp.patGmSUD + fhxssp.matGpSUD + fhxssp.matGmSUD;

% FH SUD
fhxssp.typeSUD = nan(height(fhxssp),1);
for s = 1:height(fhxssp)
    fhxssp.FHSUD(s) = calculateFH(fhxssp.parentSUD(s), fhxssp.grandparentsSUD(s), fhxssp.missingSUD(s), 1, 2);
    fhxssp.FHAUD(s) = calculateFH(fhxssp.parentAlc(s), fhxssp.grandparentAlc(s), fhxssp.missingAUD(s), 1, 2);
    fhxssp.FHDUD(s) = calculateFH(fhxssp.parentDrug(s), fhxssp.grandparentDrug(s), fhxssp.missingDUD(s), 1, 2);

    if fhxssp.FHSUD(s) == "FH-"
        fhxssp.typeSUD(s) = categorical("FH-");
    else
        if fhxssp.FHDUD(s) == "FH+" && fhxssp.FHAUD(s)  == "FH+"
            fhxssp.typeSUD(s) = categorical("FH+-Both");
        elseif fhxssp.FHDUD(s) == "FH+" && fhxssp.FHAUD(s) ~= "FH+"
            fhxssp.typeSUD(s) = categorical("FH+-DUD");
        elseif fhxssp.FHDUD(s) ~= "FH+" && fhxssp.FHAUD(s) == "FH+"
            fhxssp.typeSUD(s) = categorical("FH+-AUD");
        else
            fhxssp.typeSUD(s) = categorical("Multi");
        end
    end
end

% Family density of SUD (grandpart = 0.5, parent = 1) 
grandparentSUD = fhxssp.grandparentsSUD * 0.5; 
parentSUD = fhxssp.parentSUD; 
fhxssp.familydensitySUD = parentSUD + grandparentSUD; 

%% Parental history of MI (binary):
% Suicide, depression, mania, seen a counselor for mental health, nervous breakdown,
% hospitalized for emotional/mental problem, visions of plotting, trouble with job/police/fights 
missingParentMH = find(isnan(fhxssp.parentSCD)|isnan(fhxssp.parentDEP)|...
    isnan(fhxssp.parentMA)| isnan(fhxssp.parentPRF) | ...
    isnan(fhxssp.parentNRV)| isnan(fhxssp.parentHSPD)|...
    isnan(fhxssp.parentVS)| isnan(fhxssp.parentTRB));

fhxssp.parentMH = (fhxssp.parentSCD)>0 | (fhxssp.parentDEP) >0 | ...
    (fhxssp.parentMA) > 0 | (fhxssp.parentPRF) >0 | (fhxssp.parentNRV) >0 | ...
    (fhxssp.parentHSPD) >0 | (fhxssp.parentVS) >0 | (fhxssp.parentTRB) > 0';
fhxssp.parentMH = double(fhxssp.parentMH); 
for i = 1:height(fhxssp)
    if fhxssp.parentMH(i) == 0
        if ismember(i,missingParentMH)
            fhxssp.parentMH(i) = NaN; 
        end 
    end 
end 

%% In utero substance use - knowing of pregnancy 
% other drugs categories: 0 = None; 
% 1 = Amphetamines/ methamphetamine 
% 2 = Benzodiazepines; 
% 3 = Caffeine; 
% 4 = Cathinones (bath salts); 5 = Fake/synthetic marijuana;
% 6 = GHB; 7 = Hallucinogens (LSD or acid); 
% 8 = Inhalants; 
% 9 = Ketamine (special K); 
% 10 = MDMA (ecstasy; 
% 11 = Opioids; 
% 12 = Other;
% 13 = Barbituates; 

dhxs_var = {'src_subject_id','devhx_9_tobacco','devhx_9_alcohol' 'devhx_9_marijuana',...
    'devhx_9_her_morph', 'devhx_9_oxycont','devhx_9_other_drugs','devhx_9_other1_name_2',...
    'devhx_9_other2_name_2','devhx_9_other3_name_2','devhx_9_other4_name_2','devhx_9_other5_name_2','devhx_ss_9_alcohol_max_p',...
    'devhx_ss_9_alcohol_avg_p','devhx_caffeine_11'}; 
dhxs_var_names = {'subjectkey','pregtobacco','pregalcohol','pregmarijuana', 'pregmorph','pregoxy','pregother',...
    'pregother_name1','pregother_name2','pregother_name3','pregother_name4',...
    'pregother_name5','pregalcmax', 'pregalcweek','pregcaffeine'};
dhxs = readtable([basedir '/physical-health/ph_p_dhx.csv']);
dhxs = extractInstrument(dhxs, dhxs_var, dhxs_var_names,event,subjs);

dhxs.missingPregSub = isnan(dhxs.pregtobacco)|isnan(dhxs.pregalcohol)|...
    isnan(dhxs.pregmarijuana)|isnan(dhxs.pregmorph)|...
    isnan(dhxs.pregoxy)|isnan(dhxs.pregother);

dhxs.subDuringPreg = (dhxs.pregtobacco) == 1 | (dhxs.pregalcohol == 1) | ...
    (dhxs.pregmarijuana) == 1 | (dhxs.pregmorph) == 1 | ...
    (dhxs.pregoxy) == 1 | (dhxs.pregother) == 1; 

%% Child substance use  

% KSADS SUD  - parent reported substance use 
ksads_sud_var = {'src_subject_id','ksads_aud_raw_581_p','ksads_aud_raw_2115_p',...
    'ksads_aud_raw_587_p','ksads_dud_raw_614_p','ksads_dud_raw_615_p',...
    'ksads_dud_raw_616_p','ksads_dud_raw_617_p','ksads_dud_raw_618_p','ksads_dud_raw_619_p',...
    'ksads_dud_raw_620_p','ksads_dud_raw_621_p','ksads_dud_raw_622_p'};
ksads_sud_var_names = {'subjectkey','threedrinksday','twodrinks12months','drinkproblems',...
    'marijuana','stimulants','sedatives','cocaine','opioids','hallucinogens','tobacco',...
    'solvents','otherdrugs'}; 
ksads_sud = readtable([basedir '/substance-use/su_p_ksads_sud.csv']);
ksads_sud = extractInstrument(ksads_sud, ksads_sud_var, ksads_sud_var_names, event, subjs);

% Calculate the number of substances used
any_sub_p = nansum(ksads_sud{:, 2:end}, 2); % Sum substance use indicators for each row, ignoring NaNs
ksads_sud.num_subs = any_sub_p; % Add the total number of substances as a new column

% Create a binary variable for any substance use
ksads_sud.any_sub_parent = any_sub_p > 0; % Logical array indicating if any substance was used

% SU - child reported substance use 
su_var = {'src_subject_id','isip_2_2','tlfb_alc_use','tlfb_cig_use','tlfb_ecig_use',...
    'tlfb_chew_use','tlfb_cigar_use','tlfb_hookah_use','tlfb_pipes_use',...
    'tlfb_nicotine_use','tlfb_mj_use','tlfb_edible_use','tlfb_mj_conc_use',...
    'tlfb_mj_drink_use','tlfb_tincture_use','tlfb_mj_synth_use','tlfb_coc_use','tlfb_meth_use',...
    'tlfb_mdma_use','tlfb_ket_use','tlfb_ghb_use','tlfb_opi_use','tlfb_hall_use','tlfb_shrooms_use',...
    'tlfb_salvia_use','tlfb_steroids_use','tlfb_bitta_use',...
    'tlfb_sniff_use','tlfb_inhalant_use','tlfb_amp_use','tlfb_tranq_use',...
    'tlfb_vicodin_use','tlfb_cough_use','tlfb_other_use'};

su_var_names = {'subjectkey','cont_drinking','full_drink','cig_more_puff',...
    'ecig_more_puff','tlfb_chew_use','tlfb_cigar_use','tlfb_hookah_use','tlfb_pipes_use',...
    'tlfb_nicotine_use','tlfb_mj_use','tlfb_edible_use','tlfb_mj_conc_use','tlfb_mj_drink_use',...
    'tlfb_tincture_use','tlfb_mj_synth_use','tlfb_coc_use','tlfb_meth_use','tlfb_mdma_use','tlfb_ket_use',...
    'tlfb_ghb_use','tlfb_opi_use','tlfb_hall_use','tlfb_shrooms_use','tlfb_salvia_use','tlfb_steroids_use','tlfb_bitta_use',...
    'tlfb_sniff_use','tlfb_inhalant_use','tlfb_amp_use','tlfb_tranq_use','tlfb_vicodin_use','tlfb_cough_use','tlfb_other_use'};

su = readtable([basedir '/substance-use/su_y_sui.csv']);
su = extractInstrument(su, su_var, su_var_names,event,subjs);
any_sub = nansum(su(:,2:end),2);
any_sub = table2array(any_sub);
su.num_subs = any_sub; 
su.any_sub = any_sub > 0 ;

%% Behavioral/environmental variables of potential interest 

% CBCL 
cbcl_var = {'src_subject_id','cbcl_scr_syn_anxdep_r', 'cbcl_scr_syn_withdep_r','cbcl_scr_syn_somatic_r','cbcl_scr_dsm5_somaticpr_r',...
    'cbcl_scr_syn_social_r', 'cbcl_scr_syn_thought_r', 'cbcl_scr_syn_attention_r'...
    'cbcl_scr_syn_rulebreak_r', 'cbcl_scr_syn_aggressive_r', 'cbcl_scr_dsm5_opposit_r',...
    'cbcl_scr_07_sct_r', 'cbcl_scr_07_ocd_r','cbcl_scr_07_stress_r','cbcl_scr_dsm5_adhd_r',...
    'cbcl_scr_syn_totprob_r','cbcl_scr_dsm5_depress_r','cbcl_scr_dsm5_anxdisord_r',...
    'cbcl_scr_dsm5_conduct_r','cbcl_scr_syn_internal_r','cbcl_scr_syn_external_r','cbcl_scr_syn_internal_t', 'cbcl_scr_syn_external_t',...
    'cbcl_q56a_p','cbcl_q56b_p','cbcl_q56f_p'};
cbcl_var_names = {'subjectkey', 'anxdep', 'withdep', 'somatic','somatic_dsm', 'social', 'thought',...
    'attention','rulebreak','aggressive','opposit', 'sct', 'ocd', 'stress','adhd',...
    'totprob', 'depress','anxdisord','conduct', 'internal', 'external','internal_t','external_t','somatic_aches','somatic_headaches', 'somatic_stomach'};
cbcl = readtable([basedir '/mental-health/mh_p_cbcl.csv']);
cbcl = extractInstrument(cbcl, cbcl_var, cbcl_var_names, event, subjs);

% CBCL - follow up 
cbcl_var = {'src_subject_id','cbcl_scr_syn_somatic_r','cbcl_scr_dsm5_somaticpr_r',...
    'cbcl_q56a_p','cbcl_q56b_p','cbcl_q56f_p'};
cbcl_var_names = {'subjectkey','somatic_1yr','somatic_dsm_1yr',...
    'somatic_aches_1yr','somatic_headaches_1yr', 'somatic_stomach_1yr'};
cbcl2 = readtable([basedir '/mental-health/mh_p_cbcl.csv']);
cbcl2 = extractInstrument(cbcl2, cbcl_var, cbcl_var_names, '1_year_follow_up_y_arm_1', subjs);

cbcl_var = {'src_subject_id','cbcl_scr_syn_somatic_r','cbcl_scr_dsm5_somaticpr_r',...
    'cbcl_q56a_p','cbcl_q56b_p','cbcl_q56f_p'};
cbcl_var_names = {'subjectkey','somatic_2yr','somatic_dsm_2yr', ...
    'somatic_aches_2yr','somatic_headaches_2yr', 'somatic_stomach_2yr'};
cbcl3 = readtable([basedir '/mental-health/mh_p_cbcl.csv']);
cbcl3 = extractInstrument(cbcl3, cbcl_var, cbcl_var_names, '2_year_follow_up_y_arm_1', subjs);

cbcl_var = {'src_subject_id','cbcl_scr_syn_somatic_r','cbcl_scr_dsm5_somaticpr_r',...
    'cbcl_q56a_p','cbcl_q56b_p','cbcl_q56f_p'};
cbcl_var_names = {'subjectkey','somatic_3yr','somatic_dsm_3yr', ...
    'somatic_aches_3yr','somatic_headaches_3yr', 'somatic_stomach_3yr'};
cbcl4 = readtable([basedir '/mental-health/mh_p_cbcl.csv']);
cbcl4 = extractInstrument(cbcl4, cbcl_var, cbcl_var_names, '3_year_follow_up_y_arm_1', subjs);

% EHIS - handedness
ehis_var = {'src_subject_id','ehi_y_ss_scoreb'}; 
ehis_var_names = {'subjectkey','handedness'};
ehis = readtable([basedir '/neurocognition/nc_y_ehis.csv']);
ehis = extractInstrument(ehis, ehis_var, ehis_var_names, event, subjs);

% Screen time (parent reported)
st_p_var = {'src_subject_id','screentime1_p_hours','screentime1_p_minutes','screentime2_p_hours','screentime2_p_minutes'}; 
st_p_var_names = {'subjectkey','screen_wkday_hours','screen_wkday_mins','screen_wknd_hours','screen_wknd_mins' }; 
st_p = readtable([basedir, '/novel-technologies/nt_p_stq.csv']);
st_p = extractInstrument(st_p, st_p_var, st_p_var_names, event, subjs);
st_p.screen_week = (st_p.screen_wkday_hours*60) +  st_p.screen_wkday_mins +... 
(st_p.screen_wknd_hours*60)+ st_p.screen_wknd_mins; 
lowThreshold = prctile(st_p.screen_week, 33);
highThreshold = prctile(st_p.screen_week, 66);
edges = [0, lowThreshold, highThreshold, Inf]; % Define the edges of the bins
categories = ["1", "2", "3"]; % Category labels
st_p.screentime_category = discretize(st_p.screen_week, edges, 'Categorical', categories);

% UPPS 
upps_var = {'src_subject_id','upps_y_ss_negative_urgency','upps_y_ss_lack_of_planning',...
'upps_y_ss_sensation_seeking','upps_y_ss_positive_urgency','upps_y_ss_lack_of_perseverance'}; 
upps_var_names = {'subjectkey','upps_neg_urgency', 'upps_lack_plan',...
    'upps_sensation_seek', 'upps_pos_urgency','upps_lack_preservance'};
upps = readtable([basedir, '/mental-health/mh_y_upps.csv']);
upps = extractInstrument(upps, upps_var, upps_var_names, event, subjs);
upps.upps_sum = upps.upps_sensation_seek + upps.upps_lack_preservance + upps.upps_neg_urgency + upps.upps_pos_urgency + upps.upps_lack_plan;

% RES 
res_coi_var = {'src_subject_id','reshist_addr1_coi_z_ed_nat','reshist_addr1_coi_z_he_nat',...
'reshist_addr1_coi_z_se_nat'}; 
res_coi_var_names = {'subjectkey','z_edu','z_health_env', 'z_social_eco'};
res_coi = readtable([basedir, '/linked-external-data/led_l_coi.csv']);
res_coi = extractInstrument(res_coi, res_coi_var, res_coi_var_names, event, subjs);

% RES 
res_adi_var = {'src_subject_id','reshist_addr1_adi_edu_h','reshist_addr1_adi_income'}; 
res_adi_var_names = {'subjectkey', 'high_school_dip','median_income'};
res_adi = readtable([basedir, '/linked-external-data/led_l_adi.csv']);
res_adi = extractInstrument(res_adi, res_adi_var, res_adi_var_names, event, subjs);

% NSC 
nsc_var = {'src_subject_id','neighborhood_crime_y'}; 
nsc_var_names = {'subjectkey','neigh_crime'}; 
nsc = readtable([basedir, '/culture-environment/ce_y_nsc.csv']);
nsc = extractInstrument(nsc, nsc_var, nsc_var_names, event, subjs);

% PM child  
pm_y_var = {'src_subject_id','pmq_y_ss_mean'}; 
pm_y_var_names = {'subjectkey','parent_mon_y'}; 
pm_y = readtable([basedir, '/culture-environment/ce_y_pm.csv']);
pm_y = extractInstrument(pm_y, pm_y_var,pm_y_var_names, event, subjs);

% PM parent  
pm_p_var = {'src_subject_id','parental_monitor_ss_mean'}; 
pm_p_var_names = {'subjectkey','parent_mon_p'}; 
pm_p = readtable([basedir, '/culture-environment/ce_p_pm.csv']);
pm_p = extractInstrument(pm_p, pm_p_var,pm_p_var_names, event, subjs);

% FES 
fes_var = {'src_subject_id','fes_y_ss_fc','fes_y_ss_fc_pr'}; 
fes_var_names = {'subjectkey','fam_env_conf_raw', 'fam_env_conf_pr'}; 
fes = readtable([basedir, '/culture-environment/ce_y_fes.csv']);
fes = extractInstrument(fes, fes_var, fes_var_names, event, subjs);

% PSB Y  
psb_y_var = {'src_subject_id','psb_y_ss_mean'}; 
psb_y_var_names = {'subjectkey','prosoc_mean_y'}; 
psb_y = readtable([basedir, '/culture-environment/ce_y_psb.csv']);
psb_y = extractInstrument(psb_y, psb_y_var, psb_y_var_names, event, subjs);

% PSB P 
psb_p_var = {'src_subject_id','psb_p_ss_mean'}; 
psb_p_var_names = {'subjectkey','prosoc_mean_p'}; 
psb_p = readtable([basedir, '/culture-environment/ce_p_psb.csv']);
psb_p = extractInstrument(psb_p, psb_p_var, psb_p_var_names, event, subjs);

% CRPBI 
crpbi_var = {'src_subject_id','crpbi_y_ss_parent'}; 
crpbi_var_names = {'subjectkey','parent_ac_mean'}; 
crpbi = readtable([basedir, '/culture-environment/ce_y_crpbi.csv']);
crpbi = extractInstrument(crpbi, crpbi_var,crpbi_var_names, event, subjs);

% SRPF 
srpf_var = {'src_subject_id','srpf_y_ss_ses','srpf_y_ss_iiss','srpf_y_ss_dfs'}; 
srpf_var_names = {'subjectkey','school_env','school_inv','school_diseng'}; 
srpf = readtable([basedir, '/culture-environment/ce_y_srpf.csv']);
srpf = extractInstrument(srpf,srpf_var, srpf_var_names, event, subjs);

% PPS 
pps_var = {'src_subject_id','pps_y_ss_number','pps_y_ss_severity_score','pps_ss_mean_severity','pps_y_ss_number_nm'};
pps_var_names = {'subjectkey','total_psychosis', 'psychosis_severity','mean_pscyhosis_severity','missing_psychosis'};
pps = readtable([basedir '/mental-health/mh_y_pps.csv']);
pps = extractInstrument(pps, pps_var, pps_var_names, event, subjs);
pps.total_psychosis(pps.missing_psychosis > 6) = nan; 
pps.psychosis_severity(pps.missing_psychosis > 6) = nan; 

% PPS 
pps_var = {'src_subject_id','pps_y_ss_number','pps_y_ss_severity_score'};
pps_var_names = {'subjectkey','total_psychosis', 'psychosis_severity'};
pps = readtable([basedir '/mental-health/mh_y_pps.csv']);
pps = extractInstrument(pps, pps_var, pps_var_names, event, subjs);

% BIS BAS 
bisbas_var = {'src_subject_id','bis_y_ss_bis_sum','bis_y_ss_bas_rr',...
    'bis_y_ss_bas_drive','bis_y_ss_bas_fs','bis_y_ss_bism_sum','bis_y_ss_basm_rr',...
    'bis_y_ss_basm_drive'};
bisbas_var_names = {'subjectkey','bis_sum','bas_reward','bas_drive',...
    'bas_fun','bis_mod','bas_reward_mod','bas_drive_mod'};
bisbas = readtable([basedir '/mental-health/mh_y_bisbas.csv']);
bisbas = extractInstrument(bisbas, bisbas_var, bisbas_var_names, event, subjs);

% GBI
gbi_var = {'src_subject_id','pgbi_p_ss_score'};
gbi_var_names = {'subjectkey','mania'};
gbi = readtable([basedir '/mental-health/mh_p_gbi.csv']);
gbi = extractInstrument(gbi, gbi_var, gbi_var_names, event, subjs);

% BDEF 
bdef_var = {'src_subject_id', 'bdefs_calm_down_p','bdefs_consequences_p','bdefs_distract_upset_p',...
    'bdefs_explain_idea_p','bdefs_explain_pt_p','bdefs_explain_seq_p',...
    'bdefs_impulsive_action_p','bdefs_inconsistant_p','bdefs_lazy_p',...
    'bdefs_process_info_p','bdefs_rechannel_p','bdefs_sense_time_p','bdefs_shortcuts_p',...
    'bdefs_stop_think_p'};
bdef_var_names = {'subjectkey','calm','consequences','distract_upset','explain_idea','explain_pt',...
    'explain_seq','impulsive_action','inconsistant','bdefs_lazy_p','bdefs_process_info_p',...
    'rechannel','sense_time','shortcuts','stop_think'};
bdef = readtable([basedir, '/neurocognition/nc_p_bdef.csv']);
bdef = extractInstrument(bdef, bdef_var, bdef_var_names,'3_year_follow_up_y_arm_1', subjs);

% SST 
sst_var = {'src_subject_id', 'tfmri_sst_all_beh_crs_rt'};
sst_var_names = {'subjectkey', 'rate_stop_correct'}; 
sst = readtable([basedir '/imaging/mri_y_tfmr_sst_beh.csv']);
sst = extractInstrument(sst, sst_var, sst_var_names,event,subjs);

%% Puberty
pub_var = {'src_subject_id', 'pds_p_ss_male_category','pds_p_ss_female_category'};
pub_var_names = {'subjectkey', 'puberty_male', 'puberty_female'}; 
pub = readtable([basedir '/physical-health/ph_p_pds.csv']);
pub = extractInstrument(pub, pub_var, pub_var_names,event,subjs);
pub.puberty_combined = pub.puberty_male; % Initialize with male scores
nan_male = isnan(pub.puberty_male); % Find rows where male scores are NaN
pub.puberty_combined(nan_male) = pub.puberty_female(nan_male);

%% SPH: HORMONES 
sph_var = {'src_subject_id','hormone_sal_sex','hormon_sal_notes_y___1','hormon_sal_notes_y___2','hormon_sal_notes_y___3',...
    'hormon_sal_notes_y___4','hormon_sal_notes_y___5','hormon_sal_notes_y___6','hormone_scr_dhea_mean','hormone_scr_dhea_rep1','hormone_scr_dhea_rep1_ll','hormone_scr_dhea_rep1_qns',...
    'hormone_scr_dhea_rep1_nd','hormone_scr_dhea_rep2','hormone_scr_dhea_rep2_ll','hormone_scr_dhea_rep2_qns',...
    'hormone_scr_dhea_rep2_nd','hormone_scr_hse_mean','hormone_scr_hse_rep1','hormone_scr_hse_rep1_ll','hormone_scr_hse_rep1_qns',...
    'hormone_scr_hse_rep1_nd','hormone_scr_hse_rep2','hormone_scr_hse_rep2_ll','hormone_scr_hse_rep2_qns',...
    'hormone_scr_hse_rep2_nd','hormone_scr_ert_mean','hormone_scr_ert_rep1','hormone_scr_ert_rep1_ll','hormone_scr_ert_rep1_qns',...
    'hormone_scr_ert_rep1_nd','hormone_scr_ert_rep2','hormone_scr_ert_rep2_ll','hormone_scr_ert_rep2_qns',...
    'hormone_scr_ert_rep2_nd'}; 
sph_var_names = {'subjectkey','hormone_sex','hormon_sal_notes_y___1','hormon_sal_notes_y___2','hormon_sal_notes_y___3',...
    'hormon_sal_notes_y___4','hormon_sal_notes_y___5','hormon_sal_notes_y___6','hormone_scr_dhea_mean','hormone_scr_dhea_rep1','hormone_scr_dhea_rep1_ll','hormone_scr_dhea_rep1_qns',...
    'hormone_scr_dhea_rep1_nd','hormone_scr_dhea_rep2','hormone_scr_dhea_rep2_ll','hormone_scr_dhea_rep2_qns',...
    'hormone_scr_dhea_rep2_nd','hormone_scr_hse_mean','hormone_scr_hse_rep1','hormone_scr_hse_rep1_ll','hormone_scr_hse_rep1_qns',...
    'hormone_scr_hse_rep1_nd','hormone_scr_hse_rep2','hormone_scr_hse_rep2_ll','hormone_scr_hse_rep2_qns',...
    'hormone_scr_hse_rep2_nd','hormone_scr_ert_mean','hormone_scr_ert_rep1','hormone_scr_ert_rep1_ll','hormone_scr_ert_rep1_qns',...
    'hormone_scr_ert_rep1_nd','hormone_scr_ert_rep2','hormone_scr_ert_rep2_ll','hormone_scr_ert_rep2_qns',...
    'hormone_scr_ert_rep2_nd'}; 
sph = readtable([basedir '/physical-health/ph_y_sal_horm.csv']);
sph = extractInstrument(sph, sph_var, sph_var_names, event, subjs);
for i = 1:height(sph)
    if sph.hormone_sex(i) == 1; sph.horm_sex(i) = categorical("F");
    elseif sph.hormone_sex(i) == 2; sph.horm_sex(i) = categorical("M");
    end 
end 

% Check sample notes and set to nan if any issues with saliva sample 
inds = sph.hormon_sal_notes_y___1 == 0 |(sph.hormon_sal_notes_y___2 + sph.hormon_sal_notes_y___3 + sph.hormon_sal_notes_y___4 +...
     sph.hormon_sal_notes_y___5 + sph.hormon_sal_notes_y___6) > 0;
for i= 3:size(sph,2)-1; sph{inds, i} = nan; end 
horm = table;  horm.subjectkey = sph.subjectkey; horm.horm_sex = sph.horm_sex;

% DHEA 
horm.dhea = NaN(height(horm),1);
inds_rep1 = (sph.hormone_scr_dhea_rep1_ll + sph.hormone_scr_dhea_rep1_nd + sph.hormone_scr_dhea_rep1_qns) == 0; 
inds_rep2 = (sph.hormone_scr_dhea_rep2_ll + sph.hormone_scr_dhea_rep2_nd + sph.hormone_scr_dhea_rep2_qns) == 0; 
horm.dhea(inds_rep1 & inds_rep2) = sph.hormone_scr_dhea_mean(inds_rep1 & inds_rep2); 
horm.dhea((~inds_rep1) & inds_rep2) = sph.hormone_scr_dhea_rep2((~inds_rep1) & inds_rep2);
horm.dhea((~inds_rep2) & inds_rep1) = sph.hormone_scr_dhea_rep1((~inds_rep2) & inds_rep1);

% Estradiol
horm.hse = NaN(height(horm),1);
inds_rep1 = (sph.hormone_scr_hse_rep1_ll + sph.hormone_scr_hse_rep1_nd + sph.hormone_scr_hse_rep1_qns) == 0; 
inds_rep2 = (sph.hormone_scr_hse_rep2_ll + sph.hormone_scr_hse_rep2_nd + sph.hormone_scr_hse_rep2_qns) == 0; 
horm.hse(inds_rep1 & inds_rep2) = sph.hormone_scr_hse_mean(inds_rep1 & inds_rep2); 
horm.hse((~inds_rep1) & inds_rep2) = sph.hormone_scr_hse_rep2((~inds_rep1) & inds_rep2);
horm.hse((~inds_rep2) & inds_rep1) = sph.hormone_scr_hse_rep1((~inds_rep2) & inds_rep1);

% testosterone 
horm.ert = NaN(height(horm),1);
inds_rep1 = (sph.hormone_scr_ert_rep1_ll + sph.hormone_scr_ert_rep1_nd + sph.hormone_scr_ert_rep1_qns) == 0; 
inds_rep2 = (sph.hormone_scr_ert_rep2_ll + sph.hormone_scr_ert_rep2_nd + sph.hormone_scr_ert_rep2_qns) == 0; 
horm.ert(inds_rep1 & inds_rep2) = sph.hormone_scr_ert_mean(inds_rep1 & inds_rep2); 
horm.ert((~inds_rep1) & inds_rep2) = sph.hormone_scr_ert_rep2((~inds_rep1) & inds_rep2);
horm.ert((~inds_rep2) & inds_rep1) = sph.hormone_scr_ert_rep1((~inds_rep2) & inds_rep1);

%% Compile all data tables into mega subjectInfo table 
tabs = {site, mri, sex, dem,ehis,res_coi,res_adi,nsc, pm_y, pm_p, fes, ...
    psb_y, psb_p, crpbi, srpf, cbcl, cbcl2, cbcl3, cbcl4, pps, bisbas, gbi, bdef, fhxssp,...
    dhxs,st_p,su,pub,horm,ksads_sud,upps,sst}; 

subjectInfo = mergeTables(tabs);

sexFHSUD = strcat(char(subjectInfo.sex),'-',char(subjectInfo.FHSUD)); 
subjectInfo.sexFHSUD = categorical(cellstr(sexFHSUD)); 

%% SAVE subjectInfo table 
nsubj = height(subjectInfo); 
savedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; mkdir(savedir); cd(savedir);
save(fullfile(savedir,'nonimagingInfo_all_subjs_release51.mat'),'subjectInfo'); 
disp('Saved subjectInfo for all subjects.')
%% Helper functions 

function fhCategory = calculateFH(parentCount, grandparentCount, missing, thresholdParents, thresholdGrandparents)
    if missing
        if parentCount >= thresholdParents || grandparentCount >= thresholdGrandparents
            fhCategory = categorical({'FH+'});
        else
            fhCategory = categorical({'NaN'});
        end
    else
        if parentCount >= thresholdParents || grandparentCount >= thresholdGrandparents
            fhCategory = categorical({'FH+'});
        elseif parentCount == 0 && grandparentCount == 0
            fhCategory = categorical({'FH-'});
        elseif parentCount == 0 && grandparentCount == 1
            fhCategory = categorical({'FH+/-'});
        else
            fhCategory = categorical({'NaN'});
        end
    end
end
