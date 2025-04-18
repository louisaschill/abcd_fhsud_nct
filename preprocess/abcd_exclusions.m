%% Louisa Schilling -- updated December 2024 
% Exclusions for ABCD FH SUD Project
% This script applies various exclusion criteria to the ABCD dataset. 
clear all; close all;

%% Load All Subjects Information
basedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data'; cd(basedir);
load('nonimagingInfo_all_subjs_release51.mat'); % Loads `subjectInfo`
nTotal = size(subjectInfo, 1);
disp(['Total cohort N = ' num2str(nTotal)]);

% Initialize exclusion counter
nExcluded = 0;

%% Exclusion Criteria
% 1. Exclude subjects missing rsfMRI data
cd('/Users/louisaschilling/Desktop/Datasets/ABCD/Data');
load('subjects.mat'); % Load `subjkeys` with rsfMRI data availability
missingRsfmri = ~ismember(subjectInfo.subjectkey, subjkeys);
nExcluded = nExcluded + sum(missingRsfmri);
subjectInfo = subjectInfo(~missingRsfmri, :);
disp(['Excluded subjects missing rsfMRI data: n = ' num2str(sum(missingRsfmri))]);

% 2. Exclude Philips scanner subjects or undefined scanner type
notPhilips = subjectInfo.manufacturer ~= 'Philips Medical Systems';
notUndefined = ~isundefined(subjectInfo.manufacturer);
excludePhilips = ~(notPhilips & notUndefined);
nExcluded = nExcluded + sum(excludePhilips);
subjectInfo = subjectInfo(~excludePhilips, :);
disp(['Excluded Philips or undefined scanner type: n = ' num2str(sum(excludePhilips))]);

% 3. Exclude subjects missing FH SUD info (not FH+, FH-, nor FH+/-)
validFHSUD = ismember(subjectInfo.FHSUD, {'FH-', 'FH+', 'FH+/-'});
nExcluded = nExcluded + sum(~validFHSUD);
subjectInfo = subjectInfo(validFHSUD, :);
disp(['Excluded subjects missing FH SUD info: n = ' num2str(sum(~validFHSUD))]);

% 4. Exclude missing maternal SU 
nExcluded = nExcluded + sum(subjectInfo.missingPregSub);
subjectInfo = subjectInfo(~subjectInfo.missingPregSub, :);
disp(['Excluded missing info about in utero SU: n = ' num2str(sum(subjectInfo.missingPregSub))]);

% 5. Exclude FH- kids whose moms drank >7 drinks/week or >4 in one sitting
highMaternalDrinking = (subjectInfo.FHSUD == 'FH-' | subjectInfo.FHSUD == 'FH+/-') & ...
                       (subjectInfo.pregalcweek > 7 | subjectInfo.pregalcmax >= 4);
nExcluded = nExcluded + sum(highMaternalDrinking);
subjectInfo = subjectInfo(~highMaternalDrinking, :);
disp(['Excluded FH- kids with high maternal drinking: n = ' num2str(sum(highMaternalDrinking))]);

% 6. Exclude adopted children
isAdopted = subjectInfo.prim_id == 3; %| ~isnan(subjectInfo.age_adopted); % Primary caregiver is adoptive parent or report age of adoption 
nExcluded = nExcluded + sum(isAdopted);
subjectInfo = subjectInfo(~isAdopted, :);
disp(['Excluded adopted children: n = ' num2str(sum(isAdopted))]);

% 7. Exclude children using drugs/alcohol (parent-reported)
usingSubstancesParent = subjectInfo.any_sub_parent ~= 0;
nExcluded = nExcluded + sum(usingSubstancesParent);
subjectInfo = subjectInfo(~usingSubstancesParent, :);
disp(['Excluded children using drugs/alcohol (parent-reported): n = ' num2str(sum(usingSubstancesParent))]);

% 8. Exclude children using substances (child-reported)
usingSubstancesChild = subjectInfo.any_sub == 1;
nExcluded = nExcluded + sum(usingSubstancesChild);
subjectInfo = subjectInfo(~usingSubstancesChild, :);
disp(['Excluded children using substances (child-reported): n = ' num2str(sum(usingSubstancesChild))]);

% 9. Exclude subjects with mismatched saliva and biological sex
mismatchedSex = subjectInfo.sex ~= subjectInfo.horm_sex;
nExcluded = nExcluded + sum(mismatchedSex);
subjectInfo = subjectInfo(~mismatchedSex, :);
disp(['Excluded subjects with mismatched saliva and biological sex: n = ' num2str(sum(mismatchedSex))]);

% 10. Exclude missing parent income
missingIncome = ismissing(subjectInfo.income_cat);
nExcluded = nExcluded + sum(missingIncome);
subjectInfo = subjectInfo(~missingIncome, :);
disp(['Excluded subjects missing income info: n = ' num2str(sum(missingIncome))]);

% 11. Exclude missing parent education 
missingParentEd = ismissing(subjectInfo.parentEd_cat);
nExcluded = nExcluded + sum(missingParentEd);
subjectInfo = subjectInfo(~missingParentEd, :);
disp(['Excluded subjects missing parental education info: n = ' num2str(sum(missingParentEd))]);

% 12. Exclude missing parent mental health
missingParentMH = ismissing(subjectInfo.parentMH);
nExcluded = nExcluded + sum(missingParentMH);
subjectInfo = subjectInfo(~missingParentMH, :);
disp(['Excluded subjects missing parental mental health info: n = ' num2str(sum(missingParentMH))]);

% 13. Exclude missing puberty info 
missingPuberty = ismissing(subjectInfo.puberty_combined);
nExcluded = nExcluded + sum(missingPuberty);
subjectInfo = subjectInfo(~missingPuberty, :);
disp(['Excluded subjects missing puberty info: n = ' num2str(sum(missingPuberty))]);

%% Add extra SUD Info 
nsubjs = height(subjectInfo);

grandparentSUD = subjectInfo.grandparentsSUD * 0.5; 
parentSUD = subjectInfo.parentSUD; 
subjectInfo.familydensitySUD = parentSUD + grandparentSUD; 

subjectInfo.sex = categorical(string(subjectInfo.sex)); 
subjectInfo.FHSUD = categorical(string(subjectInfo.FHSUD)); 
subjectInfo.sexFHSUD = categorical(string(subjectInfo.sexFHSUD));
subjectInfo.manufacturer = categorical(string(subjectInfo.manufacturer));
subjectInfo.site = categorical(string(subjectInfo.site));

%% Final Summary
nIncluded = height(subjectInfo);
disp(['Total excluded subjects: n = ' num2str(nExcluded)]);
disp(['Total included subjects: n = ' num2str(nIncluded)]);

% Display counts for specific subgroups
categories = {'F-FH+', 'F-FH-', 'M-FH+', 'M-FH-', 'F-FH+/-', 'M-FH+/-'};
for i = 1:length(categories)
    count = height(subjectInfo(subjectInfo.sexFHSUD == categories{i}, :));
    disp(['Total ' categories{i} ' subjects: n = ' num2str(count)]);
end
%% Save Updated Subject Information
savedir = '/Users/louisaschilling/Desktop/FINAL CODE PAPER/Data';
cd(savedir);
save('subjectInfo_SUDcohort.mat', 'subjectInfo', 'nIncluded');
disp('subjectInfo saved.');
