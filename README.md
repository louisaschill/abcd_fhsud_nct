# abcd_fhsud_nct

Code to reproduce analysis in: Schilling et al. 2024 ("Sex-specific differences in brain activity dynamics of youth with a family history of substance use disorder", Bioxriv). 

This code relies on code from repos Singleton et al., Nature Comms, 2022 (https://github.com/singlesp/energy_landscape), which in turn relies on the repo from Cornblath et al. Comms Bio, 2020 (https://github.com/ejcorn/brain_states). Code from these have been added to this repo with some functions slightly modified.

**Requirements** 
MATLAB R2017a or later, with packages: 
- gramm

To replicate the analysis, use the following scripts for each step: 

**Data extraction and pre-processing**
1. Extract non-imaging (demographic, SU, fam hx, etc) information: *abcd_get_nonimaging.m*
2. Exclude subjects based on exclusion criteria: *abcd_exclusions.m*
3. Extract imaging info and pre-process (normalize by mean GM signal and remove outlier frames): *abcd_get_imaging.m*
4. Concatenate timeseries across all subjects for clustering: *abcd_concTS.m*

***k*-means clustering**
1. Repeat *k*-means across full range of k's to identify optimal *k* range: *abcd_repkmeans.m*
2. Elbow plot for optimal k identification: 
3. Get final clustering from partition with most adjusted mutual information across runs: *abcd_kmeans_ami.m*
4. Name clusters by cosine similarity and radar plots (Fig. 2):   

**Network control theory**
1. Calculate subject specific centroids:
2. Calculate optimal T (T with highest negative correlation between transition energy and probability): 
3. Calculate transition energies (can be done with group avg or individual SC): 

**Analyses**
1. Global TE Fig. 3: *abcd_globalTE.m*
2. Network TE Fig. 4: *abcd_networkTE.m*
3. Regional TE Fig. 5: *abcd_regionalTE.m*

**Supplementary analyses**
1. NCANDA analysis
2. Within-site
3. Within-income
4. Within-MRI model
5. k = 5
6. Individual SC 

**Review analyses**
1. Follow-up SU and behavioral/environmental correlations
2. Clusters by MRI model
3. By type of FH SUD
4. FO/AR/DT/TP by group
5. By puberty status 

**Helper functions**
- ejc_bs_code_mod: code from EJ Cornbalth modified by Parker Singleton (source: https://github.com/singlesp/energy_landscape) 
- extractInstument.m: extracts instruments from ABCD non-imaging information 
- plotMatricesDiff.m: plots difference between two matrices 
- mergeTables.m: merges tables based on subjeckeys 
- runANCOVA.m: runs ancovas on univariate or multi-variate data 
- extract_ancova_results.m: gets results from specific variables within ancovas

**Raw Data**

Processed ABCD neuroimaging data by Ooi et al. (2022) used in this study were uploaded by the original authors to the NDA. Researchers with access to the ABCD data will be able to download the data: https://nda.nih.gov/study.html?id=1368. The ABCD imaging data used in this report came from https://doi.org/10.15154/1504041 and non-imaging data was from the 5.1 release (http://dx.doi.org/10.15154/z563-zd24). These data were used in the analyses described in https://doi.org/10.15154/c9z7-ng36. The ABCD data are publicly available via the NIMH Data Archive (NDA). Collection and distribution of the NCANDA data were supported by NIH funding AA021681, AA021690, AA021691, AA021692, AA021695, AA021696, AA021697. Researchers with access to the NCANDA data will be able to download the data via https://nda.nih.gov/study.html?id=4513.

Please contact Louisa Schilling Singleton (lss4002@med.cornell.edu) with any questions regarding this code.
