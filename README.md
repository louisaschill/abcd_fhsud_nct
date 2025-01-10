# abcd_fhsud_nct
Code for paper: Sex-specific differences in brain activity dynamics of youth with a family history of substance use disorder

**Data extraction and pre-processing**
1. Extract non-imaging (demographic, SU, fam hx, etc) information: *abcd_get_nonimaging.m*
2. Exclude subjects based on exclusion criteria: *abcd_exclusions.m*
3. Extract imaging info and pre-process (normalize by mean GM signal and remove outlier frames): *abcd_get_imaging.m*
4. Concatenate timeseries across all subjects for clustering: *abcd_concTS.m*


***k*-means clustering**
1. Repeat *k*-means across full range of k's to identify optimal *k* range: *abcd_repkmeans.m*
2. Elbow plot for optimal k identification: 
3. Get final clustering from partition with most adjusted mutual information across runs: *abcd_kmeans_ami.m*
4. Name clusters by cosine similarity and radar plots (Fig. 1):  

**Network control theory**
1. Calculate subject specific centroids:
2. Calculate optimal T (T with highest negative correlation between transition energy and probability): 
3. Calculate transition energies (can be done with group avg or individual SC): 

**Analyses**
1. Global TE Fig. 2: *abcd_globalTE.m*
2. Network TE Fig. 3: *abcd_networkTE.m*
3. Regional TE Fig. 4: *abcd_regionalTE.m*

**Supplementary analyses**
1. NCANDA analysis
2. Within-site
3. Within-income
4. Within-MRI model
5. k = 5
6. Individual SC 

**Reviewer analyses**
1. Follow-up SU and behavioral/environmental correlations
2. Clusters by MRI model
3. By type of FH SUD
4. FO/AR/DT/TP by group
5. By puberty status 

**Helper functions**
- run_kmeans.m
- GET_CENTROIDS.m
- NAME_CLUSTERS_ANGLE.m
- extractInstument.m
- runANCOVA.m
