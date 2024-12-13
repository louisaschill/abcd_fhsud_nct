# abcd_fhsud_nct
Code for paper: Sex-specific differences in brain activity dynamics of youth with a family history of substance use disorder

**Data extraction and pre-processing**
1. Extract non-imaging (demographic, SU, fam hx, etc) information: *abcd_get_nonimaging.m*
2. Exclude subjects based on exclusion criteria: *abcd_exclusions.m*
3. Extract imaging info and pre-process (normalize by mean GM signal and remove outlier frames): *abcd_get_imaging.m*

**K-means clustering**
1. Concatenate timeseries across all subjects for clustering: *abcd_concTS.m*
2. Repeat k-means across full range of k's to identify optimal k range: *abcd_repkmeans.m*
3. Elbow plot for optimal k identification: 
4. Get final clustering from partition with most adjusted mutual information across runs: *abcd_kmeans_ami.m*
5. Name clusters by cosine similarity and radar plots (Fig. 1):  

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
5. Follow-up SU and behavioral/environmental correlations
6. Clusters by MRI model
7. By type of FH SUD
8. FO/AR/DT/TP by group
9. By puberty status 

**Helper functions**
- run_kmeans.m
- GET_CENTROIDS.m
- NAME_CLUSTERS_ANGLE.m
- extractInstument.m
- runANCOVA.m
