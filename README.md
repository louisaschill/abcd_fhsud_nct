# Sex-specific differences in brain activity dynamics of youth with a family history of substance use disorder

This repository contains the code to reproduce the analysis in Schilling et al. 2024, "Sex-specific differences in brain activity dynamics of youth with a family history of substance use disorder", bioRxiv.

This repo relies on that of previous work (adapted version of relevant code has been added to this repo): 
- Singleton et al., *Nature Communications*, 2022 ([Repo](https://github.com/singlesp/energy_landscape))
- Cornblath et al., *Communications Biology*, 2020 ([Repo](https://github.com/ejcorn/brain_states))

---

## Requirements
- **MATLAB R2017a or later**
- Required package:
  - `gramm` plotting toolbox (https://github.com/piermorel/gramm)
- Hardware: No non-standard hardware required; analyses were run locally on a MacBook Pro 2022.

## Installation
- Clone or download this repository.
- Make sure you have MATLAB R2017a or later.
- Install the gramm toolbox if not already installed.
- Add all scripts and data paths to your MATLAB path.
- Typical install time: Less than 5 minutes on a standard desktop or laptop.
  
---

## Usage: Replicating the Analysis

For ~6000 5 min fMRI scans and 86 ROIs, I was able to perform all of this analysis locally on a 2022 MacBook Pro with a total run time of <1 day, with the most time needed by k-means clustering. 

The following scripts should be used for each step of the analysis:

### **Data Extraction and Pre-Processing**
1. **Extract non-imaging data** (demographic, substance use, family history, etc.): `abcd_get_nonimaging.m`
2. **Apply exclusion criteria:** `abcd_exclusions.m`
3. **Extract and preprocess imaging data** (normalize by mean GM signal, remove outliers): `abcd_get_imaging.m`
4. **Concatenate timeseries across all subjects for clustering:** `abcd_concTS.m`
5. **Get mean framewise displacement for all subjects:** `abcd_get_meanFD.m`

### ***k*-means clustering**
1. **Run repeated *k*-means across a range of *k* values** to determine the optimal range: `abcd_repkmeans.m` and `run_kmeans.m`
2. **Generate an elbow plot** for *k* selection: `abcd_elbow.m`
3. **Obtain final clustering** using the partition with the highest adjusted mutual information across runs: `abcd_kmeans_ami.m`
4. **Name clusters** using cosine similarity and radar plots (*Fig. 2*): `abcd_radar_plots.m`
   
### **Network Control Theory (NCT) Analysis**
1. **Calculate subject-specific centroids.**: `abcd_subj_centroids.m` and `abcd_generate_subcentroids.m`
2. **Determine optimal T** (highest negative correlation between transition energy and probability): `abcd_tsweep.m`
3. **Calculate transition energies** (using either group-average or individual structural connectivity): `abcd_calculate_TE.m`

### **Main Analyses**
1. **Global Transition Energy (TE) Analysis** (*Fig. 3*): `abcd_globalTE.m`
2. **Network TE Analysis** (*Fig. 4*): `abcd_networkTE.m`
3. **Regional TE Analysis** (*Fig. 5*): `abcd_regionalTE.m`

### **Supplementary Analyses**
1. **NCANDA analysis**: `ncanda_supp.m` 
2. **Within-site**: `abcd_supp_singlesite.m`
3. **Within-income**: `abcd_supp_income.m`
4. **Within-MRI model**: `abcd_supp_model.m`
5. **Analysis with *k* = 5**: re-run clustering and above analysis with *k*=5
6. **Individual SC analysis**: `abcd_supp_indivSC.m`
7. **Follow-up substance use and behavioral/environmental correlations**: `abcd_supp_followup_SU.m` and `abcd_cbcl_supp.m`
8. **Clusters by MRI model**: `abcd_mri_centroids.m`  
9. **Fractional occupancy group-wise analysis**: `fractional_occupancy.m`

---

## Raw Data

### **ABCD Data**
Processed ABCD neuroimaging data used in this study were uploaded by Ooi et al. (2022) and are available via NDA Study 1368 (https://nda.nih.gov/study.html?id=1368). Researchers with access to the ABCD data will be able to download the data: https://nda.nih.gov/study.html?id=1368. The ABCD imaging data used in this report came from
https://doi.org/10.15154/1504041 and non-imaging data was from the 5.1 release (http://dx.doi.org/10.15154/z563-zd24). These data were used in the analyses described in https://doi.org/10.15154/c9z7-ng36. The ABCD data are publicly available via the NIMH Data Archive (NDA).

### **NCANDA Data**
Researchers with access to NCANDA data can download it via https://nda.nih.gov/study.html?id=4513. 

---

## Contact
For any questions regarding this code, please contact:

**Louisa Schilling**  
Email: **lss4002@med.cornell.edu**

---
