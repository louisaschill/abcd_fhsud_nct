# ABCD_FHSUD_NCT

This repository contains the code to reproduce the analysis in Schilling et al. 2024, "Sex-specific differences in brain activity dynamics of youth with a family history of substance use disorder", bioRxiv.

This code is based on previous work by (added here with modifications):
- Singleton et al., *Nature Communications*, 2022 ([Repo](https://github.com/singlesp/energy_landscape))
- Cornblath et al., *Communications Biology*, 2020 ([Repo](https://github.com/ejcorn/brain_states))

---

## Requirements
- **MATLAB R2017a or later**
- Required package:
  - `gramm`

---

## Usage: Replicating the Analysis
The following scripts should be used for each step of the analysis:

### **Data Extraction and Pre-Processing**
1. **Extract non-imaging data** (demographic, substance use, family history, etc.):  
   `abcd_get_nonimaging.m`
2. **Apply exclusion criteria:**  
   `abcd_exclusions.m`
3. **Extract and preprocess imaging data** (normalize by mean GM signal, remove outliers):  
   `abcd_get_imaging.m`
4. **Concatenate timeseries across all subjects for clustering:**  
   `abcd_concTS.m`

### ***k*-means clustering**
1. **Run repeated *k*-means across a range of *k* values** to determine the optimal range:  
   `abcd_repkmeans.m`
2. **Generate an elbow plot** for *k* selection.
3. **Obtain final clustering** using the partition with the highest adjusted mutual information across runs:  
   `abcd_kmeans_ami.m`
4. **Name clusters** using cosine similarity and radar plots (*Fig. 2*).

### **Network Control Theory (NCT) Analysis**
1. **Calculate subject-specific centroids.**
2. **Determine optimal T** (highest negative correlation between transition energy and probability).
3. **Calculate transition energies** (using either group-average or individual structural connectivity).

### **Main Analyses**
1. **Global Transition Energy (TE) Analysis** (*Fig. 3*):  
   `abcd_globalTE.m`
2. **Network TE Analysis** (*Fig. 4*):  
   `abcd_networkTE.m`
3. **Regional TE Analysis** (*Fig. 5*):  
   `abcd_regionalTE.m`

### **Supplementary Analyses**
1. **NCANDA analysis**
2. **Within-site effects**
3. **Within-income effects**
4. **Within-MRI model effects**
5. **Analysis with *k* = 5**
6. **Individual SC analysis**

### **Review Analyses**
1. **Follow-up substance use and behavioral/environmental correlations**
2. **Clusters by MRI model**
3. **Effects by type of family history of SUD**
4. **FO/AR/DT/TP group-wise analysis**
5. **Analysis by puberty status**

---

## Helper Functions
- `ejc_bs_code_mod.m` - Modified code from EJ Cornblath by Parker Singleton ([Source](https://github.com/singlesp/energy_landscape))
- `extractInstrument.m` - Extracts instruments from ABCD non-imaging data
- `plotMatricesDiff.m` - Plots differences between two matrices
- `mergeTables.m` - Merges tables based on subject keys
- `runANCOVA.m` - Runs ANCOVAs on univariate or multivariate data
- `extract_ancova_results.m` - Extracts specific results from ANCOVA outputs

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
