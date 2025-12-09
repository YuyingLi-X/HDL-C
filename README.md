# High-definition likelihood inference of genetic colocalization (HDL-C)

## Introduction
HDL-C is a likelihood-based framework for genetic colocalization analysis that extends the high-definition likelihood local genetic correlation method [HDL-L](https://www.nature.com/articles/s41588-025-02123-3). HDL-C performs inference on whether two traits share sufficiently strong local genetic effects to warrant evidence of colocalization. 

Unlike traditional Bayesian colocalization approaches that rely on restrictive assumptions (e.g., a single causal variant per trait), HDL-C leverages the multivariate Gaussian structure of GWAS z-scores under high-dimensional LD and formulates a constrained likelihood ratio test to evaluate whether the local genetic correlation exceeds a biologically meaningful threshold.

## System Requirements
HDL-C requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements

**OS Requirements**
This package is supported for macOS and Linux. 

**R Requirements**
R Dependencies: 
```R
dplyr
data.table
tidyverse
```

## Installation 
HDL-C can be easily installed from GitHub using the `remotes` package. If you don't already have `remotes` installed, the following commands will manage the installation for you:
```R
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("YuyingLi-X/HDL-C")
library(HDLC)  #Notice: The name of the R package without `-`
```
This installation process only takes a few minutes.


## Tutorial
Here, we provide a step-by-step tutorial for HDL-C and a real data example at the end. Before you begin your analysis, ensure that you have the necessary resources downloaded, so you need to download the reference panel and LD at first: 

### Step 1: Reference panel and local region definition
As with HDL-L, we have already prepared the pre-computed reference panel and LD for each protein in UKB-PPP. You can download it from [Zenodo](https://zenodo.org/records/14209926).

In the "LD.path", it includes LD files, eigenvectors, and eigen matrices for all local regions, ending with "_LDSVD.rda"

In the "bim.path", it includes bim files for local regions, which helps to clean the summary statistics data and check if there are multiallelic or duplicated SNPs

If you want to build a reference panel for your predefined loci, you can use compute_ld_eigen.R

### Step 2: The format of summary statistics
To analyze your data using HDL-L, it is crucial to format your summary statistics correctly. Below are the required columns that your input data file must include:

- `SNP`: SNP ID  
- `A1`: Effect allele  
- `A2`: Reference allele  
- `N`: Sample size  
- `Z`: Z-score  

If `Z` is not available, alternatively, you may provide:  
- `b`: Estimate of marginal effect in GWAS  
- `se`: Standard error of the estimates of marginal effects in GWAS  

If the GWAS is based on logistic regression, `b` should be the logarithm of OR (odds ratio), and `se` is the standard error of log(OR). 

The summary statistics should look like (b and se can be absent in this example since Z is available):

```R
##          SNP A1 A2      N        b       se      Z
## 1  rs3131962  G  A 205475 0.001004 0.004590 0.2187
## 2 rs12562034  A  G 205475 0.005382 0.005011 1.0740
## 3 rs11240779  A  G 205475 0.002259 0.003691 0.6119
## 4 rs57181708  G  A 205475 0.005401 0.005114 1.0562
## 5  rs4422948  G  A 205475 0.005368 0.003604 1.4893
## 6  rs4970383  A  C 205475 0.004685 0.003582 1.3080
```

You can use `HDL.data.wrangling.R` in HDL git repository to do data wrangling for data from [the Neale Lab round 2 GWAS of UK Biobank](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291) 

You can follow the [instruction](https://github.com/zhenin/HDL/wiki/Format-of-summary-statistics) to format raw GWAS summary statistics into HDL and HDL-C input.

For example 

```R
system("Rscript /Path/to/HDL.data.wrangling.R gwas.file=/Path/to/I25.gwas.imputed_v3.male.tsv.bgz LD.path=/Path/to/UKB_imputed_SVD_eigen99_extraction GWAS.type=UKB.Neale output.file=/Path/to/output/file/I25 log.file=/Path/to/output")

```

### Step 3: Running HDL.C
**For one protein**
```r
# Load required libraries
library(HDLC)
library(dplyr)
library(data.table)

# Load example data in the package
data(gwas1.df)  # For "Protein CELSR2"
data(gwas2.df)  # For "Trait B"

# Paths to the LD reference and bim files
LD.path  <- "/Path/to/LD/LD.path/"
bim.path <- "/Path/to/bim/bim.path/"
chr = 1
proname = "CELSR2"
res <- HDL.C(
  gwas1      = gwas1.df,     # Replace with your loaded data
  gwas2      = gwas2.df,     # Replace with your loaded data
  Trait1name = "CELSR2",
  Trait2name = "I25",
  LD.path    = LD.path,
  bim.path   = bim.path,
  chr        = chr,
  proname      = proname,
  N0          = 0
)
print(res)
```
The output of the HDL.C is:

```R
Analysis starts on Fri May 23 04:06:28 2025 
0 SNPs were removed in GWAS 1 due to missing N or missing test statistic.  
0 SNPs were removed in GWAS 2 due to missing N or missing test statistic.  
689 out of 689 (100%) SNPs in reference panel are available in GWAS 1.  
689 out of 689 (100%) SNPs in reference panel are available in GWAS 2.  


Estimates: 
Heritability of phenotype 1: 0.1167, P = 2.34e-38 
Heritability of phenotype 2: 6e-04, P = 4e-04 
Genetic Correlation: -0.9364, 95% confidence interval (-1, -0.7126), P = 5.86e-09 
Evidence for colocalization: *****

Codes:
      |rG| 0 ' ' 0.5 '+' 0.7 '++' 0.9 '+++' 1 
P    
1
' '          ' '     ' '     ' '      ' '  
0.05
'+'          ' '     '*'     '**'     '***' 
0.005
'++'         ' '     '**'    '***'    '****' 
0.0005
'+++'        ' '     '***'   '****'   '*****' 
0 

Analysis finished at Fri May 23 04:06:30 2025 
```

The level of evidence for colocalization is reported on a five-star scale ranging from "None" (no evidence) to "*****" (very strong evidence), providing an interpretable summary of the statistical and genetic support for shared causal variants.

The “prolist.rda” file includes all protein names and chromosome annotations.

The HDL.C function takes several arguments as outlined below:

1. **gwas1.df**
The first formatted summary statistics data. The input data frame should include the following columns: 
- `SNP`: SNP ID
- `A1`: Effect allele
- `A2`: Reference allele
- `N`: Sample size
- `Z`: Z-score
Alternatively, if `Z` is not provided, you may include:
- `b`: Estimate of marginal effect in GWAS
- `se`: Standard error of the estimates of marginal effects in GWAS.

2. **gwas2.df**
For the second formatted summary statistics data, the columns should be the same as `gwas1.df`.

3. **Trait1name**
The trait name for **gwas1.df**.

4. **Trait2name**
The trait name for **gwas2.df**.

5. **LD.path**
Path to the directory where the decompressed LD.path.zip file is stored.

6. **bim.path**
Path to the directory where the decompressed bim.path.zip file is stored.

7. **Nref**
Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335272.

8. **N0**
Number of individuals included in both cohorts.

9. **output.file**
Location where the log and results should be written. If you do not specify a file, the log will be printed on the console. The default is NULL.

10. **eigen.cut**
Specifies which eigenvalues and eigenvectors in each LD score matrix should be used for HDL. The default value is 0.99. Users are allowed to specify a numeric value between 0 and 1 for eigen.cut.

11. **intercept.output**
Logical, `FALSE` by default. Determines whether the intercept terms are included in the result `estimates.df` or not.

12. **fill.missing.N**
If `NULL` (default), SNPs with missing `N` are removed. You can specify "median", "min", or "max" so that the missing `N` will be filled accordingly. For example, "median" means the missing `N` are filled with the median `N` of the SNPs with available `N`.

13. **lim**
Tolerance limitation, default `lim = exp(-18)`.

14. **chr**
The chromosome to which the region belongs.

15. **piece**
The piece of the genome to which the region belongs. The whole genome is divided into 2,476 smaller, semi-independent blocks, each defined by LD calculated by Plink. The SNP information in each local region is included in these two data sets: "UKB_snp_counter_imputed.RData" and "UKB_snp_list_imputed_vector.RData".


## Citation
If you use the HDL-C software, please cite:
- Li Y, Zhai R, Yang Z, Li T, Pawitan Y, Shen X (2025). High-definition likelihood inference of colocalization reveals protein biomarkers for human complex diseases
- Li, Y., Pawitan, Y. & Shen, X. An enhanced framework for local genetic correlation analysis. Nat Genet 57, 1053–1058 (2025). https://doi.org/10.1038/s41588-025-02123-3
- Ning, Z., Pawitan, Y. & Shen, X. *High-definition likelihood inference of genetic correlations across human complex traits*. Nat Genet (2020).

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
