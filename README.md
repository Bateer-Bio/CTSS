# FJSD: Fold Change and Jensen-Shannon Divergence for Identifying Master Transcription Factors with Specific and High Expression
## 1.	Introduction to FJSD
### 1.1	Overview of FJSD
The FJSD (Fold Change and Jensen-Shannon Divergence) method is designed to identify master transcription factors (TFs) that are highly expressed and specific to different cell types in mouse lung data. This method combines both the fold change (FC) and Jensen-Shannon divergence (JSD) scores to select TFs with high and specific expression levels.
<img width="1974" height="929" alt="overview" src="https://github.com/user-attachments/assets/171b59f2-9e69-4958-b933-38acafcde8c1" />
### 1.2 Installing R package
```
library(devtools)
devtools::install_github("Bateer-Bio/FJSD")
```
## 2.	Main functions of FJSD
### 2.1 Load data
```
library(FJSD)
library(dplyr)
bulk_protein <- read.csv("./data/proteomeTF_data.csv",row.names = 1)
cell_type <- read.csv("./data/cell_type_sorted_proteome.csv")
cell_type_4 <- data.frame(sample = cell_type$SampleID, celltype = cell_type$cell_class)
```
**Input Requirements:​​**

- **​​Expression Matrix​​:** A raw protein/gene expression matrix with features (genes/proteins) as rows and samples as columns (rownames should be gene/protein identifiers)
- **​​Cell Type Metadata​​:** A data frame containing two required columns:  
sample: Sample identifiers matching the column names of the expression matrix  
celltype: Corresponding cell type classifications for each sample  
