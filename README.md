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

- **​​Expression Matrix​​:** A raw protein/gene expression matrix with features (genes/proteins) as rows and samples as columns 
- **​​Cell Type Metadata​​:** A data frame containing two required columns:  
sample: Sample identifiers matching the column names of the expression matrix  
celltype: Corresponding cell type classifications for each sample  
### 2.2 Calculates FJSD scores for all cell types (4 cell types)
```
FJSD_4 <- FJSD(bulk_protein, cell_type_4)
```
### 2.3 Calculates optimal FJSD score cutoff by maximizing cluster separation (measured by silhouette coefficient)
```
FJSD_4_Cutoff <- Optimal_FJSD_Cutoff(bulk_protein, cell_type_4, FJSD_4)
```
The optimal FJSD score cutoff of FJSD_4
```
FJSD_4_Cutoff$Optimal_FJSD
FJSD_cutoff  SC_mean
0.8          0.177725
```
```
FJSD_4_Cutoff$Plot
```
<img width="894" height="529" alt="截屏2025-08-06 16 45 27" src="https://github.com/user-attachments/assets/cdcc63a9-a655-48f6-87ed-7f92dbe81a40" />
### 2.4 Volcano Plot Visualization for FJSD Results
```
FJSD_4_Top <- Top_FJSD(FJSD_4$FJSD_list, Top = Inf, FJSD_cutoff = FJSD_4_Cutoff$Optimal_FJSD$FJSD_cutoff, p_val_cutoff = 0.05, p_adj_cutoff = 0.05)
FJSD_Epi <- FJSD_geom_point(FJSD_4$FJSD_list$Epithelial_cell, FJSD_Cutoff = FJSD_4_Cutoff$Optimal_FJSD$FJSD_cutoff, title = "Epithelial_cell Master TFs Analysis")
FJSD_Epi
```
<img width="742" height="588" alt="截屏2025-08-06 16 46 31" src="https://github.com/user-attachments/assets/ded08614-89db-402a-bf01-cfc216676f2f" />
### 2.5 Visualize Gene Expression Patterns Across Samples
```
FJSD_geom_line <- FJSD_geom_line(bulk_protein, FJSD_4_Top$Epithelial_cell$gene[1:5], 
                                 cell_type = cell_type_4, title = "Gene Expression")
FJSD_geom_line
```
<img width="1070" height="533" alt="截屏2025-08-06 16 46 00" src="https://github.com/user-attachments/assets/e9766e97-f60e-43be-b9af-9211730bcbb5" />















