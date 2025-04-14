# Drug Sensitivity Scoring Pipeline Based on PPI Topology and Expression Change

This repository contains the core analytical scripts used in our study titled:  
**"[Title of the Manuscript]"**, submitted to *npj Precision Oncology*.

It integrates transcriptomic data, protein-protein interaction (PPI) network topology, and ChEMBL-based drug-target information to prioritize therapeutic compounds in CDK4/6 inhibitor-resistant cancers.

---

## 🔍 Overview of Workflow

1. **Differential Gene Expression**  
   Using either DESeq2 or limma depending on dataset format.

2. **PPI Network Construction**  
   Filtered high-confidence STRING interactions mapped to DEGs.

3. **Network Topology Analysis**  
   Calculating centrality scores (Degree, Betweenness, Eigenvector, MCC, EPC) and computing PCA-based weights.

4. **Target Integration**  
   Mapping genes to UniProt, merging logFC and centrality, computing per-target importance.

5. **Drug Sensitivity Scoring**  
   Matching with ChEMBL drug-targets and scoring drug efficacy using expression-weighted topology.

6. **Statistical Validation of Scoring System**  
   Binomial test and graphical visualization for evaluating prediction accuracy across datasets.

---

## 📁 Repository Contents

| File Name | Description |
|-----------|-------------|
| `01_DESeq2_DifferentialExpression.R` | Differential gene expression analysis using DESeq2 or limma based on `counts.csv` and `metadata.csv` |
| `02_PPI_NetworkConstruction.R` | Construct STRING-based protein-protein interaction networks |
| `03_PPI_Topology_Analysis_and_Weighting.R` | Analyze network topology and calculate centrality-based scores |
| `04_DEG_Mapping_and_Uniprot_Merge.R` | Merge DEG results with UniProt mapping and interaction partners |
| `05_Drug_Target_Matching_and_Scoring.R` | Match with ChEMBL targets and compute weighted drug scores |
| `06_Validation_Pvalue_ScoringSystem.R` | Validate prediction accuracy using binomial test and rank-based p-values (input: `Summary of Datasets accuracy .csv`) |
| `counts.csv` / `metadata.csv` | Input expression matrix and sample metadata |
| `GSE200029_parental_vs_Resistant_filteredDEGs.xlsx` | Processed differential gene expression results (example dataset) |
| `gse20029_target_importance.csv` | Target-level topology-integrated importance scores |
| `drug_scores_gse20029.csv` | Final ranked drug sensitivity scores |

---

## 🚀 Reproducibility Instructions

```r
# Clone the repository
git clone https://github.com/qiankey759/drug_sensitivity_scoring.git

# Run scripts in sequence
source("scripts/01_DESeq2_DifferentialExpression.R")
source("scripts/02_PPI_NetworkConstruction.R")
source("scripts/03_PPI_Topology_Analysis_and_Weighting.R")
source("scripts/04_DEG_Mapping_and_Uniprot_Merge.R")
source("scripts/05_Drug_Target_Matching_and_Scoring.R")
source("scripts/06_Validation_Pvalue_ScoringSystem.R")
```

All scripts are self-contained and annotated. Output files are saved in-place under `results/` and `figures/`.

---

## 🔐 Data Availability

Due to the size and proprietary nature of the original datasets, only representative files are provided in this repository.  
Full data and intermediate files are available upon reasonable request from the corresponding author.

---

## 📩 Contact

For academic correspondence or access to full datasets and drug scoring results, please contact:

**Dr. Qian Keyang**  
Email: *qiankey759@163.com*  
Affiliation: Wuxi Medical College of Jiangnan University
