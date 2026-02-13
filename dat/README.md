# Data folder contents

Download the data from [Zenodo (DOI: 10.5281/zenodo.14902834)](https://zenodo.org/records/14902834). Unzip and add the contents to this folder.

Below is a description of each subfolder and file, including the units and normalization applied.

---

## 01_Sample_metadata/

Sample-level metadata for each of the 4 preparations. Used during preprocessing (Script 01) to link MS runs to cell identities.

| Subfolder | File | Description |
|-----------|------|-------------|
| Rep1_5day_male/ | `5day_male_linker.csv` | Maps MS run files to mTRAQ label and cell identity |
| | `5day_male_old_isolated.xls` | cellenONE isolation log for 24-month-old cells (cell diameter in micrometers) |
| | `5day_male_young_isolated.xls` | cellenONE isolation log for 2-month-old cells (cell diameter in micrometers) |
| Rep2_5day_female/ | `5day_female_linker.csv` | Maps MS run files to mTRAQ label and cell identity |
| | `5day_female_old_isolated.xls` | cellenONE isolation log for 24-month-old cells |
| | `5day_female_young_isolated.xls` | cellenONE isolation log for 2-month-old cells |
| Rep3_10day_male/ | `10day_male_linker.csv` | Maps MS run files to mTRAQ label and cell identity |
| | `10day_male_young_isolated.xls` | cellenONE isolation log for 2-month-old cells |
| Rep4_10day_female/ | `10day_female_linker.csv` | Maps MS run files to mTRAQ label and cell identity |
| | `10day_young_isolated.xls` | cellenONE isolation log for 2-month-old cells |

---

## 02_raw_peptide_X_singleCell/

Raw peptide-level quantification matrices (peptide x single cell). Produced by Script 01.

| File | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `r1_peptide.csv` | Peptide intensities, prep 1 (5-day male) | | |
| `r2_peptide.csv` | Peptide intensities, prep 2 (5-day female) | | |
| `r3_peptide.csv` | Peptide intensities, prep 3 (10-day male) | | |
| `r4_peptide.csv` | Peptide intensities, prep 4 (10-day female) | | |

- **Rows:** Peptide sequence-charge identifiers
- **Columns:** Single cell identifiers

---

## 03_QuantQC_objects/

QuantQC R objects containing fully processed single-cell proteomics data. Produced by Scripts 01 and updated by Script 03 (recycling correction).

| File | Description |
|------|-------------|
| `r1_5day_male.RData` | QuantQC object, prep 1 (5-day male) |
| `r2_5day_female.RData` | QuantQC object, prep 2 (5-day female) |
| `r3_10day_male.RData` | QuantQC object, prep 3 (10-day male) |
| `r4_10day_female.RData` | QuantQC object, prep 4 (10-day female) |

Each object contains:

| Slot | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `@matricies@peptide` | Peptide-level abundance matrix | | |
| `@matricies@protein` | Protein-level abundance matrix (no imputation) | | |
| `@matricies@protein.imputed` | Protein-level abundance matrix (KNN imputed) | | |
| `@matricies@protein_abs` | Absolute protein concentrations | | |
| `@miceotopes@Raw_H` | Heavy isotope peptide intensities | | |
| `@miceotopes@Raw_L` | Light isotope peptide intensities | | |
| `@miceotopes@Alpha_pep` | Recycling-corrected clearance rates (peptide level) | | |
| `@miceotopes@Alpha_prot` | Recycling-corrected clearance rates (protein level) | | |
| `@meta.data` | Cell metadata (ID, Order, label, Run, prot_total, diameter) | | |

---

## 04_Gene_X_SingleCell_and_annotations/

Integrated gene-by-single-cell matrices and cell annotations across all 4 preparations. Produced by Scripts 02 and 03.

### Annotation files

| File | Description |
|------|-------------|
| `sc_protein_annotations.csv` | Per-cell metadata for proteomics data |
| `sc_mRNA_annotations.csv` | Per-cell metadata for mRNA data |
| `MS_Run_annotations.csv` | MS run-level annotations |

**sc_protein_annotations.csv columns:**

| Column | Description | Units |
|--------|-------------|-------|
| `umap_1`, `umap_2` | UMAP coordinates | |
| `sample` | Preparation number (1-4) | |
| `cluster` | Louvain cluster assignment | |
| `Cell_Type` | Annotated cell type (Basal, Secretory, Ciliated, Chondrocyte, Fibroblast, Smooth muscle, Immune) | |
| `age` | Mouse age (Old = 24 months, Young = 2 months) | |
| `prot_total` | Total MS protein intensity per cell | |
| `diameter` | Cell diameter from cellenONE imaging | |
| `ID` | Unique cell identifier (cellID_prepN) | |
| `label` | mTRAQ label (0 or 4) | |
| `Run` | MS run identifier | |
| `Order` | Injection order within the MS run | |
| `prep` | Preparation name (one, two, three, four) | |
| `prot_total_adj` | Median-centered total protein intensity per prep | |

### Data matrices (gene x single cell)

| File | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `sc_protein_relative.csv` | Relative protein abundance | | |
| `sc_protein_absolute.csv` | Absolute protein concentrations | | |
| `clearance_relative.csv` | Relative protein clearance rates | | |
| `clearance_absolute.csv` | Absolute protein clearance rates | | |
| `Missingness_results.csv` | Missing data model output (predicted biological missingness per protein per cell type) | | |

- **Rows:** UniProt protein accession IDs
- **Columns:** Single cell identifiers (cellID_prepN)

---

## 05_Stan_model_input_data/

Formatted input data for the Bayesian regulation model (Script 07). Produced by Script 04.

| File | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `protein_input_stan.csv` | Protein abundance observations for Stan | | |
| `rna_input_stan.csv` | mRNA abundance observations for Stan | | |
| `clearance_stan.csv` | Protein clearance rate observations for Stan | | |

**protein_input_stan.csv columns:**

| Column | Description |
|--------|-------------|
| `Peptide` | Peptide sequence-charge identifier |
| `cell_type` | Cell type |
| `Protein` | UniProt accession |
| `split_gene` | Gene symbol |
| `Data_set` | Preparation (prep1-prep4) |
| `regressed_p_obs` | Residual peptide abundance after spline correction for run order |
| `n_obs` | Number of non-missing observations (used as precision weight) |
| `gene_id_pro` | Integer gene index for Stan |
| `celltype_id_pro` | Integer cell type index for Stan |

---

## 06_Gene_X_CellType/

Gene-by-cell-type summary matrices for each modality. Produced by Scripts 05, 06, and 07.

### Relative_abundance/

Relative abundances across cell types (fold changes relative to the average across cell types).

**Frequentist estimates (Script 05):**

| File | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `protein_freq.csv` | Relative protein abundance | | |
| `mRNA_freq.csv` | Relative mRNA abundance | | |
| `clearance_freq.csv` | Relative protein clearance rates | | |
| `translation_freq.csv` | Relative translation efficiency (computed as: protein + clearance - mRNA) | | |

**Bayesian posterior estimates (Script 07):**

| File | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `protein_bayes.csv` | Posterior mean protein abundance | | |
| `protein_bayes_95CI.csv` | 95% credible interval width for protein | | |
| `mRNA_bayes.csv` | Posterior mean mRNA abundance | | |
| `mRNA_bayes_95CI.csv` | 95% credible interval width for mRNA | | |
| `clearance_bayes.csv` | Posterior mean clearance rate | | |
| `clearance_bayes_95CI.csv` | 95% credible interval width for clearance | | |
| `translation_bayes.csv` | Posterior mean translation efficiency | | |
| `translation_bayes_95CI.csv` | 95% credible interval width for translation | | |

- **Rows:** Gene symbols
- **Columns:** Cell types (Basal, Secretory, Ciliated, Chondrocyte, Fibroblast, Smooth muscle, Immune)

### Absolute_abundance/

Absolute abundances per cell type. Produced by Script 06.

| File | Description | Units | Normalization |
|------|-------------|-------|---------------|
| `abundance.csv` | Absolute protein concentration per cell type | | |
| `mRNA.csv` | Absolute mRNA abundance per cell type | | |
| `clearance.csv` | Absolute protein clearance rate per cell type | | |
| `translation.csv` | Absolute translation efficiency per cell type | | |

- **Rows:** Gene symbols
- **Columns:** Cell types

---

## 07_Output_tables/

Additional output tables from downstream analysis.

---

## Reference files

| File | Description |
|------|-------------|
| `Mouse.fasta` | Mouse proteome FASTA (used for UniProt accession to gene symbol mapping) |
| `GO_Human.txt` | Gene Ontology annotation file |

---

## Supplemental tables

| File | Description |
|------|-------------|
| `Supplemental_table1.xlsx` | Supplemental Table 1 from the manuscript |
| `Supplemental_table2.csv` | GO enrichment results for regulation-classified genes (Script 08) |
| `Supplemental_table3.csv` | Single-cell covariation analysis results (Script 08) |
