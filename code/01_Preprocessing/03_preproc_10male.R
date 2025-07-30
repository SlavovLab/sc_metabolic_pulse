#### Required libraries
library(QuantQC)
library(dplyr)
library(ggpubr)
library(Seurat)


#### Set paths
path_raw <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/raw/01_Regular_searches/'
path_meta <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


#### paths for relevant raw and meta data
data_path <-  paste0(path_raw,'Rep3_10day_male/Report.tsv') 

linker_path <- paste0(path_meta,"01_Sample_metadata/Rep3_10day_male/10day_male_linker.csv")

one <- paste0(path_meta,"01_Sample_metadata/Rep3_10day_male/10day_male_young_isolated.xls")
all_cells <- list(Young = one)


#### Initialize QQC object and convert report table to peptide matrix 
#### with linked meta data
r3_10day_male <- DIANN_to_QQC(data_path,linker_path, plex = 2, carrier = F)

r3_10day_male <- Miceotope_cellXpeptide(r3_10day_male,chQVal = 1,t = 10)

r3_10day_male <- link_cellenONE_Raw(r3_10day_male,all_cells)



#### Some QC checks for LC batch effects and nPOP slide viz (optional)

#r3_10day_male <- Calculate_run_order_statistics(r3_10day_male)

# Plot MS Intensity Statistics
#PlotIntensityDrift(r3_10day_male)

#Plot Retention time Statistics
#PlotRTDrift(r3_10day_male)

#PlotSlideLayout_celltype(r3_10day_male)
#PlotSlideLayout_label(r3_10day_male)


#### Evaluating single cell quality compared to negative controls
#### and filtering out of failed cells

r3_10day_male <- EvaluateNegativeControls(r3_10day_male)

PlotNegCtrl(r3_10day_male)

# log10 intensity as failed cell filter
r3_10day_male <- FilterBadCells(r3_10day_male, min_intens = 7.1)



#### Remove extra peptides that are lowly abundant for proteins that already have
#### 5 good peptides 

nrow(r3_10day_male@matricies@peptide)
r3_10day_male <- Trim_extra_peptides(r3_10day_male)
nrow(r3_10day_male@matricies@peptide)

nrow(r3_10day_male@miceotopes@Raw_H)
r3_10day_male <- Trim_extra_peptides_miceotopes(r3_10day_male)
nrow(r3_10day_male@miceotopes@Raw_H)


#### Here we save the raw peptide abundance table for later modeling (04_Missing_data_proc.R)
write.csv(r3_10day_male@matricies@peptide,paste0(path_meta,'02_raw_reptide_X_singleCell/r3_peptide.csv'))


#### Plot cell size vs MS intensity
PlotCellSizeVsIntensity(r3_10day_male, type = 'sample')


#### Normalize peptide data and collapse to protein level
r3_10day_male@ms_type <- 'miceotopes'
r3_10day_male <- CollapseToProtein(r3_10day_male, 1,LC_correct = T)

# Viz for LC related batch effects
#PlotLC_Deviations(r3_10day_male,global_trend = T)

#### Plot number of proteins peptides and protein level data missingness 
r3_10day_male@ms_type <- 'DDA'
PlotProtAndPep(r3_10day_male)
PlotDataComplete(r3_10day_male)

#### Compute and plot correlations between relative peptide levels for peptides
#### mapping to the same protein
r3_10day_male <- SharedPeptideCor(r3_10day_male)
PlotPepCor(r3_10day_male)

#PlotMS1vMS2(r3_10day_male)

#### KNN protein level imputation
r3_10day_male <- KNN_impute(r3_10day_male)


#### Protein level batch correction for mTRAQ label bias
r3_10day_male <- BatchCorrect(r3_10day_male,run = F,labels = T)

#### Compute PCA and overlay some cell type markers and potential artifacts
r3_10day_male <- ComputePCA(r3_10day_male, impute = T)

FeaturePCA(r3_10day_male, prot = 'Q922U2', imputed = F) #bas
FeaturePCA(r3_10day_male, prot = 'Q99NB1', imputed = F) #sec
FeaturePCA(r3_10day_male, prot = 'P27573', imputed = F) #cil

FeaturePCA(r3_10day_male, prot = 'Q61233', imputed = F) # Immune
FeaturePCA(r3_10day_male, prot = 'O55226', imputed = F) # Chond
FeaturePCA(r3_10day_male, prot = 'P20152', imputed = F) # Fib 

PlotPCA(r3_10day_male, by = "Run order")
PlotPCA(r3_10day_male, by = "Total protein")
PlotPCA(r3_10day_male, by = "Label")
PlotPCA(r3_10day_male, by = "Condition")

#### Compute UMAP and overlay some cell type markers and potential artifacts
r3_10day_male <- ComputeUMAP(r3_10day_male)

PlotUMAP(r3_10day_male)
PlotUMAP(r3_10day_male, by = 'Total protein')
PlotUMAP(r3_10day_male, by = 'Run order')
PlotUMAP(r3_10day_male, by = 'Label')
PlotUMAP(r3_10day_male, by = 'Condition')

FeatureUMAP(r3_10day_male, prot = 'Q61233', imputed = F)
FeatureUMAP(r3_10day_male, prot = 'P20152', imputed = F)


#r3_10day_male <- Miceotope_protein_collapse(r3_10day_male)
save(r3_10day_male, file = paste0(path_meta,"03_QuantQC_objects/r3_10day_male.RData"))

