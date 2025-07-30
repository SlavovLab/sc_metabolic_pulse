#### Required libraries
library(QuantQC)
library(dplyr)
library(ggpubr)
library(Seurat)


#### Set paths
path_raw <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/raw/01_Regular_searches/'
path_meta <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'



#### paths for relevant raw and meta data
data_path <- paste0(path_raw,'Rep2_5day_female/')

linker_path <- paste0(path_meta,"Rep2_5day_female/5day_female_linker.csv")


one <- paste0(path_meta,"01_Sample_metadata/Rep2_5day_female/5day_female_old_isolated.xls")
two <- paste0(path_meta,"01_Sample_metadata/Rep2_5day_female/5day_female_young_isolated.xls")
all_cells <- list(Old = one,
                  Young = two)


#### Initialize QQC object and convert report table to peptide matrix 
#### with linked meta data
r2_5day_female <- DIANN_to_QQC(data_path,linker_path, plex = 2, carrier = F)

r2_5day_female <- Miceotope_cellXpeptide(r2_5day_female,chQVal = 1,t = 5)

r2_5day_female <- link_cellenONE_Raw(r2_5day_female,all_cells)



#### Some QC checks for LC batch effects and nPOP slide viz (optional)

#r2_5day_female <- Calculate_run_order_statistics(r2_5day_female)

# Plot MS Intensity Statistics

#PlotIntensityDrift(r2_5day_female)

#Plot Retention time Statistics
#PlotRTDrift(r2_5day_female)

#PlotSlideLayout_celltype(r2_5day_female)
#PlotSlideLayout_label(r2_5day_female)


#### Evaluating single cell quality compared to negative controls
#### and filtering out of failed cells

r2_5day_female <- EvaluateNegativeControls(r2_5day_female)

PlotNegCtrl(r2_5day_female)

# log10 intensity as failed cell filter
r2_5day_female <- FilterBadCells(r2_5day_female, min_intens = 7.0)


#### Remove extra peptides that are lowly abundant for proteins that already have
#### 5 good peptides 

nrow(r2_5day_female@matricies@peptide)
r2_5day_female <- Trim_extra_peptides(r2_5day_female)
nrow(r2_5day_female@matricies@peptide)

dim(r2_5day_female@miceotopes@Raw_H)
r2_5day_female <- Trim_extra_peptides_miceotopes(r2_5day_female)
dim(r2_5day_female@miceotopes@Raw_H)


#### Here we save the raw peptide quant table for later modeling (04_Missing_data_proc.R)
write.csv(r2_5day_female@matricies@peptide,paste0(path_meta,'02_raw_reptide_X_singleCell/r2_peptide.csv'))


#### Plot cell size vs MS intensity
PlotCellSizeVsIntensity(r2_5day_female, type = 'sample')


#### Normalize peptide data and collapse to protein level
r2_5day_female@ms_type <- 'miceotopes'
r2_5day_female <- CollapseToProtein(r2_5day_female, 1, LC_correct = T)

# Viz for LC related batch effects
#PlotLC_Deviations(r2_5day_female,global_trend = T)



#### Plot number of proteins peptides and protein level data missingness 
r2_5day_female@ms_type <- 'DDA'
PlotProtAndPep(r2_5day_female)
PlotDataComplete(r2_5day_female)


#### Compute and plot correlations between relative peptide levels for peptides
#### mapping to the same protein

r2_5day_female <- SharedPeptideCor(r2_5day_female)

PlotPepCor(r2_5day_female)

#PlotMS1vMS2(r2_5day_female)


#### KNN protein level imputation
r2_5day_female <- KNN_impute(r2_5day_female)


#### Protein level batch correction for mTRAQ label bias
r2_5day_female <- BatchCorrect(r2_5day_female,run = F,labels = T)

#### Compute PCA and overlay some cell type markers and potential artifacts
r2_5day_female <- ComputePCA(r2_5day_female, impute = T)


FeaturePCA(r2_5day_female, prot = 'Q922U2', imputed = F) #bas
FeaturePCA(r2_5day_female, prot = 'Q99NB1', imputed = F) #sec
FeaturePCA(r2_5day_female, prot = 'P68369', imputed = F) #cil


FeaturePCA(r2_5day_female, prot = 'Q61233', imputed = F) # Immune
FeaturePCA(r2_5day_female, prot = 'O55226', imputed = F) # Chond
FeaturePCA(r2_5day_female, prot = 'P20152', imputed = F) # Fib 

# plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(r2_5day_female, by = "Run order")
PlotPCA(r2_5day_female, by = "Total protein")
PlotPCA(r2_5day_female, by = "Label")


#### Compute UMAP and overlay some cell type markers and potential artifacts
r2_5day_female <- ComputeUMAP(r2_5day_female)

PlotUMAP(r2_5day_female)
PlotUMAP(r2_5day_female, by = 'Total protein')
PlotUMAP(r2_5day_female, by = 'Run order')
PlotUMAP(r2_5day_female, by = 'Label')
PlotUMAP(r2_5day_female, by = 'Condition')


#r2_5day_female <- Miceotope_protein_collapse(r2_5day_female)
save(r2_5day_female, file = paste0(path_meta,"03_QuantQC_objects/r2_5day_female.RData"))



