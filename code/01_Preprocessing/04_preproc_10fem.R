#### Required libraries
library(QuantQC)
library(dplyr)
library(ggpubr)
library(Seurat)


#### Set paths
path_raw <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/raw/01_Regular_searches/'
path_meta <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


#### paths for relevant raw and meta data
data_path <- paste0(path_raw,'Rep4_10day_female/Report.tsv')

linker_path <- paste0(path_meta,'01_Sample_metadata/Rep4_10day_female/10day_female_linker.csv')

one <- paste0(path_meta,'01_Sample_metadata/Rep4_10day_female/10day_young_isolated.xls')
all_cells <- list(Young = one)


#### Initialize QQC object and convert report table to peptide matrix 
#### with linked meta data
r4_10day_female <- DIANN_to_QQC(data_path,linker_path, plex = 2, carrier = F)

r4_10day_female <- Miceotope_cellXpeptide(r4_10day_female,chQVal = 1,t = 10)

r4_10day_female <- link_cellenONE_Raw(r4_10day_female,all_cells)



#### Some QC checks for LC batch effects and nPOP slide viz (optional)

#r4_10day_female <- Calculate_run_order_statistics(r4_10day_female)

# Plot MS Intensity Statistics
#PlotIntensityDrift(r4_10day_female)

#Plot Retention time Statistics
#PlotRTDrift(r4_10day_female)

#PlotSlideLayout_celltype(MiceLab_exp3)
#PlotSlideLayout_label(MiceLab_exp3)


#### Evaluating single cell quality compared to negative controls
#### and filtering out of failed cells

r4_10day_female <- EvaluateNegativeControls(r4_10day_female)

PlotNegCtrl(r4_10day_female)

# log10 intensity as failed cell filter
r4_10day_female <- FilterBadCells(r4_10day_female, min_intens = 6.5)


#### Remove extra peptides that are lowly abundant for proteins that already have
#### 5 good peptides 

nrow(r4_10day_female@matricies@peptide)
r4_10day_female <- Trim_extra_peptides(r4_10day_female)
nrow(r4_10day_female@matricies@peptide)

nrow(r4_10day_female@miceotopes@Raw_H)
r4_10day_female <- Trim_extra_peptides_miceotopes(r4_10day_female)
nrow(r4_10day_female@miceotopes@Raw_H)


#### Here we save the raw peptide quant table for later modeling (04_Missing_data_proc.R)
write.csv(r4_10day_female@matricies@peptide,paste0(path_meta,'02_raw_reptide_X_singleCell/r4_peptide.csv'))


#### Plot cell size vs MS intensity
PlotCellSizeVsIntensity(r4_10day_female, type = 'sample')


#### Normalize peptide data and collapse to protein level
r4_10day_female@ms_type <- 'miceotopes'
r4_10day_female <- CollapseToProtein(r4_10day_female, 1, LC_correct = T)

# Viz for LC related batch effects
#PlotLC_Deviations(r1_5day_male,global_trend = T)


#### Plot number of proteins peptides and protein level data missingness 
r4_10day_female@ms_type <- 'DDA'
PlotProtAndPep(r4_10day_female)
PlotDataComplete(r4_10day_female)



#### Compute and plot correlations between relative peptide levels for peptides
#### mapping to the same protein
r4_10day_female <- SharedPeptideCor(r4_10day_female)

PlotPepCor(r4_10day_female)
median(r4_10day_female@pep.cor[[1]]$Cor)

#PlotMS1vMS2(r4_10day_female)

#### KNN protein level imputation
r4_10day_female <- KNN_impute(r4_10day_female)


#### Protein level batch correction for mTRAQ label bias
r4_10day_female <- BatchCorrect(r4_10day_female,run = F,labels = T)


#### Compute PCA and overlay some cell type markers and potential artifacts
r4_10day_female <- ComputePCA(r4_10day_female, impute = T)


FeaturePCA(r4_10day_female, prot = 'Q922U2', imputed = F) #bas
FeaturePCA(r4_10day_female, prot = 'Q99NB1', imputed = F) #sec
FeaturePCA(r4_10day_female, prot = 'P27573', imputed = F) #cil

FeaturePCA(r4_10day_female, prot = 'Q61233', imputed = F) # Immune
FeaturePCA(r4_10day_female, prot = 'O55226', imputed = F) # Chond
FeaturePCA(r4_10day_female, prot = 'P20152', imputed = F) # Fib P63158

PlotPCA(r4_10day_female, by = "Run order")
PlotPCA(r4_10day_female, by = "Total protein")
PlotPCA(r4_10day_female, by = "Label")


#### Compute UMAP and overlay some cell type markers and potential artifacts
r4_10day_female <- ComputeUMAP(r4_10day_female)

PlotUMAP(r4_10day_female)
PlotUMAP(r4_10day_female, by = 'Total protein')
PlotUMAP(r4_10day_female, by = 'Run order')
PlotUMAP(r4_10day_female, by = 'Label')
PlotUMAP(r4_10day_female, by = 'Condition')


#r4_10day_female <- Miceotope_protein_collapse(r4_10day_female)
save(r4_10day_female, file = paste0(path_meta,"03_QuantQC_objects/r4_10day_female.RData"))

