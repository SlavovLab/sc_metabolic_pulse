#### Required libraries
library(QuantQC)
library(dplyr)
library(ggpubr)
library(Seurat)


#### Set paths
path_raw <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/raw/01_Regular_searches/'
path_meta <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'



#### paths for relevant raw and meta data
data_path <- paste0(path_raw,'Rep1_5day_male/Report.tsv')

linker_path <- paste0(path_meta,"01_Sample_metadata/Rep1_5day_male/5day_male_linker.csv")


one <- paste0(path_meta,"01_Sample_metadata/Rep1_5day_male/5day_male_old_isolated.xls")
two <- paste0(path_meta,"01_Sample_metadata/Rep1_5day_male/5day_male_young_isolated.xls")
all_cells <- list(Old = one,
                  Young = two)



#### Initialize QQC object and convert report table to peptide matrix 
#### with linked meta data
r1_5day_male <- DIANN_to_QQC(data_path,linker_path, plex = 2, carrier = F)

r1_5day_male <- Miceotope_cellXpeptide(r1_5day_male,chQVal = 1,t = 5)

r1_5day_male <- link_cellenONE_Raw(r1_5day_male,all_cells)



#### Some QC checks for LC batch effects and nPOP slide viz (optional)

#r1_5day_male <- Calculate_run_order_statistics(r1_5day_male)

# Plot MS Intensity Statistics
#PlotIntensityDrift(r1_5day_male)

#Plot Retention time Statistics
#PlotRTDrift(r1_5day_male)

#PlotSlideLayout_celltype(r1_5day_male)
#PlotSlideLayout_label(r1_5day_male)



#### Evaluating single cell quality compared to negative controls
#### and filtering out of failed cells

r1_5day_male <- EvaluateNegativeControls(r1_5day_male)

PlotNegCtrl(r1_5day_male)

# log10 intensity as failed cell filter
r1_5day_male <- FilterBadCells(r1_5day_male, min_intens = 7.0)



#### Remove extra peptides that are lowly abundant for proteins that already have
#### 5 good peptides 

# for abundance
nrow(r1_5day_male@matricies@peptide)
r1_5day_male <- Trim_extra_peptides(r1_5day_male)
nrow(r1_5day_male@matricies@peptide)

# for heavy and light peptides
nrow(r1_5day_male@miceotopes@Raw_H)
r1_5day_male <- Trim_extra_peptides_miceotopes(r1_5day_male)
nrow(r1_5day_male@miceotopes@Raw_H)


#### Here we save the raw peptide quant table for later modeling (04_Missing_data_proc.R)
write.csv(r1_5day_male@matricies@peptide,paste0(path_meta,'02_raw_reptide_X_singleCell/r1_peptide.csv'))



#### Plot cell size vs MS intensity
PlotCellSizeVsIntensity(r1_5day_male, type = 'sample')



#### Normalize peptide data and collapse to protein level
r1_5day_male@ms_type <- 'miceotopes'
r1_5day_male <- CollapseToProtein(r1_5day_male, 1,LC_correct = T)

# Viz for LC related batch effects
PlotLC_Deviations(r1_5day_male,global_trend = F)

################################################ 
#### Below is code for supplemental figure pannels on LC batch correction procedure
################################################ 
plot_LC_correct <- r1_5day_male@LC_batch_deviations[[2]] %>% filter(score == '0.547435706446752')

ggplot(plot_LC_correct, aes(x = order)) +
  geom_point(aes(y = data), colour = "black", size = 2, alpha = 0.7) +
  geom_line(aes(y = model), colour = "red", linewidth = 1) +
  labs(
    x = "Order",
    y = "Intensity",
    title = "Observed vs model-predicted LC signal",
    colour = NULL
  ) +
  theme_classic(base_size = 14)

ggplot(plot_LC_correct, aes(x = order)) +
  geom_point(aes(y = data-model), colour = "black", size = 2, alpha = 0.7) +
  labs(
    x = "Order",
    y = "Corrected intensity",
    title = "",
    colour = NULL
  ) +
  theme_classic(base_size = 14)


ggplot() +
  geom_histogram(
    data = data.frame(dev = r1_5day_male@LC_batch_deviations[[1]]),
    aes(x = dev),
    bins = 30,               # adjust as needed
    fill = "steelblue",
    colour = "white",
    alpha = 0.85
  ) + geom_vline(xintercept = .1,color = 'red')+
  labs(
    title = "LC batch 1 deviations",
    x = "Rsq of spline fit",
    y = "# Peptides"
  ) +
  theme_classic(base_size = 14)



#### Plot number of proteins peptides and protein level data missingness 
r1_5day_male@ms_type <- 'DDA'
PlotProtAndPep(r1_5day_male)
PlotDataComplete(r1_5day_male)



#### Compute and plot correlations between relative peptide levels for peptides
#### mapping to the same protein
r1_5day_male <- SharedPeptideCor(r1_5day_male)

PlotPepCor(r1_5day_male)
median(r1_5day_male@pep.cor[[1]]$Cor)

#PlotMS1vMS2(r1_5day_male)


#### KNN protein level imputation
r1_5day_male <- KNN_impute(r1_5day_male)


#### Protein level batch correction for mTRAQ label bias
r1_5day_male <- BatchCorrect(r1_5day_male,run = F,labels = T)


#### Compute PCA and overlay some cell type markers and potential artifacts
r1_5day_male <- ComputePCA(r1_5day_male, impute = T)

FeaturePCA(r1_5day_male, prot = 'Q922U2', imputed = F) #bas
FeaturePCA(r1_5day_male, prot = 'Q99NB1', imputed = F) #sec
FeaturePCA(r1_5day_male, prot = 'P27573', imputed = F) #cil

FeaturePCA(r1_5day_male, prot = 'Q61233', imputed = F) # Immune
FeaturePCA(r1_5day_male, prot = 'O55226', imputed = F) # Chond
FeaturePCA(r1_5day_male, prot = 'P20152', imputed = F) # Fib 

# plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(r1_5day_male, by = "Run order")
PlotPCA(r1_5day_male, by = "Total protein")
PlotPCA(r1_5day_male, by = "Label")


#### Compute UMAP and overlay some cell type markers and potential artifacts
r1_5day_male <- ComputeUMAP(r1_5day_male)

PlotUMAP(r1_5day_male)
PlotUMAP(r1_5day_male, by = 'Total protein')
PlotUMAP(r1_5day_male, by = 'Run order')
PlotUMAP(r1_5day_male, by = 'Label')
PlotUMAP(r1_5day_male, by = 'Condition')


#r1_5day_male <- Miceotope_protein_collapse(r1_5day_male)
save(r1_5day_male, file = paste0(path_meta,"03_QuantQC_objects/r1_5day_male.RData"))

