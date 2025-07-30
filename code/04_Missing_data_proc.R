#### Initial model
library(ggplot2)
library(QuantQC)
library(dplyr)
library(reshape2)
library(stringr)
library(tidyr)
library(purrr)
library(data.table)
library(doParallel)
library(foreach)
library(splines)
library(pROC)
library(gridExtra)
library(seqinr)






##################
# data path
##################

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


#########################
# Functions for missingness classifier
#########################

technical_effect_model_vec <- function(data_conc, data_miss) {

  # Initialize predicted columns
  data_conc$predicted_p_obs <- NA_real_
  data_miss$predicted_p_obs <- NA_real_

  # Group and nest data_conc
  nested_conc <- data_conc %>%
    group_by(Peptide, Data_set, Label) %>%
    nest()

  # Group and nest data_miss
  nested_miss <- data_miss %>%
    group_by(Peptide, Data_set, Label) %>%
    nest()

  # Combine into a single table so that each group from conc can
  # be matched with the corresponding group in miss
  combined_nested <- full_join(
    nested_conc, nested_miss,
    by = c("Peptide", "Data_set", "Label"),
    suffix = c("_conc", "_miss")
  )

  # For each group (row), build a spline using data_conc and then
  # predict for both data_conc and data_miss
  results <- combined_nested %>%
    rowwise() %>%
    mutate(
      result = list({

        # Extract the data frames for this group
        d_conc <- data_conc
        d_miss <- data_miss

        # Only build a spline if we have enough data in conc
        if (!is.null(d_conc) && nrow(d_conc) > 5) {

          # Rough logic to set df
          set_df <- floor(nrow(d_conc) / 11)
          if (set_df > 15) set_df <- 15
          if (set_df < 2)  set_df <- 2

          # Fit the smooth spline on (Order, p_obs)
          fit <- stats::smooth.spline(d_conc$Order, d_conc$p_obs, df = set_df)

          # Predict for both data_conc and data_miss
          d_conc$predicted_p_obs <- predict(fit, d_conc$Order)$y

          if (!is.null(d_miss) && nrow(d_miss) > 0) {
            d_miss$predicted_p_obs <- predict(fit, d_miss$Order)$y
          }
        }

        # Return both updated data frames in a list
        list(d_conc = d_conc, d_miss = d_miss)
      })
    ) %>%
    ungroup()

  # Separate the nested predictions back into two data frames
  results_expanded <- results %>%
    mutate(
      d_conc = map(result, "d_conc"),
      d_miss = map(result, "d_miss")
    ) %>%
    dplyr::select(-result)

  # Unnest and restore columns for data_conc
  final_conc <- results_expanded %>%
    dplyr::select(Peptide, Data_set, Label, d_conc) %>%
    unnest(cols = d_conc)

  # Unnest and restore columns for data_miss
  final_miss <- results_expanded %>%
    dplyr::select(Peptide, Data_set, Label, d_miss) %>%
    unnest(cols = d_miss)

  # Ungroup to remove grouping structure
  final_conc <- ungroup(final_conc)
  final_miss <- ungroup(final_miss)

  # Return list with updated data frames
  list(final_conc, final_miss)
}

Same_set_dt <- function(data_miss) {
  # Convert to data.table (if not already)
  setDT(data_miss)

  # 1) Summarize to find if a group has '0'=1 or '4'=1
  group_info <- data_miss[
    ,
    .(
      has_label0_1 = any(Label == "0" & p_obs_present == 1),
      has_label4_1 = any(Label == "4" & p_obs_present == 1)
    ),
    by = .(Peptide, Data_set, Order)
  ]

  # 2) Key both tables for fast joins
  setkey(data_miss, Peptide, Data_set, Order)
  setkey(group_info, Peptide, Data_set, Order)

  # 3) Update 'Same_Set' in data_miss by joining in group_info
  #    fcase() lets you efficiently handle multiple conditions in a single pass
  data_miss[group_info,
            Same_Set := fcase(
              i.has_label0_1 & Label == "4", 1,
              i.has_label4_1 & Label == "0", 1,
              default = NA_real_
            )
  ]

  # Return updated data
  return(as.data.frame(data_miss))
}

Eff_missingness_model_parallel_statusbar <- function(data_miss) {
  library(doSNOW)
  library(foreach)
  library(dplyr)

  # Set up parallel backend using available cores (reserve one core)
  cores <- parallel::detectCores() - 1
  cl <- makeCluster(cores)
  registerDoSNOW(cl)

  # Get unique peptides and progress bar settings
  unique_peptides <- unique(data_miss$Peptide)
  n_peptides <- length(unique_peptides)
  pb <- txtProgressBar(max = n_peptides, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Initialize prediction columns
  data_miss$predicted_p_obs_all  <- NA_real_
  data_miss$predicted_p_obs_tech <- NA_real_

  # Process each peptide in parallel with progress bar
  results <- foreach(pep = unique_peptides, .packages = c("dplyr"), .options.snow = opts) %dopar% {
    # Subset data for the current peptide
    data_miss_pep <- data_miss %>% filter(Peptide == pep)

    # If the factor 'Same_Set' or 'Data_set' has only one level, skip modeling
    if (length(unique(data_miss_pep$Same_Set)) == 1 ||
        length(unique(data_miss_pep$Data_set)) == 1) {
      return(list(
        indices = which(data_miss$Peptide == pep),
        pred_all = rep(NA_real_, nrow(data_miss_pep)),
        pred_tech = rep(NA_real_, nrow(data_miss_pep))
      ))
    }

    # Fit full model (including cell_type)
    missingness_model <- glm(
      p_obs_present ~ Size + Same_Set + predicted_p_obs + Data_set + cell_type,
      family = binomial(link = "logit"),
      data = data_miss_pep
    )

    pred_all_pep <- predict(missingness_model, type = "response")

    # Fit partial model (without cell_type) to isolate its effect
    model_without_cell_type <- glm(
      p_obs_present ~ Size + Same_Set + predicted_p_obs + Data_set,
      family = binomial(link = "logit"),
      data = data_miss_pep
    )

    pred_tech_pep <- predict(model_without_cell_type, type = "response")

    # Return the row indices and predictions for this peptide subset
    list(
      indices = which(data_miss$Peptide == pep),
      pred_all = pred_all_pep,
      pred_tech = pred_tech_pep
    )
  }

  close(pb)
  stopCluster(cl)

  # Update the original data frame with the predictions from each peptide
  for (res in results) {
    data_miss$predicted_p_obs_all[res$indices]  <- res$pred_all
    data_miss$predicted_p_obs_tech[res$indices] <- res$pred_tech
  }

  return(data_miss)
}

Proc_fasta <- function(path){
  convert_mouse <- read.fasta(path,set.attributes = T,whole.header = T)
  convert_mouse <- names(convert_mouse)
  parse_row<-grep("GN=",convert_mouse, fixed=T)
  split_prot<-str_split(convert_mouse[parse_row], pattern = fixed("GN="))
  gene<-unlist(split_prot)[seq(2,2*length(split_prot),2)]
  prot <- unlist(split_prot)[seq(1,2*length(split_prot),2)]
  prot_parse <- grep("|",prot, fixed=T)
  gene_parse <- grep(" ",gene, fixed=T)
  split_gene<-str_split(gene[parse_row], pattern = fixed(" "))
  split_gene<-unlist(split_gene)[seq(1,3*length(split_gene),3)]
  split_prot<-str_split(prot[parse_row], pattern = fixed("|"))
  split_prot<-unlist(split_prot)[seq(2,3*length(split_prot),3)]
  convert_mouse  <- as.data.frame(cbind(split_prot,split_gene))
  
  return(convert_mouse)
}

#########################
# Functions for missingness classifier
#########################


# Load in the data (QuantQC objects) from each prep
load(paste0(path_dat,'03_QuantQC_objects/r1_5day_male.RData'))
load(paste0(path_dat,'03_QuantQC_objects/r2_5day_female.RData'))
load(paste0(path_dat,'03_QuantQC_objects/r3_10day_male.RData'))
load(paste0(path_dat,'03_QuantQC_objects/r4_10day_female.RData'))


# Load in the raw peptide data from each prep (this gets normalized in the QQC objects so is not stored)
data1_abs <- read.csv(paste0(path_dat,'02_raw_reptide_X_singleCell/r1_peptide.csv'),row.names = 1)
data2_abs <- read.csv(paste0(path_dat,'02_raw_reptide_X_singleCell/r2_peptide.csv'),row.names = 1)
data3_abs <- read.csv(paste0(path_dat,'02_raw_reptide_X_singleCell/r3_peptide.csv'),row.names = 1)
data4_abs <- read.csv(paste0(path_dat,'02_raw_reptide_X_singleCell/r4_peptide.csv'),row.names = 1)
# Append prep number to IDs
colnames(data1_abs) <- paste0(colnames(data1_abs),'_prep1')
colnames(data2_abs) <- paste0(colnames(data2_abs),'_prep2')
colnames(data3_abs) <- paste0(colnames(data3_abs),'_prep3')
colnames(data4_abs) <- paste0(colnames(data4_abs),'_prep4')

# Normalize Raw peptide data (here we are not accounting for LC variation or label bias
  # which we will model subsequently) which are corrected for in the QQC object stored raw peptide data
data1 <- Normalize_reference_vector(data1_abs,log = T)
data2 <- Normalize_reference_vector(data2_abs,log = T)
data3 <- Normalize_reference_vector(data3_abs,log = T)
data4 <- Normalize_reference_vector(data4_abs,log = T)

# Remove peptides with only a few observations and melt data
data1 <- data1[rowSums(is.na(data1)==F) > 10,]
data1 <- reshape2::melt(data1)
data2 <- data2[rowSums(is.na(data2)==F) > 10,]
data2 <- reshape2::melt(data2)
data3 <- data3[rowSums(is.na(data3)==F) > 10,]
data3 <- reshape2::melt(data3)
data4 <- data4[rowSums(is.na(data4)==F) > 10,]
data4 <- reshape2::melt(data4)


# Read in meta data with cell type assignments
meta_data <- read.csv(paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),row.names = 1)

# Read in additional metadata needed about batch effect information
meta_data_mini1 <- r1_5day_male@meta.data %>% dplyr::select('ID','Order','label','Run')
meta_data_mini1$ID <- paste0(meta_data_mini1$ID,'_prep1')
meta_data_mini2 <- r2_5day_female@meta.data %>% dplyr::select('ID','Order','label','Run')
meta_data_mini2$ID <- paste0(meta_data_mini2$ID,'_prep2')
meta_data_mini3 <- r3_10day_male@meta.data %>% dplyr::select('ID','Order','label','Run')
meta_data_mini3$ID <- paste0(meta_data_mini3$ID,'_prep3')
meta_data_mini4 <- r4_10day_female@meta.data %>% dplyr::select('ID','Order','label','Run')
meta_data_mini4$ID <- paste0(meta_data_mini4$ID,'_prep4')

# Join all relevant meta data to raw peptide level data
data1 <- data1 %>% left_join(meta_data_mini1, by = c('Var2'='ID'))
data1$Data_set <- 'prep1'
data1 <- data1 %>% left_join(r1_5day_male@matricies@peptide_protein_map, by = c('Var1'='seqcharge'))

data2 <- data2 %>% left_join(meta_data_mini2, by = c('Var2'='ID'))
data2$Data_set <- 'prep2'
data2 <- data2 %>% left_join(r2_5day_female@matricies@peptide_protein_map, by = c('Var1'='seqcharge'))

data3 <- data3 %>% left_join(meta_data_mini3, by = c('Var2'='ID'))
data3$Data_set <- 'prep3'
data3 <- data3 %>% left_join(r3_10day_male@matricies@peptide_protein_map, by = c('Var1'='seqcharge'))

data4 <- data4 %>% left_join(meta_data_mini4, by = c('Var2'='ID'))
data4$Data_set <- 'prep4'
data4 <- data4 %>% left_join(r4_10day_female@matricies@peptide_protein_map, by = c('Var1'='seqcharge'))

meta_data$Order <- NULL
meta_data$label <- NULL
meta_data$Run <- NULL


data <- rbind(data1,data2,data3,data4)
data <- data %>% left_join(meta_data, by = c('Var2'='ID'))

## Read in deg rate data

data_alpha1 <- Normalize_reference_vector(r1_5day_male@miceotopes@Alpha_pep,log = T)
data_alpha2 <- Normalize_reference_vector(r2_5day_female@miceotopes@Alpha_pep,log = T)
data_alpha3 <- Normalize_reference_vector(r3_10day_male@miceotopes@Alpha_pep,log = T)
data_alpha4 <- Normalize_reference_vector(r4_10day_female@miceotopes@Alpha_pep,log = T)

colnames(data_alpha1) <- paste0(colnames(data_alpha1),'_prep1')
colnames(data_alpha2) <- paste0(colnames(data_alpha2),'_prep2')
colnames(data_alpha3) <- paste0(colnames(data_alpha3),'_prep3')
colnames(data_alpha4) <- paste0(colnames(data_alpha4),'_prep4')

data_alpha1 <- data_alpha1[rowSums(is.na(data_alpha1)==F) > 10,]
data_alpha1 <- reshape2::melt(data_alpha1)
data_alpha2 <- data_alpha2[rowSums(is.na(data_alpha2)==F) > 10,]
data_alpha2 <- reshape2::melt(data_alpha2)
data_alpha3 <- data_alpha3[rowSums(is.na(data_alpha3)==F) > 10,]
data_alpha3 <- reshape2::melt(data_alpha3)
data_alpha4 <- data_alpha4[rowSums(is.na(data_alpha4)==F) > 10,]
data_alpha4 <- reshape2::melt(data_alpha4)

data_alpha1$Data_set <- 'prep1'
data_alpha1 <- data_alpha1 %>% left_join(r1_5day_male@miceotopes@peptide_protein_map, by = c('Var1'='seqcharge'))

data_alpha2$Data_set <- 'prep2'
data_alpha2 <- data_alpha2 %>% left_join(r2_5day_female@miceotopes@peptide_protein_map, by = c('Var1'='seqcharge'))

data_alpha3$Data_set <- 'prep3'
data_alpha3 <- data_alpha3 %>% left_join(r3_10day_male@miceotopes@peptide_protein_map, by = c('Var1'='seqcharge'))

data_alpha4$Data_set <- 'prep4'
data_alpha4 <- data_alpha4 %>% left_join(r4_10day_female@miceotopes@peptide_protein_map, by = c('Var1'='seqcharge'))


data_alpha <- rbind(data_alpha1,data_alpha2,data_alpha3,data_alpha4)
data_alpha <- data_alpha %>% left_join(meta_data, by = c('Var2'='ID'))






# First pass model should have no NAs

model_data <- data



# Fit the abundance model with fixed effects for run order and label, and random effects for cell type and protein



  data <- model_data %>% dplyr::select('Var1','Var2','value','Cell_Type','Protein','prot_total','diameter','Order','label','Data_set')
  colnames(data) <- c('Peptide','Run','p_obs','cell_type','Protein','Size','Diameter','Order','Label','Data_set')

  data <- data %>% filter(cell_type != "Shwann")

  data_conc <- data %>% filter(is.na(p_obs) == F)
  data_conc$Peptide <- as.factor(data_conc$Peptide)
  data_conc$Label <- as.factor(data_conc$Label)
  data_conc$cell_type <- as.factor(data_conc$cell_type)

  data_miss <- data %>%
    mutate(p_obs_present = ifelse(is.na(p_obs), 0, 1))
  
  

  ## Running technical effects model 

  out_dat2 <- technical_effect_model_vec(data_conc,data_miss)
  data_conc <- out_dat2[[1]]
  data_miss <- out_dat2[[2]]

  print('Regressing out techincal effects')

  data_conc <- data_conc %>% filter(is.na(predicted_p_obs) == F)
  data_miss <- data_miss %>% filter(is.na(predicted_p_obs) == F)


  data_miss <- Same_set_dt(data_miss)
  data_miss$Same_Set[is.na(data_miss$Same_Set)] <- 0

  data_conc$regressed_p_obs <- data_conc$p_obs - data_conc$predicted_p_obs


  #### Plotting feature and how they relate to missing data
  
    data_miss_bas <- data_miss %>% filter( Data_set == 'prep1')#cell_type == 'Basal' &
    data_miss_bas_size <- data_miss_bas %>% group_by(Run) %>% summarise(Size = mean(Size,na.rm=T),
                                                                        Diameter = mean(Diameter,na.rm=T),
                                                                        p_obs_present = mean(p_obs_present,na.rm=T))
    
    ggplot(data_miss_bas_size,aes(x = Size,y = 100-100*p_obs_present)) +
      geom_point() + theme_classic(base_size = 18) + ylab('% missing peptide intensities per cell')+xlab('Total cell MS intensity, log2')
    ggplot(data_miss_bas_size,aes(x = Diameter,y = 100-100*p_obs_present)) + 
      geom_point() + theme_classic(base_size = 18) + ylab('% missing peptide intensities')
    
    #ggplot(data_miss_bas_size,aes(x = Diameter,y = Size)) +
    #  geom_point() + theme_classic(base_size = 18) 
    
    cor(data_miss_bas_size$Size,data_miss_bas_size$p_obs_present,method = 'pearson')^2
    cor(data_miss_bas_size$Diameter,data_miss_bas_size$p_obs_present,method = 'pearson')^2
    
    data_miss_bas <- data_miss_bas[sample(nrow(data_miss_bas),5000),]
    ggplot(data_miss_bas,aes(x = as.character(p_obs_present),y = predicted_p_obs-median(predicted_p_obs,na.rm=T))) +
      geom_hline(yintercept = 0)+
      geom_boxplot() + coord_cartesian(ylim = c(-1.5,1.5))+theme_classic(base_size = 18) +
      ylab('Predicted intensity (Retention time, label)') + xlab('') + theme_bw()
    
  
  
  ### Running effective missing data model (takes a second)
  data_miss <- Eff_missingness_model_parallel_statusbar(data_miss)

  #### Find cell type with lowest average expression each protein
  
  top_ct_per_protein <- data_miss_dat %>%                     # original data
    group_by(Protein, cell_type) %>%                           # 1. group by Prot+Cell
    summarise(mean_intensity = mean(p_obs, na.rm = TRUE),  # 2. avg intensity
              .groups = "drop_last") %>%                       # keep Protein grouping
    slice_min(mean_intensity, n = 1, with_ties = FALSE) %>%    # 3. pick top Cell_type
    ungroup() 
  
  top_ct_per_protein$prot_cell <- paste0(top_ct_per_protein$Protein,top_ct_per_protein$cell_type)
  
  
  ### Plot ROC curve for results
  data_miss_plot <- data_miss #$cell_type
  data_miss_plot$Cell_prot <- paste0(data_miss_plot$Protein,data_miss_plot$cell_type)
  data_miss_plot <- data_miss_plot %>% filter(Cell_prot %in% top_ct_per_protein$prot_cell)
  data_miss_plot <- data_miss_plot %>% filter(cell_type == 'Smooth muscle')
  
  roc_obj <- roc(response = data_miss_plot$p_obs_present,
                 predictor = data_miss_plot$predicted_p_obs_all,
                 ci = TRUE)                    # 95% CI by default
  auc_val <- auc(roc_obj)                     # numeric AUC
  print(auc_val)
  
  roc_obj_tech <- roc(response = data_miss_plot$p_obs_present,
                 predictor = data_miss_plot$predicted_p_obs_tech,
                 ci = TRUE)  
  
  auc_val <- auc(roc_obj_tech)                     # numeric AUC
  print(auc_val)
  
  
  
  roc_obj <- data.frame(sensitivities = roc_obj$sensitivities,specificities = roc_obj$specificities)
  roc_obj$type = 'Cell type'
  roc_obj_tech <- data.frame(sensitivities = roc_obj_tech$sensitivities,specificities = roc_obj_tech$specificities)
  roc_obj_tech$type = 'Technical only'
  roc_obj <- rbind(roc_obj,roc_obj_tech)
  
  roc_obj <- roc_obj[sample(nrow(roc_obj),500),]
  
  ggplot(roc_obj,aes(x = specificities,y = sensitivities,color = type)) + geom_line(size = 1.5)+ theme_classic(base_size = 18) +
    scale_x_reverse() + scale_color_manual(values = c('black','grey50')) + geom_abline(intercept = 1,slope = 1)
  
 



  ### Continue subsetting data

  data_miss_post <- data_miss %>% filter(p_obs_present == 0)
  
  data_miss_dat <- data_miss %>% filter(p_obs_present == 1)

  data_miss_post <- data_miss_post %>% filter(is.na(predicted_p_obs_all) == F)

  data_miss_post$predicted_p_obs_all <- 1- data_miss_post$predicted_p_obs_all
  data_miss_post$predicted_p_obs_tech <- 1- data_miss_post$predicted_p_obs_tech
  hist(data_miss_post$predicted_p_obs_all)
  hist(data_miss_post$predicted_p_obs_tech)

  data_miss_post$diff <- data_miss_post$predicted_p_obs_all - data_miss_post$predicted_p_obs_tech

  

  
  #### Plot model difference protein specific way
  
  diff_tech <- c()
  diff_Bio <- c()
  prot <- c()
  ct <- c()
  for(i in unique(data_miss_post$cell_type)){
    data_miss_post_hold <- data_miss_post %>% filter(cell_type == i)
    print(i)
    for(j in unique(data_miss$Protein)){
      diff_tech <- c(diff_tech,mean(data_miss_post_hold$predicted_p_obs_tech[data_miss_post_hold$Protein == j ],na.rm=T))
      diff_Bio <- c(diff_Bio,mean(data_miss_post_hold$predicted_p_obs_all[data_miss_post_hold$Protein == j ],na.rm=T))
      
      prot <- c(prot,j)
      ct <- c(ct,i)
    }
  }
  df_learn <- data.frame(prot  = prot,ct = ct, diff_tech = diff_tech,diff_b = diff_Bio)
  df_learn$prot_Cell <- paste0(df_learn$prot,df_learn$ct)
  df_learn <- df_learn %>% filter(prot_Cell %in% top_ct_per_protein$prot_cell)

  df_learn$diff <- df_learn$diff_b -df_learn$diff_tech
  
  df_learn$prot <- factor(df_learn$prot, levels = df_learn$prot[order(df_learn$diff_tech)])

  df_learn <- df_learn %>% filter(ct == 'Smooth muscle')

  df_learn$diff_b <- df_learn$diff_b*100
  df_learn$diff_tech <- df_learn$diff_tech*100
  
  ggplot(df_learn, aes(x = prot)) +
    ## one grey line per protein
    geom_segment(aes(xend = prot, y = diff_tech, yend = diff_b),
                 colour = "grey60", linewidth = 0.4) +
    ## the two points
    geom_point(aes(y = diff_tech), colour = "blue") +
    geom_point(aes(y = diff_b),    colour = "red")  +
    theme_classic(base_size = 18) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    ylab('% of missing values predicted')+ xlab('Proteins')
  
 
   
  
  data_miss_post$diff[data_miss_post$diff < 0] <- NA
  data_miss_post <- data_miss_post %>% filter(is.na(diff) == F)

  data_miss_post_collapse1 <- data_miss_post %>% group_by(cell_type,Run,Protein) %>% summarise(diff=median(diff,na.rm=T))
  data_miss_post_collapse2 <- data_miss_post_collapse1 %>% group_by(cell_type,Protein) %>% summarise(diff=sum(diff,na.rm=T))

  
  write.csv(data_miss_post_collapse2,paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/Missingness_results.csv'))


  # 1) Summarize how many times (diff) to replicate per (cell_type, Peptide)
  data_miss_post_collapse1 <- data_miss_post %>%
    group_by(cell_type, Data_set, Peptide) %>%
    summarise(diff = round(sum(diff, na.rm = TRUE)), .groups = "drop")

  # 2) For each (cell_type, Peptide), find the minimum regressed_p_obs in data_miss_obs
  min_vals <- data_conc %>%
    group_by(cell_type, Data_set,Peptide) %>%
    summarise(min_val = min(regressed_p_obs, na.rm = TRUE), .groups = "drop")

  # 3) Join the min_vals with the diff counts
  min_vals_diff <- min_vals %>%
    left_join(data_miss_post_collapse1, by = c("Data_set","cell_type", "Peptide"))

  # 4) Replicate each min_val 'diff' times
  expanded_min_vals <- min_vals_diff %>%
    filter(!is.na(diff), diff > 0) %>%
    uncount(weights = diff)

  # 5) Bring back the Protein info (if each Peptide maps to a single Protein in each cell type)
  #    We'll join to a distinct list of (cell_type, Peptide, Protein).
  protein_map <- data_conc %>%
    distinct(Data_set,cell_type, Peptide, Protein)

  expanded_min_vals <- expanded_min_vals %>%
    left_join(protein_map, by = c("Data_set","cell_type", "Peptide")) %>%
    mutate(regressed_p_obs = min_val) %>%
    dplyr::select(Data_set,cell_type, Protein, Peptide, regressed_p_obs)




  # 6) Combine with the original data
  combined_data <- bind_rows(
    data_conc %>% dplyr::select(Data_set,cell_type, Protein, Peptide, regressed_p_obs),
    expanded_min_vals
  )
  
  
  ###### Protein model input
  
  convert = Proc_fasta(paste0(path_dat,'Mouse.fasta'))

  data_protein_model_input <- combined_data %>%
    group_by(cell_type, Protein, Data_set, Peptide) %>%
    summarise( n_obs = sum(!is.na(regressed_p_obs)) ,
              regressed_p_obs = mean(regressed_p_obs, na.rm = TRUE), .groups = "drop")

  data_protein_model_input <- data_protein_model_input %>% left_join(convert, by = c('Protein'='split_prot'))
  data_protein_model_input <- data_protein_model_input %>% filter(is.na(split_gene)==F)
  
  data_protein_model_input$CT_prot <- paste0(data_protein_model_input$cell_type,data_protein_model_input$split_gene)
  
  get_prots_prot <- data_protein_model_input %>%
    group_by(Protein) %>%
    dplyr::summarise(n_cell_types = n_distinct(cell_type))
  
  get_prots_prot <- get_prots_prot %>% filter(n_cell_types > 1)

  nrow(data_protein_model_input)
  data_protein_model_input <- data_protein_model_input %>% filter(Protein %in% get_prots_prot$Protein)
  nrow(data_protein_model_input)
  
  
  ###### Curating deg rates
  
  colnames(data_alpha)[colnames(data_alpha) == 'Var1'] <- 'Peptide'

  data_degradation_model_input <- data_alpha %>%
    group_by(Cell_Type, Protein, Data_set, Peptide) %>%
    summarise( n_obs = sum(!is.na(value)) ,
               value = mean(value, na.rm = TRUE), .groups = "drop")
  
  data_degradation_model_input <- data_degradation_model_input %>% filter(n_obs > 0)
  
  
  data_degradation_model_input <- data_degradation_model_input %>% left_join(convert, by = c('Protein'='split_prot'))
  data_degradation_model_input <- data_degradation_model_input %>% filter(is.na(split_gene)==F)
  
  data_degradation_model_input$CT_prot <- paste0(data_degradation_model_input$Cell_Type,data_degradation_model_input$split_gene)
  
  get_prots_prot <- data_degradation_model_input %>%
    group_by(Protein) %>%
    dplyr::summarise(n_cell_types = n_distinct(Cell_Type))
  get_prots_prot <- get_prots_prot %>% filter(n_cell_types > 1)
  
  nrow(data_degradation_model_input)
  data_degradation_model_input <- data_degradation_model_input %>% filter(Protein %in% get_prots_prot$Protein)
  nrow(data_degradation_model_input)


  ###### Curating RNA rates
  
  rna_seq <- readRDS('/Users/andrewleduc/Desktop/Senescense/seurat_integrated_filered_700_named.rds')
  
  
  sect <- intersect(data_degradation_model_input$split_gene,data_protein_model_input$split_gene)
  sect <- intersect(sect,rownames(rna_seq@assays$RNA@counts))

  Basal = c('Basal 3','Basal 4','Basal 1','Hillock basal')
  
  Secratory = c('Secretory 1 + Goblet','Hillock luminal', 'Secretory 2')
  
  Fibroblast = c('Fibroblast 1', 'Fibroblast 2','Fibroblast 3','Fibroblast 4','Fibroblast 5','Fibroblast 6','Fibroblast 7')
  
  Chondrocyte = c('Chondrocyte 1', 'Chondrocyte 2', 'Chondrocyte 3')
  
  Cilliated = c('Ciliated', 'Ciliated + Sec->Cil')
  
  Smooth_muscle = 'Smooth muscle'
  
  Immune = c('Lymphoid (T and ILCs)','IAMs + IMs + mast','cDC 2 + B cell', 'cDC 1') #
  
  Shwann = c('Schwann')
  
  mRNA_CTs <- list(Basal,Chondrocyte,Secratory,Fibroblast,Immune,Cilliated,Smooth_muscle)
  names(mRNA_CTs) <- c('Basal','Chondrocyte','Secratory','Fibroblast','Immune',
                       'Cilliated','Smooth muscle')
  
  #intersect(names(mRNA_CTs),data_degradation_model_input$Cell_Type)
  
  df_sep <- data.frame(id = names(rna_seq@active.ident))
  df_sep$Data_set <- str_replace(df_sep$id, "^(.*?)_(.*?)_.*$", "\\1_\\2")

  mrna_mat_1 <- matrix(data = NA,nrow = nrow(rna_seq@assays$RNA@counts),ncol = length(mRNA_CTs))
  mrna_mat_2 <- matrix(data = NA,nrow = nrow(rna_seq@assays$RNA@counts),ncol = length(mRNA_CTs))
  mrna_mat_3 <- matrix(data = NA,nrow = nrow(rna_seq@assays$RNA@counts),ncol = length(mRNA_CTs))
  mrna_mat_4 <- matrix(data = NA,nrow = nrow(rna_seq@assays$RNA@counts),ncol = length(mRNA_CTs))

  mrna_data <- list(mrna_mat_1,mrna_mat_2,mrna_mat_3,mrna_mat_4)
  names(mrna_data) <- unique(df_sep$Data_set)

  count_mat <-as.matrix(rna_seq@assays$RNA@counts)
  for(i in 1:length(mRNA_CTs)){

    for(j in 1:length(mrna_data)){

      df_sep_hold <- df_sep %>% filter(Data_set == names(mrna_data)[j])

      s1 <- names(rna_seq@active.ident)[rna_seq@active.ident %in% mRNA_CTs[[i]]]

      s1 <- intersect(s1,df_sep_hold$id)

      mrna_data[[j]][,i] = (.1+ rowMeans(count_mat[,s1]))
      

    }

  }

  for(i in 1:length(mrna_data)){
    rownames(mrna_data[[i]]) <- rownames(rna_seq@assays$RNA@counts)
    colnames(mrna_data[[i]]) <- names(mRNA_CTs)
    mrna_data[[i]] <- Normalize_reference_vector(mrna_data[[i]],log = T)
    mrna_data[[i]] <- Normalize_reference_vector_log(mrna_data[[i]])
    mrna_data[[i]] <- mrna_data[[i]][sect,]

  }
  

  plot(mrna_data[[1]][,5],mrna_data[[2]][,5]);abline(a=0,b=1)
  plot(mrna_data[[1]][,5],mrna_data[[3]][,5]);abline(a=0,b=1)
  plot(mrna_data[[1]][,5],mrna_data[[4]][,5]);abline(a=0,b=1)


  mrna_prep1 <- melt(mrna_data[[1]]);colnames(mrna_prep1)[1:2] <- c('gene','Cell_type')
  mrna_prep2 <- melt(mrna_data[[2]]);colnames(mrna_prep2)[1:2] <- c('gene','Cell_type')
  mrna_prep3 <- melt(mrna_data[[3]]);colnames(mrna_prep3)[1:2] <- c('gene','Cell_type')
  mrna_prep4 <- melt(mrna_data[[4]]);colnames(mrna_prep4)[1:2] <- c('gene','Cell_type')

  rna_dat <- rbind(mrna_prep1,mrna_prep2,mrna_prep3,mrna_prep4)

  
  rna_dat$CT_prot <- paste0(rna_dat$Cell_type,rna_dat$gene)
    
  sect_CT_prot <- intersect(intersect(rna_dat$CT_prot,data_protein_model_input$CT_prot),data_degradation_model_input$CT_prot)
  
  rna_dat <- rna_dat %>% filter(CT_prot %in% sect_CT_prot)
  data_protein_model_input <- data_protein_model_input %>% filter(CT_prot %in% sect_CT_prot)
  data_degradation_model_input <- data_degradation_model_input %>% filter(CT_prot %in% sect_CT_prot)

  length(unique(rna_dat$gene))
  length(unique(data_protein_model_input$split_gene))
  length(unique(data_degradation_model_input$split_gene))
  
  
  
  all_genes     <- sort(unique(data_degradation_model_input$split_gene))
  all_celltypes <- sort(unique(data_degradation_model_input$Cell_Type))
  
  # 2. apply the same factor levels in each DF
  
  
  data_protein_model_input <- data_protein_model_input %>%
    mutate(
      split_gene        = factor(split_gene,        levels = all_genes),
      cell_type   = factor(cell_type,   levels = all_celltypes),
      gene_id_pro = as.integer(split_gene),
      celltype_id_pro = as.integer(cell_type)
    )
  data_protein_model_input2 <- data_protein_model_input %>% group_by(Data_set,Peptide) %>%
    mutate(regressed_p_obs = regressed_p_obs-mean(regressed_p_obs,na.rm=T))
  
  write.csv(data_protein_model_input2,paste0(path_dat,'05_Stan_model_input_data/protein_input_stan.csv'))
  
  rna_dat <- rna_dat %>%
    mutate(
      gene          = factor(gene,        levels = all_genes),
      Cell_type     = factor(Cell_type,   levels = all_celltypes),
      gene_id_mrna  = as.integer(gene),
      celltype_id_mrna = as.integer(Cell_type)
    )
  
  
  rna_dat2 <- rna_dat %>% group_by(Cell_type) %>%
    mutate(value = value-mean(value,na.rm=T))
  rna_dat2 <- rna_dat2 %>% group_by(gene) %>%
    mutate(value = value-mean(value,na.rm=T))
  
  write.csv(rna_dat2,paste0(path_dat,'05_Stan_model_input_data/rna_input_stan.csv'))
  
  data_degradation_model_input <- data_degradation_model_input %>%
    mutate(
      split_gene           = factor(split_gene,        levels = all_genes),
      Cell_Type      = factor(Cell_Type,   levels = all_celltypes),
      gene_id_deg    = as.integer(split_gene),
      celltype_id_deg = as.integer(Cell_Type)
    )
  data_degradation_model_input2 <- data_degradation_model_input %>% group_by(Data_set,Peptide) %>%
    mutate(value = value-mean(value,na.rm=T))
  
  
  write.csv(data_degradation_model_input2,paste0(path_dat,'05_Stan_model_input_data/clearance_stan.csv'))
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#########################
  
  
  i='Pgk1'
  data_degradation_plot <- data_degradation_model_input %>% filter(split_gene == i)
  data_degradation_plot <- data_degradation_plot %>% group_by(Data_set,Peptide) %>%
    mutate(value = value-mean(value,na.rm=T))
  data_degradation_plot <- data_degradation_plot %>% filter(n_obs > 2)
  ggplot(data_degradation_plot,aes(x = Cell_Type,y = value)) + geom_point() + ylim(c(-2,2))
  data_degradation_plot <- data_degradation_plot %>% dplyr::select(cell_type,value)
  data_degradation_plot$type = 'Deg'
  
  data_protein_plot <- data_protein_model_input %>% filter(split_gene == i)
  data_protein_plot <- data_protein_plot %>% group_by(Data_set,Peptide) %>%
    mutate(value = regressed_p_obs-mean(regressed_p_obs,na.rm=T))
  data_protein_plot <- data_protein_plot %>% filter(n_obs > 2)
  ggplot(data_protein_plot,aes(x = cell_type,y = regressed_p_obs)) + geom_point() + ylim(c(-2,2))
  data_protein_plot$Cell_Type = data_protein_plot$cell_type
  data_protein_plot <- data_protein_plot %>% dplyr::select(Cell_Type,value)
  data_protein_plot$type = 'Abundance'
  
  
  
  ggplot(data_degradation_plot,aes(x = Cell_Type,y = value)) + geom_point() + ylim(c(-2,2))
  rna_dat_plot <- rna_dat %>% filter(gene == i)
  ggplot(rna_dat_plot,aes(x = Cell_type,y = value)) + geom_point() +  ylim(c(-2,2))
  rna_dat_plot$Cell_Type <- rna_dat_plot$Cell_type
  rna_dat_plot <- rna_dat_plot %>% dplyr::select(Cell_Type,value)
  
  
  
  plot(mRNA_mat[,],protein_mat[,])
  plot(deg_mat[i,],protein_mat[i,])
  plot(trans_mat[i,],protein_mat[i,])
  
  cor(mRNA_mat[i,],protein_mat[i,],use = 'pairwise.complete.obs')
  cor(deg_mat[i,],protein_mat[i,],use = 'pairwise.complete.obs')
  cor(trans_mat[i,],protein_mat[i,],use = 'pairwise.complete.obs')
  
  library(patchwork)
  plot_gene_relationships('Nutf2')
  
  
  
  

###################### For validating with other analysis pipeline

  #Prot_mat_verify
  data_miss_obs_test <-  data_protein_model_input2 %>%
    group_by(cell_type, split_gene) %>%
    summarise(
      regressed_p_obs = sum(regressed_p_obs * n_obs, na.rm = TRUE) / sum(n_obs, na.rm = TRUE),
      .groups = "drop"
    )

  data_miss_obs_test <- dcast(data_miss_obs_test,split_gene ~ cell_type, value.var = 'regressed_p_obs')
  data_miss_obs_test <- data_miss_obs_test %>% filter(is.na(split_gene)==F)
  rownames(data_miss_obs_test) <- data_miss_obs_test$split_gene
  data_miss_obs_test$split_gene <- NULL
  #data_miss_obs_test <- Normalize_reference_vector_log(data_miss_obs_test)
  data_miss_obs_test <- as.matrix(data_miss_obs_test)
  
  #Deg_mat_verify
  data_miss_obs_test_deg <-  data_degradation_model_input %>%
    group_by(Cell_Type, split_gene) %>%
    summarise(
      regressed_p_obs = sum(value * n_obs, na.rm = TRUE) / sum(n_obs, na.rm = TRUE),
      .groups = "drop"
    )
  
  data_miss_obs_test_deg <- dcast(data_miss_obs_test_deg,split_gene ~ Cell_Type, value.var = 'regressed_p_obs')
  data_miss_obs_test_deg <- data_miss_obs_test_deg %>% filter(is.na(split_gene)==F)
  rownames(data_miss_obs_test_deg) <- data_miss_obs_test_deg$split_gene
  data_miss_obs_test_deg$Shwann <- NULL
  data_miss_obs_test_deg$split_gene <- NULL
  data_miss_obs_test_deg <- Normalize_reference_vector_log(data_miss_obs_test_deg)
  data_miss_obs_test_deg <- as.matrix(data_miss_obs_test_deg)
  
  #RNA_mat verify
  rna_dat_test  <-   rna_dat2 %>%
    group_by(Cell_type, gene) %>%
    summarise(
      value = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  rna_dat_test <- dcast(rna_dat_test,gene ~ Cell_type, value.var = 'value')
  rna_dat_test <- rna_dat_test %>% filter(is.na(gene)==F)
  rownames(rna_dat_test) <- rna_dat_test$gene
  rna_dat_test$gene <- NULL
  rna_dat_test[is.na(data_miss_obs_test)] <- NA
  rna_dat_test <- as.matrix(rna_dat_test)

  
  colnames(rna_dat_test) == colnames(data_miss_obs_test)
  rownames(rna_dat_test) == rownames(data_miss_obs_test)
  
  plot(rna_dat_test,data_miss_obs_test)
  abline(a = 0,b=1)
  cor(c(rna_dat_test),c(data_miss_obs_test),use = 'pairwise.complete.obs',method = 'spearman')
  
  data_miss_obs_test_deg[rowSums(is.na(data_miss_obs_test_deg)==F) < 4 ,] <- NA
  deg_mat[cor_store >.3 ,] <- NA
  View(deg_mat)
  cor_store
  cor_mat <- cor(protein_mat,deg_mat,use = 'pairwise.complete.obs',method = 'pearson')
  cor_store <- c()
  slope_st <- c()
  gene <- c()
  for(i in 1:nrow(protein_mat)){
    cor_store <- c(cor_store,cor(protein_mat[i,],deg_mat[i,],use = 'pairwise.complete.obs'))
    if(sum(is.na(deg_mat[i,])==F)>3){
      slope_st <- c(slope_st,TLS(deg_mat[i,],protein_mat[i,])[[1]])
    }else{
      slope_st <- c(slope_st,0)
    }
    
    gene <- c(gene,rownames(protein_mat)[i])
  }
  
 df_make <- data.frame(gene = gene,cor=cor_store,slope = slope_st)
 df_make <- df_make %>% filter(is.na(cor_store)==F)
 df_make <- df_make %>% filter(cor < -.4)
 
 View(df_make)
 
 hist(log2(-df_make$slope))
 
 plot(protein_mat['Psmd14',],deg_mat['Psmd14',])
 abline(a=0,b=-1)
  
  
Heatmap(
  cor_mat,
  name = "corr",
  cluster_rows = F,
  cluster_columns = F,
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # i, j = row/column indices in cor_mat
    # x, y = center coordinates of the cell
    grid.text(
      round(cor_mat[i, j], 2),  # round to 2 decimal places
      x = x, y = y,
      gp = gpar(fontsize = 10)
    )
  }
)
  


plot_gene_relationships <- function(gene_name) {
  # 1) Find the row index for this gene
  if (!gene_name %in% rownames(mRNA_mat)) {
    stop("Gene '", gene_name, "' not found in mRNA_mat rownames.")
  }
  i <- match(gene_name, rownames(mRNA_mat))
  
  # 2) Build a small data.frame for plotting
  df <- data.frame(
    cell_type = colnames(mRNA_mat),
    mRNA      = as.numeric(mRNA_mat[i, ]),
    protein   = as.numeric(protein_mat[i, ]),
    deg       = as.numeric(deg_mat[i, ]),
    trans     = as.numeric(trans_mat[i, ])
  )
  
  # 3) Create the three scatterplots
  p1 <- ggplot(df, aes(x = mRNA, y = protein)) +
    geom_point() +
    theme_classic() +
    labs(
      title = paste0("mRNA vs Protein"),
      x     = "log2 mRNA (fold-change)",
      y     = "log2 Protein (fold-change)"
    )
  
  p2 <- ggplot(df, aes(x = deg, y = protein)) +
    geom_point() + ylab('')+
    theme_classic() +
    labs(
      title = paste0("Degradation vs Protein"),
      x     = "log2 Degradation rate"
    )
  
  p3 <- ggplot(df, aes(x = trans, y = protein)) +
    geom_point() +ylab('')+
    theme_classic() +
    labs(
      title = paste0("Translation vs Protein"),
      x     = "log2 Translation rate"
    )
  
  # 4) Stack vertically
  combined <- p1 + p2 + p3 +
    plot_layout(ncol = 3) &
    plot_annotation(
      title = paste("Gene:", gene_name),
      theme = theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
      )
    )
  
  return(combined)
}
  
  
  
  

  
  
  
  
  
  
  
  
  

