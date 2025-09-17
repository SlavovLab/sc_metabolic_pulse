# Absolute level analysis, Figures 1f and Figure 2
library(ggplot2)
library(dplyr)
library(seqinr)
library(stringr)
library(reshape2)
library(scales)
library(tidyr)
library(forcats)


palette <- c(
  Immune       = "#F90303",
  Basal        = "#B50202",
  Secratory    = "#B606C4",
  Cilliated    = "#9B70F9",
  Fibroblast   = "#2C78FF",
  `Smooth muscle` = "#0498BA",
  Chondrocyte  = "#03C1ED"
)


##################
# data path
##################

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


TLS <- function(vect1,vect2){
  
  vect1[vect1 == -Inf] <- NA
  vect1[vect1 == Inf] <- NA
  
  
  
  vect2[vect2 == -Inf] <- NA
  vect2[vect2 == Inf] <- NA
  
  int_x <- mean(vect1,na.rm=T)
  int_y <- mean(vect2,na.rm=T)
  
  vect1 <- vect1-int_x
  vect2 <- vect2-int_y
  
  mat <- cbind(vect1,vect2)
  
  mat <- mat[complete.cases(mat),]
  
  TLS_mat <- svd(mat)$v
  
  slope <- TLS_mat[1,1]/TLS_mat[2,1]
  
  int <- c(int_x,int_y)
  
  return(list(slope,int))
  
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

Compute_cell_growth_rate <- function(Cell_Type){
  
  Cell_Type_5day <- Cell_Type %>% filter(sample %in% c('1','2'))
  Cell_Type_10day <- Cell_Type %>% filter(sample %in% c('3','4'))
  
  N_0 = sum(is.na(exp(p_all_alpha_abs['Q64524',Cell_Type_5day$ID]*5)-1) ==F)-
    sum(log2(exp(p_all_alpha_abs['Q64524',Cell_Type_5day$ID]*5)-1) > -1,na.rm=T)/2
  
  N_t = sum(is.na(exp(p_all_alpha_abs['Q64524',Cell_Type_5day$ID]*5)-1) ==F)
  
  print((log2(N_t) - log2(N_0))/5)
  
  N_0 = sum(is.na(exp(p_all_alpha_abs['Q64524',Cell_Type_10day$ID]*10)-1) ==F)-
    sum(log2(exp(p_all_alpha_abs['Q64524',Cell_Type_10day$ID]*10)-1) > -1,na.rm=T)/2
  
  N_t = sum(is.na(exp(p_all_alpha_abs['Q64524',Cell_Type_10day$ID]*10)-1) ==F)
  
  print((log2(N_t) - log2(N_0))/10)
  
  
}


#########################
# Read in data for analysis
#########################

rna_seq <- readRDS(paste0(mRNA_raw_path,'seurat_integrated_filered_700_named.rds'))


p_all_abs <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_absolute.csv'),row.names = 1)
p_all_abs <- as.matrix(p_all_abs)
p_all_alpha_abs <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/clearance_absolute.csv'),row.names = 1)
p_all_alpha_abs <- as.matrix(p_all_alpha_abs)

meta_data <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),row.names = 1)
mRNA_meta <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_mRNA_annotations.csv'),row.names = 1)
mRNA_meta <- mRNA_meta %>% filter(cell_type != 'Shwann')


### Compute cell division rates


# Get cell type IDs (two replicates for error bars)
Basal <- meta_data %>% dplyr::filter(Cell_Type == 'Basal')
Chondrocyte <- meta_data %>% filter(Cell_Type == 'Chondrocyte')
Secratory <- meta_data %>% filter(Cell_Type == 'Secratory')
Fibroblast <- meta_data %>% filter(Cell_Type == 'Fibroblast')
Immune <- meta_data %>% filter(Cell_Type == 'Immune')
Cilliated <- meta_data %>% filter(Cell_Type == 'Cilliated')
Muscle <- meta_data %>% filter(Cell_Type == 'Smooth muscle')



# plot the heavy to light for just one cell type to illustrate the point (fig2 A)

Basal_5dat <- Basal %>% filter(sample %in% c('1','2'))
Basal_10dat <- Basal %>% filter(sample %in% c('3','4'))
hist(log2(exp(p_all_alpha_abs['Q64524',Basal_5dat$ID]*5)-1),30,xlab = 'Theoretical H/L recycling adjusted')
abline(v=0,col='red')
abline(v=log2(.75/.25)+.15,col='red')


df_hist_plot_f <- data.frame(H_ov_L = log2(exp(p_all_alpha_abs['Q64524',Basal_5dat$ID]*5)-1),days = ' 5 day')
df_hist_plot_t <- data.frame(H_ov_L = log2(exp(p_all_alpha_abs['Q64524',Basal_10dat$ID]*10)-1),days = '10 day')
df_hist_plot <- rbind(df_hist_plot_f,df_hist_plot_t)

ggplot(df_hist_plot, aes(x = H_ov_L)) +
  geom_histogram() +
  facet_wrap(
    ~ days,
    ncol          = 1,           # one column → vertical stack
    strip.position = "right"     # facet labels on the right
  ) +
  xlim(-5, 4) +
  dot_plot +                     # your custom theme
  theme(
    strip.placement = "outside"  # put strips outside panels (optional)
  ) + xlab('Heavy/Light ratio, recycling adjusted')+ylab('# Single cells')



# Get growth rates

Compute_cell_growth_rate(Basal)
Compute_cell_growth_rate(Chondrocyte)
Compute_cell_growth_rate(Secratory)
Compute_cell_growth_rate(Fibroblast)
Compute_cell_growth_rate(Immune)
Compute_cell_growth_rate(Cilliated)
Compute_cell_growth_rate(Muscle)

cell_type = c('7Chondrocyte','6Smooth muscle','5Fibroblast','3Secratory','4Cilliated','2Basal','1Immune')
growth_rate <- c(0.001,         0.001,              0.01,        0.0179065, 0.01255309, 0.03,   0.06807215)

df_order_div <- data.frame(cell_type=cell_type,growth_rate=growth_rate)
ggplot(df_order_div, aes(x = cell_type,y = growth_rate)) + geom_point(size = 5)+
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))+xlab('')




### relation of regulatory steps to protein abundance

convert = Proc_fasta(paste0(path_dat,'Mouse.fasta'))
convert = convert %>% filter(split_prot %in% rownames(p_all_abs))
convert = convert %>% filter(split_gene %in% rownames(rna_seq@assays$RNA$counts))

convert_good <- convert %>% filter(split_prot %in% rownames(p_all_alpha_abs))



protein_abs_mat = matrix(data = NA,nrow = nrow(convert_good),ncol = 7)
deg_abs_mat = matrix(data = NA,nrow = nrow(convert_good),ncol = 7)
mRNA_mat = matrix(data = NA,nrow = nrow(convert_good),ncol = 7)

rely_deg <- c()
rely_abs <- c()
rely_mRNA <- c()
rely_trans <- c()

# rownames(deg_abs_mat) <- convert_good$split_prot
# colnames(deg_abs_mat) <- unique(mRNA_meta$cell_type)
# write.csv(log2(deg_abs_mat),'/Users/andrewleduc/Desktop/Projects/Miceotopes/SingleCell/Proc_files/deg_mat.csv')



for(i in 1:length(unique(mRNA_meta$cell_type))){
  
  # Get Cell type IDs
  um_plot_hold =  meta_data %>% filter(Cell_Type == unique(mRNA_meta$cell_type)[i])
  um_plot_hold <- um_plot_hold %>% filter(um_plot_hold$ID %in% colnames(p_all_abs))

  
  ## Account for missingness the two cell types that only have a few cells
  x = 1
  xx = 1
  if(unique(mRNA_meta$cell_type)[i] == 'Cilliated'){
    xx = 5
  }
  if(unique(mRNA_meta$cell_type)[i] =='Immune'){
    xx = 5
  }
  if(unique(mRNA_meta$cell_type)[i] =='Smooth muscle'){
    xx = 5
  }
  
  # Store protein deg (subsets)
  cells_hold <- intersect(colnames(p_all_alpha_abs),um_plot_hold$ID)
  cells_sub1 <- cells_hold[sample(1:length(cells_hold),length(cells_hold)/2)]
  cells_sub2 <- cells_hold[!cells_hold %in% cells_sub1]

  prot_i_sub1 = rowMeans(p_all_alpha_abs[convert_good$split_prot,cells_sub1],na.rm = T)
  prot_na <- prot_i_sub1
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_alpha_abs[convert_good$split_prot,cells_sub1])==F)  < 7] <- NA
  prot_i_sub1_deg = prot_i_sub1*prot_na

  prot_i_sub2 = rowMeans(p_all_alpha_abs[convert_good$split_prot,cells_sub2],na.rm = T)
  prot_na <- prot_i_sub2
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_alpha_abs[convert_good$split_prot,cells_sub2])==F)  < 7] <- NA
  prot_i_sub2_deg = prot_i_sub2*prot_na
  
  
  rely_deg <- c(rely_deg,cor(prot_i_sub1_deg,prot_i_sub2_deg,use = 'pairwise.complete.obs'))

  
  
  
  # Store protein deg (all cells)
  prot_i = rowMeans(p_all_alpha_abs[convert_good$split_prot,intersect(colnames(p_all_alpha_abs),um_plot_hold$ID)],na.rm = T)
  prot_na <- prot_i
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_alpha_abs[convert_good$split_prot,intersect(colnames(p_all_alpha_abs),um_plot_hold$ID)])==F)  < xx] <- NA
  prot_i = prot_i*prot_na
  deg_abs_mat[,i] = prot_i


  
  # Store protein Abundance (subsets)
  #cells_hold <- intersect(colnames(p_all_abs),rownames(um_plot_hold))
  #cells_sub1 <- cells_hold[sample(1:length(cells_hold),length(cells_hold)/2)]
  #cells_sub2 <- cells_hold[!cells_hold %in% cells_sub1]

  prot_i_sub1 = rowMeans(p_all_abs[convert_good$split_prot,cells_sub1],na.rm = T)
  # prot_na <- prot_i_sub1
  # prot_na[is.na(prot_na)==F]  <- 1
  # prot_na[rowSums(is.na(p_all_abs[convert_good$split_prot,cells_sub1])==F)  < x] <- NA
  # prot_i_sub1_abs = prot_i_sub1*prot_na

  prot_i_sub2 = rowMeans(p_all_abs[convert_good$split_prot,cells_sub2],na.rm = T)
  # prot_na <- prot_i_sub2
  # prot_na[is.na(prot_na)==F]  <- 1
  # prot_na[rowSums(is.na(p_all_abs[convert_good$split_prot,cells_sub2])==F)  < x] <- NA
  # prot_i_sub2_abs = prot_i_sub2*prot_na

  rely_abs <- c(rely_abs,cor(log2(prot_i_sub1_abs),log2(prot_i_sub2_abs),use = 'pairwise.complete.obs'))
  
  
  # Store protein abundance (all cells)
  prot_i = rowMeans(p_all_abs[convert_good$split_prot,intersect(colnames(p_all_abs),um_plot_hold$ID)],na.rm = T)
  prot_na <- prot_i
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_abs[convert_good$split_prot,intersect(colnames(p_all_abs),um_plot_hold$ID)])==F) < xx] <- NA
  prot_i = prot_i*prot_na
  protein_abs_mat[,i] = prot_i

  
  # Store mRNA abundance (all cells)
  mrna_meta_hold <- mRNA_meta %>% filter(cell_type == unique(mRNA_meta$cell_type)[i])
  
  rna_mat <- as.matrix(rna_seq@assays$RNA@counts[,rownames(mrna_meta_hold)])
  mRNA_mat[,i] = (.1+ rowMeans(rna_mat[convert_good$split_gene,]))


  
  # Store mRNA abundance (subsets)
  mRNA_sub1 <- (.1+ rowMeans(rna_mat[convert_good$split_gene,sample(1:ncol(rna_mat),50)]))
  mRNA_sub2 <- (.1+ rowMeans(rna_mat[convert_good$split_gene,sample(1:ncol(rna_mat),50)]))

  rely_mRNA <- c(rely_mRNA,cor(mRNA_sub1,mRNA_sub2))
  
  
  # Store translation (subsets)
  
  prot_i_sub1_trans <- prot_i_sub1_abs * prot_i_sub1_deg / mRNA_sub1
  prot_i_sub2_trans <- prot_i_sub2_abs * prot_i_sub2_deg / mRNA_sub2
  
  rely_trans <- c(rely_trans,cor(log2(prot_i_sub1_trans),log2(prot_i_sub2_trans),use = 'pairwise.complete.obs'))
  
  
  # 
  # df_plot <- data.frame(
  #   prot1 = (.1+ rowMeans(rna_mat[convert_good$split_gene,sample(1:ncol(rna_mat),100)])),
  #   prot2 = (.1+ rowMeans(rna_mat[convert_good$split_gene,sample(1:ncol(rna_mat),100)]))
  # )
  # 
  # # 2. Draw the plot
  # mRNA_plot <- ggplot(df_plot, aes(x = prot1, y = prot2)) +
  #   geom_point(size = 2, alpha = 0.7) +                 # points
  #   geom_smooth(method = "lm", se = FALSE,              # optional trend line
  #               colour = "red", linetype = "dashed") +
  #   scale_x_log10() +                                   # log10 x‑axis
  #   scale_y_log10() +                                   # log10 y‑axis
  #   labs(
  #     x = "Technical replicate 1",
  #     y = "Technical replicate 2",
  #     title = "# mRNA"
  #   ) +
  #   dot_plot
  # 
  # 
  # 
  # 
  # df_plot <- data.frame(
  #   prot1 = prot_i_sub1_deg,
  #   prot2 = prot_i_sub2_deg
  # )
  # 
  # # 2. Draw the plot
  # deg_plot <- ggplot(df_plot, aes(x = prot1, y = prot2)) +
  #   geom_point(size = 2, alpha = 0.7) +                 # points
  #   geom_smooth(method = "lm", se = FALSE,              # optional trend line
  #               colour = "red", linetype = "dashed") +
  #   scale_x_log10() +                                   # log10 x‑axis
  #   scale_y_log10() +                                   # log10 y‑axis
  #   labs(
  #     x = "Technical replicate 1",
  #     y = "Technical replicate 2",
  #     title = "Clearance rate, Kr (1/days)"
  #   ) +
  #   dot_plot
  # 
  # 
  # df_plot <- data.frame(
  #   prot1 = prot_i_sub1_abs,
  #   prot2 = prot_i_sub2_abs
  # )
  # 
  # # 2. Draw the plot
  # prot_plot <- ggplot(df_plot, aes(x = prot1, y = prot2)) +
  #   geom_point(size = 2, alpha = 0.7) +                 # points
  #   geom_smooth(method = "lm", se = FALSE,              # optional trend line
  #               colour = "red", linetype = "dashed") +
  #   scale_x_log10() +                                   # log10 x‑axis
  #   scale_y_log10() +                                   # log10 y‑axis
  #   labs(
  #     x = "Technical replicate 1",
  #     y = "Technical replicate 2",
  #     title = "Protein abundance, copies"
  #   ) +
  #   dot_plot
  # 
  # 
  # df_plot <- data.frame(
  #   prot1 = prot_i_sub1_trans,
  #   prot2 = prot_i_sub2_trans
  # )
  # 
  # # 2. Draw the plot
  # trans_plot <- ggplot(df_plot, aes(x = prot1, y = prot2)) +
  #   geom_point(size = 2, alpha = 0.7) +                 # points
  #   geom_smooth(method = "lm", se = FALSE,              # optional trend line
  #               colour = "red", linetype = "dashed") +
  #   scale_x_log10() +                                   # log10 x‑axis
  #   scale_y_log10() +                                   # log10 y‑axis
  #   labs(
  #     x = "Technical replicate 1",
  #     y = "Technical replicate 2",
  #     title = "Translation rate, Kt, copies protein * Day/copies mRNA"
  #   ) +
  #   dot_plot

  
}


## Plot Reliabilities
df_rely <- data.frame(clear = rely_deg,
                      abs  = rely_abs,
                      mRNA = rely_mRNA,
                      trans = rely_trans,
                      cell_type = unique(mRNA_meta$cell_type))

df_rely <- melt(df_rely,id.vars = 'cell_type')

ggplot(df_rely,aes(x = cell_type,y = value,color = variable))+ 
  ggbeeswarm::geom_beeswarm(size = 2) + theme_classic(base_size = 18) +
  ylim(c(0,1))+ ylab('Reliability') + theme(axis.text.x = element_text(angle=45,hjust=1))+
  xlab('') + scale_color_manual(values = c('black','grey50','orange','blue'))


trans_mat = protein_abs_mat*deg_abs_mat/mRNA_mat

rownames(trans_mat)<- convert_good$split_gene
colnames(trans_mat) <- unique(mRNA_meta$cell_type)

rownames(protein_abs_mat)<- convert_good$split_gene
colnames(protein_abs_mat) <- unique(mRNA_meta$cell_type)

rownames(mRNA_mat)<- convert_good$split_gene
colnames(mRNA_mat) <- unique(mRNA_meta$cell_type)

rownames(deg_abs_mat)<- convert_good$split_gene
colnames(deg_abs_mat) <- unique(mRNA_meta$cell_type)



#write.csv(trans_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/translation.csv'))
#write.csv(protein_abs_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/abundance.csv'))
#write.csv(mRNA_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/mRNA.csv'))
#write.csv(deg_abs_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/clearance.csv'))


### Make all the compare plots for Figure 2b

# Sec vs basal
rownames(deg_abs_mat)<- convert_good$split_gene
colnames(deg_abs_mat) <- unique(mRNA_meta$cell_type)
df_plot <- tibble(
  Basal   = deg_abs_mat[,1],
  Sec = deg_abs_mat[,3]
)

ggplot(df_plot, aes(x = Sec, y = Basal)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    x = "Sec",
    y = "Bas",
  ) +
  scale_x_log10(limits = c(.01,1))+
  scale_y_log10(limits = c(.01,1))+
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  ) + geom_abline(slope = 1,intercept = 0) 


TLS(log2(df_plot$Basal),log2(df_plot$Sec))[[1]]


# Chondrocyte vs basal
df_plot <- tibble(
  basal_log2   = rowMeans(p_all_alpha_abs[, Basal$ID], na.rm = TRUE),
  chondro_log2 = rowMeans(p_all_alpha_abs[, Chondrocyte$ID], na.rm = TRUE)
)

TLS(log2(df_plot$basal_log2),log2(df_plot$chondro_log2))[[1]]

cor(log2(df_plot$basal_log2),log2(df_plot$chondro_log2),use = 'pairwise.complete.obs')

ggplot(df_plot, aes(y = basal_log2, x = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  
  # reference line y = x
  geom_abline(intercept = 0, slope = 1,
              colour = "black", linetype = "solid") +
  
  # red dashed trend line (use geom_abline for fixed slope or stat_smooth for LM)
  geom_abline(intercept = 0.0, slope =  0.65,
              colour = "red", linetype = "dashed", size = 0.9) +
  
  scale_x_log10(limits = c(0.01, 1)) +   # x‑axis from 10⁻² to 1
  scale_y_log10(limits = c(0.01, 1))+
  labs(
    x = "Basal",
    y = "Chondrocyte",
    title = "Clearance rates"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )



# Chondrocyte vs basal (degredation not clearance, we are subtracting out the growth rate, Kg)
df_plot <- tibble(
  basal_log2   = rowMeans(p_all_alpha_comp[, Basal$ID]-.07, na.rm = TRUE),
  chondro_log2 = rowMeans(p_all_alpha_comp[, Chondrocyte$ID], na.rm = TRUE)
)

TLS(log2(df_plot$basal_log2),log2(df_plot$chondro_log2))[[1]]

ggplot(df_plot, aes(y = basal_log2, x = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  
  # reference line y = x
  geom_abline(intercept = 0, slope = 1,
              colour = "black", linetype = "solid") +
  
  # red dashed trend line (use geom_abline for fixed slope or stat_smooth for LM)
  geom_abline(intercept = 0.0, slope =  1,
              colour = "red", linetype = "dashed", size = 0.9) +
  
  scale_x_log10(limits = c(0.01, 1)) +   # x‑axis from 10⁻² to 1
  scale_y_log10(limits = c(0.01, 1))+
  labs(
    x = "Basal",
    y = "Chondrocyte",
    title = "Clearance rates"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )










mRNA_mat[log2(mRNA_mat) < -2] <- NA

colnames(deg_abs_mat) <- unique(mRNA_meta$cell_type)

hm_plot <- cor(log2(deg_abs_mat),use = 'pairwise.complete.obs')
diag(hm_plot) <- rely_deg
Heatmap(
  hm_plot,
  name = "corr",
  cluster_rows = T,
  cluster_columns = T,
  col = colorRamp2(c( 0, 1), c( "white", "red")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # i, j = row/column indices in cor_mat
    # x, y = center coordinates of the cell
    grid.text(
      round(hm_plot[i, j], 2),  # round to 2 decimal places
      x = x, y = y,
      gp = gpar(fontsize = 10)
    )
  }
)



rely_abs <- rely_abs * 0.7 # From literature estimate of trysin vs LysC similarity
rely_trans <- rely_trans * 0.7

rely_mRNA <- rely_mRNA * 0.7 # From literature estimate of 10x vs smart-seq bias
rely_deg <- rely_deg

deg_var <- c()
trans_var <- c()
mrna_var <- c()

for(i in 1:7){
  #deg_var <- c(deg_var,cor(log2(deg_abs_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson'))
  mrna_var <- c(mrna_var,cor(log2(mRNA_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson')^2/(rely_abs[i]*rely_mRNA[i]))
  deg_var <- c(deg_var,cor(log2(deg_abs_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson')^2/(rely_abs[i]*rely_deg[i]))
  trans_var <- c(trans_var,cor(log2(trans_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson')^2/(rely_abs[i]*rely_trans[i]))
  
}

(trans_var+mrna_var+deg_var)
trans_var_dif = 1-mrna_var-deg_var

plot(trans_var_dif,trans_var)

cell_type = c('Chondrocyte','Smooth muscle','Fibroblast','Secratory','Cilliated','Basal','Immune')

df_growth <- data.frame(ct = cell_type,kg = growth_rate)
rownames(df_growth) <- df_growth$ct
df_growth <- df_growth[unique(mRNA_meta$cell_type),]


plot(df_growth$kg,deg_var*100,ylab = 'Raw data % variance', xlab = 'Kg (growth rate)',main='pearson = -0.78')
cor(df_growth$kg,deg_var,method = 'pearson')

plot(df_growth$kg,trans_var_dif)
cor(df_growth$kg,trans_var_dif,method = 'pearson')

plot(df_growth$kg,mrna_var) 
cor(df_growth$kg,mrna_var,method = 'pearson')



df_growth$zdeg <- deg_var
df_growth$atrans <- trans_var_dif
df_growth$trans_reg <- trans_var
df_growth$mrna <- mrna_var

df_growth$ct <- factor(df_growth$ct,   levels = names(palette))

ggplot(df_growth,aes(x = kg,y = zdeg,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(axis.line      = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+scale_color_manual(values = palette)+
  ylab('Fraction protein abundance explained') + xlab('Growth rate') + ggtitle('Degradation')


ggplot(df_growth,aes(x = kg,y = atrans,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(axis.line      = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+scale_color_manual(values = palette)+
  ylab('Fraction protein abundance explained') + xlab('Growth rate')  + ggtitle('Translation')+
  scale_y_continuous() 


ggplot(df_growth,aes(x = kg,y = mrna,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(axis.line      = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+scale_color_manual(values = palette)+
  ylab('Fraction protein abundance explained') + xlab('Growth rate') + ggtitle('mRNA')

ggplot(df_growth,aes(x = trans_reg,y = atrans,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(axis.line      = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+scale_color_manual(values = palette)+
  ylab('Remaining variance estimate') + xlab('Regression estimate') + ggtitle('Translation')+
  scale_y_continuous(
   # 0.0, 0.1, 0.2, …
    # or, if you prefer the older helper:
    # labels = number_format(accuracy = 0.1)
  ) 

ggplot(df_growth,aes(x = kg,y = trans_reg,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(axis.line      = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+scale_color_manual(values = palette)+
  ylab('Regression estimate') + xlab('') + ggtitle('Translation')+
  scale_y_continuous(
     # 0.0, 0.1, 0.2, …
    # or, if you prefer the older helper:
    # labels = number_format(accuracy = 0.1)
  ) 


cor(df_growth$kg,df_growth$trans_reg)




################
# mRNA scatter plot for supplemental fig
################

df_plot <- tibble(
  basal_log2   = mRNA_mat[,3],
  chondro_log2 = (protein_abs_mat[,3])
)

mrna_var[7]
ggplot(df_plot, aes(x = basal_log2, y = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_x_log10() +   # x‑axis from 10⁻² to 1
  scale_y_log10()+ #limits = c(0.01, 1)
  labs(
    x = "mRNA copies",
    y = "Abundance",
    title = "Secratory"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )




################
# percent variance of clearnce fig 2 c
################


p_all_alpha_comp <- p_all_alpha_abs[rowSums(is.na(p_all_alpha_abs)==F) > 1000,]


test <- p_all_alpha_comp[,Cilliated$ID]
test[rowSums(is.na(test)==F) < 20] <- NA

dd <- c(sd(log2(rowMeans(p_all_alpha_comp[,Chondrocyte$ID],na.rm=T)),na.rm=T),
        sd(log2(rowMeans(p_all_alpha_comp[,Muscle$ID],na.rm=T)),na.rm=T),
        sd(log2(rowMeans(p_all_alpha_comp[,Fibroblast$ID],na.rm=T)),na.rm=T),
        sd(log2(rowMeans(p_all_alpha_comp[,Secratory$ID],na.rm=T)),na.rm=T),
        sd(log2(rowMeans(test,na.rm=T)),na.rm=T),
        sd(log2(rowMeans(p_all_alpha_comp[,Basal$ID],na.rm=T)),na.rm=T),
        sd(log2(rowMeans(p_all_alpha_comp[,Immune$ID],na.rm=T)),na.rm=T))

ff <- c(0,0,.005,.01,.01,.015,.02)

cell_type = c('Chondrocyte','Muscle','Fibroblast','Secratory','Cilliated','Basal','Immune')

df_scat <- data.frame(sds = dd, Kg = ff, cell_type = cell_type)

ggplot(df_scat,aes(x = Kg, y = sds,color = cell_type)) + geom_point(size = 5) +xlab('Growth rate, Kg') +
  ylab('Cor(Clearance,Abundance)')+
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )





