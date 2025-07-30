# Absolute level analysis, Figures 1f and Figure 2
library(ggplot2)
library(dplyr)
library(seqinr)
library(stringr)
library(reshape2)


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


#### Read in data

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


rna_seq <- readRDS('/Users/andrewleduc/Desktop/Senescense/seurat_integrated_filered_700_named.rds')


p_all_abs <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_absolute.csv'),row.names = 1)
p_all_abs <- as.matrix(p_all_abs)
p_all_alpha_abs <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/clearance_absolute.csv'),row.names = 1)
p_all_alpha_abs <- as.matrix(p_all_alpha_abs)

um_plot <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),row.names = 1)
mRNA_meta <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_mRNA_annotations.csv'),row.names = 1)

convert = Proc_fasta(paste0(path_dat,'Mouse.fasta'))

#### Read in data

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



mRNA_meta <- mRNA_meta %>% filter(cell_type != 'Shwann')
for(i in 1:length(unique(mRNA_meta$cell_type))){
  
  # Get Cell type IDs
  um_plot_hold = um_plot %>% filter(Cell_Type == unique(mRNA_meta$cell_type)[i])
  um_plot_hold <- um_plot_hold %>% filter(um_plot_hold$ID %in% colnames(p_all_abs))

  
  x = 1
  xx = 15
  if(unique(mRNA_meta$cell_type)[i] == 'Smooth muscle'){
    #x = 100
    xx = 5
  }
  if(unique(mRNA_meta$cell_type)[i] == 'Cilliated'){
    print('here')
    xx = 5
  }
  if(unique(mRNA_meta$cell_type)[i] =='Immune'){
    print('here')
    #x = 
    xx = 5
  }
  
  # Store protein deg (subsets)
  cells_hold <- intersect(colnames(p_all_alpha_abs),um_plot_hold$ID)
  cells_sub1 <- cells_hold[sample(1:length(cells_hold),length(cells_hold)/2)]
  cells_sub2 <- cells_hold[!cells_hold %in% cells_sub1]

  prot_i_sub1 = rowMeans(p_all_alpha_abs[convert_good$split_prot,cells_sub1],na.rm = T)
  prot_na <- prot_i_sub1
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_alpha_abs[convert_good$split_prot,cells_sub1])==F)  < 5] <- NA
  prot_i_sub1_deg = prot_i_sub1*prot_na

  prot_i_sub2 = rowMeans(p_all_alpha_abs[convert_good$split_prot,cells_sub2],na.rm = T)
  prot_na <- prot_i_sub2
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_alpha_abs[convert_good$split_prot,cells_sub2])==F)  < 5] <- NA
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
  prot_na <- prot_i_sub1
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_abs[convert_good$split_prot,cells_sub1])==F)  < x] <- NA
  prot_i_sub1_abs = prot_i_sub1*prot_na

  prot_i_sub2 = rowMeans(p_all_abs[convert_good$split_prot,cells_sub2],na.rm = T)
  prot_na <- prot_i_sub2
  prot_na[is.na(prot_na)==F]  <- 1
  prot_na[rowSums(is.na(p_all_abs[convert_good$split_prot,cells_sub2])==F)  < x] <- NA
  prot_i_sub2_abs = prot_i_sub2*prot_na

  
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



write.csv(trans_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/translation.csv'))
write.csv(protein_abs_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/abundance.csv'))
write.csv(mRNA_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/mRNA.csv'))
write.csv(deg_abs_mat,paste0(path_dat,'06_Gene_X_CellType/Absolute_abundance/clearance.csv'))



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


colnames(deg_abs_mat)

plot(log2(deg_abs_mat[,2]),log2(protein_abs_mat[,2]))


mRNA_mat[log2(mRNA_mat) < -2] <- NA
synth_mat = protein_abs_mat*deg_abs_mat



plot(log2(deg_abs_mat[,1]/deg_abs_mat[,2]) ,log2(protein_abs_mat[,1]/protein_abs_mat[,2]))
cor(log2(deg_abs_mat[,1]/deg_abs_mat[,2]) ,log2(protein_abs_mat[,1]/protein_abs_mat[,2]),use = 'pairwise.complete.obs',method = 'pearson')


unique(mRNA_meta$cell_type)
plot(log2(deg_abs_mat[,6]) ,log2(protein_abs_mat[,6]))

plot(log2(deg_abs_mat[,1]/deg_abs_mat[,2]) ,log2(protein_abs_mat[,1]/protein_abs_mat[,2]))

cor(log2(deg_abs_mat[,6])  ,log2(protein_abs_mat[,6]),use = 'pairwise.complete.obs',method = 'pearson')

plot(log2(synth_mat[,5]) ,log2(protein_abs_mat[,5]))
cor(log2(synth_mat[,5])  ,log2(protein_abs_mat[,5]),use = 'pairwise.complete.obs',method = 'pearson')^2

plot( log2(mRNA_mat[,5]),log2(protein_abs_mat[,5]))
cor( log2(mRNA_mat[,5]) ,log2(protein_abs_mat[,5]),use = 'pairwise.complete.obs',method = 'pearson')

plot(log2(synth_mat[,5]) - log2(mRNA_mat[,5]),log2(mRNA_mat[,5]))
cor(log2(synth_mat[,5]) - log2(mRNA_mat[,5]) ,log2(protein_abs_mat[,5]),use = 'pairwise.complete.obs',method = 'pearson')

plot(log2(trans_mat[,5]),log2(protein_abs_mat[,5]))
cor(log2(trans_mat[,5]) ,log2(protein_abs_mat[,5]),use = 'pairwise.complete.obs',method = 'pearson')




colnames(deg_abs_mat) <- unique(mRNA_meta$cell_type)
Heatmap(
  cor(deg_abs_mat,use = 'pairwise.complete.obs'),
  name = "corr",
  cluster_rows = T,
  cluster_columns = T,
  col = colorRamp2(c( 0, 1), c( "white", "red")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # i, j = row/column indices in cor_mat
    # x, y = center coordinates of the cell
    grid.text(
      round(cor(deg_abs_mat,use = 'pairwise.complete.obs')[i, j], 2),  # round to 2 decimal places
      x = x, y = y,
      gp = gpar(fontsize = 10)
    )
  }
)



rely_abs <- rely_abs-.15
rely_mRNA <- rely_mRNA * 0.75
rely_deg

deg_var <- c()
trans_var <- c()
mrna_var <- c()
for(i in 1:7){

  mrna_var <- c(mrna_var,cor(log2(mRNA_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson')^2/(rely_abs[i]*rely_mRNA[i]))
  deg_var <- c(deg_var,cor(log2(deg_abs_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson')^2/(rely_abs[i]*rely_deg[i]))
  trans_var <- c(trans_var,cor(log2(trans_mat[,i]),log2(protein_abs_mat[,i]),use = 'pairwise.complete.obs',method = 'pearson')^2/(rely_abs[i]*rely_trans[i]))
  
}

(trans_var+mrna_var+deg_var)
trans_var = 1-mrna_var-deg_var

plot(df_growth$kg,deg_var*100,ylab = 'Raw data % variance', xlab = 'Kg (growth rate)',main='pearson = -0.78')
cor(df_growth$kg,deg_var,method = 'pearson')

plot(df_growth$kg,trans_var)
cor(df_growth$kg,trans_var,method = 'pearson')

plot(df_growth$kg,mrna_var) 
cor(df_growth$kg,mrna_var,method = 'pearson')




cell_type = c('Chondrocyte','Smooth muscle','Fibroblast','Secratory','Cilliated','Basal','Immune')

df_growth <- data.frame(ct = cell_type,kg = growth_rate)
rownames(df_growth) <- df_growth$ct
df_growth <- df_growth[unique(mRNA_meta$cell_type),]

df_growth$zdeg <- deg_var
df_growth$atrans <- trans_var
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
  scale_y_continuous(
    labels = label_number(accuracy = 0.1)   # 0.0, 0.1, 0.2, …
    # or, if you prefer the older helper:
    # labels = number_format(accuracy = 0.1)
  ) 


ggplot(df_growth,aes(x = kg,y = mrna,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(axis.line      = element_blank(), 
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+scale_color_manual(values = palette)+
  ylab('Fraction protein abundance explained') + xlab('Growth rate') + ggtitle('mRNA')+
  scale_y_continuous(
    labels = label_number(accuracy = 0.1)   # 0.0, 0.1, 0.2, …
    # or, if you prefer the older helper:
    # labels = number_format(accuracy = 0.1)
  ) 




################
# mRNA 
################
unique(mRNA_meta$cell_type)
df_plot <- tibble(
  basal_log2   = mRNA_mat[,7],
  chondro_log2 = (protein_abs_mat[,7])
)

mrna_var[7]
ggplot(df_plot, aes(x = basal_log2, y = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_x_log10() +   # x‑axis from 10⁻² to 1
  scale_y_log10()+ #limits = c(0.01, 1)
  labs(
    x = "mRNA copies",
    y = "Abundance",
    title = "Smooth Muscle"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )



df_plot <- tibble(
  basal_log2   = mRNA_mat[,3],
  chondro_log2 = (protein_abs_mat[,3])
)

mrna_var[1]

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


df_plot <- tibble(
  basal_log2   = log2(mRNA_mat[,4]/mRNA_mat[,7]) - median(log2(mRNA_mat[,4]/mRNA_mat[,7]),na.rm=T),
  chondro_log2 = log2(protein_abs_mat[,4]/protein_abs_mat[,7]) - median(log2(protein_abs_mat[,4]/protein_abs_mat[,7]),na.rm=T)
)


ggplot(df_plot, aes(x = basal_log2, y = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    x = "log2 mRNA fold change",
    y = "log2 Abundance fold change",
    title = "Smooth muscle / Fibroblast"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  ) + geom_abline(slope = 1,intercept = 0) +
  xlim(c(-5,5))+ylim(c(-5,5))


cor(df_plot$basal_log2,df_plot$chondro_log2,use = 'pairwise.complete.obs',method = 'spearman')




##########

df_plot <- tibble(
  basal_log2   = log2(deg_abs_mat[,1]/deg_abs_mat[,2]) - median(log2(deg_abs_mat[,1]/deg_abs_mat[,2]),na.rm=T),
  chondro_log2 = log2(protein_abs_mat[,1]/protein_abs_mat[,2]) - median(log2(protein_abs_mat[,1]/protein_abs_mat[,2]),na.rm=T)
)


ggplot(df_plot, aes(x = basal_log2, y = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  labs(
    x = "log2 Degradation fold change",
    y = "log2 Abundance fold change",
    title = "Smooth muscle / Fibroblast"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  ) + geom_abline(slope = -1,intercept = 0) +
  xlim(c(-5,5))+ylim(c(-5,5))



df_growth$sds <-  colSds(log2(deg_abs_mat),na.rm=T)

ggplot(df_growth,aes(x = (kg),y = sds,color = ct)) + geom_point(size = 5) + 
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+
  ylab('Clearance rate variability (std.)') + xlab('Growth rate')
cor(df_growth$kg,df_growth$sds)




library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(scales)
df_long <- df_growth %>% 
  pivot_longer(
    cols      = c(zdeg, atrans, mrna),
    names_to  = "variable",
    values_to = "value"
  ) %>% 
  
  ## ── order x‑axis by kg ────────────────────────────────────────────
  mutate(ct = fct_reorder(ct, kg, .desc = FALSE))

ggplot(df_long, aes(x = ct, y = value, fill = variable)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(zdeg = "#D27628", atrans = "#F6AA4C", mrna = "gray20")) +
  labs(x = "Cell type (ordered by growth rate)", y = "Fraction of variance explained") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8),
    axis.text.x  = element_text(angle = 45, hjust = 1)
  )

p_all['P18242',]
rownames(p_all)
## Compute cell division rates

colnames(meta_data)[1] <- 'ID'

Basal <- meta_data %>% filter(Cell_Type == 'Basal')
Basal_1 <- Basal[sample(1:nrow(Basal),nrow(Basal)/2),]
Basal_2 <- Basal %>% filter(!ID %in% Basal_1$ID)
Chondrocyte <- meta_data %>% filter(Cell_Type == 'Chondrocyte')
Chondrocyte_1 <- Chondrocyte[sample(1:nrow(Chondrocyte),nrow(Chondrocyte)/2),]
Chondrocyte_2 <- Chondrocyte %>% filter(!ID %in% Chondrocyte_1$ID)
Secratory <- meta_data %>% filter(Cell_Type == 'Secratory')
Secratory_1 <- Secratory[sample(1:nrow(Secratory),nrow(Secratory)/2),]
Secratory_2 <- Secratory %>% filter(!ID %in% Secratory_1$ID)
Fibroblast <- meta_data %>% filter(Cell_Type == 'Fibroblast')
Fibroblast_1 <- Fibroblast[sample(1:nrow(Fibroblast),nrow(Fibroblast)/2),]
Fibroblast_2 <- Fibroblast %>% filter(!ID %in% Fibroblast_1$ID)
Immune <- meta_data %>% filter(Cell_Type == 'Immune')
Immune_1 <- Immune[sample(1:nrow(Immune),nrow(Immune)/2),]
Immune_2 <- Immune %>% filter(!ID %in% Immune_1$ID)
Cilliated <- meta_data %>% filter(Cell_Type == 'Cilliated')

Muscle <- meta_data %>% filter(Cell_Type == 'Smooth muscle')



prot_1 <- r1_5day_male@miceotopes@Alpha_prot
colnames(prot_1) <- paste0(colnames(prot_1),'_prep1')

prot_2 <- r2_5day_female@miceotopes@Alpha_prot
colnames(prot_2) <- paste0(colnames(prot_2),'_prep2')

prot_3 <- r3_10day_male@miceotopes@Alpha_prot
colnames(prot_3) <- paste0(colnames(prot_3),'_prep3')

prot_4 <- r4_10day_female@miceotopes@Alpha_prot
colnames(prot_4) <- paste0(colnames(prot_4),'_prep4')


p1_alpha <- as.data.frame(prot_1)
p1_alpha$prot <- rownames(p1_alpha)
p2_alpha <- as.data.frame(prot_2)
p2_alpha$prot <- rownames(p2_alpha)
p3_alpha <- as.data.frame(prot_3)
p3_alpha$prot <- rownames(p3_alpha)
p4_alpha <- as.data.frame(prot_4)
p4_alpha$prot <- rownames(p4_alpha)
p_1and2 <- p1_alpha %>% merge(p2_alpha, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_alpha, by = 'prot',all = TRUE)
p_all_alpha <- p_123 %>% merge(p4_alpha, by = 'prot',all = TRUE)
rownames(p_all_alpha) <- p_all_alpha$prot
p_all_alpha$prot <- NULL
p_all_alpha <- as.matrix(p_all_alpha)
p_all_alpha_abs <- p_all_alpha[,colnames(protein_Data)]


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

cell_type = c('7Chondrocyte','6Smooth muscle','5Fibroblast','3Secratory','4Cilliated','2Basal','1Immune')
growth_rate <- c(0.001,         0.001,              0.01,        0.0179065, 0.01255309, 0.03,   0.06807215)

df_order_div <- data.frame(cell_type=cell_type,growth_rate=growth_rate)
ggplot(df_order_div, aes(x = cell_type,y = growth_rate)) + geom_point(size = 5)+
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))+xlab('')


1032-305/2
994-409/2

(log(1032) - log(879.5))/5
(log(994) - log(789.5-50))/10

hist(log2(exp(p_all_alpha_abs['Q64524',Basal_10dat$ID]*10)-1),30,xlab = 'Theoretical H/L recycling adjusted')
abline(v=0,col='red')
abline(v=log2(.75/.25)+.15,col='red')







Fibroblast_5dat <- Secratory %>% filter(sample %in% c('1','2'))
Fibroblast_10dat <- Secratory %>% filter(sample %in% c('3','4'))
hist(log2(exp(p_all_alpha_abs['Q64524',Fibroblast_5dat$ID]*5)-1),30,xlab = 'Theoretical H/L recycling adjusted')
abline(v=0,col='red')

hist(log2(exp(p_all_alpha_abs['Q64524',Fibroblast_10dat$ID]*10)-1),30,xlab = 'Theoretical H/L recycling adjusted')
abline(v=0,col='red')
abline(v=log2(.75/.25)+.15,col='red')


N_0 = sum(is.na(exp(p_all_alpha_abs['Q64524',Fibroblast_5dat$ID]*5)-1) ==F)-
  sum(log2(exp(p_all_alpha_abs['Q64524',Fibroblast_5dat$ID]*5)-1) > -1,na.rm=T)/2

N_t = sum(is.na(exp(p_all_alpha_abs['Q64524',Fibroblast_5dat$ID]*5)-1) ==F)


N_0 = sum(is.na(exp(p_all_alpha_abs['Q64524',Fibroblast_10dat$ID]*10)-1) ==F)-
  sum(log2(exp(p_all_alpha_abs['Q64524',Fibroblast_10dat$ID]*10)-1) > -1,na.rm=T)/2

N_t = sum(is.na(exp(p_all_alpha_abs['Q64524',Fibroblast_10dat$ID]*10)-1) ==F)

(log2(N_t) - log2(N_0))/5

# Full clearance distribution
dist_plot <- melt(p_all_alpha)
dist_plot$value <- log(2)/dist_plot$value
dist_plot$value[dist_plot$value > 60] <- NA
dist_plot <- dist_plot %>% left_join(meta_data,by = c('Var2'='ID'))
dist_plot <- dist_plot %>% filter(Cell_Type != 'Shwann')
ggplot(dist_plot,aes(x = Cell_Type,y = value)) + geom_violin()+ geom_boxplot(width = .2,outlier.shape = NA)+
  scale_y_log10() + ylab('Clearance half life (Days)')+ theme_minimal(base_size = 15)+
  xlab('') + coord_cartesian(ylim = c(.5,30)) +theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))



p_all_alpha_comp <- p_all_alpha_abs[rowSums(is.na(p_all_alpha_abs)==F) > 1000,]
plot(log2(rowMeans(p_all_alpha_comp[,Basal$ID],na.rm=T)),log2(rowMeans(p_all_alpha_comp[,Chondrocyte$ID],na.rm=T)),
     xlim = c(-5,1), ylim = c(-5,1))

abline(a = 0,b=1)
abline(a = 0,b=1/.65)

df_plot <- data.frame(
  Chond   = rowMeans(p_all_alpha_comp[, Chondrocyte$ID],      na.rm = TRUE),
  Muscle = rowMeans(p_all_alpha_comp[, Muscle$ID], na.rm = TRUE),
  Fib = rowMeans(p_all_alpha_comp[,Fibroblast$ID],na.rm=T),
  Secratory = rowMeans(p_all_alpha_comp[,Secratory$ID],na.rm=T),
  Cil = rowMeans(test,na.rm=T),
  Bas = rowMeans(p_all_alpha_comp[,Basal$ID],na.rm=T),
  Immune = rowMeans(p_all_alpha_comp[,Immune$ID],na.rm=T)
)

df_plot <- data.frame(
  Chond   = rowMeans(p_all_alpha_comp[, Chondrocyte$ID],      na.rm = TRUE),
  Secratory = rowMeans(p_all_alpha_comp[,Secratory$ID],na.rm=T),
  Bas = rowMeans(p_all_alpha_comp[,Basal$ID],na.rm=T)
)

df_long <- df_plot %>% 
  # (optional) keep original row names as an identifier
  tibble::rownames_to_column("row_id") %>% 
  
  # ── 1.  melt wide → long ───────────────────────────────────────────
  pivot_longer(
    cols       = -row_id,          # every column except the row id
    names_to   = "cell_type",
    values_to  = "value"
  ) %>% 
  
  # (optional) drop missing measurements
  filter(!is.na(value)) %>%        
  
  # ── 2–4.  rank values within each cell type ───────────────────────
  group_by(cell_type) %>% 
  arrange(desc(value), .by_group = TRUE) %>%   # use `value` for ascending
  mutate(order = row_number()) %>%            # 1N in the chosen order
  ungroup()


ggplot(df_long,aes(x = -order,y = value,color = cell_type)) + geom_point() + 
  scale_y_log10(limits = c(0.05, .8))+
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+
  ylab('Clearance rates') + xlab('Rank')



df_plot <- tibble(
  basal_log2   = rowMeans(p_all_alpha_comp[, Basal$ID]-.05, na.rm = TRUE),
  chondro_log2 = rowMeans(p_all_alpha_comp[, Chondrocyte$ID], na.rm = TRUE)
)

TLS(log2(df_plot$basal_log2),log2(df_plot$chondro_log2))[[1]]

ggplot(df_plot, aes(x = basal_log2, y = chondro_log2)) +
  geom_point(alpha = 0.8, size = 2) +
  
  # reference line y = x
  geom_abline(intercept = 0, slope = 1,
              colour = "black", linetype = "solid") +
  
  # red dashed trend line (use geom_abline for fixed slope or stat_smooth for LM)
  geom_abline(intercept = 0.0, slope = 1 / 0.65,
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

df_scat <- data.frame(sds = cors, Kg = ff, cell_type = cell_type)

ggplot(df_scat,aes(x = Kg, y = cors,color = cell_type)) + geom_point(size = 5) +xlab('Growth rate, Kg') +
  ylab('Cor(Clearance,Abundance)')+
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )


prot_1 <- r1_5day_male@matricies@protein_abs
colnames(prot_1) <- paste0(colnames(prot_1),'_prep1')

prot_2 <- r2_5day_female@matricies@protein_abs
colnames(prot_2) <- paste0(colnames(prot_2),'_prep2')

prot_3 <- r3_10day_male@matricies@protein_abs
colnames(prot_3) <- paste0(colnames(prot_3),'_prep3')

prot_4 <- r4_10day_female@matricies@protein_abs
colnames(prot_4) <- paste0(colnames(prot_4),'_prep4')


p1_alpha <- as.data.frame(prot_1)
p1_alpha$prot <- rownames(p1_alpha)
p2_alpha <- as.data.frame(prot_2)
p2_alpha$prot <- rownames(p2_alpha)
p3_alpha <- as.data.frame(prot_3)
p3_alpha$prot <- rownames(p3_alpha)
p4_alpha <- as.data.frame(prot_4)
p4_alpha$prot <- rownames(p4_alpha)
p_1and2 <- p1_alpha %>% merge(p2_alpha, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_alpha, by = 'prot',all = TRUE)
p_all_prot <- p_123 %>% merge(p4_alpha, by = 'prot',all = TRUE)
rownames(p_all_prot) <- p_all_prot$prot
p_all_prot$prot <- NULL
p_all_prot <- as.matrix(p_all_prot)
p_all_prot <- p_all_prot[,colnames(protein_Data)]


p_all_alpha_comp <- p_all_alpha_abs[rowSums(is.na(p_all_alpha_abs)==F) > 200,]

plot(rowMeans(log2(p_all_alpha_comp[, Chondrocyte$ID]),na.rm = TRUE),
     rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Chondrocyte$ID]),na.rm = TRUE))

cor(rowMeans(log2(p_all_alpha_comp[, Chondrocyte$ID]),na.rm = TRUE),
    rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Chondrocyte$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
    method = 'spearman')


plot(rowMeans(log2(p_all_alpha_comp[, Basal$ID]),na.rm = TRUE),
     rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Basal$ID]),na.rm = TRUE))

cor(rowMeans(log2(p_all_alpha_comp[, Basal$ID]),na.rm = TRUE),
    rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Basal$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
    method = 'spearman')

df_plot <- tibble(
  deg   = rowMeans(p_all_alpha_comp[, Chondrocyte$ID]-.05, na.rm = TRUE),
  Abs = rowMeans(p_all_prot[rownames(p_all_alpha_comp), Chondrocyte$ID], na.rm = TRUE)
)

TLS(log2(df_plot$basal_log2),log2(df_plot$chondro_log2))[[1]]

ggplot(df_plot, aes(x = deg, y = Abs)) +
  geom_point(alpha = 0.8, size = 2) +
  
  # reference line y = x
  geom_abline(intercept = 0, slope = 1,
              colour = "black", linetype = "solid") +
  
  # red dashed trend line (use geom_abline for fixed slope or stat_smooth for LM)
  #geom_abline(intercept = 0.0, slope = 1 / 0.65,
  #           colour = "red", linetype = "dashed", size = 0.9) +
  
  scale_x_log10() +   # x‑axis from 10⁻² to 1
  scale_y_log10()+
  labs(
    x = "Clearance rate",
    y = "Abundance",
    title = "Basal"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )









plot(ff,cors)
plot(dd,cors)

cor(ff,cors)


p_all_alpha_comp <- p_all_alpha_abs[rowSums(is.na(p_all_alpha_abs)==F) > 1000,]
nrow(p_all_alpha_comp)

cors <- c(cor(rowMeans(log2(p_all_alpha_comp[, Chondrocyte$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Chondrocyte$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'),
          cor(rowMeans(log2(p_all_alpha_comp[, Muscle$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Muscle$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'),
          cor(rowMeans(log2(p_all_alpha_comp[, Fibroblast$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Fibroblast$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'),
          cor(rowMeans(log2(p_all_alpha_comp[, Secratory$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Secratory$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'),
          cor(rowMeans(log2(p_all_alpha_comp[, Cilliated$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Cilliated$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'),
          
          cor(rowMeans(log2(p_all_alpha_comp[, Basal$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Basal$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'),
          
          cor(rowMeans(log2(p_all_alpha_comp[, Immune$ID]),na.rm = TRUE),
              rowMeans(log2(p_all_prot[rownames(p_all_alpha_comp), Immune$ID]),na.rm = TRUE),use = 'pairwise.complete.obs',
              method = 'pearson'))


Muscle = rowMeans(p_all_alpha_comp[, Muscle$ID], na.rm = TRUE)
Fib = rowMeans(p_all_alpha_comp[,Fibroblast$ID],na.rm=T)
Secratory = rowMeans(p_all_alpha_comp[,Secratory$ID],na.rm=T)
Cil = rowMeans(test,na.rm=T)
Bas = rowMeans(p_all_alpha_comp[,Basal$ID],na.rm=T)
Immune = rowMeans(p_all_alpha_comp[,Immune$ID],na.rm=T)



p_all_alpha_comp <- p_all_alpha_abs[rowSums(is.na(p_all_alpha_abs)==F) > 400,]
cdif <- rowMeans(p_all_alpha_comp[, Chondrocyte$ID],na.rm = TRUE)/rowMeans(p_all_alpha_comp[, Basal$ID],na.rm = TRUE)
ddif <- rowMeans(p_all_alpha_comp[, Chondrocyte$ID],na.rm = TRUE)/rowMeans(p_all_alpha_comp[, Basal$ID]-.05,na.rm = TRUE)
pdif <- rowMeans(p_all_prot[rownames(p_all_alpha_comp), Chondrocyte$ID],na.rm = TRUE)/rowMeans(p_all_prot[rownames(p_all_alpha_comp), Basal$ID],na.rm = TRUE)

ddif[ddif >2] <- NA
sd(log2(cdif),na.rm=T)
sd(log2(ddif),na.rm=T)

plot(log2(rowMeans(p_all_alpha_comp[, Basal$ID],na.rm = TRUE)),log2(rowMeans(p_all_alpha_comp[, Chondrocyte$ID],na.rm = TRUE)))
abline(a = 0,b=1)
plot(log2(rowMeans(p_all_alpha_comp[, Basal$ID]-.05,na.rm = TRUE)),log2(rowMeans(p_all_alpha_comp[, Chondrocyte$ID],na.rm = TRUE)))
abline(a = 0,b=1)


two <- ((log2(rowMeans(p_all_alpha_comp[, Chondrocyte$ID],na.rm = TRUE))-
           log2(rowMeans(p_all_alpha_comp[, Basal$ID],na.rm = TRUE))))


df_plot <- tibble(
  deg   = two,
  Abs = pdif
)
df_plot$Abs <- log2(df_plot$Abs/ median(df_plot$Abs,na.rm=T))
df_plot$deg <- df_plot$deg - median(df_plot$deg,na.rm=T)


TLS(log2(df_plot$deg),log2(df_plot$Abs))[[1]]
cor(log2(df_plot$deg),log2(df_plot$Abs),use = 'pairwise.complete.obs')^2

mean(abs(log2(df_plot$deg)- log2(df_plot$Abs)),na.rm=T)

ggplot(df_plot, aes(x = deg, y = Abs)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_abline(intercept = 0,slope = -1)+
  labs(
    x = "Clearance (Basal / Chondrocyte)",
    y = "Abundance (Basal / Chondrocyte)",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 0.8)  # add border
  )+
  annotate(
    "text",
    x      = 4,               # choose a position within the x‑limits
    y      = 4,            # choose a y‑value > 0 (log scale)
    label  = expression(R^2 == 0.25),
    hjust  = 1,               # right‑align
    vjust  = 0,               # top‑align
    size   = 5
  )


