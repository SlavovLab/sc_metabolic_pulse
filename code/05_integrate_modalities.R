library(dplyr)
library(ggplot2)
library(stringr)
library(seqinr)
library(QuantQC)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)

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


#### Paths

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'



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


rna_seq <- readRDS('/Users/andrewleduc/Desktop/Senescense/seurat_integrated_filered_700_named.rds')


p_all <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv'),row.names = 1)
p_all <- as.matrix(p_all)
p_all_alpha <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/clearance_relative.csv'),row.names = 1)
p_all_alpha <- as.matrix(p_all_alpha)

um_plot <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),row.names = 1)
mRNA_meta <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_mRNA_annotations.csv'),row.names = 1)

prot_imp_meta <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/Missingness_results.csv'),row.names = 1)

convert = Proc_fasta(paste0(path_dat,'Mouse.fasta'))


#################################################################################
# Compare mRNA vs prot vs deg across Cell types
#################################################################################
convert = convert %>% filter(split_prot %in% rownames(p_all))
convert = convert %>% filter(split_gene %in% rownames(rna_seq@assays$RNA$counts))

convert_good <- convert %>% filter(split_prot %in% rownames(p_all_alpha))


protein_mat = matrix(data = NA,nrow = nrow(convert_good),ncol = 7)
deg_mat = matrix(data = NA,nrow = nrow(convert_good),ncol = 7)
mRNA_mat = matrix(data = NA,nrow = nrow(rna_seq@assays$RNA@counts),ncol = 7)



p_all_min <- p_all

p_all_min <- Normalize_reference_vector_log(p_all_min)

data_miss_post_collapse3 <- prot_imp_meta %>% filter(Protein %in% convert_good$split_prot)

mRNA_counts <- as.matrix(rna_seq@assays$RNA@counts)


mRNA_meta <- mRNA_meta %>% filter(cell_type != 'Shwann')
for(i in 1:length(unique(mRNA_meta$cell_type))){

    um_plot_hold = um_plot %>% filter(Cell_Type == unique(mRNA_meta$cell_type)[i])
    mRNA_meta_hold = mRNA_meta %>% filter(cell_type == unique(mRNA_meta$cell_type)[i])
    
    mRNA_mat[,i] = (.1+ rowMeans(mRNA_counts[,rownames(mRNA_meta_hold)]))
    
    
    data_miss_post_collapse_bas <- data_miss_post_collapse3 %>% filter(cell_type == unique(mRNA_meta$cell_type)[i])
    data_miss_post_collapse_bas <- data_miss_post_collapse_bas %>% filter(diff > 1)
    data_miss_post_collapse_bas <- data_miss_post_collapse_bas %>% filter(Protein %in% rownames(p_all_min))

    for(j in 1:length(data_miss_post_collapse_bas$Protein)){
       prot_v <- p_all_min[data_miss_post_collapse_bas$Protein[j],intersect(colnames(p_all_min),um_plot_hold$ID)]
       prot_v[is.na(prot_v)][1:round(data_miss_post_collapse_bas$diff[j])] <- min(prot_v,na.rm=T)
       p_all_min[data_miss_post_collapse_bas$Protein[j],intersect(colnames(p_all_min),um_plot_hold$ID)] <- prot_v

    }

    
    um_plot_hold <- um_plot_hold %>% filter(ID %in% colnames(p_all))
    prot_i <- rowMedians(p_all_min[convert_good$split_prot,um_plot_hold$ID],na.rm = T)
    prot_na <- prot_i
    prot_na[is.na(prot_na)==F]  <- 1
    prot_na[rowSums(is.na(p_all_min[convert_good$split_prot,um_plot_hold$ID])==F) < 4] <- NA
    #prot_na_vect <- rowSums(is.na(p_all_min[convert_good$split_prot,rownames(um_plot_hold)])==F)#/nrow(um_plot_hold)
    prot_i = prot_i*prot_na
    #protein_mat_NA[,i] = prot_na_vect
    protein_mat[,i] = prot_i
    
    
    prot_i <- rowMedians(p_all_alpha[convert_good$split_prot,um_plot_hold$ID],na.rm = T)
    prot_na <- prot_i
    prot_na[is.na(prot_na)==F]  <- 1
    prot_na[rowSums(is.na(p_all_alpha[convert_good$split_prot,um_plot_hold$ID])==F) < 10] <- NA
    prot_i = prot_i*prot_na
    #protein_mat_NA[,i] = prot_na_vect
    deg_mat[,i] = prot_i

    
}

colnames(mRNA_mat) = unique(mRNA_meta$cell_type)
colnames(protein_mat) = unique(mRNA_meta$cell_type)
colnames(deg_mat) = unique(mRNA_meta$cell_type)

rownames(protein_mat) = convert_good$split_gene
rownames(mRNA_mat) <- rownames(rna_seq@assays$RNA@counts)
rownames(deg_mat) = convert_good$split_gene

ordering = c('Basal','Secratory','Cilliated','Chondrocyte','Fibroblast','Smooth muscle','Immune')
mRNA_mat = mRNA_mat[,ordering]
protein_mat = protein_mat[,ordering]
deg_mat = deg_mat[,ordering]

protein_mat[is.na(deg_mat)==T] <- NA

mRNA_mat_ <- mRNA_mat
mRNA_mat <- mRNA_mat[rownames(protein_mat),]

mRNA_mat[is.na(protein_mat)==T] <- NA

mRNA_mat <- Normalize_reference_vector(mRNA_mat,log = T)
mRNA_mat <- Normalize_reference_vector_log(mRNA_mat)
protein_mat[protein_mat==-Inf] <- NA
protein_mat <- Normalize_reference_vector_log(protein_mat)
protein_mat[protein_mat==-Inf] <- NA

deg_mat <- Normalize_reference_vector_log(deg_mat)

trans_mat <- protein_mat + deg_mat - mRNA_mat



write.csv(protein_mat,paste0(path_dat,'06_Gene_X_CellType/Relative_abundance/protein_freq.csv'))
write.csv(mRNA_mat,paste0(path_dat,'06_Gene_X_CellType/Relative_abundance/mRNA_freq.csv'))
write.csv(trans_mat,paste0(path_dat,'06_Gene_X_CellType/Relative_abundance/translation_freq.csv'))
write.csv(deg_mat,paste0(path_dat,'06_Gene_X_CellType/Relative_abundance/clearance_freq.csv'))




cor_mat <- cor(
  protein_mat, mRNA_mat,
  use = "pairwise.complete.obs",
  method = "pearson"
)




# 2) Create the heatmap, printing correlation values in each square
Heatmap(
  cor_mat,
  name = "corr",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
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

### Making scatter plots

### Basal cells
df_look <- data.frame(mrna = mRNA_mat[,1],prot = protein_mat[,1],deg = deg_mat[,2],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,1],prot = data_miss_obs_test[sect,1]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle('Basal') + xlim(c(-5,5))  + ylim(c(-5,5))

### Secratory cells
df_look <- data.frame(mrna = mRNA_mat[,2],prot = protein_mat[,2],deg = deg_mat[,2],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,2],prot = data_miss_obs_test[sect,2]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle(colnames(mRNA_mat)[2]) + xlim(c(-5,5))  + ylim(c(-5,5))

### Cilliated cells
df_look <- data.frame(mrna = mRNA_mat[,3],prot = protein_mat[,3],deg = deg_mat[,3],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,3],prot = data_miss_obs_test[sect,3]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  #scale_color_gradient2(midpoint = median(df_look$deg,na.rm=T),high='red',low = 'blue')+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle(colnames(mRNA_mat)[3]) + xlim(c(-5,5))  + ylim(c(-5,5))


### Chondrocyte cells
df_look <- data.frame(mrna = mRNA_mat[,4],prot = protein_mat[,4],deg = deg_mat[,4],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,4],prot = data_miss_obs_test[sect,4]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle(colnames(mRNA_mat)[4]) + xlim(c(-5,5))  + ylim(c(-5,5))


### Fibroblast cells
df_look <- data.frame(mrna = mRNA_mat[,5],prot = protein_mat[,5],deg = deg_mat[,5],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,5],prot = data_miss_obs_test[sect,5]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  #scale_color_gradient2(midpoint = median(df_look$deg,na.rm=T),high='red',low = 'blue')+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle(colnames(mRNA_mat)[5]) + xlim(c(-5,5))  + ylim(c(-5,5))


### Smooth muscle cells
df_look <- data.frame(mrna = mRNA_mat[,6],prot = protein_mat[,6],deg = deg_mat[,6],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,6],prot = data_miss_obs_test[sect,6]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle(colnames(mRNA_mat)[6]) + xlim(c(-5,5))  + ylim(c(-5,5))


### Immune cells
df_look <- data.frame(mrna = mRNA_mat[,7],prot = protein_mat[,7],deg = deg_mat[,7],
                      gene = convert_good$split_prot)
#df_look <- data.frame(mrna = mRNA_mat[sect,7],prot = data_miss_obs_test[sect,7]) 

ggplot(df_look,aes(x = mrna,y = prot)) + geom_point()+
  geom_abline(intercept = 0,slope = 1) + dot_plot + xlab('mRNA') + ylab('Protein')+
  ggtitle(colnames(mRNA_mat)[7]) + xlim(c(-5,5))  + ylim(c(-5,5))





#### Comparing Cell type slopes

slopes <- c()
cors <- c()
for(i in 1:7){
  slope_ <- TLS(mRNA_mat[,i],protein_mat[,i])[[1]]
  #mRNA_mat[,i] = mRNA_mat[,i]/slope_
  slopes <- c(slopes,slope_)
}

#before <- slopes
after <- slopes

df_plot <- data.frame(slopes = c(after,1/before),correct = c(rep('Post',7),rep(' Pre',7)))
ggplot(df_plot, aes(x = correct,y = slopes))+ geom_boxplot(outliers = F)+geom_jitter(width = .2)  + ylim(c(.25,1.25))+
  ggtitle('Across genes') + dot_plot + ylab('# Cell types') + xlab('')






#### Comparing gene slopes


protein_mat_ <- protein_mat[rowSums(is.na(protein_mat)==F) > 4,]
mRNA_mat__ <- mRNA_mat[rowSums(is.na(mRNA_mat)==F) > 4,]
deg_mat_ <- deg_mat[rowSums(is.na(deg_mat)==F) > 4,]
dim(protein_mat_)
sect <- intersect(intersect(rownames(mRNA_mat__),rownames(protein_mat_)),rownames(deg_mat_))
mRNA_mat__ <- mRNA_mat[sect,]
protein_mat_ <- protein_mat_[sect,]
deg_mat_ <- deg_mat[sect,]


trans_eff_mat = protein_mat_ + deg_mat_  - mRNA_mat__

cor_store <- c()
cor_store_prot_deg <- c()
cor_store_trans_eff <- c()
var_store <- c()
slopes_gene <- c()
na_deg <- c()
for(i in 1:nrow(mRNA_mat__)){
  cor_store <- c(cor_store,cor(mRNA_mat__[i,],protein_mat_[i,],use = 'pairwise.complete.obs'))
  slopes_gene <- c(slopes_gene,TLS(mRNA_mat__[i,],protein_mat_[i,])[[1]])
  #cor_store_prot_deg <- c(cor_store_prot_deg,cor(deg_mat_[i,],protein_mat_[i,],use = 'pairwise.complete.obs'))
  #cor_store_trans_eff <- c(cor_store_trans_eff,cor(trans_eff_mat[i,],protein_mat_[i,],use = 'pairwise.complete.obs'))
  #var_avg = (mean(abs(mRNA_mat__[i,]),na.rm=T)+mean(abs(protein_mat_[i,]),na.rm=T))/2
  var_avg = mean(abs(deg_mat_[i,]),na.rm=T)
  var_store <- c(var_store,var_avg)

  na_deg <- c(na_deg,sum(is.na(deg_mat_[i,])==F))
}

plot(log2(before),log2(after))
abline(a = 0,b = 1)

before <- slopes_gene
#after <- slopes_gene

df_plot_gene <- data.frame(slopes = c(after,1/before),correct = c(rep('Post',length(slopes_gene)),rep(' Pre',length(slopes_gene))))
ggplot(df_plot_gene, aes(fill = correct,x = slopes))+ geom_density(alpha = .2)  +
  ggtitle('Across genes') + dot_plot + ylab('Genes') + xlab('Slopes') +scale_x_log10(limits = c(0.1, 3))


View(df_plot_gene)

hist(log2(slopes_gene))








plt_gene <- data.frame(prot = protein_mat['Cavin1',],deg = -deg_mat['Cavin1',],mRNA = mRNA_mat['Cavin1',],cell_type = colnames(mRNA_mat))
plt_long <- plt_gene %>%
  pivot_longer(
    cols      = c(prot, deg, mRNA),
    names_to  = "modality",
    values_to = "abundance"
  ) %>%
  mutate(
    modality   = recode(modality, prot = "Protein", deg = "Deg",mRNA = 'mRNA'),
    cell_type  = factor(cell_type, levels = unique(cell_type))   # preserve order
  )

# 2) dot‑and‑line plot
ggplot(plt_long, aes(x = cell_type,
                     y = abundance,
                     group = modality,
                     colour = modality)) +
  geom_line() +
  geom_point(size = 4) +
  labs(
    x = "Cell type",
    y = "Relative abundance, log2",
    colour = NULL,
    title = "Hspa9"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+scale_color_manual(values = c('blue','#000','#ff7f0e'))+xlab('')

View(cor_store_df)
View(cor_store_df)
cor_store_df = as.data.frame(cor_store)
cor_store_df$slopes <- slopes_gene
cor_store_df$prot_deg <- cor_store_prot_deg
cor_store_df$trans_eff <- cor_store_trans_eff
cor_store_df$var = var_store
cor_store_df$deg_n <- na_deg
cor_store_df$prot = rownames(mRNA_mat__)
ggplot(cor_store_df, aes(x = cor_store)) + geom_histogram()+ dot_plot+
  xlab('Protein X mRNA') +ylab('# Genes') + ggtitle('Accross cluster correlations')+
  geom_vline(xintercept = median(cor_store_df$cor_store))

median(cor_store_df$cor_store)

View(cor_store_df)

ggplot(cor_store_df, aes(x = log2(slopes_gene))) + geom_histogram()+ dot_plot+
  xlab('Protein X mRNA') +ylab('# Genes') + ggtitle('Accross cluster log2(slopes)')


cor_store_df <- cor_store_df %>% filter(deg_n > 2)
#cor_store_df <- cor_store_df %>% filter(var > .5)

ggplot(cor_store_df, aes(x = prot_deg)) + geom_histogram()+ dot_plot+
  xlab('Protein abundance X Degradation') +ylab('# Genes') +  ggtitle('Accross cluster correlations')+
  geom_vline(xintercept = median(cor_store_df$prot_deg)-.1)

mean(cor_store_df$prot_deg)

ggplot(cor_store_df, aes(x = trans_eff)) + geom_histogram()+ dot_plot+
  xlab('Protein X Translation/mRNA') +ylab('# Genes') +  ggtitle('Accross cluster correlations')


View(cor_store_df)

cor_store_df <- cor_store_df %>% filter(is.na(prot_deg) == F)
cor_store_df$prot_deg <- round(cor_store_df$prot_deg/4,digits = 1)*4
cor_store_df$prot_deg <- factor(cor_store_df$prot_deg, levels = sort(unique(cor_store_df$prot_deg)))

#cor_store_df <- cor_store_df %>% filter(var > .5)
ggplot(cor_store_df, aes(x = prot_deg, y = cor_store)) +
  geom_boxplot() +
  stat_summary(fun.data = function(x) data.frame(y = max(x), label = length(x)),
               geom = "text", vjust = -0.5, color = "blue") +
  xlab('Protein X Deg. rate') +
  ylab('Protein X mRNA') +
  dot_plot  + ggtitle('Correlations')




cor_store_df$trans_eff <- round(cor_store_df$trans_eff/4,digits = 1)*4
cor_store_df$trans_eff <- factor(cor_store_df$trans_eff, levels = sort(unique(cor_store_df$trans_eff)))

ggplot(cor_store_df, aes(x = trans_eff, y = cor_store)) +
  geom_boxplot() +
  stat_summary(fun.data = function(x) data.frame(y = max(x), label = length(x)),
               geom = "text", vjust = -0.5, color = "blue") +
  xlab('Protein X Translation efficiency') +
  ylab('Protein X mRNA') +
  dot_plot  + ggtitle('Correlations')



