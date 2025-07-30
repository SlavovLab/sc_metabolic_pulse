# Required libraries
library(ggplot2)
library(dplyr)
library(Hmisc)
library(sva)
library(Seurat)
library(stringr)
library(QuantQC)
library(ComplexHeatmap)
library(circlize)

### Functions for analysis
findMarker_LikeX <- function(X,mat){
  
  cor_mat = cor(t(mat),use = 'pairwise.complete.obs')
  obs_mat = psych::pairwiseCount(t(mat))
  
  numb_obs <- rowSums(is.na(mat)==F)
  numb_obs <- outer(numb_obs, numb_obs, function(x, y) (x + y) / 2)
  
  new_metric = obs_mat/numb_obs
  
  df_out <- as.data.frame(cor_mat[,X])
  colnames(df_out) = 'Cor'
  df_out$Count = new_metric[,X]
  return(df_out)
  
}

Plot_markers <- function(prot_n,gene_n){
  
  um_plot$prot <- log10(p_all_abs[prot_n,])
  
  um_plot$prot[is.na(um_plot$prot)] <- min(um_plot$prot,na.rm=T)
  ggplot(um_plot, aes(x = Cell_Type,y = prot,color = prot)) + xlab('')+ ylab('Protein conc. per cell')+ 
    geom_boxplot(outliers = F)+ theme_classic()  + theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))
  
  
  mRNA_umap$gene <- rna_seq@assays$RNA@counts[gene_n,rownames(mRNA_umap)]
  ggplot(mRNA_umap, aes(x = cell_type, y = gene)) + ylab('# mRNA count per cell')+xlab('')+
    geom_boxplot(outliers = F)+ theme_classic() + theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))
  

  ## ── 1.  Protein plot (no x‑axis ticks or labels) ──────────────────────
  p_prot <- ggplot(um_plot,
                   aes(x = Cell_Type, y = prot, colour = prot)) +
    geom_boxplot(outliers = F) +              # hide outlier points
    ylab("Protein conc. per cell") +
    xlab(NULL) +
    theme_classic() +
    theme(
      axis.text.x  = element_blank(),   # remove tick text
      axis.ticks.x = element_blank()    # remove tick marks
    )+
    ggtitle(gene_n)
  
  ## ── 2.  mRNA plot (keep x‑axis ticks & angled labels) ─────────────────
  p_mrna <- ggplot(mRNA_umap,
                   aes(x = cell_type, y = gene)) +
    geom_boxplot(outliers = F) +
    ylab("# mRNA count per cell") +
    xlab(NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  
  ## ── 3.  Stack the two plots vertically ───────────────────────────────
  return(p_prot / p_mrna)  
  
}

GetCorVect <- function(prot,mRNA){
  prot_cor <- rcorr(t(prot))$r
  rna_cor <- rcorr(t(mRNA))$r
  
  diag(prot_cor) <- NA
  diag(rna_cor) <- NA
  
  genes <- c()
  cors <- c()
  sd_diff <- c()
  
  for(i in 1:ncol(rna_cor)){
    genes <- c(genes,rownames(rna_cor)[i])
    cors <- c(cors,cor(prot_cor[,i],rna_cor[,i], use = 'pairwise.complete.obs'))
    sd_diff <- c(sd_diff,(sd(rna_cor[,i],na.rm = T)-sd(prot_cor,na.rm = T)))
  }
  
  df <- as.data.frame(cors)
  df$genes <- genes
  
  return(df)
  
  
}


### load protein data sets
path <- "/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/"

load(paste0(path,'03_QuantQC_objects/r1_5day_male.RData'))
load(paste0(path,'03_QuantQC_objects/r2_5day_female.RData'))
load(paste0(path,'03_QuantQC_objects/r3_10day_male.RData'))
load(paste0(path,'03_QuantQC_objects/r4_10day_female.RData'))

#### Add the path of the seurat object for mRNA seq data to the 
rna_seq <- readRDS('/Users/andrewleduc/Desktop/Senescense/seurat_integrated_filered_700_named.rds')



# mRNA data QC 

df_qc <- rna_seq@meta.data[, c("nFeature_RNA", "nCount_RNA")] |>
  mutate(cluster = Idents(rna_seq)) 

df_qc <- df_qc %>% filter(cluster %in% names(mRNA_CTs))
ggplot(df_qc, aes(x = cluster, y = nCount_RNA, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = .1, outlier.shape = NA) +
  theme_classic() + xlab('')+ ylab('# Unique transcripts quantified per cell')+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_y_log10()






##########################################################################################
## Show quant for different data sets
##########################################################################################


Day1_prot <- colSums(is.na(r1_5day_male@matricies@protein)==F)
Day1_prot <- as.data.frame(Day1_prot)
colnames(Day1_prot) <- 'Number_Proteins'
Day1_prot$cond <- 'Day1'

Day2_prot <- colSums(is.na(r2_5day_female@matricies@protein)==F)
Day2_prot <- as.data.frame(Day2_prot)
colnames(Day2_prot) <- 'Number_Proteins'
Day2_prot$cond <- 'Day2'

Day3_prot <- colSums(is.na(r3_10day_male@matricies@protein)==F)
Day3_prot <- as.data.frame(Day3_prot)
colnames(Day3_prot) <- 'Number_Proteins'
Day3_prot$cond <- 'Day3'

Day4_prot <- colSums(is.na(r4_10day_female@matricies@protein)==F)
Day4_prot <- as.data.frame(Day4_prot)
colnames(Day4_prot) <- 'Number_Proteins'
Day4_prot$cond <- 'Day4'

prot_all = rbind(Day1_prot,Day2_prot,Day3_prot)#Day4_prot)

ggplot(prot_all,aes(x = cond,y = Number_Proteins)) +  geom_boxplot() +
  dot_plot + ylab('# proteins / cell') + xlab('')



get_ids <- r1_5day_male@matricies@peptide_protein_map %>% filter(Protein == 'P21981')

df_pep_scatter <- tibble(
  peptide_1 = r1_5day_male@matricies@peptide[get_ids$seqcharge[1], ],
  peptide_2 = r1_5day_male@matricies@peptide[get_ids$seqcharge[2], ]
)

ggplot(df_pep_scatter, aes(x = peptide_1, y = peptide_2)) +
  geom_point(size = 2, alpha = 0.7) +                 # plain dots
  labs(
    x = paste0("Concentration: ", get_ids$seqcharge[1]),
    y = paste0("Concentration: ", get_ids$seqcharge[2]),
    title = "Peptides from protein Tgm2"
  ) +
  theme_classic(base_size = 14)

col_fun = colorRamp2(c(0, 1), c("white", "red"))
Heatmap(cor(t(r1_5day_male@matricies@peptide[get_ids$seqcharge, ]),use = 'pairwise.complete.obs'),
        col = col_fun,show_row_dend = F,show_column_dend = F)


Quant1 <- r1_5day_male@pep.cor[[1]]
Quant1$cond = 'Day1'
Quant2 <- r2_5day_female@pep.cor[[1]]
Quant2$cond = 'Day2'
Quant3 <- r3_10day_male@pep.cor[[1]]
Quant3$cond = 'Day3'
Quant4 <- r4_10day_female@pep.cor[[1]]
Quant4$cond = 'Day4'

Quant_all <- rbind(Quant1,Quant2,Quant3,Quant4)

ggplot(Quant_all,aes(y = Cor, x = FC,fill = cond)) + geom_boxplot()+ dot_plot +
  ylab('Cor, peptides from same protein') + xlab('Fold change across single cells')





##########################################################################################
## Save matricies and integrate
##########################################################################################



p1 <- r1_5day_male@matricies@protein.imputed
p2 <- r2_5day_female@matricies@protein.imputed
p3 <- r3_10day_male@matricies@protein.imputed
p4 <- r4_10day_female@matricies@protein.imputed

p1_no <- r1_5day_male@matricies@protein
p2_no <- r2_5day_female@matricies@protein
p3_no <- r3_10day_male@matricies@protein
p4_no <- r4_10day_female@matricies@protein


p1_prot_abs <-  r1_5day_male@matricies@protein_abs
p2_prot_abs  <-  r2_5day_female@matricies@protein_abs
p3_prot_abs  <-  r3_10day_male@matricies@protein_abs
p4_prot_abs  <-  r4_10day_female@matricies@protein_abs


sect <- intersect(intersect(intersect(rownames(p1),rownames(p2)),rownames(p3)),rownames(p4))

p1 <- p1[sect,]
colnames(p1) <- paste0(colnames(p1),'_prep1')
p2 <- p2[sect,]
colnames(p2) <- paste0(colnames(p2),'_prep2')
p3 <- p3[sect,]
colnames(p3) <- paste0(colnames(p3),'_prep3')
p4 <- p4[sect,]
colnames(p4) <- paste0(colnames(p4),'_prep4')



## Make cor vects
imputed_cor_vect2 <- GetCorVect(p1,p2)
imputed_cor_vect2 <- imputed_cor_vect2 %>% filter(cors > .6)

# imputed_cor_vect1_2 <- GetCorVect(p1,p2)
# imputed_cor_vect1_2$comp = '1 vs 2'
# imputed_cor_vect1_3 <- GetCorVect(p1,p3)
# imputed_cor_vect1_3$comp = '1 vs 3'
# imputed_cor_vect1_4 <- GetCorVect(p1,p4)
# imputed_cor_vect1_4$comp = '1 vs 4'
# imputed_cor_vect2_3 <- GetCorVect(p2,p3)
# imputed_cor_vect2_3$comp = '2 vs 3'
# imputed_cor_vect2_4 <- GetCorVect(p2,p4)
# imputed_cor_vect2_4$comp = '2 vs 4'
# imputed_cor_vect3_4 <- GetCorVect(p3,p4)
# imputed_cor_vect3_4$comp = '3 vs 4'
# 
# comp_plot = rbind(imputed_cor_vect1_2,imputed_cor_vect1_3,imputed_cor_vect1_4,
#                   imputed_cor_vect2_3,imputed_cor_vect2_4,imputed_cor_vect3_4)
# 
# ggplot(comp_plot, aes(x = comp,y = cors))+geom_boxplot()+dot_plot+
#   ylim(c(-.2,1)) + ylab('Correlations')+ xlab('Pairwise dataset comparisons')



# take best cor vect proteins for integration

sect2 <- intersect(sect,imputed_cor_vect2$genes)

p1 <- as.data.frame(p1)
p1$prot <- rownames(p1)
p2 <- as.data.frame(p2)
p2$prot <- rownames(p2)
p3 <- as.data.frame(p3)
p3$prot <- rownames(p3)
p4 <- as.data.frame(p4)
p4$prot <- rownames(p4)

p_1and2 <- p1 %>% merge(p2, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3, by = 'prot',all = TRUE)
protein_Data <- p_123 %>% merge(p4, by = 'prot',all = TRUE)
rownames(protein_Data) <- protein_Data$prot
protein_Data$prot <- NULL
protein_Data <- as.matrix(protein_Data)
protein_Data <- protein_Data[sect2,]
dim(protein_Data)




### Save protein expression for all IDed proteins
colnames(p1_no) <- paste0(colnames(p1_no),'_prep1')
colnames(p2_no) <- paste0(colnames(p2_no),'_prep2')
colnames(p3_no) <- paste0(colnames(p3_no),'_prep3')
colnames(p4_no) <- paste0(colnames(p4_no),'_prep4')
p1_no <- as.data.frame(p1_no)
p1_no$prot <- rownames(p1_no)
p2_no <- as.data.frame(p2_no)
p2_no$prot <- rownames(p2_no)
p3_no <- as.data.frame(p3_no)
p3_no$prot <- rownames(p3_no)
p4_no <- as.data.frame(p4_no)
p4_no$prot <- rownames(p4_no)

p_1and2 <- p1_no %>% merge(p2_no, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_no, by = 'prot',all = TRUE)
p_all <- p_123 %>% merge(p4_no, by = 'prot',all = TRUE)
rownames(p_all) <- p_all$prot
p_all$prot <- NULL
p_all <- as.matrix(p_all)
p_all <- p_all[,colnames(protein_Data)]

write.csv(p_all,paste0(path,'04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv'))




### Save absolute abundances for all IDed proteins
colnames(p1_prot_abs) <- paste0(colnames(p1_prot_abs),'_prep1')
colnames(p2_prot_abs) <- paste0(colnames(p2_prot_abs),'_prep2')
colnames(p3_prot_abs) <- paste0(colnames(p3_prot_abs),'_prep3')
colnames(p4_prot_abs) <- paste0(colnames(p4_prot_abs),'_prep4')
p1_prot_abs <- as.data.frame(p1_prot_abs)
p1_prot_abs$prot <- rownames(p1_prot_abs)
p2_prot_abs <- as.data.frame(p2_prot_abs)
p2_prot_abs$prot <- rownames(p2_prot_abs)
p3_prot_abs <- as.data.frame(p3_prot_abs)
p3_prot_abs$prot <- rownames(p3_prot_abs)
p4_prot_abs <- as.data.frame(p4_prot_abs)
p4_prot_abs$prot <- rownames(p4_prot_abs)
p_1and2 <- p1_prot_abs %>% merge(p2_prot_abs, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_prot_abs, by = 'prot',all = TRUE)
p_all_abs <- p_123 %>% merge(p4_prot_abs, by = 'prot',all = TRUE)
rownames(p_all_abs) <- p_all_abs$prot
p_all_abs$prot <- NULL
p_all_abs <- as.matrix(p_all_abs)
p_all_abs <- p_all_abs[,colnames(protein_Data)]


write.csv(p_all_abs,paste0(path,'04_Gene_X_SingleCell_and_annotations/sc_protein_absolute.csv'))


##########################################################################################
## Integrate and assign cell types
##########################################################################################


protein_Data <- sva::ComBat(protein_Data,batch = (str_sub(colnames(protein_Data),-1)))
protein_Data <- Normalize_reference_vector_log(protein_Data)

prot_umap <- CreateSeuratObject(counts = protein_Data, project = "prot_mat")
prot_umap <- NormalizeData(prot_umap, normalization.method = "LogNormalize", scale.factor = 10000)
prot_umap@assays$RNA@layers$data <- protein_Data

all.genes <- rownames(protein_Data)
prot_umap <- ScaleData(prot_umap, features = all.genes)
prot_umap@assays$RNA@layers$data <- protein_Data
prot_umap <- Seurat::RunPCA(prot_umap, features = all.genes)
prot_umap <- FindNeighbors(prot_umap, dims = 1:6)
prot_umap <- FindClusters(prot_umap, resolution = 0.5)
prot_umap <- RunUMAP(prot_umap, dims = 1:6)

um_plot <- as.data.frame(prot_umap@reductions[["umap"]]@cell.embeddings)

um_plot$sample <- str_sub(colnames(protein_Data),-1)
um_plot$cluster <- prot_umap@meta.data[["RNA_snn_res.0.5"]]


# Plot by Batch
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = sample)) +
  geom_point(size = 1)+ theme_classic() +
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title = element_text(size = 16))

# Plot by Louvain cluster
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = cluster)) +
  geom_point(size = 1)+ theme_classic()


um_plot$Cell_Type = NA
um_plot$Cell_Type[um_plot$cluster == '12'] <- 'Smooth muscle'
um_plot$Cell_Type[um_plot$cluster %in% c('3','13')] <- 'Chondrocyte'
um_plot$Cell_Type[um_plot$cluster %in% c('10','4','6')] <- 'Fibroblast'
um_plot$Cell_Type[um_plot$cluster == '8'] <- 'Immune'
um_plot$Cell_Type[um_plot$cluster %in% c('9','5','1','0')] <- 'Basal'
um_plot$Cell_Type[um_plot$cluster %in% c('2','14')] <- 'Secratory'
um_plot$Cell_Type[um_plot$cluster %in% c('11')] <- 'Cilliated'
um_plot$Cell_Type[um_plot$cluster %in% c('7')] <- 'Shwann'
um_plot$Cell_Type[um_plot$Cell_Type == 'Cilliated' & um_plot$umap_1 < 1] <- 'Chondrocyte'


# plot by cell type now
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = Cell_Type)) +
  geom_point(size = 1)+ theme_classic() +
  theme(axis.text = element_blank(),axis.ticks = element_blank(),
        axis.title = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 3)))


##########################################################################################
## Get meta data
##########################################################################################

meta1 <- r1_5day_male@meta.data
meta1$prep = 'one'
meta1$prot_total_adj = meta1$prot_total - median(meta1$prot_total,na.rm=T)
meta1$ID <- paste0(meta1$ID,'_prep1')

meta2 <- r2_5day_female@meta.data
meta2$prep = 'two'
meta2$prot_total_adj = meta2$prot_total - median(meta2$prot_total,na.rm=T)
meta2$ID <- paste0(meta2$ID,'_prep2')

meta3 <- r3_10day_male@meta.data
meta3$prep = 'three'
meta3$prot_total_adj = meta3$prot_total - median(meta3$prot_total,na.rm=T)
meta3$ID <- paste0(meta3$ID,'_prep3')

meta4 <- r4_10day_female@meta.data
meta4$prep = 'four'
meta4$prot_total_adj = meta4$prot_total - median(meta4$prot_total,na.rm=T)
meta4$ID <- paste0(meta4$ID,'_prep4')

meta_mega <- rbind(meta1,meta2,meta3,meta4)

rownames(meta_mega) <- meta_mega$ID
meta_mega <- meta_mega[colnames(p_all),]

um_plot$age = meta_mega$sample

ggplot(meta_mega, aes(x = prep,fill = sample)) + geom_bar(position = 'dodge')


# Inspect diameters/size

um_plot$prot_total <- meta_mega$prot_total
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot_total)) +
   geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 25, low = 'blue',mid = 'white', high = 'red')

um_plot$diameter <- meta_mega$diameter
um_plot2 <- um_plot %>% group_by(sample) %>% mutate(prot_total_correct = prot_total - median(prot_total,na.rm=T))


ggplot(um_plot, aes(x = umap_1,y = umap_2,color = diameter)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = median(um_plot$diameter), low = 'blue',mid = 'white', high = 'red')

um_plot2 <- um_plot2 %>% filter(Cell_Type != 'Shwann')

palette <- c(
  Immune       = "#F90303",
  Basal        = "#B50202",
  Secratory    = "#B606C4",
  Cilliated    = "#9B70F9",
  Fibroblast   = "#2C78FF",
  `Smooth muscle` = "#0498BA",
  Chondrocyte  = "#03C1ED"
)

um_plot2$Cell_Type <- factor(um_plot2$Cell_Type, levels = names(palette))
ggplot(um_plot2, aes(y = prot_total_correct,x = log2(diameter^3),color = Cell_Type)) +
  geom_point(size = 1)+ theme_classic(base_size = 18) +ggtitle('pearson = 0.75')+
  xlab('Volume, um^3') +ylab('Total protein intensity')+
  scale_colour_manual(values = palette,guide  = guide_legend(override.aes = list(size = 5))) 

round(cor(um_plot2$prot_total_correct,log2(um_plot2$diameter^3),
    use = 'pairwise.complete.obs'),2)

meta_mega$diameter <- NULL
meta_mega$sample <- NULL
meta_mega$prot_total <- NULL

um_plot$ID <- rownames(um_plot)
um_plot <- um_plot %>% left_join(meta_mega,by = 'ID')

write.csv(um_plot,paste0(path,'04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'))

#####################################################################
# Collapse mRNA annotations
#####################################################################

unique(rna_seq@active.ident)

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

mRNA_umap = as.data.frame(rna_seq@reductions$umap@cell.embeddings)
mRNA_umap$cell_type = NA
mRNA_umap$CT = rna_seq@active.ident
mRNA_umap$cell_type[mRNA_umap$CT %in% Fibroblast] = 'Fibroblast'
mRNA_umap$cell_type[mRNA_umap$CT %in% Immune] = 'Immune'
mRNA_umap$cell_type[mRNA_umap$CT %in% Smooth_muscle] = 'Smooth muscle'
mRNA_umap$cell_type[mRNA_umap$CT %in% Cilliated] = 'Cilliated'
mRNA_umap$cell_type[mRNA_umap$CT %in% Chondrocyte] = 'Chondrocyte'
mRNA_umap$cell_type[mRNA_umap$CT %in% Secratory] = 'Secratory'
mRNA_umap$cell_type[mRNA_umap$CT %in% Basal] = 'Basal'
mRNA_umap$cell_type[mRNA_umap$CT %in% Shwann] = 'Shwann'


write.csv(mRNA_umap,paste0(path,'04_Gene_X_SingleCell_and_annotations/sc_mRNA_annotations.csv'))


#####################################################################
# # Final UMAPs
#####################################################################


mRNA_umap$cell_type <- factor(mRNA_umap$cell_type, levels = names(palette))
mRNA_umap <- mRNA_umap %>% filter(is.na(cell_type)==F)

# for the protein UMAP (note the capitalisation difference)
um_plot$Cell_Type   <- factor(um_plot$Cell_Type,   levels = names(palette))
um_plot <- um_plot %>% filter(is.na(Cell_Type)==F)


p1 <- ggplot(mRNA_umap,
             aes(UMAP_1, UMAP_2, colour = cell_type)) +
  geom_point(size = 0.7) +
  scale_colour_manual(values = palette,guide  = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.key.size = unit(5, "mm"))

p2 <- ggplot(um_plot,
             aes(umap_1, umap_2, colour = Cell_Type)) +
  geom_point(size = 0.7) +
  scale_colour_manual(values = palette,
                      guide  = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        legend.key.size = unit(5, "mm"))

p1
p2 




#####################################################################
# Marker plots /exploration  used to assign clusters 
#####################################################################



#Immune
um_plot$prot <- p_all['O89053',] # Q61233 O89053
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Coro1a')

Plot_markers('O89053','Coro1a')


#Basal
um_plot$prot <- p_all['Q922U2',]
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Krt5')

Plot_markers('Q922U2','Krt5')


#Sec cells
um_plot$prot <- p_all['O88312',]
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Agr2')

Plot_markers('O88312','Agr2')


# Cilliated
um_plot$prot <- p_all['O35685',rownames(um_plot)] # Q9WTM5, P60122
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

Cil_Markers  = findMarker_LikeX('O35685',p_all)
Cil_Markers <- Cil_Markers %>% filter(Cor > .4)
Cil_Markers <- Cil_Markers %>% filter(Count > .2)
um_plot$prot <- p_all['P60122',]
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Crocc2')
FeaturePlot(rna_seq,features = 'Ruvbl1')
FeaturePlot(rna_seq,features = 'Ruvbl2')

Plot_markers('P60122','Ruvbl1')


#Shwann cells
um_plot$prot <- p_all['P97361',] # O55103 P40240 P16330
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

um_plot$prot <- p_all['O89023',] # O55103 P40240 P16330
um_plot$prot <- p_all_alpha['O89023',] # O55103 P40240 P16330

ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Banf1')
FeaturePlot(rna_seq,features = 'Cnp')
FeaturePlot(rna_seq,features = 'Prx')
FeaturePlot(rna_seq,features = 'Lipa')


#Adipocytes
um_plot$prot <- p_all['P05784',]
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')
FeaturePlot(rna_seq,features = 'Fabp4')


# Chond
um_plot$prot <- p_all['P06802',] # Q60604, Q64669
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Enpp1')

Plot_markers('P06802','Enpp1')


#Smooth muscle
um_plot$prot <- p_all['O08638',] # O08638 P37804 A2ARA8
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Itga8')

Plot_markers('O08638','Itga8')

#Fibroblasts
um_plot$prot <- p_all['P40936',]  # P20152 P40936 P20152 O08553 P16045
ggplot(um_plot, aes(x = umap_1,y = umap_2,color = prot)) +
  geom_point(size = 1)+ theme_classic() +scale_color_gradient2(midpoint = 0, low = 'blue',mid = 'white', high = 'red')

FeaturePlot(rna_seq,features = 'Inmt')

Plot_markers('P40936','Inmt')


 

