library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(stringr)
library(seqinr)
library(Seurat)
library(QuantQC)
library(psych)
library(patchwork)
library(circlize)
library(reshape2)
library(STRINGdb)
library(igraph)
library(purrr)
library(matrixStats)
library(tidyr)
library(clusterProfiler) 
library(msigdbr)
library(org.Mm.eg.db);
library(AnnotationDbi)
library(GO.db)
library(Hmisc)


palette <- c(
  Immune       = "#F90303",
  Basal        = "#B50202",
  Secratory    = "#B606C4",
  Cilliated    = "#9B70F9",
  Fibroblast   = "#2C78FF",
  `Smooth muscle` = "#0498BA",
  Chondrocyte  = "#03C1ED"
)


## Minimal enrichment with GO BP + Hallmark + Reactome


##################
# data path
##################

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'

mRNA_raw_path <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/'

#########################
# Functions for analysis
#########################

get_interact <- function(my_genes){
  
  ############################################################
  ## 0 ·  Packages
  ############################################################
  
  
  
  # <--- replace with all 834
  species_id <- 10090                                   # 9606 = human, 10090 = mouse
  
  ## initialise STRING at “high confidence” (700 = 0.7)
  string <- STRINGdb$new(version = "12",
                         species = species_id,
                         score_threshold = 500)
  
  ############################  map & pull edges  ######################
  mapped <- string$map(data.frame(gene = my_genes),
                       "gene", removeUnmappedRows = TRUE)
  
  ids <- mapped$STRING_id
  edges_raw <- string$get_interactions(ids)
  
  ############################  keep only  (A) our genes AND (B) score≥700  ##
  edges_filt <- edges_raw |>
    filter(from %in% ids,
           to   %in% ids,
           combined_score >= 750)
  
  ############################  make clean, non-duplicated pair list  ########
  edges_df <- edges_filt |>
    mutate(gene1 = mapped$gene[match(from, mapped$STRING_id)],
           gene2 = mapped$gene[match(to,   mapped$STRING_id)]) |>
    ## sort the pair so “A–B” == “B–A”
    rowwise() |>
    mutate(pair = list(sort(c(gene1, gene2)))) |>
    mutate(gene1 = pair[[1]], gene2 = pair[[2]]) |>
    ungroup() |>
    dplyr::select(gene1, gene2) |>
    distinct()
  
  return(edges_df)
  
}

pair_df <- function(gene_vec, complex){
  gene_vec <- sort(unique(gene_vec))
  if(length(gene_vec) < 3) return(NULL)
  tibble(
    gene1   = combn(gene_vec, 2L)[1, ],
    gene2   = combn(gene_vec, 2L)[2, ],
    complex = complex
  )
}

normalize_within_prep <- function(mat,meta){
  #meta$ID <- rownames(meta)
  for(i in unique(meta$sample)){
    meta_hold <- meta %>% filter(sample == i)
    mat[,meta_hold$ID] <- Normalize_reference_vector_log(mat[,meta_hold$ID])
  }
  
  mat <- mat[rowSums(is.na(mat)==F)>100,]
  
  return(mat)
  
}

interaction_plot <- function(prot_mat,mrna_mat){
  ## ---------------------------------------------------------------------
  ## 1 · helper that returns all unique, order-standardised pairs
  keep_syms <- rownames(prot_mat)[!rownames(prot_mat) %in% c(ribo)]
  upper2mouse <- setNames(keep_syms, toupper(keep_syms))   # names = UPPER, value = mouse
  
  
  ## ---------------------------------------------------------------
  ## 1 ·  helper that produces all pairs (unchanged)
  ## ---------------------------------------------------------------
  
  
  ## ---------------------------------------------------------------
  ## 2 ·  loop over complexes — same as before
  ## ---------------------------------------------------------------
  pair_table <- complexes %>% 
    split(.$ComplexName) %>% 
    map_dfr(function(df){
      genes <- str_split(df$subunits.Gene.name., ";") %>% 
        unlist() %>% str_trim() %>% toupper()
      in_mat <- intersect(genes, names(upper2mouse))
      pair_df(in_mat, df$ComplexName[1])
    })
  
  ## ---------------------------------------------------------------
  ## 3 ·  convert the UPPER symbols back to mouse style
  ##     (NAs will appear only if a symbol was not present in the lookup)
  ## ---------------------------------------------------------------
  pair_table_mouse <- pair_table %>% 
    mutate(
      gene1 = upper2mouse[gene1],
      gene2 = upper2mouse[gene2]
    )
  
  
  prot_row_syms <- rownames(prot_mat)[!rownames(prot_mat) %in% c(ribo)]    # row-names to caps
  string_interact <- get_interact(prot_row_syms)
  
  string_interact$both <- paste0(string_interact$gene1,string_interact$gene2)
  string_interact <- string_interact %>% filter(!both %in% paste0(pair_table_mouse$gene1,pair_table_mouse$gene2))
  
  
  interactors <- rbind(string_interact[,1:2],pair_table_mouse[,1:2])
  interactors$both <- paste0(interactors$gene1,interactors$gene2)
  interactors <- interactors %>% distinct(both,.keep_all = T)
  
  #c('Aldh2','Tst','Ivd','Acsf2','Suclg2','Pycr2','Pcca','Pccb','Idh2','Echs1','Aldh4a1','Aldh6a1')
  
  cor_prot <- c()
  cor_mRNA <- c()
  pairwise_obs <- c()
  for(i in 1:nrow(interactors)){
    cor_prot <- c(cor_prot,cor(prot_mat[interactors$gene1[i],],prot_mat[interactors$gene2[i],],use = 'pairwise.complete.obs'))
    cor_mRNA <- c(cor_mRNA,cor(mrna_mat[interactors$gene1[i],],mrna_mat[interactors$gene2[i],]))
    pairwise_obs <- c(pairwise_obs, pairwiseCount(prot_mat[interactors$gene1[i],],prot_mat[interactors$gene2[i],]))
  }
  
  df_cov_string <- data.frame(prot = cor_prot,rna = cor_mRNA, count = pairwise_obs,gene1 =interactors$gene1,gene2 = interactors$gene2 )
  df_cov_string <-df_cov_string %>% dplyr::filter(count > 100)
  
  #ggplot(df_cov_string,aes(x = prot)) + geom_histogram() + geom_vline(xintercept = median(df_cov_string$prot))
  #ggplot(df_cov_string,aes(x = rna)) + geom_histogram() + geom_vline(xintercept = median(df_cov_string$rna))
  df_cov_string$count <- NULL
  df_cov_string2 <- melt(df_cov_string,ids = c('gene1','gene2'))
  
  return(df_cov_string2)
}

functional_heatmap_annotation <- function(mat1, mat2,gene_list){
  
  
  hm_name = 'test'

  
  mat1 <- mat1[gene_list,]
  mat2 <- mat2[gene_list,]
  
  cor_mat1 <- cor(t(mat1), method = "pearson",use = 'pairwise.complete.obs')
  cor_mat1[is.na(cor_mat1)]<-0
  cor_mat2 <- cor(t(mat2), method = "pearson",use = 'pairwise.complete.obs')
  cor_mat2[is.na(cor_mat2)]<-0
  
  k <- 10  # number of clusters you want
  
  # 1) Cluster and make a non-overlapping order
  hc <- hclust(as.dist(1 - cor_mat1), method = "average")
  cl <- cutree(hc, k = k)                                # names = genes
  ord_hc   <- hc$order
  genes_hc <- rownames(cor_mat1)[ord_hc]
  idx_by_cl <- split(seq_along(genes_hc), cl[genes_hc])
  clus_ids  <- sort(as.integer(names(idx_by_cl)))
  genes_ord <- unlist(lapply(clus_ids, function(i) genes_hc[idx_by_cl[[as.character(i)]]]),
                      use.names = FALSE)
  
  cor_ord <- cor_mat1[genes_ord, genes_ord, drop = FALSE]
  cl_fac  <- factor(cl[genes_ord], levels = clus_ids)
  

  
  # 3) Compute block indices and outline with thin black lines
  align_to <- split(seq_along(genes_ord), cl_fac)
  
  
  
  #### Get mRNA/Protein half half correlation heat map and plot boxes
  cor_mat1 <- cor_mat1[genes_ord,genes_ord]
  cor_mat2 <- cor_mat2[genes_ord,genes_ord]
  diag(cor_mat1) <- NA
  diag(cor_mat2) <- NA
  
  cor_mat1[upper.tri(cor_mat1,diag = T)] <- 0
  cor_mat2[lower.tri(cor_mat2,diag = T)] <- 0
  
  ht <- Heatmap(cor_mat1+cor_mat2,
          name = hm_name,
          cluster_rows = FALSE, cluster_columns = FALSE,   # keep our non-overlap order
          show_row_names = FALSE, show_column_names = FALSE)
  

  
  ### Enrichment for gene cluster modules
  
  genes_by_cluster <- split(genes_ord, cl_fac)
  

  
  ## Background/universe
  bg <- if (exists("bg_symbols")) bg_symbols else unique(unlist(genes_by_cluster))
  univ <- bitr(bg, "SYMBOL", "ENTREZID", org.Mm.eg.db) |> pull(ENTREZID) |> unique()
  
  ## Gene sets
  H  <- msigdbr("Mus musculus", category = "H") |>
    dplyr::select(term = gs_name, gene = entrez_gene)
  RE <- msigdbr("Mus musculus", category = "C2", subcategory = "REACTOME") |>
    dplyr::select(term = gs_name, gene = entrez_gene)
  
  ## GO BP TERM2GENE / TERM2NAME from OrgDb
  go_map <- AnnotationDbi::select(org.Mm.eg.db, keys = univ, keytype = "ENTREZID",
                                  columns = c("GOALL","ONTOLOGYALL")) |>
    filter(ONTOLOGYALL == "BP") |>
    transmute(term = GOALL, gene = ENTREZID) |>
    distinct()
  
  go_name <- AnnotationDbi::select(GO.db::GO.db, keys = unique(go_map$term),
                                   keytype = "GOID", columns = "TERM") |>
    transmute(term = GOID, name = TERM) |>
    distinct()
  
  enrich_one <- function(sym) {
    g <- bitr(sym, "SYMBOL", "ENTREZID", org.Mm.eg.db) |> pull(ENTREZID) |> unique()
    if (length(g) < 5) return(NULL)
    list(
      GO_BP = enricher(g, TERM2GENE = go_map, TERM2NAME = go_name,
                       universe = univ, pAdjustMethod = "BH",
                       minGSSize = 5, maxGSSize = 5000),
      Hallmark = enricher(g, TERM2GENE = H,
                          universe = univ, pAdjustMethod = "BH",
                          minGSSize = 5, maxGSSize = 5000),
      Reactome = enricher(g, TERM2GENE = RE,
                          universe = univ, pAdjustMethod = "BH",
                          minGSSize = 5, maxGSSize = 5000)
    )
  }
  
  res <- lapply(genes_by_cluster, enrich_one)
  
  top_terms <- lapply(names(res), function(cid) {
    x <- res[[cid]]; if (is.null(x)) return(NULL)
    bind_rows(
      if (!is.null(x$GO_BP))   as_tibble(x$GO_BP@result)   |> mutate(source = "GO_BP")   else NULL,
      if (!is.null(x$Hallmark))as_tibble(x$Hallmark@result)|> mutate(source = "Hallmark")else NULL,
      if (!is.null(x$Reactome))as_tibble(x$Reactome@result)|> mutate(source = "Reactome")else NULL
    ) |>
      filter(!is.na(p.adjust)) |>
      arrange(p.adjust) |>
      slice_head(n = 5) |>
      mutate(cluster = as.integer(cid))
  }) |> bind_rows()
  
  return(list(top_terms,ht,hm_name,clus_ids,align_to,genes_ord))
}

add_heatmap_clusters <- function(hm_info){
  
  hm_name <- hm_info[[3]]
  clus_ids <- hm_info[[4]]
  align_to <- hm_info[[5]]
  genes_ord <- hm_info[[6]]
  
  decorate_heatmap_body(hm_name, {
    n <- length(genes_ord)
    for (cid in clus_ids) {
      idx <- align_to[[cid]]
      y1 <- 1 - max(idx)/n; y2 <- 1 - (min(idx)-1)/n
      x1 <- (min(idx)-1)/n; x2 <- max(idx)/n
      # box (already drawn earlier, keep or remove)
      grid.rect(x = unit((x1 + x2)/2, "npc"),
                y = unit((y1 + y2)/2, "npc"),
                width = unit(x2 - x1, "npc"),
                height = unit(y2 - y1, "npc"),
                gp = gpar(fill = NA, col = "black", lwd = 0.8))
      # numeric label
      
    }
  })
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

heat_map_clust_avg <- function(mat1,CT){
  mat_make_fib <- matrix(data = NA,ncol = 2,nrow = 2)
  colnames(mat_make_fib) <- c('Purple','Green')
  rownames(mat_make_fib) <- c('Top','Bottom')
  
  mat_make_fib[1,1] <- cor(colMeans(mat1[convert_set1$split_gene,],na.rm = T),
                           colMeans(mat1[intersect(rownames(a2)[1:158],rownames(mat1)),],na.rm = T),use = 'pairwise.complete.obs')
  
  mat_make_fib[2,1] <- cor(colMeans(mat1[convert_set1$split_gene,],na.rm = T),
                           colMeans(mat1[intersect(rownames(a2)[159:nrow(a2)],rownames(mat1)),],na.rm = T),use = 'pairwise.complete.obs')
  
  
  mat_make_fib[1,2] <- cor(colMeans(mat1[convert_set2$split_gene,],na.rm = T),
                           colMeans(mat1[intersect(rownames(a2)[1:158],rownames(mat1)),],na.rm = T),use = 'pairwise.complete.obs')
  
  mat_make_fib[2,2] <- cor(colMeans(mat1[convert_set2$split_gene,],na.rm = T),
                           colMeans(mat1[intersect(rownames(a2)[159:nrow(a2)],rownames(mat1)),],na.rm = T),use = 'pairwise.complete.obs')
  
  Heatmap(mat_make_fib,col=col_fun,cluster_rows = F,cluster_columns = F,column_title = CT)
}




# Single cell analysis
ribo <- c("Rpl13" ,  "Rps11"  ,  "Rpl24"   ,  "Rps18"  ,  "Rpl6" ,
          "Rpl28" ,  "Rps26" ,  "Rps3" ,     "Rpl4"  ,  "Rps4x"  ,"Rpl21" ,   "Rpl10a",  "Rps23"  ,
          "Rpl18a",  "Rps21"  , "Rpl7"  ,  "Rps7"  ,  "Rps25"  , "Rps17"  , "Rpl35" ,  "Rpl3" ,   "Rpl5"  ,  "Rplp2" ,  "Rplp0"  ,
          "Rpl31" ,  "Rpl34"  ,  "Rps2"  ,  "Rpl23" ,  "Rpl23a" , "Rpl11" ,  "Rps8"  ,  "Rps6"  ,  "Rpl7a" ,  "Rpl38"  ,
          "Rpl27"  , "Rps16"  , "Rps19"  , "Rps9" ,   "Rpsa"  ,  "Rpl36a",  "Rpl9"  ,  "Rps10" ,  "Rpl13a",  "Rpl12" ,  "Rpl29" , 
          "Rpl22"  , "Rps15a" , "Rps20"  , "Rps13" ,  "Rps27"  , "Fau"  ,   "Rps12"  , "Rps24"  , "Rpl8" ,  
          "Rpl17"  , "Rps14",'Rpl26','Rpl27a','Rps5','Rpl30','Rpl19','Rpl35a','Rps28','Rpl18','Rpl14','Rps28')

mito <- c("Hadhb",   "Acadvl",  "Hadha",   "Acadm",   "Acads",   "Acadl",   "Acaa2",
          "Aco2",    "Cs",      "Idh3a" ,  "Idh3g"  , "Ogdh"  ,  "Dlst" ,   "Sucla2"  , "Suclg2",
          "Pdhb" ,    "Pdha1" ,"Mpc2", "Aldh6a1", 'Aldh4a1' ,"Pcca",    "Pccb",    "Etfa",
          "Etfb"  ,  "Etfdh"  , "Uqcrc1" , "Uqcrc2" ,  "Uqcrq"  , "Uqcrfs1", "Uqcrb",   "Sdhb",
          "Cox4i1",  "Cox6c",   "Cox7a1",  "Ndufa4",  "Ndufa5",  "Ndufs1",  "Slc25a20", "Letmd1",
          "Ucp1",    "Decr1",   "Ech1",    "Eci1",    "Gpd2",    "Prdx3",   "Park7",'Mdh2',
          "SOD2","Sod2", "Prdx5", "Hspa9", "Hspd1", "Hadh", "Echs1",
          "Dlat", "Sdha", "Acat1", "Idh2", "Ivd", "Tst",
          "Pycr2", "Acsf2","Coq9","Aldh2",'Cox7a2',
          "Ndufa10", "Ndufa12", "Ndufa2" , "Ndufa6" , "Ndufa8",  "Ndufa9" , "Ndufb10" ,"Ndufb11", "Ndufb5"  ,"Ndufc2" ,
          "Ndufs2" , "Ndufs3" , "Ndufs4" , "Ndufs5" , "Ndufs6" , "Ndufs7" , "Ndufs8" , "Ndufv2" ,
          "Coq8a", "Mtch2", "Cyc1", "Cox7c", "Cox5b", "Cox5a", "Ndufa13", "Ndufb4", "Ndufv1", "Mpc1", 
          "Slc25a1", "Slc25a4", "Slc25a5", "Slc25a3", "Slc25a11", "Cpt2", "Chchd3", "Immt", "Tufm", "Hibadh", "Dnajb7","Slc25a22", 
          "Lrpprc",  "Micos13", "Letm1",  "Lonp1", "Iars2",  "Cbr4",
          "Acad9",    "Vwa8",    "Cox6b1",  "Suclg1", "Mecr",  "Uqcrh",  "Hspe1",
          "Pdhx",     "Dld",     "Gcdh",    "Acad8",  "Hibch", "Phb2",  "Uqcr10",
          "Grpel1",   "Crat",    "Ndufv3",  "Clybl",  "Acot13", "Sepsecs", "Naxd",
          "Hint2",    "Bcat2",   "Slc25a10","Acadsb", "Mccc1", "Mmut",  "Auh",
          "Sfxn1",    "Eci2",    "Acad10",  "Nipsnap2","Hsd17b10","Acss1","Slc25a35",
          "Oat")




#########################
# Read in data for analysis
#########################


rna_seq <- readRDS(paste0(mRNA_raw_path,'seurat_integrated_filered_700_named.rds'))

p_all <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_relative.csv'),row.names = 1)
p_all <- as.matrix(p_all)
p_all_alpha <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/clearance_relative.csv'),row.names = 1)
p_all_alpha <- as.matrix(p_all_alpha)

p_all_alpha_abs <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/clearance_absolute.csv'),row.names = 1)
p_all_alpha_abs <- as.matrix(p_all_alpha_abs)

um_plot <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'),row.names = 1)
mRNA_meta <- read.csv(paste0(path_dat,'/04_Gene_X_SingleCell_and_annotations/sc_mRNA_annotations.csv'),row.names = 1)
mRNA_meta$ID <- rownames(mRNA_meta)
mRNA_meta <- mRNA_meta %>% filter(is.na(cell_type)==F)


convert = Proc_fasta(paste0(path_dat,'Mouse.fasta'))

convert = convert %>% filter(split_prot %in% rownames(p_all))
convert = convert %>% filter(split_gene %in% rownames(rna_seq@assays$RNA$counts))

convert_good$split_gene[1001:1455]

convert_good <- convert %>% filter(split_prot %in% rownames(p_all_alpha))


rna_seq <- NormalizeData(rna_seq)
rna_seq <- ScaleData(rna_seq,features = convert$split_gene,do.scale = F,
                     do.center = TRUE)

mRNA_norm <- as.matrix(rna_seq@assays$RNA@scale.data)

#### Extract data for the individual cell types 

# MRNA

mRNA_mat_bas = mRNA_norm[convert$split_gene,mRNA_meta$ID[mRNA_meta$cell_type == 'Basal']]

mRNA_mat_bas_count <- as.matrix( rna_seq@assays$RNA@counts)
mRNA_mat_bas_count = mRNA_mat_bas_count[convert$split_gene,mRNA_meta$ID[mRNA_meta$cell_type == 'Basal']]

mRNA_mat_fib = mRNA_norm[convert$split_gene,mRNA_meta$ID[mRNA_meta$cell_type == 'Fibroblast']]
mRNA_mat_chond = mRNA_norm[convert$split_gene,mRNA_meta$ID[mRNA_meta$cell_type == 'Chondrocyte']]
mRNA_mat_sec = mRNA_norm[convert$split_gene,mRNA_meta$ID[mRNA_meta$cell_type == 'Secratory']]

# Protein

um_plot_hold <- um_plot %>% filter(ID %in% colnames(p_all))


um_plot_bas = um_plot_hold %>% filter(Cell_Type == 'Basal')
um_plot_fib = um_plot_hold %>% filter(Cell_Type == 'Fibroblast')
um_plot_chond = um_plot_hold %>% filter(Cell_Type == 'Chondrocyte')
um_plot_sec = um_plot_hold %>% filter(Cell_Type == 'Secratory')

prot_mat_bas <- p_all[convert$split_prot,um_plot_bas$ID]
rownames(prot_mat_bas) <- convert$split_gene
prot_mat_fib <- p_all[convert$split_prot,um_plot_fib$ID]
rownames(prot_mat_fib) <- convert$split_gene
prot_mat_chond <- p_all[convert$split_prot,um_plot_chond$ID]
rownames(prot_mat_chond) <- convert$split_gene
prot_mat_sec <- p_all[convert$split_prot,um_plot_sec$ID]
rownames(prot_mat_sec) <- convert$split_gene


convert_deg <- convert %>% filter(split_prot %in% rownames(p_all_alpha))
deg_mat_bas <- p_all_alpha[convert_deg$split_prot,um_plot_bas$ID]
rownames(deg_mat_bas) <- convert_deg$split_gene
deg_mat_fib <- p_all_alpha[convert_deg$split_prot,um_plot_fib$ID]
rownames(deg_mat_fib) <- convert_deg$split_gene
deg_mat_chond <- p_all_alpha[convert_deg$split_prot,um_plot_chond$ID]
rownames(deg_mat_chond) <- convert_deg$split_gene
deg_mat_sec <- p_all_alpha[convert_deg$split_prot,um_plot_sec$ID]
rownames(deg_mat_sec) <- convert_deg$split_gene



prot_mat_bas <- normalize_within_prep(prot_mat_bas,um_plot_bas)
mRNA_mat_bas <- mRNA_mat_bas[rownames(prot_mat_bas),]
prot_mat_fib <- normalize_within_prep(prot_mat_fib,um_plot_fib)
mRNA_mat_fib <- mRNA_mat_fib[rownames(prot_mat_fib),]
prot_mat_chond <- normalize_within_prep(prot_mat_chond,um_plot_chond)
mRNA_mat_chond <- mRNA_mat_chond[rownames(prot_mat_chond),]
prot_mat_sec <- normalize_within_prep(prot_mat_sec,um_plot_sec)
mRNA_mat_sec <- mRNA_mat_sec[rownames(prot_mat_sec),]



df_sd <- data.frame(mRNA = rowSds(mRNA_mat_bas,na.rm=T),protein = rowSds(prot_mat_bas,na.rm=T))
ggplot(df_sd,aes(x = mRNA,y = protein)) + geom_point()+ xlim(c(0,2.2))+ ylim(c(0,2.2))+
  theme_classic(base_size = 18) + geom_abline(intercept = 0,slope = 1)+
  theme(
    axis.line       = element_blank(),                 # strip the default axes
    panel.border    = element_rect(                    # draw a new, even border
      colour = "black", fill = NA, linewidth = 2
    )
  ) + xlab('Proteins') + ylab('mRNAs') + ggtitle('Dynamic range (Sd)')

deg_mat_bas <- normalize_within_prep(deg_mat_bas,um_plot_bas)
deg_mat_fib <- normalize_within_prep(deg_mat_fib,um_plot_fib)
deg_mat_chond <- normalize_within_prep(deg_mat_chond,um_plot_chond)
deg_mat_sec <- normalize_within_prep(deg_mat_sec,um_plot_sec)




############--------############--------############--------############--------############--------
# Figures a-c

# Basal comparison
cor_mat_bas_prot <- cor(t(prot_mat_bas),use = 'pairwise.complete.obs')
cor_mat_bas_mRNA <- cor(t(mRNA_mat_bas))
obs_mat_bas_prot <- pairwiseCount(t(prot_mat_bas))
idx <- which(lower.tri(cor_mat_bas_prot, diag = FALSE), arr.ind = TRUE)
gene1 <- rownames(cor_mat_bas_prot)[idx[, 1]]   # first gene in the pair (row index)
gene2 <- colnames(cor_mat_bas_prot)[idx[, 2]]   # second gene (column index)

df_cor_basal <- data.frame(cor_prot = cor_mat_bas_prot[lower.tri(cor_mat_bas_prot)],
                           cor_mRNA = cor_mat_bas_mRNA[lower.tri(cor_mat_bas_mRNA)],
                           obs = obs_mat_bas_prot[lower.tri(obs_mat_bas_prot)],
                           gene1 = gene1,gene2 = gene2)
df_cor_basal <-df_cor_basal %>% filter(obs > 100 )
df_cor_basal <-df_cor_basal %>% filter(!gene1 %in% ribo)
df_cor_basal <-df_cor_basal %>% filter(!gene2 %in% ribo)   
df_cor_basal <-df_cor_basal %>% filter(!gene1 %in% mito)
df_cor_basal <-df_cor_basal %>% filter(!gene2 %in% mito)  

sum(df_cor_basal$cor_mRNA > .3)


df_cor_basal_sub <- df_cor_basal[sample(nrow(df_cor_basal),60000),]

cor(df_cor_basal_sub$cor_mRNA,df_cor_basal_sub$cor_prot,use = 'pairwise.complete.obs')

ggplot(df_cor_basal_sub,aes(x = cor_mRNA,y = cor_prot)) + ggpointdensity::geom_pointdensity()+
  theme_bw(base_size = 18) + geom_abline(intercept = 0,slope = 1) + xlim(c(-.6,.8))+ ylim(c(-.6,.8))



#### Explaining more significant protein protein correlations

df_cor_basal2 <- df_cor_basal %>% filter(abs(cor_prot) > .35)

df_cor_basal2_rna <- df_cor_basal %>% filter(abs(cor_mRNA) > .2)
length(unique(c(df_cor_basal2_rna$gene1,df_cor_basal2_rna$gene2)))

df_cor_basal2$deg_cor <- NA
df_cor_basal2$mRNA_avg <- NA


df_cor_basal2 <- df_cor_basal2 %>% filter(gene1 %in% rownames(deg_mat_bas))
df_cor_basal2 <- df_cor_basal2 %>% filter(gene2 %in% rownames(deg_mat_bas))

bas_count_df <- as.data.frame(rowSums(mRNA_mat_bas_count))
colnames(bas_count_df) <- 'Count'

for(i in 1:nrow(df_cor_basal2)){
  f1 <- NA;f2 <- NA
  f3 <- NA;f4 <- NA
  f5 <- NA;f6 <- NA
  #f7 <- NA;f8 <- NA
  
  if(pairwiseCount(prot_mat_bas[df_cor_basal2$gene1[i],],deg_mat_bas[df_cor_basal2$gene1[i],]) > 99){
    f1 <- cor(prot_mat_bas[df_cor_basal2$gene1[i],],deg_mat_bas[df_cor_basal2$gene1[i],],use = 'pairwise.complete.obs')
    f3 <- bas_count_df$Count[rownames(bas_count_df) == df_cor_basal2$gene1[i]]
    #f5 <- #df_prot_noise$sds[df_prot_noise$gene== gene1[i]]
    #f7
  }
  if(pairwiseCount(prot_mat_bas[df_cor_basal2$gene2[i],],deg_mat_bas[df_cor_basal2$gene2[i],]) > 99){
    f2 <- cor(prot_mat_bas[df_cor_basal2$gene2[i],],deg_mat_bas[df_cor_basal2$gene2[i],],use = 'pairwise.complete.obs')
    f4 <- bas_count_df$Count[rownames(bas_count_df) == df_cor_basal2$gene2[i]]
    #f6 <- #df_prot_noise$sds[df_prot_noise$gene== gene1[i]]
    #f8
  }

  df_cor_basal2$deg_cor[i] <-  mean(c(f1,f2),na.rm = T)
  df_cor_basal2$mRNA_avg[i] <- mean(c(f3,f4),na.rm = T)
  df_cor_basal2$prot_quant[i] <- mean(c(f5,f6),na.rm = T)
  
}
plot(df_cor_basal2$deg_cor,abs(df_cor_basal2$cor_prot -df_cor_basal2$cor_mRNA),
     ylab = "Protein mRNA correlation difference", xlab = 'Deg protein correlation',main = 'Cor = 0.43')
cor(abs(df_cor_basal2$cor_prot) - abs(df_cor_basal2$cor_mRNA),log2(df_cor_basal2$mRNA_avg),
    use = 'pairwise.complete.obs',method = 'pearson')

cor(df_cor_basal2$obs,abs(df_cor_basal2$cor_prot) - abs(df_cor_basal2$cor_mRNA))


#### Plot example indivudual proteins

plot(prot_mat_bas['Psma1',],deg_mat_bas['Psma1',])
plot(prot_mat_bas['Hmga1',],prot_mat_bas['Hmgn2',])
plot(sanity_basal['Hmga1',],sanity_basal['Hmgn2',])


p1 <- ggplot() +
  geom_point(aes(deg_mat_bas['Hmga1', ],prot_mat_bas['Hmga1', ]),alpha = .5) +
  labs(x = "Protein Hmga1", y = "Clearance Hmga1")+ theme_classic() +
  xlim(c(-3,2))

p2 <- ggplot() +
  geom_point(aes(prot_mat_bas['Hmga1', ], prot_mat_bas['Hmgn2', ]),alpha = .5) +
  labs(x = "Protein Hmga1", y = "Protein Hmgn2")+ theme_classic()+xlim(c(-4,2))+ylim(c(-2.5,2))

p3 <- ggplot() +
  geom_point(aes(sanity_basal['Hmga1', ], sanity_basal['Hmgn2', ]),alpha = .5) +
  labs(x = "Sanity Hmga1", y = "Sanity Hmgn2")+ theme_classic()

p4 <- ggplot() +
  geom_point(aes(deg_mat_bas['Hmgn2', ],prot_mat_bas['Hmgn2', ],),alpha = .5) +
  labs(y = "Clearance Hmgn2", x = "Protein Hmgn2")+ theme_classic()


(p2 + p1) / (p3 +p4)

plot(df_cor_basal2$HL1 + df_cor_basal2$HL2,abs(df_cor_basal2$cor_prot) - abs(df_cor_basal2$cor_mRNA),
     xlab = "Protein quant consistency, shared pep cor", ylab = 'Protein pairwise cor minus mRNA pairwise cor',main = 'Cor = 0')

plot(log2(df_cor_basal2$mRNA_avg),abs(df_cor_basal2$cor_prot) - abs(df_cor_basal2$cor_mRNA),
     ylab = "Protein mRNA correlation difference", xlab = 'Deg protein correlation',main = 'Cor = 0.35')

df_mod <- df_cor_basal2 %>%
  mutate(
    prot_mRNA_cor_diff = (abs(cor_prot) - abs(cor_mRNA)),   # response
    log2_mRNA_avg      = log2(mRNA_avg)                   # second predicto
  )

## 2.  Fit the multiple-linear model -----------------------------------
lm_fit <- lm(prot_mRNA_cor_diff ~  deg_cor + log2_mRNA_avg, data = df_mod)

## 3.  Inspect results --------------------------------------------------
summary(lm_fit)

df_mod$pred <- predict(lm_fit, newdata = df_mod)

# 2.  Scatter-plot predicted versus observed
ggplot(df_mod, aes(x = (prot_mRNA_cor_diff), y = pred)) +
  ggpointdensity::geom_pointdensity()+scale_color_viridis()+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(
    x = "Actual  (prot_mRNA_cor_diff)",
    y = "Predicted",
    title = "Linear model: predicting correlation difference",
    subtitle = 'Rsq = 0.30,R = 0.55'
  ) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none")



Bas_interactions<- interaction_plot(prot_mat_bas,sanity_basal)
Bas_interactions$cell_type <- 'Basal'


null_df <- melt(df_cor_basal %>% dplyr::select(cor_mRNA,cor_prot))
ggplot(null_df,aes(x = variable,y = value)) + geom_boxplot() + ylab('Correlations') + xlab('')+ 
  ggtitle('') + theme_classic(base_size = 18) + geom_hline(yintercept = 0)+
  ylim(c(-.5,.8))



#### Making heatmap mRNA protein correlation modules

df_cor_basal2 <- df_cor_basal %>% filter(abs(cor_prot) > .35)

df_cor_basal2_rna <- df_cor_basal %>% filter(abs(cor_mRNA) > .2)
length(unique(c(df_cor_basal2_rna$gene1,df_cor_basal2_rna$gene2)))

look <- cor(t(prot_mat_bas[unique(c(df_cor_basal2$gene1,df_cor_basal2$gene2)),]),use = 'pairwise.complete.obs')
cvm <- c()
for(i in 1:ncol(look)){
  cvm <- c(cvm,sum(abs(look[,i]),na.rm=T))
}
hist(cvm)
df_look <- data.frame(prot = colnames(look),mag = cvm)
df_look <- df_look %>% filter(mag > 50)
prot_mat_bas2 <- prot_mat_bas[rowSums(is.na(prot_mat_bas)==F) > 200,]
df_look <- df_look %>% filter(prot %in% rownames(prot_mat_bas2))


prot_cor_mod <- functional_heatmap_annotation(prot_mat_bas,mRNA_mat_bas,df_look$prot)
prot_cor_mod[[2]]
add_heatmap_clusters(prot_cor_mod)



rna_cor_mod <- functional_heatmap_annotation(mRNA_mat_bas,prot_mat_bas,unique(c(df_cor_basal2_rna$gene1,df_cor_basal2_rna$gene2)))
rna_cor_mod[[2]]
add_heatmap_clusters(rna_cor_mod)





############
## Secratory comparison
############

cor_mat_sec_prot <- cor(t(prot_mat_sec),use = 'pairwise.complete.obs')
cor_mat_sec_mRNA <- cor(t(mRNA_mat_sec))
obs_mat_sec_prot <- pairwiseCount(t(prot_mat_sec))
idx <- which(lower.tri(cor_mat_sec_prot, diag = FALSE), arr.ind = TRUE)
gene1 <- rownames(cor_mat_sec_prot)[idx[, 1]]   # first gene in the pair (row index)
gene2 <- colnames(cor_mat_sec_prot)[idx[, 2]]   # second gene (column index)

df_cor_sec <- data.frame(cor_prot = cor_mat_sec_prot[lower.tri(cor_mat_sec_prot)],
                           cor_mRNA = cor_mat_sec_mRNA[lower.tri(cor_mat_sec_mRNA)],
                           obs = obs_mat_sec_prot[lower.tri(obs_mat_sec_prot)],
                           gene1 = gene1,gene2 = gene2)
df_cor_sec <-df_cor_sec %>% filter(obs > 100 )
df_cor_sec <-df_cor_sec %>% filter(!gene1 %in% ribo)
df_cor_sec <-df_cor_sec %>% filter(!gene2 %in% ribo)  
df_cor_sec <-df_cor_sec %>% filter(!gene1 %in% mito)
df_cor_sec <-df_cor_sec %>% filter(!gene2 %in% mito) 

df_cor_sec2 <- df_cor_sec %>% filter(abs(cor_prot) > .25)
df_cor_sec2$deg_cor <- NA
df_cor_sec2 <- df_cor_sec2 %>% filter(gene1 %in% rownames(deg_mat_sec))
df_cor_sec2 <- df_cor_sec2 %>% filter(gene2 %in% rownames(deg_mat_sec))

for(i in 1:nrow(df_cor_sec2)){
  f1 <- NA;f2 <- NA
  if(pairwiseCount(prot_mat_sec[df_cor_sec2$gene1[i],],deg_mat_sec[df_cor_sec2$gene1[i],]) > 99){
    f1 <- cor(prot_mat_sec[df_cor_sec2$gene1[i],],deg_mat_sec[df_cor_sec2$gene1[i],],use = 'pairwise.complete.obs')
  }
  if(pairwiseCount(prot_mat_sec[df_cor_sec2$gene2[i],],deg_mat_sec[df_cor_sec2$gene2[i],]) > 99){
    f2 <- cor(prot_mat_sec[df_cor_sec2$gene2[i],],deg_mat_sec[df_cor_sec2$gene2[i],],use = 'pairwise.complete.obs')
  }
  
  df_cor_sec2$deg_cor[i] <-  mean(c(f1,f2),na.rm = T)
  
}
plot(df_cor_sec2$deg_cor,abs(df_cor_sec2$cor_prot) - abs(df_cor_sec2$cor_mRNA),
     ylab = "Protein mRNA correlation difference", xlab = 'Deg protein correlation',main = 'Cor = 0.21')
cor(abs(df_cor_sec2$cor_prot) - abs(df_cor_sec2$cor_mRNA),df_cor_sec2$deg_cor,
    use = 'pairwise.complete.obs',method = 'pearson')


cor(df_cor_sec$cor_mRNA,df_cor_sec$cor_prot)

Sec_interactions <- interaction_plot(prot_mat_sec,mRNA_mat_sec)
Sec_interactions$cell_type = 'Secretory'

deg_sec <- data.frame(HL1 = deg_abs_mat[,3],gene1 = rownames(deg_abs_mat))
df_cor_sec <- df_cor_sec %>% left_join(deg_sec, by = c('gene1'))
deg_sec <- data.frame(HL2 = deg_abs_mat[,3],gene2 = rownames(deg_abs_mat))
df_cor_sec <- df_cor_sec %>% left_join(deg_sec, by = c('gene2'))


####################################
# Fib comparison
####################################
cor_mat_fib_prot <- cor(t(prot_mat_fib),use = 'pairwise.complete.obs')
cor_mat_fib_mRNA <- cor(t(mRNA_mat_fib)) 
obs_mat_fib_prot <- pairwiseCount(t(prot_mat_fib))
idx <- which(lower.tri(cor_mat_fib_prot, diag = FALSE), arr.ind = TRUE)
gene1 <- rownames(cor_mat_fib_prot)[idx[, 1]]   # first gene in the pair (row index)
gene2 <- colnames(cor_mat_fib_prot)[idx[, 2]]   # second gene (column index)

df_cor_fib <- data.frame(cor_prot = cor_mat_fib_prot[lower.tri(cor_mat_fib_prot)],
                           cor_mRNA = cor_mat_fib_mRNA[lower.tri(cor_mat_fib_mRNA)],
                           obs = obs_mat_fib_prot[lower.tri(obs_mat_fib_prot)],
                           gene1 = gene1,gene2 = gene2)
df_cor_fib <- df_cor_fib %>% filter(obs > 100 )
df_cor_fib <- df_cor_fib %>% filter(!gene1 %in% ribo)
df_cor_fib <- df_cor_fib %>% filter(!gene2 %in% ribo)
df_cor_fib <- df_cor_fib %>% filter(!gene1 %in% mito)
df_cor_fib <- df_cor_fib %>% filter(!gene2 %in% mito)






df_cor_fib$dif <- df_cor_fib$cor_prot - df_cor_fib$cor_mRNA

cor(df_cor_fib$cor_mRNA,df_cor_fib$cor_prot,use = 'pairwise.complete.obs')

Fib_interactions <- interaction_plot(prot_mat_fib,mRNA_mat_fib)
Fib_interactions$cell_type = 'Fibroblast'

nrow(df_cor_fib %>% filter(cor_prot > .25 & dif > .1))
nrow(df_cor_fib %>% filter(cor_mRNA > .25 & dif < -.1))


deg_fib <- data.frame(HL1 = deg_abs_mat[,4],gene1 = rownames(deg_abs_mat))
df_cor_fib <- df_cor_fib %>% left_join(deg_fib, by = c('gene1'))
deg_fib <- data.frame(HL2 = deg_abs_mat[,4],gene2 = rownames(deg_abs_mat))
df_cor_fib <- df_cor_fib %>% left_join(deg_fib, by = c('gene2'))

df_cor_fib2 <- df_cor_fib %>% filter(abs(cor_prot) > .25)
df_cor_fib2$deg_cor <- NA
df_cor_fib2 <- df_cor_fib2 %>% filter(gene1 %in% rownames(deg_mat_fib))
df_cor_fib2 <- df_cor_fib2 %>% filter(gene2 %in% rownames(deg_mat_fib))

for(i in 1:nrow(df_cor_fib2)){
  f1 <- NA;f2 <- NA
  if(pairwiseCount(prot_mat_fib[df_cor_fib2$gene1[i],],deg_mat_fib[df_cor_fib2$gene1[i],]) > 99){
    f1 <- cor(prot_mat_fib[df_cor_fib2$gene1[i],],deg_mat_fib[df_cor_fib2$gene1[i],],use = 'pairwise.complete.obs',
              method = 'pearson')
  }
  if(pairwiseCount(prot_mat_fib[df_cor_fib2$gene2[i],],deg_mat_fib[df_cor_fib2$gene2[i],]) > 99){
    f2 <- cor(prot_mat_fib[df_cor_fib2$gene2[i],],deg_mat_fib[df_cor_fib2$gene2[i],],use = 'pairwise.complete.obs',
              method = 'pearson')
  }
  
  df_cor_fib2$deg_cor[i] <-  mean(c(f1,f2),na.rm = T)
  
}
plot(df_cor_fib2$deg_cor,abs(df_cor_fib2$cor_prot) - abs(df_cor_fib2$cor_mRNA),
     ylab = "Protein mRNA correlation difference", xlab = 'Deg protein correlation',main = 'Cor = 0.21')
cor(abs(df_cor_fib2$cor_prot) - abs(df_cor_fib2$cor_mRNA),df_cor_fib2$deg_cor,
    use = 'pairwise.complete.obs',method = 'pearson')



df_cor_fib <- df_cor_fib %>% filter(cor_prot > .3)

plot(df_cor_fib$dif,log2(df_cor_fib$HL1+df_cor_fib$HL2))

cor(df_cor_fib$dif,log2(df_cor_fib$HL1+df_cor_fib$HL2),use = 'pairwise.complete.obs')

################################################
# Chond comparison
################################################
cor_mat_chond_prot <- cor(t(prot_mat_chond),use = 'pairwise.complete.obs')
cor_mat_chond_mRNA <- cor(t(mRNA_mat_chond))
obs_mat_chond_prot <- pairwiseCount(t(prot_mat_chond))
idx <- which(lower.tri(cor_mat_chond_prot, diag = FALSE), arr.ind = TRUE)
gene1 <- rownames(cor_mat_chond_prot)[idx[, 1]]   # first gene in the pair (row index)
gene2 <- colnames(cor_mat_chond_prot)[idx[, 2]]   # second gene (column index)

df_cor_chond <- data.frame(cor_prot = cor_mat_chond_prot[lower.tri(cor_mat_chond_prot)],
                           cor_mRNA = cor_mat_chond_mRNA[lower.tri(cor_mat_chond_mRNA)],
                           obs = obs_mat_chond_prot[lower.tri(obs_mat_chond_prot)],
                           gene1 = gene1,gene2 = gene2)
df_cor_chond <-df_cor_chond %>% filter(obs > 100 )
df_cor_chond <-df_cor_chond %>% filter(!gene1 %in% ribo)
df_cor_chond <-df_cor_chond %>% filter(!gene2 %in% ribo)
df_cor_chond <-df_cor_chond %>% filter(!gene1 %in% mito)
df_cor_chond <-df_cor_chond %>% filter(!gene2 %in% mito)


Chond_interactions <- interaction_plot(prot_mat_chond,mRNA_mat_chond)
Chond_interactions$cell_type = 'Chondrocyte'

deg_mat_chond <- deg_mat_chond[rowSums(is.na(deg_mat_chond)==F) >100,]
cor_mat_chond_deg <- cor(t(deg_mat_chond), use = 'pairwise.complete.obs')
obs_mat_chond_deg <- pairwiseCount(t(deg_mat_chond))
idx <- which(lower.tri(cor_mat_chond_deg, diag = FALSE), arr.ind = TRUE)
gene1 <- rownames(cor_mat_chond_deg)[idx[, 1]]   # first gene in the pair (row index)
gene2 <- colnames(cor_mat_chond_deg)[idx[, 2]]   # second gene (column index)

df_chond_deg <- data.frame(deg_cor = cor_mat_chond_deg[lower.tri(cor_mat_chond_deg)],
                           obs_deg = obs_mat_chond_deg[lower.tri(obs_mat_chond_deg)],
                           gene1 = gene1,gene2 = gene2)

df_chond_deg <- df_chond_deg %>% filter(obs_deg > 100)

df_chond_deg <- df_chond_deg %>% left_join(df_cor_chond, by = c('gene1','gene2'))

cor(df_chond_deg$deg_cor,df_chond_deg$cor_prot,use = 'pairwise.complete.obs')

df_chond_deg_sub <- df_chond_deg[sample(nrow(df_chond_deg),50000),]

cor(df_chond_deg_sub$deg_cor,df_chond_deg_sub$cor_prot,use = 'pairwise.complete.obs')

ggplot(df_chond_deg_sub,aes(x = deg_cor,y = cor_prot)) + ggpointdensity::geom_pointdensity()+
  theme_bw(base_size = 18)




sanity_basal <- read.delim('/Users/andrewleduc/Desktop/Github/Sanity/bin/sanity_out/log_transcription_quotients_vmax.txt')
rownames(sanity_basal) <- rownames(prot_mat_bas)
sanity_basal$GeneID <- NULL
sanity_basal <- as.matrix(sanity_basal)
dim(sanity_basal)



colnames(test_t) <- paste0('Cell ', 1:ncol(test_t))

rownames(test_t) <-  paste0('Gene ', 1:nrow(test_t))


test_t <- rna_seq@assays$RNA@counts[convert$split_gene,rna_seq@active.ident %in% mRNA_CTs[['Basal']]]



df_rna_plot <- data.frame(
  Cav1 = sanity_basal["Cav1", ],
  Cavin1  = sanity_basal["Cavin1",  ]
)

a <- ggplot(df_rna_plot, aes(x = Cav1, y = Cavin1)) + geom_point(alpha = .3)+
  theme_classic(base_size = 15) + ggtitle('mRNA')

df_prot_plot <- data.frame(
  Cav1 = prot_mat_bas["Cav1", ],
  Cavin1  = prot_mat_bas["Cavin1",  ]
)

b <- ggplot(df_prot_plot, aes(x = Cav1, y = Cavin1)) + geom_point(alpha = .5)+
  theme_classic(base_size = 15) + ggtitle('Protein')




sd(mRNA_mat_bas["Cavin1", ] - mean(mRNA_mat_bas["Cavin1", ]))
sd(prot_mat_bas["Cavin1", ]-mean(prot_mat_bas["Cavin1", ],na.rm=T),na.rm=T)

df_rna_plot <- data.frame(
  value = c(mRNA_mat_bas["Cavin1", ] - mean(mRNA_mat_bas["Cavin1", ]),prot_mat_bas["Cavin1", ]-mean(prot_mat_bas["Cavin1", ],na.rm=T)) ,
  type  = c(rep('mRNA',ncol(mRNA_mat_bas)),rep('prot',ncol(prot_mat_bas))))

ggplot(df_rna_plot, aes(fill = type, x = value)) + geom_density(alpha = .3)+
  theme_classic(base_size = 15) + ggtitle('Cavin1')+ scale_fill_manual(values = c('#F7941D','#00A651'))+
  xlab('Concentrations') + ylab('Density, single cells')

library(dplyr)
library(ggplot2)

# build centered vectors (drop NAs)
v_m <- as.numeric(sanity_basal["Cavin1", ])
v_p <- as.numeric(prot_mat_bas["Cavin1", ])
df_rna_plot <- tibble(
  value = c(v_m - mean(v_m, na.rm = TRUE),
            v_p - mean(v_p, na.rm = TRUE)),
  type  = c(rep("mRNA", length(v_m)),
            rep("prot", length(v_p)))
) %>% filter(!is.na(value))

# rank within each type
df_rank <- df_rna_plot %>%
  group_by(type) %>%
  arrange(value, .by_group = TRUE) %>%
  mutate(rank = row_number(),
         rank_pct = rank / n()) %>%
  ungroup()

ggplot(df_rank, aes(x = rank_pct, y = value, color = type, group = type)) +
  geom_line(size = 1) +
  geom_point(size = 1, alpha = 0.7) +
  theme_classic(base_size = 15) +
  ggtitle("Cavin1") +
  scale_color_manual(values = c(mRNA = "#F7941D", prot = "#00A651")) +
  xlab("Rank percentile") + ylab("Log2 fold change")



plot(c(.27,.15,.21,.18),c(.03,0,.011,.016))

df_cor_comp <- data.frame(cor = c(.27,.15,.21,.18),
                          growth = c(.03,0,.017,.011),
                          tissue = c('Basal','Chondrocyte',
                                     'Secratory','Fibroblast'))

ggplot(df_cor_comp,aes(y = cor, x = growth,color = tissue)) + geom_point(size = 5) + theme_bw(base_size  = 16)



###### Correlations for degradation rates

deg_cor_comp <- function(deg_mat,prot_mat){
  sect <- intersect(rownames(deg_mat),rownames(prot_mat))
  deg_mat <- deg_mat[sect,]
  prot_mat <- prot_mat[sect,]
  
  cor_prot <- cor(t(prot_mat),use = 'pairwise.complete.obs')
  obs_prot <- pairwiseCount(t(prot_mat))
  cor_deg <- cor(t(deg_mat),use = 'pairwise.complete.obs')
  obs_deg <- pairwiseCount(t(deg_mat))
  idx <- which(lower.tri(cor_prot, diag = FALSE), arr.ind = TRUE)
  gene1 <- rownames(cor_prot)[idx[, 1]]   # first gene in the pair (row index)
  gene2 <- colnames(cor_prot)[idx[, 2]]   # second gene (column index)
  
  df_comp_cor <- data.frame(cor_deg = cor_deg[lower.tri(cor_deg)],
                            obs_deg = obs_deg[lower.tri(obs_deg)],
                            cor_prot = cor_prot[lower.tri(cor_prot)],
                            obs_prot = obs_prot[lower.tri(obs_prot)],
                             gene1 = gene1,gene2 = gene2)
  
  df_comp_cor <- df_comp_cor %>% filter(obs_prot > 100)
  df_comp_cor <- df_comp_cor %>% filter(obs_deg > 100)
  
  return(df_comp_cor)
  
}

deg_mat_chond <- deg_mat_chond[rowSums(is.na(deg_mat_chond)==F) >100,]
df_chond_deg <- deg_cor_comp(deg_mat_chond,prot_mat_chond)

ggplot(df_chond_deg,aes(x = cor_deg,y = cor_prot)) + ggpointdensity::geom_pointdensity(alpha = .2)+
  theme_bw(base_size = 18) + geom_abline(intercept = 0,slope = 1) + xlim(c(-.7,1)) +  ylim(c(-.7,1))

cor(df_chond_deg$cor_deg,df_chond_deg$cor_prot,use = 'pairwise.complete.obs')

interact_deg_chond <- interaction_plot_deg(deg_mat_chond)
interact_deg_chond$cell_type <- 'Chondrocyte'

plot(df_chond_deg$cor_deg,df_chond_deg$cor_prot,xlim = c(-.7,1),ylim = c(-.7,1))

deg_mat_fib <- deg_mat_fib[rowSums(is.na(deg_mat_fib)==F) >100,]
df_fib_deg <- deg_cor_comp(deg_mat_fib,prot_mat_fib)

cor(df_fib_deg$cor_deg,df_fib_deg$cor_prot,use = 'pairwise.complete.obs')

interact_deg_fib <- interaction_plot_deg(deg_mat_fib)
interact_deg_fib$cell_type <- 'Fibroblast'

deg_mat_sec <- deg_mat_sec[rowSums(is.na(deg_mat_sec)==F) >100,]
df_sec_deg <- deg_cor_comp(deg_mat_sec,prot_mat_sec)

cor(df_sec_deg$cor_deg,df_sec_deg$cor_prot,use = 'pairwise.complete.obs')

interact_deg_sec <- interaction_plot_deg(deg_mat_sec)
interact_deg_sec$cell_type <- 'Secratory'

deg_mat_bas <- deg_mat_bas[rowSums(is.na(deg_mat_bas)==F) >100,]
df_bas_deg <- deg_cor_comp(deg_mat_bas,prot_mat_bas)

cor(df_bas_deg$cor_deg,df_bas_deg$cor_prot,use = 'pairwise.complete.obs')


interact_deg_bas <- interaction_plot_deg(deg_mat_bas)
interact_deg_bas$cell_type <- 'Basal'

df_cor_comp <- data.frame(cor = c(0.08,0.2504645,0.133848,0.1226482),
                          growth = c(.03,0,.017,.011),
                          tissue = c('Basal','Chondrocyte','Secratory','Fibroblast'))

ggplot(df_cor_comp,aes(y = cor, x = growth,color = tissue)) + 
  geom_point(size = 5) + theme_bw(base_size  = 16) + xlab('Growth rate, days^-1')+
  ylab('cor(Clearance,Protein) pairwise correlations')




#### Make interaction plot for all cell types


interact_all <- rbind(Bas_interactions,Sec_interactions,Chond_interactions,Fib_interactions,
      interact_deg_fib,interact_deg_chond,interact_deg_sec,interact_deg_bas )


interact_all <- interact_all %>% filter(variable != 'deg')

ggplot(interact_all,aes(x = variable,y = value)) + geom_boxplot()+
  facet_wrap(~cell_type,ncol=4) + theme_classic(base_size = 18)+
  xlab('') + geom_hline(yintercept = 0) + ylab('Correlations') + ylim(c(-.7,.7))


Null_dist <- rbind(data.frame(cor = df_cor_basal$cor_prot,type = 'protein'),
      data.frame(cor = df_cor_basal$cor_mRNA,type = 'rna'),
      data.frame(cor = df_bas_deg$cor_deg,type = 'zdeg'))

Null_dist <- Null_dist %>% filter(type != 'zdeg')

ggplot(Null_dist,aes(x = type,y = cor)) + geom_boxplot()+
  theme_classic(base_size = 18)+
  xlab('') + geom_hline(yintercept = 0) + ylab('Correlations')+ ylim(c(-.7,.7))





####################################################################################
# Go term plotting for dot and heatmap
####################################################################################


go_list <- read.delim('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/Gene_sets/GO_Human.txt',sep = ' ')

go_term_celltype_compute <- function(prot_mat_in,mRNA_mat_in,deg_mat_in = NULL){
  p_val_prot <- c()
  p_val_mRNA <- c()
  terms <- c()
  med_prot <- c()
  med_rna <- c()
  med_deg <- c()
  numb <- c()
  
  prot_mat_in <- prot_mat_in[!rownames(prot_mat_in) %in% ribo,]
  mRNA_mat_in <- mRNA_mat_in[!rownames(mRNA_mat_in) %in% ribo,]
  sect <- intersect(rownames(prot_mat_in),rownames(mRNA_mat_in))
  prot_mat_in <- prot_mat_in[sect,]
  mRNA_mat_in <- mRNA_mat_in[sect,]
  
  prot_null_dist <- cor(t(prot_mat_in),use = 'pairwise.complete.obs')
  prot_null_dist <- prot_null_dist[lower.tri(prot_null_dist,diag = F)]
  
  rna_null_dist <- cor(t(mRNA_mat_in),use = 'pairwise.complete.obs')
  rna_null_dist <- rna_null_dist[lower.tri(rna_null_dist,diag = F)]
  
  for(i in unique(go_list$GO_term_name)){
    
    go_list_hold <- go_list %>% filter(GO_term_name == i)
    if(length(intersect(go_list_hold$Gene,toupper(rownames(prot_mat_in)))) > 2){
      
      prot_mat_in_hold <- prot_mat_in[toupper(rownames(prot_mat_in)) %in% go_list_hold$Gene,]
      mRNA_mat_in_hold <- mRNA_mat_in[toupper(rownames(prot_mat_in)) %in% go_list_hold$Gene,]
      #deg_mat_basal_hold <- deg_mat_in[toupper(rownames(deg_mat_basal)) %in% go_list_hold$Gene,]
      
      
      cor_prot <- cor(t(prot_mat_in_hold),use = 'pairwise.complete.obs')
      cor_prot <- cor_prot[lower.tri(cor_prot,diag = F)]
      
      cor_rna <- cor(t(mRNA_mat_in_hold),use = 'pairwise.complete.obs')
      cor_rna <- cor_rna[lower.tri(cor_rna,diag = F)]
      
      #cor_deg <- cor(t(deg_mat_basal_hold),use = 'pairwise.complete.obs')
      #cor_deg <- cor_deg[lower.tri(cor_deg,diag = F)]
      
      
      numb <- c(numb,length(intersect(go_list_hold$Gene,toupper(rownames(prot_mat_in)))) )
      if(sum(is.na(cor_prot))== length(cor_prot)){
        p_val_prot <- c(p_val_prot,NA)
      }else{
        p_val_prot <- c(p_val_prot,t.test(prot_null_dist,cor_rna)$p.value)
      }
      
      p_val_mRNA <- c(p_val_mRNA,t.test(rna_null_dist,cor_rna)$p.value)
      terms <- c(terms,i)
      med_prot <- c(med_prot,median(cor_prot,na.rm = T))
      med_rna <- c(med_rna,median(cor_rna,na.rm = T))
      #med_deg <- c(med_deg,median(cor_deg,na.rm = T))
      
      
    }
    
  }
  
  df_go <- data.frame(numb = numb,p_val_prot = p_val_prot,p_val_mRNA=p_val_mRNA,
                      go = terms, prot = med_prot, rna = med_rna) # deg = med_deg
  
  df_go <- df_go %>% filter(is.na(p_val_prot)==F)
  df_go$qval_prot <- p.adjust(df_go$p_val_prot,method = 'BH')
  df_go$qval_rna <- p.adjust(df_go$p_val_mRNA,method = 'BH')
  #df_go <- df_go %>% filter(qval_prot < .01 | qval_rna < .01 )
  return(df_go)
}


df_go_bas <- go_term_celltype_compute(prot_mat_in = prot_mat_bas,
                                    mRNA_mat_in = mRNA_mat_bas)
df_go_bas$type <- 'Bas'


df_go_bas2 <- df_go_bas %>% filter(qval_prot < .01 | qval_rna < .01 )
df_go_bas2 <- df_go_bas2 %>% filter(prot > .1 | rna > .1)
df_go_bas <- df_go_bas %>% dplyr::select('go','prot','rna','type')

df_go_fib <- go_term_celltype_compute(prot_mat_in = prot_mat_fib,
                                    mRNA_mat_in = mRNA_mat_fib)
df_go_fib$type <- 'Fib'
df_go_fib2 <- df_go_fib %>% filter(qval_prot < .01 | qval_rna < .01 )
df_go_fib2 <- df_go_fib2 %>% filter(prot > .1 | rna > .1)
df_go_fib <- df_go_fib %>% dplyr::select('go','prot','rna','type')


df_go_chond <- go_term_celltype_compute(prot_mat_in = prot_mat_chond,
                                      mRNA_mat_in = mRNA_mat_chond)
df_go_chond$type <- 'Chond'

df_go_chond2 <- df_go_chond %>% filter(qval_prot < .01 | qval_rna < .01 )
df_go_chond2 <- df_go_chond2 %>% filter(prot > .1 | rna > .1)
df_go_chond <- df_go_chond %>% dplyr::select('go','prot','rna','type')



df_go_sec <- go_term_celltype_compute(prot_mat_in = prot_mat_sec,
                                        mRNA_mat_in = mRNA_mat_sec)
df_go_sec$type <- 'Sec'
df_go_sec2 <- df_go_sec%>% filter(qval_prot < .01 | qval_rna < .01 )
df_go_sec2 <- df_go_sec2 %>% filter(prot > .1 | rna > .1)
df_go_sec <- df_go_sec %>% dplyr::select('go','prot','rna','type')


sect_go <- intersect(intersect(intersect(df_go_sec$go,df_go_chond$go),df_go_fib$go),df_go_bas$go)
sect_go <- intersect(sect_go,unique(c(df_go_sec2$go,df_go_chond2$go,df_go_fib2$go,df_go_bas2$go)))

View(df_go_chond2)
df_go_all <- rbind(df_go_fib,df_go_sec,df_go_chond,df_go_bas)
df_go_all2 <- rbind(df_go_fib2,df_go_sec2,df_go_chond2,df_go_bas2)
write.csv(df_go_all2,'/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/Supplemental_table2.csv')


df_go_all <- df_go_all %>% filter(go %in% sect_go)
df_go_all <- melt(df_go_all, ids = c('go','type'))
df_go_all$tv <- paste0(df_go_all$type,'_',df_go_all$variable)
df_go_all <- dcast(df_go_all,go~tv,value.var = 'value')
rownames(df_go_all) <- df_go_all$go
df_go_all$go <- NULL

Heatmap(as.matrix(df_go_all),cluster_rows = T,cluster_columns = F,col=col_fun)

colnames(df_go_all)
a <- Heatmap(as.matrix(df_go_all)[terms,c(1,7,3,5)],cluster_rows = T,cluster_columns = F,col=col_fun)
a 
df_rnaterm <- as.matrix(df_go_all)[terms,c(2,8,4,6)]
Heatmap(df_rnaterm[row_order(a),],cluster_rows = F,cluster_columns = F,col=col_fun)

terms <- terms[!terms %in% c('positive regulation of ubiquitin-protein ligase activity involved in mitotic cell cycle',
                             'regulation of ubiquitin-protein ligase activity involved in mitotic cell cycle',
                             'olichyl-diphosphooligosaccharide-protein glycotransferase activity',
                             'COPI coating of Golgi vesicle',
                             'SRP-dependent cotranslational protein targeting to membrane',
                             'ER to Golgi vesicle-mediated transport',
                             'nuclear-transcribed mRNA poly(A) tail shortening',
                             'COPII vesicle coating',
                             'Golgi organization',
                             'keratin filament',
                             'cell-cell adherens junction',
                             'endocytic vesicle lumen',
                             'nuclear speck',
                             'regulation of cellular amino acid metabolic process',
                             'oligosaccharyltransferase complex',
                             'proteasome core complex, alpha-subunit complex',
                             'intra-Golgi vesicle-mediated transport',
                             'endoplasmic reticulum-Golgi intermediate compartment membrane',
                             'protein N-linked glycosylation via asparagine',
                             'retrograde vesicle-mediated transport, Golgi to ER',
                             'oxidative phosphorylation',
                             'tricarboxylic acid cycle','fatty acid beta-oxidation','peroxidase activity','cytochrome-c oxidase activity')]


Heatmap(as.matrix(df_go_all)[,c(1,7,3,5)],cluster_rows = T,cluster_columns = T,col=col_fun)

# Picking a limited set of go terms to display
df_go_all2 <- as.matrix(df_go_all )
df_go_all2 <- as.data.frame(df_go_all2)
df_go_all2 <- df_go_all2 %>% filter(Bas_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Chond_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Fib_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Sec_prot > .15)

terms <- rownames(df_go_all2)

df_go_all2 <- as.matrix(df_go_all )
df_go_all2 <- as.data.frame(df_go_all2)
df_go_all2 <- df_go_all2 %>% filter(Bas_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Chond_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Fib_prot > .1)
#df_go_all2 <- df_go_all2 %>% filter(Sec_prot < .1)

terms <- c(terms,rownames(df_go_all2))
terms <- unique(terms)

df_go_all2 <- as.matrix(df_go_all )
df_go_all2 <- as.data.frame(df_go_all2)
df_go_all2 <- df_go_all2 %>% filter(Bas_prot >.2)
#df_go_all2 <- df_go_all2 %>% filter(Chond_prot > .15)
#df_go_all2 <- df_go_all2 %>% filter(Fib_prot > .3)
#df_go_all2 <- df_go_all2 %>% filter(Sec_prot > .3)

terms <- c(terms,rownames(df_go_all2))
terms <- unique(terms)

df_go_all2 <- as.matrix(df_go_all )
df_go_all2 <- as.data.frame(df_go_all2)
df_go_all2 <- df_go_all2 %>% filter(Bas_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Chond_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Fib_prot < .1)
df_go_all2 <- df_go_all2 %>% filter(Sec_prot > .2)

terms <- c(terms,rownames(df_go_all2))
terms <- unique(terms)

Heatmap(as.matrix(df_go_all)[terms,c(1,7,3,5)],cluster_rows = T,cluster_columns = T,col=col_fun)



# 
# 
# go_list_hold <- go_list %>% filter(GO_term_name == 'COPI vesicle coat') 
# intersect(go_list_hold$Gene,toupper(rownames(prot_mat_basal)))
# 
# plot_terms <- function(term){
#   
#   go_list_hold <- go_list %>% filter(GO_term_name == term)
#   
#   prot_mat_basal_hold <- prot_mat_basal[toupper(rownames(deg_mat_basal)) %in% go_list_hold$Gene,]
#   mRNA_mat_basal_hold <- mRNA_mat_basal[toupper(rownames(mRNA_mat_basal)) %in% go_list_hold$Gene,]
#   deg_mat_basal_hold <- deg_mat_basal[toupper(rownames(deg_mat_basal)) %in% go_list_hold$Gene,]
#   
#   
#   cor_prot <- cor(t(prot_mat_basal_hold),use = 'pairwise.complete.obs')
#   cor_prot <- cor_prot[lower.tri(cor_prot,diag = F)]
#   
#   cor_rna <- cor(t(mRNA_mat_basal_hold),use = 'pairwise.complete.obs')
#   cor_rna <- cor_rna[lower.tri(cor_rna,diag = F)]
#   
#   cor_deg <- cor(t(deg_mat_basal_hold),use = 'pairwise.complete.obs')
#   cor_deg <- cor_deg[lower.tri(cor_deg,diag = F)]
#   
#   
#   
#   df_plot_term <- data.frame(prot = cor_prot,rna = cor_rna,cor_deg = cor_deg)
#   df_plot_term <- reshape2::melt(df_plot_term)
#   
#   ggplot(df_plot_term,aes(x = variable,y = value))+ geom_boxplot() + dot_plot+
#     ylab('GeneXGene correlations') + xlab('')
# 
#   
#   
# }
# 
# 
# plot_term_heatmap <- function(term){
#   go_list_hold <- go_list %>% filter(GO_term_name == term)
#   go_list_hold <- go_list_hold %>% filter(Gene != 'ANXA1')
#   
#   prot_mat_basal_hold <- prot_mat_basal[toupper(rownames(deg_mat_basal)) %in% go_list_hold$Gene,]
#   mRNA_mat_basal_hold <- mRNA_mat_basal[toupper(rownames(mRNA_mat_basal)) %in% go_list_hold$Gene,]
#   deg_mat_basal_hold <- deg_mat_basal[toupper(rownames(deg_mat_basal)) %in% go_list_hold$Gene,]
#   
#   dim(prot_mat_basal_hold)
#   
#   cor_prot <- cor(t(prot_mat_basal_hold),use = 'pairwise.complete.obs')
#   cor_prot[lower.tri(cor_prot,diag = T)] <- 0
#   
#   a = Heatmap(cor_prot)
#   ht = row_order(a)
#   
#   cor_rna <- cor(t(mRNA_mat_basal_hold),use = 'pairwise.complete.obs')
#   cor_rna[upper.tri(cor_rna,diag = T)] <- 0
#   
#   
#   
#   cor_deg <- cor(t(deg_mat_basal_hold),use = 'pairwise.complete.obs')
#   cor_deg[lower.tri(cor_deg,diag = F)] <- 0
#   
#   
#   prot_rna <- cor_prot+cor_rna
#   prot_deg <- cor_rna+cor_deg
#   
#   diag(prot_rna) <- 0
#   diag(prot_deg) <- 0
#   
#   col_fun = colorRamp2(c(-.8,0,.8), c("blue","white", "red"))
#   
#   b = Heatmap(prot_rna[ht,ht],cluster_rows = F,cluster_columns = F,col = col_fun)
#   c = Heatmap(prot_deg[ht,ht],cluster_rows = F,cluster_columns = F,col = col_fun)
#   
#   b
#   
# }
# 
# 
# 
# complexes <- read.csv('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/Gene_sets/coreComplexes.csv')
# 
# 
# 
# p_val <- c()
# terms <- c()
# med_prot <- c()
# med_rna <- c()
# med_deg <- c()
# n_prot<-c()
# 
# for(i in unique(complexes$ComplexName)){
#   
#   comp_list_hold <- complexes %>% filter(ComplexName == i)
#   comp_genes <- unlist(str_split(comp_list_hold$subunits.Gene.name.,';'))
#   
#   if(length(intersect(comp_genes,toupper(rownames(prot_mat_basal)))) > 2){
#     
#     prot_mat_basal_hold <- prot_mat_basal[toupper(rownames(deg_mat_basal)) %in% comp_genes,]
#     mRNA_mat_basal_hold <- mRNA_mat_basal[toupper(rownames(mRNA_mat_basal)) %in% comp_genes,]
#     deg_mat_basal_hold <- deg_mat_basal[toupper(rownames(deg_mat_basal)) %in% comp_genes,]
#     
#     
#     cor_prot <- cor(t(prot_mat_basal_hold),use = 'pairwise.complete.obs')
#     cor_prot <- cor_prot[lower.tri(cor_prot,diag = F)]
#     
#     cor_rna <- cor(t(mRNA_mat_basal_hold),use = 'pairwise.complete.obs')
#     cor_rna <- cor_rna[lower.tri(cor_rna,diag = F)]
#     
#     cor_deg <- cor(t(deg_mat_basal_hold),use = 'pairwise.complete.obs')
#     cor_deg <- cor_deg[lower.tri(cor_deg,diag = F)]
#     
#     
#     
#     p_val <- c(p_val,t.test(cor_prot,cor_rna)$p.value)
#     terms <- c(terms,i)
#     med_prot <- c(med_prot,median(cor_prot,na.rm = T))
#     med_rna <- c(med_rna,median(cor_rna,na.rm = T))
#     med_deg <- c(med_deg,median(cor_deg,na.rm = T))
#     n_prot <- c(n_prot,length(intersect(comp_genes,toupper(rownames(prot_mat_basal)))))
#     
#   }
#   
# }
# 
# df_go <- data.frame(pval = p_val,go = terms,numb = n_prot, prot = med_prot, rna = med_rna, deg = med_deg)
# df_go$qval <- p.adjust(df_go$pval,method = 'BH')
# df_go <- df_go %>% filter(qval < .05)
# 
# hist(df_go$prot-df_go$rna)
# 






####################  KNN pseudo cell algo ###################### ######################  


library(Matrix)
library(FNN)          # fast k-nearest neighbours
library(Matrix)
library(irlba)     # fast PCA
library(ClusterR)  # kmeans++ (scales to >10k cells)


k_clusters   <- 50          # number of pseudo-cell clusters

X_top   <- t(mRNA_mat_bas)        # cells × genes (dense rows)


set.seed(42)
pca_res <- irlba::prcomp_irlba(X_top, n = 30, center = TRUE, scale. = FALSE)
pc_mat  <- pca_res$x                      # cells × 30


kfit <- ClusterR::KMeans_rcpp(
  data         = pc_mat,
  clusters     = k_clusters,
  num_init     = 5,              # multiple restarts
  max_iters    = 300,
  initializer  = 'kmeans++',
  verbose      = FALSE
)

cell_cluster <- kfit$clusters            # length = #cells (1…50)


pseudo_mat <- vapply(
  1:k_clusters,
  function(cl){
    idx <- which(cell_cluster == cl)
    Matrix::rowMeans(mRNA_mat_bas[, idx, drop = FALSE])
  },
  numeric(nrow(mRNA_mat_bas))
)

colnames(pseudo_mat) <- paste0("pseudo_", sprintf("%02d", 1:k_clusters))
rownames(pseudo_mat) <- rownames(mRNA_mat_bas)


cor_mat_bas_prot_test <- cor(t(prot_mat_bas),use = 'pairwise.complete.obs')
cor_mat_bas_mRNA_test <- cor(t(pseudo_mat))
obs_mat_bas_prot_test <- pairwiseCount(t(prot_mat_bas))
idx <- which(lower.tri(cor_mat_bas_prot_test, diag = FALSE), arr.ind = TRUE)
gene1 <- rownames(cor_mat_bas_prot_test)[idx[, 1]]   # first gene in the pair (row index)
gene2 <- colnames(cor_mat_bas_prot_test)[idx[, 2]]   # second gene (column index)

df_cor_basal <- data.frame(cor_prot = cor_mat_bas_prot_test[lower.tri(cor_mat_bas_prot_test)],
                           cor_mRNA = cor_mat_bas_mRNA_test[lower.tri(cor_mat_bas_mRNA_test)],
                           obs = obs_mat_bas_prot_test[lower.tri(obs_mat_bas_prot_test)],
                           gene1 = gene1,gene2 = gene2)
df_cor_basal <-df_cor_basal %>% filter(obs > 100 )
df_cor_basal <-df_cor_basal %>% filter(!gene1 %in% ribo)
df_cor_basal <-df_cor_basal %>% filter(!gene2 %in% ribo)   
df_cor_basal <-df_cor_basal %>% filter(!gene1 %in% mito)
df_cor_basal <-df_cor_basal %>% filter(!gene2 %in% mito)  

cor(df_cor_basal$cor_prot,df_cor_basal$cor_mRNA)





#################### Figure 5 #################### 



sect <- intersect(rownames(prot_mat_bas),rownames(prot_mat_fib))
#sect <- sect[!sect %in% c(mito,ribo)]




cor_mat_bas_prot2 <- cor_mat_bas_prot[sect,sect]
cor_mat_fib_prot2 <- cor_mat_fib_prot[sect,sect]





diag(cor_mat_fib_prot2) <- NA
diag(cor_mat_bas_prot2) <- NA

cors <- c()
cor_null <- c()
for(i in sect){
  
    cors <- c(cors,cor(cor_mat_fib_prot2[i,],cor_mat_bas_prot2[i,],use = 'pairwise.complete.obs'))
    v1 <- cor_mat_fib_prot2[sample(nrow(cor_mat_fib_prot2),1),]
    v1_shuffled <- sample(v1)
    v2 <- cor_mat_bas_prot2[sample(nrow(cor_mat_bas_prot2),1),]
    cor_null <- c(cor_null,cor(v1_shuffled,v2,use = 'pairwise.complete.obs'))
}

df_cva <- data.frame(cor = cors,null_dist = cor_null,prot = sect)
hist(cor_null)


# reshape to long --------------------------------------------------------------
df_long <- df_cva %>%
  dplyr::select(cor, null_dist) %>%
  pivot_longer(cols = everything(),
               names_to  = "type",
               values_to = "corr") %>%
  mutate(type = factor(type,
                       levels = c("null_dist", "cor"),
                       labels = c("Null", "Observed")))


ggplot(df_long, aes(x = type, y = corr, fill = type)) +
  geom_violin(alpha = .2,trim = FALSE,
              bw   = 0.015,       # ↓ bandwidth -> finer detail (try 0.01–0.02)
              scale = "width",
              adjust = 1) +       # adjust multiplies bw; keep at 1
  scale_fill_manual(values = c("Null" = "#ff7f0e", "Observed" = "#1f77b4"),
                    name = NULL) +
  labs(x = NULL,
       y = 'Correlation vector corr',
       title = "") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") 


df_cva_bad <- df_cva %>% filter(cor < .1)
df_cva_good <- df_cva %>% filter(cor > .4)

df_cva$set <- NA
df_cva$set[df_cva$prot %in% df_cva_bad$prot] <- 'Set2'
df_cva$set[df_cva$prot %in% df_cva_good$prot] <- 'Set1'
df_cva$null_dist <- NA
write.csv(df_cva,'/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/Supplemental_table3.csv')


cm1 <- cor_mat_bas_prot2[df_cva_good$prot,df_cva_good$prot]
cm2 <- cor_mat_fib_prot2[df_cva_good$prot,df_cva_good$prot]
a <- Heatmap(cor_mat_bas_prot2[df_cva_good$prot,df_cva_good$prot])
cm2 <- cm2[row_order(a),row_order(a)]
cm1 <- cm1[row_order(a),row_order(a)]
Heatmap(cm1,cluster_rows = F,cluster_columns = F)

cm1[upper.tri(cm1,diag = T)] <- 0
cm2[lower.tri(cm2,diag = T)] <- 0

Heatmap(cm1+cm2,cluster_rows = F,cluster_columns = F)


test <- functional_heatmap_annotation(cor_mat_bas_prot2,cor_mat_fib_prot2,df_cva_good$prot)
add_heatmap_clusters(test)
View(test)


G1 <- df_cva_good$prot        # proteins that form the fingerprint (rows)
G2 <- df_cva_bad$prot    # proteins you compare across cell types (cols)

cor_mat_bas_prot2[is.na(cor_mat_bas_prot2)] <- 0
cor_mat_fib_prot2[is.na(cor_mat_fib_prot2)] <- 0

cor_mat_bas_prot3 <- cor_mat_bas_prot2[G1,G2]
cor_mat_fib_prot3 <- cor_mat_fib_prot2[G1,G2]


a <- Heatmap(cor_mat_fib_prot3,show_row_names = F,show_column_names = F,show_row_dend = F,show_column_dend = F)
a2 <- cor_mat_fib_prot3[row_order(a),column_order(a)]
b <- Heatmap(cor_mat_bas_prot3[row_order(a),column_order(a)],cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F)
b2 <- cor_mat_bas_prot3[row_order(a),column_order(a)]


top_half_bas <- c()
top_half_fib <- c()
bot_half_fib <- c()
bot_half_fib <- c()
for(i in 1:ncol(a2)){
  
  top_half_bas <- c(top_half_bas,median(a2[1:158,i],na.rm = T)-median(a2[159:nrow(a2),i],na.rm = T))
  top_half_fib <- c(top_half_fib,median(b2[1:158,i],na.rm = T)-median(b2[159:nrow(a2),i],na.rm = T))
  #bot_half_fib <- c(bot_half_fib,)
  #bot_half_fib <- c(bot_half_fib,)
  
}

df_look <- data.frame(bas = top_half_bas,fib = top_half_fib,gene = colnames(a2))
df_look <- df_look %>%
  filter(abs(fib) > .05, abs(bas) > .05) %>%
  filter(sign(fib) != sign(bas))


a2 <- a2[,colnames(a2) %in% df_look$gene]
b2 <- b2[,colnames(b2) %in% df_look$gene]

a2_rna <- cor_mat_bas_mRNA[rownames(a2),colnames(a2)]
b2_rna <- cor_mat_fib_mRNA[rownames(a2),colnames(a2)]

cor_mat_bas_deg <- cor(t(deg_mat_bas),use = 'pairwise.complete.obs')
cor_mat_fib_deg <- cor(t(deg_mat_fib),use = 'pairwise.complete.obs')

a2_row <- rownames(a2)[rownames(a2) %in% rownames(cor_mat_bas_deg)]
a2_col <- colnames(a2)[colnames(a2) %in% rownames(cor_mat_bas_deg)]
a2_row <- a2_row[a2_row %in% rownames(cor_mat_fib_deg)]
a2_col <- a2_col[a2_col %in% rownames(cor_mat_fib_deg)]

a2_deg <- cor_mat_bas_deg[a2_row,a2_col]
b2_deg <- cor_mat_fib_deg[a2_row,a2_col]

col_fun = colorRamp2(c(-.4,0,.4), c("blue","white", "red"))
Heatmap(a2,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,col = col_fun)

Heatmap(b2,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,col = col_fun)

Heatmap(a2_rna,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,col = col_fun)

Heatmap(b2_rna,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,col = col_fun)

Heatmap(a2_deg,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,col = col_fun)

Heatmap(b2_deg,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,col = col_fun)



look <- cor(cor_mat_bas_prot2,cor_mat_bas_prot2,use = 'pairwise.complete.obs')

a <- Heatmap(cor(cor_mat_fib_prot3,cor_mat_fib_prot3,use = 'pairwise.complete.obs'),show_row_names = F,show_column_names = F,show_row_dend = F,show_column_dend = F)
a
Heatmap(cor(cor_mat_bas_prot3,cor_mat_bas_prot3,use = 'pairwise.complete.obs')[row_order(a),column_order(a)],cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F)



enrich_terms_debug <- function(genes_up_symbols,
                               background_symbols = NULL,
                               species = c("mouse","human"),
                               min_set_size = 5,
                               p_adj_cutoff = 1,   # start lenient; filter later
                               verbose = TRUE) {
  species <- match.arg(species)
  if (species == "mouse") {
    OrgDb <- org.Mm.eg.db; kegg_org <- "mmu"; react_org <- "mouse"
  } else {
    OrgDb <- org.Hs.eg.db; kegg_org <- "hsa"; react_org <- "human"
  }
  
  map_syms <- function(x) {
    suppressMessages(
      AnnotationDbi::select(OrgDb, keys = unique(x),
                            keytype = "SYMBOL", columns = "ENTREZID") |>
        distinct(SYMBOL, ENTREZID) |>
        filter(!is.na(ENTREZID))
    )
  }
  
  up_map <- map_syms(genes_up_symbols)
  up_entrez <- unique(up_map$ENTREZID)
  
  if (is.null(background_symbols)) {
    bg_syms <- keys(OrgDb, keytype = "SYMBOL")
  } else {
    bg_syms <- unique(background_symbols)
  }
  bg_map <- map_syms(bg_syms)
  bg_entrez <- unique(bg_map$ENTREZID)
  
  if (verbose) {
    message(sprintf("Mapped %d/%d up-regulated symbols to Entrez.",
                    length(up_entrez), length(unique(genes_up_symbols))))
    message(sprintf("Background mapped: %d/%d symbols.",
                    length(bg_entrez), length(bg_syms)))
  }
  if (length(up_entrez) < min_set_size && verbose) {
    message("Warning: up list smaller than min_set_size; consider lowering min_set_size.")
  }
  
  tidy_enrich <- function(x) {
    if (inherits(x, "enrichResult") || inherits(x, "gseaResult")) {
      df <- as.data.frame(x) |> tibble::as_tibble()
      # Some versions don’t populate 'qvalue'; avoid filtering on missing column.
      if (!"qvalue" %in% names(df)) df$qvalue <- NA_real_
      df
    } else tibble::tibble()
  }
  
  maybe_simplify <- function(er) {
    df <- tidy_enrich(er)
    if (nrow(df) == 0) return(df)
    # only simplify if there are rows
    er2 <- tryCatch(simplify(er, OrgDb = OrgDb, by = "p.adjust",
                             cutoff = 0.7, select_fun = min),
                    error = function(e) er)
    tidy_enrich(er2)
  }
  
  do_go <- function(ont) {
    er <- tryCatch(
      enrichGO(gene = up_entrez, universe = bg_entrez,
               OrgDb = OrgDb, keyType = "ENTREZID", ont = ont,
               pAdjustMethod = "BH",
               pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = min_set_size, readable = TRUE),
      error = function(e) NULL
    )
    if (is.null(er)) return(tibble::tibble())
    maybe_simplify(er) |>
      filter(Count >= min_set_size, p.adjust <= p_adj_cutoff)
  }
  
  go_bp <- do_go("BP"); go_mf <- do_go("MF"); go_cc <- do_go("CC")
  
  kegg <- tryCatch({
    enrichKEGG(gene = up_entrez, organism = kegg_org,
               pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1) |>
      setReadable(OrgDb, keyType = "ENTREZID") |>
      tidy_enrich() |>
      filter(Count >= min_set_size, p.adjust <= p_adj_cutoff)
  }, error = function(e) tibble::tibble())
  
  react <- tryCatch({
    enrichPathway(gene = up_entrez, organism = react_org,
                  pvalueCutoff = 1, pAdjustMethod = "BH",
                  qvalueCutoff = 1, minGSSize = min_set_size,
                  readable = TRUE) |>
      tidy_enrich() |>
      filter(Count >= min_set_size, p.adjust <= p_adj_cutoff)
  }, error = function(e) tibble::tibble())
  
  if (verbose) {
    message(sprintf("GO_BP rows: %d | GO_MF: %d | GO_CC: %d | KEGG: %d | Reactome: %d",
                    nrow(go_bp), nrow(go_mf), nrow(go_cc), nrow(kegg), nrow(react)))
  }
  
  list(GO_BP = go_bp, GO_MF = go_mf, GO_CC = go_cc,
       KEGG = kegg, Reactome = react,
       mapped_up = up_map, mapped_bg = bg_map)
}



test <- enrich_terms_from_set(colnames(a2)[25:ncol(a2)],background_symbols = rownames(cor_mat_bas_prot),
                              species = 'mouse')


# ----- your two sets -----
set_fib <- c("Glud1","Sfxn3","Dpysl3","Ptgis","Ugdh","Ugp2","Uap1l1","Rras","Ephx1","Parva",
             "Tram1","Hsd17b4","Rsu1","Epb41l2","Cnpy2","Ehd1","Uso1","Lman1","Erp44","Dad1",
             "Rcn3","Prkcsh","Sparc","Lypla1","Flot2","Sntb2","Lmnb2","Dpm1","Actn1","Galk1","Ftl1")

set_bas <- c("Dsp","Sec14l2","Ybx1","Rpl14","Trim28","Hnrnpa1","Agr2","Snu13","Fkbp4","Hmgb2",
             "Rnpep","Ces2c","Gsto1","Aldh1a1","Acot2","Fasn","Crocc","Tpm1","Krt15","Krt5",
             "Cndp2","Cotl1","Npepps","Dbi")

convert_set1 <- convert %>% filter(split_gene %in% set_fib)
convert_set1 <- convert_set1 %>% filter(split_gene %in% rownames(deg_mat_bas))
convert_set1 <- convert_set1 %>% filter(split_gene %in% rownames(deg_mat_fib))

convert_set2 <- convert %>% filter(split_gene %in% set_bas)
convert_set2 <- convert_set2 %>% filter(split_gene %in% rownames(deg_mat_bas))
convert_set2 <- convert_set2 %>% filter(split_gene %in% rownames(deg_mat_fib))

deg_mat_bas_norm <- Normalize_reference_vector(deg_mat_bas,log = T)
deg_mat_fib_norm <- Normalize_reference_vector(deg_mat_fib,log = T)

cor_bas_set1 <- c()
cor_fib_set1 <- c()
for(i in 1:nrow(convert_set1)){
  prot <- convert_set1$split_gene[i]
  cor_bas_set1 <- c(cor_bas_set1,cor(deg_mat_bas[prot,],prot_mat_bas[prot,],use = 'pairwise.complete.obs'))
  cor_fib_set1 <- c(cor_fib_set1,cor(deg_mat_fib[prot,],prot_mat_fib[prot,],use = 'pairwise.complete.obs'))
}
cor_bas_set2 <- c()
cor_fib_set2 <- c()
for(i in 1:nrow(convert_set2)){
  prot <- convert_set2$split_gene[i]
  cor_bas_set2 <- c(cor_bas_set2,cor(deg_mat_bas[prot,],prot_mat_bas[prot,],use = 'pairwise.complete.obs'))
  cor_fib_set2 <- c(cor_fib_set2,cor(deg_mat_fib[prot,],prot_mat_fib[prot,],use = 'pairwise.complete.obs'))
}
df_fib_make <- data.frame(cors = c(cor_fib_set1,cor_fib_set2),
                          set = c(rep('set1',length(cor_fib_set1)),rep('set2',length(cor_fib_set2))))

ggplot(df_fib_make,aes(x = set,y = cors))+ geom_boxplot() + coord_cartesian(ylim = c(-.25,.25))

df_bas_make <- data.frame(cors = c(cor_bas_set1,cor_bas_set2),
                          set = c(rep('set1',length(cor_fib_set1)),rep('set2',length(cor_fib_set2))))

ggplot(df_bas_make,aes(x = set,y = cors))+ geom_boxplot()+ coord_cartesian(ylim = c(-.25,.25))


mat_make_bas <- matrix(data = NA,ncol = 2,nrow = 2)
colnames(mat_make_bas) <- c('Purple','Green')
rownames(mat_make_bas) <- c('Top','Bottom')

mRNA_mat_bas

mat_make_bas[1,1] <- cor(colMeans(prot_mat_bas[convert_set1$split_gene,],na.rm = T),
     colMeans(prot_mat_bas[intersect(rownames(a2)[1:158],rownames(deg_mat_bas)),],na.rm = T),use = 'pairwise.complete.obs')

mat_make_bas[2,1] <- cor(colMeans(prot_mat_bas[convert_set1$split_gene,],na.rm = T),
     colMeans(prot_mat_bas[intersect(rownames(a2)[159:nrow(a2)],rownames(deg_mat_bas)),],na.rm = T),use = 'pairwise.complete.obs')

mat_make_bas[1,2] <- cor(colMeans(prot_mat_bas[convert_set2$split_gene,],na.rm = T),
     colMeans(prot_mat_bas[intersect(rownames(a2)[1:158],rownames(deg_mat_bas)),],na.rm = T),use = 'pairwise.complete.obs')

cor.test(colMeans(mRNA_mat_bas[convert_set2$split_gene,],na.rm = T),
         colMeans(mRNA_mat_bas[intersect(rownames(a2)[1:158],rownames(deg_mat_bas)),],na.rm = T))

mat_make_bas[2,2] <- cor(colMeans(prot_mat_bas[convert_set2$split_gene,],na.rm = T),
     colMeans(prot_mat_bas[intersect(rownames(a2)[159:nrow(a2)],rownames(deg_mat_bas)),],na.rm = T),use = 'pairwise.complete.obs')

cor.test(colMeans(mRNA_mat_bas[convert_set2$split_gene,],na.rm = T),
         colMeans(mRNA_mat_bas[intersect(rownames(a2)[159:nrow(a2)],rownames(deg_mat_bas)),],na.rm = T))



heat_map_clust_avg(prot_mat_fib,'Fibroblast')
heat_map_clust_avg(deg_mat_fib,'Fibroblast')
heat_map_clust_avg(mRNA_mat_fib,'Fibroblast')

heat_map_clust_avg(prot_mat_bas,'Basal')
heat_map_clust_avg(mRNA_mat_bas,'Basal')
heat_map_clust_avg(deg_mat_bas,'Basal')






# ----- tiny helpers -----
to_entrez <- function(sym) {
  bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) |>
    distinct(ENTREZID) |>
    pull(ENTREZID)
}

go_bp <- function(sym) {
  enrichGO(gene = to_entrez(sym),
           OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
           ont = "BP", pAdjustMethod = "BH",
           pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE) |>
    as.data.frame() |>
    select(ID, Description, p.adjust, GeneRatio, Count) |>
    arrange(p.adjust)
}

# ----- run -----
res_fib_BP <- go_bp(set_fib)
res_bas_BP <- go_bp(set_bas)

res_fib_BP_set1 <- go_bp(rownames(a2)[1:158])
res_bas_BP_set1 <- go_bp(rownames(a2)[159:nrow(a2)])

# view top terms
head(res_fib_BP, 20)
head(res_bas_BP, 20)
head(res_fib_BP_set1, 20)
head(res_bas_BP_set1, 20)




##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Fibroblast epithelial vim - krt5 comparison plots Fig5 pan a, b
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 


prot_vec <- df_cva_good$prot

## build one data-frame with the four scatter pairs ---------------------
df_plot <- bind_rows(
  tibble(
    panel = "Basal Krt5  vs  Fib Vim",
    x     = cor_mat_bas_prot2[prot_vec,"Krt5"],
    y     = cor_mat_fib_prot2[prot_vec,"Vim"]
  ),
  tibble(
    panel = "Basal Krt5  vs  Fib Krt5",
    x     = cor_mat_bas_prot2[prot_vec,"Krt5"],
    y     = cor_mat_fib_prot2[prot_vec,"Krt5"]
  ),
  tibble(
    panel = "Basal Vim   vs  Fib Vim",
    x     = cor_mat_bas_prot2[prot_vec,"Vim"],
    y     = cor_mat_fib_prot2[prot_vec,"Vim"]
  ),
  tibble(
    panel = "Basal Vim   vs  Fib Krt5",
    x     = cor_mat_bas_prot2[prot_vec,"Vim"],
    y     = cor_mat_fib_prot2[prot_vec,"Krt5"]
  )
)

## compute correlation for subtitle (optional) --------------------------
corr_tbl <- df_plot %>%
  group_by(panel) %>%
  summarise(r = cor(x, y, use = "pairwise.complete.obs")) %>%
  mutate(label = sprintf("r = %.2f", r))

## 2×2 facet plot --------------------------------------------------------
ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(-.4,.6) + ylim(-.4,.6)+
  facet_wrap(~ panel, ncol = 2) +
  labs(x = "Basal correlation",
       y = "Fibroblast correlation",
       title = "") +
  geom_text(data = corr_tbl, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 4, inherit.aes = FALSE) +
  theme_classic(base_size = 14)

ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(-0.4, 0.6) + ylim(-0.4, 0.6) +
  facet_wrap(~ panel, ncol = 2) +
  labs(x = "Basal correlation",
       y = "Fibroblast correlation",
       title = "") +
  geom_text(data = corr_tbl, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 4, inherit.aes = FALSE) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_blank(),                       # remove facet headers
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # border
    panel.spacing = unit(1, "lines"),                  # space between facets
    axis.line = element_blank()                        # remove original x/y axis lines
  )

df_plot_cor <- data.frame(Krt5 = cor_mat_bas_prot2['Krt5',],Vim = cor_mat_fib_prot2['Vim',])
df_plot_cor <- df_plot_cor %>% filter(abs(Krt5) > .2)
df_plot_cor <- df_plot_cor %>% filter(abs(Vim) > .25)

hm <- prot_mat_bas[c(rownames(df_plot_cor),'Krt5','Krt15'),] #'Aldh2','Cyp2f2','Cbr2','Pccb','Atp1b1'
#hm <- mRNA_mat_basal[c(rownames(df_plot_cor),'Krt5','Krt15','Aldh2','Cyp2f2'),sample(ncol(mRNA_mat_basal),2000,)]

hm <- hm[rowSums(is.na(hm)==F) > 1000,]
hm <- hm[,colSums(is.na(hm)==F) > 35]

dim(hm)

for(i in 1:nrow(hm)){
  hm[i,] <- hm[i,] - mean(hm[i,],na.rm=T)
}

a <- Heatmap(hm,show_row_dend = F,show_column_dend = F,show_column_names = F)
hm <- hm[row_order(a),]
Heatmap(hm,show_row_dend = F,show_column_dend = F,show_column_names = F,cluster_rows = F)

rownames(hm)

joint_anti <- c('Prdx6','Idh2','Aldh6a1','Cbr2','Hmgn2','Gstm2','Gsta3','Pccb') #''
joint_pro <- c('Lmna','Anxa2','Anxa1','Ehd2','Cavin1','Cav1','Eif4a1','Tuba4a','Plec','Tagln2')

unique_anti_fib <- c()
unique_pro_fib <- c('Vim','Dpysl2','Cyb5r3','Rtn4','Glud1','S100a10')

unique_anti_bas <- c('Sec14l2','Cbr2')
unique_pro_bas <- c('Krt5','Krt15','Sfn')

hm <- prot_mat_fib[c(unique_pro_fib,joint_pro,joint_anti,unique_anti_fib),]
dim(hm)
hm <- hm[,colSums(is.na(hm)==F) > 15]
dim(hm)

hm <- prot_mat_bas[c(unique_pro_bas,joint_pro,joint_anti,unique_anti_bas),]
dim(hm)
hm <- hm[,colSums(is.na(hm)==F) > 15]
dim(hm)

for(i in 1:nrow(hm)){
  hm[i,] <- hm[i,] - mean(hm[i,],na.rm=T)
}


a <- Heatmap(hm,show_row_dend = F,show_column_dend = F,show_column_names = F)
hm <- hm[row_order(a),]
Heatmap(hm,show_row_dend = F,show_column_dend = F,show_column_names = F,cluster_rows = F)

Heatmap(cor(t(hm),use = 'pairwise.complete.obs'))

comp <-c('Krt5',
         'Krt15',
         'Krt14',
         'Itga6',
         'Itgb4',
         'Plec',
         'Dsp',
         'Jup',
         'Dst',
         'Eppk1',
         'Sfn', 
         'Serpinb5',  
         'Lmna','Vim','Des','S100a4','Lgals1','Tln1','Msn') 

'Vimentin, desmin, Tgm2, Rab10 and Ywhaz'






