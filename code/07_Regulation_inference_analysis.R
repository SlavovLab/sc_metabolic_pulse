library(cmdstanr)
library(posterior)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(circlize)
library(ggbeeswarm)
library(AnnotationDbi)
library(purrr)
library(ComplexHeatmap)     # for row_order()
library(clusterProfiler)
library(org.Mm.eg.db)
library(circlize)



##################
# data path
##################

path <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'
path_code <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/code/'


#########################
# Functions for analysis
#########################


mrna_plot <- function(gene_in,df_means,df_widths,palette){
  df_means <- df_means %>% filter(gene == gene_in)
  df_widths <- df_widths %>% filter(gene == gene_in)
  
  df_plot <- df_means %>%                       # means
    left_join(df_widths,                   # widths
              by = c("gene", "cell_type")) %>%
    ## keep only the two modalities we want to show
    dplyr::select(cell_type,
                  protein_mean       = protein_abundance,
                  mRNA_mean          = mRNA,
                  protein_width,
                  mRNA_width)
  
  ## ── 2. Long format: one row per (cell type × modality) ─────────────
  df_plot_long <- df_plot %>%
    pivot_longer(
      cols      = c(protein_mean, mRNA_mean,
                    protein_width, mRNA_width),
      names_to  = c("modality", ".value"),
      names_sep = "_"
    ) %>%                         # .value makes “mean” and “width” real columns
    mutate(
      modality = recode(modality,
                        protein = "Protein",
                        mRNA    = "mRNA"),
      colour   = if_else(modality == "Protein", "orange", "#008000")
    )
  
  ## Decide what an “error bar” represents
  ## If width := SD  → use ±1 SD            (change to 2*width for ±2 SD)
  ## If width := CI width → use ±width / 2  (half‑width of the CI)
  df_plot_long <- df_plot_long %>%
    mutate(ymin = mean - width/2,
           ymax = mean + width/2)
  wanted_order <- names(palette)          # character vector of cell-type names
  
  df_plot_long <- df_plot_long %>% 
    mutate(
      cell_type = factor(cell_type, levels = wanted_order)
    )
  
  
  ## ── 3.Draw the figure ─────────────────────────────────────────────
  ggplot(df_plot_long,
         aes(x        = cell_type,
             y        = mean,
             group    = modality,
             colour   = modality)) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax),
                  width = 0.2) +
    scale_colour_manual(values = c("Protein" = "orange",
                                   "mRNA"    = "#008000")) +
    labs(x = "Cell type",
         y = "Posterior mean (log2 fold change)",
         colour = NULL,
         title = paste0(gene_in," protein vs mRNA across cell types")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

trans_plot <- function(gene_in,df_posterior,df_posteriors_width) {
  ## ── 1Filter the two master frames to this gene ──────────────────────────
  
  df_mean  <- df_posterior        %>% filter(gene == gene_in)
  df_width <- df_posteriors_width %>% filter(gene == gene_in)
  
  ## ── 2Merge means + widths and keep the three modalities we want ─────────
  df_plot <- df_mean %>%
    left_join(df_width, by = c("gene", "cell_type")) %>%
    dplyr::select(cell_type,
                  protein_mean      = protein_abundance,
                  translation_mean  = translation,
                  mRNA_mean         = mRNA,
                  protein_width,
                  translation_width,
                  mRNA_width)
  
  ## ── 3 Long format (pivot) ────────────────────────────────────────────────
  df_long <- df_plot %>%
    pivot_longer(
      cols      = c(protein_mean, translation_mean, mRNA_mean,
                    protein_width, translation_width, mRNA_width),
      names_to  = c("modality", ".value"),
      names_sep = "_"
    ) %>%
    mutate(
      modality = recode(modality,
                        protein     = "Protein",
                        translation = "Translation",
                        mRNA        = "mRNA"),
      ymin = mean - width / 2,   # use 2*width for ±2 SD if width==SD
      ymax = mean + width / 2
    )
  wanted_order <- names(palette)          # character vector of cell-type names
  
  df_long <- df_long %>% 
    mutate(
      cell_type = factor(cell_type, levels = wanted_order)
    )
  ## ── 4 Plot ───────────────────────────────────────────────────────────────
  ggplot(df_long,
         aes(x = cell_type,
             y = mean,
             group  = modality,
             colour = modality)) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
    scale_colour_manual(values = c("Protein"     = "orange",
                                   "Translation" = "blue",
                                   "mRNA"        = "#008000")) +
    labs(x = "Cell type",
         y = "Posterior mean (log2 fold change)",
         colour = NULL,
         title = paste0(gene_in, ": protein, translation and mRNA across cell types")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

deg_plot <- function(gene_in,df_posterior,df_posteriors_width) {
  ## ── 1Filter the two master frames to this gene ──────────────────────────
  
  df_mean  <- df_posterior        %>% filter(gene == gene_in)
  df_width <- df_posteriors_width %>% filter(gene == gene_in)
  
  ## ── 2Merge means + widths and keep the three modalities we want ─────────
  df_plot <- df_mean %>%
    left_join(df_width, by = c("gene", "cell_type")) %>%
    dplyr::select(cell_type,
                  protein_mean      = protein_abundance,
                  deg_mean  = deg,
                  mRNA_mean         = mRNA,
                  protein_width,
                  deg_width,
                  mRNA_width)
  
  df_plot$deg_mean <- -df_plot$deg_mean
  
  ## ── 3 Long format (pivot) ────────────────────────────────────────────────
  df_long <- df_plot %>%
    pivot_longer(
      cols      = c(protein_mean, deg_mean, mRNA_mean,
                    protein_width, deg_width, mRNA_width),
      names_to  = c("modality", ".value"),
      names_sep = "_"
    ) %>%
    mutate(
      modality = recode(modality,
                        protein     = "Protein",
                        deg = "Clearance",
                        mRNA        = "mRNA"),
      ymin = mean - width / 2,   # use 2*width for ±2 SD if width==SD
      ymax = mean + width / 2
    )
  wanted_order <- names(palette)          # character vector of cell-type names
  
  df_long <- df_long %>% 
    mutate(
      cell_type = factor(cell_type, levels = wanted_order)
    )
  ## ── 4 Plot ───────────────────────────────────────────────────────────────
  ggplot(df_long,
         aes(x = cell_type,
             y = mean,
             group  = modality,
             colour = modality)) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
    scale_colour_manual(values = c("Protein"     = "orange",
                                   "Clearance" = "#AE6828",
                                   "mRNA"        = "#008000")) +
    labs(x = "Cell type",
         y = "Posterior mean (log2 fold change)",
         colour = NULL,
         title = paste0(gene_in, ": protein, clearance, and mRNA across cell types")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


reg_heatmap <- function(i){
  
  
  df_clust <- cbind(df_posterior_mRNA[,i],df_posterior_trans[,i],-df_posterior_deg[,i],df_posterior_prot[,i]) 
  colnames(df_clust) <- c('mRNA','Trans','Deg','Protein')
  df_clust <- df_clust[rowSums(is.na(df_clust)==F)==4,]
  
  
  df_clust_mRNA <- df_clust#[df_clust[,1]>0,]
  df_clust_mRNA <- df_clust_mRNA[abs(df_clust_mRNA[,1])>.5,]
  df_clust_mRNA <- df_clust_mRNA[abs(df_clust_mRNA[,2])<.5,]
  df_clust_mRNA <- df_clust_mRNA[abs(df_clust_mRNA[,3])<.5,]
  
  #df_clust1 <- df_clust1[df_clust1[,4]>1,]
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  a <- Heatmap(df_clust_mRNA,cluster_columns = F,col = col_fun,show_row_dend = F,
               cluster_rows = T,show_row_names = T)
  
  df_clust_mRNA <- df_clust_mRNA[row_order(a),]
  
  
  df_clust_trans <- df_clust#[df_clust[,1]>0,]
  df_clust_trans <- df_clust_trans[abs(df_clust_trans[,1])<.5,]
  df_clust_trans <- df_clust_trans[abs(df_clust_trans[,2])>.5,]
  df_clust_trans <- df_clust_trans[abs(df_clust_trans[,3])<.5,]
  
  
  b <- Heatmap(df_clust_trans,cluster_columns = F,col = col_fun,show_row_dend = F,
               cluster_rows = T,show_row_names = T)
  
  df_clust_trans <- df_clust_trans[row_order(b),]
  
  df_clust_clear <- df_clust#[df_clust[,1]>0,]
  df_clust_clear <- df_clust_clear[abs(df_clust_clear[,1])<.5,]
  df_clust_clear <- df_clust_clear[abs(df_clust_clear[,2])<.5,]
  df_clust_clear <- df_clust_clear[abs(df_clust_clear[,3])>.5,]
  
  c <- Heatmap(df_clust_clear,cluster_columns = F,col = col_fun,show_row_dend = F,
               cluster_rows = T,show_row_names = T)
  
  df_clust_clear <- df_clust_clear[row_order(c),]
  
  df_all <- rbind(df_clust_mRNA,df_clust_trans,df_clust_clear)
  
  Heatmap(df_all,cluster_columns = F,col = col_fun,show_row_dend = F,
          cluster_rows = F,show_row_names = T)
  
  
  go_ids <- c("GO:0006412", "GO:0006414")        # example: translation terms
  go_genes <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys     = go_ids,
    columns  = "SYMBOL",
    keytype  = "GOALL"
  )$SYMBOL |> unique()
  membership <- ifelse(rownames(df_all) %in% go_genes, "in_term", "other")
  
  ha <- rowAnnotation(
    GO = membership,
    col = list(GO = c(in_term = "black", other = "white")),
    show_annotation_name = FALSE,
    width = unit(5, "mm")
  )
  
  Heatmap(
    as.matrix(df_all[, 1:4]),           # original four columns only
    cluster_columns = FALSE,
    show_row_dend   = FALSE,
    cluster_rows    = FALSE,
    show_row_names  = F,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  ) + ha                                     # “+” combines annotation with heatmap
  
  
}

## Convenience: SYMBOL → Entrez, first hit only ------------------------------
sym2entrez <- function(sym_vec) {
  mapIds(
    org.Mm.eg.db,
    keys      = sym_vec,
    keytype   = "SYMBOL",
    column    = "ENTREZID",
    multiVals = "first"
  ) |>
    na.omit() |>
    unname()
}

## ---- helper that returns a list(mRNA = …, Translation = …, Clearance = …) --
get_modality_sets <- function(i) {
  mat <- cbind(df_posterior_mRNA[, i],
               df_posterior_trans[, i],
               -df_posterior_deg[, i],      # clearance = −deg
               df_posterior_prot[, i])
  
  colnames(mat) <- c("mRNA", "Trans", "Clear", "Protein")
  mat <- mat[rowSums(!is.na(mat)) == 4, ]            # complete rows only
  
  ## mRNA-driven --------------------------------------------------------------
  mRNA_set <- rownames(mat)[
    abs(mat[, "mRNA" ]) >  thr_high &
      abs(mat[, "Trans" ]) <  thr_low  &
      abs(mat[, "Clear" ]) <  thr_low
  ]
  
  ## Translation-driven -------------------------------------------------------
  trans_set <- rownames(mat)[
    abs(mat[, "mRNA" ]) <  thr_low  &
      abs(mat[, "Trans" ]) >  thr_high &
      abs(mat[, "Clear" ]) <  thr_low
  ]
  
  ## Clearance-driven ---------------------------------------------------------
  clear_set <- rownames(mat)[
    abs(mat[, "mRNA" ]) <  thr_low  &
      abs(mat[, "Trans" ]) <  thr_low &
      abs(mat[, "Clear" ]) >  thr_high
  ]
  
  list(mRNA = mRNA_set,
       Translation = trans_set,
       Clearance   = clear_set)
}



plot_GO_dot_clustered <- function(ego_obj, title,
                                  modality_matrix,
                                  top_n = 10,
                                  dist_method = "euclidean",
                                  hclust_method = "complete") {
  
  x_levels <- colnames(modality_matrix)  # << fixed X order
  
  # 1) keep top_n per cluster; size = -log10(adj p)
  res <- ego_obj@compareClusterResult %>%
    as_tibble() %>%
    group_by(Cluster) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  score_matrix <- res %>%
    mutate(score = -log10(p.adjust)) %>%
    select(Description, Cluster, score) %>%
    pivot_wider(names_from = Cluster, values_from = score, values_fill = 0) %>%
    column_to_rownames("Description") %>%
    as.matrix()
  
  # ---- FIX X ORDER (no column clustering) -----------------------------------
  # ensure all x_levels exist; add 0 columns for any missing clusters
  missing_cols <- setdiff(x_levels, colnames(score_matrix))
  if (length(missing_cols) > 0) {
    add <- matrix(0, nrow = nrow(score_matrix), ncol = length(missing_cols),
                  dimnames = list(rownames(score_matrix), missing_cols))
    score_matrix <- cbind(score_matrix, add)
  }
  # reorder columns to fixed x_levels
  score_matrix <- score_matrix[, x_levels, drop = FALSE]
  
  # 2) cluster rows ONLY
  row_ord <- if (nrow(score_matrix) > 1)
    hclust(dist(score_matrix, method = dist_method), method = hclust_method)$order else 1
  score_matrix <- score_matrix[row_ord, , drop = FALSE]
  
  # 3) long format for plotting
  res_long <- as.data.frame(score_matrix) %>%
    rownames_to_column("Description") %>%
    pivot_longer(-Description, names_to = "Cluster", values_to = "score")
  
  # 4) union-based term genes (assumes readable=TRUE so geneID are SYMBOLs separated by "/")
  term_genes <- ego_obj@compareClusterResult %>%
    as_tibble() %>%
    select(Description, geneID) %>%
    mutate(gene_list = strsplit(geneID, "/")) %>%
    select(Description, gene_list) %>%
    group_by(Description) %>%
    summarise(genes = list(unique(unlist(gene_list))), .groups = "drop")
  
  # 5) compute per-term, per-cluster median effect from modality_matrix
  all_clusters <- x_levels
  med_tbl <- expand_grid(Description = term_genes$Description,
                         Cluster = all_clusters) %>%
    left_join(term_genes, by = "Description") %>%
    rowwise() %>%
    mutate(
      median_effect = {
        g <- intersect(genes, rownames(modality_matrix))
        if (length(g) == 0) NA_real_
        else median(as.numeric(modality_matrix[g, Cluster, drop = TRUE]), na.rm = TRUE)
      }
    ) %>%
    ungroup() %>%
    select(Description, Cluster, median_effect) %>%
    tidyr::replace_na(list(median_effect = 0))
  
  # 6) join medians to plotting frame
  res_long <- res_long %>% left_join(med_tbl, by = c("Description","Cluster"))
  
  ggplot(res_long,
         aes(x = factor(Cluster, levels = x_levels),
             y = factor(Description, levels = rev(rownames(score_matrix))),
             size = score,
             colour = median_effect)) +
    geom_point() +
    scale_size_continuous(range = c(1, 6),
                          name  = expression(-log[10](adj~p))) +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank()) +
    scale_colour_gradient2(
      name = "Median effect",
      midpoint = 0, low = "blue", mid = "white", high = "red",
      limits = c(-1, 1),
      oob = scales::squish,
      breaks = c(-1, -0.5, 0, 0.5, 1)
    )
}




waterfall_arrow_plot <- function(cell_type,
                                 protein_name,
                                 delta_mrna,  se_mrna,
                                 delta_translation, se_translation,
                                 delta_degradation, se_degradation,
                                 baseline = 0) {
  
  # --------------------------------------------------------------------------
  contrib_tbl <- tibble(
    effect = factor(c("mRNA", "Translation", "Degradation"),
                    levels = c("mRNA", "Translation", "Degradation")),
    delta = c(delta_mrna, delta_translation, delta_degradation),
    se    = c(se_mrna,    se_translation,    se_degradation)
  ) |>
    mutate(x     = as.numeric(effect),
           start = baseline + c(0, head(cumsum(delta), -1)),
           end   = start + delta,
           tip   = end,                         # where the arrow ends
           ymin  = tip - se,                   # error bar (± se)
           ymax  = tip + se)
  
  total <- baseline + sum(contrib_tbl$delta)
  total_x <- length(levels(contrib_tbl$effect)) + 1
  
  ggplot(contrib_tbl) +
    geom_hline(yintercept = baseline, linetype = "dashed") +
    
    ## 1 ⃣  error bar first (so the arrow sits on top of it)
    geom_linerange(aes(x = x, ymin = ymin, ymax = ymax),
                   linewidth = .8, colour = "grey30") +
    
    ## 2 ⃣  contribution arrow on top
    geom_segment(aes(x = x, xend = x, y = start, yend = end, colour = effect),
                 linewidth = 3,
                 arrow = arrow(type   = "closed",
                               length = unit(0.35, "cm"),
                               angle  = 30)) +
    
    ## 3 ⃣  total bar
    annotate("segment", x = total_x, xend = total_x,
             y = baseline, yend = total, linewidth = 1.4) +
    annotate("text", x = total_x, y = total,
             label = sprintf("%.2f", total), vjust = -0.6, fontface = "bold") +
    
    ## 4 ⃣  labels (optional)
    geom_text(aes(x = x, y = tip, label = sprintf("%.2f", delta)),
              vjust = ifelse(contrib_tbl$delta >= 0, -0.4, 1.2), size = 4) +
    
    scale_colour_manual(values = c(mRNA = "#008000",
                                   Translation = "blue",
                                   Degradation = "purple")) +
    scale_x_continuous(breaks = c(1:3, total_x),
                       labels = c(levels(contrib_tbl$effect), "")) +
    labs(x = NULL, y = "Log2 protein fold change",
         title = protein_name, subtitle = cell_type) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}


param_to_long <- function(draws, par_regex, value_name,
                          genes, celltypes) {
  m <- as_draws_matrix(draws)              # iterations × (G·C)
  means <- colMeans(m)
  
  tibble(variable = names(means),
         value    = means) %>%
    extract(variable,
            into = c("gene_idx", "ct_idx"),
            regex = par_regex,              # e.g. "logP\\[([0-9]+),([0-9]+)\\]"
            convert = TRUE) %>%
    mutate(
      gene      = genes[gene_idx],
      cell_type = celltypes[ct_idx],
      !!value_name := value
    ) %>%
    dplyr::select(gene, cell_type, !!value_name)
}

param_to_long_width <- function(draws, par_regex, value_name,
                                genes, celltypes,
                                type = c("sd", "ci"),  # choose: "sd" or "ci"
                                prob = 0.95) {         # CI mass, if type == "ci"
  type <- match.arg(type)
  m <- posterior::as_draws_matrix(draws)      # iterations × (G·C)
  
  width <- switch(type,
                  sd = apply(m, 2, sd),
                  ci = {
                    q_lo <- (1 - prob) / 2
                    q_hi <- 1 - q_lo
                    qs   <- apply(m, 2, stats::quantile, probs = c(q_lo, q_hi))
                    qs[2, ] - qs[1, ]                          # upper − lower
                  }
  )
  
  tibble(variable = names(width),
         value    = as.numeric(width)) %>%
    tidyr::extract(
      variable,
      into   = c("gene_idx", "ct_idx"),
      regex  = par_regex,                      # e.g. "logP\\[([0-9]+),([0-9]+)\\]"
      convert = TRUE
    ) %>%
    mutate(
      gene      = genes[gene_idx],
      cell_type = celltypes[ct_idx],
      !!value_name := value
    ) %>%
    dplyr::select(gene, cell_type, !!value_name)
}


palette <- c(
  Immune       = "#F90303",
  Basal        = "#B50202",
  Secratory    = "#B606C4",
  Cilliated    = "#9B70F9",
  Fibroblast   = "#2C78FF",
  `Smooth muscle` = "#0498BA",
  Chondrocyte  = "#03C1ED"
)




#########################
# Read in data for analysis
#########################
model_file <- file.path(paste0(path_code,"07_regulation_model.stan"))
model_file <- file.path("/Users/andrewleduc/Desktop/Projects/Miceotopes/SingleCell/Model.stan")



mRNA_dat <- read.csv(paste0(path,'05_Stan_model_input_data/rna_input_stan.csv')) 
protien_dat <- read.csv(paste0(path,'05_Stan_model_input_data/protein_input_stan.csv')) 
clearance_dat <- read.csv(paste0(path,'05_Stan_model_input_data/clearance_stan.csv')) 


prot_abs <- read.csv(paste0(path,'06_Gene_X_CellType/Absolute_abundance/abundance.csv'),row.names = 1) 
deg_abs <- read.csv(paste0(path,'06_Gene_X_CellType/Absolute_abundance/clearance.csv'),row.names = 1) 


go_list <- read.delim(paste0(path,'GO_Human.txt'),sep = ' ')



#########################
# Compare quickly relative levels explained by deg for two different cell type comparisons
# Figure 3d
#########################

df_plot <- tibble(
  diff_prot = log2(prot_abs$Basal) - log2(prot_abs$Chondrocyte),
  diff_deg  = log2(deg_abs$Basal)  - log2(deg_abs$Chondrocyte)
)

ggplot(df_plot, aes(x = diff_prot, y = diff_deg)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 14) + ylim(c(-5,5)) + xlim(c(-5,5)) + geom_abline(intercept = 0, slope = -1)

mat_dif <- matrix(data = NA,nrow = ncol(prot_abs),ncol= ncol(prot_abs))
colnames(mat_dif) <- colnames(prot_abs)
rownames(mat_dif) <- colnames(prot_abs)

df_rely <- data.frame(cell_type = c('Basal','Fibroblast','Secratory','Cilliated','Chondrocyte','Immune','Smooth.muscle'),
                      rely = c(0.7064013,0.5578751,0.6354048,0.5404068,0.7265930,0.4617232,0.6183131))


for(i in 1:(ncol(prot_abs)-1)){
  for(j in (i+1):ncol(prot_abs)){
    cor_val = cor(log2(prot_abs[,i]) - log2(prot_abs[,j]),
                        log2(deg_abs[,i]) - log2(deg_abs[,j]),use = 'pairwise.complete.obs')
    
    rely_val <- df_rely$rely[df_rely$cell_type == colnames(prot_abs)[i]]*df_rely$rely[df_rely$cell_type == colnames(prot_abs)[j]]
    if(cor_val < 0){
      cor_val <- cor_val^2 / rely_val
    }else{
      cor_val <- .01
    }
    
    mat_dif[j,i] <- cor_val
  }
}


diag(mat_dif) <- 0

mat_dif[upper.tri(mat_dif)] <- t(mat_dif)[upper.tri(mat_dif)] 
col_fun <- colorRamp2(c(0, .55), c("white", "red"))
Heatmap(mat_dif,col = col_fun)
Heatmap(
  mat_dif,
  col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- mat_dif[i, j]
    if (!is.na(v)) {
      grid.text(sprintf("%.2f", v), x, y, gp = gpar(col = "black", fontsize = 8))
    }
  }
)

#########################
# Stan model 
#########################

# input data
stan_data <- list(
  G = length(unique(mRNA_dat$gene)),
  C = length(unique(mRNA_dat$Cell_type)),
  
  # Protein
  M_pro = nrow(protien_dat),
  gene_id_pro = protien_dat$gene_id_pro,
  celltype_id_pro = protien_dat$celltype_id_pro,
  bar_pro = protien_dat$regressed_p_obs,
  n_pro = protien_dat$n_obs,
  
  # mRNA
  M_mrna = nrow(mRNA_dat),
  gene_id_mrna = mRNA_dat$gene_id_mrna,
  celltype_id_mrna = mRNA_dat$celltype_id_mrna,
  bar_mrna = mRNA_dat$value,
  #n_mRNA = mRNA_dat$n_obs,
  
  # Deg
  M_deg = nrow(clearance_dat),
  gene_id_deg = clearance_dat$gene_id_deg,
  celltype_id_deg = clearance_dat$celltype_id_deg,
  bar_deg = clearance_dat$value,
  n_deg = clearance_dat$n_obs
)



#------------------------------------------------------------------------------
# Fit the main model in ~/model.stan
#------------------------------------------------------------------------------

mod_main <- cmdstan_model(model_file)

# #### Run MCMC 
# fit_main <- mod_main$sample(
#   data = stan_data,
#   seed = 121,
#   chains = 4,
#   parallel_chains = 4,
#   iter_warmup = 1000,
#   iter_sampling = 2000
# )

#### Run VI model 


fit_vi <- mod_main$variational(
  data = stan_data,
  seed = 123,
  algorithm = "meanfield"
)



################ ################ ################ ################ 
####  Get gene and cell type anno from model
################ ################ ################ ################ 

all_genes <- mRNA_dat |>
  arrange(gene_id_mrna) |>         # gene_id_mrna is 1-based index in CSV
  distinct(gene_id_mrna, gene) |>
  pull(gene)

# cell types in the order 1 … C
all_celltypes <- mRNA_dat |>
  arrange(celltype_id_mrna) |>
  distinct(celltype_id_mrna, Cell_type) |>
  pull(Cell_type)




################ ################ ################ ################ 
####  Get matricies
################ ################ ################ ################ 



df_prot <- param_to_long(fit_vi$draws("logP"),
                         "logP\\[([0-9]+),([0-9]+)\\]",
                         "protein_abundance",
                         all_genes, all_celltypes)

df_deg  <- param_to_long(fit_vi$draws("logDeg"),
                         "logDeg\\[([0-9]+),([0-9]+)\\]",
                         "deg",
                         all_genes, all_celltypes)

df_trans<- param_to_long(fit_vi$draws("T"),
                         "T\\[([0-9]+),([0-9]+)\\]",
                         "translation",
                         all_genes, all_celltypes)

df_mrna <- param_to_long(fit_vi$draws("logM"),
                         "logM\\[([0-9]+),([0-9]+)\\]",
                         "mRNA",
                         all_genes, all_celltypes)





df_prot_w <- param_to_long_width(
  fit_vi$draws("logP"),  "logP\\[([0-9]+),([0-9]+)\\]",
  "protein_width",       # name of the new column
  all_genes, all_celltypes,
  type = "ci"            # or type = "ci", prob = 0.80, …
)

df_deg_w <- param_to_long_width(
  fit_vi$draws("logDeg"),"logDeg\\[([0-9]+),([0-9]+)\\]",
  "deg_width",
  all_genes, all_celltypes,
  type = "ci"
)

df_trans_w <- param_to_long_width(
  fit_vi$draws("T"),     "T\\[([0-9]+),([0-9]+)\\]",
  "translation_width",
  all_genes, all_celltypes,
  type = "ci"
)

df_mrna_w <- param_to_long_width(
  fit_vi$draws("logM"),  "logM\\[([0-9]+),([0-9]+)\\]",
  "mRNA_width",
  all_genes, all_celltypes,
  type = "ci"
)



df_mrna$mRNA[df_prot_w$protein_width > 20] <- NA
df_trans$translation[df_prot_w$protein_width > 20] <- NA
df_deg$deg[df_prot_w$protein_width > 20] <- NA
df_prot$protein_abundance[df_prot_w$protein_width > 20] <- NA

df_posterior <- df_prot  %>%
  left_join(df_deg,   by = c("gene", "cell_type")) %>%
  left_join(df_trans, by = c("gene", "cell_type")) %>%
  left_join(df_mrna,  by = c("gene", "cell_type"))

df_posteriors_width <- df_prot_w %>%
  left_join(df_deg_w,   by = c("gene", "cell_type")) %>%
  left_join(df_trans_w, by = c("gene", "cell_type")) %>%
  left_join(df_mrna_w,  by = c("gene", "cell_type"))



####  making the individual matricies for each mode of regulation


df_posterior_mRNA <- df_posterior %>% dplyr::select(gene,cell_type,mRNA)
df_posterior_prot <- df_posterior %>% dplyr::select(gene,cell_type,protein_abundance)
df_posterior_deg <- df_posterior %>% dplyr::select(gene,cell_type,deg)
df_posterior_trans <- df_posterior %>% dplyr::select(gene,cell_type,translation)

df_posterior_mRNA <- dcast(df_posterior_mRNA,gene~cell_type,value.var = 'mRNA')
df_posterior_prot <- dcast(df_posterior_prot,gene~cell_type,value.var = 'protein_abundance')
df_posterior_deg <- dcast(df_posterior_deg,gene~cell_type,value.var = 'deg')
df_posterior_trans <- dcast(df_posterior_trans,gene~cell_type,value.var = 'translation')


rownames(df_posterior_mRNA) <- df_posterior_mRNA$gene
df_posterior_mRNA$gene <- NULL
df_posterior_mRNA <- as.matrix(df_posterior_mRNA)
rownames(df_posterior_prot) <- df_posterior_prot$gene
df_posterior_prot$gene <- NULL
df_posterior_prot <- as.matrix(df_posterior_prot)
rownames(df_posterior_deg) <- df_posterior_deg$gene
df_posterior_deg$gene <- NULL
df_posterior_deg <- as.matrix(df_posterior_deg)
rownames(df_posterior_trans) <- df_posterior_trans$gene
df_posterior_trans$gene <- NULL
df_posterior_trans <- as.matrix(df_posterior_trans)


df_posterior_mRNA_w <- df_posteriors_width %>% dplyr::select(gene,cell_type,mRNA_width)
df_posterior_prot_w <- df_posteriors_width %>% dplyr::select(gene,cell_type,protein_width)
df_posterior_deg_w <- df_posteriors_width %>% dplyr::select(gene,cell_type,deg_width)
df_posterior_trans_w <- df_posteriors_width %>% dplyr::select(gene,cell_type,translation_width)

df_posterior_mRNA_w <- dcast(df_posterior_mRNA_w,gene~cell_type,value.var = 'mRNA_width')
df_posterior_prot_w <- dcast(df_posterior_prot_w,gene~cell_type,value.var = 'protein_width')
df_posterior_deg_w <- dcast(df_posterior_deg_w,gene~cell_type,value.var = 'deg_width')
df_posterior_trans_w <- dcast(df_posterior_trans_w,gene~cell_type,value.var = 'translation_width')


rownames(df_posterior_mRNA_w) <- df_posterior_mRNA_w$gene
df_posterior_mRNA_w$gene <- NULL
df_posterior_mRNA_w <- as.matrix(df_posterior_mRNA_w)
rownames(df_posterior_prot_w) <- df_posterior_prot_w$gene
df_posterior_prot_w$gene <- NULL
df_posterior_prot_w <- as.matrix(df_posterior_prot_w)
rownames(df_posterior_deg_w) <- df_posterior_deg_w$gene
df_posterior_deg_w$gene <- NULL
df_posterior_deg_w <- as.matrix(df_posterior_deg_w)
rownames(df_posterior_trans_w) <- df_posterior_trans_w$gene
df_posterior_trans_w$gene <- NULL
df_posterior_trans_w <- as.matrix(df_posterior_trans_w)


# write.csv(df_posterior_trans,paste0(path,'06_Gene_X_CellType/Relative_abundance/translation_bayes.csv'))
# write.csv(df_posterior_trans_w,paste0(path,'06_Gene_X_CellType/Relative_abundance/translation_bayes_95CI.csv'))
# 
# write.csv(df_posterior_deg,paste0(path,'06_Gene_X_CellType/Relative_abundance/clearance_bayes.csv'))
# write.csv(df_posterior_deg_w,paste0(path,'06_Gene_X_CellType/Relative_abundance/clearance_bayes_95CI.csv'))
# 
# write.csv(df_posterior_prot,paste0(path,'06_Gene_X_CellType/Relative_abundance/protein_bayes.csv'))
# write.csv(df_posterior_prot_w,paste0(path,'06_Gene_X_CellType/Relative_abundance/protein_bayes_95CI.csv'))
# 
# write.csv(df_posterior_mRNA,paste0(path,'06_Gene_X_CellType/Relative_abundance/mRNA_bayes.csv'))
# write.csv(df_posterior_mRNA_w,paste0(path,'06_Gene_X_CellType/Relative_abundance/mRNA_bayes_95CI.csv'))

################ ################ ################ ################ 
####  Correlation distributions to look at regulation for each mode
################ ################ ################ ################ 


cors_rna <- c()
cors_trans <- c()
cors_deg <- c()
abs_prot <- c()
numb <- c()

deg_null <- c()
rna_null <- c()
trans_null<- c()

for(i in 1:nrow(df_posterior_mRNA)){
  cors_rna <- c(cors_rna,cor(df_posterior_mRNA[i,],df_posterior_prot[i,],use = 'pairwise.complete.obs'))
  rna_null <- c(rna_null,cor(df_posterior_mRNA[i,],df_posterior_prot[sample(nrow(df_posterior_prot),1),],use = 'pairwise.complete.obs'))
  cors_trans <- c(cors_trans,cor(df_posterior_trans[i,],df_posterior_prot[i,],use = 'pairwise.complete.obs'))
  trans_null <- c(trans_null,cor(df_posterior_trans[i,],df_posterior_prot[sample(nrow(df_posterior_prot),1),],use = 'pairwise.complete.obs'))
  cors_deg <- c(cors_deg,cor(df_posterior_deg[i,],df_posterior_prot[i,],use = 'pairwise.complete.obs'))
  deg_null <- c(deg_null,cor(df_posterior_deg[i,],df_posterior_prot[sample(nrow(df_posterior_prot),1),],use = 'pairwise.complete.obs'))
  
  abs_prot <- c(abs_prot,sum(abs(df_posterior_prot[i,]),na.rm=T))
  numb <- c(numb,sum(is.na(df_posterior_prot[i,])==F))
  
  
}


df_cond <- data.frame(mRNA = cors_rna,trans = cors_trans,deg = -cors_deg,var = abs_prot,numb = numb)
df_cond$gene <- rownames(df_posterior_mRNA)

df_cond <- df_cond %>% filter(numb > 4)

##############################
### Plot individual proteins and their regulation across cell types
##############################
mrna_plot('Sult1d1',df_posterior,df_posteriors_width,palette)

trans_plot('Sec14l2',df_posterior,df_posteriors_width)

deg_plot('Tuba4a',df_posterior,df_posteriors_width)

trans_plot('St13',df_posterior,df_posteriors_width)




df_binned <- df_cond %>%
  mutate(
    deg_bin_num = round(deg/2,1)*2 ,                     # nearest 0.2
  )


df_binned$deg_bin  = factor(df_binned$deg_bin_num, levels = sort(unique(df_binned$deg_bin_num)))

ggplot(df_binned, aes(x = deg_bin, y = mRNA)) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.2) +
  labs(
    x = "Clearance half life correlation",
    y = "mRNA correlation"
  ) +
  theme_classic(base_size = 14)

###

df_binned <- df_cond %>%
  mutate(
    deg_bin_num = round(trans/2,1)*2 ,                     # nearest 0.2
  )


df_binned$deg_bin  = factor(df_binned$deg_bin_num, levels = sort(unique(df_binned$deg_bin_num)))

ggplot(df_binned, aes(x = deg_bin, y = mRNA)) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.2) +
  labs(
    x = "Translation correlation",
    y = "mRNA correlation"
  ) +
  theme_classic(base_size = 14)

###




#### Histogram for all the modes correlating to protein abundance for each gene

df_cond$var <- NULL
df_cond$numb <- NULL
df_cond.m <- melt(df_cond)

med_df <- df_cond.m %>%
  group_by(variable) %>%
  summarise(med = median(value, na.rm = TRUE), .groups = "drop")

## 2) histogram + median line -----------------------------------------
ggplot(df_cond.m, aes(y = value)) +                             # x = value
  geom_histogram(binwidth = 0.05, colour = "black", fill = "grey80") +
  geom_hline(data = med_df, aes(yintercept = med),            # red median
             colour = "red", linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(y = "Correlation to protein abundance",
       x = "# Genes") +
  theme_bw(base_size = 15) + ylim(c(-1,1))



################ ################ ################ ################ 
####  Correlation go term analysis
################ ################ ################ ################ 

pvals_rna <- c()
pvals_trans <- c()
pvals_deg <- c()
med_cor_rna <- c()
med_cor_trans <- c()
med_cor_deg <- c()
term <- c()
for(i in unique(go_list$GO_term_name)){
  go_hold <- go_list %>% filter(GO_term_name == i)
  if(length(intersect(df_cond$gene,go_hold$Gene))>3){
    
    df_hold <- df_cond %>% filter(gene %in% go_hold$Gene)
    
    pvals_rna <- c(pvals_rna,t.test(df_hold$mRNA,rna_null)$p.value)
    pvals_trans <- c(pvals_trans,t.test(df_hold$trans,trans_null)$p.value)
    pvals_deg <- c(pvals_deg,t.test(df_hold$deg,deg_null)$p.value)
    
    med_cor_rna <- c(med_cor_rna,median(df_hold$mRNA))
    med_cor_trans <- c(med_cor_trans,median(df_hold$trans))
    med_cor_deg <- c(med_cor_deg,median(df_hold$deg))
    
    term <- c(term,i)
    
  }
}

df_use_mRNA <- data.frame(p_val = pvals_rna,cor =med_cor_rna,go =  term)
df_use_trans <- data.frame(p_val = pvals_trans,cor =med_cor_trans,go =  term)
df_use_deg <- data.frame(p_val = pvals_deg,cor =med_cor_deg,go =  term)

df_use_mRNA$qval <- p.adjust(df_use_mRNA$p_val,method = 'BH')
df_use_trans$qval <- p.adjust(df_use_trans$p_val,method = 'BH')
df_use_deg$qval <- p.adjust(df_use_deg$p_val,method = 'BH')

df_use_mRNA <- df_use_mRNA %>% filter(qval < .05 & cor > 0)
df_use_trans <- df_use_trans %>% filter(qval < .05 & cor > 0)
df_use_deg <- df_use_deg %>% filter(qval < .05 & cor < 0)
df_use_deg$cor <- -df_use_deg$cor

df_use_mRNA$type <- 'mRNA'
df_use_trans$type <- 'Translation'
df_use_deg$type <- 'Clearance'



df_go_plot <- df_cond
# Deg terms
deg_terms <- c('respiratory electron transport chain',
               'cytosolic small ribosomal subunit',
               'mitotic nuclear envelope reassembly',
               'glucose metabolic process',
               'proteasome accessory complex')

deg_terms <- c('glycogen catabolic process',
               'activation of MAPKK activity',
               'mRNA transport')




count = 0
for(i in deg_terms){
  if(count == 0){
    go_hold <- go_list %>% filter(GO_term_name == i)
    df_go_plot <- df_cond %>% filter(gene %in% go_hold$Gene)
    df_go_plot$term = i
    df_go_plot$type = 'Deg'
    df_go_plot$mRNA <- NULL
    df_go_plot$trans <- NULL
    colnames(df_go_plot)[colnames(df_go_plot) == 'deg'] <- 'Value'
    count = 1
  }else{
    go_hold <- go_list %>% filter(GO_term_name == i)
    df_go_plot2 <- df_cond %>% filter(gene %in% go_hold$Gene)
    df_go_plot2$term = i
    df_go_plot2$type = 'Deg'
    df_go_plot2$mRNA <- NULL
    df_go_plot2$trans <- NULL
    colnames(df_go_plot2)[colnames(df_go_plot2) == 'deg'] <- 'Value'
    
    df_go_plot <- rbind(df_go_plot,df_go_plot2)
  }
  
}

df_go_plot$Value = - df_go_plot$Value



# mRNA terms
mRNA_terms <- c('glutathione metabolic process',
                'Rho protein signal transduction',
                'actin filament capping',
                'stress fiber')


for(i in mRNA_terms){
  go_hold <- go_list %>% filter(GO_term_name == i)
  df_go_plot2 <- df_cond %>% filter(gene %in% go_hold$Gene)
  df_go_plot2$term = i
  df_go_plot2$type = 'mRNA'
  df_go_plot2$deg <- NULL
  df_go_plot2$trans <- NULL
  colnames(df_go_plot2)[colnames(df_go_plot2) == 'mRNA'] <- 'Value'
  
  df_go_plot <- rbind(df_go_plot,df_go_plot2)
  
  
}

# Translation terms
trans_terms <- c('regulation of protein localization',
                 'kinase binding',
                 'intrinsic apoptotic signaling pathway',
                 'nuclear speck',
                 'toxin metabolic process')


for(i in trans_terms){
  go_hold <- go_list %>% filter(GO_term_name == i)
  df_go_plot2 <- df_cond %>% filter(gene %in% go_hold$Gene)
  df_go_plot2$term = i
  df_go_plot2$type = 'Trans'
  df_go_plot2$deg <- NULL
  df_go_plot2$mRNA <- NULL
  colnames(df_go_plot2)[colnames(df_go_plot2) == 'trans'] <- 'Value'
  
  df_go_plot <- rbind(df_go_plot,df_go_plot2)
  
}


ggplot(df_go_plot,aes(y = Value,x = term)) + ggbeeswarm::geom_beeswarm(cex = 2,method = 'swarm') +
  facet_wrap(~type,ncol = 3,scales = 'free_x') + xlab('Correlation')+
  ylab('') + theme_bw(base_size = 15) + ylim(c(-1,1)) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))


########## ########## ########## ########## ########## ##########
# Heatmap by important protein groups, Fig 3 C
########## ########## ########## ########## ########## ##########
e3_ligase <- c("Nedd4","Trim23","Trim28","Trim42","Prpf19","Ufl1" ,"Fbxo2","Ddb1","Kctd12")

tf_genes <- c(
  "Stat3", "Nfix", "Nfib", "Tsc22d1", "Gtf2i", "Myef2", "Tfam"
)

chromatin_regulators <- c(
  "Chd4", "Smarcc2", "Trim28", "Psip1",
  "Cbx1", "Cbx5",
  "Hmga1", "Hmgb1", "Hmgb2", "Hmgb3", "Hmgn1", "Hmgn2", "Hmgn5",
  "Mecp2", "Mta2",
  "Nap1l1", "Anp32a", "Anp32e",
  "Ruvbl1", "Ruvbl2",
  "Rad21", "Set", "Hdgf", "Zmym4"
)

chromatin_regulators <- c(tf_genes,chromatin_regulators)


proteosome_core_complex <- c("Psma1"   ,  "Psma2"   ,  "Psma3"   ,  "Psma4"  ,   "Psma5"  ,   "Psma6"  ,   "Psma7"  ,   "Psmb3"  ,   "Psmb4"  ,   "Psmb6" ,   
                             "Psmb7"    , "Psmc1"  ,   "Psmc2"   ,  "Psmc3"  ,   "Psmc4"  ,   "Psmc5"  ,   "Psmc6"  ,   "Psmd1"  ,   "Psmd11"  ,  "Psmd12" ,  
                             "Psmd13" ,   "Psmd14"  ,  "Psmd2"  ,   "Psmd3"   ,  "Psmd4"  ,   "Psmd7"   ,  "Psmd8"  ,   "Psme1"  ,   "Psme2" )




glycolysis <- c(
  "Hk1",
  "Pfkl","Pfkm","Pfkp",
  "Aldoa","Aldoc",
  "Tpi1",
  "Gapdh",
  "Pgk1",
  "Pgam1",
  "Eno1",
  "Pkm",
  "Ldha","Ldhb"
)
membrane_pm <- c(
  # Channels & pumps
  "Aqp1","Aqp5",
  "Atp1a1","Atp1a2","Atp1b1","Atp1b3","Atp2b1",
  # Transporters (PM)
  "Slc1a5","Slc3a2","Slc12a2",
  # Adhesion / receptors
  "Cdh1","Dag1","Epcam","Icam1",
  "Itga6","Itgav","Itgb1","Itgb4","Itgb5","Itgb6","Itgb7",
  # CD / tetraspanins / junction
  "Cd9","Cd34","Cd36","Cd44","Cd47","Cd200","Cldn3",
  # Rafts / caveolae
  "Cav1","Flot1","Flot2",
  # GPI-anchored / cell-surface glycoproteins
  "Lypd2","Meltf","Msln","Mpz","Emb",
  # ECM/cell-surface enzyme
  "Entpd1","Plpp3",
  # Collagen XVII (hemidesmosomes)
  "Col17a1",
  # Tetraspanins
  "Tspan6","Tspan8",
  # Ferlin family (C-term TM; PM/endosome)
  "Myof"
)


ribo <- c("Rpl10a","Rpl11","Rpl12","Rpl13","Rpl13a","Rpl14","Rpl15","Rpl17","Rpl18","Rpl18a",
"Rpl21","Rpl22","Rpl23","Rpl23a","Rpl24","Rpl27","Rpl28","Rpl29","Rpl3","Rpl31",
"Rpl34","Rpl35","Rpl36a","Rpl38","Rpl4","Rpl5","Rpl6","Rpl7","Rpl7a","Rpl8",
"Rpl9","Rplp0","Rplp2","Rpn1","Rpn2","Rprd1b","Rps10","Rps11","Rps12","Rps13",
"Rps14","Rps15a","Rps16","Rps17","Rps18","Rps19","Rps2","Rps20","Rps21","Rps23",
"Rps24","Rps25","Rps26","Rps27","Rps3","Rps4x","Rps5","Rps6","Rps6ka3","Rps7",
"Rps8","Rps9","Rpsa")

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

mrna_binding_genes <- c(
  "Elavl1",
  "Hnrnpa0","Hnrnpa1","Hnrnpa2b1","Hnrnpa3","Hnrnpab","Hnrnpc","Hnrnpd",
  "Hnrnpf","Hnrnph1","Hnrnpk","Hnrnpl","Hnrnpll","Hnrnpm","Hnrnpu","Hnrnpul2",
  "Ptbp1","Syncrip","Ybx1","Snd1","Caprin1","G3bp1","Cirbp","Ncl","Nono",
  "Srsf1","Srsf3","Srsf5","Srsf6","Srsf7","Tra2b","Sfpq",
  "Tardbp","Khdrbs1","Khsrp","Rbm39","Rbmxl1","Pabpc1","Pabpn1","Pcbp1","Pcbp2",
  "Serbp1","Pura","Purb","Srrt","Rtcb","Raly",
  "Ddx1","Ddx17","Ddx19a","Ddx21","Ddx39b","Ddx3x","Ddx5","Ddx6",
  "Dhx8","Dhx9","Dhx15",
  "Sf3a1","Sf3a3","Sf3b1","Sf3b3","U2af2","U2surp",
  "Snrnp200","Snrnp70","Snrpd1","Snrpd2","Snu13",
  "Prpf6","Prpf8","Prpf19"
)

membrane <- read.delim('/Users/andrewleduc/Desktop/Projects/AliveDead/SubCellLoc/Membrane.tsv')

df_cond_plot <- df_cond
df_cond_plot$var <- NULL
df_cond_plot$numb <- NULL
rownames(df_cond_plot) <- df_cond_plot$gene
df_cond_plot$gene <- NULL

df_cond_plot <- as.matrix(df_cond_plot)




df_core_bio <- rbind(colMedians(df_cond_plot[proteosome_core_complex,],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,chromatin_regulators),],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,membrane$Gene.Names),],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,ribo),],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,mito),],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,mrna_binding_genes),],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,e3_ligase),],na.rm=T),
                     colMedians(df_cond_plot[intersect(df_cond$gene,glycolysis),],na.rm=T)
)




rownames(df_core_bio) <- c('Proteosome','Chromatin binding', 'Membrane','Ribosome','Mitochondrial','mRNA binding','E3 ligase', 'Glycolysis')
colnames(df_core_bio) <- c('mRNA','Translation','Clearance half life')

col_fun <- colorRamp2(c( 0, 1), c("white", "red"))
Heatmap(df_core_bio,col = col_fun,cluster_columns = F)


################ ################ ################ ################ 
####  Regulation for individual cell type heatmap 
################ ################ ################ ################ 



#reg_heatmap('Basal')
reg_heatmap('Chondrocyte')
#reg_heatmap('Fibroblast')
colnames(df_posterior_mRNA)

i=2
df_clust <- cbind(df_posterior_mRNA[,i],df_posterior_trans[,i],-df_posterior_deg[,i],df_posterior_prot[,i]) 
colnames(df_clust) <- c('mRNA','Trans','Deg','Protein')
df_clust <- df_clust[rowSums(is.na(df_clust)==F)==4,]
df_clust_mRNA <- df_clust#[df_clust[,1]>0,]
df_clust_mRNA <- df_clust_mRNA[abs(df_clust_mRNA[,1])<.5,]
df_clust_mRNA <- df_clust_mRNA[abs(df_clust_mRNA[,2])>.5,]
df_clust_mRNA <- df_clust_mRNA[abs(df_clust_mRNA[,3])<.5,]
gene_clust <- rownames(df_clust_mRNA)
i=4
df_clust <- cbind(df_posterior_mRNA[,i],df_posterior_trans[,i],-df_posterior_deg[,i],df_posterior_prot[,i]) 
colnames(df_clust) <- c('mRNA','Trans','Deg','Protein')
df_clust <- df_clust[gene_clust,]
df_clust <- df_clust[rowSums(is.na(df_clust)==F)==4,]
Heatmap(df_clust,          # your 4‑column matrix
        cluster_rows = TRUE, cluster_columns = FALSE)

ht <- Heatmap(df_clust,          # your 4‑column matrix
              cluster_rows = TRUE, cluster_columns = FALSE)

# draw() is needed to materialise the clustering
ht <- draw(ht)

# extract the hclust object for rows
row_hc <- row_dend(ht) |> as.hclust()

# cut into k clusters (choose k by eye or algorithmically)
k <- 5
cl <- cutree(row_hc, k = k)          # named vector: gene → cluster
cl[names(cl)==3] <- 0
cl[names(cl) %in% rownames(df_clust1)] = 4

cl <- cl[cl!=0]
# 1 = mrna
# 2 = deg
# 3 = trans
# 4 = trans+mRNA

df_clust2 <- df_clust[names(cl),]
Heatmap(df_clust2,cluster_columns = F,col = col_fun,show_row_dend = F,
        cluster_rows = T,show_row_names = T)

################ ################ ################ ################ 
#### Go terms analysis cell type specific regulation
################ ################ ################ ################ 


## ---- parameters ------------------------------------------------------------
cell_types <- colnames(df_posterior_mRNA)   # adjust if needed
thr_high   <- 0.5                           # |Δ| above this → “high”
thr_low    <- 0.5                           # |Δ| below this → “low”


## ---- iterate over cell types ----------------------------------------------
modal_sets_by_ct <- setNames(
  lapply(seq_along(cell_types), get_modality_sets),
  cell_types
)

## Split into three “big” lists: one per modality -----------------------------
mRNA_clusters  <- lapply(modal_sets_by_ct, `[[`, "mRNA")
trans_clusters <- lapply(modal_sets_by_ct, `[[`, "Translation")
clear_clusters <- lapply(modal_sets_by_ct, `[[`, "Clearance")


mRNA_entrez  <- lapply(mRNA_clusters,  sym2entrez)
trans_entrez <- lapply(trans_clusters, sym2entrez)
clear_entrez <- lapply(clear_clusters, sym2entrez)


ego_mRNA <- compareCluster(
  geneClusters = mRNA_entrez,
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE            # converts Entrez back to SYMBOL in output
)

ego_trans <- compareCluster(
  geneClusters = trans_entrez,
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_clear <- compareCluster(
  geneClusters = clear_entrez,
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)




plot_GO_dot_clustered(ego_mRNA,  "mRNA-driven genes",
                       df_posterior_mRNA)

plot_GO_dot_clustered(ego_trans,  "Trans-driven genes",
                       df_posterior_trans)

plot_GO_dot_clustered(ego_clear,  "Clearance-driven genes",
                       df_posterior_deg)




############################################
# Doing single gene plots now
############################################


# Krt5, Krt18, Krt75, Krt15, Krt14, Krt19, Krt8 Krt4, Krt13 Krt7 Krt19
View(df_posterior_prot)
View(df_posterior_deg)

colnames(df_posterior_prot)
ct = colnames(df_posterior_prot)[1]
gene_add = 'S100a1' 

waterfall_arrow_plot(
  cell_type          = ct,
  protein_name       = gene_add,
  delta_mrna         = df_posterior_mRNA[gene_add,ct],df_posterior_mRNA_w[gene_add,ct],
  delta_translation  = df_posterior_trans[gene_add,ct], df_posterior_trans_w[gene_add,ct],
  delta_degradation  = -df_posterior_deg[gene_add,ct],df_posterior_deg_w[gene_add,ct]
)

df_posterior_prot[gene_add,ct]

View(df_posterior_prot)

df_look <- data.frame(prot = df_posterior_prot[,1],trans = df_posterior_trans[,1], gene = rownames(df_posterior_prot))
df_look <- df_look %>% filter(prot > .5 & trans > .5)

cors_rna <- c()
cors_trans <- c()
cors_deg <- c()
for(i in 1:nrow(df_posterior_mRNA)){
  cors_rna <- c(cors_rna,cor(df_posterior_mRNA[i,],df_posterior_prot[sample(nrow(df_posterior_prot),1),]))
  cors_trans <- c(cors_trans,cor(df_posterior_trans[i,],df_posterior_prot[sample(nrow(df_posterior_prot),1),]))
  cors_deg <- c(cors_deg,cor(df_posterior_deg[i,],df_posterior_prot[sample(nrow(df_posterior_prot),1),]))
}






################################################################################################################################################################

################ ################ ################ ################ 
####  plot model vs empirical for sanity check
################ ################ ################ ################ 


##### mRNA epirical vs model
df_emp <- mRNA_dat %>%
  dplyr::group_by(gene, Cell_type) %>%
  dplyr::summarize(emp_mean = mean(value), .groups = "drop")

rna_draws   <- fit_vi$draws("logM")     # draws × (G*C)
pm_rna      <- colMeans(as.matrix(rna_draws))

df_post_rna <- tibble(variable = names(pm_rna),
                      post_mean = pm_rna) %>%
  extract(
    variable,
    into = c("gene_idx", "ct_idx"),
    regex = "logM\\[([0-9]+),([0-9]+)\\]",
    convert = TRUE
  ) %>%
  mutate(
    gene      = all_genes[gene_idx],
    cell_type = all_celltypes[ct_idx]
  ) %>%
  dplyr::select(gene, cell_type, post_mean)


df_emp$cell_type <- df_emp$Cell_type
# 3) join empirical & posterior
df_cmp <- df_emp %>%
  inner_join(df_post_rna, by = c("gene", "cell_type"))


# 4) plot
ggplot(df_cmp, aes(x = emp_mean, y = post_mean)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  labs(
    x = "Empirical mean UMI count",
    y = "Posterior mean logM",
    title = "Empirical vs. Posterior mRNA levels"
  ) +
  theme_classic(base_size = 14)


mat_rna_model <- dcast(df_cmp,gene~Cell_type,value.var = 'post_mean')
mat_rna_emp <- dcast(df_cmp,gene~Cell_type,value.var = 'emp_mean')




##### protein epirical vs model
df_emp_prot <- protien_dat %>%
  dplyr::group_by(split_gene, cell_type) %>%
  dplyr::summarize(emp_mean_prot = sum(regressed_p_obs * n_obs, na.rm = TRUE) / sum(n_obs, na.rm = TRUE), .groups = "drop")


prot_draws <- fit_vi$draws("logP")
pm_prot    <- colMeans(as.matrix(prot_draws))
sd_prot    <- colSds(as.matrix(prot_draws))
#off_draws  <- fit_vi$draws("gene_offset_prot")
#pm_off     <- colMeans(as.matrix(off_draws))


# 3) Build df_post_prot and subtract the offset by gene
df_post_prot <- tibble(
  variable  = names(pm_prot),
  post_mean = pm_prot,
  post_width = sd_prot
) %>%
  extract(
    variable,
    into = c("gene_idx", "ct_idx"),
    regex   = "logP\\[([0-9]+),([0-9]+)\\]",
    convert = TRUE
  ) %>%
  mutate(
    gene        = all_genes[gene_idx],
    cell_type   = all_celltypes[ct_idx],         # pick the matching offset
    post_adj    = post_mean,           # subtract it here
    post_width = post_width
  ) %>%
  dplyr::select(gene, cell_type, post_adj,post_width)

# 4) Join empirical & offset-corrected posterior
df_cmp_prot <- df_emp_prot %>%
  inner_join(
    df_post_prot,
    by = c("split_gene" = "gene", "cell_type")
  )

# 5) Plot: empirical vs. offset-corrected posterior
ggplot(df_cmp_prot, aes(x = emp_mean_prot, y = post_adj)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  labs(
    x = "Empirical mean protein (centered)",
    y = "Posterior mean logP",
    title = "Empirical vs. Offset-Corrected Posterior Protein"
  ) +
  theme_classic(base_size = 14)

df_cmp_prot$emp_mean_prot[df_cmp_prot$post_width > .8] <- NA
df_cmp_prot$post_adj[df_cmp_prot$post_width > .8] <- NA

mat_prot_model <- dcast(df_cmp_prot,split_gene~cell_type,value.var = 'post_adj')
mat_prot_emp <- dcast(df_cmp_prot,split_gene~cell_type,value.var = 'emp_mean_prot')


##### clearance epirical vs model
df_emp_deg <- clearance_dat %>%
  group_by(split_gene, Cell_Type) %>%
  summarise(emp_mean_prot = sum(value * n_obs, na.rm = TRUE) /
              sum(n_obs, na.rm = TRUE),
            .groups = "drop")

#Posterior means *and* widths
draws_mat <- as_draws_matrix(fit_vi$draws("logDeg"))

post_mean  <- colMeans(draws_mat)
post_width <- apply(draws_mat, 2, \(x)
                    quantile(x, 0.975) - quantile(x, 0.025))          # 95 % width

df_post_deg <- tibble(
  variable   = names(post_mean),
  post_mean  = post_mean,
  post_width = post_width
) %>%
  extract(variable,
          into = c("gene_idx", "ct_idx"),
          regex   = "logDeg\\[([0-9]+),([0-9]+)\\]",
          convert = TRUE) %>%
  mutate(
    gene      = all_genes[gene_idx],
    cell_type = all_celltypes[ct_idx]
  ) %>%
  dplyr::select(gene, cell_type, post_mean, post_width)

#3Join empirical & posterior
df_cmp_deg <- df_emp_deg %>%
  inner_join(df_post_deg,
             by = c("split_gene" = "gene", "Cell_Type" = "cell_type"))

#4Scatter, coloured by posterior width
ggplot(df_cmp_deg,
       aes(x = emp_mean_prot,
           y = post_mean,
           color = post_width)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50") +
  scale_colour_viridis_c(name = "95% width",
                         option = "plasma") +
  labs(
    x     = "Empirical degradation fold change",
    y     = "Posterior mean degradation fold change",
    title = "Posterior means vs. empirical means\n(Color = posterior 95% width)"
  ) +
  theme_classic(base_size = 14)

df_cmp_deg$post_mean[df_cmp_deg$post_width > 2] <- NA
df_cmp_deg$emp_mean_prot[df_cmp_deg$post_width > 2] <- NA

mat_deg_model <- dcast(df_cmp_deg,split_gene~Cell_Type,value.var = 'post_mean')
mat_deg_emp <- dcast(df_cmp_deg,split_gene~Cell_Type,value.var = 'emp_mean_prot')


mat_prot_model[mat_prot_model > 2] <- NA
mat_prot_model[mat_prot_model < -2] <- NA

plot(mat_prot_model$Secratory , 
     mat_deg_model$Secratory)



