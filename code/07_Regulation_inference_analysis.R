library(cmdstanr)
library(posterior)


######################################################
# Simulate
######################################################

set.seed(123)

G <- 10   # number of genes
C <- 5   # number of cell types

# We'll simulate logM, logDeg, T, but since we have no mRNA data in the final,
# the model might still have logM parameters. We'll fill them in for completeness.

logM <- matrix(rnorm(G*C, 0, 0.5), nrow=G, ncol=C)    # log-mRNA (unused if M_mrna=0)
logDeg <- matrix(rnorm(G*C, -1, 0.4), nrow=G, ncol=C) # log-deg
Tmat   <- matrix(rnorm(G*C, 0, 0.3), nrow=G, ncol=C)  # offset
logP   <- logM + Tmat - logDeg   # exact relationship: logP = logM + T - logDeg




# We'll have gene×celltype-specific noise for protein and deg measurements
sigma_pro <- matrix(abs(rnorm(G*C, 0.3, 0.1)), nrow=G, ncol=C)
sigma_deg <- matrix(abs(rnorm(G*C, 0.2, 0.05)), nrow=G, ncol=C)
sigma_M <- matrix(abs(rnorm(G*C, 0.2, 0.05)), nrow=G, ncol=C)

# Now let's build the "observed" data sets:

XX = 7

# 1) Protein
# We'll measure M_pro = G*C groups (one measurement per gene×celltype).
M_pro <- G*C*XX
gene_id_pro     <- integer(M_pro)
celltype_id_pro <- integer(M_pro)
bar_pro         <- numeric(M_pro)
n_pro           <- integer(M_pro)  # number of observations used for each group

idx <- 1
for(i in 1:XX){
  for(g in 1:G) {
    for(c in 1:C) {
      gene_id_pro[idx]     <- g
      celltype_id_pro[idx] <- c
      # simulate one group measurement:
      n_pro[idx]           <- sample(2:5, 1)  # e.g. random group size 2..5
      # The group mean is roughly Normal(logP[g,c], sigma_pro[g,c]/sqrt(n_pro[idx]))
      # but let's just do a single measurement as if n_pro=1 => we can do:
      # For realism, let's do a single "collapsed" data point:
      x = 0
      if(g == 1){
        x = 2
      }
      bar_pro[idx] <- rnorm(1, mean=logP[g,c], sd=sigma_pro[g,c]) + x
      idx <- idx + 1
    }
  }
}


# 2) Degradation
# We'll measure M_deg = G*C as well for simplicity
M_deg <- G*C*XX
gene_id_deg     <- integer(M_deg)
celltype_id_deg <- integer(M_deg)
bar_deg         <- numeric(M_deg)
n_deg           <- integer(M_deg)

idx <- 1
for(i in 1:XX){
  for(g in 1:G) {
    for(c in 1:C) {
      gene_id_deg[idx]     <- g
      celltype_id_deg[idx] <- c
      n_deg[idx]           <- sample(2:5, 1)
      bar_deg[idx] <- rnorm(1, mean=logDeg[g,c], sd=sigma_deg[g,c])
      idx <- idx + 1
    }
  }
}


# mRNA data is not used => M_mrna=0
M_mrna <- G*C*XX
gene_id_mrna     <- integer(M_mrna)
celltype_id_mrna <- integer(M_mrna)
bar_mrna         <- numeric(M_mrna)

idx <- 1
for(i in 1:XX){
  for(g in 1:G) {
    for(c in 1:C) {
      gene_id_mrna[idx]     <- g
      celltype_id_mrna[idx] <- c
      bar_mrna[idx] <- rnorm(1, mean=logM[g,c], sd=sigma_deg[g,c])
      idx <- idx + 1
    }
  }
}



stan_data <- list(
  G = G,
  C = C,

  # Protein
  M_pro = M_pro,
  gene_id_pro = gene_id_pro,
  celltype_id_pro = celltype_id_pro,
  bar_pro = bar_pro,
  n_pro = n_pro,

  # mRNA
  M_mrna = M_mrna,
  gene_id_mrna = gene_id_mrna,
  celltype_id_mrna = celltype_id_mrna,
  bar_mrna = bar_mrna,

  # Deg
  M_deg = M_deg,
  gene_id_deg = gene_id_deg,
  celltype_id_deg = celltype_id_deg,
  bar_deg = bar_deg,
  n_deg = n_deg
)




######################################################
# Plug real data
######################################################

stan_data <- list(
  G = length(unique(rna_dat$gene)),
  C = length(unique(rna_dat$Cell_type)),
  
  # Protein
  M_pro = nrow(data_protein_model_input),
  gene_id_pro = data_protein_model_input$gene_id_pro,
  celltype_id_pro = data_protein_model_input$celltype_id_pro,
  bar_pro = data_protein_model_input2$regressed_p_obs,
  n_pro = data_protein_model_input$n_obs,
  
  # mRNA
  M_mrna = nrow(rna_dat),
  gene_id_mrna = rna_dat$gene_id_mrna,
  celltype_id_mrna = rna_dat$celltype_id_mrna,
  bar_mrna = rna_dat2$value,
  
  # Deg
  M_deg = nrow(data_degradation_model_input),
  gene_id_deg = data_degradation_model_input$gene_id_deg,
  celltype_id_deg = data_degradation_model_input$celltype_id_deg,
  bar_deg = data_degradation_model_input2$value,
  n_deg = data_degradation_model_input$n_obs
)









#------------------------------------------------------------------------------
# 3) Fit the main model in ~/model.stan
#------------------------------------------------------------------------------
# We'll assume your final model has the correct logic for no mRNA, etc.
# E.g. "cmdstan_model('~/model.stan')"
# But for demonstration, let's say the file is actually at:
model_file <- file.path("~", "/Desktop/Projects/Miceotopes/SingleCell/Model.stan")

mod_main <- cmdstan_model(model_file)

# Run MCMC
fit_main <- mod_main$sample(
  data = stan_data,
  seed = 121,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000
)
print(fit_main$summary())

look <- fit_main$draws("gene_offset")
summarise_draws(look, "median")

draws_mcmc <- fit_main$draws("logP")



# Use .width = 0.95 to get 95% intervals (which will yield columns q2.5 and q97.5)
logP_summary_mcmc <- summarise_draws(
  draws_mcmc,
  ~quantile2(.x, probs = c(0.025, 0.975))
)

# Extract parameter names and indices using stringr
params <- logP_summary_mcmc$variable
matches <- str_match(params, "logP\\[(\\d+),(\\d+)\\]")  # capture gene and cell type indices

# Create a data frame with gene, celltype, median, and 95% CI bounds:
logP_df_mcmc <- as.data.frame(logP_summary_mcmc) %>%
  mutate(
    gene = as.integer(matches[,2]),
    celltype = as.integer(matches[,3])
  ) %>%
  dplyr::select(gene, celltype,
                lower_logP = q2.5, upper_logP = q97.5)

logP_df_mcmc$median <- as.data.frame(summarise_draws(draws_mcmc, "median"))$median


# For comparison, suppose you have empirical protein means from your data:
emp_df <- data.frame(
  gene = stan_data$gene_id_pro,
  celltype = stan_data$celltype_id_pro,
  empirical_logP = stan_data$bar_pro
) %>%
  group_by(gene, celltype) %>%
  summarise(empirical_logP = mean(empirical_logP), .groups = "drop")

plot_df_mcmc <- inner_join(emp_df, logP_df_mcmc, by = c("gene", "celltype"))

true_logP_df <- expand.grid(gene = 1:G, celltype = 1:C) %>%
  mutate(true_logP = as.vector(logP))


# --- Join the posterior summary with the true values by gene and cell type.
plot_true_df_MCMC <- inner_join(true_logP_df, logP_df_mcmc, by = c("gene", "celltype"))

# Join the VI estimates with empirical data by gene and cell type:

# Plot VI estimates vs. empirical means with 95% credible intervals as error bars:
ggplot(plot_true_df_MCMC, aes(x = true_logP, y = median)) +
  geom_point(size = 2, color = "blue") +
  geom_errorbar(aes(ymin = lower_logP, ymax = upper_logP), width = 0.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +

  labs(x = "Simulated logP",
       y = "Full MCMC Estimated logP (median with 95% CI)",
       title = "Comparison of Simulated vs. MCMC Estimated Protein Means") +
  theme_minimal()



mod_main <- cmdstan_model(model_file)
fit_vi <- mod_main$variational(
  data = stan_data,
  seed = 123,
  algorithm = "meanfield"
)


fit_vi <- mod_main$variational(
  data = stan_data,
  seed = 123,
  algorithm = "meanfield",
  iter = 20000,        # increase total iterations
  grad_samples = 3,    # number of gradient evals for each iteration
  elbo_samples = 100,  # number of samples for ELBO evaluation
  tol_rel_obj = 1e-3   # convergence tolerance
)




View(fit_vi$draws)

draws_vi_prot <- (fit_vi$draws("logP"))
draws_vi_rna <- (fit_vi$draws("logM"))
draws_vi_deg <- (fit_vi$draws("logDeg"))
draws_vi_T <- (fit_vi$draws("T"))

hist(draws_vi_prot[,3])

df <- data.frame(prot = colSds(draws_vi_prot),mRNA = colSds(draws_vi_rna),
           deg = colSds(draws_vi_deg),trans = colSds(draws_vi_T))
df <- melt(df)

ggplot(df, aes(x = variable,y = log10(value))) + geom_boxplot() + dot_plot+
  ylab('Posterior widths, log10(sd)') +xlab('')


plot(colMeans(draws_vi_rna),colMeans(draws_vi_prot),main = 'Pearson, 0.74')
abline(a=0,b=1,col='red')

cor(colMeans(draws_vi_prot)-colMeans(draws_vi_rna),colMeans(draws_vi_deg))
plot(colMeans(draws_vi_prot)-colMeans(draws_vi_rna),colMeans(draws_vi_T), main = 'Translation, cor = 0.8')
0.4479305 * cor(colMeans(draws_vi_prot)-colMeans(draws_vi_rna)-colMeans(draws_vi_T),colMeans(draws_vi_deg))^2


plot(colMeans(draws_vi_prot)-colMeans(draws_vi_rna),colMeans(draws_vi_deg))










df_emp <- rna_dat2 %>%
  dplyr::group_by(gene, Cell_type) %>%
  dplyr::summarize(emp_mean = mean(value), .groups = "drop")

# 2) posterior means of logM
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

hist(df_cmp$post_mean)

df_cmp$diff <- df_cmp$emp_mean - df_cmp$post_mean

View(df_cmp)

rna_dat_test <- rna_dat %>% filter(Cell_type == 'Chondrocyte' & gene == 'Cbr2')

View(rna_dat_test)



df_Df <- df_emp_prot %>% left_join(df_emp, by = c('split_gene'='gene','cell_type'= 'Cell_type')) 
plot(df_Df$emp_mean_prot,df_Df$emp_mean)
abline(a=0,b=1)
cor(df_Df$emp_mean_prot,df_Df$emp_mean,method='spearman')





# 1) Empirical protein means
df_emp_prot <- data_protein_model_input2 %>%
  dplyr::group_by(split_gene, cell_type) %>%
  dplyr::summarize(emp_mean_prot = sum(regressed_p_obs * n_obs, na.rm = TRUE) / sum(n_obs, na.rm = TRUE), .groups = "drop")


prot_draws <- fit_vi$draws("logP")
pm_prot    <- colMeans(as.matrix(prot_draws))

#off_draws  <- fit_vi$draws("gene_offset_prot")
#pm_off     <- colMeans(as.matrix(off_draws))


# 3) Build df_post_prot and subtract the offset by gene
df_post_prot <- tibble(
  variable  = names(pm_prot),
  post_mean = pm_prot
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
    post_adj    = post_mean           # subtract it here
  ) %>%
  dplyr::select(gene, cell_type, post_adj)

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


sd(df_cmp_prot$post_adj)




df_emp_prot <- data_degradation_model_input2 %>%
  group_by(split_gene, Cell_Type) %>%
  summarise(emp_mean_prot = sum(value * n_obs, na.rm = TRUE) /
              sum(n_obs, na.rm = TRUE),
            .groups = "drop")

#Posterior means *and* widths
draws_mat <- as_draws_matrix(fit_vi$draws("logDeg"))

post_mean  <- colMeans(draws_mat)
post_width <- apply(draws_mat, 2, \(x)
                    quantile(x, 0.975) - quantile(x, 0.025))          # 95 % width

df_post_prot <- tibble(
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
df_cmp_prot <- df_emp_prot %>%
  inner_join(df_post_prot,
             by = c("split_gene" = "gene", "Cell_Type" = "cell_type"))

#4Scatter, coloured by posterior width
ggplot(df_cmp_prot,
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

cor(df_cmp_prot$post_mean,df_cmp_prot$emp_mean_prot)











library(posterior)   # as_draws_matrix()
library(tidyr)       # pivot_longer()
library(dplyr)       # joins & pipes
library(tibble)      # rownames_to_column()

# helper ──────────────────────────────────────────────────────────────
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

# make sure you have vectors 'all_genes' and 'all_celltypes' used
# when you built stan_data
# --------------------------------------------------------------------

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

# join them all together ─────────────────────────────────────────────
df_posterior <- df_prot  %>%
  left_join(df_deg,   by = c("gene", "cell_type")) %>%
  left_join(df_trans, by = c("gene", "cell_type")) %>%
  left_join(df_mrna,  by = c("gene", "cell_type"))

# inspect



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

## ---------- build one width‑column for each modality --------------------------
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

## ---------- join them together (same pattern you used) ------------------------
df_posteriors_width <- df_prot_w %>%
  left_join(df_deg_w,   by = c("gene", "cell_type")) %>%
  left_join(df_trans_w, by = c("gene", "cell_type")) %>%
  left_join(df_mrna_w,  by = c("gene", "cell_type"))

df_posteriors_width_Cyp2b10 <- df_posteriors_width %>% filter(gene == 'Sult1d1')
df_posteriors_Cyp2b10 <- df_posterior %>% filter(gene == 'Sult1d1')

df_posterior <- df_posterior[rownames(protein_mat),colnames(rownames(protein_mat))]
df_posterior[is.na(protein_mat)] <- NA


mrna_plot('Sult1d1')
trans_plot2('Sec14l2')
prot_mrna_deg_plot('Tuba4a')


mrna_plot <- function(gene_in){
  df_plot <- df_posteriors_Cyp2b10 %>%                       # means
    left_join(df_posteriors_width_Cyp2b10,                   # widths
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

names(palette)

deg_plot = function(){
df_plot <- df_posteriors_Cyp2b10 %>%                       # means
  left_join(df_posteriors_width_Cyp2b10,                   # widths
            by = c("gene", "cell_type")) %>%
  ## keep only the two modalities we want to show
  dplyr::select(cell_type,
         protein_mean       = protein_abundance,
         deg_mean          = deg,
         protein_width,
         deg_width)

df_plot$deg_mean <- -df_plot$deg_mean

## ── 2.Long format: one row per 
df_plot_long <- df_plot %>%
  pivot_longer(
    cols      = c(protein_mean, deg_mean,
                  protein_width, deg_width),
    names_to  = c("modality", ".value"),
    names_sep = "_"
  ) %>%                         # .value makes “mean” and “width” real columns
  mutate(
    modality = recode(modality,
                      protein = "Protein",
                      deg    = "deg"),
    colour   = if_else(modality == "Protein", "orange", "green")
  )

## Decide what an “error bar” represents

df_plot_long <- df_plot_long %>%
  mutate(ymin = mean - width/2,
         ymax = mean + width/2)

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
                                 "deg"    = "green")) +
  labs(x = "Cell type",
       y = "Posterior mean (log‑scale)",
       colour = NULL,
       title = "Cyp2b10: protein vs deg across cell types") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

}


trans_plot2 <- function(gene_in) {
  ## ── 1 Filter the two master frames to this gene ──────────────────────────

  df_mean  <- df_posterior        %>% filter(gene == gene_in)
  df_width <- df_posteriors_width %>% filter(gene == gene_in)
  
  ## ── 2 Merge means + widths and keep the three modalities we want ─────────
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


prot_mrna_deg_plot <- function(gene_in) {
  
  ## ── 1 Filter to the chosen gene ──────────────────────────────────────────
  df_mean  <- df_posterior        %>% filter(gene == gene_in)
  df_width <- df_posteriors_width %>% filter(gene == gene_in)
  
  ## ── 2 Merge means + widths; keep the three modalities ────────────────────
  df_plot <- df_mean %>%
    left_join(df_width, by = c("gene", "cell_type")) %>%
    dplyr::select(cell_type,
           protein_mean  = protein_abundance,
           mRNA_mean     = mRNA,
           deg_mean      = deg,
           protein_width,
           mRNA_width,
           deg_width)
  
  df_plot$deg_mean <- -df_plot$deg_mean
  
  ## ── 3 Long format (pivot) ────────────────────────────────────────────────
  df_long <- df_plot %>%
    pivot_longer(
      cols      = c(protein_mean, mRNA_mean, deg_mean,
                    protein_width, mRNA_width, deg_width),
      names_to  = c("modality", ".value"),
      names_sep = "_"
    ) %>%
    mutate(
      modality = recode(modality,
                        protein = "Protein",
                        mRNA    = "mRNA",
                        deg     = "Degradation"),
      ymin = mean - width / 2,      # if width==SD use 2*width for ±2 SD
      ymax = mean + width / 2
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
                                   "mRNA"        = "#008000",
                                   "Degradation" = "purple")) +
    labs(x = "Cell type",
         y = "Posterior mean (log2 scale)",
         colour = NULL,
         title = paste0(gene_in, ": protein, mRNA and degradation across cell types")) +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

## Example:
# prot_mrna_deg_plot("Cyp2b10")

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

colnames(mRNA_mat)
i=7
df_clust <- cbind(mRNA_mat[,i],trans_mat[,i],-deg_mat[,i],protein_mat[,i]) 
colnames(df_clust) <- c('mRNA','Trans','Deg','Protein')
df_clust <- df_clust[rowSums(is.na(df_clust)==F)==4,]

d1 <- c(d1,df_clust1[,3])
d2 <- c(d2,rep('Immune',nrow(df_clust1)))

dfdfdf <- data.frame(deg = d1,ct = d2)
ggplot(dfdfdf,aes(x = ct,y = deg))+geom_boxplot()

df_clust1 <- df_clust#[df_clust[,1]>0,]
df_clust1 <- df_clust1[abs(df_clust1[,3])>.5,]
df_clust1 <- df_clust1[abs(df_clust1[,2])>.5,]
df_clust1 <- df_clust1[abs(df_clust1[,1])>.4,]

df_clust1 <- df_clust1[df_clust1[,4]>1,]



df_clust1 <- df_clust1[abs(df_clust1[,3])> abs(df_clust1[,1]),]

Heatmap(df_clust1,cluster_columns = F,col = col_fun,show_row_dend = F,
        cluster_rows = T,show_row_names = T)

df_clust1[,4] <- df_clust1[,4] - df_clust1[,1]
df_clust1 <- df_clust1[,2:4]

colnames(df_clust1)[3] <- 'Protein/mRNA ratio'

df_clust1 <- df_clust1[abs(df_clust1[,1]) > .2,]
df_clust1 <- df_clust1[abs(df_clust1[,3]) < .2,]


library(ComplexHeatmap)

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

gene_clusters <- split(names(cl), cl)  # list( cluster1 = c("GeneA","GeneB",…),
#       cluster2 = … )

# Optionally convert to Entrez if you started with symbols

gene_clusters <- lapply(gene_clusters, \(sym) {
  mapIds(org.Mm.eg.db,
         keys         = sym,
         keytype      = "SYMBOL",
         column       = "ENTREZID",
         multiVals    = "first") |>
    na.omit()
})

all_syms <- unique(unlist(gene_clusters, use.names = FALSE))
library(tibble)
## 2.  Map → Entrez and build a data‑frame
gene_map <- AnnotationDbi::select(org.Mm.eg.db,
                                  keys     = all_syms,
                                  keytype  = "ENTREZID",
                                  columns  = c("SYMBOL", "ENTREZID")) 




ego <- compareCluster(
  geneClusters = gene_clusters,
  fun          = "enrichGO",
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

head(ego@compareClusterResult[, c("ID", "Description", "p.adjust")], 100)

View(ego@compareClusterResult)

nrow(ego@compareClusterResult)

term <- "GO:0006749"                  # or "translation"

genes_term <- ego@compareClusterResult |>
  subset(ID == term | Description == term) |>
  (\(x) x[["geneID"]])() |>          # <- this is now a real function
  strsplit("/") |>
  unlist() |>
  unique()

gene_map2 <- gene_map %>% filter(ENTREZID %in% genes_term)

df_plot <- as.data.frame(df_clust[gene_map2$SYMBOL,])
df_plot$gene <- rownames(df_plot)
df_plot <- melt(df_plot,id.vars = 'gene')

ggplot(df_plot,aes(x = variable,y = value))+geom_point() + theme_classic()+
  facet_wrap(~gene)

df_plot <- df_plot %>% filter(variable != 'Trans')
df_plot$variable <- factor(df_plot$variable,
                           levels = c("mRNA", "Trans", "Protein"))

## 2 ─── Order panels by descending protein value ───────────────────────────
gene_order <- df_plot %>%
  filter(variable == "Protein") %>%
  arrange(desc(value)) %>%
  pull(gene)

df_plot$gene <- factor(df_plot$gene, levels = gene_order)

## 3 ─── Lollipop plot with a 0-baseline ────────────────────────────────────

ggplot(df_plot,aes(x = variable,y= value))+ geom_boxplot() + ylim(c(0,2.5))

ggplot(df_plot,
       aes(y = variable,
           x = value,
           fill = variable)) +
  
  ## baseline at 0 (drawn first so the stick sits on top)
  geom_vline(xintercept = 0, colour = "black", linewidth = 1) +
  
  geom_bar(stat = 'identity')+
  ## stick
  # geom_segment(aes(yend = variable, xend = 0),
  #              colour = "grey70", linewidth = 0.6) +
  # 
  # ## lollipop head
  # geom_point(size = 2.2) +
  
  facet_wrap(~gene, ncol = 11) +
  coord_flip() +                          # keeps the panels wider than tall
  labs(x = NULL, y = "Value") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        axis.text.x      = element_text(angle = 45, hjust = 1))+xlab('log2 fold change')+
  ylab('')


dotplot(ego, showCategory = 10)
new_names <- c("1" = "Transcription",
               "3" = "Translation",
               "4" = "Transcription + translation")

ego2 <- ego                                      # keep original intact
ego2@compareClusterResult$Cluster <-
  plyr::revalue(ego2@compareClusterResult$Cluster,
                new_names)

cnetplot(ego2, foldChange=, showCategory = 10)
ego2 = setReadable(ego2, 'org.Mm.eg.db', 'ENTREZID')

enrichplot::dotplot(ego2, showCategory = 5) +                       # your renamed object
  ggtitle("Basal") +                                     # add “Basal” to title
  theme(axis.text.x = element_text(                      # tilt x-axis labels
    angle = 45, hjust = 1, vjust = 1))

ekk <- compareCluster(
  geneClusters = gene_clusters,
  fun          = "enrichKEGG",
  organism     = "mmu",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dotplot(ekk, showCategory = 10)




df_posterior_trans <- Normalize_reference_vector_log(df_posterior_trans)

df_try <- data.frame(rna = df_posterior_mRNA[,1],synth = df_posterior_trans[,1])
df_try$col <- ''
df_try$col[rownames(df_try) %in% c('Krt5', 'Krt18', 'Krt75', 'Krt15', 'Krt14', 'Krt19', 'Krt8', 'Krt4', 'Krt13', 'Krt7', 'Krt19')] <- 'Krt'
ggplot(df_try,aes(x = rna,y = synth,colour = col,alpha = col))+ geom_point() + geom_abline(intercept = 0,slope = 0)+
  scale_alpha_manual(values = c(.3,1)) + ylab('Protein synthesis')+ xlab('mRNA abundance')+ dot_plot



View(df_posterior_prot)

cors_rna <- c()
cors_trans <- c()
cors_deg <- c()
for(i in 1:nrow(df_posterior_mRNA)){
  cors_rna <- c(cors_rna,cor(df_posterior_mRNA[i,],df_posterior_prot[i,]))
  cors_trans <- c(cors_trans,cor(df_posterior_trans[i,],df_posterior_prot[i,]))
  cors_deg <- c(cors_deg,cor(df_posterior_deg[i,],df_posterior_prot[i,]))
}

df_cond <- data.frame(mRNA = cors_rna,trans = cors_trans,deg = cors_deg)
df_cond$gene <- rownames(df_posterior_deg)

df_cond <- melt(df_cond)
ggplot(df_cond,aes(y = value))+ geom_histogram() + facet_wrap(~variable) +
  ylab('Correlation to protein abundance') + xlab('# Genes') + theme_bw(base_size = 15)

df_cond <- data.frame(mRNA = cors_rna,trans = cors_trans,deg = cors_deg,gene = toupper(rownames(df_posterior_mRNA)))
go_list <- read.delim('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/Gene_sets/GO_Human.txt',sep = ' ')

pvals <- c()
med_cor <- c()
term <- c()
for(i in unique(go_list$GO_term_name)){
  go_hold <- go_list %>% filter(GO_term_name == i)
  if(length(intersect(df_cond$gene,go_hold$Gene))>3){
    
    df_hold <- df_cond %>% filter(gene %in% go_hold$Gene)
    
    pvals <- c(pvals,t.test(df_hold$trans,cors_rna)$p.value)
    med_cor <- c(med_cor,median(df_hold$trans))
    term <- c(term,i)
    
  }
}

df_use <- data.frame(p_val = pvals,cor =med_cor,go =  term)



df_go_plot <- df_cond
# Deg terms
deg_terms <- c('respiratory electron transport chain',
'cytosolic small ribosomal subunit',
'mitotic nuclear envelope reassembly',
'glucose metabolic process',
'proteasome accessory complex')
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



View(df_go_plot)


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

ggplot(df_go_plot,aes(x = Value,y = term)) + ggbeeswarm::geom_beeswarm() +
  facet_wrap(~type,ncol = 1,scales = 'free_y') + xlab('Correlation')+
  ylab('') + theme_bw()



go_hold <- go_list %>% filter(GO_term_name == 'toxin metabolic process')
df_hold <- df_cond %>% filter(gene %in% go_hold$Gene)
ggplot(df_hold,aes(x = trans,y = 'cytosolic small ribosomal subunit')) + xlim(-1,1)+
  ggbeeswarm::geom_beeswarm()



install.packages('ggbeeswarm')
library(ggbeeswarm)
View(df_posterior_prot)
plot(df_posterior_prot[,7],df_posterior_mRNA[,7])
abline(a=0,b=1)
colnames(df_posterior_prot)
df_posterior_prot['Acan',]

mRNA_mat['Krt5',]

dim(df_posterior_prot)

# Vim, Des, Col6a2

# Acan -> secreded chondrocytes

#Vcl

# Krt5, Krt18, Krt75, Krt15, Krt14, Krt19, Krt8 Krt4, Krt13 Krt7 Krt19

ct = colnames(df_posterior_prot)[2]
gene_add = 'Pgls' 

waterfall_arrow_plot(
  cell_type          = ct,
  protein_name       = gene_add,
  delta_mrna         = df_posterior_mRNA[gene_add,ct],
  delta_translation  = df_posterior_trans[gene_add,ct],
  delta_degradation  = -df_posterior_deg[gene_add,ct]
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


# View(df_posterior)
# 
# library(qusage)
# go_list <- read.delim('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/Gene_sets/GO_Human.txt',sep = ' ')
# df_posterior_bas$set <- 'No'
# df_posterior_bas$gene_upper <- toupper(df_posterior_bas$gene)
# df_posterior_bas$set[df_posterior_bas$gene_upper %in% mito$Gene] <- 'Mito'
# 
# 
# unique(df_posterior$cell_type)
# df_posterior_bas <- df_posterior %>% filter(cell_type == 'Chondrocyte')
# df_posterior_bas$total <- abs(df_posterior_bas$deg) + abs(df_posterior_bas$translation)+
#   abs(df_posterior_bas$mRNA)
# 
# 

########

# parameters {
#   //--------------------------------------------
#     // 1) True log-levels per (gene, cell type)
# //--------------------------------------------
#   matrix[G, C] logM;    // log(mRNA)
# matrix[G, C] logDeg;  // log(degRate)
# 
# // Offset T[g,c], controlling Protein = mRNA + T - deg
# matrix[G, C] T;
# 
# // 2) Hierarchical noise parameters for each modality
# //    (protein, mRNA, deg), partially pooled across gene×celltype
# //    with cell-type hyperparameters
# 
# // (A) Cell-type baseline noise
# vector<lower=0>[C] alpha_pro;
# vector<lower=0>[C] alpha_mrna;
# vector<lower=0>[C] alpha_deg;
# 
# // (B) Global scale of gene-level deviations
# real<lower=0> scale_pro;
# real<lower=0> scale_mrna;
# real<lower=0> scale_deg;
# 
# // offset mean protein normalization
# vector[G] gene_offset_prot;
# vector[G] gene_offset_deg;
# 
# // (C) Actual gene×celltype noise
# matrix<lower=0>[G, C] sigma_pro;
# matrix<lower=0>[G, C] sigma_mrna;
# matrix<lower=0>[G, C] sigma_deg;
# }
# 
# transformed parameters {
#   // Exact relationship: logProtein[g,c] = logM[g,c] + T[g,c] - logDeg[g,c]
#   matrix[G, C] logP;
#   
#   
#   
#   for (g in 1:G) {
#     logDeg[g,] = logDeg[g,] + gene_offset_deg[g]
#     for (c in 1:C) {
#       logP[g,c] = gene_offset_prot[g] + (logM[g,c] + T[g,c] - logDeg[g,c]);
#     }
#   }
# }
# 
# model {
#   //----------------------------------------------------
#     // (A) Priors for the log-levels and offsets
#   //----------------------------------------------------
#     to_vector(logM)  ~ normal(0, 1.5);
#   to_vector(logDeg) ~ normal(0, 1);
#   to_vector(T)     ~ normal(0, 0.5);  // encourages T ≈ 0
#   
#   
#   gene_offset ~ normal(0, 1);
#   //----------------------------------------------------
#     // (B) Priors for hierarchical noise structures
#   //----------------------------------------------------
#     // Cell-type baseline
#   alpha_pro  ~ cauchy(0, 1);
#   alpha_mrna ~ cauchy(0, 1);
#   alpha_deg  ~ cauchy(0, 1);
#   
#   // Global scale for gene-level deviations
#   scale_pro  ~ cauchy(0, 1);
#   scale_mrna ~ cauchy(0, 1);
#   scale_deg  ~ cauchy(0, 1);
#   
#   // Partial pooling for sigma_*[g,c]
#   for (c in 1:C) {
#     for (g in 1:G) {
#       sigma_pro[g,c]  ~ normal(alpha_pro[c],  scale_pro)  T[0,];
#       sigma_mrna[g,c] ~ normal(alpha_mrna[c], scale_mrna) T[0,];
#       sigma_deg[g,c]  ~ normal(alpha_deg[c],  scale_deg)  T[0,];
#     }
#   }
#   
#   //----------------------------------------------------
#     // (C) Likelihood for each data stream
#   //----------------------------------------------------
#     
#     // Protein
#   for (m in 1:M_pro) {
#     int g = gene_id_pro[m];
#     int c = celltype_id_pro[m];
#     real se = sigma_pro[g,c] / sqrt(n_pro[m]);
#     bar_pro[m] ~ normal(logP[g,c], se);
#   }
#   
#   // mRNA
#   for (m in 1:M_mrna) {
#     int g = gene_id_mrna[m];
#     int c = celltype_id_mrna[m];
#     //real se = sigma_mrna[g,c] / sqrt(n_mrna[m]);
#     bar_mrna[m] ~ normal(logM[g,c], sigma_mrna[g,c]);
#   }
#   
#   // Degradation
#   for (m in 1:M_deg) {
#     int g = gene_id_deg[m];
#     int c = celltype_id_deg[m];
#     real se = sigma_deg[g,c] / sqrt(n_deg[m]);
#     bar_deg[m] ~ normal(logDeg[g,c], se);
#   }
# }
# 
# 
# 
# 
# 


library(ggplot2)
library(dplyr)

# --------------------------------------------------------------------
# Waterfall plot  (ONE protein × ONE cell type)
# --------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(grid)   # for arrow()

waterfall_arrow_plot <- function(cell_type,
                                 protein_name,
                                 delta_mrna,
                                 delta_translation,
                                 delta_degradation,
                                 baseline = 0) {
  
  # 1. contribution table with cumulative positions ---------------------------
  contrib_tbl <- tibble(
    effect = factor(c("mRNA", "Translation", "Degradation"),
                    levels = c("mRNA", "Translation", "Degradation")),
    delta  = c(delta_mrna, delta_translation, delta_degradation)
  ) %>% 
    mutate(x     = as.numeric(effect),                  # x ticks 1‒3
           start = baseline + c(0, head(cumsum(delta), -1)),
           end   = start + delta)
  
  total   <- baseline + sum(contrib_tbl$delta)
  total_x <- length(levels(contrib_tbl$effect)) + 1     # x‑pos for black line
  
  # 2. plot -------------------------------------------------------------------
  ggplot(contrib_tbl) +
    geom_hline(yintercept = baseline, linetype = "dashed") +
    
    # coloured arrows for each contribution

    geom_segment(
      aes(x     = x,
          xend  = x,
          y     = start,
          yend  = end,
          colour = effect),
      linewidth = 3,
      lineend   = "butt",      # no rounded caps
      linejoin  = "mitre",     # sharp corners where shaft meets head
      arrow = arrow(type   = "closed",
                    length = unit(0.35, "cm"),  # bigger head
                    angle  = 30)                # smaller angle → pointier tip
    )+


    # black line from baseline to final total
    annotate("segment",
             x = total_x, xend = total_x,
             y = baseline, yend = total,
             linewidth = 1.4, colour = "black") +
    annotate("text",
             x = total_x, y = total,
             label = sprintf("%.2f", total),
             vjust = -0.6, fontface = "bold") +
    
    # contribution labels near arrow tips (optional)
    geom_text(aes(x = x,
                  y = if_else(delta >= 0, end, start),
                  label = sprintf("%.2f", delta),
                  vjust = if_else(delta >= 0, -0.4, 1.2)),
              size = 4) +
    
    scale_colour_manual(values = c("mRNA"        = "#008000",
                                   "Translation" = "blue",
                                   "Degradation" = "purple")) +
    scale_x_continuous(breaks = c(1:3, total_x),
                       labels = c(levels(contrib_tbl$effect), "")) +
    labs(x = NULL,
         y = "Log2 protein fold change",
         title    = protein_name,
         subtitle = cell_type,
         colour   = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45,hjust=1,vjust=1))
}

## Example call:
waterfall_arrow_plot(
  cell_type          = "Chondrocytes",
  protein_name       = "Glycogen phosphorylase",
  delta_mrna         = 0.12,
  delta_translation  = 0.35,
  delta_degradation  = -0.27
)




