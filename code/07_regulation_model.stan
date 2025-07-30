data {
  //--- Protein data
  int<lower=0> M_pro;
  int<lower=1> G;               // number of genes
  int<lower=1> C;               // number of cell types
  array[M_pro] int<lower=1, upper=G> gene_id_pro;
  array[M_pro] int<lower=1, upper=C> celltype_id_pro;
  vector[M_pro] bar_pro;
  array[M_pro] int<lower=1> n_pro;


  //--- mRNA data
  int<lower=0> M_mrna;
  array[M_mrna] int<lower=1, upper=G> gene_id_mrna;
  array[M_mrna] int<lower=1, upper=C> celltype_id_mrna;
  vector[M_mrna] bar_mrna;
  //array[M_mrna] int<lower=1> n_rna;

  //--- Degradation-rate data
  int<lower=0> M_deg;
  array[M_deg] int<lower=1, upper=G> gene_id_deg;
  array[M_deg] int<lower=1, upper=C> celltype_id_deg;
  vector[M_deg] bar_deg;
  array[M_deg] int<lower=1> n_deg;
}

parameters {
  matrix[G, C] logM;
  matrix[G, C] logDeg;
  matrix[G, C] T;

  // noise parms …
  //vector<lower=0>[C] alpha_pro;
  //vector<lower=0>[C] alpha_mrna;
  //vector<lower=0>[C] alpha_deg;
  //real<lower=0> scale_pro;
  //real<lower=0> scale_mrna;
  //real<lower=0> scale_deg;

  // **two gene-level offsets**
  //vector[G] gene_offset_prot;
  //vector[G] gene_offset_deg;
  //vector[G] gene_offset_M;

  //matrix<lower=0>[G,C] sigma_pro;
  //matrix<lower=0>[G,C] sigma_mrna;
  //matrix<lower=0>[G,C] sigma_deg;
}
transformed parameters {
  matrix[G, C] logP;
  //matrix[G,C] logDeg_true;
  //matrix[G,C] logM_true;
  
  
  
  for (g in 1:G) {
    for (c in 1:C) {
      // incorporate both offsets here, without clobbering logDeg
      //logDeg_true[g,c] = logDeg[g,c] + gene_offset_deg[g];
      // gene_offset_prot[g]
      logP[g,c] =  logM[g,c]
                +      T[g,c]
                - logDeg[g,c];
    }
  }
}
model {
  // 1) Priors on levels
  to_vector(logM)   ~ normal(0, 5);
  to_vector(logDeg) ~ normal(0, 5);
  to_vector(T)      ~ normal(0, 5);

  // 2) Priors on your two offsets
  //gene_offset_prot  ~ normal(0, 1);
  //gene_offset_deg   ~ normal(0, 1);

  // 3) Hierarchical noise priors (unchanged)
  //alpha_pro  ~ cauchy(0, 1);
  //alpha_mrna ~ cauchy(0, 1);
  //alpha_deg  ~ cauchy(0, 1);

  //scale_pro  ~ cauchy(0, 1);
  //scale_mrna ~ cauchy(0, 1);
  //scale_deg  ~ cauchy(0, 1);

  //for (g in 1:G) {
  //  for (c in 1:C) {
  //    sigma_pro[g,c]  ~ normal(alpha_pro[c],  scale_pro)  T[0,];
  //    sigma_mrna[g,c] ~ normal(alpha_mrna[c], scale_mrna) T[0,];
  //    sigma_deg[g,c]  ~ normal(alpha_deg[c],  scale_deg)  T[0,];
  //  }
  //}

  // 4) Likelihood (unchanged, except uses our new logP)
  for (m in 1:M_pro) {
    int g = gene_id_pro[m];
    int c = celltype_id_pro[m];
    real se = 3 / sqrt(n_pro[m]);
    bar_pro[m] ~ normal(logP[g,c], se);
  }
  for (m in 1:M_mrna) {
    int g = gene_id_mrna[m];
    int c = celltype_id_mrna[m];
    real se = 3.0 / sqrt(20); // Can scale proportionally to the number of cells without using the absolute number
                              // can make "3" a parameter for each modality 
    bar_mrna[m] ~ normal(logM[g,c], se);
  }
  for (m in 1:M_deg) {
    int g = gene_id_deg[m];
    int c = celltype_id_deg[m];
    real se = 3 / sqrt(n_deg[m]);
    bar_deg[m] ~ normal(logDeg[g,c], se);
  }
}

