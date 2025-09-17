data {
  
  int<lower=1> G;   // number of genes
  int<lower=1> C;   // number of cells
  
  //--- Protein data
  int<lower=0> M_pro;        
  array[M_pro] int<lower=1, upper=G> gene_id_pro;
  array[M_pro] int<lower=1, upper=C> celltype_id_pro;
  vector[M_pro] bar_pro;
  array[M_pro] int<lower=1> n_pro; // Number non NA observations


  //--- mRNA data
  int<lower=0> M_mrna;
  array[M_mrna] int<lower=1, upper=G> gene_id_mrna;
  array[M_mrna] int<lower=1, upper=C> celltype_id_mrna;
  vector[M_mrna] bar_mrna;
  array[M_mrna] int<lower=1> n_rna; //// *** Number non 0 observations

  //--- Clearance rate data
  int<lower=0> M_deg;
  array[M_deg] int<lower=1, upper=G> gene_id_deg;
  array[M_deg] int<lower=1, upper=C> celltype_id_deg;
  vector[M_deg] bar_deg;
  array[M_deg] int<lower=1> n_deg; // Number non NA observations
  
}

parameters {
  matrix[G, C] logM; // average mRNA level for a cell type and cluster
  matrix[G, C] logDeg; // average clearance level for a cell type and cluster
  matrix[G, C] Trans; // average translation level for a cell type and cluster

}

transformed parameters {
  matrix[G, C] logP; // average protein level for a cell type and cluster
  
  
  // The protein data point
  for (g in 1:G) {
    for (c in 1:C) {
      
      logP[g,c] =  logM[g,c]
                +      Trans[g,c]
                - logDeg[g,c];
    }
  }
  
}

model {
  // 1) Priors on levels
  
  // Set these as parameters, the normal std... (see what they end up being)
  
  to_vector(logM)   ~ normal(0, 5);
  to_vector(logDeg) ~ normal(0, 3);
  to_vector(Trans)  ~ normal(0, 3);

  
  for (m in 1:M_pro) {
    int g = gene_id_pro[m];
    int c = celltype_id_pro[m];
    real se = 1 / sqrt(n_pro[m]);
    bar_pro[m] ~ normal(logP[g,c], se);
  }
  
  for (m in 1:M_mrna) {
    int g = gene_id_mrna[m];
    int c = celltype_id_mrna[m];
    real se = 1 / sqrt(n_mRNA[m]); //
    bar_mrna[m] ~ normal(logM[g,c], se);
  }
  
  for (m in 1:M_deg) {
    int g = gene_id_deg[m];
    int c = celltype_id_deg[m];
    real se = 1 / sqrt(n_deg[m]);
    bar_deg[m] ~ normal(logDeg[g,c], se);
  }
  

}

