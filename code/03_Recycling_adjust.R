#### Required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(QuantQC)
library(stringr)
library(foreach)
library(doParallel)
library(reshape2)

library(lme4)
library(splines)
library(gridExtra)
library(mgcv)




##################
# raw data path
##################

path_raw <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/raw/02_Recycle_searches/'

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


##################
# Functions for recycle correction analysis
##################
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

proc_recycle <- function(data){

  
  data <- data %>% filter(Ms1.Area != 0)

  data$kcount <- str_count(data$Stripped.Sequence, "K")
  data <- data %>% filter(kcount == 2)
  data$hcount <- str_count(data$Precursor.Id,"Label")

  data$seqcharge <- paste0(data$Stripped.Sequence,data$Precursor.Charge)

  data$plex <- substr(data$Precursor.Id[1:nrow(data)], 10, 10)
  data$ID <- paste0(data$Run,data$plex)

  HL_dat <- data %>% filter(hcount == 1)
  HH_dat <- data %>% filter(hcount == 2)
  sect <- intersect(HL_dat$seqcharge,HH_dat$seqcharge)
  HL_dat <- HL_dat %>% filter(seqcharge %in% sect)
  HH_dat <- HH_dat %>% filter(seqcharge %in% sect)

  HL_dat <- HL_dat %>% group_by(seqcharge,ID) %>% summarise(Ms1.Area = median(Ms1.Area,na.rm=T))

  HL_dat <- reshape2::dcast(HL_dat,seqcharge ~ ID, value.var = 'Ms1.Area')
  HH_dat <- reshape2::dcast(HH_dat,seqcharge ~ ID, value.var = 'Ms1.Area')
  sect2 <- intersect(colnames(HL_dat),colnames(HH_dat))

  HL_dat <- HL_dat[,sect2]
  HH_dat <- HH_dat[,sect2]
  rownames(HL_dat) <- HL_dat$seqcharge
  rownames(HH_dat) <- HH_dat$seqcharge
  HL_dat$seqcharge <- NULL
  HH_dat$seqcharge <- NULL
  HL_dat <- as.matrix(HL_dat)
  HH_dat <- as.matrix(HH_dat)

  #pct_Heavy <- (HL_dat/HH_dat)/(HL_dat/HH_dat+2) # Marko
  pct_Heavy <- 2/(HL_dat/HH_dat+2)
  #pct_Heavy <- HL_dat/HH_dat

  return(pct_Heavy)
}

meta_data_proc <- function(meta,data,prep){

  meta$label[meta$label == "4"] <- "8"
  meta$ID_mod <- paste0(meta$Run,meta$label)
  meta <- meta %>% filter(is.na(Run) == F)
  rownames(meta) <- meta$ID_mod
  sect <- intersect(colnames(data),rownames(meta))

  meta <- meta[sect,]
  data <- data[,sect]
  colnames(data) <- meta$ID

  data <- colMeans(data,na.rm = T)
  data <- as.data.frame(data)
  data$ID <- paste0(rownames(data),prep)

  return(data)
}

fit_sum_exp_model <- function(recycle_df){


  # IC for guess
  a_init <- log(2)/0.7
  b_init <- log(2)/20


  # Save parameters
  age_save <- c()
  tissue_save <- c()
  a_save <- c()
  b_save <- c()



  for(i in unique(recycle_df$age)){
    print(i)
    for (j in unique(recycle_df$Recycle_type)){
      print(j)

      recycle_df_temp <- recycle_df %>% filter(age == i)
      recycle_df_temp <- recycle_df_temp %>% filter(Recycle_type == j)

      if(nrow(recycle_df_temp) >0){

        time_obs <- c(0, recycle_df_temp$Time)
        recycle_value_obs <- c(0,recycle_df_temp$avg)


        # Use nls() to fit the model and tune a and b
        fit <- nls(recycle_value_obs ~ 1 - .5 * exp(-a * time_obs) - .5 * exp(-b * time_obs),
                   start = list(a = a_init, b = b_init),
                   control = nls.control(maxiter = 1000),
                   algorithm = "port")



        # Extract the fitted parameters
        params <- coef(fit)
        a_fitted <- as.numeric(params["a"])
        b_fitted <- as.numeric(params["b"])

        age_save <- c(age_save,i)
        tissue_save <- c(tissue_save,j)
        a_save <- c(a_save,a_fitted)
        b_save <- c(b_save,b_fitted)

      }
    }
  }


  recycle_df_final <- as.data.frame(cbind(age_save,tissue_save,a_save,b_save))
  recycle_df_final$a_save <- as.numeric(recycle_df_final$a_save)
  recycle_df_final$b_save <- as.numeric(recycle_df_final$b_save)

  return(recycle_df_final)
}

NA_false_ratios = function(ratio_mat, t, anno, recycle_params,prep){

  recycle_params <- recycle_params %>% filter(Time == t)

  colnames(ratio_mat) <- paste0(colnames(ratio_mat),prep)

  ratio_mat <- ratio_mat[,intersect(colnames(ratio_mat),meta_data$ID)]

  for(i in 1:ncol(ratio_mat)){

    meta_temp <- meta_data %>% filter(ID == colnames(ratio_mat)[i])
    Final_values_temp <- recycle_params %>% filter(age == meta_temp$age)
    if(meta_temp$Cell_Type %in% c('Basal','Secratory','Cilliated','Immune')){
      Final_values_temp <- Final_values_temp %>% filter(tissue_save =='Epithelial')
    }else{
      Final_values_temp <- Final_values_temp %>% filter(tissue_save =='NotEpithelial')
    }

    for (j in 1:nrow(ratio_mat)) {

      # Get the current L and L0 values
      L_measured <- (ratio_mat[j, i])

      if(is.na(L_measured) == F){
        if(L_measured < Final_values_temp$avg){
          ratio_mat[j, i] <- NA
        }


      }
    }

  }
  return(ratio_mat)
}

# Function that solves for Kr (alpha is variable used) math is in methods or can read from function
L_function <- function(alpha, L0, a, b, t) {
  c = .5
  beta = L0*alpha
  #Lt <- (exp(-(a*t) - b*t)*(-(a*beta*c*exp(a*t)) - b*beta*c*exp(b*t) +
  #                            a*beta*c*exp(a*t + b*t) + b*beta*c*exp(a*t + b*t) + a*b*exp(a*t + b*t)*L0))/(a*b)


  Lt <- (a*beta*c - 2*alpha*beta*c + b*beta*c + alpha*beta*c*exp((-a + alpha)*t) - b*beta*c*exp((-a + alpha)*t) -
           a*beta*c*exp((alpha - b)*t) + alpha*beta*c*exp((alpha - b)*t) - a*alpha*L0 + alpha^2*L0 + a*b*L0 - alpha*b*L0)/
    ((-a + alpha)*(alpha - b)*exp(alpha*t))

  return(Lt)
}

# RSS objective function for optim to solve for Kr
objective_function <- function(alpha, L_measured, L0, a, b, t) {

  L_function <- function(alpha, L0, a, b, t) {
    c = .5
    beta = L0*alpha
    #Lt <- (exp(-(a*t) - b*t)*(-(a*beta*c*exp(a*t)) - b*beta*c*exp(b*t) +
    #                            a*beta*c*exp(a*t + b*t) + b*beta*c*exp(a*t + b*t) + a*b*exp(a*t + b*t)*L0))/(a*b)


    Lt <- (a*beta*c - 2*alpha*beta*c + b*beta*c + alpha*beta*c*exp((-a + alpha)*t) - b*beta*c*exp((-a + alpha)*t) -
             a*beta*c*exp((alpha - b)*t) + alpha*beta*c*exp((alpha - b)*t) - a*alpha*L0 + alpha^2*L0 + a*b*L0 - alpha*b*L0)/
      ((-a + alpha)*(alpha - b)*exp(alpha*t))

    return(Lt)
  }

  L_predicted <- L_function(alpha, L0, a, b, t)
  residuals <- L_predicted - (L_measured)
  return(sum(residuals^2))
}

NA_false_ratios = function(ratio_mat, t, anno, recycle_params,prep){

  recycle_params <- recycle_params %>% filter(Time == t)

  colnames(ratio_mat) <- paste0(colnames(ratio_mat),prep)

  ratio_mat <- ratio_mat[,intersect(colnames(ratio_mat),meta_data$ID)]

  for(i in 1:ncol(ratio_mat)){

    meta_temp <- anno %>% filter(ID == colnames(ratio_mat)[i])
    Final_values_temp <- recycle_params %>% filter(age == meta_temp$age)

    if(meta_temp$Cell_Type %in% c('Basal','Secratory','Cilliated','Immune')){
      Final_values_temp <- Final_values_temp %>% filter(Recycle_type =='Epithelial')
    }else{
      Final_values_temp <- Final_values_temp %>% filter(Recycle_type =='NotEpithelial')
    }

    for (j in 1:nrow(ratio_mat)) {

      # Get the current L and L0 values
      L_measured <- (ratio_mat[j, i])

      if(is.na(L_measured) == F){

        if((1-L_measured) > Final_values_temp$avg){
          ratio_mat[j, i] <- NA
        }


      }
    }

  }
  return(ratio_mat)
}

Recycle_adjust_par_bar <- function(L_mat, L0_mat, t, anno, recycle_params, prep, meta_data, objective_function, ncores = 2) {
  # Load required packages
  library(dplyr)
  library(foreach)
  library(doSNOW)

  # Adjust column names by appending 'prep'
  colnames(L_mat)  <- paste0(colnames(L_mat), prep)
  colnames(L0_mat) <- paste0(colnames(L0_mat), prep)

  # Subset to the common columns in meta_data
  common_IDs <- intersect(colnames(L_mat), meta_data$ID)
  L_mat  <- L_mat[, common_IDs, drop = FALSE]
  L0_mat <- L0_mat[, common_IDs, drop = FALSE]

  # Prepare an output matrix with appropriate dimnames
  out_mat <- matrix(NA, nrow = nrow(L_mat), ncol = ncol(L_mat),
                    dimnames = list(rownames(L_mat), colnames(L_mat)))

  # Set up a SNOW cluster and register it
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)

  # Set up a progress bar that tracks progress over columns
  pb <- txtProgressBar(max = ncol(L_mat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Parallelize the outer loop over columns using foreach:
  # Each iteration returns a vector of optimized alpha values for one column.
  result_list <- foreach(i = seq_len(ncol(L_mat)),
                         .packages = "dplyr",
                         .options.snow = opts) %dopar% {
                           col_res <- rep(NA, nrow(L_mat))  # initialize result vector for column i
                           current_ID <- colnames(L_mat)[i]

                           # Get metadata for the current column/sample
                           meta_temp <- anno %>% filter(ID == current_ID)

                           # Filter recycle_params based on age and tissue type.
                           if (meta_temp$Cell_Type %in% c('Basal', 'Secratory', 'Cilliated', 'Immune')) {
                             Final_values_temp <- recycle_params %>%
                               filter(age_save == meta_temp$age, tissue_save == 'Epithelial')
                           } else {
                             Final_values_temp <- recycle_params %>%
                               filter(age_save == meta_temp$age, tissue_save == 'NotEpithelial')
                           }

                           # Process each row (peptide) in the current column
                           for (j in seq_len(nrow(L_mat))) {
                             L_measured <- L_mat[j, i]
                             L0_value   <- L0_mat[j, i]
                             if (!is.na(L_measured) && !is.na(L0_value)) {
                               # Run the 1D optimization (using L-BFGS-B here)
                               result <- optim(
                                 par = 0.3,
                                 fn = objective_function,
                                 L_measured = L_measured,
                                 L0 = L0_value,
                                 a = Final_values_temp$a_save,
                                 b = Final_values_temp$b_save,
                                 t = t,
                                 method = "L-BFGS-B",
                                 lower = 0,
                                 upper = 10
                               )
                               col_res[j] <- result$par
                             }
                           }
                           # Return the vector for this column
                           col_res
                         }

  # Close the progress bar and stop the cluster
  close(pb)
  stopCluster(cl)

  # Combine the results from each column into a matrix.
  # result_list is a list of column vectors; cbind them together.
  out_mat <- do.call(cbind, result_list)
  colnames(out_mat) <- colnames(L_mat)
  rownames(out_mat) <- rownames(L_mat)

  return(out_mat)
}


##################
# Read in raw data
##################


columns_to_read <-c('Genes','Run','Lib.PG.Q.Value','RT','Precursor.Id','Stripped.Sequence','Precursor.Mz',
                    'Precursor.Charge','Precursor.Quantity','Ms1.Area','Protein.Group','Translated.Q.Value','Channel.Q.Value')

male5day_1 <- paste0(path_raw,'Rep1_5day_male/Report_plate1.tsv')
male5day_2 <- paste0(path_raw,'Rep1_5day_male/Report_plate2.tsv')
male5day_1 <- data.table::fread(male5day_1,select = columns_to_read)
male5day_2 <- data.table::fread(male5day_2,select = columns_to_read)
male5day <- rbind(male5day_1,male5day_2)
rep1_5day <- proc_recycle(male5day)
load(paste0(path_dat,'03_QuantQC_objects/r1_5day_male.RData'))
meta_data_mini <- r1_5day_male@meta.data %>% dplyr::select('ID','Order','label','Run')
rep1_5day <- meta_data_proc(meta_data_mini,rep1_5day,'_prep1')
rep1_5day$Time <- 5

female5day_1 <- paste0(path_raw,'Rep2_5day_female/Report_plate1.tsv')
female5day_2 <- paste0(path_raw,'Rep2_5day_female/Report_plate2.tsv')
female5day_1 <- data.table::fread(female5day_1,select = columns_to_read)
female5day_2 <- data.table::fread(female5day_2,select = columns_to_read)
female5day <- rbind(female5day_1,female5day_2)
rep2_5day <- proc_recycle(female5day)
load(paste0(path_dat,'03_QuantQC_objects/r2_5day_female.RData'))
meta_data_mini <- r2_5day_female@meta.data %>% dplyr::select('ID','Order','label','Run')
rep2_5day <- meta_data_proc(meta_data_mini,rep2_5day,'_prep2')
rep2_5day$Time <- 5

male10day <- paste0(path_raw,'Rep3_10day_male/Report.tsv')
male10day <- data.table::fread(male10day,select = columns_to_read)
rep1_10day <- proc_recycle(male10day)
load(paste0(path_dat,'03_QuantQC_objects/r3_10day_male.RData'))
meta_data_mini <- r3_10day_male@meta.data %>% dplyr::select('ID','Order','label','Run')
rep1_10day <- meta_data_proc(meta_data_mini,rep1_10day,'_prep3')
rep1_10day$Time <- 10

female10day <- paste0(path_raw,'Rep4_10day_female/Report.tsv')
female10day <- data.table::fread(female10day,select = columns_to_read)
rep2_10day <- proc_recycle(female10day)
load(paste0(path_dat,'03_QuantQC_objects/r4_10day_female.RData'))
meta_data_mini <- r4_10day_female@meta.data %>% dplyr::select('ID','Order','label','Run')
rep2_10day <- meta_data_proc(meta_data_mini,rep2_10day,'_prep4')
rep2_10day$Time <- 10


all_data <- rbind(rep1_5day,rep2_5day,rep1_10day,rep2_10day)
meta_data <- read.csv(paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/sc_protein_annotations.csv'))
meta_data$Order <- NULL
meta_data$label <- NULL
meta_data$Run <- NULL

all_data <- all_data %>% left_join(meta_data, by = c('ID'))
all_data <- all_data %>% filter(is.na(age)==F)
all_data <- all_data %>% filter(Cell_Type!='Shwann')


##################
# Solve for cell type specific free light AA pool
##################

count_df <- all_data %>%
  dplyr::group_by(Time, Cell_Type, age) %>%
  dplyr::summarize(
    n = n(),                      # number of data points
    y_max = max(data, na.rm = TRUE)  # max y-value in this group
  )

# 2. Plot: boxplot + facet + text for counts
ggplot(all_data, aes(x = Cell_Type, y = data, color = age)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~Time) + ylab('% heavy AA')+
  # Add geom_text for the counts
  geom_text(
    data = count_df,
    aes(
      x = Cell_Type,
      y = y_max + 0.1 * y_max,  # place text 10% above the max
      label = n,
      color = age
    ),
    position = position_dodge(width = 0.8),
    vjust = 0
  )

all_data$Recycle_type <- NA
all_data$Recycle_type[all_data$Cell_Type %in% c('Basal','Secratory','Cilliated','Immune')] <- 'Epithelial'
all_data$Recycle_type[!all_data$Cell_Type %in% c('Basal','Secratory','Cilliated','Immune')] <- 'NotEpithelial'

Final_values <- all_data %>%
  dplyr::group_by(Time, age,Recycle_type) %>%
  dplyr::summarize(
    avg = median(data, na.rm = TRUE)  # max y-value in this group
  )

Final_values[(nrow(Final_values)+1):(nrow(Final_values)+2),] <- NA
Final_values$Time[7:8] <- 10
Final_values$age[7:8] <- 'Old'
Final_values$Recycle_type[7:8] <- c("Epithelial","NotEpithelial")
Final_values$avg[7:8] <- c(0.6797455-.03,0.6600058-.03)



Final_values$avg[Final_values$Time==10] <- Final_values$avg[Final_values$Time==10] + .1
Final_values$avg[Final_values$Time==5] <- Final_values$avg[Final_values$Time==5]  + .06


Fitted_curves <- fit_sum_exp_model(Final_values)


t <- seq(0,20,by = .1)
n = 3
L_aa_fitted <-  .5 * exp(-Fitted_curves$a_save[n] * t) + .5 * exp(-Fitted_curves$b_save[n] * t)
df_AA_fit <- data.frame(time = t,H_aa_free = 1- L_aa_fitted)


ggplot(df_AA_fit,aes(x = t,y = H_aa_free)) + geom_line() +
  theme_classic(base_size = 15) + xlab('Days') +ylab('% Heavy amino acid in free pool')+
  ylim(c(0,1)) + xlim(c(0,20))



##################
# Apply recycling corrections to single cells
##################


r1 <- r1_5day_male@miceotopes@Raw_L/(r1_5day_male@miceotopes@Raw_L+r1_5day_male@miceotopes@Raw_H)
r2 <- r2_5day_female@miceotopes@Raw_L/(r2_5day_female@miceotopes@Raw_L+r2_5day_female@miceotopes@Raw_H)
r3 <- r3_10day_male@miceotopes@Raw_L/(r3_10day_male@miceotopes@Raw_L+r3_10day_male@miceotopes@Raw_H)
r4 <- r4_10day_female@miceotopes@Raw_L/(r4_10day_female@miceotopes@Raw_L+r4_10day_female@miceotopes@Raw_H)



na_r1 <- NA_false_ratios(r1,5,meta_data,Final_values,'_prep1')
na_r2 <- NA_false_ratios(r2,5,meta_data,Final_values,'_prep2')
na_r3 <- NA_false_ratios(r3,10,meta_data,Final_values,'_prep3')
na_r4 <- NA_false_ratios(r4,10,meta_data,Final_values,'_prep4')


colnames(r1)<-paste0(colnames(r1),'_prep1')
sum(is.na(na_r1)==F)/sum(is.na(r1[,colnames(na_r1)])==F)

colnames(r2)<-paste0(colnames(r2),'_prep2')
sum(is.na(na_r2)==F)/sum(is.na(r2[,colnames(na_r2)])==F)

colnames(r3)<-paste0(colnames(r3),'_prep3')
sum(is.na(na_r3)==F)/sum(is.na(r3[,colnames(na_r3)])==F)

colnames(r4)<-paste0(colnames(r4),'_prep4')
sum(is.na(na_r4)==F)/sum(is.na(r4[,colnames(na_r4)])==F)


adj_r1 <- Recycle_adjust_par_bar(r1_5day_male@miceotopes@Raw_L,(r1_5day_male@miceotopes@Raw_L+r1_5day_male@miceotopes@Raw_H)
                       ,5,meta_data,Fitted_curves,'_prep1',meta_data = meta_data,objective_function)

adj_r1[adj_r1==10] <- NA

adj_r2 <- Recycle_adjust_par_bar(r2_5day_female@miceotopes@Raw_L,(r2_5day_female@miceotopes@Raw_L+r2_5day_female@miceotopes@Raw_H),
                         5,meta_data,Fitted_curves,'_prep2',meta_data = meta_data,objective_function)
adj_r2[adj_r2==10] <- NA

adj_r3 <- Recycle_adjust_par_bar(r3_10day_male@miceotopes@Raw_L,(r3_10day_male@miceotopes@Raw_L+r3_10day_male@miceotopes@Raw_H),
                         10,meta_data,Fitted_curves,'_prep3',meta_data = meta_data,objective_function)
adj_r3[adj_r3==10] <- NA

adj_r4 <- Recycle_adjust_par_bar(r4_10day_female@miceotopes@Raw_L,(r4_10day_female@miceotopes@Raw_L+r4_10day_female@miceotopes@Raw_H),
                         10,meta_data,Fitted_curves,'_prep4',meta_data = meta_data,objective_function)
adj_r4[adj_r4==10] <- NA


#### Looking at not corrected values for comparison
r1 <- -log(r1[,colnames(na_r1)])/5
r2 <- -log(r2[,colnames(na_r2)])/5
r3 <- -log(r3[,colnames(na_r3)])/10
r4 <- -log(r4[,colnames(na_r4)])/10

median(log(2)/r1,na.rm=T)
median(log(2)/r2,na.rm=T)
median(log(2)/r3,na.rm=T)
median(log(2)/r4,na.rm=T)

median((log(2)/adj_r1),na.rm=T)
median((log(2)/adj_r2),na.rm=T)
median((log(2)/adj_r3),na.rm=T)
median((log(2)/adj_r4),na.rm=T)




adj_r1_ <- adj_r1[rowSums(is.na(adj_r1)==F) > 20,]
adj_r2_ <- adj_r2[rowSums(is.na(adj_r2)==F) > 20,]
adj_r3_ <- adj_r3[rowSums(is.na(adj_r3)==F) > 20,]
adj_r4_ <- adj_r4[rowSums(is.na(adj_r4)==F) > 20,]

ss <- intersect(intersect(intersect(rownames(adj_r1_),rownames(adj_r2_)),rownames(adj_r3_)),rownames(adj_r4_))


adj_r1[adj_r1==0] <- NA
adj_r2[adj_r2==0] <- NA
adj_r3[adj_r3==0] <- NA
adj_r4[adj_r4==0] <- NA

young_p1 <- meta_data %>% filter(age == 'Young' & sample == '1' & Cell_Type == 'Basal')
young_p2 <- meta_data %>% filter(age == 'Young' & sample == '2' & Cell_Type == 'Basal')

young_p3 <- meta_data %>% filter(age == 'Young' & sample == '3' & Cell_Type == 'Basal')
young_p4 <- meta_data %>% filter(age == 'Young' & sample == '4' & Cell_Type == 'Basal')

old_p1 <- meta_data %>% filter(age == 'Old' & sample == '1')



#### Compare rates for 5 and 10 day mice with corrected data
df_corrected_compare <- data.frame(
  logx = rowMeans(log2(cbind(adj_r2[ss,young_p2$ID],adj_r1[ss,young_p1$ID])), na.rm = TRUE),
  logy = rowMeans(log2(cbind(adj_r3[ss,young_p3$ID],adj_r4[ss,young_p4$ID])), na.rm = TRUE)
)

cor(df_corrected_compare$logx,df_corrected_compare$logy,use = 'pairwise.complete.obs')

df_corrected_compare$X <- 2 ^ df_corrected_compare$logx                
df_corrected_compare$Y <- 2 ^ df_corrected_compare$logy                

ggplot(df_corrected_compare, aes(X, Y)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  geom_abline(intercept = -.07, slope = 1/1.1,color = 'red')+
  scale_x_log10()+scale_y_log10()+
  coord_cartesian(xlim = c(.05,1),ylim = c(.05,1))+
  labs(x = "5 day degradation rate, Days^-1",
       y = "10 day degradation rate, Days^-1") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1), 
    axis.line     = element_blank()                                       
  )

1/TLS(df$logx,df$logy)[[1]]


#### Compare rates for 5 and 10 day mice without correction applied for reference
df_uncor <- data.frame(
  logx = rowMeans(log2(cbind(r1[ss,young_p1$ID],r2[ss,young_p2$ID])), na.rm = TRUE),
  logy = rowMeans(log2(cbind(r3[ss,young_p3$ID],r4[ss,young_p4$ID])), na.rm = TRUE)
)

cor(df_uncor$logx,df_uncor$logy,use = 'pairwise.complete.obs')

df_uncor$X <- 2 ^ df_uncor$logx                
df_uncor$Y <- 2 ^ df_uncor$logy                

ggplot(df_uncor, aes(X, Y)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  geom_abline(intercept = -.38, slope = 1/1.6, linetype = "solid",color = 'red')+
  scale_x_log10()+scale_y_log10()+
  coord_cartesian(xlim = c(.03,.2),ylim = c(.03,.2))+
  labs(x = "5 day degradation rate, Days^-1",
       y = "10 day degradation rate, Days^-1") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # add frame
    axis.line     = element_blank()                                       # hide dup. axes
  )

1/TLS(df$logx,df$logy)[[1]]


#### Compare magnitudes of corrected and uncorrected half life

df_compare <- data.frame(
  uncorrected = rowMeans(r1[ss,], na.rm = TRUE),
  corrected = rowMeans(adj_r1[ss,], na.rm = TRUE)
)
df_compare <- melt(df_compare)

ggplot(df_compare, aes(y = log(2)/value, fill = variable)) +
  geom_histogram(alpha = 0.5, bins = 30, colour = "black",
                 orientation = "y") +           # <- horizontal bars
  facet_wrap(~ variable, ncol = 2) +
  labs(
    x = "Count in each bin",
    y = "Value",
    title = ""
  ) +
  dot_plot +
  theme(legend.position = "none")+
  scale_y_log10() + ylab('Clearance half life (days)')+ xlab('# of proteins')+
  scale_fill_manual(values = c('grey80','black'))


#### Compare old and young clearance rates

df_old_vs_young <- data.frame(
  logx = rowMeans((cbind(adj_r1[ss,young_p1$ID])), na.rm = TRUE),
  logy = rowMeans((cbind(adj_r1[ss,old_p1$ID])), na.rm = TRUE)
)

cor(df_old_vs_young$logx,df_old_vs_young$logy,use = 'pairwise.complete.obs')

ggplot(df_old_vs_young, aes(logx, logy)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_abline(intercept = 0,  slope = 1,  linetype = "dashed") +
  geom_abline(intercept = -.1, slope = 1,  colour = "red") +
  scale_x_log10() + 
  scale_y_log10() +
  coord_cartesian(xlim = c(.05, 1), ylim = c(.05, 1)) +
  labs(x = "2 months", y = "24 months") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # add frame
    axis.line     = element_blank()                                       # hide dup. axes
  )



#### Save data in QQC objects 

colnames(adj_r1) <- str_remove(colnames(adj_r1), "_prep1")
colnames(adj_r2) <- str_remove(colnames(adj_r2), "_prep2")
colnames(adj_r3) <- str_remove(colnames(adj_r3), "_prep3")
colnames(adj_r4) <- str_remove(colnames(adj_r4), "_prep4")

r1_5day_male@miceotopes@Alpha_pep <- adj_r1
r2_5day_female@miceotopes@Alpha_pep <- adj_r2
r3_10day_male@miceotopes@Alpha_pep <- adj_r3
r4_10day_female@miceotopes@Alpha_pep <- adj_r4


r1_5day_male <- Miceotope_protein_collapse(r1_5day_male)
r2_5day_female <- Miceotope_protein_collapse(r2_5day_female)
r3_10day_male <- Miceotope_protein_collapse(r3_10day_male)
r4_10day_female <- Miceotope_protein_collapse(r4_10day_female)


##################
# Save updated QQC objects
##################

save(r1_5day_male, file = paste0(path_dat,'03_QuantQC_objects/r1_5day_male.RData'))
save(r2_5day_female, file = paste0(path_dat,'03_QuantQC_objects/r2_5day_female.RData'))
save(r3_10day_male, file = paste0(path_dat,'03_QuantQC_objects/r3_10day_male.RData'))
save(r4_10day_female, file = paste0(path_dat,'03_QuantQC_objects/r4_10day_female.RData'))

##################
# Save Single cell matricies for deg
##################

p1_alpha <-  QuantQC::Normalize_reference_vector(r1_5day_male@miceotopes@Alpha_prot, log = T)
p2_alpha <-  QuantQC::Normalize_reference_vector(r2_5day_female@miceotopes@Alpha_prot , log = T)
p3_alpha <-  QuantQC::Normalize_reference_vector(r3_10day_male@miceotopes@Alpha_prot , log = T)# QuantQC::Normalize_reference_vector
p4_alpha <-  QuantQC::Normalize_reference_vector(r4_10day_female@miceotopes@Alpha_prot , log = T)# QuantQC::Normalize_reference_vector

p1_alpha_abs <-  r1_5day_male@miceotopes@Alpha_prot
p2_alpha_abs <-  r2_5day_female@miceotopes@Alpha_prot
p3_alpha_abs <-  r3_10day_male@miceotopes@Alpha_prot
p4_alpha_abs <-  r4_10day_female@miceotopes@Alpha_prot

### Save deg rates for all IDed proteins
colnames(p1_alpha) <- paste0(colnames(p1_alpha),'_prep1')
colnames(p2_alpha) <- paste0(colnames(p2_alpha),'_prep2')
colnames(p3_alpha) <- paste0(colnames(p3_alpha),'_prep3')
colnames(p4_alpha) <- paste0(colnames(p4_alpha),'_prep4')
p1_alpha <- as.data.frame(p1_alpha)
p1_alpha$prot <- rownames(p1_alpha)
p2_alpha <- as.data.frame(p2_alpha)
p2_alpha$prot <- rownames(p2_alpha)
p3_alpha <- as.data.frame(p3_alpha)
p3_alpha$prot <- rownames(p3_alpha)
p4_alpha <- as.data.frame(p4_alpha)
p4_alpha$prot <- rownames(p4_alpha)
p_1and2 <- p1_alpha %>% merge(p2_alpha, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_alpha, by = 'prot',all = TRUE)
p_all_alpha <- p_123 %>% merge(p4_alpha, by = 'prot',all = TRUE)
rownames(p_all_alpha) <- p_all_alpha$prot
p_all_alpha$prot <- NULL
p_all_alpha <- as.matrix(p_all_alpha)

write.csv(p_all_alpha,paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/clearance_relative.csv'))


colnames(p1_alpha_abs) <- paste0(colnames(p1_alpha_abs),'_prep1')
colnames(p2_alpha_abs) <- paste0(colnames(p2_alpha_abs),'_prep2')
colnames(p3_alpha_abs) <- paste0(colnames(p3_alpha_abs),'_prep3')
colnames(p4_alpha_abs) <- paste0(colnames(p4_alpha_abs),'_prep4')
p1_alpha_abs <- as.data.frame(p1_alpha_abs)
p1_alpha_abs$prot <- rownames(p1_alpha_abs)
p2_alpha_abs <- as.data.frame(p2_alpha_abs)
p2_alpha_abs$prot <- rownames(p2_alpha_abs)
p3_alpha_abs <- as.data.frame(p3_alpha_abs)
p3_alpha_abs$prot <- rownames(p3_alpha_abs)
p4_alpha_abs <- as.data.frame(p4_alpha_abs)
p4_alpha_abs$prot <- rownames(p4_alpha_abs)
p_1and2 <- p1_alpha_abs %>% merge(p2_alpha_abs, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_alpha_abs, by = 'prot',all = TRUE)
p_all_alpha_abs <- p_123 %>% merge(p4_alpha_abs, by = 'prot',all = TRUE)
rownames(p_all_alpha_abs) <- p_all_alpha_abs$prot
p_all_alpha_abs$prot <- NULL
p_all_alpha_abs <- as.matrix(p_all_alpha_abs)

write.csv(p_all_alpha_abs,paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/clearance_absolute.csv'))





