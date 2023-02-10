#####################################
# Sorghum BCNAM analytical pipeline #
#####################################

# library and ad-hoc functions

library(mppR)
library(dplyr)
library(data.table)
library(ggplot2)
library(igraph)

# set working directory
# location of sorghum_BCNAM_analysis folder
setwd('~/sorghum_BCNAM_analysis')

# load ad-hoc functions
source('./scripts/BCNAM_analysis_fct.R')

#### overview of the data and results folder structure ####

# ./software: .tar.gz file of the mppR package used for the analysis

# ./scripts: R script of the pipeline and of the ad-hoc functions

# ./data/mppData : mppData object containing processed genotypes and phenotype data
# The mppData objects are too voluminous to be stored on the github repository
# they can be download from ... and stored in the ./data/mppData folder.

# ./data/threshold: QTL detection threshold estimated using Li and Ji method

# ./data/climate : climatic data of the testing sites (Sotuba, Cinzana in 2012-13)
# and of the grid of points covering part of Mali.
# .data/EC: list of EC variables per testing site and point grid

# ./results/QTL_detection: list of detected QTLs per population x trait configuration
# ./results/pheno_x_EC: list of the 5 best covariates with window information
# .results/QTL_x_EC: list with the QTL_main vs QTL_x_E effect and the QTLxEC detail res
# for 5 best ECs.
# ./results/database: database with all relevant QTL results (position, effect, etc.)
# 

#### QTL detection ####

load('./data/threshold/d_thr.RData')
pop_id <- c('KK2012', 'KK2013', 'GR2012', 'GR2013', 'Lata')
n_pop <- length(pop_id)
tr_id <- c("FLAGDD", "PH", "NODE_N", "NODE_L", "PED_L", "PAN_L", "G1000_WGH", "G_YIELD")
n_trait <- length(tr_id)

for(p in 1:n_pop){
  
  load(file = file.path('./data/mppData', pop_id[p], 'mppData.RData'))
  thr_i <- d_thr[rownames(d_thr) == pop_id[p], 2]
  
  for(t in 1:n_trait){
    
    tr_pos <- grep(pattern = tr_id[t], x = colnames(mppData$pheno))
    
    # subset the mppData object: select the genotypes with at least
    # phenotypic information in one environment.
    
    pheno_i <- mppData$pheno[, tr_pos]
    g_sel_i <- apply(pheno_i, 1, function(x) any(!is.na(x)))
    mppData_i <- subset(mppData, gen.list = g_sel_i)
    env_empty <- apply(pheno_i, 2, function(x) all(is.na(x)))
    
    tr_s <- colnames(pheno_i)[!env_empty]
    
    res_fold <- file.path('./results/QTL_analysis', pop_id[p])
    
    t1 <- Sys.time()
    QTL <- mppGE_proc(trait.name = tr_id[t], mppData = mppData_i,
                      trait = tr_s, VCOV_data = "minus_cof", thre.cof = thr_i,
                      win.cof = 500, thre.QTL = thr_i, win.QTL = 30, verbose = FALSE,
                      output.loc = res_fold, n.cores = 5)
    t2 <- Sys.time()
    t_diff <- t2 - t1
    print(tr_id[t])
    print(t_diff)
    
  }
  
}

#### phenotype x EC analysis ####

tr_id <- c("FLAGDD", "PH", 'NODE_N', 'NODE_L', "PED_L", "PAN_L", "G_YIELD")
n_trait <- length(tr_id)
pop_id <- c('GR2012', 'GR2013', 'KK2012', 'KK2013')
n_pop <- length(pop_id)
year_id <- c('2012', '2013', '2012', '2013')
env_id <- c('SB1', 'SB2', 'CZ1', 'CZ2')

EC_id <- c('rain', 'hum', 'VPD', 'SVP', 'ETP', 'PETP', 'Tmin', 'Tmax',
           'Trange', 'DD', 'FRUE', 'hSun', 'photoperiod', 'solarRad', 'photothermal')

type <- c('sum', 'mean', 'mean', 'mean', 'mean', 'mean', 'mean', 'mean', 'mean',
          'sum', 'mean', 'sum', 'mean', 'sum', 'sum')

# EC list
load(file = './data/climate/EC_2012.RData'); EC_list_2012 <- EC_list
load(file = './data/climate/EC_2013.RData'); EC_list_2013 <- EC_list
rm(EC_list)
EC_list <- list(`2012` = EC_list_2012, `2013` = EC_list_2013)

# space to store the results
EC_best_list <- vector(mode = 'list', length = n_pop)
names(EC_best_list) <- pop_id

res_list <- vector(mode = 'list', length = n_pop)
names(res_list) <- pop_id

for(p in 1:n_pop){

  load(file = file.path('./data/mppData', pop_id[p], 'mppData.RData'))
  
  pheno <- data.frame(mppData$pheno)
  pheno <- pheno %>% select(contains('FLAG'))
  pheno <- pheno[, 1:4]
  FLAG_mean <- colMeans(x = pheno, na.rm = TRUE)
  
  crop_duration <- ceiling(mean(FLAG_mean)) + 40
  flag_time <- ceiling(mean(FLAG_mean)) + 10
  
  EC_list_p <- EC_list[year_id[p]][[1]]
  plot_dir <- file.path('results', 'Pheno_EC_plot', pop_id[p])
  
  EC_best_list[[p]] <- vector(mode = 'list', length = n_trait)
  names(EC_best_list[[p]]) <- tr_id
  
  # loop over trait
  for(t in 1:n_trait){
    
    pheno <- data.frame(mppData$pheno)
    pheno <- pheno %>% select(contains(tr_id[t]))
    trait_env_mean <- colMeans(x = pheno, na.rm = TRUE)
    
    if(tr_id[t] == 'FLAGDD'){ crop_dur_e <- flag_time
    } else {crop_dur_e <- crop_duration}
    
    EC_effect_pt <- EC_effect(trait_env_mean = trait_env_mean,
                              crop_duration = crop_dur_e,
                              EC_list = EC_list_p, type = type,
                              plot = TRUE, plot_dir = plot_dir,
                              p_title = tr_id[t], env_nm = env_id)
    
    EC_best_list[[p]][[t]] <- EC_effect_pt[-1, -4] %>% slice_max(R2_glb, n=5)
    
  }
}

# Save the results
save(EC_best_list, file = './results/Pheno_EC/EC_best_list.RData')

#### QTLxEC analysis ####

pop_id <- c('GR2012', 'GR2013', 'KK2012', 'KK2013', 'Lata')
n_pop <- length(pop_id)

tr_list <- c('FLAGDD', 'PH', 'NODE_N', 'NODE_L', 'PED_L', 'PAN_L', 'G1000_WGH', 'G_YIELD')
tr_list2 <- c('FLAGDD', 'PH', 'NODE_N', 'NODE_L', 'PED_L', 'PAN_L', 'G100_WGH', 'G_YIELD')
n_trait <- length(tr_list)

env_id <- c('SB1', 'SB2', 'CZ1', 'CZ2')
env_id2 <- c('LP', 'HP', 'KOL')

# QTL_file_loc <- 'C:/Users/vince/OneDrive/Documents/WD/ICRISAT/BCNAM/Results/QTL_analysis'
# res_folder <- 'C:/Users/vince/OneDrive/Documents/WD/ICRISAT/BCNAM/Results/QTLxEC'

# List of environment covariates
load(file = './results/Pheno_EC/EC_best_list.RData')

# space to store the results
res_list <- vector(mode = 'list', length = n_pop)
names(res_list) <- pop_id
for(i in 1:n_pop){
  res_list[[i]] <- vector(mode = 'list', length = n_trait)
  names(res_list[[i]]) <- tr_list
}

for(p in 1:n_pop){
  
  load(file = file.path('./data/mppData', pop_id[p], 'mppData.RData'))
  
  # loop over trait
  # 1:n_trait
  for(t in 1:n_trait){
    
    file_t <- file.path('./results/QTL_analysis', pop_id[p],
                        paste0('QTLGE_MPP_', tr_list[t], '_UN'), 'QTLs.RData')
    
    if(file.exists(file_t)){
      
      load(file = file_t)
      
      if(pop_id[p] == 'Lata'){
        
        if(tr_list[t] %in% c('PED_L', 'PAN_L')){ env_pt <- env_id2[1:2] } else { env_pt <- env_id2 }
        
        trait <- paste0(env_pt, '_', tr_list2[t])
        
        Qeff <- tryCatch(QTL_effect_main_QxE(mppData = mppData,  trait = trait, QTL = QTL,
                                             env_id = env_pt), error = function(x) NULL)
        
      } else {
        
        if(tr_list[t] == 'G1000_WGH'){
          
          trait <- paste0(env_id[1:2], '_', tr_list[t])
          
          Qeff <- tryCatch(QTL_effect_main_QxE(mppData = mppData,  trait = trait, QTL = QTL,
                                               env_id = env_id[1:2]), error = function(x) NULL)
          
        } else {
          
          trait <- paste0(env_id, '_', tr_list[t])
          
          Qmain_QxE <- tryCatch(QTL_effect_main_QxE(mppData = mppData,  trait = trait, QTL = QTL,
                                                    env_id = env_id), error = function(x) NULL)
          
          # form EC matrix
          EC_mat <- EC_best_list[[pop_id[p]]][[tr_list[t]]][, 5:8]
          EC_vect <- rownames(EC_mat)
          n_EC <- length(EC_vect)
          
          Qeff <- vector(mode = 'list', length = n_EC)
          names(Qeff) <- EC_vect
          EC_mat <- as.matrix(t(EC_mat))
          
          for(i in 1:n_EC){
            
            EC_mat_i <- EC_mat[, i, drop = FALSE]
            
            system.time(
              Qeff[[i]] <- tryCatch(QTL_effect_QxEC(mppData = mppData,  trait = trait,
                                                    QTL = QTL, env_id = env_id,
                                                    EC = EC_mat_i, thre_QTL = 1.301,
                                                    Qmain_QxE = Qmain_QxE),
                                    error = function(x) NULL)
            )
            
          }
          
        }
        
      }
      
      
      if(!is.null(Qeff)){
        
        res_list[[p]][[t]] <- Qeff
        print(paste(pop_id[p], tr_list[t]))
        print(Sys.time())
        
      } else {
        
        print(paste(pop_id[p], tr_list[t]))
        print('Failed')
        
      }
      
    }
    
  }
  
}

# Save the results
QTLxEC_res <- res_list
save(QTLxEC_res, file = './results/QTLxEC/QTLxEC_res.RData')

#### QTL database: QTL positions ####

pop_id <- c('GR2012', 'GR2013','KK2012', 'KK2013', 'Lata')
n_pop <- length(pop_id)
RP <- c('Grinkan', 'Grinkan', 'Kenin-Keni', 'Kenin-Keni', 'Lata')
year <- c('2012', '2013', '2012', '2013', '2013')

tr_list <- c('FLAGDD', 'PH', "NODE_N", "NODE_L", 'PAN_L', 'PED_L', 'G1000_WGH', 'G_YIELD')
n_trait <- length(tr_list)
c_nm <- c('RP', 'year', 'pop', 'trait', 'N_QTL', 'R2_glb', 'QTL', 'marker', 'chr',
          'bp', 'cM', 'log10pval', 'R2')

DB <- c()

for(p in 1:n_pop){
  
  RP_p <- RP[p]
  year_p <- year[p]
  pop_p <- pop_id[p]
  
  for(t in 1:n_trait){
    
    tr <- tr_list[t]
    
    file_pt <- file.path('./results/QTL_analysis', pop_p, paste0('QTLGE_MPP_', tr, '_UN'))
    
    file_t <- file.path(file_pt, 'QTLs.RData')
    G_res <- file.path(file_pt, 'Glb_res.RData')
    R2_res <- file.path(file_pt, 'QTL_R2.RData')
    
    if(file.exists(file_t)){
      
      load(file = file_t)
      load(file = G_res)
      load(file = R2_res)
      
      Q_id <- paste0('Q', 1:nrow(QTL))
      bp_pos <- as.numeric(sapply(strsplit(x = QTL$mk.names, split = '_'), `[[`, 2))
      
      DB_pt <- data.frame(RP_p, year_p, pop_p, tr, glb_res$N_QTL, glb_res$adj_R2,
                          Q_id, QTL$mk.names, QTL$chr, bp_pos, QTL$pos.cM,
                          QTL$log10pval, QTL.R2$adj.R2.diff)
      
      colnames(DB_pt) <- c_nm
      DB <- rbind(DB, DB_pt)
      
      
    }
    
  }
  
}

#### QTL database: aggregate QTLs at same position ####

DB <- data.table(DB)

cM_min <- 5

# add a unique indicator to the QTL
DB$Q_id <- paste0(DB$pop, '_', DB$trait, '_', DB$QTL)

tr_list <- c('FLAGDD', 'PH', "NODE_N", "NODE_L", 'PAN_L', 'PED_L', 'G1000_WGH', 'G_YIELD')
tr_sh_id <- c('FLAG', 'PH', "NODE_N", "NODE_L", 'PAN', 'PED', 'GWGH', 'YIELD')
n_trait <- length(tr_list)

QTL_un_vect <- c()
QTL_un_id <- c()

for(i in 1:n_trait){
  
  chr_list <- DB[trait == tr_list[i]][, .(chr_list = sort(unique(chr)))]
  chr_list <- unlist(chr_list[, 1])
  n_chr <- length(chr_list)
  
  for(c in 1:n_chr){
    
    DB_tr_chr <- DB[trait == tr_list[i] &  chr == chr_list[c]]
    
    pos_cM <- DB_tr_chr$cM
    n_pos <- length(pos_cM)
    pos_vect <- paste0('p_', round(pos_cM, 2))
    pos_num <- 1:n_pos
    
    if(n_pos > 1){
      
      r_vect <- rep(1:(n_pos-1), times = (n_pos-1):1)
      c_vect <- c()
      for(v in 2:n_pos){c_vect <- c(c_vect, v:n_pos)}
      
      vertex_pl <- c()
      vertex <- c()
      
      for(p in 1:n_pos){
        vertex_pl <- rbind(vertex_pl, rep(pos_vect[p], 2))
        vertex <- rbind(vertex, rep(pos_num[p], 2))
      }
      
      for(n in 1:length(r_vect)){
        
        if(abs(pos_cM[r_vect[n]] - pos_cM[c_vect[n]]) < cM_min){
          vertex_pl <- rbind(vertex_pl, c(pos_vect[c_vect[n]], pos_vect[r_vect[n]]))
          vertex_pl <- rbind(vertex_pl, c(pos_vect[r_vect[n]], pos_vect[c_vect[n]]))
          vertex <- rbind(vertex, c(pos_num[c_vect[n]], pos_num[r_vect[n]]))
          vertex <- rbind(vertex, c(pos_num[r_vect[n]], pos_num[c_vect[n]]))
        }
      }
      
      vertex <- apply(X = vertex, MARGIN = 1, FUN = function(x) x)
      vertex_pl <- apply(X = vertex_pl, MARGIN = 1, FUN = function(x) x)
      g_QTL <- graph(c(vertex))
      
      clu <- components(graph = g_QTL, mode = 'weak')
      grp <- igraph::groups(clu)
      n_qtl <- sapply(grp, length)
      
      Q_id <- rep(paste0('QTL_', tr_sh_id[i], '_', chr_list[c], '_'), n_pos)
      
      for(u in 1:length(n_qtl)){
        Q_id[grp[[u]]] <- paste0(Q_id[grp[[u]]], u)
      }
      
      QTL_un_vect <- c(QTL_un_vect, Q_id)
      QTL_un_id <- c(QTL_un_id, DB_tr_chr$Q_id)
      
    } else {
      
      QTL_un_vect <- c(QTL_un_vect, paste0('QTL_', tr_sh_id[i], '_', chr_list[c], '_1'))
      QTL_un_id <- c(QTL_un_id, DB_tr_chr$Q_id)
      
    }
    
  }
  
}

names(QTL_un_vect) <- QTL_un_id
DB$QTL_un_id <- QTL_un_vect[DB$Q_id]

DB$trait <- factor(DB$trait, levels = tr_list)
DB <- DB[order(trait, chr)]

R2 <- DB[, .(mean(R2)), by = .(QTL_un_id)]
Q_un_R2 <- R2$V1
names(Q_un_R2) <- R2$QTL_un_id
DB$QTL_un_R2 <- Q_un_R2[DB$QTL_un_id]

# change the name of the unique QTL id by adding the average cM position

cM <- DB[, .(mean(cM)), by = .(QTL_un_id)]
Q_un_cM <- cM$V1
names(Q_un_cM) <- cM$QTL_un_id
cM_vect <- as.character(round(Q_un_cM[DB$QTL_un_id]))
QTL_un_id_list <- strsplit(x = DB$QTL_un_id, split = '_')
for(i in 1:length(QTL_un_id_list)){
  QTL_un_id_list[[i]][length(QTL_un_id_list[[i]])] <- cM_vect[i]}
QTL_un_id <- lapply(X = QTL_un_id_list, FUN = function(x) paste0(x, collapse = '_'))
DB$QTL_un_id <- unlist(QTL_un_id)

#### QTL database: Add the QTL effects ####

load(file = './results/QTLxEC/QTLxEC_res.RData')

# get the information structure
pop_id <- c('GR2012', 'GR2013','KK2012', 'KK2013', 'Lata')
n_pop <- length(pop_id)
RP <- c('Grinkan', 'Grinkan', 'Kenin-Keni', 'Kenin-Keni', 'Lata')
year <- c('2012', '2013', '2012', '2013', '2013')

tr_list <- c('FLAGDD', 'PH', "NODE_N", "NODE_L", 'PAN_L', 'PED_L', 'G1000_WGH', 'G_YIELD')
n_trait <- length(tr_list)

env_lk <- c('SB1', 'SB2', 'CZ1', 'CZ2')
names(env_lk) <- paste0('E', 1:4)
env_lk2 <- c('LP', 'HP', 'KOL')
names(env_lk2) <- paste0('E', 1:3)

DB_eff <- c()

options(warn = 2)

for(p in 1:n_pop){
  
  RP_p <- RP[p]
  year_p <- year[p]
  pop_p <- pop_id[p]
  
  if(RP_p == 'Lata'){env_lk_p <- env_lk2} else {env_lk_p <- env_lk}
  
  for(t in 1:n_trait){
    
    tr <- tr_list[t]
    file_pt <- file.path('./results/QTL_analysis', pop_p, paste0('QTLGE_MPP_', tr, '_UN'))
    Qeff_f <- file.path(file_pt, 'QTL_effects.RData')
    
    if(file.exists(Qeff_f)){
      
      load(file = Qeff_f)
      n_QTL <- length(Q_eff)
      
      for(q in 1:n_QTL){
        
        # subset the specific part of the DB
        DB_ptq <- DB[pop == pop_p & trait == tr & QTL == paste0('Q', q)]
        Q_eff_q <- Q_eff[[q]]
        
        # replace Short_Kaur by ShortKaur
        rownames(Q_eff_q) <- gsub(pattern = 'Short_Kaur', replacement = 'ShortKaur',
                                  x = rownames(Q_eff_q))
        
        rownames(Q_eff_q) <- gsub(pattern = 'IS23645', replacement = 'Hafijeka',
                                  x = rownames(Q_eff_q))
        
        rownames(Q_eff_q) <- gsub(pattern = 'FaraFara', replacement = 'Fara',
                                  x = rownames(Q_eff_q))
        
        # organise effect identifiers
        eff_id <- rownames(Q_eff_q)
        eff_id <- strsplit(x = eff_id, split = '_')
        par_id <- sapply(X = eff_id, `[[`, 3)
        env_id <- env_lk_p[sapply(X = eff_id, `[[`, 4)]
        
        Q_eff_q <- data.table(par = par_id, env = env_id,
                              effect = Q_eff_q$Effect,
                              log10pval = -log10(Q_eff_q$p.val))
        
        Q_main_QxE <- QTLxEC_res[[pop_p]][[tr]]
        
        if(!is.null(Q_main_QxE)){
          
          if(tr == 'G1000_WGH' | pop_p == 'Lata'){
            
            Q_main_QxE_ptq <- Q_main_QxE[[q]]
            
          } else {
            
            Q_main_QxE_ptq <- Q_main_QxE[[1]]$Qeff_main_QxE[[q]]
            
          }
          
          Q_main_QxE_ptq <- data.table(par = rownames(Q_main_QxE_ptq),
                                       main_eff = Q_main_QxE_ptq$Effect_main,
                                       logp_main = Q_main_QxE_ptq$logP_main,
                                       logp_QxE = Q_main_QxE_ptq$logP_QxE)
          
          # insert the effect in the matrix
          par_id <- unique(par_id)
          n_par <- length(par_id)
          
          Q_eff_q2 <- c()
          
          for(p_n in 1:n_par){
            
            Q_eff_q_p <- Q_eff_q[par == par_id[p_n]]
            Q_main_p <- Q_main_QxE_ptq[par ==  par_id[p_n]]
            Q_main_p <- data.table(par = Q_main_p$par, env = 'main',
                                   effect = Q_main_p$main_eff,
                                   log10pval = Q_main_p$logp_main)
            
            Q_eff_q_p <- rbind(Q_main_p, Q_eff_q_p)
            Q_eff_q_p$log10_main <- Q_main_QxE_ptq[par ==  par_id[p_n]]$logp_main
            Q_eff_q_p$log10_QxE <- Q_main_QxE_ptq[par ==  par_id[p_n]]$logp_QxE
            Q_eff_q2 <- rbind(Q_eff_q2, Q_eff_q_p)
            
          }
          
          Q_eff_q <- Q_eff_q2
          
        } else {
          
          Q_eff_q$log10_main <- NA
          Q_eff_q$log10_QxE <- NA
          
        }
        
        # order per parents and environment
        Q_eff_q$par <- factor(Q_eff_q$par, levels = unique(par_id))
        Q_eff_q$env <- factor(Q_eff_q$env, levels = c('main', unique(env_id)))
        Q_eff_q <- Q_eff_q[order(par, env)]
        
        DB_ptq <- data.table(DB_ptq, Q_eff_q)
        DB_eff <- rbind(DB_eff, DB_ptq)
        
      }
      
    }
    
  }
  
}

# change the QTL parent allele sign name
colnames(DB_eff)[20] <- 'log10pval_Qp'

# Correct the names Hafijega for Hafijeka
DB_eff$par <- as.character(DB_eff$par)
DB_eff$par[DB_eff$par == 'Hafijega'] <-  'Hafijeka'
DB_eff$par <- factor(DB_eff$par)

#### QTL database: Add the QTLxEC effects ####

# extend the DB with empty columns for EC effects

EC_id <- c('rain', 'hum', 'VPD', 'SVP', 'ETP', 'PETP', 'Tmin', 'Tmax',
           'Trange', 'DD', 'FRUE', 'hSun', 'photoperiod', 'solarRad', 'photothermal')

n_EC <- length(EC_id)

ext_col_nm <- paste0(rep(EC_id, 3), '_', rep(c('main', 'B', 'sign'), each = n_EC))
ext_inf_mat <- data.table(matrix(NA, nrow = nrow(DB_eff), ncol = length(ext_col_nm)))
colnames(ext_inf_mat) <- ext_col_nm
ext_inf_mat <- apply(X = ext_inf_mat, 2 , FUN = as.numeric)
n_var_DB_org <- ncol(DB_eff)

DB <- data.table(DB_eff, ext_inf_mat)
rm(DB_eff)

# load QTLxEC results
load(file = './results/QTLxEC/QTLxEC_res.RData')

# get the information structure
pop_id <- c('GR2012', 'GR2013','KK2012', 'KK2013')
n_pop <- length(pop_id)

trait <- c('FLAGDD', 'PH', "NODE_N", "NODE_L", 'PAN_L', 'PED_L', 'G_YIELD')
n_trait <- length(trait)

EC_v_main <- EC_v_B <- EC_v_sign <- rep(NA, n_EC)
names(EC_v_main) <- names(EC_v_B) <- names(EC_v_sign) <- EC_id

for(p in 1:n_pop){
  
  pop <- pop_id[p]
  
  for(t in 1:n_trait){
    
    tr <- trait[t]
    
    Q_res_pt <- QTLxEC_res[[pop]][[tr]]
    EC_pt <- names(Q_res_pt)
    
    if(any(!(EC_pt %in% EC_id))){ stop('EC id not present') }
    
    for(ec in 1:length(Q_res_pt)){
      
      EC_res <- Q_res_pt[[ec]][[2]]
      
      for(q in 1:length(EC_res)){
        
        EC_res_q <- EC_res[[q]]
        
        # restrict to the parent with significant effect
        EC_res_q <- EC_res_q[!is.na(EC_res_q[, 3]), ]
        
        p_id <- rownames(EC_res_q)
        
        # change the problematic p_nm
        p_id[p_id == 'Short_Kaur'] <- 'ShortKaur'
        p_id[p_id == 'IS23645'] <- 'Hafijeka'
        p_id[p_id == 'FaraFara'] <- 'Fara'
        
        # p_nb <- 1
        
        for(p_nb in 1:length(p_id)){
          
          # find the line in the DB that correspond to those results
          line_s <- which(DB$pop == pop & DB$trait == tr & DB$QTL == paste0('Q', q) &
                            DB$par == p_id[p_nb] & DB$env == 'main')
          
          pos_int <- n_var_DB_org + which(ext_col_nm == paste0(EC_pt[ec], '_main'))
          pos_B <- n_var_DB_org + which(ext_col_nm == paste0(EC_pt[ec], '_B'))
          pos_sign <- n_var_DB_org + which(ext_col_nm == paste0(EC_pt[ec], '_sign'))
          
          DB[line_s, pos_int] <- EC_res_q[p_nb, 1] # intercept
          DB[line_s, pos_B] <- EC_res_q[p_nb, 3] # Beta
          DB[line_s, pos_sign] <- EC_res_q[p_nb, 4] # sign
          
        }
        
      }
      
    }
    
  }
  
}

# change the trait name
tr_lk <- c('FLAG', 'PH', 'NODE_N', 'NODE_L', 'PED', 'PAN', 'GWGH', 'YIELD')
names(tr_lk) <- c('FLAGDD', 'PH', 'NODE_N', 'NODE_L', 'PED_L', 'PAN_L', 'G1000_WGH', 'G_YIELD')
DB$trait <- as.character(DB$trait)
DB$trait <- tr_lk[DB$trait]
DB$trait <- factor(DB$trait, levels = c('FLAG', 'PH', 'NODE_N', 'NODE_L', 'PED', 'PAN', 'GWGH', 'YIELD'))

# save DB
save(DB, file = './results/database/database.RData')