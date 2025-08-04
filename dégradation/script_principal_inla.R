# Script principal pour faire de la validation croisée en intégrant la dégradation 

#L'idée est de partir des données ponctuelles (contenant les covaribles et les valeurs de la propriété cible).
#Suivant les différentes voies de validation croisée (classique ou spatiale), on définit des folds et à chaque itérartion,
#on calivre sur les k-1 folds et on teste sur le fold restant. Mais la calibration ne se fait pas sur l'intégralité des points
#des k-1 folds mais plutut sur un échantillon ce qui permet de simuler la dégradation de l'information.

#Pour chaque échantillon on implémente d'abord le Random Forest pour avoir les prédictions, les metriques, les identifiant' numéro folds,
#id, numéro repetition) pour ensuite les utiliser et eviter la perte d'information afin d'implémenter le KO et le KED.

#1. Chargement des packages----

library(sf)
# library(tmap)
library(readxl)
library(tidyr)
library(dplyr)
library(foreach)
library(raster) # pourquoi Raster????
library(purrr)
# library(ggpubr)
# library(ggnewscale)
# library(doParallel)
library(tuneRanger)
library(INLA)
library(inlabru)
library(iml)
library(mlr)
library(caret)

# setwd("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R")

#2. Chargement des fonctions RF, INLA  Myeval et dataINLA----

source("dégradation/fonction_RF.R")

source("dégradation/fonction_inla.R")
source("dégradation/fonction_geomasking.R")
#Myeval
Myeval <- function(x, y){
  ME <- mean(y - x, na.rm = TRUE)
  RMSE <- sqrt(mean((y - x)^2, na.rm = TRUE))
  MAE <- mean(abs(y - x), na.rm = TRUE)
  r2 <- (cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2)
  SSE <- sum((y - x)^2, na.rm = T)
  SST <- sum((y - mean(y, na.rm = T))^2, na.rm = T)
  NSE <- 1 - SSE/SST
  r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
  v <- sd(x, na.rm = T) / sd(y, na.rm = T)
  sx2 <- var(x, na.rm = T) * (length(x) - 1) / length(x)
  sy2 <- var(y, na.rm = T) * (length(y) - 1) / length(y)
  u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
  Cb <- ((v + 1/v + u^2)/2)^-1
  rCb <- r * Cb
  CCC <- rCb
  data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE, CCC = CCC, Cb = Cb)
}


# dataINLA

# Paramètres

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

name <- args[1]  #arg
kmax <- 30
ntree <- 350
NomsCoord <- c("x", "y")
sample_sizes <- args[2] # c(500,1000 ) # c(600,800,1000,1200,1300,1400,1600,1800,2000,3000,4000,5000,6000,7000,7600)
repets <- args[3]
types_validation <- c("Classique", "Spatiale")
drive = "/media/communs_infosol/" # ou "Y:/"
DistanceGeomasking = 0

#3. Chargement des données---- 

com <- st_read( paste0(drive, "BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp") )
centroides_communes <- st_centroid(com) %>% dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID) %>% st_drop_geometry()

datacov <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/donnees_ponctuelles", name, ".rds"))
moyenne_covariable <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/moyenne_covariable", name, ".rds"))

cov_brt <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/", name, "_cov_brt.rds"))

#4. loop of the ----


cat("\n==============  Starting loops ===============\n")

bru_safe_inla(multicore = FALSE)

pred_RF_full <-  foreach(
  
  n = sample_sizes,
  .combine = rbind.data.frame ,
  .errorhandling='pass'
  
  ) %do% {
    cat("\n============== Taille d'échantillon :", n, "===============\n")
    
    k <- ceiling(nrow(datacov) / 1000)
    foreach (
      
      rep = 1:repets,
      .combine = rbind.data.frame,
      
      .errorhandling='pass'
      
             ) %do%  {
      cat("\n---- Répétition :", rep, "----\n")
      
      set.seed(1001 + rep)
      
      foreach (type_val = types_validation,
               .combine = rbind.data.frame,
               .errorhandling='pass'
               ) %do% {
        cat(">> Validation :", type_val, "\n")
        
        folds <- if (type_val == "Classique") {
          createFolds(datacov[[name]], k = k, list = TRUE)
        } else {
          communes <- unique(datacov$INSEE_COM)
          split(communes, sample(rep(1:k, length.out = length(communes))))
        }
        
        foreach (fold_idx = 1:k,
                 .combine = rbind.data.frame,
                 .errorhandling='pass'
                 ) %do%  {
          
          cat("   > Fold", fold_idx, "sur", k, "\n")
          
          if (type_val == "Classique") {
            idx_test <- folds[[fold_idx]]
            idx_calib <- unlist(folds[-fold_idx])
            test_data <- datacov[idx_test, ]
            calib_pool <- datacov[idx_calib, ]
          } else {
            com_test <- folds[[fold_idx]]
            idx_test <- which(datacov$INSEE_COM %in% com_test)
            idx_calib <- which(!datacov$INSEE_COM %in% com_test)
            test_data <- datacov[idx_test, ]
            calib_pool <- datacov[idx_calib, ]
          }
          
          if (nrow(calib_pool) > n) calib_points <- calib_pool %>% sample_n(n) else calib_points <- calib_pool
          
          if ( DistanceGeomasking >0 ) {
            calib_points_geomasked <- geomasking(calib_points,geomasking)
          }
          
          
          # RF Ponctuelle
          res_rf_p <- run_rf(
            "Ponctuelle",
            type_val,
            calib_points,
            test_data,
            cov_brt,
            moyenne_covariable,
            name,
            ntree,
            kmax,
            NomsCoord
          )
          
          calib_points$pred  <- res_rf_p$predCal
          
          
          # agrgégation de la variable cible pour la méthode Centroide
          agg_target <- calib_points %>%
            dplyr::select(INSEE_COM, all_of(name)) %>%
            group_by(INSEE_COM) %>%
            dplyr::summarise(!!name := mean(.data[[name]], na.rm = TRUE),
                      .groups = "drop")
          
          train_aggr <- moyenne_covariable %>%
            inner_join(agg_target, 
                       by = "INSEE_COM"
                       )
          
          res_rf_c <- run_rf(
            "Centroide",
            type_val,
            train_aggr,
            test_data,
            cov_brt,
            moyenne_covariable,
            name,
            ntree,
            kmax,
            NomsCoord
          )
          
          
          train_aggr$pred = res_rf_c$predCal
          
          
          
          cat("=> Exécution INLA ", type_val, "\n")
          
          res_ko_p <- run_inla_spde_core(calib_points,
                                         test_data,
                                         name,
                                         "KO",
                                         "Ponctuelle",
                                         type_val)
          
          
          train_aggr <- train_aggr  %>%
            inner_join(centroides_communes, by = "INSEE_COM") %>%
            dplyr::rename(x = X, y = Y) %>%
            mutate(id = INSEE_COM)
          
          
          res_ko_c <- run_inla_spde_core(as.data.frame(train_aggr),
                                         test_data,
                                         name,
                                         "KO",
                                         "Centroide",
                                         type_val)
          
          
          res_ked_p <- run_inla_spde_core(calib_points,
                                          res_rf_p$detail,
                                          name,
                                          "KED",
                                          "Ponctuelle",
                                          type_val)
          
          gc()
          
          res_ked_c <- run_inla_spde_core(as.data.frame(train_aggr),
                                          res_rf_c$detail ,
                                          
                                          name,
                                          "KED",
                                          "Centroide",
                                          type_val)
          
          gc()
          
          
          
          # Stockage des prédictions RF
          # results_rf_all[[length(results_rf_all) + 1]] <-
            
            bind_rows(
            res_rf_p$detail %>% 
              mutate(approach = "Ponctuelle", 
                     type_val = type_val,
                     sample_size = n,
                     rep = rep,
                     fold = fold_idx,
                     predKO = res_ko_p,
                     predKED = res_ked_p
                     ) %>%
              dplyr::rename( obs := !!name  ) %>%
              dplyr::select(id,approach ,type_val,
                            sample_size,rep,fold,
                            obs,pred ,predKO,predKED),
            res_rf_c$detail %>% 
              mutate(approach = "Centroide",
                     type_val = type_val, 
                     sample_size = n,
                     rep = rep,
                     fold = fold_idx,
                     predKO = res_ko_c,
                     predKED = res_ked_c
              )%>%
              dplyr::rename( obs := !!name ) %>%
              dplyr::select(id,approach ,type_val,
                            sample_size,rep,fold,
                            obs,pred,predKO,predKED)
          )
          
          
        }
      }
      
      
    }
  }


cat("FIN DES CALCULS--------------------")

  
  # Fusion de toutes les prédictions RF
saveRDS(pred_RF_full, 
        paste0(
          "output/Xval_",
          name,
          paste0(sample_sizes,collapse = "_") ,
          ".rds")
        )
  
