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
library(raster)
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
prepare_dataINLA <- function(type, approach, data_train, data_test, name) {
  pred_rf_col <- paste0("predRF_", substr(approach, 1, 1), substr(type_val, 1, 1))
  
  if (approach == "Centroide") {
    agg <- data_train %>%
      group_by(INSEE_COM) %>%
      summarise(
        activ = mean(.data[[name]], na.rm = TRUE),
        pred = if (type == "KED") mean(.data[[pred_rf_col]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      )
    train_aggr <- centroides_communes %>%
      inner_join(agg, by = "INSEE_COM") %>%
      rename(x = X, y = Y)%>%
      mutate(id = INSEE_COM)
  } else {
    train_aggr <- data_train %>%
      transmute(
        id = id,
        INSEE_COM = INSEE_COM,
        x = .data[[NomsCoord[1]]],
        y = .data[[NomsCoord[2]]],
        activ = .data[[name]],
        pred = if (type == "KED") .data[[pred_rf_col]] else NA_real_
      ) %>%
      filter(!is.na(x), !is.na(y), !is.na(activ), if (type == "KED") !is.na(pred) else TRUE)
  }
  
  coords_test <- data_test[, NomsCoord, drop= FALSE]
  x_test <- coords_test[[NomsCoord[1]]]
  y_test <- coords_test[[NomsCoord[2]]]
  
  elt_test <- rep(NA_real_, nrow(coords_test))
  pred_test <- if (type == "KED") data_test[[pred_rf_col]] else NULL
  
  if (!"id" %in% names(test_data)) {
    test_data <- test_data %>% mutate(id = INSEE_COM)
  }
  if (!"id" %in% names(calib_data_rf)) {
    calib_data_rf <- calib_data_rf %>% mutate(id = INSEE_COM)
  }
  
  #création datINLA
  dataINLA <- data.frame(
    id = c(train_aggr$id, data_test$id),
    x = c(train_aggr$x, x_test),
    y = c(train_aggr$y, y_test),
    elt = c(train_aggr$activ, elt_test),
    id_point = c(rep(NA, nrow(train_aggr)), data_test$id)
  )
  
  if (type == "KED") {
    dataINLA$pred <- c(train_aggr$pred, pred_test)
  }
  
  return(dataINLA)
}



# Paramètres
name <- "arg"
kmax <- 10
ntree <- 350
NomsCoord <- c("x", "y")
sample_sizes <- c(600,7600 ) # c(600,800,1000,1200,1300,1400,1600,1800,2000,3000,4000,5000,6000,7000,7600)
repets <- 2
types_validation <- c("Classique", "Spatiale")
drive = "/media/communs_infosol/" # ou "Y:/"

#3. Chargement des données---- 

com <- st_read( paste0(drive, "BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp") )
centroides_communes <- st_centroid(com) %>% dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID) %>% st_drop_geometry()

datacov <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/donnees_ponctuelles", name, ".rds"))
moyenne_covariable <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/moyenne_covariable", name, ".rds"))

cov_brt <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/", name, "_cov_brt.rds"))

#4. RF pour la validation croisée avec dégradation----

#initialisation
results_rf_all <- list()
results_rf_all_metrics <- list()

for (n in sample_sizes) {
    cat("\n============== Taille d'échantillon :", n, "===============\n")
    
  k <- ceiling(nrow(datacov) / 1000)
  
  for (rep in 1:repets) {
    cat("\n---- Répétition :", rep, "----\n")
    
    set.seed(1000 + rep)
    
    for (type_val in types_validation) {
      cat(">> Validation :", type_val, "\n")
      
      folds <- if (type_val == "Classique") {
        createFolds(datacov[[name]], k = k, list = TRUE)
      } else {
        communes <- unique(datacov$INSEE_COM)
        split(communes, sample(rep(1:k, length.out = length(communes))))
      }
      
      for (fold_idx in 1:k) {
        
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
        
        if (nrow(calib_pool) < n) next
        calib_points <- calib_pool %>% sample_n(n)
        
        # RF Ponctuelle
        res_rf_p <- run_rf("Ponctuelle", type_val, calib_points, test_data, cov_brt, moyenne_covariable, name, ntree, kmax,NomsCoord)
        
        # agrgégation de la variable cible pour la méthode Centroide
        agg_target <- calib_points %>%
          dplyr::select(INSEE_COM, all_of(name)) %>%
          group_by(INSEE_COM) %>%
          summarise(!!name := mean(.data[[name]], na.rm = TRUE), .groups = "drop")
        
        train_aggr <- moyenne_covariable %>%
          inner_join(agg_target, by = "INSEE_COM")
        
        res_rf_c <- run_rf("Centroide", type_val, train_aggr, test_data, cov_brt, moyenne_covariable, name, ntree, kmax,NomsCoord)
        
        # Stockage des prédictions RF
        results_rf_all[[length(results_rf_all) + 1]] <- bind_rows(
          res_rf_p$detail %>% mutate(approach = "Ponctuelle", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_rf_c$detail %>% mutate(approach = "Centroide", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx)
        )
        
        # Stockage des métriques RF
        results_rf_all_metrics[[length(results_rf_all_metrics) + 1]] <- bind_rows(
          res_rf_p$evaluation %>% mutate(approach = "Ponctuelle", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_rf_c$evaluation %>% mutate(approach = "Centroide", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx)
        )
      }
    }
    
    
  }
}

# Fusion de toutes les prédictions RF
pred_RF_full <- bind_rows(results_rf_all)
saveRDS(pred_RF_full, "output/pred_RF_full.rds")

# Fusion de toutes les métriques RF
metrics_RF_full <- bind_rows(results_rf_all_metrics)
saveRDS(metrics_RF_full, "output/metrics_RF_full.rds")

#5. INLA pour la validation croisée avec dégradation----

#Initialisation
results_inla_all_preds <- list()
results_inla_all_metrics <- list()

for (n in sample_sizes) {
  cat("\n============== Taille d'échantillon :", n, "===============\n")
  k <- ceiling(nrow(datacov) / 1000)
  
  for (rep in 1:repets) {
    cat("\n---- Répétition :", rep, "----\n")
    
    
    for (type_val in types_validation) {
      
      cat(">> Validation :", type_val, "\n")
      
      set.seed(1000 + rep)
      
      folds <- if (type_val == "Classique") {
        createFolds(datacov[[name]], k = k, list = TRUE)
      } else {
        communes <- unique(datacov$INSEE_COM)
        split(communes, sample(rep(1:k, length.out = length(communes))))
      }
      
      for (fold_idx in 1:k) {
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
        
        if (nrow(calib_pool) < n) next
        calib_points <- calib_pool %>% sample_n(n)
        
        # Extraction des prédictions RF
        pred_rf_fold <- pred_RF_full %>% filter(sample_size == n, rep == rep, type_val == type_val, fold == fold_idx)
        
        pred_col_pc <- if (type_val == "Classique") "predRF_PC" else "predRF_PS"
        pred_col_cc <- if (type_val == "Classique") "predRF_CC" else "predRF_CS"
        
        test_data_rf <- left_join(
          test_data,
          pred_rf_fold %>% dplyr::select(all_of(c("id", pred_col_pc, pred_col_cc))),
          by = "id"
        )
        
        calib_data_rf <- left_join(calib_points, pred_rf_fold[, c("id", pred_col_pc, pred_col_cc)], by = "id")
        
      
        # Exécution INLA
        
        cat("=> Exécution INLA", type, "-", approach, "\n")
        res_ko_p <- run_inla_spde_core(prepare_dataINLA("KO", "Ponctuelle", calib_points, test_data, name), test_data, name, "KO", "Ponctuelle", type_val)
        
        cat("=> Exécution INLA", type, "-", approach, "\n")
        res_ko_c <- run_inla_spde_core(prepare_dataINLA("KO", "Centroide", calib_points, test_data, name), test_data, name, "KO", "Centroide", type_val)
        
        cat("=> Exécution INLA", type, "-", approach, "\n")
        res_ked_p <- run_inla_spde_core(prepare_dataINLA("KED", "Ponctuelle", calib_data_rf, test_data_rf, name), test_data_rf, name, "KED", "Ponctuelle", type_val)
        
        cat("=> Exécution INLA", type, "-", approach, "\n")
        res_ked_c <- run_inla_spde_core(prepare_dataINLA("KED", "Centroide", calib_data_rf, test_data_rf, name), test_data_rf, name, "KED", "Centroide", type_val)
      
        
        # Stockage des prédictions INLA
        results_inla_all_preds[[length(results_inla_all_preds) + 1]] <- bind_rows(
          res_ko_p$detail %>% mutate(approach = "Ponctuelle", method = "INLA_KO", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_ko_c$detail %>% mutate(approach = "Centroide", method = "INLA_KO", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_ked_p$detail %>% mutate(approach = "Ponctuelle", method = "INLA_KED", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_ked_c$detail %>% mutate(approach = "Centroide", method = "INLA_KED", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx)
        )
        
        # Stockage des métriques INLA
        results_inla_all_metrics[[length(results_inla_all_metrics) + 1]] <- bind_rows(
          res_ko_p$evaluation %>% mutate(approach = "Ponctuelle", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_ko_c$evaluation %>% mutate(approach = "Centroide", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_ked_p$evaluation %>% mutate(approach = "Ponctuelle", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx),
          res_ked_c$evaluation %>% mutate(approach = "Centroide", type_val = type_val, sample_size = n, rep = rep, fold = fold_idx)
        )
      }
    }
  }
}

# Fusion des métriques INLA
metrics_INLA_full <- bind_rows(results_inla_all_metrics)
saveRDS(metrics_INLA_full, "output/metrics_INLA_full.rds")

pred_INLA_full <- bind_rows(results_inla_all_preds)
saveRDS(pred_INLA_full, "output/pred_INLA_full.rds")
