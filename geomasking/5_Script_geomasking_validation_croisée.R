#===========================================================================================================#
# Script : Script principal geomasking

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la cartographie des propriétés des sols en utilisant la méthode Random Forest 
#               Krigeage ordinaire et Krigeage avec Dérive Externe avec un geomasking

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    14-10-2025

# Entrees :     Données ponctuelles de la BDAT et IGCS floutées sur les propriétés des sols

# Sorties :     Carte des propriétés des sols sur la zone agricole

# Modification : 
#===========================================================================================================#

#================================================DEBUT DU SCRIPT====================================================#


#Liste des packages utilisés -----

library(sf) # manipuler les datas SIG vector
library (sp)
library(terra) # manipuler des données raster-pixels
library(raster) # manipuler des données raster
library(tidyr) # reroganisation des données
library(dplyr)  # grammaire de manipulation de données
library(randomForest) # randomForest pour utiliser tuneRandomForest
library(Boruta)  # implementation  de selection de variables de randomForest
library(ranger) # quantile randomForest, version optimée de randomForest
library(caret) # modelisation creation de modeles prédictifs avec tuning validation-croisée
library(tuneRanger) # pour le tuning de randomForest
library(readxl) # lecture de fichier excel
library(ggplot2) # graphiques
library(tmap)  # outil d'édution de carte
library(ggnewscale) # pour les graphiques avec plusieurs échelles de couleurs
library(ggpubr) # pour les graphiques
library(purrr)
library(foreach) # boucles optimisées de R plus facile à coder voir option combine
library(doParallel) # pour le calcul parallèle

# INLA SPDE avec inla bru, approche bayesienne de la géostatistique
library(INLA)
library(inlabru)
library(forcats)

#Liste des Fonctions----

##Indicateurs de qualité de prédictions----

#Ces indicateurs sont utilisés pour évaluer la qualité des modèles de prédiction et sont calculés en comparant les valeurs prédites (x) aux valeurs observées (y).
#lors du processus de validation croisée.

Myeval <- function(x, y){
  
  # Erreur Moyenne (EM)
  EM <- round(mean(y - x, na.rm = TRUE), digits = 2)
  
  # Racine de l'erreur quadratique moyenne (EQM)
  REQM <-   round(sqrt(mean((y - x)^2, na.rm = TRUE)), digits = 2)
  
  # Carré de la corrélation de Pearson (r2) https://fr.wikipedia.org/wiki/Coefficient_de_d%C3%A9termination
  r2 <-  round((cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2), digits = 2)
  
  # Coefficient d'efficacité du modèle de Nash-Sutcliffe (NSE) https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
  SSE <- sum((y - x) ^ 2, na.rm = T)
  SST <- sum((y - mean(y, na.rm = T)) ^ 2, na.rm = T)
  NSE <- round((1 - SSE/SST), digits = 2)
  
  # Coefficient de corrélation de concordance (CCC) https://thedatascientist.com/concordance-correlation-coefficient/#:~:text=Definition%20and%20Purpose,and%200%20signifies%20no%20agreement.
  n <- length(x)
  sdx <- sd(x, na.rm = T)
  sdy <- sd(y, na.rm = T)
  r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
  # scale shift
  v <- sdx / sdy
  sx2 <- var(x, na.rm = T) * (n - 1) / n
  sy2 <- var(y, na.rm = T) * (n - 1) / n
  # location shift relative to scale
  u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
  Cb <- ((v + 1 / v + u^2)/2)^-1
  rCb <- r * Cb
  CCC <- round(rCb, digits = 2)
  
  Cb <- round(Cb, digits = 2)
  r <- round(r, digits = 2)
  
  # retourner les résultats dans un datarame
  evalRes <- data.frame(EM = EM, REQM = REQM, r = r, r2 = r2, NSE = NSE, CCC = CCC, Cb = Cb)
  
  return(evalRes)
}

##Fonction pour extraire les covariables à partir des coordonnées spatiales----
extract_covs <- function(tab, coord_names) {
  v <- terra::vect(tab, geom = coord_names, crs = crs_epsg)
  vals <- terra::extract(st, v)
  vals <- vals[, -1, drop = FALSE]       # retirer col. id interne terra
  dplyr::bind_cols(tab, tibble::as_tibble(vals))
}

##Fonction pour la calibration du RF-----
run_rf_geomask <- function(data_train, data_test, cov_brt, name, ntree, kmax,
                           NomsCoord_train = c("x_moved","y_moved"),
                           NomsCoord_test  = c("x","y")) {
  
  rf_task <- mlr::makeRegrTask(data = data_train[, c(name, cov_brt)], target = name)
  res_tune <- tuneRanger::tuneRanger(rf_task, num.trees = ntree, iters = 10,
                                     num.threads = kmax, show.info = FALSE)
  
  rf_model <- ranger::ranger(
    formula = stats::as.formula(paste0(name, " ~ .")),
    data = data_train[, c(name, cov_brt)],
    num.trees = ntree,
    mtry = res_tune$recommended.pars$mtry,
    min.node.size = res_tune$recommended.pars$min.node.size,
    quantreg = FALSE,
    max.depth = 15,
    importance = "permutation",
    scale.permutation.importance = TRUE,
    keep.inbag = FALSE
  )
  
  predCal <- rf_model$predictions
  preds   <- predict(rf_model, data = data_test[, cov_brt], num.threads = kmax)$predictions
  preds   <- round(preds, 2)
  
  test_out <- data_test[, c("id", name, "INSEE_COM", NomsCoord_test, cov_brt)]
  test_out$pred <- preds
  
  eval <- Myeval(preds, test_out[[name]]) %>% dplyr::mutate(method = "RF_PC")
  list(
    rf_model   = rf_model,
    best_pars  = res_tune$recommended.pars,
    predCal    = predCal,
    evaluation = eval,
    detail     = test_out[, c("id", name, "INSEE_COM", NomsCoord_test, "pred")]
  )
}

##Fonction pour la calibration de KO et KED----
run_inla_spde_core_geomask <- function(dataINLA, data_test, name,
                                       type = c("KO","KED"),
                                       NomsCoord_train = c("x_moved","y_moved"),
                                       NomsCoord_test  = c("x","y"),
                                       crs_str = "epsg:2154") {
  type <- match.arg(type)
  dataINLA$elt <- dataINLA[, name]
  
  coords_train <- as.matrix(dataINLA[, NomsCoord_train])
  coords_test  <- as.matrix(data_test[,  NomsCoord_test])
  
  sp::coordinates(dataINLA) <- NomsCoord_train
  sp::proj4string(dataINLA) <- sp::CRS(crs_str)
  sp::coordinates(data_test) <- NomsCoord_test
  sp::proj4string(data_test) <- sp::CRS(crs_str)
  
  max.edge    <- diff(range(coords_train[,1])) / (3 * 5)
  bound.outer <- diff(range(coords_train[,1])) / 3
  mesh3 <- INLA::inla.mesh.2d(
    loc      = coords_train,
    max.edge = c(1, 2) * max.edge,
    offset   = c(max.edge, bound.outer),
    cutoff   = max.edge / 10
  )
  matern <- INLA::inla.spde2.pcmatern(
    mesh3, alpha = 2, prior.sigma = c(10, 0.01), prior.range = c(500, 0.01)
  )
  
  cmp <- if (type == "KO") {
    elt ~ Intercept(1) + field(coordinates, model = matern)
  } else {
    if (!"pred" %in% names(dataINLA)) stop("KED: 'pred' manquant dans dataINLA (train).")
    elt ~ Intercept(1) + rfpred(pred, model = "linear") + field(coordinates, model = matern)
  }
  cmpPred <- if (type == "KO") { ~ Intercept + field } else {
    if (!"pred" %in% names(data_test)) stop("KED: 'pred' manquant dans data_test (test).")
    ~ Intercept + rfpred + field
  }
  
  fit <- inlabru::bru(cmp, data = dataINLA, family = "Gaussian",
                      options = list(control.inla = list(int.strategy = "eb"),
                                     verbose = FALSE))
  gc()
  predObj <- predict(fit, data_test, cmpPred)
  predObj$median
}

# 1. Définition des variables ------

name="arg" # changer la variable d'intérêt au besoin 
kmax= 23 # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
k=10 # pour la validation croisée
nsim=100 # for bayesian inla simulation
NomsCoord1 <- c("x", "y")              
NomsCoord2 <- c("x_moved", "y_moved") 
crs_epsg   <- "EPSG:2154"
out_dir <- "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/"
predict_script <- "geomasking/geomasking_mapping.R" 

# Séquence des distances
d_seq <- seq(100, 2500, by = 100)

#Preparation des données pour la spatialisation----
##Importation des données----

rast_za <- rast(out_dir,"rast_za.tif") # raster de la zone agricole

cov_brt <- readRDS(paste0(out_dir, name, "_cov_brt.rds"))


##Extraction des matrices de covariables pour les données ponctuelles

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

#liste des fichiers de covariables
l<- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)

#stack des covariables
st <- rast(l)

gXY <- as.data.frame(st , xy=TRUE) %>%
  na.omit( )

# Harmonisation / clés / contrôle coords
metrics_all_d <- foreach::foreach(
  d = d_seq,
  .packages = c("dplyr","tidyr","terra","sp","INLA","inlabru","ranger","mlr","tuneRanger","caret","tibble")
) %do% {
  
  message("\n===================== d = ", d, " m =====================")
  
  #charger les données géomasquées pour cette distance 
  dtTB_path <- paste0(out_dir,"geomasked", name, "_dist", d, ".rds")
  stopifnot(file.exists(dtTB_path))
  dtTB <- readRDS(dtTB_path)
  
  #harmonisation / clés / coords
  dt <- dtTB %>%
    dplyr::mutate(
      id = dplyr::coalesce(as.numeric(ID), as.numeric(id_profil), dplyr::row_number()),
      INSEE_COM = as.character(INSEE_COM)
    ) %>%
    dplyr::select(id, INSEE_COM, source, annee = annee, !!rlang::sym(name), x, y, x_moved, y_moved) %>%
    dplyr::filter(is.finite(x), is.finite(y), is.finite(x_moved), is.finite(y_moved))
  
  #xtraction covariables aux deux jeux de coords 
  orig_df <- dt %>%
    dplyr::select(id, INSEE_COM, source, annee, !!rlang::sym(name), dplyr::all_of(NomsCoord1)) %>%
    extract_covs(NomsCoord1)
  
  mask_df <- dt %>%
    dplyr::select(id, INSEE_COM, source, annee, !!rlang::sym(name), dplyr::all_of(NomsCoord2)) %>%
    extract_covs(NomsCoord2)
  
  # filtre NA synchronisé 
  ok_orig <- stats::complete.cases(orig_df[, c(name, cov_brt), drop = FALSE])
  ok_mask <- stats::complete.cases(mask_df[, c(name, cov_brt), drop = FALSE])
  ok_both <- ok_orig & ok_mask
  orig_df <- orig_df[ok_both, , drop = FALSE]
  mask_df <- mask_df[ok_both, , drop = FALSE]
  
  stopifnot(nrow(orig_df) == nrow(mask_df))
  stopifnot(all(orig_df$id == mask_df$id))
  
  # garder colonnes utiles 
  orig_df <- orig_df %>% dplyr::select(id, INSEE_COM, source, annee, !!rlang::sym(name),
                                       dplyr::all_of(NomsCoord1), dplyr::all_of(cov_brt))
  mask_df <- mask_df %>% dplyr::select(id, INSEE_COM, source, annee, !!rlang::sym(name),
                                       dplyr::all_of(NomsCoord2), dplyr::all_of(cov_brt))
  
  # folds identiques
  set.seed(1001)  # même folds par distance
  folds <- caret::createFolds(orig_df$id, k = k, list = TRUE, returnTrain = FALSE)
  
  # CV 
  results_all <- vector("list", length = k)
  
  for (fold_idx in seq_len(k)) {
    message("  Fold ", fold_idx, "/", k)
    test_idx  <- folds[[fold_idx]]
    train_idx <- unlist(folds[-fold_idx])
    
    data_train <- mask_df[train_idx, ]
    data_test  <- orig_df[test_idx,  ]
    
    rf_res <- run_rf_geomask(
      data_train = data_train,
      data_test  = data_test,
      cov_brt    = cov_brt,
      name       = name,
      ntree      = ntree,
      kmax       = kmax,
      NomsCoord_train = c("x_moved","y_moved"),
      NomsCoord_test  = c("x","y")
    )
    test_with_rf <- rf_res$detail
    
    ko_pred <- run_inla_spde_core_geomask(
      dataINLA = data_train, data_test = data_test, name = name, type = "KO",
      NomsCoord_train = c("x_moved","y_moved"), NomsCoord_test = c("x","y"),
      crs_str = "epsg:2154"
    )
    
    data_train_ked <- data_train; data_train_ked$pred <- rf_res$predCal
    ked_pred <- run_inla_spde_core_geomask(
      dataINLA = data_train_ked, data_test = test_with_rf, name = name, type = "KED",
      NomsCoord_train = c("x_moved","y_moved"), NomsCoord_test = c("x","y"),
      crs_str = "epsg:2154"
    )
    
    results_all[[fold_idx]] <- rf_res$detail %>%
      dplyr::mutate(approach = "Ponctuelle", type_val = "Classique", fold = fold_idx,
                    predKO = ko_pred, predKED = ked_pred) %>%
      dplyr::rename(obs = !!rlang::sym(name)) %>%
      dplyr::select(id, approach, type_val, fold, obs, pred, predKO, predKED, INSEE_COM, x, y)
  }
  
  resuXval_geomask <- dplyr::bind_rows(results_all)
  
  # métriques par distance
  resuXvalQRF <-  Myeval(resuXval_geomask$pred,    resuXval_geomask$obs )
  resuXvalKO  <-  Myeval(resuXval_geomask$predKO,  resuXval_geomask$obs )
  resuXvalKED <-  Myeval(resuXval_geomask$predKED, resuXval_geomask$obs )
  
  metrics_overall <- dplyr::bind_rows(
    dplyr::mutate(resuXvalQRF, method = "RF"),
    dplyr::mutate(resuXvalKO,  method = "KO"),
    dplyr::mutate(resuXvalKED, method = "KED")
  ) %>% dplyr::mutate(dist = d, .before = 1)
  
  #sauvegardes fichiers de cette distance
  saveRDS(resuXval_geomask, file = paste0(out_dir, "geomasked_Xval_", name, "_dist", d, ".rds"))
  write.csv(metrics_overall, file = paste0(out_dir, "metrics_overall_", name, "_dist", d, ".csv"),
            row.names = FALSE)
  
  metrics_overall  # valeur renvoyée pour cette distance
  
  rf_fit <- run_rf_geomask(
    data_train = mask_df,
    data_test  = mask_df[0, ],
    cov_brt    = cov_brt,                  
    name       = name,
    ntree      = ntree,
    kmax       = kmax,
    NomsCoord_train = c("x_moved","y_moved"),
    NomsCoord_test  = c("x","y")
  )
  
  rf_model <- rf_fit$rf_model  
  best_pars <- rf_fit$best_pars
  
  source(predict_script, local = TRUE)
}

metrics_all <- dplyr::bind_rows(metrics_all_d)


metrics_long <- metrics_all %>%
  dplyr::select(dist, method, NSE, CCC, REQM) %>%
  tidyr::pivot_longer(cols = c(NSE, CCC, REQM),
                      names_to = "metric", values_to = "value") %>%
  dplyr::mutate(
    method = factor(method, levels = c("RF","KO","KED")),
    metric = factor(metric, levels = c("NSE","CCC","REQM")),
    dist   = as.numeric(dist)
  )

p_line <- ggplot(metrics_long, aes(x = dist, y = value, group = 1)) +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(metric ~ method, scales = "free_y") +
  labs(x = "Distance de géomasking (m)", y = "Valeur de l'indicateur",
       title = paste0("Evolution des indicateurs en fontion de la distance selon les modèles pour le ", name)) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))
print(p_line)


#================================================FIN DU SCRIPT====================================================#
