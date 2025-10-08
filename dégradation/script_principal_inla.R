#=============================================================================================================================#
# Script :      Script principale pour la validation croisée intégrant un processus de dégradation de l'information ponctuelle

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : La validation croisée est une technique utilisée pour évaluer la performance prédictive d'un modèle en le 
#               testant sur des sous-ensembles de données . Elle pemet d'évaluer la capacité de généralisation du modèle.
#               Pour effectuer la validation croisée, il est nécessaire de diviser le jeu de données en k groupes. k itérations sont 
#               effectuées où à chaque itération un groupe est utilisé comme jeu de test et les k-1 autres groupes comme jeu de calibration.
#               Cette séparation des données en groupes peut se faire de façon classique comme spatiale
#               La forme classique consiste à diviser aléatoirement le jeu de données en k groupes alors que celle 
#               spatiale divise les données en k groupes en fonction des communes donc les points utilisés pour la calibration appartiennent au 
#               même groupe de communes.
#               Dans un soucis d'évaluer l’influence de la densité d’échantillonnage sur la qualité des prédictions et la sensibilité 
#               des modèles utilisées face à une diminution de la densité d’échantillonnage du jeu de calibration, la dégradation est intégrée 
#               dans le processus de validation croisée. L'idée est de  réduire la taille de l'échantillon
#               utilisé pour la calibration du modèle. Ainsi à chaque itération un échantillon du pool de calibration (k-1 groupes) est prélevé
#               aléatoirement pour calibrer le modèle. Ce processus est repeté 30 fois pour chaque taille d'échantillon et pour chaque forme 
#               de validation croisée (classique et spatiale).


# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    21-07-2025

# Entrees :     Jeu de données sur les propriétés des sols, scripts utilitaireS pour le RF, KO et KED

# Sorties :     Tableu avec les valeurs des indicateurs de performance pour chaque modèle et selon la taille d'échantillon la forme de la CV 
#               et l4approche de CSMS

# Modification : 08-10-2025
#===========================================================================================================================#

#================================================DEBUT DU SCRIPT============================================================#



#1. Chargement des packages----

library(sf)
library(readxl)
library(tidyr)
library(dplyr)
library(foreach)
library(raster) 
library(purrr)
library(tuneRanger)
library(INLA)
library(inlabru)
library(iml)
library(mlr)
library(caret)

# setwd("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R")

#2. Fonctions RF, INLA, geomasking,  Myeval----

source("dégradation/fonction_RF.R")

source("dégradation/fonction_inla.R")

source("dégradation/fonction_geomasking.R")

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


#3. Definir les paramètres----

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

name <- args[1]  #nom de la variable cible (ex. "arg" ou "pH")

sample_sizes <- args[2] #tailles d'échantillon de calibration pour la dégradation

repets <- args[3] #nombre de répétitions par taille d’échantillon (entier)

kmax <- 30 #nombre de coeurs pour la parallelisation

ntree <- 350  #nombre d'arbres pour la RF

NomsCoord <- c("x", "y") #noms des colonnes des coordonnées

types_validation <- c("Classique", "Spatiale") #types de validation croisée

drive <- "/media/communs_infosol/" # ou "Y:/" #chemin d'accès aux données

DistanceGeomasking = 0


#4. Chargement des données---- 


com <- st_read( paste0(drive, "BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp") )

datacov <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/donnees_ponctuelles", name, ".rds"))

moyenne_covariable <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/moyenne_covariable", name, ".rds"))

cov_brt <- readRDS(paste0(drive, "BDAT/traitement_donnees/MameGadiaga/resultats/", name, "_cov_brt.rds"))

centroides_communes <- st_centroid(com) %>% dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID) %>% st_drop_geometry()



#5. Boucle de la dégradation dans le processus de validation croisée ----


cat("\n==============  Starting loops ===============\n")

bru_safe_inla(multicore = FALSE)

results_rf_all <- list()  #liste qui va accumuler les data.frames de résultats.

resLoop <- foreach( # Boucle de la dégradation sur les tailles d'échantillon n
  
  n = sample_sizes,
  .combine = rbind.data.frame ,
  .errorhandling='pass'
  
  ) %do% {
    
    n <- as.numeric(n)# test pour voir si cela corrige le pb de test degradation
    
    cat("\n============== Taille d'échantillon :", n, "===============\n")
    
    k <- 10 #  ceiling(nrow(datacov) / 1000)
    
    foreach (# Boucle sur les répétitions pour chaque taille d'échantillon n
      
      rep = 1:repets,
      .combine = rbind.data.frame,
      .errorhandling='pass'
      
             ) %do%  {
               
      cat("\n---- Répétition :", rep, "----\n")
      
      set.seed(1001 + rep)
      
      foreach (# Boucle sur les types de validation croisée
        
        type_val = types_validation,
               .combine = rbind.data.frame,
               .errorhandling='pass'
        
               ) %do% {
                 
        cat(">> Validation :", type_val, "\n")
        
        # Création des groupes ou folds en fonction du type de validation croisée
                 
        folds <- if (type_val == "Classique") {
          createFolds(datacov[[name]], k = k, list = TRUE)
        } else {
          communes <- unique(datacov$INSEE_COM)
          split(communes, sample(rep(1:k, length.out = length(communes))))
        }
        
        # Boucle sur les folds
        
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
          
           # Echantillonnage du jeu de calibration à partir du pole de calibration (k-1 groupes)  
                   
          NCalibTotal <- nrow(calib_pool)
          
          cat("-----taille du groupe de la kfold  ", NCalibTotal ,"\net n = ", n ,"et le test est ",  NCalibTotal>=n ," \n>>>>>")
                   
          if ( NCalibTotal>=n ) {
            
            calib_points <- calib_pool %>% sample_n(n)  
            
            cat("----- resample  ",n," >>>>>")
            
          }  else {
            
            cat("-----Pas de resample  !!!",n," >>>>>\n")
            
            calib_points <- calib_pool
            
          }
            
          cat("-----taille du groupe resample de la kfold  ",nrow(calib_points)," >>>>>")
          
          
          # Application des modèles RF, KO et KED selon l'appoche de CSMS (Données ponctuelles vs Désagrégation)
    
          # RF Ponctuelle
          
          cat("=> Exécution RF ; Approche : Données Ponctuelle ", type_val, "\n")
          
            ## Appel de la fonction RF
          
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
                NomsCoord)
                       
          
            ## Ajout des prédictions RF aux jeux de calibration pour KED
          
            calib_points$pred  <- res_rf_p$predCal
          
          # RF Désagrégation
          
          ## Agrgégation de la variable cible paà l'échelle communale
          
            agg_target <- calib_points %>%
              dplyr::select(INSEE_COM, all_of(name)) %>%
              group_by(INSEE_COM) %>%
              dplyr::summarise(!!name := mean(.data[[name]], na.rm = TRUE),
                      .groups = "drop")
          
              train_aggr <- moyenne_covariable %>%
                            inner_join(agg_target, 
                             by = "INSEE_COM" )
          
          cat("=> Exécution RF ; Approche : Désagrégation" , type_val, "\n")
          
          ## Appel de la fonction RF
          
            res_rf_c <- run_rf(
              "Désagrégation",
              type_val,
              train_aggr,
              test_data,
              cov_brt,
              moyenne_covariable,
              name,
              ntree,
              kmax,
              NomsCoord )
          
          ## Ajout des prédictions RF aux jeux de calibration pour KED
            
          train_aggr$pred = res_rf_c$predCal
          
          
          # KO Ponctuelle
          
          cat("=> Exécution KO ; Approche : Données ponctuelles ", type_val, "\n")
          
            ## Appel de la fonction INLA SPDE
              res_ko_p <- run_inla_spde_core(calib_points,
                                         test_data,
                                         name,
                                         "KO",
                                         "Ponctuelle",
                                         type_val)
          
              
           #KO Désagrégation
              
            ## Aggrégation de la variable cible à l'échelle communale
              
          train_aggr <- train_aggr  %>%
            inner_join(centroides_communes, by = "INSEE_COM") %>%
            dplyr::rename(x = X, y = Y) %>%
            mutate(id = INSEE_COM)
          
            ## Appel de la fonction INLA SPDE
          
          cat("=> Exécution KO ; Approche : Désagrégation ", type_val, "\n")
          
          res_ko_c <- run_inla_spde_core(as.data.frame(train_aggr),
                                         test_data,
                                         name,
                                         "KO",
                                         "Désagrégation",
                                         type_val)
          
          # KED Ponctuelle
          
          cat("=> Exécution KED ; Approche : Données ponctuelles ", type_val, "\n")
          
          ## Appel de la fonction INLA SPDE
          res_ked_p <- run_inla_spde_core(calib_points,
                                          res_rf_p$detail,
                                          name,
                                          "KED",
                                          "Ponctuelle",
                                          type_val)
          
          gc()
          
          # KED Désagrégation
          
          cat("=> Exécution KED ; Approche : Désagrégation ", type_val, "\n")
          
          ## Appel de la fonction INLA SPDE
          res_ked_c <- run_inla_spde_core(as.data.frame(train_aggr),
                                          res_rf_c$detail ,
                                          
                                          name,
                                          "KED",
                                          "Désagrégation",
                                          type_val)
          
          gc()
          
          
          
          # Stockage des prédictions
          
          
# Le bloc suivant construit un data.frame en empilant les deux approches : Ponctuelle et Désagrégation.
# Ajoute les métadonnées (type_val, sample_size, rep, fold), renomme la cible en 'obs',
# puis ne conserve que : id, approach, type_val, sample_size, rep, fold, obs, pred, predKO, predKED.
# Le tableau obtenu est ajouté à la liste 'results_rf_all' (un élément par fold)          
          
          
      results_rf_all[[length(results_rf_all) + 1]] <- bind_rows(
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
              mutate(approach = "Désagrégation",
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

pred_RF_full <- bind_rows(results_rf_all)

# 6. Save file -------------------- 

if(drive == "Y:/") {
  Myfile = paste0(drive, 
                  "BDAT/traitement_donnees/MameGadiaga/resultats/",
                  paste0("Xval_",
                         name,
                         sample_sizes,
                                collapse = "_") ,
                         ".rds")
  } else {
    
    Myfile = paste0( 
                    "output/",
                    paste0("Xval_",
                           name,
                           sample_sizes,
                           collapse = "_") ,
                    ".rds")
    
    Myfile2 = paste0( 
      "output/",
      paste0("Xval2_",
             name,
             sample_sizes,
             collapse = "_") ,
      ".rds")
    
  }
  


                       

saveRDS(pred_RF_full, 
        Myfile)
  
saveRDS(resLoop, 
        Myfile2)

