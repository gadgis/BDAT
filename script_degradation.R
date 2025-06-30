library(sf) # manipuler les datas SIG vector
library(tmap)  # outil d'édution de carte
library(readxl) # lecture de fichier excel
library(tidyr) # reroganisation des données
library(dplyr)  # grammaire de manipulation de données
library(terra) # manipuler des données raster-pixels
library(ggplot2) # graphiques
library(sp) # ancetre de sf. Pourquoi l'utiliser ?
library(randomForest) # randomForest pour utiliser tuneRandomForest
library(Boruta)  # implementation  de selection de variables de randomForest
library(ranger) # quantile randomForest, version optimée de randomForest
library(caret) # modelisation creation de modeles prédictifs avec tuning validation-croisée
library(foreach) # boucles optimisées de R plus facile à coder voir option combine
library(raster) # manipuler des données raster
library(purrr)
library(ggpubr) # pour les graphiques
library(ggnewscale)
library(doParallel)
library(raster) # copie de terra mais pour le package OGC
library(tuneRanger)
# INLA SPDE avec inla bru, approche bayesienne de la géostatistique
library(INLA)
library(inlabru)
library(iml)
library(mlr)


## Fonction utilsée pour la validation --------------

Myeval <- function(x, y){
  
  # mean error
  ME <- round(mean(y - x, na.rm = TRUE), digits = 2)
  
  # root mean square error
  RMSE <-   round(sqrt(mean((y - x)^2, na.rm = TRUE)), digits = 2)
  
  # root mean absolute error
  MAE <-   round(mean(abs(y - x), na.rm = TRUE), digits = 2)
  
  # Pearson's correlation squared
  r2 <-  round((cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2), digits = 2)
  
  # MEC
  SSE <- sum((y - x) ^ 2, na.rm = T)
  SST <- sum((y - mean(y, na.rm = T)) ^ 2, na.rm = T)
  NSE <- round((1 - SSE/SST), digits = 2)
  
  # concordance correlation coefficient
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
  rhoC <- round(rCb, digits = 2)
  
  Cb <- round(Cb, digits = 2)
  r <- round(r, digits = 2)
  
  # return the results
  evalRes <- data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE, rhoC = rhoC, Cb = Cb)
  
  return(evalRes)
}

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R")

## Définition des variables ------

name="pH"
kmax= 23 # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
k=10 # pour la k fold cross validation
nsim=100 # for bayesian inla simulation
NomsCoord <- c("x","y")
sample_sizes <- c(600,800,1000,1200,1400,1600,1800,2000,3000,4000,5000,6000,7000,8000,9000,10000) 
repets <- 30

# 1 Preparation des données pour la spatialisation
#Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

##liste des fichiers de covariables----
l<- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)

st <- rast(l)


dtTB <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")

# Extraction covariables
datacov <- terra::extract( st , 
                           dtTB %>%
                             st_as_sf(coords = NomsCoord ,
                                      crs = 2154)
                           ) %>% 
  bind_cols(dtTB %>%
              dplyr::select(all_of(c( name,NomsCoord,"source")  )
                            )
            ) %>%

  mutate(id = row_number()) %>%
  na.omit()

# spécifier les numéro des colonnes pour les covariables et la variable
idcovs = 2:65
idvar = 66

# Initialisation du stockage
results_all <- list()

#Boucle principale
for (n in sample_sizes) {
  cat(">> Taille d'échantillon :", n, "\n")
  k <- ifelse(n <= 2000, ceiling(n / 200), ceiling(n / 500))
  
  for (rep in 1:repets) {
    cat("  → Répétition :", rep, "\n")
    set.seed(rep)
    
    # Sous-échantillonnage
    data_sample <- datacov %>% dplyr::sample_n(n)
    data_sample$id <- 1:nrow(data_sample)
    
    # Création des folds
    fold <- split(data_sample$id, rep(1:k, length.out = n))
    
    # Random Forest 
    source("Random_Forest_D.R")  
    
    # Préparation unique de dataINLA pour les deux krigeages
    dataINLA <- data_sample[, c(NomsCoord, name,"predRF")]
    dataINLA$activ <- dataINLA[[name]]
     
  
    coordinates(dataINLA) <- NomsCoord
    proj4string(dataINLA) <- CRS("epsg:2154")
    
    coords <- data_sample[, NomsCoord]  # pour le maillage
    
    # INLA KO
    source("KO_INLASPDE_D.R") 
    
    # INLA KED
    source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/KED_INLASPDE_D.R")  
    
    # Ajout des métadonnées
    res_rf$method <- "RF"
    res_ko$method <- "INLA_KO"
    res_ked$method <- "INLA_KED"
    
    # Compilation
    res_all <- bind_rows(res_rf, res_ko, res_ked) %>%
      mutate(rep = rep, sample_size = n)
    
    results_all[[paste0(n, "_", rep)]] <- res_all
  }
}

results_final <- bind_rows(results_all)

results_summary <- results_final %>%
  group_by(sample_size, method) %>%
  summarise(across(c(ME, MAE, RMSE, r, r2, NSE, rhoC, Cb), mean, na.rm = TRUE), .groups = "drop")

results_summary<-results_summary %>%
  mutate(across(c(ME, MAE, RMSE, r, r2, NSE, rhoC, Cb), round, digits = 2))

saveRDS(results_summary, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/resultats_degradation_pH.rds")


