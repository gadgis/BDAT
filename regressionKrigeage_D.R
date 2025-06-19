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

library(gstlearn) # géostatistique
library(ggpubr) # pour les graphiques
library(ggnewscale)
library(doParallel)

library(raster) # copie de terra mais pour le package OGC
library(tuneRanger)

# INLA SPDE avec inla bru, approche bayesienne de la géostatistique
library(INLA)
library(inlabru)

#Fonction de cal des indicarteurs
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

name="pH"
kmax= detectCores()-1  # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
nbOGC = 5 # le nombre de pseudo covariables oblique


nsim=100 # for bayesian inla simulation


NomsCoord <- c("x","y")

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R")

#Chargement des données

# Preparation des données pour la spatialisation
#Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

##liste des fichiers de covariables----
l<- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)

st <- rast(l)
rast_za <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_za.tif")
plot(st)
rmqs<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/RMQS_pH.rds")

# prepare covar into a table from the stack r1
gXY <- as.data.frame(st , xy=TRUE) %>%
  na.omit( )


dtTB <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")


# attribution d'un id par sites (pour les doublons analytique possible)
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

colnames(datacov)

# spécifier les numéro des colonnes pour les covariables et la variable
idcovs = 2:65
idvar = 66

# Paramètres RF
X <- datacov[, idcovs]
Y <- datacov[, idvar]

# Chargement des covariables sélectionnées (Boruta)
cov_brt <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_cov_brt.rds")
cov_brt <- c(name, cov_brt)

# Paramètres globaux
repets <- 3
sample_sizes <- c(600, 800)

results_all <- list()

source("Random_Forest_D.R")

write.csv(results_df, "resultats_degradation_validation.csv", row.names = FALSE)
