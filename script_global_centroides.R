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

library(gstlearn) # géostatistique
library(ggpubr) # pour les graphiques
library(ggnewscale)
library(doParallel)

library(raster) # copie de terra mais pour le package OGC
#library(OGC) # Pour le calcul des coordonnées obliques
library(tuneRanger)

# INLA SPDE avec inla bru, approche bayesienne de la géostatistique
library(INLA)
library(inlabru)
library(RColorBrewer) 
library(patchwork) 

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


## Définition des variables ------

name="moy_pH"
kmax= 22 # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
nbOGC = 5 # le nombre de pseudo covariables oblique

k=10 # pour la k fold cross validation

nsim=100 # for bayesian inla simulation


NomsCoord <- c("X","Y")

# 1 Preparation des données pour la spatialisation
#Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

##liste des fichiers de covariables----
l<- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)

st <- rast(l)
rast_za <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_za.tif")
plot(st)
# 
# writeRaster(st, file = "output/covariables.tiff" , overwrite = T)

# prepare covar into a table from the stack r1
gXY <- as.data.frame(st , xy=TRUE) %>%
  na.omit( )


dtTB <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/centroides_covariables.rds")


# attribution d'un id par sites (pour les doublons analytique possible)
datacov <-  dtTB %>%
  st_as_sf(coords = NomsCoord ,  crs = 2154, remove = FALSE) %>%
  mutate(id = row_number()) %>%
  na.omit()


fold = createFolds(y = datacov$id, k = k)


# 3 Modélisation par RF -----

colnames(datacov)

# spécifier les numéro des colonnes pour les covariables et la variable
idcovs = 6:69
idvar = 4

source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/RandomForest_centroides.R")

resuXvalQRF

saveRDS(resuXvalQRF, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_QRF_centroides.rds") 

# 2 krigeage ordinaire --------
# https://inlabru-org.github.io/inlabru/articles/random_fields_2d.html

# prepare sp data for inlabru

dataINLA <- datacov[,c(NomsCoord,name)]
dataINLA<-st_drop_geometry(dataINLA)

dataINLA$activ <- dataINLA[,name]

coords <- datacov[,NomsCoord]
coords<-st_drop_geometry(coords)

coordinates(dataINLA) <- NomsCoord
proj4string(dataINLA) <-  "epsg:2154"


# creation d'un tableau avec les prédictions random forest

r <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/moy_pHqrf_centroides.tif")

dataINLA$qrf <-  terra::extract(  r , vect(dataINLA)  )$QRF_Median

source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/KO_INLASPDE_centroides.R")


resuXvalTKO
saveRDS(resuXvalTKO, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_KO_centroides.rds")

# 4 Krigeage avec dérive externe -------------



source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/KED_INLASPDE_centroides.R")
resuXvalpredINLAKED
saveRDS(resuXvalpredINLAKED, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_KED_centroides.rds")

## Modelisation et spatilisation



# 5 Synthèse -----

## Validation croisée

resuXvalQRF
resuXvalTKO
resuXvalpredINLAKED

## Carte

qrf =   rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pHqrf_centroides.tif")
koINLA =   rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKOINLA_centroides.tif")
kedINLA =   rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKEDINLA_centroides.tif")

qrf = terra::resample(qrf,kedINLA)
koINLA = terra::resample(koINLA,kedINLA)

#apllication du mask

crs(rast_za) <- "EPSG:2154"
crs(qrf) <- "EPSG:2154"
crs(koINLA) <- "EPSG:2154"
crs(kedINLA) <- "EPSG:2154"


qrf<- extend(qrf, rast_za, snap = "near")
koINLA<- extend(koINLA, rast_za, snap = "near")
kedINLA<- extend(kedINLA, rast_za, snap = "near")

qrf = mask(qrf, rast_za)
koINLA = mask(koINLA, rast_za)
kedINLA = mask(kedINLA, rast_za)

predstack <- c(koINLA,qrf,kedINLA)
names(predstack) <- c("Krigeage Ordi. INLA","QRF","KED-INLA")
plot(predstack)

# tm_shape(predstack) +
#   tm_raster(
#     col.scale = tm_scale(values = terrain.colors(10) ,
#                          style = "quantile",n = 10),
#     col.legend = 
#       tm_legend(
#         position = c(0,0.3),
#         item.height = .6,
#         item.width = .5,
#         item.r = 0, 
#         text.size = .3,
#         item.space = 0.05, 
#         item.na.space = .51, 
#         title.align = "pH")
#   )

# tmap_mode("view")
# tm_shape(predstack[[3]]) + tm_raster(style="quantile" , n=8)
# tmap_mode("plot")


