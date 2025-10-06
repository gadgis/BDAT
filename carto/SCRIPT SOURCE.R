#===========================================================================================================#
# Script : Script principal 

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la cartographie des propriétés des sols en utilisant la méthode Random Forest 
#               Krigeage ordinaire et Krigeage avec Dérive Externe 

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Données ponctuelles de la BDAT et IGCS sur les propriétés des sols

# Sorties :     Carte des propriétés des sols sur la zone agricole

# Modification : 06-10-2025
#===========================================================================================================#


#DEBUT DU SCRIPT------------------------------------------------------------------------------------------

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


# 1. Définition des variables ------

name="arg" # changer la variable d'intérêt au besoin 
kmax= 23 # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
k=10 # pour la validation croisée
nsim=100 # for bayesian inla simulation
NomsCoord <- c("x","y")

#2.Preparation des données pour la spatialisation----
##2.1. Importation des données----
datacov <- readRDS(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/donnees_ponctuelles", name, ".rds")) # jeu de données ponctuelles avec les covariables

moyenne_covariable <- readRDS(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/moyenne_covariable", name, ".rds")) # jeu de données avec les moyennes des covariables par commune

com <- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp") # shapefile des communes

rast_za <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_za.tif") # raster de la zone agricole

df_vars <- read.csv("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/COV_RF.csv", sep = ";", stringsAsFactors = FALSE  ) # Nommenclature des covariables

##2.2. Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

#liste des fichiers de covariables
l<- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)

#stack des covariables
st <- rast(l)
plot(st)


# Préparer les covariables sous forme de tableau à partir du stack
gXY <- as.data.frame(st , xy=TRUE) %>%
  na.omit( )

centroides_communes <- st_centroid(com) %>% dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID) %>% st_drop_geometry()

#Enlever le commentaire et exécuter les lignes suivantes si vous vouler utilisiser l'approche centroides
# Changer la variable d'intérêt au besoin

# y_agg<-datacov %>%
#   group_by(INSEE_COM) %>%
#   summarise(arg = mean(arg, na.rm = TRUE)) %>%
#   ungroup()
# 
# datacov<-moyenne_covariable %>%
#   left_join(y_agg, by = "INSEE_COM") %>%
#   left_join(centroides_communes, by = "INSEE_COM")%>%
#   rename(x=X, y=Y)
# 
# 
# datacov$id<- 1:nrow(datacov)

fold = createFolds(y = datacov$id, k = k)

# 3 Modélisation par Random Forest -----

colnames(datacov)

# spécifier les numéro des colonnes pour les covariables et la variable
idcovs = 2:65
idvar = 66

source("carto/RandomForest.R")

resuXvalQRF


# 4. Krigeage Ordinaire -----

#https://inlabru-org.github.io/inlabru/articles/random_fields_2d.html

##4.1.  préparation du fichier sp pour inlabru----

dataINLA <- datacov[,c(NomsCoord,name,"predRF")]
dataINLA$activ <- dataINLA[,name]
coords <- datacov[,NomsCoord]

coordinates(dataINLA) <- NomsCoord
proj4string(dataINLA) <-  "epsg:2154"


## 4.2. Création d'un tableau avec les prédictions random forest----

r <- rast(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "qrf.tif"))


dataINLA$qrf <-  terra::extract(  r , vect(dataINLA)  )$QRF_Median

source("carto/KO_INLASPDE.R")

resuXvalTKO


# 5. Krigeage avec dérive externe -------------


source("carto//KED_INLASPDE.R")
resuXvalpredINLAKED

# 6. Cartographie des résultats-----

# Reéchantillonnage des rasters 

qrf =   rast(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"qrf.tif"))
koINLA =   rast(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"predKOINLA.tif"))
kedINLA =   rast(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"predKEDINLA.tif"))

qrf = terra::resample(qrf,kedINLA)
koINLA = terra::resample(koINLA,kedINLA)

#Affection du même système de projection

crs(rast_za) <- "EPSG:2154"
crs(qrf) <- "EPSG:2154"
crs(koINLA) <- "EPSG:2154"
crs(kedINLA) <- "EPSG:2154"

# Extension et masquage des rasters à la zone agricole
qrf<- extend(qrf, rast_za, snap = "near")
koINLA<- extend(koINLA, rast_za, snap = "near")
kedINLA<- extend(kedINLA, rast_za, snap = "near")

qrf_agri = mask(qrf, rast_za)
koINLA_agri = mask(koINLA, rast_za)
kedINLA_agri = mask(kedINLA, rast_za)

# 10. Sauvegarde des rasters finaux-----

writeRaster(qrf_agri, file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"qrf_final_cent.tif"), overwrite = TRUE)
writeRaster(koINLA_agri, file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"predKOINLA_final_cent.tif"), overwrite = TRUE)
writeRaster(kedINLA_agri, file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"predKEDINLA_final_cent.tif"), overwrite = TRUE)

# 11. Visualisation des résultats-----
predstack <- c(koINLA_agri,qrf_agri,kedINLA_agri)
names(predstack) <- c("Krigeage Ordi. INLA","QRF","KED-INLA")
plot(predstack)

tm_shape(predstack) +
  tm_raster(
    col.scale = tm_scale(values = terrain.colors(10) ,
                         style = "quantile",n = 10),
    col.legend = 
      tm_legend(
        position = c(0,0.3),
        item.height = .6,
        item.width = .5,
        item.r = 0, 
        text.size = .3,
        item.space = 0.05, 
        item.na.space = .51, 
        title.align = "Carbone")
  )

#FIN DU SCRIPT------------------------------------------------------------------------------------------