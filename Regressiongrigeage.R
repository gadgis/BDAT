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
library(OGC) # Pour le calcul des coordonnées obliques
library(tuneRanger)

# INLA SPDE avec inla bru, approche bayesienne de la géostatistique
library(INLA)
library(inlabru)



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

name="pH"
kmax= 22 # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
nbOGC = 5 # le nombre de pseudo covariables oblique

k=10 # pour la k fold cross validation

nsim=100 # for bayesian inla simulation


NomsCoord <- c("x","y")

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



fold = createFolds(y = datacov$id, k = k)


# 3 Modélisation par RF -----

colnames(datacov)

# spécifier les numéro des colonnes pour les covariables et la variable
idcovs = 2:65
idvar = 66

source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/RandomForest.R")

resuXvalQRF
saveRDS(resuXvalQRF, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_qrf.rds")

# 2 krigeage ordinaire --------
# https://inlabru-org.github.io/inlabru/articles/random_fields_2d.html

# prepare sp data for inlabru

dataINLA <- datacov[,c(NomsCoord,name)]
dataINLA$activ <- dataINLA[,name]
coords <- datacov[,NomsCoord]

coordinates(dataINLA) <- NomsCoord
proj4string(dataINLA) <-  "epsg:2154"


# creation d'un tableau avec les prédictions random forest

r <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pHqrf.tif")

dataINLA$qrf <-  terra::extract(  r , vect(dataINLA)  )$QRF_Median

source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/KO_INLASPDE.R")

resuXvalTKO
saveRDS(resuXvalTKO, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_KO.rds")


# 4 Krigeage avec dérive externe -------------


source("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R/KED_INLASPDE.R")
resuXvalpredINLAKED
resuXvalpredINLAKEDTotal

saveRDS(resuXvalpredINLAKED, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_KED.rds")
saveRDS(resuXvalpredINLAKEDTotal, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrique_KED_total.rds")

datacov <- datacov %>%
  mutate(predRF=round(predRF,1),
         predINLAKO=round(predINLAKO,1),
         predINLAKED=round(predINLAKED,1),
         predINLAKEDTotal=round(predINLAKEDTotal,1))
  
saveRDS(datacov, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/results_pts.rds") 
  
  
  
## Modelisation et spatilisation



# 5 Synthèse -----

## Validation croisée

resuXvalQRF
resuXvalTKO
resuXvalpredINLAKED
resuXvalpredINLAKEDTotal

## Carte

qrf =   rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pHqrf.tif")
koINLA =   rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKOINLA.tif")
kedINLA =   rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKEDINLA.tif")

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

qrf_agri = mask(qrf, rast_za)
koINLA_agri = mask(koINLA, rast_za)
kedINLA_agri = mask(kedINLA, rast_za)

writeRaster(qrf_agri, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pHqrf_final.tif", overwrite = T)
writeRaster(koINLA_agri, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKOINLA_final.tif", overwrite = T)
writeRaster(kedINLA_agri, file = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKEDINLA_final.tif", overwrite = T)

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

# tmap_mode("view")
# tm_shape(predstack[[3]]) + tm_raster(style="quantile" , n=8)
# tmap_mode("plot")



