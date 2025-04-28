#Chargement des packages----
library(sf)
library(terra)
library(dplyr)
library(randomForest)
library(ranger)

#Chargement des données----
donnnees <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")
ocsol<- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/OCCUPATION_SOL.shp")

#Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

##liste des fichiers de covariables----
fichiers_cov <- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)
print(fichiers_cov)

##Création d'une liste de rasters----
rasters_cov <- lapply(fichiers_cov, function(x) {
  rast(x)
})

##Spécification du type de variable pour le re-échantillonnage----
types_covariables <- c(
  # bdforet 0-6 et 30-50 : CATEGORICAL
  rep("categorical", 10),
  
  # Clay_eu23 : CONTINUOUS
  "continuous",
  
  # cti : CONTINUOUS
  "continuous",
  
  # curvature : CONTINUOUS
  "continuous",
  
  # elevation_BDalti : CONTINUOUS
  "continuous",
  
  # exposition : CONTINUOUS
  "continuous",
  
  # Geology1 - Geology14 : CATEGORICAL
  rep("categorical", 14),
  
  # idpr : CONTINUOUS
  "continuous",
  
  # Kgamma : CONTINUOUS
  "continuous",
  
  # OCS201611 à OCS201651 : CATEGORICAL
  rep("categorical", 16),
  
  # PC1_NDVI, PC2_NDVI, PC3_NDVI : CONTINUOUS
  rep("continuous", 3),
  
  # roughness : CONTINUOUS
  "continuous",
  
  # Sand_eu23 : CONTINUOUS
  "continuous",
  
  # scale_pos : CONTINUOUS
  "continuous",
  
  # Silt_eu23 : CONTINUOUS
  "continuous",
  
  # slope : CONTINUOUS
  "continuous",
  
  # slopeascos : CONTINUOUS
  "continuous",
  
  # Thgamma : CONTINUOUS
  "continuous",
  
  # ThK_ratio : CONTINUOUS
  "continuous",
  
  # ThU_ratio : CONTINUOUS
  "continuous",
  
  # typo2 - typo5 : CATEGORICAL
  rep("categorical", 4),
  
  # Ugamma : CONTINUOUS
  "continuous"
)

##Définir une projection cible pour les raster
projection_cible <- st_crs(ocsol)
projection_cible_txt <- as.character(crs(ocsol))

##Reprojetion des rasters----
rasters_cov <- lapply(rasters_cov, function(r) {
  if (as.character(crs(r)) != projection_cible_txt) {
    project(r, projection_cible_txt)
  } else {
    r
  }
})

# # Tous les CRS sont-ils identiques ?
# all(sapply(rasters_cov, function(r) as.character(crs(r)) == as.character(crs(ocsol))))

##Re-échantillonnage des rasters----

#Création d'un modèle de raster à 90 m
template_90m <- rast(ext(ocsol),
                     resolution = 90,
                     crs = projection_cible_txt)

#re-échantillonnage au pas de 90 m

rasters_cov90m <- mapply(function(r, type) {
  methode <- ifelse(type == "continuous", "bilinear", "near")
  resample(r, template_90m, method = methode)
}, rasters_cov, types_covariables, SIMPLIFY = FALSE)

# Vérification de la résolution
# Afficher la résolution de tous les rasters
# lapply(rasters_cov90m, res)

