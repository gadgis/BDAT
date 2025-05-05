#Chargement des packages----
library(sf)
library(terra)
library(dplyr)
library(randomForest)
library(ranger)
library(exactextractr)
library(tmap)

#Chargement des données----
donnnees <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")
ocsol<-vect("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/OCCUPATION_SOL.shp")
communes <- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\communes53_2014.shp")
ph_moyen<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/ph_moyen.rds")
dept<- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\dept53_rep.shp")

#Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

##liste des fichiers de covariables----
fichiers_cov <- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)
print(fichiers_cov)

##Création d'une liste de rasters----
rasters_cov <- lapply(fichiers_cov, function(x) {
  rast(x)
})
print(rasters_cov)

#définir le type de covariables
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


##Empilage des rasters----
rst <- do.call(c, rasters_cov)
# Séparer les rasters continus et catégoriels
rasters_continuous <- subset(rst, which(types_covariables == "continuous"))
rasters_categorical <- subset(rst, which(types_covariables == "categorical"))
print(rasters_continuous)
print(rasters_categorical)


##Définir un masque par la zone agricole----
zone_agricole <- ocsol[ocsol$CODE_US == "US1.1", ]

# Rasteriser cette zone agricole
rast_za <- rasterize(zone_agricole, rst, field=1, background=NA)
print(rast_za)

## Appliction du masque au stack----
# rst_za <- mask(rst, rast_za)
# print(rst_za)


# conversion en df du stack en zone agricole
gXY <- as.data.frame(rst, xy = TRUE, na.rm = TRUE) %>%  
  mutate(across(
    .cols = names(rasters_categorical),  
    .fns = as.factor                    
  ))

str(gXY)


##Extraction des valeurs de covariables pour les points----
donnees_pH<-st_as_sf(donnnees, coords = c("x", "y"), crs = "2154")

matrice_cov <- terra::extract(rst, donnees_pH)%>%  
  mutate(across(
    .cols = names(rasters_categorical),  
    .fns = as.factor                    
  ))

summary(matrice_cov)

#Joindre les attributs des points avec les covariables extraites

donnees_extraits <- cbind(as.data.frame(donnees_pH), matrice_cov)

summary(donnees_extraits)

#Sauvegarde de la matrice de covariables
saveRDS(donnees_extraits, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/matrice_covariables.rds")

#Extraction des matrices de covariables pour les centroides----

##Création des centroides des communes----
ph_moyen <- ph_moyen %>%
  mutate(INSEE_COM = as.character(INSEE_COM)) 

##Création des centroïdes avec gestion des types
centroides <- communes %>%
  st_centroid() %>%
  select(INSEE_COM) %>%
  mutate(INSEE_COM = as.character(INSEE_COM)) %>%  
  left_join(ph_moyen, by = "INSEE_COM") %>%
  mutate(
    X = st_coordinates(.)[,1],  
    Y = st_coordinates(.)[,2]   
  ) %>%
  select(INSEE_COM, X, Y, moy_ph)%>%
  st_drop_geometry()

##Extraction des valeurs----

# Fonction mode (valeur la plus fréquente)
mode_function <- function(values, coverage_fraction) {
  values <- na.omit(values)
  if (length(values) == 0) {
    return(NA)
  } else {
    return(as.numeric(names(which.max(table(values)))))
  }
}

# Extraire les moyennes des variables continues
moyennes_cov <- exact_extract(rasters_continuous, communes,
                              fun = function(values, coverage_fraction) mean(values, na.rm = TRUE),
                              stack_apply = TRUE)

names(moyennes_cov) <- gsub("^fun\\.", "", names(moyennes_cov))

# Extraire les modes des variables catégorielles
mode_cov <- exact_extract(rasters_categorical, communes,
                              fun = mode_function,
                              stack_apply = TRUE)
names(mode_cov) <- gsub("^fun\\.", "", names(mode_cov))

#Fusionner les résultats
df_covariables <- cbind(moyennes_cov, mode_cov)
df_covariables <- df_covariables %>%
  mutate(INSEE_COM = communes$INSEE_COM)

#jointure avec les centroides
centroides <- centroides %>%
  left_join(df_covariables, by = "INSEE_COM")

saveRDS(centroides, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/centroides_covariables.rds")

