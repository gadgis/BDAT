#=============================================================================================================================#
# Script :      Script d'extraction des covariables 

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Ce script permet d'extraire, depuis un stack de raster de covariables les valeurs de ces dernières 
#               au niveau des observations. Il calcule également des statistiques zonalees par commune à partir du 
#               stack raster :
#               - moyenne pour les couches numériques,
#               - mode (valeur la plus fréquente) pour les couches catégorielles.

#               Les résultats sont sauvegardés en deux fichiers .rds :
#               - données ponctuelles enrichies par les covariables,
#               - tableau des covariables agrégées par commune.

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    21-07-2025

# Entrees :     Shapefile des communes, observations ponctuelles, les covariables

# Sorties :     data.frame avec les valeurs des covariables pour les données ponctuelles et un data.frame
#               avec la moyenne des covariables par commune.

# Modification : 09-10-2025
#===========================================================================================================================#

#================================================DEBUT DU SCRIPT============================================================#


#Fonction pour extraire la moyenne des covariables par commune


zonal_covariates_by_commune <- function(st, com, cat_vars, num_vars,
                                        id_col = "INSEE_COM", exact = FALSE) {
  # Harmoniser CRS et passer en SpatVector
  com_aligned <- sf::st_transform(com, crs = terra::crs(st))
  vcom <- terra::vect(com_aligned)
  
  if (!id_col %in% names(vcom)) {
    stop(sprintf("La colonne '%s' n'existe pas dans l'objet des communes.", id_col))
  }
  
  # Ne garder que les couches présentes dans le stack
  all_layers <- names(st)
  cat_layers <- intersect(cat_vars, all_layers)
  num_layers <- intersect(num_vars, all_layers)
  
  # Moyennes pour les couches numériques
  res_num <- if (length(num_layers)) {
    terra::extract(
      st[[num_layers]], vcom,
      fun = mean, na.rm = TRUE, exact = exact, bind = TRUE
    ) |>
      as.data.frame() |>
      dplyr::select(dplyr::all_of(c(id_col, num_layers)))
  } else NULL
  
  # Modes pour les couches catégorielles
  res_cat <- if (length(cat_layers)) {
    terra::extract(
      st[[cat_layers]], vcom,
      fun = terra::modal, na.rm = TRUE, exact = exact, bind = TRUE
    ) |>
      as.data.frame() |>
      dplyr::select(dplyr::all_of(c(id_col, cat_layers)))
  } else NULL
  
  #  Fusion
  datacov_mean <- if (!is.null(res_num) && !is.null(res_cat)) {
    dplyr::left_join(res_num, res_cat, by = id_col)
  } else if (!is.null(res_num)) {
    res_num
  } else if (!is.null(res_cat)) {
    res_cat
  } else {
    stop("Aucune couche de 'cat_vars' ou 'num_vars' trouvée dans le stack.")
  }
  
  #  Remettre les catégorielles en facteur
  if (length(cat_layers)) {
    datacov_mean <- datacov_mean |>
      dplyr::mutate(dplyr::across(dplyr::all_of(cat_layers), as.factor))
  }
  
  datacov_mean
}



#1. Chargement des packages----
library(sf)
library(tidyr)
library(dplyr)
library(terra)
library(sp)
library(raster)
library(DescTools)

#2.Définir les paramètres----

name <- "pH"

NomsCoord <- c("x", "y")


#3.Chargement des données et des covariables----

com<-st_read("Y:/BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp")

dtTB<-readRDS(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_", name, ".rds"))


#4. Stack des covariables et rxtraction des centroides----

chemin_cov <- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall" #chemin pour les covariables

st <- rast(list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)) #stack des covariables

centroides_communes<-st_centroid(com)%>%
  dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID)%>%
  st_drop_geometry()  #extraction des centroides


#5. Extraction des covariables sur les points----

datacov <- terra::extract(st, st_as_sf(dtTB, coords = NomsCoord, crs = 2154)) %>%
  bind_cols(dtTB %>% 
              dplyr::select(all_of(c(name, NomsCoord, "source", "INSEE_COM")))) %>%
  mutate(id = row_number()) %>%
  na.omit()

datacov <- datacov %>%
  mutate(INSEE_COM = as.character(INSEE_COM))


#6. Définir les types des covariables----


all_vars<-colnames(datacov[2:65])

cat_vars <- c("bdforet0", "bdforet1", "bdforet2", "bdforet3", "bdforet30", "bdforet4", "bdforet40", "bdforet5",
                  "bdforet50", "bdforet6", "Geology1", "Geology10", "Geology11", "Geology12", "Geology13", "Geology14",
                  "Geology2", "Geology3", "Geology4", "Geology5", "Geology6", "Geology7", "Geology8", "Geology9",
                  "OCS201611", "OCS201612", "OCS2016211", "OCS2016221", "OCS2016222", "OCS201631", "OCS201632",
                  "OCS201634", "OCS201636", "OCS201641", "OCS201642", "OCS201643", "OCS201644", "OCS201645",
                  "OCS201646", "OCS201651", "typo2", "typo3", "typo4", "typo5") #variables catégorielles

num_vars <- setdiff(all_vars, cat_vars) #variables numériques

#transformations en facteur des variables catégorielles

datacov <- datacov %>%
  mutate(across(all_of(cat_vars), as.factor))

#7. Moyenne des covariables par commune----

datacov_mean <- zonal_covariates_by_commune(
  st = st,
  com = com,
  cat_vars = cat_vars,
  num_vars = num_vars,
  id_col = "INSEE_COM",
  exact = FALSE 
)


#sauvegarde des résultats

saveRDS(datacov, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/donnees_ponctuelles", name, ".rds"))
saveRDS(datacov_mean, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/moyenne_covariable", name, ".rds"))

#=================================================FIN DU SCRIPT==============================================================#