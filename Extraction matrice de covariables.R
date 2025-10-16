
install.packages("DescTools")

library(sf)
library(tidyr)
library(dplyr)
library(terra)
library(sp)
library(raster)
library(DescTools)

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

name <- "pH"
NomsCoord <- c("x", "y")

#1.Chargement des données et des covariables----

chemin_cov <- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"

st <- rast(list.files(chemin_cov, pattern = ".tif$", full.names = TRUE))

dtTB<-readRDS(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_", name, ".rds"))

com<-st_read("Y:/BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp") 

centroides_communes<-st_centroid(com)%>%
  dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID)%>%
  st_drop_geometry()


#2. Extraction des covariables sur les points----

datacov <- terra::extract(st, st_as_sf(dtTB, coords = NomsCoord, crs = 2154)) %>%
  bind_cols(dtTB %>% 
              dplyr::select(all_of(c(name, NomsCoord, "source", "INSEE_COM")))) %>%
  mutate(id = row_number()) %>%
  na.omit()

datacov <- datacov %>%
  mutate(INSEE_COM = as.character(INSEE_COM))


#3. Définition des types----

#Toutes les variables catégorielles
all_vars<-colnames(datacov[2:65])

cat_vars <- c("bdforet0", "bdforet1", "bdforet2", "bdforet3", "bdforet30", "bdforet4", "bdforet40", "bdforet5",
                  "bdforet50", "bdforet6", "Geology1", "Geology10", "Geology11", "Geology12", "Geology13", "Geology14",
                  "Geology2", "Geology3", "Geology4", "Geology5", "Geology6", "Geology7", "Geology8", "Geology9",
                  "OCS201611", "OCS201612", "OCS2016211", "OCS2016221", "OCS2016222", "OCS201631", "OCS201632",
                  "OCS201634", "OCS201636", "OCS201641", "OCS201642", "OCS201643", "OCS201644", "OCS201645",
                  "OCS201646", "OCS201651", "typo2", "typo3", "typo4", "typo5")

num_vars <- setdiff(all_vars, cat_vars)

#transformations en facteur
datacov <- datacov %>%
  mutate(across(all_of(cat_vars), as.factor))

#4. Moyenne des covariables par commune----

#selection des covariables et de l'INSEE_COM
datacov_mean <- datacov %>%
  dplyr::select(all_vars, INSEE_COM)

#calcul de la moyenne (var_num) et du mode (cat_vars) par commune
datacov_mean[cat_vars] <- lapply(datacov_mean[cat_vars], as.character)

# calculer les moyennes et les modes
datacov_mean <- datacov_mean %>%
  group_by(INSEE_COM) %>%
  summarise(
    across(all_of(num_vars), ~mean(.x, na.rm = TRUE)),
    across(all_of(cat_vars), ~as.character(names(which.max(table(.))))),
    .groups = "drop"
  )
# Convertir les colonnes catégorielles en facteurs
datacov_mean <- datacov_mean %>%
  mutate(across(all_of(cat_vars), as.factor))

#sauvegarde des résultats

saveRDS(datacov, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/donnees_ponctuelles", name, ".rds"))
saveRDS(datacov_mean, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/moyenne_covariable", name, ".rds"))
