#===========================================================================================================#
# Script : Script prétraitement et extraction des données pour la propriété cible

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description :
# Ce scrpt prépare automatiquement une base d’apprentissage pour une propriété de sol cible (`name`)
# en combinant BDAT et IGCS sur la Mayenne (53). Le script :
# - nettoie/reprojette les coordonnées, vérifie les codes INSEE et garde les points agricoles (OCS-GE) ;
# - convertit les mesures BDAT et récupère la variable cible via un mapping générique ;
# - harmonise IGCS à 0–30 cm (cas 1/2 horizons : moyenne pondérée ; ≥3 horizons : spline mpspline2) ;
# - exclut les profils composites RMQS ;
# - fusionne IGCS + BDAT, supprime les valeurs manquantes et exporte : BDAT ciblé, RMQS,
#   et la base finale `igcs_bdat_<name>.rds`.

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Shapefile (OCS_GE, communes), observations ponctuelles (BDAT et IGCS),

# Sorties :     base de données final IGCS et BDAT qui servira pour extraire les covariables 
#               et entrainer les modèles de prédiction

# Modification : 10-10-2025
#===========================================================================================================#

#================================================DEBUT DU SCRIPT====================================================#


#1.Chargement de package ----
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(tibble)
library(mpspline2)
library(rlang)

# Paramètre : variable cible 

name <- "pH"   # ex. "pH","CEC","C","clay","sand","silt"

# Correspondance BDAT -> nom standard
map_bdat <- c(
  pH   = "pho",
  CEC  = "cecmet",
  C    = "corgox",
  arg = "argi",
  sand = "sabt",
  silt = "limt"
)
stopifnot(name %in% names(map_bdat))

#2.Importation des couches ----

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees")

BDAT <- read.csv("bdat_53_x_y.csv", sep=",", header=TRUE)
igcs <- st_read("DonneesIGCSStage.shp")
dpt  <- st_read("dept_53.shp")
com  <- st_read("COMMUNE.SHP")
ocsol <- st_read("OCCUPATION_SOL.shp")
zone_urb <- st_read("ZONE_CONSTRUITE.shp")
nomenclature <- read.csv("NomenclatureOCSGE.csv", sep=";", header=TRUE)

#3.Pretraitements sur la BDAT ----

##3.1 Correction de la base ----

### 3.1.1. Reprojection et selection des points se trouvant en Mayenne----

#### Selection des colonnes d'intéret
bdat <- BDAT %>%
  dplyr::select(annee,insee,x,y,x_commune,y_commune,bdatid)

#### Suppression des lignes avec un X=NULL ou un Y=NULL
bdat <- bdat %>%
  filter(x!="NULL", y!="NULL")

#### Identification des différents SCR 
bdat_L93  <- bdat %>% filter(y>6000000)                     # Lambert 93
bdat_L2E  <- bdat %>% filter(y>1600000 & y<2700000)         # Lambert 2 Etendu
bdat_wgs84<- bdat %>% filter(y>40 & y<50)                   # WGS84
bdat_invrs<- bdat %>% filter(y>=-2.45 & y<=7.23) %>%        # coords inversées lon/lat
  rename(y=x, x=y)

#### Transformation en sf et reprojection en EPSG:2154 
bdat_L93_sf   <- st_as_sf(bdat_L93,   coords = c("x","y"), crs=2154)
bdat_L2E_sf   <- st_as_sf(bdat_L2E,   coords = c("x","y"), crs=27572) %>% st_transform(crs=2154)
bdat_wgs84_sf <- st_as_sf(bdat_wgs84, coords = c("x","y"), crs=4326)  %>% st_transform(crs=2154)
bdat_invrs_sf <- st_as_sf(bdat_invrs, coords = c("x","y"), crs=4326)  %>% st_transform(crs=2154)

#### Fusion des sous bases 
bdat_sf <- rbind(bdat_wgs84_sf,bdat_L93_sf,bdat_invrs_sf,bdat_L2E_sf)

#### Extraction des points se trouvant dans le département de la Mayenne
bdat_53 <- st_filter(bdat_sf, dpt, .predicate=st_within)

### 3.1.2. Verification de l'exactitude des codes INSEE----

# Jointure spatiale bdat et communes (codes INSEE conformes lorsque verif=0)
bdat_join <- bdat_53 %>%
  st_join(com, join = st_intersects) %>%
  mutate(INSEE_COM=as.numeric(INSEE_COM),
         x_commune=as.numeric(x_commune),
         y_commune=as.numeric(y_commune),
         verif= (insee-INSEE_COM))

coords <- st_coordinates(bdat_join)
bdat_join <- bdat_join %>%
  mutate(X=coords[,1],
         Y=coords[,2])

# Deux lots de points (1: INSEE conforme ; 2: non conforme)
bdat_INSEE_OK  <- bdat_join %>% filter(verif==0) %>% mutate(codification = 1)
bdat_INSEE_ncf <- bdat_join %>% filter(verif!=0) %>% mutate(codification = 2)

# Distance minimale à la commune attendue pour les non conformes (<= 10 km gardés)
bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  rowwise() %>%
  mutate(
    min_dist = min(as.numeric(
      st_distance(
        geometry,
        st_boundary(st_geometry(com[com$INSEE_COM == insee, ]))
      )
    ))
  ) %>%
  mutate(min_dist = min_dist / 1000) %>%     # km
  filter(min_dist <= 10) %>%
  ungroup() %>%
  dplyr::select(-min_dist)

# Fusion avec les conformes
bdat_ok <- rbind(bdat_INSEE_OK, bdat_INSEE_ncf) %>%
  dplyr::select(annee,insee,X,Y,bdatid,INSEE_COM,codification, geometry)

### 3.1.3. Identification des points en zone agricole ----
bdat_ocsol_sf <- st_join(bdat_ok, ocsol, join = st_intersects)

bdat_agri_sf <- bdat_ocsol_sf %>%
  filter(CODE_US=="US1.1") %>%
  dplyr::select(annee, insee, X, Y, bdatid, INSEE_COM, CODE_US, codification, geometry)

## 3.2. Remobilisation du tableau avec les propriétés des sols (BDAT) ----
BDAT <- BDAT %>%
  inner_join(
    bdat_agri_sf %>% dplyr::select(bdatid, INSEE_COM, CODE_US, codification,X,Y),
    by = "bdatid"
  )

# Conversion numérique des colonnes mesures BDAT (pho/argi/limt/sabt/corgox/cecmet)
BDAT <- BDAT %>%
  mutate(
    pho    = suppressWarnings(as.numeric(pho)),
    argi   = suppressWarnings(as.numeric(argi)),
    limt   = suppressWarnings(as.numeric(limt)),
    sabt   = suppressWarnings(as.numeric(sabt)),
    corgox = suppressWarnings(as.numeric(corgox)),
    cecmet = suppressWarnings(as.numeric(cecmet))
  )


# Table BDAT générique pour chaque propriété cible

BDAT_var <- BDAT %>%
  transmute(
    bdatid, annee, insee, INSEE_COM, codification, X, Y,
    !!name := .data[[ map_bdat[[name]] ]],
    source = "BDAT"
  ) %>%
  filter(!is.na(.data[[name]]))%>%
  rename(
    x = X,
    y = Y
  ) %>%
  mutate(
  annee     = as.factor(annee),
  INSEE_COM = as.character(INSEE_COM),
  x         = as.numeric(x),
  y         = as.numeric(y),
  !!name    := as.numeric(.data[[name]])
) %>%
  filter(!is.na(.data[[name]]))

saveRDS(BDAT_var, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/BDAT_53", name, ".rds"))

# 4.Pretraitement sur IGCS ----

##4.1. Identification des points en zone agricole ----
igcs$dt_cmp_ <- as.Date(igcs$dt_cmp_, format="%d/%m/%Y")

igcs_sf <- igcs %>%
  dplyr::select(id_prfl, dt_cmp_) %>%
  group_by(id_prfl, dt_cmp_) %>%
  summarise() %>%
  mutate(annee=as.factor(format(dt_cmp_, "%Y"))) %>%
  st_transform(crs=2154)

igcs_ocsol_sf <- st_join(igcs_sf, ocsol, join = st_intersects)

igcs_agri_sf <- igcs_ocsol_sf %>%
  filter(CODE_US=="US1.1") %>%
  dplyr::select(annee, id_prfl, CODE_US, geometry, dt_cmp_)

##4.2. Remobilisation du tableau avec les propriétés des sols ----
igcs_agri <- igcs %>%
  filter(id_prfl %in% igcs_agri_sf$id_prfl) %>%
  mutate(annee=as.factor(format(dt_cmp_, "%Y")))

##4.3. Harmonisation IGCS à 0-30 cm ----
# Renommer colonnes 
igcs_agri <- igcs_agri %>%
  mutate(limon=1000-(argile+sand)) %>%
  rename(
    pH   = ph_eau,
    CEC  = cec,
    C    = carbone,
    clay = argile,
    silt = limon
  ) %>%
  st_transform(crs=2154)

# Récupération INSEE de IGCS
igcs_agri <- igcs_agri %>%
  st_join(com %>% dplyr::select(INSEE_COM), join = st_intersects)

cds <- st_coordinates(igcs_agri)
igcs_agri <- igcs_agri %>% mutate(x=cds[,1], y=cds[,2])

###4.3.1. Transformation en dataframe ----
igcs_agri_df <- igcs_agri %>% 
  rename(id_profil = id_prfl,
         top       = prof_sp,
         bottom    = prof_nf,
  ) %>%
  st_drop_geometry() %>% 
  mutate(
    top = case_when(
      !is.na(top)               ~ top,
      is.na(top) & no_hrzn == 1 ~ 0,
      TRUE                      ~ lag(bottom)
    )
  )

###4.3.2. Filtrage des horizons sans limites définies + exclusion RMQS ----
igcs_agri_df <- igcs_agri_df %>%
  filter(!is.na(top), !is.na(bottom)) %>%
  filter(is.na(typ_pr_))   # exclut composites RMQS

# Sauvegarde RMQS si besoin
RMQS <- igcs_agri_df %>% filter( typ_pr_=="C" | typ_pr_=="F" )
saveRDS(RMQS, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/RMQS.rds")

###4.3.4. Harmonisation des horizons (générique pour 'name') ----
# Profils à horizon unique
pf_uniq <- igcs_agri_df %>%
  filter(top < bottom) %>% arrange(id_profil, top) %>%
  group_by(id_profil) %>% filter(n() == 1) %>% ungroup()

# Profils à 2 horizons avec 1er horizon 30-35 cm
P2hrzn1_30_35 <- igcs_agri_df %>%
  filter(top < bottom) %>% arrange(id_profil, top) %>%
  group_by(id_profil) %>% filter(n() == 2) %>% ungroup() %>%
  filter(no_hrzn==1, bottom >= 30 & bottom <= 35)

# Profils à 2 horizons à moyenne pondérée (hors cas précédent)
pf_mp <- igcs_agri_df %>%
  filter(top < bottom) %>% arrange(id_profil, top) %>%
  group_by(id_profil) %>% filter(n() == 2) %>% ungroup() %>%
  anti_join(P2hrzn1_30_35 %>% distinct(id_profil), by = "id_profil")

p_mp <- pf_mp %>%
  mutate(thick = bottom - top) %>%
  group_by(id_profil, annee, INSEE_COM, x, y) %>%
  summarise(
    !!name := round(weighted.mean(.data[[name]], thick, na.rm = TRUE), 1),
    .groups = "drop"
  )

# Profils à plus de 3 horizons avec 1er horizon 30-35 cm (prise directe)
p3hrzn1_30_35 <- igcs_agri_df %>%
  filter(top < bottom) %>% arrange(id_profil, top) %>%
  group_by(id_profil) %>% filter(n() >= 3) %>% ungroup() %>%
  filter(no_hrzn==1, bottom >= 30 & bottom <= 35)

# Profils à  plus de 3 horizons à spliner (hors cas précédent)
spl_dfs <- igcs_agri_df %>%
  filter(top < bottom) %>% arrange(id_profil, top) %>%
  group_by(id_profil) %>% filter(n() >= 3) %>% ungroup() %>%
  anti_join(p3hrzn1_30_35 %>% distinct(id_profil), by = "id_profil") %>%
  dplyr::select(id_profil, top, bottom, annee, INSEE_COM, x, y, !!name)

# Spline 0-30 cm (variable cible = name)
spl_layers <- mpspline(
  obj = spl_dfs,
  var_name = name,
  d = c(0, 30),
  lam = 0.1
)

col_spline <- paste0(name, "_000_030_cm")

igcs_spline <- purrr::imap_dfr(spl_layers, ~ {
  tibble(
    id_profil = .x$id_profil,
    !!col_spline := round(.x$est_dcm[["000_030_cm"]], 1),
    RMSE = .x$est_err[["RMSE"]],
    RMSE_IQR = .x$est_err[["RMSE_IQR"]]
  )
}) %>%
  left_join(
    spl_dfs %>% distinct(id_profil, INSEE_COM, annee, x, y),
    by = "id_profil"
  ) %>%
  rename(!!name := all_of(col_spline)) %>%
  dplyr::select(id_profil, annee, INSEE_COM, x, y, !!name)

# Harmonisation colonnes pour l’assemblage
pf_uniq       <- pf_uniq       %>% 
  dplyr::select(id_profil, annee, INSEE_COM, x, y, !!name)
p3hrzn1_30_35 <- p3hrzn1_30_35 %>% 
  dplyr::select(id_profil, annee, INSEE_COM, x, y, !!name)
P2hrzn1_30_35 <- P2hrzn1_30_35 %>% 
  dplyr::select(id_profil, annee, INSEE_COM, x, y, !!name)

igcs_final_var <- bind_rows(p_mp, pf_uniq, p3hrzn1_30_35, igcs_spline, P2hrzn1_30_35) %>%
  mutate(source="IGCS")
  

#Création de la base finale----
base_final <- bind_rows(igcs_final_var, BDAT_var)%>%
  filter(!is.na(.data[[name]]))

saveRDS(base_final, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_", name, ".rds"))

#================================================FIN DU SCRIPT====================================================#