#Carte type BDAT en faisant une agrégation par commune



#Chargement de package----
library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggmap)
library(ggspatial)
library(prettymapr)

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

#Importation des couches----
communes <- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\communes53_2014.shp")
dpt <- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\dept53_rep.shp")
d_C<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_C.rds")
d_pH<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")
d_C <- d_C %>%
  select(id_profil,annee,x,y,C,source,bdatid,insee, INSEE_COM,codification)
d_pH <- d_pH %>%
  select(id_profil,annee,x,y,pH,source,bdatid,insee, INSEE_COM,codification)
#Transformation en sf----
d_C_sf<-st_as_sf(d_C, coords = c("x", "y"), crs = 2154)
d_pH_sf<-st_as_sf(d_pH, coords = c("x", "y"), crs = 2154)
#Agrégation des propriétés par commune

C_com <- d_C_sf %>%
  group_by(INSEE_COM) %>%
  summarise(
    moy_C = round(mean(C, na.rm = TRUE), 1),
    med_C = round(median(C, na.rm = TRUE), 1),
    sd_C  = round(sd(C, na.rm = TRUE), 1),
    n      = n()
  )

C_com_df<-C_com %>%
  st_drop_geometry()

pH_com <- d_pH_sf %>%
  group_by(INSEE_COM) %>%
  summarise(
    moy_pH = round(mean(pH, na.rm = TRUE), 1),
    med_pH= round(median(pH, na.rm = TRUE), 1),
    sd_pH  = round(sd(pH, na.rm = TRUE), 1),
    n      = n()
  )

pH_com_df<-pH_com %>%
  st_drop_geometry()

#Jointure de l'agrégation" aux communes

mayenne_pH <- st_join(communes, pH_com, join = st_intersects)
mayenne_C <- st_join(communes, C_com, join = st_intersects)


#Représentation graphique
#Carte de la moyenne de pH par commune
mayenne_pH <- mayenne_pH %>%
  mutate(ph_class = cut(
    moy_pH,
    breaks = c(-Inf, 6, 6.3, 6.5, 6.8, 7, Inf),
    labels = c(
      "<6",
      "[6-6.3]",
      "]6.3-6.5]",
      "]6.5-6.8]",
      "]6.8-7]",
      ">7"
    )
  ))
ggplot() +
  geom_sf(data = mayenne_pH, aes(fill = ph_class), color = "white") +
  scale_fill_manual(
    values = c(
      "<6"      = "#8c2d04",
      "[6-6.3]" = "#cc4c02",
      "]6.3-6.5]"= "#ec7014",
      "]6.5-6.8]"  = "#fdae6b",
      "]6.8-7]" = "#a1d99b",
      ">7"          = "#31a354"
    ),
    name = "Classes pH",
    na.value = "grey80"
  ) +
  theme_minimal() +
  labs(
    title = "pH moyen par commune" )
#Carte de la moyenne de C par commune
mayenne_C <- mayenne_C %>%                     
  mutate(c_org_class = cut(
    moy_C,                                  
    breaks = c(-Inf, 15, 20, 25, 30, 35, Inf), 
    labels = c(
      "<15",            
      "[15-20]",        
      "]20-25]",        
      "]25-30]",        
      "]30-35]",        
      ">35"             
    )
  ))

ggplot() +
  geom_sf(data = mayenne_C, aes(fill = c_org_class), color = "white") +
  scale_fill_manual(
    values = c(
      "<15"      = "#8c2d04",
      "[15-20]"  = "#cc4c02",
      "]20-25]"  = "#ec7014",
      "]25-30]"  = "#fdae6b",
      "]30-35]"  = "#a1d99b",
      ">35"      = "#31a354"
    ),
    name = "Classes C org (g kg-¹)",
    na.value = "grey80"
  )+
  theme_minimal() +
  labs(
    title = "Carbone moyen par commune" )

#Export des fichiers----
saveRDS(C_com_df, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/C_moyen.rds")
saveRDS(pH_com_df, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_moyen.rds")
