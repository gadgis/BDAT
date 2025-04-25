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
d_pH<-read.csv("BDAT_IGCS_pH.csv", sep=",", header=TRUE)

d_pH <- d_pH %>%
  select(id_profil,annee,x,y,pH,source,bdatid,insee, INSEE_COM,codification)

#Transformation en sf----
d_pH_sf<-st_as_sf(d_pH, coords = c("x", "y"), crs = 2154)

#Agrégation des propriétés par commune

pH_com <- d_pH_sf %>%
  group_by(INSEE_COM) %>%
  summarise(
    moy_ph = round(mean(pH, na.rm = TRUE), 1),
    med_ph = round(median(pH, na.rm = TRUE), 1),
    sd_ph  = round(sd(pH, na.rm = TRUE), 1),
    n      = n()
  )


#Jointure de l'agrégation" aux communes

mayenne_ph <- st_join(communes, pH_com, join = st_intersects)


#Représentation graphique

mayenne_ph <- mayenne_ph %>%
  mutate(ph_class = cut(
    moy_ph,
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
  geom_sf(data = mayenne_ph, aes(fill = ph_class), color = "white") +
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
