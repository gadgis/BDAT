#Chargement des packages----
library(tidyverse)
library(sf)
library(tmap)
library(tmaptools)

#Chargement des données----

dpt <- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\dept53_rep.shp")
dC<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_C.rds")
dpH<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")

#Transformation en sf----
dC_sf<-st_as_sf(dC, coords = c("x", "y"), crs = 2154)
dpH_sf<-st_as_sf(dpH, coords = c("x", "y"), crs = 2154)

#carte du pH----
pH_map <- dpH_sf %>% 
  mutate(
    ph_class = cut(
      pH,
      breaks = c(-Inf, 6, 6.3, 6.5, 6.8, 7, Inf),
      labels = c("<6", "[6–6.3]", "]6.3–6.5]", "]6.5–6.8]", "]6.8–7]", ">7"),
      right  = TRUE                
    )
  )

pal_ph <- c(
  "<6"       = "#8c2d04",
  "[6–6.3]"  = "#cc4c02",
  "]6.3–6.5]"= "#ec7014",
  "]6.5–6.8]"= "#fdae6b",
  "]6.8–7]"  = "#a1d99b",
  ">7"       = "#31a354"
)
pal_ph <- pal_ph[levels(pH_map$ph_class)]

tmap_mode("plot") 

tm_shape(dpt) +
  tm_polygons(
    col = "white",
    border.col = "black",
    lwd = 0.2
  ) +
tm_shape(pH_map) +
  tm_dots(
    col     = "ph_class",
    size = 0.45,
    palette = pal_ph,
    title   = "Classes de pH",
    border.col = "white",
    lwd       = 0.2
  ) +
  tm_layout(
    title = "Distribution du pH ",
    title.size = 1.2,
    legend.outside = TRUE,
    legend.title.size = 0.9,
    legend.text.size  = 0.8,
    frame = FALSE
  )

