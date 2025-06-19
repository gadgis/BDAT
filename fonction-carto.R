library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(ggspatial)
library(raster)
library(sf)
library(patchwork)

plot_raster_stack <- function(rasters, raster_names, breaks, labels, palette = "Spectral") {
  # Vérif
  stopifnot(length(rasters) == length(raster_names))
  
  # Stack raster
  r_stack <- rast(rasters)
  names(r_stack) <- raster_names
  
  # Convertir en data.frame long
  r_df <- as.data.frame(r_stack, xy = TRUE, na.rm = TRUE) %>%
    pivot_longer(cols = -c(x, y), names_to = "Méthode", values_to = "pH") %>%
    mutate(pH = round(pH, 2)) %>%
    mutate(pH_class = factor(cut(pH, breaks = breaks, labels = labels, include.lowest = TRUE),
                             levels = labels))
  # Forcer l'affichage de toutes les classes même absentes
  fictives <- data.frame(
    x = NA, y = NA, Méthode = raster_names[1], pH = NA, pH_class = factor(labels, levels = labels)
  )
  
  r_df <- bind_rows(r_df, fictives)
  
  # Palette
  colors <- c("#D53E4F", "#F57647", "#FDBC6C" ,"#FEEDA1" ,"#F0F9A8" ,"#BBE3A0" ,"#6FC5A4" ,"#3288BD")
  
  
  # Carte
  ggplot(r_df, aes(x = x, y = y, fill = pH_class)) +
    geom_raster() +
    coord_sf(crs = 2154, datum = NA) +
    facet_wrap(~ Méthode, ncol=length(rasters)) +
    scale_fill_manual(
      values = colors,
      name = "pH",
      drop = FALSE
    ) +
    annotation_scale(location = "bl", width_hint = 0.15,
                     line_width = 0.8, text_cex = 0.8) +
    annotation_north_arrow(location = "tl", which_north = "true",
                           height = unit(1.2, "cm"),
                           width = unit(1.2, "cm"),
                           style = north_arrow_fancy_orienteering) +
    theme_minimal() +
    labs(title = "Prédiction du pH") +
    theme(
      strip.text = element_text(face = "bold", size=12),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
}



#Chargement des données rasters

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")


#ETUDE DE LA DISTRIBUTION----
# Fichiers raster
rt<- c("pHqrf_centroides_final.tif", "predKOINLA_centroides_final.tif", "predKEDINLA_centroides_final.tif",
            "pHqrf_final.tif", "predKOINLA_final.tif","predKEDINLA_final.tif")

noms_t<- c("RF_cent", "KO_INLA_cent", "KED_INLA_cent","RF", "KO_INLA", "KED_INLA")

raster_stack <- rast(rt)
names(raster_stack) <- noms_t

#  Extraire les valeurs sous forme de data.frame
r_vals <- as.data.frame(raster_stack, na.rm = TRUE)

#  Long format pour ggplot
r_vals_long <- r_vals %>%
  pivot_longer(cols = everything(), names_to = "Méthode", values_to = "pH")
summary(r_vals_long)

# Résumé statistique par méthode
resum_stats <- r_vals_long %>%
  group_by(Méthode) %>%
  summarise(
    min = round(min(pH), 2),
    q25 = round(quantile(pH, 0.25), 2),
    median = round(median(pH), 2),
    mean = round(mean(pH), 2),
    q75 = round(quantile(pH, 0.75), 2),
    max = round(max(pH), 2),
    sd = round(sd(pH), 2),
    .groups = "drop"
  )

print(resum_stats)

#CARTOGRAPHIE-----
# Classes de pH
breaks <- c(-Inf, 6, 6.2,6.5, 6.8, 7, 8, 9, Inf)

labels <- c("[3 - 6[", "[6 – 6.2[", "[6.2 – 6.5[", "[6.5 – 6.8[", "[6.8 – 7[",

          "[7 – 8[", "[8 – 9[", "[9 - 10]")

# Appel fonction 
r <- c("predKEDINLA_final.tif")
noms <- c("KED_INLA")
plot_raster_stack(r, noms, breaks, labels)


#CARTE BDAT----
pH_med<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_median.rds")
communes <- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/COMMUNE.shp")
pH_med <- pH_med %>%
  mutate(pH_class = factor(cut(med_pH, breaks = breaks, labels = labels, include.lowest = TRUE),
                           levels = labels))

#jointure de pH_med à communes
pH_med <- pH_med %>%
  mutate(INSEE_COM= as.character(INSEE_COM))
         
communes<-left_join(communes, pH_med, by = "INSEE_COM")

communes <- communes %>%
  filter(!is.na(pH_class))

colors <- c("#D53E4F", "#F57647", "#FDBC6C" ,"#FEEDA1" ,"#F0F9A8" ,"#BBE3A0" ,"#6FC5A4" ,"#3288BD")

# Créer des lignes fictives avec géom vide pour forcer toutes les classes dans la légende
col_template <- st_drop_geometry(communes)[1, ]
col_template[] <- NA
fictives_df <- col_template[rep(1, length(labels)), ]
fictives_df$pH_class <- factor(labels, levels = labels)

# Géométrie vide (même CRS que communes)
geom_na <- st_sfc(rep(st_geometry(st_point()), length(labels)), crs = st_crs(communes))
fictives_sf <- st_sf(fictives_df, geometry = geom_na)

#  Fusionner
communes_compl <- rbind(communes, fictives_sf)

#  Carte
g_commune<-ggplot(communes_compl) +
  geom_sf(aes(fill = pH_class), color = "white", size = 0.2) +
  scale_fill_manual(values = colors, name = "pH", drop = FALSE) +
  coord_sf(crs = 2154, datum = NA) +
  annotation_scale(location = "bl", width_hint = 0.15) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  labs(title = "pH moyen par commune") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

