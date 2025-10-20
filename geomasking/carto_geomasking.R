#===================================================================================================================#
# Script :      carto_geomasking 

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour cartographier les rasters issues du geomasking 

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    20-10-2025

# Entrees :     rasters pour chaque distance 

# Sorties :     Carte en facettes des rasters

# Modification : 20-10-2025
#===================================================================================================================#

#==========================================DEBUT DU SCRIPT=========================================================#

library(terra)
library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
library(purrr)
library(tibble)

#1. Définition des Paramètres ---
name    <- "arg"
out_dir <- "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/"
couleurs_argile <- c("#0000FF","#00BFFF","#00FFBF","#FFFFBF","#FF8000","#A52A2A","#800000")
valeurs_coupures <- c(0, 50, 100, 150, 200, 400, 600)
limits_cont <- c(100, 600)
breaks_leg  <- seq(100, 600, 100)

#2. Chargement des rasters ---
pat <- paste0("^", name, "qrf_geomask(\\d+)\\.tif$") #pattern pour extraire les rasters
files <- list.files(out_dir, pattern = pat, full.names = TRUE) #liste des fichiers correspondant au pattern
stopifnot(length(files) > 0)

file_tbl <- tibble(
  file = files,
  dist = as.integer(str_match(basename(files), pat)[, 2])
) %>%
  arrange(dist) %>%
  mutate(
    panel = factor(paste0("d = ", dist, " m"), levels = paste0("d = ", dist, " m"))
  ) # Ajout de la distance pour chaque raster

#3. Table de vérification (stats sur les rasters) ---
check_tbl <- purrr::map_dfr(
  file_tbl$file,
  ~{
    r <- terra::rast(.x)
    tibble::tibble(
      file_full = .x,
      file      = basename(.x),
      mean = mean(terra::values(r), na.rm = TRUE),
      sd   = sd(terra::values(r),   na.rm = TRUE),
      min  = suppressWarnings(min(terra::values(r), na.rm = TRUE)),
      max  = suppressWarnings(max(terra::values(r), na.rm = TRUE))
    )
  }
) # création d'un tableau récapitulatif avec des statistiques pour chaque raster

# table de correspondance "basename(file)" -> dist
file_key <- file_tbl %>%
  dplyr::mutate(file = basename(file)) %>%
  dplyr::select(file, dist)

# jointure pour ajouter la distance dans le tableau de vérification
check_tbl <- check_tbl %>%
  dplyr::left_join(file_key, by = "file") %>%
  dplyr::arrange(dist)

print(check_tbl)

# 4. Définir un tableau long pour la carte en facette ---
df_all <- purrr::map2_dfr(
  file_tbl$file, file_tbl$panel,  # panel = factor with levels in sorted order
  ~{
    r  <- terra::rast(.x)
    df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
    val_col <- names(df)[3]
    dplyr::transmute(df, x, y, val = .data[[val_col]], panel = .y)
  }
)

# 5. Carte en facettes ---
p <- ggplot(df_all, aes(x = x, y = y, fill = val)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~ panel, ncol = 5, drop = FALSE) +  
  scale_fill_gradientn(
    colours = couleurs_argile,
    limits  = limits_cont,
    values  = scales::rescale(valeurs_coupures),
    breaks  = breaks_leg,
    labels  = scales::label_number(accuracy = 1),
    oob     = scales::squish,
    name    = "Argile (g/kg)"
  ) +
  labs(title = paste0(name, " — RF par distance de géomasking")) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    strip.text  = element_text(face = "bold")
  )

print(p)
