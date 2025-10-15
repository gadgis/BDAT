#=============================================================================================================================#
# Script :      Script pour la visualisation des cartes

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Ce script produit des cartes (individuelles et en facettes) à partir de rasters
#               de prédiction pour une propriété de sol. Il :
#               + charge les rasters issus des trois modèles (RF, KO-INLA, KED-INLA) et des deux
#                   approches CSMS (« Désagrégation » vs « Données ponctuelles ») ;
#               + convertit chaque raster en data.frame (x, y, valeur) pour ggplot2 ;
#               + propose deux fonctions génériques :
#                   - raster_to_df() : conversion SpatRaster -> data.frame avec étiquettes de modèle/approche ;
#                   - plot_raster_facets() / plot_raster_indiv() : génération de cartes ggplot
#                     avec barre d’échelle et flèche du nord ;
#               + gère automatiquement l’échelle de couleur selon la nature de la variable :
#                   - continue (ex. argile) : dégradé (gradient) avec bornes et coupures optionnelles ;
#                   - discrète (ex. classes de pH) : palette manuelle (scale_fill_manual).
 

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    21-07-2025

# Entrees :     raster des prédictions sur la zone agricole

# Sorties :     carte de  la propriété cible selon l'approche de CSMS

# Modification : 09-10-2025
#===========================================================================================================================#
#================================================DEBUT DU SCRIPT============================================================#

# chargement des packages----

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(purrr)
library(ggspatial)
library(scales)     
library(grid)       
library(rlang)

# 1. Dossier de travail paramètres et fonctions----

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")
name <- "pH"

#  utilitaire de calcul
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Fonction pour transformer un raster en data.frame utilisable par ggplot2

raster_to_df <- function(r, modele, type_approche, prop_name = "val") {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  val_col <- colnames(df)[3]
  df <- df %>%
    dplyr::rename(!!prop_name := all_of(val_col)) %>%
    dplyr::mutate(modèle = modele, type = type_approche)
  return(df)
}


# Fonction pour la visualisation de toutes les cartes avec facettes selon le modèle et l'approche de CSMS


plot_raster_facets <- function(
    df,
    fill_col,                     # la colonne à cartographier 
    type_col   = "type",            # colonne des lignes de facet 
    model_col  = "modèle",          # colonne des colonnes de facet
    # paramètres pour continu
    palette_cont   = NULL,        # couleurs_argile
    cuts           = NULL,        # valeurs_coupures (pour values = rescale(cuts))
    limits_cont    = NULL,        # ex. c(100,600)
    # paramètres pour discret
    palette_disc   = NULL,        # ex. couleurs_ph 
    # communs
    fill_label     = NULL,        # titre de la légende
    na_color       = "transparent",
    base_size      = 12
) {
  # tidy-eval
  fill_col  <- enquo(fill_col)
  type_col  <- enquo(type_col)
  model_col <- enquo(model_col)
  
  # détection du type
  fill_name <- as_name(fill_col)
  is_cont   <- is.numeric(df[[fill_name]])
  
  p <- ggplot(df, aes(x = x, y = y, fill = !!fill_col)) +
    geom_raster() +
    facet_grid(rows = vars(!!type_col), cols = vars(!!model_col)) +
    annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.6) +
    annotation_north_arrow(
      location = "tl",
      which_north = "true",
      style  = north_arrow_fancy_orienteering(),
      height = unit(1, "cm"),
      width  = unit(1, "cm")
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(size = 12, face = "bold"),
      legend.position  = "right",
      plot.title       = element_text(hjust = 0.5),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
  
  # échelle selon le type
  if (is_cont) {
    p <- p + scale_fill_gradientn(
      colours = palette_cont %||% viridisLite::viridis(9),
      values  = if (!is.null(cuts)) rescale(cuts) else NULL,
      limits  = limits_cont,
      name    = if (is.null(fill_label)) fill_name else fill_label,
      na.value = na_color
    )
  } else {
    # si pas de palette fournie, palette automatique (hue) nommée par niveaux
    if (is.null(palette_disc)) {
      lv <- sort(unique(df[[fill_name]]))
      palette_disc <- setNames(scales::hue_pal()(length(lv)), lv)
    }
    p <- p + scale_fill_manual(
      values       = palette_disc,
      drop         = FALSE,
      na.value     = na_color,
      na.translate = FALSE,
      name         = if (is.null(fill_label)) fill_name else fill_label
    )
  }
  
  p
}


# Fonction pour produire des carte individuelle (sans facette) pour chaque raster

plot_raster_indiv <- function(
    df,
    fill_col,                 
    palette_cont = NULL,      
    cuts = NULL,              
    limits_cont = NULL,     
    palette_disc = NULL,      
    # communs
    fill_label = NULL,
    which_north = "grid",
    base_size = 12
) {
  fill_col <- enquo(fill_col)
  fill_name <- as_name(fill_col)
  is_cont <- is.numeric(df[[fill_name]])
  
  p <- ggplot(df, aes(x = x, y = y, fill = !!fill_col)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.6) +
    annotation_north_arrow(
      location = "tl",
      which_north = which_north,
      style  = north_arrow_fancy_orienteering()
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  if (is_cont) {
    # échelle continue (argile)
    p <- p + scale_fill_gradientn(
      colours = palette_cont %||% viridisLite::viridis(9),
      values  = if (!is.null(cuts)) rescale(cuts) else NULL,
      limits  = limits_cont,
      name    = fill_label %||% fill_name,
      na.value = "transparent"
    )
  } else {
    # échelle discrète (pH classé)
    if (is.null(palette_disc)) {
      lv <- sort(unique(df[[fill_name]]))
      palette_disc <- setNames(hue_pal()(length(lv)), lv)
    }
    p <- p + scale_fill_manual(
      values       = palette_disc,
      drop         = FALSE,
      na.value     = "transparent",
      na.translate = FALSE,
      name         = fill_label %||% fill_name
    )
  }
  
  p
}

# 2.Charger les rasters----

ked_cent <- rast(paste0(name,"predKEDINLA_final_cent.tif")) # A modifier selon la propriété cible
ko_cent  <- rast(paste0(name,"predKOINLA_final_cent.tif"))
rf_cent  <- rast(paste0(name,"qrf_final_cent.tif"))

ked_ponc <- rast(paste0(name,"predKEDINLA_final.tif"))
ko_ponc  <- rast(paste0(name,"predKOINLA_final.tif"))
rf_ponc  <- rast(paste0(name,"qrf_final.tif"))

# Appliquer la fonction pour créer les data.frames
df_ked_cent <- raster_to_df(ked_cent, "KED", "Désagrégation",     prop_name = "val")
df_ko_cent  <- raster_to_df(ko_cent,  "KO",  "Désagrégation",     prop_name = "val")
df_rf_cent  <- raster_to_df(rf_cent,  "RF",  "Désagrégation",     prop_name = "val")

df_ked_ponc <- raster_to_df(ked_ponc, "KED", "Données ponctuelles", prop_name = "val")
df_ko_ponc  <- raster_to_df(ko_ponc,  "KO",  "Données ponctuelles", prop_name = "val")
df_rf_ponc  <- raster_to_df(rf_ponc,  "RF",  "Données ponctuelles", prop_name = "val")


#3. Cartographie de la propriété pour chaque modèle selon les approches----


# Fusion de tous les jeux de données

df_all <- bind_rows(
  df_ked_cent, df_ko_cent, df_rf_cent,
  df_ked_ponc, df_ko_ponc, df_rf_ponc
)                                     #df_all sera utilisé pour obtenir une seule carte avec des facette selon 
                                      #le modèle et l'approche de CSMS

# Ordonner les facteurs

df_all$modèle <- factor(df_all$modèle, levels = c("KED", "KO", "RF"))
df_all$type <- factor(df_all$type, levels = c("Désagrégation", "Données ponctuelles"))

## 3.1. Pour l'argile----

# Définir la palette 

couleurs_argile <- c(
  "#0000FF",  # bleu foncé (<100 g/kg)
  "#00BFFF",  # bleu clair (~150)
  "#00FFBF",  # turquoise (~200)
  "#FFFFBF",  # jaune pâle (~300)
  "#FF8000",  # orange (~400)
  "#A52A2A",  # brun (~500)
  "#800000"   # rouge foncé (>600)
)

# Valeurs de coupure correspondantes
valeurs_coupures <- c(0, 50, 100, 150, 200, 400, 600)

# cration du graphique
p_arg <- plot_raster_facets(
  df_all,
  fill_col    = val,
  palette_cont = couleurs_argile,
  cuts         = valeurs_coupures,
  limits_cont  = c(100, 600),
  fill_label   = "Argile (g/kg)"
)


## 3.2. Pour le pH----

# Catégorisation par classe de pH
df_all <- df_all %>%
  mutate(classe_ph = cut(val,
                         breaks = c(3, 4, 5, 6, 6.5, 7, 8, 9, 10),
                         labels = c("[3 – 4[", "[4 – 5[", "[5 – 6[", "[6 – 6.5[", "[6.5 – 7[", "[7 – 8[", "[8 – 9[", "[9 – 10]"),
                         include.lowest = TRUE))

# Palette couleur
couleurs_ph <- c(
  "[3 – 4[" = "#ca0020",
  "[4 – 5[" = "#f46d43",
  "[5 – 6[" = "#fdae61",
  "[6 – 6.5[" = "#fee08b",
  "[6.5 – 7[" = "#d9ef8b",
  "[7 – 8[" = "#a6d96a",
  "[8 – 9[" = "#66bd63",
  "[9 – 10]" = "#1a9850"
)

# Création du graphique
p_ph <- plot_raster_facets(
  df_all,
  fill_col   = classe_ph,
  palette_disc = couleurs_ph,
  fill_label   = "pH (classes)"
)


# 4. Cartographie pour chaque raster séparément----

## 4.1. Argile----

p_arg_indiv <- plot_raster_indiv(
  df      = df_ked_ponc,# A adapter ici c'est le KED en approche ponctuelle 
  fill_col = val,
  palette_cont = couleurs_argile,
  cuts        = valeurs_coupures,
  limits_cont = c(100, 600),
  fill_label  = "Argile (g/kg)"
)
print(p_arg_indiv)


## 4.2. pH----
df_ked_ponc_ph  <- df_ked_ponc %>%
  mutate(classe_ph = cut(val,
                         breaks = c(3, 4, 5, 6, 6.5, 7, 8, 9, 10),
                         labels = c("[3 – 4[", "[4 – 5[", "[5 – 6[", "[6 – 6.5[", "[6.5 – 7[", "[7 – 8[", "[8 – 9[", "[9 – 10]"),
                         include.lowest = TRUE))
p_ph_indiv <- plot_raster_indiv(
  df       = df_ked_ponc_ph,# A adapter ici c'est le KED en approche ponctuelle
  fill_col = classe_ph,
  palette_disc = couleurs_ph,
  fill_label   =  "pH (classes)"
)
print(p_ph)

#============================================FIN DU SCRIPT=================================================================#

