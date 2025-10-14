library(sf)
library(dplyr)
library(terra)

# Sous-fonction : déplacement aléatoire (angle & pv par point)
deplacer_aleatoire <- function(x, y, distance) {
  n <- length(x)
  angle <- runif(n, 0, 2*pi)        # un angle par point
  pv    <- runif(n, 0.8, 1.2) # facteur variation de la distance
  list(
    x = x + pv * distance * cos(angle),
    y = y + pv * distance * sin(angle)
  )
}

# Fonction principale 

geomasking <- function(dataset,  NomsCoord, dist_seq , seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  x0  <- dataset[[ NomsCoord[1] ]]
  y0  <- dataset[[ NomsCoord[2] ]]
  ids <- dataset[[ "ID" ]]
  
  res_list <- lapply(dist_seq, function(d) {
    # graine différente par distance pour des tirages indépendants (optionnel)
    if (!is.null(seed)) set.seed(seed + as.integer(d))
    moved <- deplacer_aleatoire(x0, y0, distance = d)
    data.frame(
      ID         = ids,
      x_origine  = x0,
      y_origine  = y0,
      d          = d,
      x_nouveau  = moved$x,
      y_nouveau  = moved$y,
      row.names  = NULL
    )
  })
  
  do.call(rbind, res_list)
}

#chargement des données
name <- "pH"
NomsCoord <- c("x", "y")
dist_seq <- seq(100, 2500, 100)
data<-readRDS(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_", name, ".rds"))
zone_agri<-rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_za.tif")

#Creation d'une colonne ID dans data
data<- data %>%
  mutate(ID = row_number())

#Déplacement aléatoire des points
res_geomask <- geomasking(
  dataset = data,
  NomsCoord,
  dist_seq 
)
res_geomask

#Verification si présence dans la zone agricole pour chaque distance

# Conversion du data.frame en sf
res_sf <- st_as_sf(
  res_geomask,
  coords = c("x_nouveau", "y_nouveau"),
  crs = 2154  
)
# Extraction des valeurs du raster à chaque point
ext <- terra::extract(zone_agri, terra::vect(res_sf), bind = TRUE)

# nom de la colonne raster (souvent "layer" ou le nom du fichier)
val_col <- names(zone_agri)[1]

# 3) Marque "en zone agricole"
ext$in_agri <- !is.na(ext[[val_col]]) & ext[[val_col]] == 1  # marche si hors-zone = NA ou 0

# 4) Résumé par distance d
effectifs_agri <- ext |>
  as.data.frame() |>
  group_by(d) |>
  summarise(
    n_total = n(),
    n_agri  = sum(in_agri),
    prop_agri = round(100 * n_agri / n_total, 2),
    .groups = "drop"
  ) |>
  arrange(d)

effectifs_agri

# garder les points en zone agricole,
res_in_agri <- res_sf[ext$in_agri, ]
