#===================================================================================================================#
# Script :      Random Forest 

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la prédiction des propriétés des sols en utilisant 
#               la méthode Random Forest et INLA SPDE pour chaque distance de geomasking

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Observations ponctuelles (geomasquées) et les covariables 

# Sorties :     Praster de prédiction sur la zone agricole

# Modification : 06-10-2025
#===================================================================================================================#

#==========================================DEBUT DU SCRIPT=========================================================#

# ============================================================


# RANDOM FOREST : fit sur mask_df, prédiction sur gXY 

# prédiction sur la grille 
testD <- gXY %>%
  dplyr::select(dplyr::all_of(cov_brt[-1]))

# Prédictions RF
QRF_Median2 <- predict(
  rf_model,
  testD,
  num.threads = kmax
)$predictions


QRF_Median50 <- dplyr::bind_cols(
  gXY %>%
    dplyr::filter_at(dplyr::vars(cov_brt[-1]), dplyr::all_vars(!is.na(.))) %>%
    dplyr::select(x, y),
  QRF_Median = QRF_Median2
)

# Rasterisation RF puis masque ZA 

r_qrf_full <- terra::rast(QRF_Median50, type = "xyz")
terra::crs(r_qrf_full) <- crs_epsg
names(r_qrf_full) <- "qrf"

r_qrf_out <- terra::extend(r_qrf_full, rast_za, snap = "near")
qrf_agri  <- terra::mask(r_qrf_out, rast_za)

terra::writeRaster(
  qrf_agri,
  filename = paste0(out_dir, name, "qrf_geomask", d, ".tif"),
  overwrite = TRUE
)

# INLA options 
INLA::inla.setOption(
  num.threads = max(1, floor(kmax / 2)),
  inla.mode   = "experimental"
)

#  KO : fit INLA sur mask_df, prédiction sur la grille 

# Préparer dataINLA (colonne 'activ' fidèle à ton style)
dataINLA <- mask_df[, c("x_moved","y_moved", name)]
names(dataINLA)[names(dataINLA) == name] <- "activ"
coords <- as.matrix(dataINLA[, c("x_moved","y_moved")])

sp::coordinates(dataINLA) <- c("x_moved","y_moved")
sp::proj4string(dataINLA) <- sp::CRS(crs_epsg)

# Mesh 
max.edge    <- diff(range(coords[,1])) / (3 * 5)
bound.outer <- diff(range(coords[,1])) / 3
mesh3 <- INLA::inla.mesh.2d(
  loc      = coords,
  max.edge = c(1, 2) * max.edge,
  offset   = c(max.edge, bound.outer),
  cutoff   = max.edge / 10
)
matern <- INLA::inla.spde2.pcmatern(
  mesh3, alpha = 2, prior.sigma = c(10, 0.01), prior.range = c(500, 0.01)
)

# Grille pxl 
pxl <- as.data.frame(r_qrf_full, xy = TRUE)
colnames(pxl)[3] <- "qrf"      # KO s'en fiche, KED l'utilise
sp::gridded(pxl) <- ~ x + y

# Modèle KO
cmp_KO <- activ ~ Intercept(1) + field(coordinates, model = matern)

fitKO <- inlabru::bru(
  components = cmp_KO,
  data       = dataINLA,
  family     = "Gaussian",
  domain     = list(coordinates = mesh3),
  options    = list(control.inla = list(int.strategy = "eb"), verbose = FALSE)
)

# Prédiction KO -> raster non masqué puis masque ZA à la sortie
predKO <- predict(
  fitKO,
  n.samples = nsim,
  pxl,
  ~ Intercept + field,
  num.threads = max(1, floor(kmax / 2))
)

r_ko_out <- terra::extend(terra::rast(predKO)[["mean"]], rast_za, snap = "near")
r_ko_out <- terra::mask(r_ko_out, rast_za)
terra::writeRaster(
  r_ko_out,
  filename = paste0(out_dir, name, "predKOINLA_geomask", d, ".tif"),
  overwrite = TRUE
)

#  KED : fit INLA avec dérive RF, prédiction sur la grille 

# Dérive RF aux points d’apprentissage (depuis le raster RF plein)
dataINLA$qrf <- terra::extract(r_qrf_full, terra::vect(dataINLA))[,1]

cmp_KED <- activ ~ Intercept(1) + rfpred(qrf, model = "linear") + field(coordinates, model = matern)

fitKED <- inlabru::bru(
  components = cmp_KED,
  data       = dataINLA,
  family     = "Gaussian",
  domain     = list(coordinates = mesh3),
  options    = list(control.inla = list(int.strategy = "eb"), verbose = FALSE)
)

# Prédiction KED -> raster non masqué puis masque ZA à la sortie
predKED <- predict(
  fitKED,
  n.samples = nsim,
  pxl,
  ~ Intercept + field + rfpred,
  num.threads = max(1, floor(kmax / 2))
)

r_ked_out <- terra::extend(terra::rast(predKED)[["mean"]], rast_za, snap = "near")
r_ked_out <- terra::mask(r_ked_out, rast_za)
terra::writeRaster(
  r_ked_out,
  filename = paste0(out_dir, name, "predKEDINLA_geomask", d, ".tif"),
  overwrite = TRUE
)
