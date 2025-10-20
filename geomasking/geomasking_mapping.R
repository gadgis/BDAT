#===================================================================================================================#
# Script :      Random Forest 

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la prédiction des propriétés des sols en utilisant 
#               la méthode Random Forest et INLA SPDE pour chaque distance de geomasking

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Observations ponctuelles (geomasquées) et les covariables 

# Sorties :     Raster de prédiction sur la zone agricole

# Modification : 06-10-2025
#===================================================================================================================#

#==========================================DEBUT DU SCRIPT=========================================================#



#1. Prédiction par le random Forest----

#Calibration du modèle sur les données masquées (mask_df) , prédiction sur la grille gXY 

covonly = c(name, cov_brt) #ajout du nom de la variable cible à la liste des covariables sélectionnées
datacov_shrt = mask_df[covonly]                   #récupération dans datacov des colonnes conservées

SOC.task = makeRegrTask(data = datacov_shrt, target = name)

# Estimation approximative du temps de réglage
estimateTimeTuneRanger(SOC.task,num.threads = 60,
                       num.trees = ntree)

#Tunning des hyperparamètres du modèle Random Forest
res = tuneRanger(SOC.task,
                 num.trees = ntree,
                 iters = 100,
                 num.threads = 60)


# Affichage des hyperparamètres optimaux

res$recommended.pars$mtry
res$recommended.pars$min.node.size

##1.1Calibration du modèle Random Forest avec les hyperparamètres optimaux----

fomula.ranger <- as.formula(paste0(name,"~."))
rf_model <- ranger(formula = fomula.ranger ,
                   data = datacov_shrt,
                   num.trees = ntree,
                   min.node.size = res$recommended.pars$min.node.size ,
                   quantreg = F,
                   max.depth = 15, 
                   mtry=res$recommended.pars$mtry ,
                   importance="permutation", 
                   scale.permutation.importance = TRUE, #division par l'écart-type de la variable (mise des permutations entre 0 et 1)
                   keep.inbag = F)

##1.2. prédiction sur la grille gXY----
testD <- gXY %>%
  dplyr::select(dplyr::all_of(covonly[-1])) #selectionne que les covariables à utiliser pour la prédiction

# Prédictions RF
QRF_Median2 <- predict(
  rf_model,
  testD,
  num.threads = kmax
)$predictions


QRF_Median50 <- dplyr::bind_cols(
  gXY %>%
    dplyr::filter_at(dplyr::vars(covonly[-1]), dplyr::all_vars(!is.na(.))) %>%
    dplyr::select(x, y),
  QRF_Median = QRF_Median2
)

# Rasterisation RF puis masque ZA 

r_qrf_full <- terra::rast(QRF_Median50, type = "xyz")
terra::crs(r_qrf_full) <- crs_epsg
names(r_qrf_full) <- "qrf"

r_qrf_out <- terra::extend(r_qrf_full, rast_za, snap = "near") #étendre le raster RF à la zone agricole
qrf_agri  <- terra::mask(r_qrf_out, rast_za) #appliquer le masque ZA

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

#2.Prédiction par KO----

##2.1. Préparation du jeu d'entrée dataINLA ----
dataINLA <- mask_df[, c("x_moved","y_moved", name)]
names(dataINLA)[names(dataINLA) == name] <- "activ"
coords <- as.matrix(dataINLA[, c("x_moved","y_moved")])

sp::coordinates(dataINLA) <- c("x_moved","y_moved")
sp::proj4string(dataINLA) <- sp::CRS(crs_epsg)

##2.2. définir la Mesh et la grille de prédiction----
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
colnames(pxl)[3] <- "qrf"      
sp::gridded(pxl) <- ~ x + y

##2.3. Définir le Modèle KO, ajustement puis prédiction----
cmp_KO <- activ ~ Intercept(1) + field(coordinates, model = matern)

fitKO <- inlabru::bru(
  components = cmp_KO,
  data       = dataINLA,
  family     = "Gaussian",
  domain     = list(coordinates = mesh3),
  options    = list(control.inla = list(int.strategy = "eb"), verbose = FALSE)
)

# Prédiction KO 
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

#3.Prédiction par KED----

##3.1. Calibration du KED----
dataINLA$qrf <- terra::extract(r_qrf_full, terra::vect(dataINLA))[,1]

cmp_KED <- activ ~ Intercept(1) + rfpred(qrf, model = "linear") + field(coordinates, model = matern)

fitKED <- inlabru::bru(
  components = cmp_KED,
  data       = dataINLA,
  family     = "Gaussian",
  domain     = list(coordinates = mesh3),
  options    = list(control.inla = list(int.strategy = "eb"), verbose = FALSE)
)

##3.2. Prédiction et masque sur la ZA
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
