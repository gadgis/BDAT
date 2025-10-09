#=============================================================================================================================#
# Script :      Script utilitaire pour le modèle de Random Forest

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Ce script fournit la fonction utilitaire `run_inla_spde_core()` qui entraîne un modèle
#               de Krigeage ordinaire et de Krigeage avec dérive externe sur un jeu d’apprentissage
#               puis prédit sur un jeu test.
#
#               D'abord les données sont convertit  en objets `sp. Ensuite un maillage (mesh) 2D est modélisé
#               avec `INLA::inla.mesh.2d` et un champ spatial Matérn est définit via `inla.spde2.pcmatern`.

#               Puis le modèle est définit puis ajusté :
#                   - KO  (Krigeage Ordinaire) : Intercept + champ spatial,
#                   - KED (Krigeage à Dérive Externe) : Intercept + dérive linéaire `pred` + champ spatial.

#               Enfin les prédictions sont effectuées sur le jeu test via `inlabru::predict()`.

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    21-07-2025

# Entrees :     data.frame d’apprentissage et de test avec des colonnes pour la variable cible `name`, 
#               les coordonnées `NomsCoord`, les predictions du Random Forest pour le KED. 

# Sorties :     data.frame du jeu test avec la prédiction (median) des modèles de KO et KED

# Modification : 09-10-2025
#===========================================================================================================================#

#================================================DEBUT DU SCRIPT============================================================#

run_inla_spde_core <- function(dataINLA, 
                               data_test,
                               name,
                               type = c("KO", "KED"),
                               approach = c("Ponctuelle", "Désagrégation"),
                               type_val = c("Classique", "Spatiale"),
                               NomsCoord = c("x", "y")
                               ) {
  
  type <- match.arg(type)
  approach <- match.arg(approach)
  type_val <- match.arg(type_val)
  
  
  # Transformation en sp pour INLA
  dataINLA$elt <- dataINLA[,name]
  
  coords <- dataINLA[, NomsCoord]
  
  coordinates(dataINLA) <- NomsCoord
  proj4string(dataINLA) <- CRS("epsg:2154")
  
  coordinates(data_test) <- NomsCoord
  proj4string(data_test) <- CRS("epsg:2154")
  
  # Création du mesh
  
  max.edge <- diff(range(coords[, 1])) / (3 * 5)
  bound.outer <- diff(range(coords[, 1])) / 3
  mesh3 <- inla.mesh.2d(
    loc = coords,
    max.edge = c(1, 2) * max.edge,
    offset = c(max.edge, bound.outer),
    cutoff = max.edge / 10
  )
  
  matern <- INLA::inla.spde2.pcmatern(
    mesh3, alpha = 2,
    prior.sigma = c(10, 0.01),
    prior.range = c(500, 0.01)
  )
  
  # Définition du modèle
  cmp <- if (type == "KO") {
    elt ~ Intercept(1) + field(coordinates, model = matern)
  } else {
    elt ~ Intercept(1) + rfpred(pred, model = "linear") + field(coordinates, model = matern)
  }
  
  cmpPred <- if (type == "KO") {
    ~ Intercept + field
  } else {
    ~Intercept + rfpred + field
  }
  
  # Ajustement du modèle
    fit <- bru(
    cmp,
    data = dataINLA,
    family = "Gaussian",
    
    options = list(control.inla = list(int.strategy = "eb"), 
                   verbose = FALSE)
  )
  
  gc()
  
  
  # Prédictions
   predObj <-  predict(fit,
                       data_test ,
                      cmpPred  
                       
                     )
  

  return(preds = predObj$median)
}