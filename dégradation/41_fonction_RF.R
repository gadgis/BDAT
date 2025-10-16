#=============================================================================================================================#
# Script :      Script utilitaire pour le modèle de Random Forest

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Ce script fournit la fonction utilitaire `run_rf()` qui entraîne un modèle
#               Random Forest sur un jeu d’apprentissage puis prédit sur un jeu test, avec
#               réglage automatique d’hyperparamètres.
#
#               D'abord on prépare une tâche de régression (mlr) à partir de `data_train` où :
#                   - cible = `name`
#                   - covariables = colonnes listées dans `cov_brt`

#               Ensuite on effectue un tuning via `tuneRanger` (10 itérations) pour estimer
#              `mtry` et `min.node.size`, 
        
#               Enfin on entraîne un modèle `ranger`(ntree = 350, max.depth = 15, importance = "permutation").

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    21-07-2025

# Entrees :     Jeu de données sur les propriétés des sols, les variables importantes retenues (cov_brt)

# Sorties :     data.frame du jeu test avec colonnes id, la cible (name), INSEE_COM, NomsCoord et
#               la prédiction RF pred

# Modification : 08-10-2025
#===========================================================================================================================#

#================================================DEBUT DU SCRIPT============================================================#


run_rf <- function(approach = c("Ponctuelle", "Désagrégation"),
                   type_val = c("Classique", "Spatiale"),
                   data_train,
                   data_test,
                   cov_brt,
                   moyenne_covariable,
                   name,
                   ntree,
                   kmax,
                   NomsCoord) {
  
  approach <- match.arg(approach)
  type_val <- match.arg(type_val)
  
  # Nom de la colonne de prédiction
  pred_col <- "pred"
  
  # Apprentissage et calibration du modèle
  
  rf_task <- makeRegrTask(data = data_train[, c(name, cov_brt)], target = name)
  res_tune <- tuneRanger(rf_task, 
                         num.trees = ntree,
                         iters = 10,
                         num.threads = kmax,
                         show.info =  FALSE)
  
  rf_model <- ranger(
    formula = as.formula(paste0(name, " ~ .")),
    data = data_train[, c(name, cov_brt)],
    num.trees = ntree,
    mtry = res_tune$recommended.pars$mtry,
    min.node.size = res_tune$recommended.pars$min.node.size,
    quantreg = FALSE,
    max.depth = 15,
    importance = "permutation",
    scale.permutation.importance = TRUE,
    keep.inbag = FALSE
  )
  
  # Identification du jeu test
  test <- data_test[, c("id", name, "INSEE_COM", NomsCoord, cov_brt)]
  
  # Prédiction
  
  predCal <- rf_model$predictions #prédiction qui vont servir dans l'approche KED
  preds <- predict(rf_model, data = test[, cov_brt], num.threads = kmax)$predictions
  preds <- round(preds, 2)
  test[[pred_col]] <- preds
  
  # Évaluation
  eval <- Myeval(preds, test[[name]]) %>%
    mutate(method = paste0("RF_", substr(approach, 1, 1), substr(type_val, 1, 1)))
  
  return(list(
    predCal = predCal, #dérive externe pour le KED
    evaluation = eval,
    detail = test[, c("id", name, "INSEE_COM", NomsCoord, pred_col)]
    )
  )
}
