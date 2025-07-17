run_rf <- function(approach = c("Ponctuelle", "Centroide"),
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
  pred_col <- paste0("predRF_", 
                     ifelse(approach == "Ponctuelle", "P", "C"),
                     ifelse(type_val == "Classique", "C", "S"))
  
  # Apprentissage
  rf_task <- makeRegrTask(data = data_train[, c(name, cov_brt)], target = name)
  res_tune <- tuneRanger(rf_task, 
                         num.trees = ntree, iters = 100,
                         num.threads = kmax)
  
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
  
  # Jeu de test
  test <- data_test[, c("id", name, "INSEE_COM", NomsCoord, cov_brt)]
  
  # Prédiction
  preds <- predict(rf_model, data = test[, cov_brt], num.threads = kmax)$predictions
  preds <- round(preds, 2)
  test[[pred_col]] <- preds
  
  # Évaluation
  eval <- Myeval(preds, test[[name]]) %>%
    mutate(method = paste0("RF_", substr(approach, 1, 1), substr(type_val, 1, 1)))
  
  return(list(
    preds = preds,
    evaluation = eval,
    detail = test[, c("id", name, "INSEE_COM", NomsCoord, pred_col)]
  ))
}
