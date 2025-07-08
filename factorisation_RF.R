# Fonction factorisée pour entraîner un Random Forest sur données agrégées

run_rf <- function(method = c("Nicolas", "Lucille"),
                   data_sample,
                   fold_structure,
                   cov_brt,
                   cat_vars,
                   num_vars,
                   name,
                   k,
                   ntree,
                   kmax) {
  
  method <- match.arg(method)
  pred_col <- ifelse(method == "Nicolas", "predRF_N", "predRF_L")
  data_sample[[pred_col]] <- NA_real_
  
  resuXval <- foreach(i = 1:k, .combine = rbind, .packages = c("dplyr", "ranger", "mlr", "tuneRanger")) %do% {
    cat("Fold:", i, "[", method, "]\n")
    
    #1.Séparation train/test----
    
    if (method == "Nicolas") {
      nblignes <- fold_structure[[i]]
      train_data <- data_sample[-nblignes, ]
      test_data  <- data_sample[nblignes, ]
    } else {
      test_communes  <- fold_structure[[i]]
      train_communes <- setdiff(unique(data_sample$INSEE_COM), test_communes)
      train_data <- data_sample %>% filter(INSEE_COM %in% train_communes)
      test_data  <- data_sample %>% filter(INSEE_COM %in% test_communes)
    }
    
    #2. Agrégation des données d'entraînement----
    
    train_aggr <- train_data %>%
      group_by(INSEE_COM) %>%
      summarise(
        across(all_of(num_vars), ~mean(.x, na.rm = TRUE)),
        across(all_of(cat_vars), ~as.character(names(which.max(table(.))))),
        .groups = "drop"
      )
    
    #3. Agrégation de la variable cible----
    y_aggr <- train_data %>%
      group_by(INSEE_COM) %>%
      summarise(value = mean(.data[[name]], na.rm = TRUE), .groups = "drop") %>%
      pull(value)
    train_aggr[[name]] <- y_aggr
    
    # Conversion des types
    train_aggr[cat_vars] <- lapply(train_aggr[cat_vars], factor)
    train_aggr[num_vars] <- lapply(train_aggr[num_vars], as.numeric)
    
    # Tuning avec mlr + tuneRanger
    rf_task <- makeRegrTask(data = train_aggr[, c(name, cov_brt)], target = name)
    res_tune <- tuneRanger(rf_task, num.trees = ntree, iters = 100, num.threads = 4)
    
    # 4. Entraînement du modèle----
    
    rf_model <- ranger(
      formula = as.formula(paste0(name, " ~ .")),
      data = train_aggr[, c(name, cov_brt)],
      num.trees = ntree,
      mtry = res_tune$recommended.pars$mtry,
      min.node.size = res_tune$recommended.pars$min.node.size,
      importance = "permutation",
      keep.inbag = FALSE
    )
    
    # Préparation des données de test
    test <- test_data[, cov_brt]
    test[cat_vars] <- lapply(test[cat_vars], factor)
    test[num_vars] <- lapply(test[num_vars], as.numeric)
    
    # Harmoniser les niveaux des facteurs pour Lucille
    
    if (method == "Lucille") {
      for (v in cat_vars) {
        test[[v]] <- factor(test[[v]], levels = levels(train_aggr[[v]]))
      }
    }
    
    # 5. Prédiction----
  
    preds <- predict(rf_model, data = test, num.threads = kmax)$predictions
    preds <- round(preds, 2)
    
    # Affectation dans la colonne prédiction
    
    if (method == "Nicolas") {
      data_sample[[pred_col]][nblignes] <- preds
    } else {
      data_sample[[pred_col]][which(data_sample$INSEE_COM %in% test_communes)] <- preds
    }
    
    # Résultat par fold
    return(data.frame(fold = i, pred = preds, obs = test_data[[name]]))
  }
  
  # Évaluation globale
  eval <- Myeval(data_sample[[pred_col]], data_sample[[name]]) %>%
    mutate(method = paste0("RF_", substr(method, 1, 1)))
  
  return(list(
    preds = data_sample[[pred_col]],
    evaluation = eval,
    detail = resuXval,
    data_sample = data_sample 
  ))
}
