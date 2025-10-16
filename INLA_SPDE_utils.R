#Script pour la factorisation du KED et du KO

run_inla_spde <- function(type = c("KO", "KED"),
                          aggregation = c("Nicolas", "Lucille"),
                          data_sample,
                          centroides_communes,
                          fold,
                          name,
                          k,
                          NomsCoord = c("x", "y")) {
  
  type <- match.arg(type)
  aggregation <- match.arg(aggregation)
  pred_col <- paste0("predINLA", type, "_", substr(aggregation, 1, 1))
  
  # 1. Définir la colonne de prédiction RF à utiliser----
  
  pred_rf_col <- if (aggregation == "Nicolas") "predRF_N" else "predRF_L"
  
  # Ajouter la colonne préd_rf_col si elle n'existe pas
  if (!pred_rf_col %in% names(data_sample)) {
    data_sample[[pred_rf_col]] <- NA_real_
  }
  
 #2. Initialiser la colonne de prédiction INLA----
  
  data_sample[[pred_col]] <- NA_real_
  
  resuXval <- foreach(i = 1:k, .combine = rbind, .packages = c("dplyr", "sp", "inlabru")) %do% {
    
#3. Séparation train/test----
    
    if (aggregation == "Nicolas") {
      nblignes <- fold[[i]]
      train_data <- data_sample[-nblignes, ]
      test_data  <- data_sample[nblignes, ]
    } else {
      test_communes <- fold[[i]]
      train_data <- data_sample[!data_sample$INSEE_COM %in% test_communes, ]
      test_data  <- data_sample[data_sample$INSEE_COM %in% test_communes, ]
    }
    
#4. Agrégation----
    
    train_aggr <- train_data %>%
      group_by(INSEE_COM) %>%
      summarise(
        activ = mean(.data[[name]], na.rm = TRUE),
        pred = if (type == "KED") mean(.data[[pred_rf_col]], na.rm = TRUE) else NA_real_,
        .groups = "drop"
      ) %>%
      left_join(centroides_communes, by = "INSEE_COM") %>%
      filter(!is.na(X), !is.na(Y),!is.na(activ), if (type == "KED") !is.na(pred) else TRUE)
    
    coords_test <- test_data[, NomsCoord]
    elt_test <- rep(NA_real_, nrow(coords_test))
    
#5. Préparation données pour INLA----
    
    if (type == "KO") {
      dataINLA <- data.frame(
        x = c(train_aggr$X, coords_test$x),
        y = c(train_aggr$Y, coords_test$y),
        elt = c(train_aggr$activ, elt_test)
      )
    } else {
      pred_test <- test_data[[pred_rf_col]]
      dataINLA <- data.frame(
        x = c(train_aggr$X, coords_test$x),
        y = c(train_aggr$Y, coords_test$y),
        elt = c(train_aggr$activ, elt_test),
        pred = c(train_aggr$pred, pred_test)
      )
    }
    
    
#6. Création du mesh----
    
    coords <- dataINLA[, NomsCoord]
    max.edge <- diff(range(coords[,1]))/(3*5)
    bound.outer <- diff(range(coords[,1])) / 3
    mesh3 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
                          offset = c(max.edge, bound.outer), cutoff = max.edge/10)
    matern <- INLA::inla.spde2.pcmatern(mesh3, alpha = 2,
                                        prior.sigma = c(10, 0.01),
                                        prior.range = c(500, 0.01))
    
    dataINLA$id_point <- c(rep(NA, nrow(train_aggr)), test_data$id)
    id_point_vect <- dataINLA$id_point
    
    coordinates(dataINLA) <- NomsCoord
    proj4string(dataINLA) <- CRS("epsg:2154")
    
    cmp <- if (type == "KO") {
      elt ~ Intercept(1) + field(coordinates, model = matern)
    } else {
      elt ~ Intercept(1) + rfpred(pred, model = "linear") + field(coordinates, model = matern)
    }
    
    fit <- bru(
      cmp,
      data = dataINLA,
      family = "Gaussian",
      options = list(control.inla = list(int.strategy = "eb"), verbose = FALSE)
    )
    
#7. Prédictions et appariement par id_point----
    
    fitted_df <- data.frame(
      id_point = id_point_vect,
      pred = fit$summary.fitted.values$mean[1:length(dataINLA$activ)]
    ) %>%
      filter(!is.na(id_point))
    
    data_sample[[pred_col]][match(fitted_df$id_point, data_sample$id)] <- fitted_df$pred
    
    return(data.frame(fold = i, pred = fitted_df$pred, obs = test_data[[name]]))
  }
  
#8. Évaluation globale----
  
  eval <- Myeval(data_sample[[pred_col]], data_sample[[name]]) %>%
    mutate(method = paste0("INLA", type, "_", substr(aggregation, 1, 1)))
  
  return(list(preds = data_sample[[pred_col]], evaluation = eval, detail = resuXval))
}
