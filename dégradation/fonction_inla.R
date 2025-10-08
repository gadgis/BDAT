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
  
  # pred_col <- paste0("predINLA", type, "_", substr(approach, 1, 1), substr(type_val, 1, 1))
  
  # id_point_vect <- dataINLA$id_point
  
  dataINLA$elt <- dataINLA[,name]
  
  coords <- dataINLA[, NomsCoord]
  
  coordinates(dataINLA) <- NomsCoord
  proj4string(dataINLA) <- CRS("epsg:2154")
  

  
  coordinates(data_test) <- NomsCoord
  proj4string(data_test) <- CRS("epsg:2154")

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
  
  
    fit <- bru(
    cmp,
    data = dataINLA,
    family = "Gaussian",
    
    options = list(control.inla = list(int.strategy = "eb"), 
                   verbose = FALSE)
  )
  
  gc()
  
  # # Sécurisation de la jointure
  # data_test_unique <- data_test[, c("id", name)] %>% distinct(id, .keep_all = TRUE)
  # 
  # fitted_df <- data.frame(
  #   id = id_point_vect,
  #   pred = fit$summary.fitted.values$mean[1:length(id_point_vect)]
  # ) %>%
  #   filter(!is.na(id)) %>%
  #   left_join(data_test_unique, by = "id") %>%
  #   rename(obs = !!name)
  # 
  

   predObj <-  predict(fit,
                       data_test ,
                      cmpPred  
                       
                     )
  

  return(preds = predObj$median)
}