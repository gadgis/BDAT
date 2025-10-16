



cmp <- activ ~ Intercept(1) +
  rfpred(qrf, model = 'linear' ) +
  field(coordinates, model = matern)

fitKED <- inlabru::bru(components = cmp,
                    data = dataINLA,
                    family = "Gaussian",
                    domain = list(coordinates = mesh3),
                    options = list(
                      control.inla = list(int.strategy = "eb"),
                      verbose = FALSE)
)

summary(fitKED)


fi2plot = fitKED
int.plot <- plot(fi2plot, "Intercept")
derive.plot <- plot(fi2plot, "rfpred")
spde.range <- spde.posterior(fi2plot, "field", what = "range")
spde.logvar <- spde.posterior(fi2plot, "field", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot,derive.plot)


pred <- predict(
  fitKED,
  n.samples = nsim,
  pxl,
  ~Intercept + field + rfpred  ,
  num.threads = 10
)

summary((rast(pred)))

 p = tm_shape(rast(pred) ) +
  tm_raster(c("mean"),
            col.scale = tm_scale(style = "quantile",n = 10))

print(p)

<<<<<<< HEAD:KED_INLASPDE.R
writeRaster(rast(pred)[["mean"]],file="Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKEDINLA.tif", overwrite = T)
=======
terra::writeRaster(
  rast(pred)[["mean"]],
  file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "predKEDINLA.tif"),
  overwrite = TRUE
)

>>>>>>> nsa-sauv-fin-stage:carto/KED_INLASPDE.R

# Validation croisée-------------
print("Validation croisée----------------")


datacov$predINLAKED = NA


dataINLA$predRF <- datacov$predRF


resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
            
            # set to na to run a cross valid with inla
            # Mettre en NA les individus pour la validation crois?e
            
<<<<<<< HEAD:KED_INLASPDE.R
            # dataINLA$elt <- dataINLA$activ
            # dataINLA$elt[ fold[[i]] ]  <- NA
            # 
            # cmp <- elt ~ Intercept(1) +
            #   rfpred(predRF, model = 'linear' ) +
            #   field(coordinates, model = matern)
            # 
            # 
            # Myfit <- bru(cmp,
            #              data = dataINLA,
            #              family = "Gaussian",
            #              options = list(
            #                control.inla = list(int.strategy = "eb"),
            #                verbose = FALSE)
            # )
            # 
            # # find the prediction in the output....
            # fitted <- Myfit$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            # 
            # mask <- is.na(dataINLA$elt)
            # 
            # datacov$predINLAKED[ fold[[i]] ] <-  fitted[mask] 
=======

>>>>>>> nsa-sauv-fin-stage:carto/KED_INLASPDE.R
            
            ### Avec le RF de la xva
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            cmp <- activ ~ Intercept(1) +
              rfpred(predRF, model = 'linear' ) +
              field(coordinates, model = matern)
            
            
            Myfit_KED <- bru(cmp,
                         data = dataINLA,
                         family = "Gaussian",
                         options = list(
                           control.inla = list(int.strategy = "eb"),
                           verbose = FALSE)
            )
            
            # find the prediction in the output....
            fitted <- Myfit_KED$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            
            mask <- is.na(dataINLA$elt)
            
            datacov$predINLAKED[ fold[[i]] ] <-  fitted[mask] 
            
          

          }


resuXvalpredINLAKED <-  Myeval(datacov$predINLAKED,    datacov[,name]  )


#Validation croisée  Méthode Nicolas----

<<<<<<< HEAD:KED_INLASPDE.R
#Validation croisee KED
registerDoParallel(cores = kmax) 
set.seed(123)
seg <- split(1:nrow(datacov), sample(rep(1:k, length.out = nrow(datacov))))
datacov<- datacov %>%
  mutate(INSEE_COM = as.character(INSEE_COM))

# Initialisation
datacov$predINLAKED_aggr <- NA
datacov$predINLAKED_aggr <- as.numeric(datacov$predINLAKED_aggr)

# Boucle de validation croisée
resuXval_aggr <- foreach(i = 1:k, .errorhandling = 'pass') %do% {
  
  print(i)
  
 
  # Séparer les points test et entraînement
  nblignes <- seg[[i]]
  train_data <- datacov[-nblignes, ]
  test_data  <- datacov[nblignes, ]
  
  # Agrégation de la variable cible par commune
  train_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(activ = mean(.data[[name]], na.rm = TRUE),
              pred = mean(predRF_aggr, na.rm = TRUE)) %>%
    ungroup() %>%
    left_join(centroides_communes, by = "INSEE_COM") %>%
    filter(!is.na(x) & !is.na(y))
  
  # Préparation des coordonnées des points test
  coords_test <- datacov[nblignes, c("x", "y")]
  pred_test   <- test_data$predRF_aggr
  
  # Création du data.frame d'entrée pour bru
  dataINLA_aggr <- data.frame(
    x = c(train_aggr$x, coords_test$x),
    y = c(train_aggr$y, coords_test$y),
    elt = c(train_aggr$activ, rep(NA, nrow(coords_test))),
    pred=c(train_aggr$pred, pred_test)
  )
  
  coordinates(dataINLA_aggr) <- NomsCoord
  proj4string(dataINLA_aggr) <-  "epsg:2154"
  
  # Formule du modèle
  cmp <- elt ~ Intercept(1) +
    rfpred(pred, model="linear")+
    field(coordinates,
          model = matern)
  
  # Calibration du modèle
  Myfit_KED <- bru(
    cmp,
    data = dataINLA_aggr,
    family = "Gaussian",
    options = list(
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
  
  # Prédictions sur tous les points
  fitted_KED <- Myfit_KED$summary.fitted.values$mean[1:nrow(dataINLA_aggr)]
  
  # Récupération des prédictions sur les points test uniquement
  mask <- is.na(dataINLA_aggr$elt)
  datacov$predINLAKED_aggr[nblignes] <- fitted_KED[mask]
  
  return(data.frame(seg = i, pred = fitted_KED[mask]))
}

resuXvalINLAKED_aggr <- Myeval(datacov$predINLAKED_aggr, datacov[[name]])

#Validation croisée méthode Lucille
#Validation croisee KED

set.seed(123)
seg <- split(1:nrow(datacov), sample(rep(1:k, length.out = nrow(datacov))))
datacov<- datacov %>%
  mutate(INSEE_COM = as.character(INSEE_COM))

communes_unique <- unique(datacov$INSEE_COM)
commune_folds <- split(communes_unique, sample(rep(1:k, length.out = length(communes_unique))))

# Initialisation
datacov$predINLAKED_aggrcom <- NA
datacov$predINLAKED_aggrcom <- as.numeric(datacov$predINLAKED_aggrcom)

resuXval_aggrcom<- foreach(i = 1:k, .errorhandling = 'pass')  %do% {
  
  cat("seg", i, "\n")
  
  # Communes test
  test_communes <- commune_folds[[i]]
  
  # Séparer entraînement et test
  train_data <- datacov[!datacov$INSEE_COM %in% test_communes, ]
  test_data  <- datacov[datacov$INSEE_COM %in% test_communes, ]
  
  # Agrégation des données par commune pour l'entraînement
  train_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(
      activ = mean(.data[[name]], na.rm = TRUE),
      pred = mean(predRF_aggrcom, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    left_join(centroides_communes, by = "INSEE_COM") %>%
    filter(!is.na(x) & !is.na(y) & !is.na(pred))
  
  # Préparation des prédicteurs pour les points test
  coords_test <- test_data[, c("x", "y")]
  pred_test   <- test_data$predRF_aggrcom
  
  # Construction du SpatialPointsDataFrame
  dataINLA_aggr <- data.frame(
    x = c(train_aggr$x, coords_test$x),
    y = c(train_aggr$y, coords_test$y),
    elt = c(train_aggr$activ, rep(NA, nrow(coords_test))),
    pred = c(train_aggr$pred, pred_test)
  )
  
  coordinates(dataINLA_aggr) <- NomsCoord
  proj4string(dataINLA_aggr) <- CRS("epsg:2154")
  
  # Formule KED
  cmp <- elt ~ Intercept(1) +
    rfpred(pred, model = "linear") +
    field(coordinates, model = matern)
  
  # Calibration
  Myfit_KED <- bru(
    cmp,
    data = dataINLA_aggr,
    family = "Gaussian",
    options = list(
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
  
  # Prédictions
  fitted_KED <- Myfit_KED$summary.fitted.values$mean[1:nrow(dataINLA_aggr)]
  mask <- is.na(dataINLA_aggr$elt)
  datacov$predINLAKED_aggrcom[which(datacov$INSEE_COM %in% test_communes)] <- fitted_KED[mask]
  
  return(data.frame(fold = i, pred = fitted_KED[mask], obs = test_data[[name]]))
}

# Fermeture du cluster
stopImplicitCluster()

# Évaluation
resuXvalINLAKED_aggr_com <- Myeval(datacov$predINLAKED_aggrcom, datacov[[name]])
=======
>>>>>>> nsa-sauv-fin-stage:carto/KED_INLASPDE.R
