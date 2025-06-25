library(INLA)
library(inlabru)


# INLA -------------
#inla.setOption(pardiso.license ="82F70E96BE7DA6A5956D4DF8F31E127ACCB33C981DE83F430BA469A2")
inla.setOption(
  num.threads = 15 ,
  inla.mode="experimental"
)
# il faut récupérer une licence sur Pardiso : inla.pardiso()






# Grille de prédictions
# transformation en sp depuis terra
pxl <- as.data.frame(r,xy=T)
colnames(pxl)[3] <- "qrf"
gridded(pxl) <- ~x+y




# Préparation du mesh pour la méthode INLA pour la discrétisation


# max.edge = diff(range(coords[,1]))/(3*5)
# bound.outer = diff(range(range(coords[,1])))/3
# 
# bndint <- inla.nonconvex.hull( points = dataINLA, convex=-.05)
# bndext <- inla.nonconvex.hull(points = dataINLA, convex=-.2)
# 
# # Use of inla.mesh.2d 
# prmesh = inla.mesh.2d(loc = coords,
#                       boundary = list(int = limit_zone,
#                                       out = bndext),
#                       max.edge = c(.5,2)*max.edge, 
#                       cutoff = cutoffValue,
#                       crs = crs(limit_zone)
# )
# 
# 
# px_mesh <- fm_mesh_2d_inla(
#   loc = dataINLA,
#   boundary = limit_zone,
#   max.edge = c(.5,2)*max.edge,
#   crs = st_crs(limit_zone)
# )


max.edge = diff(range(coords[,1]))/(3*5) # taille des triangles au cœur du maillage

bound.outer = diff(range(coords[,1]) ) / 3 # marge externe où les triangles deviennent plus grands


mesh3 = inla.mesh.2d(loc=coords,        # noeuds forcés aux points de mesure
                     max.edge = c(1,2)*max.edge,
                     offset=c(max.edge, bound.outer),
                     cutoff = max.edge/10)

ggplot() +
  # geom_sf(data=limit_zone,
  #         color='turquoise',fill='transparent')+  
  gg(mesh3) +
  geom_sf(data=as(dataINLA,"sf"),
          col='purple',
          size=1.7,alpha=0.5) 


# ggplot() +
#   geom_fm(data = px_mesh) +
#   geom_sf(
#     data = counts_df[counts_df$count > 0, ],
#     aes(color = count),
#     size = 1,
#     pch = 4
#   ) +
#   theme_minimal()
# # 
# ggplot() +
#   gg(mesh3) +
#   gg( dataINLA ) +
#   gg( dataINLA1,col=2 ) +
#   coord_equal()
# 
# # 
# mesh <- fm_mesh_2d_inla(loc=coords,boundary = limit_zone, max.edge = 50)
# ggplot() +
#   geom_fm(data = mesh)
# 

# Modelling inla

##Définit le champ gaussien Matérn
matern <-
  INLA::inla.spde2.pcmatern(mesh3,
                            alpha = 2,
                            prior.sigma = c(10, 0.01),# P(sigma > 1) = 0.5
                            prior.range = c(500, 0.01)  # P(range < 100000 m) = 0.9
  )
##Définit le modèle
cmp <- activ ~ Intercept(1) +
  field(coordinates, model = matern)

#Ajuste le modèle
fitKO <- inlabru::bru(components = cmp,
                      data = dataINLA,
                      family = "Gaussian",
                      domain = list(coordinates = mesh3),
                      options = list(
                        control.inla = list(int.strategy = "eb"),
                        verbose = FALSE)
)


summary(fitKO)


fi2plot = fitKO
#diagnostic des postrior
spde.range <- spde.posterior(fi2plot, "field", what = "range")
spde.logvar <- spde.posterior(fi2plot, "field", what = "log.variance")

int.plot <- plot(fi2plot, "Intercept")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot)



predKO <- predict(
  fitKO,
  n.samples = nsim,
  pxl,
  ~Intercept + field   ,
  num.threads = 10
)



p = tm_shape(rast(predKO) ) +
  tm_raster(c("mean"),
            col.scale = tm_scale(style = "quantile",n = 10))


print(p)

terra::writeRaster(rast(predKO)[["mean"]],file="Y:/BDAT/traitement_donnees/MameGadiaga/resultats/predKOINLA.tif", overwrite = T)

# Validation croisée-------------
print("Validation croisée----------------")


datacov$predINLAKO = NA

resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
            
            # set to na to run a cross valid with inla
            # Mettre en NA les individus pour la validation crois?e
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            
            
           cmp <- elt ~ Intercept(1) +
              field(coordinates,
                    model = matern)
            
            Myfit_KO <- bru(cmp,
                         data = dataINLA,
                         family = "Gaussian",
                         options = list(
                           control.inla = list(int.strategy = "eb"),
                           verbose = FALSE)
            )
            
            # find the prediction in the output....
            fitted <- Myfit_KO$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            
            mask <- is.na(dataINLA$elt)
            
            datacov$predINLAKO[ fold[[i]] ] <-  fitted[mask] 
            
            fitted[mask] 
          }


resuXvalTKO <-  Myeval(datacov$predINLAKO,   datacov[,name] )


#Validation croisée  Méthode Nicolas----

# Validation croisee KO
registerDoParallel(cores = kmax) 
set.seed(123)
seg <- split(1:nrow(datacov), sample(rep(1:k, length.out = nrow(datacov))))
datacov<- datacov %>%
  mutate(INSEE_COM = as.character(INSEE_COM))
# Initialisation
datacov$predINLAKO_aggr <- NA
datacov$predINLAKO_aggr <- as.numeric(datacov$predINLAKO_aggr)

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
    summarise(activ = mean(.data[[name]], na.rm = TRUE)) %>%
    ungroup() %>%
    left_join(centroides_communes, by = "INSEE_COM") %>%
    filter(!is.na(x) & !is.na(y))
  
  # Préparation des coordonnées des points test
  coords_test <- datacov[nblignes, c("x", "y")]
  
  # Création du data.frame d'entrée pour bru
  dataINLA_aggr <- data.frame(
    x = c(train_aggr$x, coords_test$x),
    y = c(train_aggr$y, coords_test$y),
    elt = c(train_aggr$activ, rep(NA, nrow(coords_test)))
  )
  
  coordinates(dataINLA_aggr) <- NomsCoord
  proj4string(dataINLA_aggr) <-  "epsg:2154"
  
  # Formule du modèle
  cmp <- elt ~ Intercept(1) +
    field(coordinates,
          model = matern)
  
  # Calibration du modèle
  Myfit_aggr <- bru(
    cmp,
    data = dataINLA_aggr,
    family = "Gaussian",
    options = list(
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
  
  # Prédictions sur tous les points
  fitted_aggr <- Myfit_aggr$summary.fitted.values$mean[1:nrow(dataINLA_aggr)]
  
  # Récupération des prédictions sur les points test uniquement
  mask <- is.na(dataINLA_aggr$elt)
  datacov$predINLAKO_aggr[nblignes] <- fitted_aggr[mask]
  
  return(data.frame(seg = i, pred = fitted_aggr[mask]))
}

resuXvalINLAKO_aggr <- Myeval(datacov$predINLAKO_aggr, datacov[[name]])


#Validation croisée Méthode Lucille----
# Liste unique des communes
communes_unique <- unique(datacov$INSEE_COM)

# Créer les folds au niveau des communes
set.seed(123)
commune_folds <- split(communes_unique, sample(rep(1:k, length.out = length(communes_unique))))

#initialisation

datacov$predINLAKO_aggr_com <- NA  

resuXval_aggr_com <- foreach(i = 1:k, .errorhandling = 'pass') %do% {
  
  print(paste("Fold", i))
  
  # Identifier les communes du fold test
  test_communes <- commune_folds[[i]]
  train_data <- datacov[!datacov$INSEE_COM %in% test_communes, ]
  test_data  <- datacov[datacov$INSEE_COM %in% test_communes, ]
  
  # Agréger les données d'entraînement par commune
  train_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(activ = mean(.data[[name]], na.rm = TRUE)) %>%
    ungroup() %>%
    left_join(centroides_communes, by = "INSEE_COM") %>%
    filter(!is.na(x) & !is.na(y))
  
  # Préparation des coordonnées des points test
  coords_test <- test_data[, c("x", "y")]
  
  # Création du data.frame d'entrée pour bru
  dataINLA_aggr <- data.frame(
    x = c(train_aggr$x, coords_test$x),
    y = c(train_aggr$y, coords_test$y),
    elt = c(train_aggr$activ, rep(NA, nrow(coords_test)))
  )
  
  coordinates(dataINLA_aggr) <- NomsCoord
  proj4string(dataINLA_aggr) <- CRS("epsg:2154")
  
  # Formule du modèle
  cmp <- elt ~ Intercept(1) +
    field(coordinates, model = matern)
  
  # Calibration du modèle
  Myfit_aggr <- bru(
    cmp,
    data = dataINLA_aggr,
    family = "Gaussian",
    options = list(
      control.inla = list(int.strategy = "eb"),
      verbose = FALSE
    )
  )
  
  # Prédictions sur tous les points
  fitted_aggr <- Myfit_aggr$summary.fitted.values$mean[1:nrow(dataINLA_aggr)]
  
  # Récupération des prédictions sur les points test uniquement
  mask <- is.na(dataINLA_aggr$elt)
  datacov$predINLAKO_aggr_com[which(datacov$INSEE_COM %in% test_communes)] <- fitted_aggr[mask]
  
  return(data.frame(fold = i, pred = fitted_aggr[mask], obs = test_data[[name]]))
}

# Évaluation
resuXvalINLAKO_aggr_com <- Myeval(datacov$predINLAKO_aggr_com, datacov[[name]])
print(resuXvalINLAKO_aggr_com)

stopImplicitCluster()