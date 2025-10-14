#===================================================================================================================#
# Script :      Krigeage avec dérive externe avec INLA SPDE

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la prédiction des propriétés des sols en utilisant 
#               la méthode Krigeage avec dérive externe avec INLA SPDE 

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Observations ponctuelles avec la localisation (x,y) et le prédictions du Random forest 

# Sorties :     Distribution a posteriori des paramètres du modèles, valeurs de la propriété cible sur la grille
#               de prédictions et les indicateurs de performance

# Modification : 06-10-2025
#===================================================================================================================#

#================================================DEBUT DU SCRIPT====================================================#

#1. Définir et ajuster le modèle-----

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

#2. Visualiser les résultats du modèle----

fi2plot = fitKED
int.plot <- plot(fi2plot, "Intercept")
derive.plot <- plot(fi2plot, "rfpred")
spde.range <- spde.posterior(fi2plot, "field", what = "range")
spde.logvar <- spde.posterior(fi2plot, "field", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot,derive.plot)

#3. Prédiction et visualisation -----

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

#4. Sauvegarde des résultats-----

terra::writeRaster(
  rast(pred)[["mean"]],
  file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "predKEDINLA_geomask",d,".tif"),
  overwrite = TRUE
)


# Validation croisée------


print("Validation croisée----------------")

#initialisation
datacov$predINLAKED = NA


dataINLA$predRF <- datacov$predRF


resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
            
            # Mettre en NA les individus pour la validation croisée
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            # Calibration du modèle avec les prédiction du RF sur les k-1 groupes
            
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
            
            # prédiction sur le groupe restant
            
            fitted <- Myfit_KED$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            
            mask <- is.na(dataINLA$elt)
            
            datacov$predINLAKED[ fold[[i]] ] <-  fitted[mask] 
            
          

          }

#Calcul des indicateurs de performance
resuXvalpredINLAKED <-  Myeval(datacov$predINLAKED,    datacov[,name]  )


#================================================FIN DU SCRIPT====================================================#
