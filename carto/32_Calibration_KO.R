#===================================================================================================================#
# Script :      Krigeage Ordinaire avec INLA SPDE

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la prédiction des propriétés des sols en utilisant 
#               la méthode Krigeage Ordinaire avec INLA SPDE 

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Observations ponctuelles avec la localisation (x,y) sous forme sp

# Sorties :     Distribution a posteriori des paramètres du modèles, valeurs de la propriété cible sur la grille
#               de prédictions et les indicateurs de performance

# Modification : 06-10-2025
#===================================================================================================================#

#================================================DEBUT DU SCRIPT====================================================#

# il faut récupérer une licence sur Pardiso : inla.pardiso()
#inla.setOption(pardiso.license ="82F70E96BE7DA6A5956D4DF8F31E127ACCB33C981DE83F430BA469A2")

inla.setOption(
  num.threads = 15 ,
  inla.mode="experimental"
)

# 1.Définir de la grille de prédictions-----

pxl <- as.data.frame(r,xy=T)
colnames(pxl)[3] <- "qrf"
gridded(pxl) <- ~x+y


# 2. Préparation du mesh (maillage) pour la discrétisation du champ spatial----

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


#3.Définir le champ gaussien Matérn----

matern <-
  INLA::inla.spde2.pcmatern(mesh3,
                            alpha = 2,
                            prior.sigma = c(10, 0.01),# P(sigma > 1) = 0.5
                            prior.range = c(500, 0.01)  # P(range < 100000 m) = 0.9
  )
#4.Définir et ajuster le modèle----

cmp <- activ ~ Intercept(1) +
  field(coordinates, model = matern)

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

#5. Diagnostic de la distribution a posteriori----

spde.range <- spde.posterior(fi2plot, "field", what = "range")
spde.logvar <- spde.posterior(fi2plot, "field", what = "log.variance")

int.plot <- plot(fi2plot, "Intercept")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot)


# 6. Prédictions sur la grille pxl-----

predKO <- predict(
  fitKO,
  n.samples = nsim,
  pxl,
  ~Intercept + field   ,
  num.threads = 10
)


#7. Visualisation et sauvegarde-----

p = tm_shape(rast(predKO) ) +
  tm_raster(c("mean"),
            col.scale = tm_scale(style = "quantile",n = 10))


print(p)


terra::writeRaster(
  rast(predKO)[["mean"]],
  file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "predKOINLA.tif"),
  overwrite = TRUE
)


# 8.Validation croisée-----

print("Validation croisée----------------")

#initialisation 

datacov$predINLAKO = NA

resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
            
            # Mettre en NA les individus pour la validation croisée
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            
            #calibration du modèle
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
            
            # predition sur le groupe mis à NA
            fitted <- Myfit_KO$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            
            mask <- is.na(dataINLA$elt)
            
            datacov$predINLAKO[ fold[[i]] ] <-  fitted[mask] 
            
            fitted[mask] 
          }

#Calcul des indicateurs de performance
resuXvalTKO <-  Myeval(datacov$predINLAKO,   datacov[,name] )

#================================================FIN DU SCRIPT====================================================#
