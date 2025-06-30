



# Préparation du mesh pour la méthode INLA pour la discrétisation

max.edge = diff(range(coords[,1]))/(3*5) # taille des triangles au cœur du maillage

bound.outer = diff(range(coords[,1]) ) / 3 # marge externe où les triangles deviennent plus grands


mesh3 = inla.mesh.2d(loc=coords,        # noeuds forcés aux points de mesure
                     max.edge = c(1,2)*max.edge,
                     offset=c(max.edge, bound.outer),
                     cutoff = max.edge/10)

# Préparation des données pour inlabru
matern <-
  INLA::inla.spde2.pcmatern(mesh3,
                            alpha = 2,
                            prior.sigma = c(10, 0.01),# P(sigma > 1) = 0.5
                            prior.range = c(500, 0.01)  # P(range < 100000 m) = 0.9
  )


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
            
            Myfit <- bru(cmp,
                         data = dataINLA,
                         family = "Gaussian",
                         options = list(
                           control.inla = list(int.strategy = "eb"),
                           verbose = FALSE)
            )
            
            # find the prediction in the output....
            fitted <- Myfit$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            
            mask <- is.na(dataINLA$elt)
            
            data_sample$predINLAKO[ fold[[i]] ] <-  fitted[mask] 
            
            fitted[mask] 
          }


res_ko  <-  Myeval(data_sample$predINLAKO,   data_sample[,name] )






