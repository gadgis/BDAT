



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


 p = tm_shape(rast(pred) ) +
  tm_raster(c("mean"),
            col.scale = tm_scale(style = "quantile",n = 10))

print(p)

writeRaster(rast(pred)[["mean"]],file="output/predKEDINLA.tif", overwrite = T)

# Validation croisée-------------
print("Validation croisée----------------")


datacov$predINLAKED = NA

  

resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
            
            # set to na to run a cross valid with inla
            # Mettre en NA les individus pour la validation crois?e
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            cmp <- activ ~ Intercept(1) +
              rfpred(qrf, model = 'linear' ) +
              field(coordinates, model = matern)
            
            
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
            
            datacov$predINLAKED[ fold[[i]] ] <-  fitted[mask] 
            

          }


resuXvalpredINLAKED <-  Myeval(datacov$predINLAKED,    datacov[,name]  )






