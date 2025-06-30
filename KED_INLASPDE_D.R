

# Validation croisée-------------
print("Validation croisée----------------")


data_sample$predINLAKED = NA

dataINLA$predRF <- data_sample$predRF


resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
          
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            cmp <- elt ~ Intercept(1) +
              rfpred(predRF, model = 'linear' ) +
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
            
            data_sample$predINLAKED[ fold[[i]] ] <-  fitted[mask] 
            

          }


res_ked  <-  Myeval(data_sample$predINLAKED,    data_sample[,name]  )






