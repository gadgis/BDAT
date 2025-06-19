
# RandomForest Adapté au protocole de dégradation ----

library(ranger)
library(mlr)
library(tuneRanger)
library(dplyr)
library(foreach)
library(doParallel)


cores <- kmax
cl <- makeCluster(cores)
registerDoParallel(cl)

clusterExport(cl, varlist = c(
  "datacov", "rmqs", "st", "name", "NomsCoord",
  "cov_brt", "kmax", "ntree", "idcovs", "idvar", "Myeval","repets","sample_sizes","results_all"
))


for (n in sample_sizes) {
  cat("Taille d'échantillon :", n, "\n")
  
  repeat_results <- foreach(rep = 1:repets, .combine = rbind, .packages = c("dplyr", "tidyr", "foreach", "randomForest", "INLA", "inlabru", "terra")) %dopar% {
    
    # Sous-échantillonnage aléatoire
    set.seed(rep)
    datacov_sample <- datacov %>%
      sample_n(n)
    
    datacov_sample$id <- 1:nrow(datacov_sample)
    
    # Définition des folds pour validation croisée interne
    k <- ifelse(n <= 2000, ceiling(n / 200), ceiling(n / 500))
    fold <- split(1:nrow(datacov_sample), rep(1:k, length.out = nrow(datacov_sample)))
    
    datacov_shrt <- datacov_sample[, cov_brt]
    
    # Création de la tâche ML
    SOC.task <- makeRegrTask(data = datacov_shrt, target = name)
    
    # Tuning des paramètres
    res <- tuneRanger(SOC.task, num.trees = ntree, iters = 100, num.threads = kmax)
    mtry_best <- res$recommended.pars$mtry
    min_nodesize_best <- res$recommended.pars$min.node.size
    
    # Validation croisée 
    
    # Initialisation
    datacov_sample$predRF <- NA
    
    resuXval <- foreach(i = 1:k, .combine = rbind, .packages = c("ranger", "dplyr")) %dopar% {
      
      nblignes <- which(datacov_sample$id %in% datacov_sample$id[ fold[[i]] ])
      RF_Mod.G <- ranger(
        formula = as.formula(paste0(name, "~ .")),
        data = datacov_shrt[-nblignes, ],
        num.trees = ntree,
        mtry = mtry_best,
        min.node.size = min_nodesize_best,
        quantreg = FALSE,
        max.depth = 15,
        importance = "permutation",
        scale.permutation.importance = FALSE,
        keep.inbag = FALSE
      )
      
      preds <- predict(RF_Mod.G, datacov_shrt[nblignes, ])$predictions
      
      datacov_sample$predRF[nblignes] <- preds
      
      return(data.frame(fold = i, pred = preds))
    }
    
  
    
    
    # Validation externe sur le RF avec RMQS 
    
    # Extraction des covariables au niveau des points RMQS
    rmqs_mod <- rmqs 
    
    vext <- terra::extract(st,
                           rmqs_mod %>%
                             st_as_sf(coords = NomsCoord, crs = 2154)) %>%
      bind_cols(rmqs_mod %>% dplyr::select(all_of(c(name, NomsCoord)))) %>%
      na.omit()
    
    # Préparation des covariables pour prédiction
    datacov_shrt_vext <- vext[cov_brt]  
    
    test <- datacov_shrt_vext %>%
      dplyr::select(all_of(cov_brt[-1]))  
    
    # Prédiction avec le modèle RF ajusté 
    datacov_shrt_vext$predRF <- predict(RF_Mod.G,
                                        test,
                                        num.threads = kmax)$predictions
    
    # Arrondir + calculer les métriques
    datacov_shrt_vext <- datacov_shrt_vext %>%
      mutate(predRF = round(predRF, 2))
    
    met_RF_ext <- Myeval(datacov_shrt_vext$predRF, datacov_shrt_vext[[name]])
    
    
    # Calcul des métriques
    met_RF <- Myeval(datacov_sample$predRF, datacov_sample[[name]])
    
    
    res_interne <- cbind(
      repet= rep,
      sample_size = n,
      model = "RF",
      type = "interne",
      met_RF
    )
    
    res_externe <- cbind(
      repet = rep,
      sample_size = n,
      model = "RF",
      type = "externe",
      met_RF_ext
    )
    
    return(rbind(res_interne, res_externe))
    
  }
  
  results_all[[as.character(n)]] <- repeat_results
}

stopCluster(cl)

# Combiner tous les résultats
results_df <- do.call(rbind, results_all)



