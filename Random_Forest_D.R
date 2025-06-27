library(sf) # manipuler les datas SIG vector
library(tidyr) # reroganisation des données
library(dplyr)  # grammaire de manipulation de données
library(terra) # manipuler des données raster-pixels
library(ggplot2) # graphiques
library(sp) # ancetre de sf. Pourquoi l'utiliser ?
library(randomForest) # randomForest pour utiliser tuneRandomForest
library(Boruta)  # implementation  de selection de variables de randomForest
library(ranger) # quantile randomForest, version optimée de randomForest
library(caret) # modelisation creation de modeles prédictifs avec tuning validation-croisée
library(foreach) # boucles optimisées de R plus facile à coder voir option combine
library(raster) # manipuler des données raster
library(purrr)
library(ggpubr) # pour les graphiques
library(doParallel)
library(raster) 
library(tuneRanger)
library(iml)
library(mlr)

## Fonction utilsée pour la validation --------------

Myeval <- function(x, y){
  
  # mean error
  ME <- round(mean(y - x, na.rm = TRUE), digits = 2)
  
  # root mean square error
  RMSE <-   round(sqrt(mean((y - x)^2, na.rm = TRUE)), digits = 2)
  
  # root mean absolute error
  MAE <-   round(mean(abs(y - x), na.rm = TRUE), digits = 2)
  
  # Pearson's correlation squared
  r2 <-  round((cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2), digits = 2)
  
  # MEC
  SSE <- sum((y - x) ^ 2, na.rm = T)
  SST <- sum((y - mean(y, na.rm = T)) ^ 2, na.rm = T)
  NSE <- round((1 - SSE/SST), digits = 2)
  
  # concordance correlation coefficient
  n <- length(x)
  sdx <- sd(x, na.rm = T)
  sdy <- sd(y, na.rm = T)
  r <- stats::cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')
  # scale shift
  v <- sdx / sdy
  sx2 <- var(x, na.rm = T) * (n - 1) / n
  sy2 <- var(y, na.rm = T) * (n - 1) / n
  # location shift relative to scale
  u <- (mean(x, na.rm = T) - mean(y, na.rm = T)) / ((sx2 * sy2)^0.25)
  Cb <- ((v + 1 / v + u^2)/2)^-1
  rCb <- r * Cb
  rhoC <- round(rCb, digits = 2)
  
  Cb <- round(Cb, digits = 2)
  r <- round(r, digits = 2)
  
  # return the results
  evalRes <- data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE, rhoC = rhoC, Cb = Cb)
  
  return(evalRes)
}



## Définition des variables ------

name="pH"
kmax= 23 # pour la parallelisation, le nb de coeurs
ntree = 350 # le nbre d'arbre de random forest
sample_sizes <- c(600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 3000, 4000, 5000,6000,7000,8000,9000,10000)
repets <- 100
NomsCoord <- c("x","y")


# 1 Preparation des données pour la spatialisation
#Extraction des matrices de covariables pour les données ponctuelles----

chemin_cov<- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"
setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

##liste des fichiers de covariables----
l<- list.files(chemin_cov, pattern = ".tif$", full.names = TRUE)

st <- rast(l)

communes <- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/COMMUNE.shp")


dtTB <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")

cov_brt<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_cov_brt.rds")


# attribution d'un id par sites (pour les doublons analytique possible)

dtTB <- dtTB %>%
  mutate(
    identifiant = ifelse(source == "IGCS", id_profil, bdatid)
  )

datacov <- terra::extract( st , 
                           dtTB %>%
                             st_as_sf(coords = NomsCoord ,
                                      crs = 2154)
) %>% 
  bind_cols(dtTB %>%
              dplyr::select(all_of(c( name,NomsCoord,"INSEE_COM","source","identifiant")  )
              )
  ) %>%
  
  mutate(id = row_number()) %>%
  na.omit()


# 3 Modélisation par RF -----
#Création du tableau des covariables séléctionnées

cov_brt = c(name, cov_brt)


# Création du cluster parallèle
registerDoParallel(cores = parallel::detectCores() - 1)

results_all <- list()
results_preds <- list()

for (n in sample_sizes) {
  cat("Échantillon de taille :", n, "\n")
  
  # Déterminer le nombre de folds 
  indiv_per_fold <- ifelse(n <= 2000, 200, 500)
  k <- ceiling(n / indiv_per_fold)
  
  tuned_params <- NULL
  
  resu_rep <- foreach(rep = 1:repets, 
                      
                      .packages = c("randomForest", "dplyr", "ranger", "tuneRanger", "iml", "mlr"),
                      .export = c("Myeval", "name", "cov_brt", "datacov", "k", "ntree", "kmax", "NomsCoord"))%dopar% {
    set.seed(rep)
                                                                        
    datacov_sample <- datacov %>% sample_n(n)
    datacov_sample_shrt <- datacov_sample[, cov_brt]
    
    
    #tune des hyperparamètres
    if (rep == 1) {
      SOC.task <- makeRegrTask(data = datacov_sample_shrt, target = name)
      res <- tuneRanger(SOC.task, num.trees = ntree, iters = 100, num.threads = 1)
      tuned_params <- res$recommended.pars
    }
    
    
    # Validation croisée
    datacov_sample$id <- 1:nrow(datacov_sample)
    
    
    folds <- split(datacov_sample$id, rep(1:k, length.out = n))
  
    
    datacov_sample$predRF <- NA
    
    print("Validation croisée----------------")
    
    #initialisation
    datacov_sample$predRF <- NA
    
    for (i in 1:k) {
      print(i)
      
      nblignes = which( datacov_sample$id %in% datacov_sample$id[ folds[[i]] ] )
    
      #Calibre le modèle
      fomula.ranger <- as.formula(paste0(name,"~."))
      RF_Mod.G <- ranger(formula = fomula.ranger ,
                         data = datacov_sample_shrt[-nblignes, ],
                         num.trees = ntree,
                         min.node.size = tuned_params$min.node.size,
                         quantreg = F,
                         max.depth = 15, 
                         mtry=tuned_params$mtry ,
                         importance="permutation", 
                         scale.permutation.importance = TRUE, 
                         keep.inbag = F)
      
     
      # Prédiction
      datacov_sample$predRF[ nblignes ] <- predict(RF_Mod.G,
                                                   datacov_sample_shrt[ nblignes , ],
                                                   num.threads = kmax )$predictions
     
    }
    
    metrics_df  <- Myeval(datacov_sample$predRF, datacov_sample[[name]]) %>%
      mutate(rep = rep, sample_size = n)
    
    pred_df <- datacov_sample %>%
      select(identifiant, all_of(NomsCoord),predRF) %>%
      mutate(rep = rep, sample_size = n)
    
    list(metrics = metrics_df, preds = pred_df)
                      }
  
  results_all[[as.character(n)]] <- resu_rep
}

stopImplicitCluster()

results_metrics <- bind_rows(lapply(results_all, function(l) bind_rows(lapply(l, `[[`, "metrics"))))
results_preds   <- bind_rows(lapply(results_all, function(l) bind_rows(lapply(l, `[[`, "preds"))))

results_summary <- results_metrics %>%
  group_by(sample_size) %>%
  summarise(across(c(ME, MAE, RMSE, r2, NSE), mean, na.rm = TRUE))

saveRDS(results_summary, "resultats_RF_Degrad.rds")

saveRDS(results_preds,   "predictions_RF_Degradation.rds")
