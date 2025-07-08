library(sf)
library(tmap)
library(readxl)
library(tidyr)
library(dplyr)
library(terra)
library(ggplot2)
library(sp)
library(randomForest)
library(Boruta)
library(ranger)
library(caret)
library(foreach)
library(raster)
library(purrr)
library(ggpubr)
library(ggnewscale)
library(doParallel)
library(tuneRanger)
library(INLA)
library(inlabru)
library(iml)
library(mlr)

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/Codes R")

Myeval <- function(x, y){
  
  # mean error
  ME <- mean(y - x, na.rm = TRUE)
  
  # root mean square error
  RMSE <-   sqrt(mean((y - x)^2, na.rm = TRUE))
  
  # root mean absolute error
  MAE <-   mean(abs(y - x), na.rm = TRUE)
  
  # Pearson's correlation squared
  r2 <-  (cor(x, y, method = 'pearson', use = 'pairwise.complete.obs')^2)
  
  # MEC
  SSE <- sum((y - x) ^ 2, na.rm = T)
  SST <- sum((y - mean(y, na.rm = T)) ^ 2, na.rm = T)
  NSE <- 1 - SSE/SST
  
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
  CCC <- rCb
  
  Cb <- Cb
  r <- r
  
  # return the results
  evalRes <- data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE, CCC = CCC, Cb = Cb)
  
  return(evalRes)
}


# Paramètres globaux

name <- "arg"
kmax <- 23
ntree <- 350
NomsCoord <- c("x", "y")
sample_sizes <- c(600,800,1000,1200,1300,1400,1600,1800,2000,3000,4000,5000,6000,7000,7706) # a modifier pour le pH par incrément de 1000
repets <- 30

#1.Chargement des données et des covariables----

chemin_cov <- "Y:/BDAT/traitement_donnees/MameGadiaga/data/Covariates_MAall"
st <- rast(list.files(chemin_cov, pattern = ".tif$", full.names = TRUE))
dtTB <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_arg.rds")
com<-st_read("Y:/BDAT/traitement_donnees/MameGadiaga/data/commune_53.shp") 

centroides_communes<-st_centroid(com)%>%
  dplyr::select(INSEE_COM, X = X_CENTROID, Y = Y_CENTROID)%>%
  st_drop_geometry()


#2. Extraction des covariables sur les points----

datacov <- terra::extract(st, st_as_sf(dtTB, coords = NomsCoord, crs = 2154)) %>%
  bind_cols(dtTB %>% 
              dplyr::select(all_of(c(name, NomsCoord, "source", "INSEE_COM")))) %>%
  mutate(id = row_number()) %>%
  na.omit()

datacov <- datacov %>%
  mutate(INSEE_COM = as.character(INSEE_COM))

#3. Variables retenues par Boruta----

cov_brt <- readRDS(paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "_cov_brt.rds"))
cat_vars_all <- c("bdforet0", "bdforet1", "bdforet2", "bdforet3", "bdforet30", "bdforet4", "bdforet40", "bdforet5",
                  "bdforet50", "bdforet6", "Geology1", "Geology10", "Geology11", "Geology12", "Geology13", "Geology14",
                  "Geology2", "Geology3", "Geology4", "Geology5", "Geology6", "Geology7", "Geology8", "Geology9",
                  "OCS201611", "OCS201612", "OCS2016211", "OCS2016221", "OCS2016222", "OCS201631", "OCS201632",
                  "OCS201634", "OCS201636", "OCS201641", "OCS201642", "OCS201643", "OCS201644", "OCS201645",
                  "OCS201646", "OCS201651", "typo2", "typo3", "typo4", "typo5")

cat_vars <- intersect(cov_brt, cat_vars_all)
num_vars <- setdiff(cov_brt, cat_vars)

#4. Source de fonctions de factorisation----

source("INLA_SPDE_utils.R")
source("factorisation_RF.R")

#5. Boucle principale sur les tailles et répétitions----

results_all <- list()

for (n in sample_sizes) {
  cat(">> Taille d'échantillon :", n, "\n")
  k <- ifelse(n <= 2000, ceiling(n / 200), ceiling(n / 500))
  
  for (rep in 1:repets) {
    cat("  → Répétition :", rep, "\n")
    set.seed(rep)
    data_sample <- datacov %>% sample_n(n)
    data_sample$id <- 1:nrow(data_sample)
    
    ## 5.a. Folds ponctuels et communaux----
    
    fold <- split(data_sample$id, rep(1:k, length.out = n))
    communes <- unique(data_sample$INSEE_COM)
    commune_folds <- split(communes, sample(rep(1:k, length.out = length(communes))))
    
    #5.b.   Appels Random Forest factorisés----
    
    res_rf_N <- run_rf("Nicolas", data_sample, fold, cov_brt, cat_vars, num_vars, name, k, ntree, kmax)
    
    res_rf_L <- run_rf("Lucille", data_sample, commune_folds, cov_brt, cat_vars, num_vars, name, k, ntree, kmax)
    
    # 5.c. Ajout des prédictions dans data_sample----
    
    # Fusion des deux colonnes de prédiction dans data_sample
    data_sample <- res_rf_L$data_sample  
    data_sample$predRF_N <- res_rf_N$preds  
    
    
    # Restriction de data_sample aux colonnes nécessaires pour INLA
    data_sample <- data_sample %>%
      dplyr::select(all_of(c(name, "INSEE_COM", "id", "predRF_N", "predRF_L",NomsCoord)))
    
    # 5.d. Appels INLA SPDE (KO / KED)----
    res_ko_N <- run_inla_spde(
      type = "KO",aggregation = "Nicolas", data_sample = data_sample,
      centroides_communes = centroides_communes,fold = fold,
      name = name,k = k, NomsCoord = c("x", "y")
    )
    res_ko_L <- run_inla_spde(
      type = "KO", aggregation = "Lucille", data_sample = data_sample,
      centroides_communes = centroides_communes, fold = commune_folds,
      name = name, k = k, NomsCoord = c("x", "y")
    )
    res_ked_N <- run_inla_spde(
      type = "KED", aggregation = "Nicolas", data_sample = data_sample,
      centroides_communes = centroides_communes, fold = fold,
      name = name, k = k, NomsCoord = c("x", "y")
    )
    res_ked_L <- run_inla_spde(
      type = "KED", aggregation = "Lucille", data_sample = data_sample,
      centroides_communes = centroides_communes, fold = commune_folds,
      name = name, k = k, NomsCoord = c("x", "y")
    )
    
    
    
    # 5.e. Résultats----
    
    res_all <- bind_rows(
      res_rf_N$evaluation,
      res_rf_L$evaluation,
      res_ko_N$evaluation,
      res_ko_L$evaluation,
      res_ked_N$evaluation,
      res_ked_L$evaluation
    ) %>% mutate(rep = rep, sample_size = n)
    
    
    results_all[[paste0(n, "_", rep)]] <- res_all
  }
}

# 6. Résumé des résultats----

results_final <- bind_rows(results_all)
results_summary <- results_final %>%
  group_by(sample_size, method) %>%
  summarise(across(c(ME, MAE, RMSE, r, r2, NSE, CCC, Cb), mean, na.rm = TRUE), .groups = "drop") 

# 7. Sauvegarde des résultats----

saveRDS(results_summary, paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/resultats_degradation_", name, ".rds"))

# 8. Graphiques----

metrics <- c("ME", "MAE", "RMSE", "r", "r2", "NSE", "CCC", "Cb")

for (metric in metrics) {
  gg <- ggplot(results_summary, aes(x = sample_size, y = .data[[metric]], color = method)) +
    geom_vline(xintercept = 2000, color = "black", linetype = "solid", size = 1) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(x = "Taille de l'échantillon", y = metric, color = "Méthode",
         title = paste("Évolution de", metric, "selon la taille d'échantillon")) +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
  print(gg)
}
