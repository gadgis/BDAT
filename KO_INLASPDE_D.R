library(INLA)
library(inlabru)
library(sp)
library(sf)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(purrr)
library(terra)

# Fonction d'évaluation----
Myeval <- function(x, y){
  ME <- round(mean(y - x, na.rm = TRUE), 2)
  RMSE <- round(sqrt(mean((y - x)^2, na.rm = TRUE)), 2)
  MAE <- round(mean(abs(y - x), na.rm = TRUE), 2)
  r2 <- round((cor(x, y, use = 'pairwise.complete.obs')^2), 2)
  SSE <- sum((y - x)^2, na.rm = TRUE)
  SST <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  NSE <- round(1 - SSE/SST, 2)
  data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r2 = r2, NSE = NSE)
}

# Paramètres----
name <- "pH"
sample_sizes <- c(600, 800)
repets <- 5
NomsCoord <- c("x", "y")

# Données----
dtTB<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")
datacov <- dtTB %>%
  dplyr::mutate(identifiant = ifelse(source == "IGCS", id_profil, bdatid)) %>%
  dplyr::select(all_of(c(name, NomsCoord, "INSEE_COM", "source", "identifiant"))) %>%
  na.omit()

# Définir le champ spde----
coords <- datacov[, NomsCoord]

mesh3 <- inla.mesh.2d(
  loc = coords,
  max.edge = c(1, 2) * diff(range(coords[,1])) / (3*5),
  offset = c(diff(range(coords[,1])) / (3*5), diff(range(coords[,1])) / 3),
  cutoff = diff(range(coords[,1])) / (3*5*10)
)

matern <- inla.spde2.pcmatern(
  mesh = mesh3,
  alpha = 2,
  prior.range = c(500, 0.01),
  prior.sigma = c(10, 0.01)
)

#Boucle sur les tailles d'échantillon----
registerDoParallel(cores = parallel::detectCores() - 1)
results_all <- list()

for (n in sample_sizes) {
  cat("\n>>> Taille d’échantillon :", n, "\n")
  indiv_per_fold <- ifelse(n <= 2000, 200, 500)
  k <- ceiling(n / indiv_per_fold)
  
  resu_rep <- foreach(rep = 1:repets, .packages = c("INLA", "inlabru", "sp", "dplyr", "sf"),
                      .export = c("Myeval","matern")) %do% {
    set.seed(rep)
    
                        ##Echantillonage----
    data_sample <- datacov %>% sample_n(n)
    data_sample$id <- 1:nrow(data_sample)
    
    ##Définit les folds----
    
    folds <- split(data_sample$id, rep(1:k, length.out = n))
    
    data_sample$elt <- data_sample[[name]]
    data_sample$predINLAKO <- NA
    
    coordinates(data_sample) <- NomsCoord
    proj4string(data_sample) <- CRS("epsg:2154")
    
    ##Validation croisée----
    
    data_sample$predINLAKO <- NA
    
    for (i in 1:k) {
      data_sample$elt[folds[[i]]] <- NA
      
      cmp <- elt ~ Intercept(1) + field(coordinates, model = matern)
      
      fit <- bru(
        cmp,
        data = data_sample,
        family = "Gaussian",
        domain = list(coordinates = mesh3),
        options = list(control.inla = list(int.strategy = "eb"), verbose = FALSE)
      )
      
      preds <- fit$summary.fitted.values$mean[1:nrow(data_sample)]
      mask <- is.na(data_sample$elt)
      data_sample$predINLAKO[mask] <- preds[mask]
    }
    
    Myeval(data_sample$predINLAKO, data_sample[[name]]) %>%
      mutate(rep = rep, sample_size = n)
    
  }
  
  results_all[[as.character(n)]] <- resu_rep
}

#Fin cluster paral----
stopImplicitCluster()

results_metrics <- bind_rows(results_all)

results_summary <- results_metrics %>%
  group_by(sample_size) %>%
  summarise(across(c(ME, MAE, RMSE, r2, NSE), mean, na.rm = TRUE))

saveRDS(results_summary, "resultats_INLAKO_Degrad.rds")
