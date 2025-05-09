#Chargement des packages----
library(sf)
library(randomForest)
library(ranger)
library(Boruta)
library(terra)
library(raster)
library(writexl)
library(ggplot2)
library(foreach)
library(dplyr)
library(tidyverse)
library(tidyr)

#Définition des fonctions ----
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
  
  # return the results
  evalRes <- data.frame(ME = ME, MAE = MAE, RMSE = RMSE, r = r, r2 = r2, NSE = NSE)
  
  return(evalRes)
}

kmax<-5
ntree=500

# Chargement des données----
datacov <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/centroides_covariables.rds")

datacov_s <- datacov %>%
  select(-c(INSEE_COM,X,Y)) %>%
  na.omit()

# définir les colonnes de covariables
covs <- datacov_s[, 2:65]  
idcovs <- names(covs)

datacov_s$id <- 1:nrow(datacov_s)

# Définir la variable réponse
idvar <- "moy_pH"
Y <- datacov_s[[idvar]]
X <- datacov_s[, idcovs]

summary(datacov_s)
  
  # Étape 1 : Tune mtry
  bestmtry <- tuneRF(
    x = X,
    y = Y,
    stepFactor = 1.3,
    mtryStart = round(sqrt(length(idcovs))),
    improve = 1e-5,
    ntree = ntree,
    plot = TRUE
  )
  mtry_opt <- bestmtry[which.min(bestmtry[, 2]), 1]
  
  # Étape 2 : Boruta
  set.seed(456)
  result_brt <- Boruta(X, Y, mtry = mtry_opt, ntree = ntree)
  result_brt_approche <- TentativeRoughFix(result_brt)
  cov_brt <- getSelectedAttributes(result_brt_approche)
  cov_brt <- c("moy_pH", cov_brt)
  
  # Étape 3 : Création du sous-jeu Boruta
  datacov_shrt <- datacov_s[, c(cov_brt, "id")]
  formule.ranger <- as.formula("moy_pH ~ .")
  
  # Étape 4 : Visualisation de l'importance des variables
  QRF_Mod.G <- ranger(
    formula = formule.ranger,
    data = datacov_shrt[, !(names(datacov_shrt) %in% c("id"))],
    num.trees = ntree,
    min.node.size = 3,
    quantreg = TRUE,
    max.depth = 15,
    mtry = mtry_opt,
    importance = "permutation",
    scale.permutation.importance = TRUE
  )
  
  Imp_QRF <- data.frame(Importance = QRF_Mod.G$variable.importance,
                        Variable = names(QRF_Mod.G$variable.importance)) %>%
    arrange(desc(Importance))
  
  plot_path <- paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/varimp_centroides_", ntree, ".png")
  png(plot_path, width = 800, height = 600)
  print(
    ggplot(Imp_QRF, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_point() +
      geom_segment(aes(xend = Variable, yend = 0)) +
      coord_flip() +
      ylab("Importance (Permutation)") +
      xlab("Covariables") +
      ggtitle(paste("Importance des covariables - ntree =", ntree))
  )
  dev.off()
  
  # Étape 5 : Validation croisée 10 plis
  k <- 10
  set.seed(789)
  folds <- cut(seq(1, nrow(datacov_shrt)), breaks = k, labels = FALSE)
  datacov_shrt$predRF <- NA
  
  for (i in 1:k) {
    cat(" Fold", i, "\n")
    test_idx <- which(folds == i)
    test_data <- datacov_shrt[test_idx, ]
    train_data <- datacov_shrt[-test_idx, ]
    
    QRF_Mod.G <- ranger(
      formula = formule.ranger,
      data = train_data[, !(names(train_data) %in% c("id", "predRF"))],
      num.trees = ntree,
      min.node.size = 3,
      quantreg = TRUE,
      max.depth = 15,
      mtry = mtry_opt,
      importance = "permutation",
      scale.permutation.importance = FALSE
    )
    
    
    preds <- predict(QRF_Mod.G,
                     data=test_data[ , cov_brt[-1]],
                     type = "quantiles",
                     quantiles = 0.5,
                     num.threads = kmax)$predictions
    
    datacov_shrt$predRF[test_idx] <- round(preds, 2)
  }
  
  # Étape 6 : Évaluation
  resuXvalQRF <- Myeval(datacov_shrt$predRF, datacov_shrt$pH)
  resuXvalQRF$ntree <- ntree
  
  
  # Étape 7 : Export Excel
  nom_fichier <- paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrics_RF_centroides_", ntree, ".xlsx")
  write_xlsx(resuXvalQRF, path = nom_fichier)


