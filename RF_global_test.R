#Chargement des packages----
library(sf)
library(randomForest)
library(ranger)
library(Boruta)
library(raster)
library(writexl)
# library(ggplot2) --> pour info, le tidyverse importe ggplot2, dplyr et tidyr
library(foreach)
# library(dplyr)
library(tidyverse)
# library(tidyr)

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

kCV <- function(f){
  require(ranger)
  QRF_Mod.G <- ranger(
    formula = as.formula("pH ~ ."),
    data = X2[-f, -1],
    num.trees = nt,
    min.node.size = 3,
    quantreg = TRUE,
    max.depth = 15,
    mtry = mtry_opt,
    importance = "permutation",
    scale.permutation.importance = FALSE,
    keep.inbag=FALSE
  )
  
  preds <- predict(QRF_Mod.G,
                   data= X2[f, -c(1,2)],  # retirer pH
                   type = "quantiles",
                   quantiles = 0.5,
                   num.threads = parallel::detectCores())$predictions
  
  preds
}

# Définition des paramètres à tester-----
liste_ntree <- seq(50,400,50)
tous_les_resultats <- list()  # stockera les résultats pour chaque ntree

# Chargement des données----
# changer le répertoire de travail
setwd("J:/traitement_donnees/MameGadiaga/resultats")

datacov <- readRDS("matrice_covariables.rds")

datacov_s <- datacov %>%
  select(-c(id_profil, bdatid, insee, codification)) %>%
  na.omit()

# définir les colonnes de covariables
covs <- datacov_s[, 7:70]  
idcovs <- names(covs)

datacov_s$id <- 1:nrow(datacov_s)

# Définir la variable réponse
idvar <- "pH"
Y <- datacov_s[[idvar]]
X <- datacov_s[, idcovs]

# Sélection d'un sous ensemble de variables -------------
# Boruta sélectionne les variables pertinentes en se basant sur l’importance de permutation
result_brt <- Boruta(X, Y)

plot(result_brt)

# Forcer une décision sur les variables incertaines
result_brt_approche <- TentativeRoughFix(result_brt)

# Liste finale des covariables confirmées
cov_brt <- getSelectedAttributes(result_brt_approche)

#enregistrer la liste des covariables
save(cov_brt, file = "liste_variables_boruta.Rds")
load("liste_variables_boruta.Rds")

# créer un sous-jeu de données avec les variables retenues par Boruta
X2 <-  datacov_s[, c("id", "pH", cov_brt)]


# Calibrage des hyperparamètres -----------------

## mettre en place la parallélisation (sera utilisée pour la kCV) ----------
#nombre de coeurs utilisés
n.cores <- parallel::detectCores() - 1

#créér le cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# définir les clusters comme hôtes de la fonction %dopar%
doParallel::registerDoParallel(cl = my.cluster)


## Calcul ------------
res_calib <- lapply(liste_ntree,
                    
                    function(nt){
                      # Étape 1 : Tune mtry
                      bestmtry <- tuneRF(
                        x = X2[,-c(1,2)],
                        y = Y,
                        stepFactor = 1.3,
                        mtryStart = round(sqrt(length(idcovs))),
                        improve = 1e-5,
                        ntree = nt,
                        plot = FALSE
                      )
                      mtry_opt <- bestmtry[which.min(bestmtry[, 2]), 1]
                      
                      # Etape 2 : Validation croisée 
                      folds <- pls::cvsegments(N = length(Y),k = 10)
                      
                      res <- foreach(f = folds) %dopar% {
                        require(ranger)
                        QRF_Mod.G <- ranger(
                          formula = as.formula("pH ~ ."),
                          data = X2[-f, -1],
                          num.trees = nt,
                          min.node.size = 3,
                          quantreg = TRUE,
                          max.depth = 15,
                          mtry = mtry_opt,
                          importance = "permutation",
                          scale.permutation.importance = FALSE,
                          keep.inbag=FALSE
                        )
                        
                        preds <- predict(QRF_Mod.G,
                                         data= X2[f, -c(1,2)],  # retirer pH
                                         type = "quantiles",
                                         quantiles = 0.5,
                                         num.threads = parallel::detectCores())$predictions
                        
                        preds
                      }
                      
                      # Etape 3 : calculer les métriques
                      
                      # recoller les résultats
                      res_df <- data.frame(obs = Y, pred = NA) 
                      for(k in 1:10){
                        res_df$pred[folds[[k]]] <- res[[k]]
                      }
                      
                      #calculer les métriques
                      metrique <- Myeval(res_df$obs, res_df$pred)
                      metrique$ntree <- nt
                      metrique$mtry <- mtry_opt
                      
                      # renvoyer les métriques
                      return(metrique)
                    })

# SUITE SCRIPT MAME -----------------------------

for (ntree_val in liste_ntree) {
  cat("\n Traitement ntree =", ntree_val, "----------------------\n")
  
  # Étape 1 : Tune mtry
  bestmtry <- tuneRF(
    x = X,
    y = Y,
    stepFactor = 1.3,
    mtryStart = round(sqrt(length(idcovs))),
    improve = 1e-5,
    ntree = ntree_val,
    plot = FALSE
  )
  mtry_opt <- bestmtry[which.min(bestmtry[, 2]), 1]
  
  # Étape 2 : Boruta
  set.seed(456)
  result_brt <- Boruta(X, Y, mtry = mtry_opt, ntree = ntree_val)
  result_brt_approche <- TentativeRoughFix(result_brt)
  cov_brt <- getSelectedAttributes(result_brt_approche)
  cov_brt <- c("pH", cov_brt)
  
  # Étape 3 : Création du sous-jeu Boruta
  datacov_shrt <- datacov_s[, c(cov_brt, "id")]
  formule.ranger <- as.formula("pH ~ .")
  
  # Étape 4 : Visualisation de l'importance des variables
  QRF_Mod.G <- ranger(
    formula = formule.ranger,
    data = datacov_shrt[, !(names(datacov_shrt) %in% c("id"))],
    num.trees = ntree_val,
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
  
  plot_path <- paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/varimp_ntree_", ntree_val, ".png")
  png(plot_path, width = 800, height = 600)
  print(
    ggplot(Imp_QRF, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_point() +
      geom_segment(aes(xend = Variable, yend = 0)) +
      coord_flip() +
      ylab("Importance (Permutation)") +
      xlab("Covariables") +
      ggtitle(paste("Importance des covariables - ntree =", ntree_val))
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
      num.trees = ntree_val,
      min.node.size = 3,
      quantreg = TRUE,
      max.depth = 15,
      mtry = mtry_opt,
      importance = "permutation",
      scale.permutation.importance = FALSE,
      keep.inbag=FALSE
    )
    
  
    preds <- predict(QRF_Mod.G,
                    data= test_data[, cov_brt[-1]],  # retirer pH
                     type = "quantiles",
                     quantiles = 0.5,
                     num.threads = parallel::detectCores())$predictions
    
    datacov_shrt$predRF[test_idx] <- round(preds, 2)
  }
  
  # Étape 6 : Évaluation
  resuXvalQRF <- Myeval(datacov_shrt$predRF, datacov_shrt$pH)
  resuXvalQRF$ntree <- ntree_val
  tous_les_resultats[[paste0("ntree_", ntree_val)]] <- resuXvalQRF
  
  # Étape 7 : Export Excel
  nom_fichier <- paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/metrics_RF_ntree_", ntree_val, ".xlsx")
  write_xlsx(resuXvalQRF, path = nom_fichier)
  cat("Export metrics dans :", nom_fichier, "\n")
}

# Résumé final de tous les résultats
df_resultats_final <- dplyr::bind_rows(tous_les_resultats)

# Export du tableau récapitulatif complet
write_xlsx(df_resultats_final, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/tableau_recapitulatif_RF.xlsx")

# Graphe des métriques selon ntree

df_long <- df_resultats_final %>%
  pivot_longer(cols = c(ME, MAE, RMSE, r2, NSE), names_to = "metrique", values_to = "valeur")

ggplot(df_long, aes(x = factor(ntree), y = valeur, fill = metrique)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Nombre d'arbres (ntree)", y = "Valeur de la métrique", title = "Performance du modèle RF selon ntree") +
  theme_minimal()

# Optionnel : sauvegarder le graphe
ggsave("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/performance_RF_par_ntree.png", width = 10, height = 6)
