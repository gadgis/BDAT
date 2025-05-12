#Chargement des packages----
library(sf)
library(randomForest)
library(ranger)
library(Boruta)
library(raster)
library(writexl)
library(foreach)
library(tidyverse)


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


# Définition des paramètres à tester-----
liste_ntree <- seq(50,500,50)

# Chargement des données----
# changer le répertoire de travail
setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

data<- readRDS("matrice_covariables_C.rds")

datacov_s <- data %>%
  select(-c(id_profil, bdatid, insee, codification)) %>%
  na.omit()

# définir les colonnes de covariables
covs <- datacov_s[, 7:70]  
idcovs <- names(covs)

datacov_s$id <- 1:nrow(datacov_s)

# Définir la variable réponse
idvar <- "C"
Y <- datacov_s[[idvar]]
X <- datacov_s[, idcovs]

# Sélection d'un sous ensemble de variables -------------
#Boruta sélectionne les variables pertinentes en se basant sur l’importance de permutation
result_brt <- Boruta(X, Y)

plot(result_brt)

# Forcer une décision sur les variables incertaines
result_brt_approche <- TentativeRoughFix(result_brt)

# Liste finale des covariables confirmées
cov_brt_C <- getSelectedAttributes(result_brt_approche)
# 
# #enregistrer la liste des covariables
save(cov_brt_C, file = "cov_brt_C.Rds")
load("cov_brt_C.Rds")

# créer un sous-jeu de données avec les variables retenues par Boruta
X2 <-  datacov_s[, c("id", "C", cov_brt_C)]


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

taille  <- ncol(X2) - 2                    # nb de variables explicatives
folds   <- pls::cvsegments(N = length(Y), k = 10)

# Boucle parallèle sur ntree -

res_calib <- foreach(
  nt        = liste_ntree,               # valeur de ntree (une par tâche)
  .combine  = rbind, 
  .packages = c("randomForest", "ranger", "pls"),  
  .export   = c("X2", "Y", "taille", "folds", "Myeval")
) %dopar% {
  
  ## -- Étape 1 : Tune mtry -----
  bestmtry <- tuneRF(
    x          = X2[, -c(1, 2)],
    y          = Y,
    mtryStart  = round(sqrt(taille)),
    stepFactor = 1.3,
    improve    = 1e-5,
    ntreeTry   = nt,            # <<< argument correct (et non ntree)
    trace      = FALSE,
    plot       = FALSE
  )
  mtry_opt <- bestmtry[which.min(bestmtry[, "OOBError"]), "mtry"]
  
  ## -- Étape 2 : 10-fold CV ----------
  preds <- numeric(length(Y))              # vecteur final de prédictions
  
  for (seg in folds) {
    mod <- ranger(
      formula      = C ~ .,
      data         = X2[-seg, -1],         # -1 supprime la colonne id
      num.trees    = nt,
      mtry         = mtry_opt,
      min.node.size = 3,
      quantreg     = TRUE,
      max.depth    = 15,
      importance   = "permutation",
      num.threads  = 1,                    # 1 thread par worker
      keep.inbag   = FALSE
    )
    
    preds[seg] <- predict(
      mod,
      data       = X2[seg, -c(1, 2)],      # set de test (sans C et id)
      type       = "quantiles",
      quantiles  = 0.5
    )$predictions
  }
  
  ## -- Étape 3 : métriques 
  metrique         <- Myeval(Y, preds)
  metrique$ntree   <- nt
  metrique$mtry    <- mtry_opt
  metrique                              # renvoyé à .combine = rbind
}

# 3. Libération des ressources
parallel::stopCluster(my.cluster)

saveRDS(res_calib, file = "resultats_calib_C.Rds")

# res1<- readRDS("resultats_calib.Rds")
# res2<- readRDS("resultats_calib2.Rds")
# 
# res_calib <- rbind(res1, res2)

# # Graphe des métriques selon ntree
# 
# res_calib <- res_calib %>%
#   pivot_longer(cols = c(ME, MAE, RMSE, r2, NSE), names_to = "metrique", values_to = "valeur")
# 
# ggplot(res_calib, aes(x = factor(ntree), y = valeur, fill = metrique)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(x = "Nombre d'arbres (ntree)", y = "Valeur de la métrique", title = "Performance du modèle RF selon ntree") +
#   theme_minimal()

 
  # Étape 4 : Visualisation de l'importance des 


#   QRF_Mod.G <- ranger(
#     formula = C ~ .,
#     data = X2[, !(names(X2) %in% c("id"))],
#     num.trees = ,
#     min.node.size = 3,
#     quantreg = TRUE,
#     max.depth = 15,
#     mtry = mtry_opt,
#     importance = "permutation",
#     scale.permutation.importance = TRUE
#   )
#   
#   Imp_QRF <- data.frame(Importance = QRF_Mod.G$variable.importance,
#                         Variable = names(QRF_Mod.G$variable.importance)) %>%
#     arrange(desc(Importance))
#   
#   plot_path <- paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/varimp_ntree_", ntree_val, ".png")
#   png(plot_path, width = 800, height = 600)
#   print(
#     ggplot(Imp_QRF, aes(x = reorder(Variable, Importance), y = Importance)) +
#       geom_point() +
#       geom_segment(aes(xend = Variable, yend = 0)) +
#       coord_flip() +
#       ylab("Importance (Permutation)") +
#       xlab("Covariables") +
#       ggtitle(paste("Importance des covariables - ntree =", ntree_val))
#   )
#   dev.off()
#   
#   
# 
# 
# 
# # Optionnel : sauvegarder le graphe
# ggsave("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/performance_RF_par_ntree.png", width = 10, height = 6)
