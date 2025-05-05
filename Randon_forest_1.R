#Chargement des packages nécessaires ----
library(sf)          # manipulation de données spatiales
library(dplyr)       # manipulation de données
library(Boruta)      # sélection des variables importantes
library(randomForest) 
library(ranger)      # random forest rapide
library(ggplot2)     # visualisation
library(terra)       # manipulation raster
library(foreach)     # boucle parallèle (optionnelle ici)
library(tidyr)       # fonctions tidyverse



#Charger les données contenant les points avec pH et covariables ----

datacov <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/matrice_covariables.rds")

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

#Définition des variables ------

kmax= 5 # pour la parallelisation, le nb de coeurs
ntree = 400 # le nbre d'arbre de random forest
k=10 # pour la k fold cross validation


#Nettoyage des données ----
datacov<-datacov %>% 
  select(-c(id_profil,bdatid,insee,codification)) %>% 
  na.omit()

covs <- datacov[, 7:70]
summary(covs)

summary(datacov)

# Échantillonnage à 5000 points pour réduire la taille du jeu de données
set.seed(123)
datacov_s <- datacov %>% sample_n(5000)
datacov_s <- datacov_s %>% 
  st_drop_geometry()

# Ajout d’un identifiant unique pour la validation croisée
datacov_s$id <- 1:nrow(datacov_s)

#Définir la variable cible et les covariables ----

idvar <- "pH"
# Sélectionner les noms des colonnes de covariables (hors pH, geometry et id)
idcovs <- names(covs)

# Créer les objets X (covariables) et Y (réponse)
X <- datacov_s[, idcovs]
Y <- datacov_s[[idvar]]

#Calibration automatique de mtry pour Random Forest ----

taille <- length(idcovs)  # nombre de covariables disponibles

# Calibration du paramètre mtry avec tuneRF
bestmtry <- tuneRF(
  x = X,
  y = Y,
  stepFactor = 1.3,
  mtryStart = round(sqrt(taille)),
  improve = 1e-5,
  ntree = ntree,
  plot = TRUE
)

# Sélection du meilleur mtry basé sur l’erreur OOB minimale
mtry_opt <- bestmtry[which.min(bestmtry[, 2]), 1]

#Sélection des variables importantes avec Boruta ----
set.seed(456)

# Boruta sélectionne les variables pertinentes en se basant sur l’importance de permutation
result_brt <- Boruta(X, Y, mtry = mtry_opt, ntree = ntree)

# Forcer une décision sur les variables incertaines
result_brt_approche <- TentativeRoughFix(result_brt)

# Liste finale des covariables confirmées
cov_brt <- getSelectedAttributes(result_brt_approche)

# Tri des variables par importance médiane
classement_brt <- attStats(result_brt_approche) %>%
  arrange(desc(medianImp))

# Ajouter la variable réponse aux covariables retenues
cov_brt <- c("pH", cov_brt)

# Création d’un sous-jeu de données avec uniquement les colonnes sélectionnées + identifiant
datacov_shrt <- datacov_s[, c(cov_brt, "id")]

#Entraînement du modèle final de Random Forest Quantile ----
# Création de la formule modèle : pH ~ toutes les autres variables
formule.ranger <- as.formula("pH ~ .")

# Entraînement du modèle avec ranger
QRF_Mod.G <- ranger(
  formula = formule.ranger,
  data = datacov_shrt[, -which(names(datacov_shrt) == "id")],
  num.trees = ntree,
  min.node.size = 3,
  quantreg = TRUE,
  max.depth = 15,
  mtry = mtry_opt,
  importance = "permutation",
  scale.permutation.importance = TRUE
)

#Visualisation de l’importance des covariables retenues ----
# Création d’un tableau trié
Imp_QRF <- data.frame(Importance = QRF_Mod.G$variable.importance,
                      Variable = names(QRF_Mod.G$variable.importance)) %>%
  arrange(desc(Importance))

# Graphe ggplot
ggplot(Imp_QRF, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_point() +
  geom_segment(aes(xend = Variable, yend = 0)) +
  coord_flip() +
  ylab("Importance (Permutation)") +
  xlab("Covariables") +
  ggtitle("Importance des covariables sélectionnées par Boruta")

# Validation croisée à 10 plis ----

cat(" Validation croisée 10 plis ----------------\n")

# Définir le nombre de plis
k <- k
set.seed(789)
folds <- cut(seq(1, nrow(datacov_shrt)), breaks = k, labels = FALSE)

# Initialisation
datacov_shrt$predRF <- NA

resuXval <- foreach(i = 1:k, .errorhandling = 'pass') %do% {
  cat("Fold", i, "\n")
  
  # Séparer les données
  test_indices <- which(folds == i, arr.ind = TRUE)
  test_data <- datacov_shrt[test_indices, ]
  train_data <- datacov_shrt[-test_indices, ]
  
  # Entraîner le modèle
  QRF_Mod.G <- ranger(
    formula = formule.ranger,
    data = train_data[, !(names(train_data) %in% c("id", "predRF"))],
    num.trees = ntree,
    min.node.size = 3,
    quantreg = TRUE,
    max.depth = 15,
    mtry = mtry_opt,
    importance = "permutation",
    scale.permutation.importance = FALSE,
    keep.inbag = FALSE
  )
  
  # Prédire sur l'ensemble de test
  preds <- predict(QRF_Mod.G,
                   data = test_data[, cov_brt[-1]],  # retirer pH
                   type = "quantiles",
                   quantiles = 0.5,
                   num.threads = parallel::detectCores())$predictions
  
  # Stocker les prédictions
  datacov_shrt$predRF[test_indices] <- round(preds, 2)
  
  # Retourner quelque chose pour chaque pli si nécessaire
  list(fold = i, indices = test_indices, predictions = preds)
}

# Calculer les statistiques de validation croisée
resuXvalQRF <-  Myeval(datacov_shrt$predRF,   datacov_shrt$pH )

resuXvalQRF
