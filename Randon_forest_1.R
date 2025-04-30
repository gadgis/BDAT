#Chargement des packages nécessaires ----
library(dplyr)       # manipulation de données
library(Boruta)      # sélection des variables importantes
library(ranger)      # random forest rapide
library(ggplot2)     # visualisation
library(terra)       # manipulation raster
library(foreach)     # boucle parallèle (optionnelle ici)
library(tidyr)       # fonctions tidyverse



#Charger les données contenant les points avec pH et covariables ----

datacov <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/matrice_covariables.rds")


#Nettoyage des données ----
datacov<-datacov %>% 
  select(-c(id_profil,bdatid,insee,codification)) %>% 
  na.omit()

covs <- datacov[, 6:70]
summary(covs)

summary(datacov)

# Échantillonnage à 5000 points pour réduire la taille du jeu de données
set.seed(123)
datacov_s <- datacov %>% sample_n(5000)

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

ntree <- 400  # nombre d’arbres à utiliser
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
set.seed(789)
k <- 10  # nombre de plis
folds <- cut(seq(1, nrow(datacov_shrt)), breaks = k, labels = FALSE)

# Initialiser la colonne de prédictions
datacov_shrt$predRF <- NA

# Boucle CV
resuXval <- foreach(i = 1:k, .combine = rbind) %do% {
  cat(" Pli CV :", i, "\n")
  
  test_idx <- which(folds == i)
  train_data <- datacov_shrt[-test_idx, ]
  test_data <- datacov_shrt[test_idx, ]
  
  rf_cv <- ranger(
    formula = formule.ranger,
    data = train_data[, !(names(train_data) %in% c("id", "predRF"))],
    num.trees = ntree,
    quantreg = TRUE,
    mtry = mtry_opt,
    max.depth = 15,
    min.node.size = 3
  )
  
  preds <- predict(rf_cv, test_data[, cov_brt[-1]], type = "quantiles", quantiles = 0.5)$predictions
  data.frame(id = test_data$id, predRF = preds)
}

# Associer les prédictions aux bonnes lignes
datacov_shrt <- left_join(datacov_shrt, resuXval, by = "id")

#Évaluation des performances du modèle ----
# Fonctions d’évaluation
RMSE <- function(obs, pred) sqrt(mean((obs - pred)^2))
R2 <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
MAE <- function(obs, pred) mean(abs(obs - pred))

cat("résultats de la validation croisée (10 plis) :\n")
cat("- RMSE :", RMSE(datacov_shrt$pH, datacov_shrt$predRF), "\n")
cat("- MAE  :", MAE(datacov_shrt$pH, datacov_shrt$predRF), "\n")
cat("- R²   :", R2(datacov_shrt$pH, datacov_shrt$predRF), "\n")

#Spatialisation finale sur la grille gXY ----
# gXY est un tableau de données rasterisé avec les mêmes covariables (voir script 1)
gXY <- readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/gXY_covariables.rds")

# Préparer les données en sélectionnant les mêmes covariables
testD <- gXY %>% dplyr::select(all_of(cov_brt[-1]))  # on enlève pH

# Filtrer uniquement les lignes complètes
pred_input <- gXY %>%
  filter_at(vars(all_of(cov_brt[-1])), all_vars(!is.na(.))) %>%
  dplyr::select(x, y)

# Prédiction avec le modèle RF
QRF_Median <- predict(QRF_Mod.G,
                      testD,
                      type = "quantiles",
                      quantiles = c(0.5),
                      num.threads = parallel::detectCores())$predictions

# Création d’un tableau spatial avec les prédictions
QRF_Median50 <- cbind(pred_input, QRF_Median = QRF_Median[, 1])

# Création raster
r <- rast(QRF_Median50, type = "xyz")

# Sauvegarde du raster final
writeRaster(r, filename = "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/prediction_pH_RF.tif", overwrite = TRUE)
