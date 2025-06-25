# RandomForest ----

##### BORUTA #####
### A CHANGER SI CHANGE NB DE COVARIABLES

X = datacov[,idcovs]                 # X: variables indépendantes -> covariables
Y = datacov[,idvar]                  # Y: Variable cible (target variable, outcome) -> ETM

taille = ncol(datacov)-3
#Tune du mtry
# bestmtry = tuneRF(X,Y,                            #tune du RF pour déterminer le meilleur mtry
#                   stepFactor = 1.3,
#                   mtry = round(sqrt(taille)),
#                   improve = 1e-5,
#                   ntree = ntree,
#                   plot = FALSE)           #utilisation du même ntree que soil 2.0 pour leur rfe
# 



library(iml)
library(mlr)




# #Boruta
# result_brt = Boruta(X, Y,                         #classification de l'importance des covariables par boruta
#                     mtry = taille ,
#                     min.node.size = 3 ,
#                     ntree = ntree)
# 
# 
# Stats_brt = attStats(result_brt)                  #résultats du boruta
# cov_brt = getSelectedAttributes(result_brt)       #sélection des covariables confirmées
# 
# 
# result_brt_approche = TentativeRoughFix(result_brt)#résultats en forçant les indécis
# cov_brt = getSelectedAttributes((result_brt_approche))      #sélection du nouvel ensemble
# 
# classement_brt = Stats_brt %>% 
#   arrange(desc(medianImp))                        #classe les covariables selon la valeur de la médiane
# 
# classement_brt_approche = attStats(result_brt_approche) %>% 
#   arrange(desc(medianImp))                        #pareil pour les covariables plus complètes
# 
# classement_brt 


cov_brt<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_cov_brt.rds")

#Création du tableau des covariables séléctionnées

cov_brt = c(name, cov_brt)
datacov_shrt = datacov[cov_brt]                   #récupération dans datacov des colonnes conservées

SOC.task = makeRegrTask(data = datacov_shrt, target = name)

# Rough Estimation of the Tuning time
estimateTimeTuneRanger(SOC.task,num.threads = 60,
                       num.trees = ntree)

# Tuning process (takes around 1 minute); Tuning measure is the multiclass brier score
res = tuneRanger(SOC.task,
                 num.trees = ntree,
                 iters = 100,
                 num.threads = 60)



res$recommended.pars$mtry
res$recommended.pars$min.node.size



# fomula.ranger <- as.formula(paste0(name,"~."))
# QRF_Mod.G <- ranger(formula = fomula.ranger ,
#                     data = datacov_shrt,
#                     num.trees = ntree,
#                     min.node.size = res$recommended.pars$mtry ,
#                     quantreg = TRUE,
#                     max.depth = 15, 
#                     mtry=res$recommended.pars$mtry ,
#                     importance="permutation", 
#                     scale.permutation.importance = TRUE, #division par l'écart-type de la variable (mise des permutations entre 0 et 1)
#                     keep.inbag = F)

fomula.ranger <- as.formula(paste0(name,"~."))
RF_Mod.G <- ranger(formula = fomula.ranger ,
                    data = datacov_shrt,
                    num.trees = ntree,
                    min.node.size = res$recommended.pars$mtry ,
                    quantreg = F,
                    max.depth = 15, 
                    mtry=res$recommended.pars$mtry ,
                    importance="permutation", 
                    scale.permutation.importance = TRUE, #division par l'écart-type de la variable (mise des permutations entre 0 et 1)
                    keep.inbag = F)
#Variable importance
Imp_RF <- data.frame(RF_Mod.G$variable.importance)
Imp_RF$Vars <- row.names(Imp_RF)
Imp_RF <- Imp_RF[order(Imp_RF$RF_Mod.G.variable.importance,decreasing = T),]


Imp_RF <- Imp_RF %>%
  left_join(df_vars, by = "Vars")

Imp_RF <- Imp_RF %>%
  rename(importance = RF_Mod.G.variable.importance)

Imp_RF$label <- iconv(Imp_RF$label, from = "", to = "UTF-8")

Imp_RF_top10 <- Imp_RF %>%
  arrange(desc(importance)) %>%
  slice(1:10)

varimp <- ggplot(Imp_RF_top10, 
                 aes(x = reorder(label, importance), 
                     y = importance,
                     color = var_group)
) + 
  geom_point(size = 1.5) +
  geom_segment(aes(x = label, xend = label, y = 0, yend = importance)) +
  scale_color_manual(
    name = "Groupe de variables",
    values = c(
      "Topographie" = "grey50",
      "Climat" = "#1f77b4",
      "Vegetation" = "#2ca02c",
      "Sol" = "#8c564b",
      "Lithologie"= "#d62728",
      "Radiometrie" = "#9467bd",
      "Occupation_Sol" = "#ffdf00" )) +
  ylab("Importance (IncNodePurity)") +
  xlab("Variables") +
  ggtitle(paste("Dix premières variables importantes -", "pH")) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

print(varimp)

# Predictions

#spatial prediction


# prepare for predction the new data table
testD <- gXY %>% 
  dplyr::select(all_of(cov_brt[-1])) 



# QRF_Median <- predict(QRF_Mod.G,
#                       testD,
#                       type = "quantiles",
#                       quantiles =  c(0.05,0.5,0.95),
#                       num.threads = kmax )$predictions

QRF_Median2 <- predict(RF_Mod.G,
                      testD,
                      num.threads = kmax )$predictions

QRF_Median50 <- bind_cols(gXY %>%
                            filter_at(vars(cov_brt[-1]),
                                      all_vars(!is.na(.))
                                      ) %>%
                            dplyr::select(x,y)     ,
                          QRF_Median = QRF_Median2)

r <- rast(QRF_Median50, type="xyz")


terra::writeRaster(r,file=paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"qrf.tif"),overwrite=TRUE)

# k fold -------------

print("Validation croisée----------------")

datacov$predRF <- NA


resuXval <-
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {

            print(i)

            # collecter les # des lignes (gérer les doublons)
            nblignes = which( datacov$id %in% datacov$id[ fold[[i]] ] )

            RF_Mod.G <- ranger(formula = fomula.ranger ,
                                data = datacov_shrt[-nblignes  , ],
                                num.trees = ntree,
                                mtry=res$recommended.pars$mtry ,
                                min.node.size = res$recommended.pars$mtry ,

                                quantreg = FALSE,
                                max.depth = 15,
                                importance="permutation",
                                scale.permutation.importance = FALSE, #division par l'écart-type de la variable (mise des permutations entre 0 et 1)
                                keep.inbag = F)


            datacov$predRF[ nblignes ] <- predict(RF_Mod.G,
                                                    datacov_shrt[ nblignes , ],
                                                    num.threads = kmax )$prediction



          }


resuXvalQRF <-  Myeval(datacov$predRF,   datacov[,name] )

#Validation croisée  Méthode Nicolas----

# Validation croisee Random Forest

registerDoParallel(cores = kmax) 
set.seed(123)
seg <- split(1:nrow(datacov), sample(rep(1:k, length.out = nrow(datacov))))

# Initialisation
datacov$predRF_aggr <- NA
datacov$predRF_aggr <- as.numeric(datacov$predRF_aggr)

# Séparer les covariables numériques et catégorielles
cat_vars <- c("Geology1", "Geology10", "Geology11", "Geology12", "Geology2", "Geology3",
              "Geology5", "Geology8", "Geology9", "OCS201611", "OCS2016211",
              "OCS2016221", "typo3", "typo4", "typo5")
num_vars <- setdiff(cov_brt[-1], cat_vars)

resuXval_aggr <- foreach(i = 1:k,  .combine = rbind, .packages = c("dplyr", "ranger", "mlr", "tuneRanger")) %do% {
  print(i)
  cat("seg:", i, "\n")
  
  nblignes <- seg[[i]] 
  train_data <- datacov[-nblignes, ]
  test_data <- datacov[nblignes, ]
  
  # Agrégation sur les données d'entraînement - covariables
  train_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(
      across(all_of(num_vars), ~mean(.x, na.rm = TRUE)),
      across(all_of(cat_vars), ~as.character(names(which.max(table(.)))))
    ) %>%
    ungroup()
  
  # Agrégation de la variable cible - variable réponse
  y_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(value = mean(.data[[name]], na.rm = TRUE)) %>%
    pull(value)
  train_aggr[[name]] <- y_aggr
  
  # Conversion des variables catégorielles
  train_aggr[cat_vars] <- lapply(train_aggr[cat_vars], factor)
  train_aggr[num_vars] <- lapply(train_aggr[num_vars], as.numeric)
  
  # Tâche mlr - définir le modèle
  rf_task <- makeRegrTask(data = train_aggr[, c(name, cov_brt[-1])], target = name)
  res_tune <- tuneRanger(rf_task, num.trees = ntree, iters = 100, num.threads = 4)
  
  # Entraînement sur les données agrégées
  rf_model <- ranger(
    formula = as.formula(paste0(name, " ~ .")),
    data = train_aggr[, c(name, cov_brt[-1])],
    num.trees = ntree,
    mtry = res_tune$recommended.pars$mtry,
    min.node.size = res_tune$recommended.pars$min.node.size,
    importance = "permutation",
    keep.inbag = FALSE
  )
  

  
  # Données de test : covariables brutes
  test <- test_data[, cov_brt[-1]]
  test[cat_vars] <- lapply(test[cat_vars], factor)
  test[num_vars] <- lapply(test[num_vars], as.numeric)
  preds <- predict(rf_model, data = test, num.threads = kmax)$predictions
  preds <- round(preds, 2)
  datacov$predRF_aggr[nblignes] <- preds
  
  # Renvoi des résultats pour le tableau final
  return(data.frame(fold = i, pred = preds, obs = test_data[[name]]))
}

resuXvalRF_commune <- Myeval(datacov$predRF_aggr, datacov[[name]])

#Méthode de Lucille----
# Créer les folds par commune
set.seed(123)
communes <- unique(datacov$INSEE_COM)
folds_com <- split(communes, sample(rep(1:k, length.out = length(communes))))

# Initialisation de la colonne des prédictions
datacov$predRF_aggrcom <- NA_real_

# Boucle séquentielle (pas dopar ici)
resuXval_aggrcom <- foreach(i = 1:k, .combine = rbind, .packages = c("dplyr", "ranger", "mlr", "tuneRanger")) %do% {
  cat("Fold commune:", i, "\n")
  
  test_communes <- folds_com[[i]]
  train_communes <- setdiff(communes, test_communes)
  
  # Séparation du jeu d'entraînement et de test
  train_data <- datacov %>% filter(INSEE_COM %in% train_communes)
  test_data  <- datacov %>% filter(INSEE_COM %in% test_communes)
  
  # Agrégation du jeu d'entraînement par commune
  train_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(
      across(all_of(num_vars), ~mean(.x, na.rm = TRUE)),
      across(all_of(cat_vars), ~as.character(names(which.max(table(.))))),
      .groups = "drop"
    )
  
  # Ajouter la variable cible agrégée
  y_aggr <- train_data %>%
    group_by(INSEE_COM) %>%
    summarise(value = mean(.data[[name]], na.rm = TRUE)) %>%
    pull(value)
  train_aggr[[name]] <- y_aggr
  
  # Conversion des types
  train_aggr[cat_vars] <- lapply(train_aggr[cat_vars], factor)
  train_aggr[num_vars] <- lapply(train_aggr[num_vars], as.numeric)
  
  # Création de la tâche mlr
  rf_task <- makeRegrTask(data = train_aggr[, c(name, cov_brt[-1])], target = name)
  res_tune <- tuneRanger(rf_task, num.trees = ntree, iters = 100, num.threads = 4)
  
  # Entraînement du modèle sur les communes agrégées
  rf_model <- ranger(
    formula = as.formula(paste0(name, " ~ .")),
    data = train_aggr[, c(name, cov_brt[-1])],
    num.trees = ntree,
    mtry = res_tune$recommended.pars$mtry,
    min.node.size = res_tune$recommended.pars$min.node.size,
    importance = "permutation",
    keep.inbag = FALSE
  )
  
  # Prédiction sur les données ponctuelles du fold
  test <- test_data[, cov_brt[-1]]
  test[cat_vars] <- lapply(test[cat_vars], factor)
  test[num_vars] <- lapply(test[num_vars], as.numeric)
  
  # Harmonisation des niveaux des facteurs
  for (v in cat_vars) {
    levels_train <- levels(train_aggr[[v]])
    test[[v]] <- factor(test[[v]], levels = levels_train)
  }
  
  preds <- predict(rf_model, data = test, num.threads = kmax)$predictions
  preds <- round(preds, 2)
  
  datacov$predRF_aggrcom[which(datacov$INSEE_COM %in% test_communes)] <- preds
  
  return(data.frame(fold = i, pred = preds, obs = test_data[[name]]))
}

# Évaluation finale
resuXvalRF_aggrcom <- Myeval(datacov$predRF_aggrcom, datacov[[name]])
resuXvalRF_aggrcom

stopImplicitCluster()
