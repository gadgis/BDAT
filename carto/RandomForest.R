#===================================================================================================================#
# Script :      Random Forest 

# Institution : UMR Infos&Sols /GisSol/BDAT

# Description : Script pour la prédiction des propriétés des sols en utilisant 
#               la méthode Random Forest avec INLA SPDE 

# Auteurs :     Mame Cheikh Gadiaga, Nicolas Saby

# Contact :     gadiagacheikh1998@gmail.com | nicolas.saby@inrae.fr

# Creation :    23-04-2025

# Entrees :     Observations ponctuelles et les covariables 

# Sorties :     Prediction de la propriété cible, indicateurs de performance, Ordre d'importance des covariables

# Modification : 06-10-2025
#===================================================================================================================#

#==========================================DEBUT DU SCRIPT=========================================================#

# Liste des packages utilisés -----

library(iml) # pour l'interprétabilité des modèles de machine learning
library(mlr) # pour la création et le tuning des modèles de machine learning

#1. Définir les types de variables----

X = datacov[,idcovs]                 # X: variables indépendantes -> covariables
Y = datacov[[idvar]]                  # Y: Variable cible (target variable, outcome) -

taille = ncol(X)

#2. Sélélction des covariables pertinantes par Boruta----

result_brt = Boruta(X, Y,                         
                    mtry = min(taille, floor(sqrt(taille))) ,
                    min.node.size = 3 ,
                    ntree = ntree)


Stats_brt = attStats(result_brt)                  #résultats du boruta
cov_brt = getSelectedAttributes(result_brt)       #sélection des covariables confirmées


result_brt_approche = TentativeRoughFix(result_brt)#résultats en forçant les indécis
cov_brt = getSelectedAttributes((result_brt_approche))      #sélection du nouvel ensemble

classement_brt = Stats_brt %>%
  arrange(desc(medianImp))                        #classe les covariables selon la valeur de la médiane

classement_brt_approche = attStats(result_brt_approche) %>%
  arrange(desc(medianImp))                        #pareil pour les covariables plus complètes

# sauvegarde des covariables sélectionnées
saveRDS(cov_brt, file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "_cov_brt.rds"))

#3. Création du jeu de données d'entrée pour la calibration du modèle

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
#                     min.node.size = res$recommended.pars$min.node.size ,
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
                    min.node.size = res$recommended.pars$min.node.size ,
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

saveRDS(Imp_RF, file=paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/",name,"_Imp_RF.rds"))

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
  ggtitle(paste("Dix premières variables importantes -",name)) +
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


terra::writeRaster(r, file = paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/", name, "qrf.tif"), overwrite = TRUE)


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
                                min.node.size = res$recommended.pars$min.node.size ,

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

