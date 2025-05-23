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

varimp <- ggplot(Imp_RF, 
       aes(x=reorder(Vars, RF_Mod.G.variable.importance),
           y=RF_Mod.G.variable.importance #,  color=as.factor(var_categ)
       )
) + 
  geom_point() +
  geom_segment(aes(x=Vars,xend=Vars,y=0,yend=RF_Mod.G.variable.importance)) +
  scale_color_discrete(name="Variable Group") +
  ylab("IncNodePurity") +
  xlab("Variable Name") +
  ggtitle(name)+
  coord_flip()

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
                                                    type = "quantiles",
                                                    quantiles =  c(0.5),
                                                    num.threads = kmax )$prediction
            
            

          }


resuXvalQRF <-  Myeval(datacov$predRF,   datacov[,name] )



