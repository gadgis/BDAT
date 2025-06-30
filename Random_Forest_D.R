library(ranger)
library(tuneRanger)
library(mlr)
library(foreach) 

cov_brt<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_cov_brt.rds")

#Création du tableau des covariables séléctionnées

cov_brt = c(name, cov_brt)
datacov_shrt = data_sample[,cov_brt]                   

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



print("Validation croisée----------------")

data_sample$predRF <- NA


resuXval <-
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {

            print(i)

            # collecter les # des lignes (gérer les doublons)
            nblignes = which( data_sample$id %in% data_sample$id[ fold[[i]] ] )
            
            fomula.ranger <- as.formula(paste0(name, "~ ."))
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


            data_sample$predRF[ nblignes ] <- predict(RF_Mod.G,
                                                    datacov_shrt[ nblignes , ],
                                                    num.threads = kmax )$predictions



          }


res_rf  <-  Myeval(data_sample$predRF,   data_sample[,name] )



