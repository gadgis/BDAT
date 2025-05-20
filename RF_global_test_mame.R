#Chargement des packages----
library(sf)
library(randomForest)
library(ranger)
library(Boruta)
library(raster)
library(writexl)
library(foreach)
library(tidyverse)
library(tmap)
library(viridisLite) 

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

#changer le répertoire de travail
setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

# Définition des paramètres à tester-----
liste_ntree <- seq(50,500,50)

# Chargement des données----

dpt <- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/dept_53.shp")
gXY <- readRDS("gXY.rds")
rast_za <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_za.tif")
data<- readRDS("matrice_covariables.rds")

datacov_s <- data %>%
  select(-c(id_profil, bdatid, insee, codification)) %>%
  na.omit()

# datacov_s <- data %>%
#   select(-c(X,Y,INSEE_COM)) %>%
#   na.omit()

# définir les colonnes de covariables
covs <- datacov_s[, 7:70]  
idcovs <- names(covs)

datacov_s$id <- 1:nrow(datacov_s)

# Définir la variable réponse
idvar <- "pH"
Y <- datacov_s[[idvar]]
X <- datacov_s[, idcovs]

# Sélection d'un sous ensemble de variables -----
#Boruta sélectionne les variables pertinentes en se basant sur l’importance de permutation
# result_brt <- Boruta(X, Y)
# 
# plot(result_brt)
# 
# # Forcer une décision sur les variables incertaines
# result_brt_approche <- TentativeRoughFix(result_brt)
# 
# # Liste finale des covariables confirmées
# cov_brt <- getSelectedAttributes(result_brt_approche)
# # 
# # #enregistrer la liste des covariables
# save(cov_brt_C, file = "cov_brt_C.Rds")
load("liste_variables_boruta.rds")
load("cov_brt_C.rds")

# créer un sous-jeu de données avec les variables retenues par Boruta
# X2 <-  datacov_s[, c("id", "moy_ph", cov_brt)]
X2 <-  datacov_s[, c("id", "pH", cov_brt)]

# Calibrage des hyperparamètres -----

## mettre en place la parallélisation (sera utilisée pour la kCV) 
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

##Calibrage du modèle
res_calib <- foreach(
  nt        = liste_ntree,               # valeur de ntree (une par tâche)
  .combine  = rbind,
  .packages = c("randomForest", "ranger", "pls"),
  .export   = c("X2", "Y", "taille", "folds", "Myeval")
) %dopar% {

  ##Étape 1 : Tune mtry 
  bestmtry <- tuneRF(
    x          = X2[, -c(1, 2)],
    y          = Y,
    mtryStart  = round(sqrt(taille)),
    stepFactor = 1.3,
    improve    = 1e-5,
    ntreeTry   = nt,
    trace      = FALSE,
    plot       = FALSE
  )
  mtry_opt <- bestmtry[which.min(bestmtry[, "OOBError"]), "mtry"]

  ## Étape 2 : 10-fold CV 
  preds <- numeric(length(Y))              # vecteur final de prédictions

  for (seg in folds) {
    mod <- ranger(
      formula      = C ~ .,
      data         = X2[-seg, -1],         
      num.trees    = nt,
      mtry         = mtry_opt,
      min.node.size = 3,
      quantreg     = TRUE,
      max.depth    = 15,
      importance   = "permutation",
      num.threads  = 1,
      keep.inbag   = FALSE
    )

    preds[seg] <- predict(
      mod,
      data       = X2[seg, -c(1, 2)],      
      type       = "quantiles",
      quantiles  = 0.5
    )$predictions
  }

  ## Étape 3 : métriques
  metrique         <- Myeval(preds, Y)
  metrique$ntree   <- nt
  metrique$mtry    <- mtry_opt
  metrique                              
}

saveRDS(res_calib, file = "resultats_calib_c.rds")

# res_calib_list <- lapply(liste_ntree, function(nt) {
#   
#   ## Étape 1 :  tuning de mtry 
#   bestmtry <- tuneRF(
#     x          = X2[, -c(1, 2)],
#     y          = Y,
#     mtryStart  = round(sqrt(taille)),
#     stepFactor = 1.3,
#     improve    = 1e-5,
#     ntreeTry   = nt,
#     trace      = FALSE,
#     plot       = FALSE
#   )
#   mtry_opt <- bestmtry[which.min(bestmtry[, "OOBError"]), "mtry"]
#   
#   ## -- Étape 2 : 10-fold CV, exécutée en parallèle 
#   res <- foreach(
#     seg      = folds,
#     .combine = rbind,
#     .packages = "ranger",
#     .export   = c("X2", "nt", "mtry_opt")
#   ) %dopar% {
#     mod <- ranger(
#       formula        = pH ~ .,
#       data           = X2[-seg, -1],     # apprentissage (sans colonne id)
#       num.trees      = nt,
#       mtry           = mtry_opt,
#       min.node.size  = 3,
#       quantreg       = TRUE,
#       max.depth      = 15,
#       importance     = "permutation",
#       num.threads    = 1
#     )
#     preds <- predict(mod,
#                      data= X2[seg, -c(1,2)],  
#                      type = "quantiles",
#                      quantiles = 0.5,
#                      num.threads = parallel::detectCores())$predictions
#     
#     preds
#   }
 # recoller les résultats
#   res_df <- data.frame(observation = Y, prediction = NA) 
#   for(k in 1:10){
#     res_df$prediction[folds[[k]]] <- res[[k]]
#   }
#   
#   ##Étape 3 : métriques 
#   metrique         <- Myeval(x = res_df$prediction, y = res_df$observation)
#   metrique$ntree   <- nt
#   metrique$mtry    <- mtry_opt
#   return(metrique)
# })
# 
# res_calib<-bind_rows(res_calib_list)

# Libération des ressources
parallel::stopCluster(my.cluster)

resC<-readRDS("resultats_calib_C.rds")
respH<-readRDS("resultats_calib_pH.rds")
respH_cent<-readRDS("result_calib_pH_centroides.rds")

# saveRDS(res_calib, file = "resultats_calib_C.rds")

 ### Le modèle retenu est celui avec ntree=500 et mtry_opt=14

# Visualisation de l'importance des variables----

##Pour le pH

  QRF_Mod.G <- ranger(
    formula = pH ~ .,
    data = X2[, !(names(X2) %in% c("id"))],
    num.trees = 350,
    min.node.size = 3,
    quantreg = TRUE,
    max.depth = 15,
    mtry = 11,
    importance = "permutation",
    scale.permutation.importance = TRUE
  )

  Imp_QRF <- data.frame(Importance = QRF_Mod.G$variable.importance,
                        Variable = names(QRF_Mod.G$variable.importance)) %>%
    arrange(desc(Importance))

  plot_path <- paste0("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/variables_importantes.png")
  png(plot_path, width = 800, height = 600)
  print(
    ggplot(Imp_QRF, aes(x = reorder(Variable, Importance), y = Importance)) +
      geom_point() +
      geom_segment(aes(xend = Variable, yend = 0)) +
      coord_flip() +
      ylab("Importance (Permutation)") +
      xlab("Covariables") +
      ggtitle(paste("Importance des covariables"))
  )
  dev.off()

#Application du modèle----
  
  #Sélection des covariables spatiales (mêmes que pour le modèle)
  testD <- gXY %>% 
    select(all_of(cov_brt))
  
  # Prédiction des quantiles (5%, médiane, 95%)
  QRF_Median <- predict(
    QRF_Mod.G,               
    data = testD,            
    type = "quantiles",
    quantiles = c(0.05, 0.5, 0.95),  
    num.threads = parallel::detectCores() - 1  
  )$prediction

  #Fusion avec les coordonnées spatiales

  QRF_Median50 <- bind_cols(
    gXY %>% 
      filter_at(vars(all_of(cov_brt)), all_vars(!is.na(.))) %>% 
      select(x, y),      
    data.frame(                 
      Q5 = QRF_Median[, 1],     
      Median = QRF_Median[, 2], 
      Q95 = QRF_Median[, 3]     
    )
  )
  
  # Création d'un raster pour chaque quantile
  rast_q5 <- rast(QRF_Median50 %>% select(x, y, Q5), type = "xyz", crs= "EPSG:2154")
  rast_median <- rast(QRF_Median50 %>% select(x, y, Median), type = "xyz",crs= "EPSG:2154")
  rast_q95 <- rast(QRF_Median50 %>% select(x, y, Q95), type = "xyz", crs= "EPSG:2154")
  crs(rast_za) <- "EPSG:2154"
  
  rast_q5 <- extend(rast_q5, rast_za, snap = "near")
  rast_median <- extend(rast_median, rast_za, snap = "near")
  rast_q95 <- extend(rast_q95, rast_za, snap = "near")
  
  #Application du masque aux rasters
  rast_q5 <- mask(rast_q5, rast_za)
  rast_median <- mask(rast_median, rast_za)
  rast_q95 <- mask(rast_q95, rast_za)
  
  # Sauvegarde des rasters
  writeRaster(rast_q5, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_q5.tif", filetype="GTiff",overwrite=TRUE)
  writeRaster(rast_median, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_median.tif", filetype="GTiff",overwrite=TRUE)
  writeRaster(rast_q95, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_q95.tif", filetype="GTiff",overwrite=TRUE)
  
#charger raster 
rast_median_cent <- rast("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/rast_median_centr.tif")
#représentation graphique----
  
  brk   <- c(-Inf, 6, 6.3, 6.5, 6.8, 6.7,  Inf)          
  labs  <- c("[4.3 - 6[",
             "[6 – 6.3[",
             "[6.3 – 6.5[",
             "[6.5 – 6.8[",
             "[6.8 – 7[",
             "[7 - 8.5]")
  
 pal <- brewer.pal(6, "Spectral")
 
#carte de la médiane avec tmap en faisant des breaks

 df_med <- as.data.frame(rast_median_cent, xy = TRUE, na.rm = TRUE)
 names(df_med)[3] <- "pH"
 
 df_med$pH_cls <- cut(df_med$pH,
                      breaks = brk,
                      labels = labs,
                      right  = FALSE)      
 

 ggplot(df_med) +
   geom_raster(aes(x, y, fill = pH_cls)) +
   scale_fill_manual(values = pal,
                     drop   = FALSE,       
                     name   = "pH") +
   coord_equal() +
   labs(title = "Prédiction de la valeur médiane du pH") +
   theme_minimal() +
   theme(
     legend.position = "right",
     legend.title    = element_text(face = "bold"),
     legend.text     = element_text(size = 8),
     plot.title      = element_text(face = "bold", hjust = 0.5)
   )
########################################
 
 df_cent <- as.data.frame(rast_median_cent, xy = TRUE, na.rm = TRUE) |>
   rename(pH = 3) |>
   mutate(pH_cls = cut(pH, breaks = brk, labels = labs, right = FALSE),
          raster  = "Modèle_centroïdes")
 
 
 df_median  <- as.data.frame(rast_median, xy = TRUE, na.rm = TRUE) |>
   rename(pH = 3) |>
   mutate(pH_cls = cut(pH, breaks = brk, labels = labs, right = FALSE),
          raster  = "Modèle_données ponctuelles ") 
 
 df_all  <- bind_rows(df_cent, df_median) 

 ggplot(df_all) +
   geom_raster(aes(x, y, fill = pH_cls)) +
   facet_wrap(~ raster, ncol = 2) +          
   scale_fill_manual(values = pal,
                     drop   = FALSE,
                     name   = "pH") +
   coord_equal() +
   labs(title = "Prédictions médianes du pH – comparaison de 2 modèles") +
   theme_minimal() +
   theme(
     legend.position = "right",
     legend.title    = element_text(face = "bold"),
     legend.text     = element_text(size = 8),
     strip.text      = element_text(face = "bold"),
     plot.title      = element_text(face = "bold", hjust = 0.5)
   )
 