#chargement des packages
library(sf) # manipuler les datas SIG vector
library(tmap)  # outil d'édution de carte
library(readxl) # lecture de fichier excel
library(tidyr) # reroganisation des données
library(dplyr)  # grammaire de manipulation de données
library(terra) # manipuler des données raster-pixels
library(ggplot2) # graphiques
library(sp) # ancetre de sf. Pourquoi l'utiliser ?
library(randomForest) # randomForest pour utiliser tuneRandomForest
library(Boruta)  # implementation  de selection de variables de randomForest
library(ranger) # quantile randomForest, version optimée de randomForest
library(caret) # modelisation creation de modeles prédictifs avec tuning validation-croisée
library(foreach) # boucles optimisées de R plus facile à coder voir option combine
library(raster) # manipuler des données raster

library(gstlearn) # géostatistique
library(ggpubr) # pour les graphiques
library(ggnewscale)
library(doParallel)

library(raster) # copie de terra mais pour le package OGC
#library(OGC) # Pour le calcul des coordonnées obliques
library(tuneRanger)

# INLA SPDE avec inla bru, approche bayesienne de la géostatistique
library(INLA)
library(inlabru)
library(RColorBrewer) 
library(patchwork) 

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

#chargement des résultats sur les metriques
metriqueRF<-readRDS("metrique_qrf.rds")
metriqueRF_cent<- readRDS("metrique_QRF_centroides.rds")
metriqueKO<- readRDS("metrique_KO.rds")
metriqueKO_cent<- readRDS("metrique_KO_centroides.rds")
metriqueKED<- readRDS("metrique_KED.rds")
metriqueKED_cent<- readRDS("metrique_KED_centroides.rds")

#Chargement des données rasters

RF<-rast("pHqrf_final.tif")
RF_cent<-rast("pHqrf_centroides_final.tif")
KO<-rast("predKOINLA_final.tif")
KO_cent<-rast("predKOINLA_centroides_final.tif")
KED<-rast("predKEDINLA_final.tif")
KED_cent<-rast("predKEDINLA_centroides_final.tif")

#metriques finales----
Metrique<-rbind(metriqueRF,metriqueRF_cent,metriqueKO,metriqueKO_cent,metriqueKED,metriqueKED_cent)

modele<-c("RF","RF_centroides","KO","KO_centroides","KED","KED_centroides","KED_total")
Metrique$modele<-modele

#Représentation des rasters----

##RF----

# Intervalles de classes
brk <- c(-Inf, 6, 6.2, 6.3, 6.4, 6.5, Inf)

# Étiquettes de classes
labs <- c("[5.7 - 6[",
          "[6 – 6.2[",
          "[6.2 – 6.3[",
          "[6.3 – 6.4[",
          "[6.4 – 6.5[",
          "[6.5 - 6.6]")

# Palette Spectral (inversée pour ressembler à ton image)
colors <- rev(RColorBrewer::brewer.pal(length(labs), "Spectral"))

# Affecter les noms de classes comme levels (important pour la légende)

#conversion raster en df pour ggplot
df_raster <- as.data.frame(RF_cent, xy = TRUE, na.rm = TRUE)
colnames(df_raster)[3] <- "pH"
df_raster <- df_raster %>%
  mutate(pH = round(pH, 2))  

# Classification
df_raster <- df_raster %>%
  mutate(pH_class = cut(pH, breaks = brk, labels = labs, include.lowest = TRUE))

# ggplot
ggplot(df_raster, aes(x = x, y = y, fill = pH_class)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values = rev(colors), name = "pH", drop = FALSE) +
  labs(title = "Prédiction pH - RF_centroÏdes") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

##KO----

# Intervalles de classes
brk <- c(-Inf, 6, 6.2, 6.3, 6.4, 6.5, Inf)

# Étiquettes de classes
labs <- c("[5.7 - 6[",
          "[6 – 6.2[",
          "[6.2 – 6.3[",
          "[6.3 – 6.4[",
          "[6.4 – 6.5[",
          "[6.5 - 6.6]")

# Palette Spectral (inversée pour ressembler à ton image)
colors <- rev(RColorBrewer::brewer.pal(length(labs), "Spectral"))

# Affecter les noms de classes comme levels (important pour la légende)

#conversion raster en df pour ggplot
df_raster <- as.data.frame(KO_cent, xy = TRUE, na.rm = TRUE)
colnames(df_raster)[3] <- "pH"
df_raster <- df_raster %>%
  mutate(pH = round(pH, 2))  

# Classification
df_raster <- df_raster %>%
  mutate(pH_class = cut(pH, breaks = brk, labels = labs, include.lowest = TRUE))

# ggplot
ggplot(df_raster, aes(x = x, y = y, fill = pH_class)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values = rev(colors), name = "pH", drop = FALSE) +
  labs(title = "Prédiction pH - KO_centroÏdes") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
##KED----
# Intervalles de classes
brk <- c(-Inf, 6, 6.2, 6.5, 6.7, 6.9, Inf)

# Étiquettes de classes
labs <- c("[5.2 - 6[",
          "[6 - 6.2[",
          "[6.2 – 6.5[",
          "[6.5 – 6.7[",
          "[6.7 – 6.9[",
          "[6.9 - 7]")

# Palette Spectral (inversée pour ressembler à ton image)
colors <- rev(RColorBrewer::brewer.pal(length(labs), "Spectral"))

# Affecter les noms de classes comme levels (important pour la légende)

#conversion raster en df pour ggplot
df_raster <- as.data.frame(KED_cent, xy = TRUE, na.rm = TRUE)
colnames(df_raster)[3] <- "pH"
df_raster <- df_raster %>%
  mutate(pH = round(pH, 2))  

# Classification
df_raster <- df_raster %>%
  mutate(pH_class = cut(pH, breaks = brk, labels = labs, include.lowest = TRUE))

# ggplot
ggplot(df_raster, aes(x = x, y = y, fill = pH_class)) +
  geom_raster() +
  coord_equal() +
  scale_fill_manual(values = rev(colors), name = "pH", drop = FALSE) +
  labs(title = "Prédiction pH - KED_centroÏdes") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

#Comparaison des prédictions----

# Charger tous les rasters
# Charger les rasters
r_stack <- rast(c("pHqrf_centroides_final.tif",
                  "predKOINLA_centroides_final.tif",
                  "predKEDINLA_centroides_final.tif"))
names(r_stack) <- c("RF_Centroides", "INLA_KO_Centroides", "INLA_KED_Centroides")

# Convertir en data.frame long
r_df <- as.data.frame(r_stack, xy = TRUE) %>%
  pivot_longer(cols = -c(x, y), names_to = "Méthode", values_to = "pH") %>%
  mutate(pH = round(pH, 2))

# Définir les classes, labels et couleurs
brk <- c(-Inf, 6, 6.2, 6.4, 6.6, 6.8,  Inf)

labs <- c("[5.2 - 6[",
          "[6 – 6.2[",
          "[6.2 – 6.4[",
          "[6.4 – 6.6[",
          "[6.6 – 6.8[",
          "[6.8 - 7]")

colors <- rev(brewer.pal(length(labs), "Spectral"))

# Créer une variable catégorisée
r_df <- r_df %>%
  mutate(pH_class = cut(pH, breaks = brk, labels = labs, include.lowest = TRUE))

# Tracer avec ggplot2
ggplot(r_df, aes(x = x, y = y, fill = pH_class)) +
  geom_raster() +
  coord_equal() +
  facet_wrap(~ Méthode) +
  scale_fill_manual(
    values = rev(colors),
    name = "pH",
    drop = FALSE
  ) +
  labs(title = "Comparaison des prédictions de pH selon les méthodes") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title = element_blank(),
    legend.position = "bottom"
  )
#Graphique sur les prédictions de pH selon les méthodes----

points<-readRDS("results_pts.rds")
centroides<-readRDS("results_centroides.rds")

##distribution des prédictions----
points_long <- points %>%
  dplyr::select(pH,predRF,predINLAKO, predINLAKED) %>%
  dplyr::rename(observed=pH,
         RF=predRF,
         INLAKO=predINLAKO,
         INLAKED=predINLAKED)%>%
  tidyr::pivot_longer(cols = everything(), names_to = "Variable", values_to = "value")

# Tracer la densité
ggplot(points_long, aes(x = value, fill = Variable, color = Variable)) +
  geom_density(alpha = 0.4,adjust = 1.5) +
  theme_minimal() +
  labs(
    title = "Densités des prédictions et des observations du pH",
    x = "pH",
    y = "Densité"
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
    
  )
centroides_long <- centroides %>%
  dplyr::select(moy_pH,predRF,predINLAKO, predINLAKED) %>%
  dplyr::rename(observed=moy_pH,
         RF=predRF,
         INLAKO=predINLAKO,
         INLAKED=predINLAKED)%>%
  tidyr::pivot_longer(cols = everything(), names_to = "Variable", values_to = "value")
# Tracer la densité
ggplot(centroides_long, aes(x = value, fill = Variable, color = Variable)) +
  geom_density(alpha = 0.4, adjust = 1.5) +
  theme_minimal() +
  labs(
    title = "Densités des prédictions et des observations du pH",
    x = "pH",
    y = "Densité"
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
    
  )

# Graphique de comparaison entre les valeurs observées et prédites du pH----
points_long <- points %>%
  dplyr::select(observed = pH, predRF, predINLAKO, predINLAKED) %>%
  dplyr::rename(
    RF = predRF,
    INLAKO = predINLAKO,
    INLAKED = predINLAKED
  ) %>%
  tidyr::pivot_longer(
    cols = c(RF, INLAKO, INLAKED),  
    names_to = "Variable",
    values_to = "value"
  )


ggplot(points_long, aes(x = value, y = observed, color = Variable)) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  facet_wrap(~ Variable) +
  theme_minimal() +
  labs(
    x = "Valeurs prédites",
    y = "Valeurs observées",
    title = "Comparaison entre valeurs observées et prédites du pH"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
