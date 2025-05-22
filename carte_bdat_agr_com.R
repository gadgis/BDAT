#Carte type BDAT en faisant une agrégation par commune



#Chargement de package----
library(sf)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggmap)
library(ggspatial)
library(prettymapr)

setwd("Y:/BDAT/traitement_donnees/MameGadiaga/resultats")

#Importation des couches----
communes <- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\communes53_2014.shp")
dpt <- st_read("Y:\\BDAT\\traitement_donnees\\MameGadiaga\\prétraitement\\data\\dept53_rep.shp")

d_pH<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat.rds")

d_pH <- d_pH %>%
  select(id_profil,annee,x,y,pH,source,bdatid,insee, INSEE_COM,codification)

#transformation des données pH en exponnentiel
d_pH <- d_pH %>%
  mutate(H = 10^(-pH))
  

#exploration des données
summary(d_pH)

# #graphe de distribution (courbe) du pH
# ggplot(d_pH, aes(x = pH)) +
#   geom_density(color = "red", size = 1) +
#   labs(title = "Distribution du pH dans les données IGCS et BDAT",
#        x = "pH",
#        y = "Densité") +
#   theme_minimal()
# 
# d_pH2 <- d_pH %>% 
#   filter(source %in% c("BDAT", "IGCS")) %>% 
#   mutate(source = factor(source, levels = c("BDAT", "IGCS")))
# 
# # 2. Distribution pH par source
# ggplot(d_pH, aes(x = pH, color = source)) +
#   geom_density(alpha = 0.25, linewidth = 1) +     
#   scale_color_manual(values = c("#1f77b4", "#d62728"),  
#                      name   = "Source") +
#   scale_fill_manual(values  = c("#1f77b4", "#d62728"), 
#                     guide   = "none") +            
#   labs(title = "Distribution du pH selon la source",
#        x     = "pH",
#        y     = "Densité") +
#   theme_minimal() +
#   theme(plot.title = element_text(face = "bold", hjust = 0.5))

#Transformation en sf

d_pH_sf<-st_as_sf(d_pH, coords = c("x", "y"), crs = 2154)


pal <- brewer.pal(8, "Spectral")

brks <- c(-Inf, 5, 6, 6.3, 6.5, 6.8, 7, 8, +Inf)
labs <- c("[4.1 - 5[","[5 - 6[", "[6 – 6.3[", "[6.3 – 6.5[", "[6.5 – 6.8[", "[6.8 – 7[", "[7 - 8[", "[8 - 9.4]")

d_pH_sf <- d_pH_sf |>
  mutate(pH_cls = cut(pH, breaks = brks,
                      labels = labs,
                      right  = FALSE))

ggplot() +
  geom_sf(data = dept, fill = NA, colour = "grey70", linewidth = 0.3) +
  geom_sf(data = d_pH_sf,
          aes(colour = pH_cls),
          size = 1.22,
          alpha=0.9) +
  scale_colour_manual(values = pal, drop=F, name = "pH") +
  coord_sf() +
  labs(title   = "Distribution du pH") +
  theme_minimal() +
  theme(legend.position   = "right",
        legend.title      = element_text(size = 10),
        legend.text       = element_text(size = 8),
        plot.title        = element_text(face = "bold", hjust = 0.5))




#Agrégation des propriétés par commune
# 
# C_com <- d_C_sf %>%
#   group_by(INSEE_COM) %>%
#   summarise(
#     moy_C = round(mean(C, na.rm = TRUE), 1),
#     med_C = round(median(C, na.rm = TRUE), 1),
#     sd_C  = round(sd(C, na.rm = TRUE), 1),
#     n      = n()
#   )
# 
# C_com_df<-C_com %>%
#   st_drop_geometry()
# 
pH_com <- d_pH_sf %>%
  group_by(INSEE_COM) %>%
  summarise(
    moy_H = mean(H, na.rm = TRUE),
    med_H= median(H, na.rm = TRUE),
    sd_H  = sd(H, na.rm = TRUE),
    n      = n()
  )
#pH par commune comme log10 de la moy_H
pH_com <- pH_com %>%
  mutate(moy_pH = round(-log10(moy_H),1))%>%
  select(-sd_pH)
summary(pH_com)

pH_com_df<-pH_com %>%
st_drop_geometry()

#Jointure de l'agrégation" aux communes

mayenne_pH <- st_join(communes, pH_com, join = st_intersects)
# mayenne_C <- st_join(communes, C_com, join = st_intersects)


#Représentation graphique
#Carte de la moyenne de pH par commune
mayenne_pH <- mayenne_pH %>%
  mutate(ph_class = cut(
    med_pH,
    breaks = c(-Inf, 6, 6.3, 6.5, 6.8, 7, Inf),
    labels = c(
      "[5.7 - 6[",
      "[6 - 6.3[",
      "[6.3 - 6.5[",
      "[6.5 - 6.8[",
      "[6.8 - 7[",
      "[7 - 7.4["
    ),
    right  = FALSE
  ))
palet <- brewer.pal(6, "Spectral")
ggplot() +
  geom_sf(data = mayenne_pH, aes(fill = ph_class), color = "white") +
  scale_fill_manual(
    values = palet, name = "Classes pH",na.value = "grey80" ) +
  theme_minimal() +
  labs(
    title = "pH médian par commune" )

# #Carte de la moyenne de C par commune
# mayenne_C <- mayenne_C %>%                     
#   mutate(c_org_class = cut(
#     moy_C,                                  
#     breaks = c(-Inf, 15, 20, 25, 30, 35, Inf), 
#     labels = c(
#       "<15",            
#       "[15-20]",        
#       "]20-25]",        
#       "]25-30]",        
#       "]30-35]",        
#       ">35"             
#     )
#   ))
# 
# ggplot() +
#   geom_sf(data = mayenne_C, aes(fill = c_org_class), color = "white") +
#   scale_fill_manual(
#     values = c(
#       "<15"      = "#8c2d04",
#       "[15-20]"  = "#cc4c02",
#       "]20-25]"  = "#ec7014",
#       "]25-30]"  = "#fdae6b",
#       "]30-35]"  = "#a1d99b",
#       ">35"      = "#31a354"
#     ),
#     name = "Classes C org (g kg-¹)",
#     na.value = "grey80"
#   )+
#   theme_minimal() +
#   labs(
#     title = "Carbone moyen par commune" )

# 2. Table avec n points par commune (pH_com vient de ton code)

# 3. Palette continue bleu-vert (viridis option "mako")
pal1 <- viridis(100, option = "rocket", direction = -1)

# 4. Carte
ggplot() +
  geom_sf(data = mayenne_pH,
          aes(fill = n),
          colour = "grey30",    # fins contours
          linewidth = 0.05) +
  scale_fill_gradientn(
    colours  = pal1,
    name     = "Nombre\nde points",
    na.value = "grey90"         # gris pâle pour communes sans point
  ) +
  labs(title = "Effectif de points pH par commune") +
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 8),
    plot.title      = element_text(face = "bold", hjust = 0.5)
  )

#Carte de la repartition des points selon la source
ggplot() +
  geom_sf(data = d_pH_sf, aes(fill = source, color=source)) +
  scale_fill_manual(values = c("#1f77b4", "#d62728"), name = "source") +
  theme_minimal() +
  labs(
    title = "Origine des mesure de pH",
    ) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 12)
  )
  



#Export des fichiers----
saveRDS(C_com_df, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/C_moyen.rds")
saveRDS(pH_com_df, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/pH_median.rds")
