
# Chargement de package ----
library(ggmap)
library(ggspatial)
library(prettymapr)
library(soiltexture)
library(transformr)
library(gifski)

#Analyse des données BDAT ----

###Répartition par annee

####points en zone agricole par année
pts_bdat_agri_sf<-bdat_agri_sf %>%
  group_by(annee) %>%
  summarise(Effectif = n())

pts_bdat_agri<-pts_bdat_agri_sf %>%
  st_drop_geometry()
####points en zone urbaine par annee
pts_bdat_urb_sf<-bdat_urb_sf %>%
  group_by(annee) %>%
  summarise(Effectif = n())

pts_bdat_urb<-pts_bdat_urb_sf %>%
  st_drop_geometry()

##Représentation graphique ----

###Hisrogrammes ----

ggplot(pts_bdat_agri, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points BDAT en zone agricole par année", x="Année", y="Effectif") +
  theme_minimal()

ggplot(pts_bdat_urb, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points BDAT en zone urbaine par année", x="Année", y="Effectif") +
  theme_minimal()

###Cartes ----
ggplot() +
  annotation_map_tile(type = "osm") + 
  geom_sf(data = bdat_agri_sf, aes(color = "agri"), size = 2, alpha = 0.7) +
  geom_sf(data = bdat_urb_sf, aes(color = "urb"), size = 2, alpha = 0.7)+
  geom_sf(data = dpt, fill = NA, color = "black", size = 1.5) +
  scale_color_viridis_d(name = "Légende") + 
  labs(
    title = "Répartition des points BDAT catégorisés en fonction de l'ocsol",
    caption = "Source: Données BDAT et OCS_GE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering())

##carte distributiondes point agri par annee
ggplot() +
  annotation_map_tile(type = "osm") + 
  geom_sf(data = bdat_agri_sf, aes(color = as.factor(annee)), size = 2, alpha = 0.7) +
  geom_sf(data = dpt, fill = NA, color = "black", size = 1.5) +
  scale_color_viridis_d(name = "Année") + 
  labs(
    title = "Répartition des points BDAT en zone agricole par année",
    caption = "Source: Données BDAT et OCS_GE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering())


##Statistiques sur les propriétés des sol ----


summary(BDAT %>% select(pho, argi, sabt,limt, corgox, cecmet))

# distribution des propriétés des sols
#pH
ggplot(BDAT, aes(x = pho)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  labs(title = "Distribution du pH", x = "pH ", y = "Effectif") +
  theme_minimal()

#CARBONE
ggplot(BDAT, aes(x = corgox)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black") +
  labs(title = "Distribution du carbone", x = "Carbone ", y = "Effectif") +
  theme_minimal()

#CEC
ggplot(BDAT, aes(x = cecmet)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black") +
  labs(title = "Distribution de la CEC", x = "CEC ", y = "Effectif") +
  theme_minimal()

#triangle texturale BDAT

# Créer le df texture
texture<-BDAT %>%
  select(argi, limt, sabt) 

colnames(texture) <- c("CLAY", "SILT", "SAND")

texture <- texture %>% 
  na.omit() %>%
  mutate(total = CLAY + SILT + SAND) %>%
  mutate(
    CLAY = CLAY / total * 100,
    SILT = SILT / total * 100,
    SAND = SAND / total * 100
  ) %>%
  select(CLAY, SILT, SAND)

# Tracer le triangle textural
texture <- TT.normalise.sum(tri.data = texture, residuals = TRUE)


TT.plot(
  tri.data    = texture,
  class.sys   = "FR.GEPPA.TT",       
  main        = "Triangle textural BDAT",
  cex         = 1.2,               
  pch         = 22,             
  col         = "black",             
  bg          = "skyblue",  
  frame.bg.col = "white",       
  class.lab.col = "black",           
  grid.show    = TRUE                
)



#Analyse des données IGCS----

##Répartition par annee

#points en zone agricole par année
pts_igcs_agri_sf<-igcs_agri_sf %>%
  group_by(annee) %>%
  summarise(Effectif = n())

pts_igcs_agri<-pts_igcs_agri_sf %>%
  st_drop_geometry()
#points en zone urbaine par annee
pts_igcs_urb_sf<-igcs_urb_sf %>%
  group_by(annee) %>%
  summarise(Effectif = n())

pts_igcs_urb<-pts_igcs_urb_sf %>%
  st_drop_geometry()

##Représentation graphique ----


###Hisrogramme ----

ggplot(pts_igcs_agri, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points IGCS en zone agricole par année", x="Année", y="Effectif") +
  theme_minimal()

ggplot(pts_igcs_urb, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points IGCS en zone urbaine par année", x="Année", y="Effectif") +
  theme_minimal()

###Cartes ----
ggplot() +
  annotation_map_tile(type = "osm") + 
  geom_sf(data = igcs_agri_sf, aes(color = "agri"), size = 2, alpha = 0.7) +
  geom_sf(data = igcs_urb_sf, aes(color = "urb"), size = 2, alpha = 0.7)+
  geom_sf(data = dpt, fill = NA, color = "black", size = 1.5) +
  scale_color_viridis_d(name = "Légende") + 
  labs(
    title = "Répartition des points IGCS catégorisés en fonction de l'ocsol",
    caption = "Source: Données IGCS et OCS_GE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering())

##carte distribution des point agri par annee
ggplot() +
  annotation_map_tile(type = "osm") + 
  geom_sf(data = igcs_agri_sf, aes(color = as.factor(annee)), size = 2, alpha = 0.7) +
  geom_sf(data = dpt, fill = NA, color = "black", size = 1.5) +
  scale_color_viridis_d(name = "Année") + 
  labs(
    title = "Répartition des points IGCS en zone agricole par année",
    caption = "Source: Données BDAT et OCS_GE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) + 
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering())

##Statistique sur les propriétés des sol----

summary(igcs_agri)

igcs_agri<-igcs_agri %>%
  mutate(annee=as.factor(format(dt_cmp_,"%Y")),
         limon=1000-(argile+sand))

summary(igcs_agri %>% select(carbone, argile, limon, sand, cec, ph_eau,sand,mat_org,k_ech, mg_ech))

# distribution des propriétés des sols
#pH
ggplot(igcs_agri, aes(x = ph_eau)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "black") +
  labs(title = "Distribution du pH", x = "pH ", y = "Effectif") +
  theme_minimal()


ggplot(igcs_agri, aes(x = ph_eau)) +
  geom_density(fill = "steelblue", color = "black", alpha = 0.6) +
  labs(
    title = "Distribution du pH",
    x = "Carbone (%)",
    y = "Densité",
    caption = "Source : Données IGCS"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(size = 10)
  )


#CARBONE
ggplot(igcs_agri, aes(x = carbone)) +
  geom_histogram(binwidth = 1.5, fill = "steelblue", color = "black") +
  labs(title = "Distribution du carbone", x = "Carbone ", y = "Effectif") +
  theme_minimal()

ggplot(igcs_agri, aes(x = carbone)) +
  geom_density(fill = "steelblue", color = "black", alpha = 0.6) +
  labs(
    title = "Distribution de la teneur en carbone",
    x = "Carbone (%)",
    y = "Densité",
    caption = "Source : Données IGCS"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(size = 10)
  )

#CEC
ggplot(igcs_agri, aes(x = cec)) +
  geom_density(fill = "steelblue", color = "black", alpha = 0.6) +
  labs(
    title = "Distribution de la CEC",
    x = "CEC",
    y = "Densité",
    caption = "Source : Données IGCS"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(size = 10))
"scale_x_continuous(breaks = seq(0, 150, by = 2), limits = c(0, 150))"

#triangle texturale 

# Créer le df texture
texture_igcs <- igcs_agri %>%
  select(argile, limon, sand)%>%
  st_drop_geometry() 

colnames(texture_igcs) <- c("CLAY", "SILT", "SAND")

texture_igcs <- texture_igcs %>% 
  na.omit() %>%
  mutate(total = CLAY + SILT + SAND) %>%
  mutate(
    CLAY = CLAY / total * 100,
    SILT = SILT / total * 100,
    SAND = SAND / total * 100
  ) %>%
  select(CLAY, SILT, SAND)

# Tracer le triangle textural
texture_igcs <- TT.normalise.sum(tri.data = texture_igcs, residuals = TRUE)

TT.plot(
  tri.data    = texture_igcs,
  class.sys   = "FR.GEPPA.TT",       
  main        = "Triangle textural IGCS",
  cex         = 1.2,               
  pch         = 22,             
  col         = "black",             
  bg          = "skyblue",     
  frame.bg.col = "white",         
  class.lab.col = "black",           
  grid.show    = TRUE                
)

