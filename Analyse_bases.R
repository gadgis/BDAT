install.packages(c("gganimate", "transformr", "gifski"))
install.packages("magick")


# Chargement de package ----
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggmap)
library(ggspatial)
library(writexl)
library(prettymapr)
library(soiltexture)
library(transformr)
library(gifski)
library(mpspline2)
library(purrr)
library(aqp)
"library(httpgd)
hgd()"
#importation des couches
setwd("C:/Users/gadiaga/Desktop/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees")
BDAT<-read.csv("bdat_53_x_y.csv", sep=",", header=TRUE)
igcs<-st_read("C:/Users/gadiaga/Desktop/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/DonneesIGCSStage.shp")
dpt <- st_read("C:/Users/gadiaga/Desktop/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/dept_53.shp")
com <- st_read("C:/Users/gadiaga/Desktop/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/COMMUNE.SHP")
ocsol<- st_read("C:/Users/gadiaga/Desktop/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/OCCUPATION_SOL.shp")
nomenclature<-read.csv("NomenclatureOCSGE.csv", sep=";", header=TRUE)

#correction des baseS

##BDAT
###selection des colonnes d'intéret
bdat<-BDAT %>%
  select(annee,insee,x,y,x_commune,y_commune,bdatid)

###Suppression dees lignes avec un X=NULL ou un Y=NULL
bdat<-bdat %>%
  filter(x!="NULL")
bdat<-bdat %>%
  filter(y!="NULL")

###Identification des différents SCR
bdat_L93<-bdat %>%
  filter(y>6000000)

bdat_L2E<-bdat %>%
  filter(y>1600000 & y<2700000)

bdat_wgs84<-bdat %>%
  filter(y>40 & y<50) 

bdat_invrs<-bdat %>%
  filter(y>=-2.45 & y<=7.23)

bdat_invrs<-bdat_invrs %>%
  rename(y=x, x=y)

######Effectifs pour les différents SCR
count(bdat_L2E)
count(bdat_wgs84)
count(bdat_L93)
count(bdat_invrs)

###Transformation en sf et reprojection en EPSG:2154
bdat_L93_sf<-st_as_sf(bdat_L93, coords = c("x","y"),crs=2154)

bdat_L2E_sf<-st_as_sf(bdat_L2E, coords = c("x","y"),crs=27572)
bdat_L2E_sf<-st_transform(bdat_L2E_sf, crs=2154)

bdat_wgs84_sf<-st_as_sf(bdat_wgs84, coords = c("x","y"),crs=4326)
bdat_wgs84_sf<-st_transform(bdat_wgs84_sf, crs=2154)

bdat_invrs_sf<-st_as_sf(bdat_invrs, coords = c("x","y"),crs=4326)
bdat_invrs_sf<-st_transform(bdat_invrs_sf, crs=2154)

bdat_sf<-rbind(bdat_wgs84_sf,bdat_L93_sf,bdat_invrs_sf,bdat_L2E_sf)

####Point BDAT dans le département de la Mayenne
bdat_53<-st_filter(bdat_sf, dpt, .predicate=st_within)
count(bdat_53)

###jointure spatiale bdat et communes
bdat_join<-bdat_53 %>%
  st_join(com, join = st_intersects)

bdat_join<-bdat_join %>%
  mutate(INSEE_COM=as.numeric(INSEE_COM),
         x_commune=as.numeric(x_commune),
         y_commune=as.numeric(y_commune),
         verif= (insee-INSEE_COM))

coords<- st_coordinates(bdat_join)
bdat_join<-bdat_join %>%
  mutate(x=coords[,1],
         y=coords[,2])

bdat_join<-bdat_join %>%
  mutate(d_com_bdat= (sqrt((x-x_commune)^2+(y-y_commune)^2)), d_com_bdat=round(d_com_bdat,0),
         d_com= (sqrt((x-X_CENTROID)^2+(y-Y_CENTROID)^2)), d_com=round(d_com,0))

bdat_join<-bdat_join %>%
  mutate(d_com_bdat=d_com_bdat/1000,
         d_com=d_com/1000)

bdat_INSEE_OK<-bdat_join %>%
  filter(verif==0)
count(bdat_INSEE_OK)

bdat_INSEE_ncf<-bdat_join %>%
  filter(verif!=0)
count(bdat_INSEE_ncf)

st_geometry(bdat_INSEE_ncf)


# Calculer la distance minimale et maximale entre chaque point et les limites de la commune qui lui est associée

bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  rowwise() %>%
  mutate(
    min_dist = min(as.numeric(
      st_distance(
        geometry,
        st_boundary(st_geometry(com[com$INSEE_COM == insee, ]))
      )
    )
  ))

#conversion en km

bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  mutate(min_dist = min_dist / 1000) 

#conservation des points avec une distance inférieure à 10 km
bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  filter(min_dist <= 10)

#fusion avec les point conforme
bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  select(-min_dist) 

bdat_ok<-rbind(bdat_INSEE_OK,bdat_INSEE_ncf)
bdat_ok<-bdat_ok %>%
  select(annee,insee,x,y,bdatid,INSEE_COM,geometry)
  

##Répartition des points suivant l'occupation du sol

bdat_ocsol_sf<- st_join(bdat_ok, ocsol, join = st_intersects)

######points en zone agricole

bdat_agri_sf <- bdat_ocsol_sf %>%
  filter(CODE_US=="US1.1")

bdat_agri_sf<- bdat_agri_sf %>%
  select(annee, insee, x, y, bdatid, INSEE_COM, CODE_US, geometry) 

count(bdat_agri_sf)

######points en zone urbaine

bdat_urb_sf<- bdat_ocsol_sf %>%
  filter(CODE_US =="US5")

bdat_urb_sf<- bdat_urb_sf %>%
  select(annee, insee, x, y, bdatid, INSEE_COM, CODE_US, geometry)

count(bdat_urb_sf)

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

##############Représentation graphique

########Hisrogramme

ggplot(pts_bdat_agri, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points BDAT en zone agricole par année", x="Année", y="Effectif") +
  theme_minimal()

ggplot(pts_bdat_urb, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points BDAT en zone urbaine par année", x="Année", y="Effectif") +
  theme_minimal()

############Cartes
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


##############Statistique sur les propriétés des sol

#####Remobilisation du tableau avec les propriétés des sols

BDAT<- BDAT %>%
  inner_join(bdat_agri_sf, by="bdatid")
#Transformation des champs en numérique
BDAT<-BDAT %>%
  mutate(
    pho = as.numeric(pho),
    argi = as.numeric(argi),
    limt= as.numeric(limt),
    sabt = as.numeric(sabt),
    corgox = as.numeric(corgox),
    cecmet = as.numeric(cecmet),
  )

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

texture <- texture %>% na.omit()

texture <- texture %>%
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
  class.sys   = "FR.GEPPA.TT",       # Système français
  main        = "Triangle textural BDAT",
  cex         = 1.2,               # Taille des points
  pch         = 22,             # Type de point (rond avec fond)
  col         = "black",             # Bordure
  bg          = "skyblue",     # Couleur de remplissage
  frame.bg.col = "white",         # Fond du graphique
  class.lab.col = "black",           # Couleur du nom des classes
  grid.show    = TRUE                # Affiche la grille
)



# IGCS ----

###selection des colonnes d'intéret

igcs$dt_cmp_ <- as.Date(igcs$dt_cmp_, format="%d/%m/%Y")

igcs_sf<-igcs %>%
  select(id_prfl, dt_cmp_)

igcs_sf<-igcs_sf %>%
  group_by(id_prfl,dt_cmp_) %>%
  summarise()

igcs_sf<-igcs_sf %>%
  mutate(annee=as.numeric(format(dt_cmp_, "%Y")))%>%
  st_transform(igcs_sf, crs=2154)



###Répartition des points suivant l'occupation du sol

igcs_ocsol_sf<- st_join(igcs_sf, ocsol, join = st_intersects)

######points en zone agricole

igcs_agri_sf <- igcs_ocsol_sf %>%
  filter(CODE_US=="US1.1")

igcs_agri_sf<- igcs_agri_sf %>%
  select(annee, id_prfl, CODE_US, geometry,dt_cmp_)

count(igcs_agri_sf)

######points en zone urbaine

igcs_urb_sf<- igcs_ocsol_sf %>%
  filter(CODE_US =="US5")

igcs_urb_sf<- igcs_urb_sf %>%
  select(annee, id_prfl, CODE_US, geometry,dt_cmp_)

count(igcs_urb_sf)

###Répartition par annee

####points en zone agricole par année
pts_igcs_agri_sf<-igcs_agri_sf %>%
  group_by(annee) %>%
  summarise(Effectif = n())

pts_igcs_agri<-pts_igcs_agri_sf %>%
  st_drop_geometry()
####points en zone urbaine par annee
pts_igcs_urb_sf<-igcs_urb_sf %>%
  group_by(annee) %>%
  summarise(Effectif = n())

pts_igcs_urb<-pts_igcs_urb_sf %>%
  st_drop_geometry()

##############Représentation graphique

########Hisrogramme

ggplot(pts_igcs_agri, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points IGCS en zone agricole par année", x="Année", y="Effectif") +
  theme_minimal()

ggplot(pts_igcs_urb, aes(x=annee, y=Effectif)) +
  geom_bar(stat="identity", fill="steelblue") +
  labs(title="Distribution des points IGCS en zone urbaine par année", x="Année", y="Effectif") +
  theme_minimal()

############Cartes
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

###carte animée


# Étape 1 : créer la colonne 'decennie'
igcs_agri_sf <- igcs_agri_sf %>%
  mutate(decennie = cut(annee, 
                        breaks = seq(floor(min(annee, na.rm = TRUE) / 10) * 10, 
                                     ceiling(max(annee, na.rm = TRUE) / 10) * 10, 
                                     by = 10),
                        include.lowest = TRUE,
                        right = FALSE,
                        labels = paste0(seq(floor(min(annee, na.rm = TRUE) / 10) * 10,
                                            ceiling(max(annee, na.rm = TRUE) / 10) * 10 - 10,
                                            by = 10),
                                        "–",
                                        seq(floor(min(annee, na.rm = TRUE) / 10) * 10 + 9,
                                            ceiling(max(annee, na.rm = TRUE) / 10) * 10 - 1,
                                            by = 10)
                        )
  ))

# Étape 2 : tracer l’animation avec gganimate
anim_plot <- ggplot() +
  annotation_map_tile(type = "osm") + 
  geom_sf(data = dpt, fill = NA, color = "black", size = 1) +
  geom_sf(data = igcs_agri_sf, aes(color = decennie), size = 2, alpha = 0.7) +
  scale_color_viridis_d(name = "Décennie") +
  labs(
    title = "Répartition des points IGCS par décennie",
    subtitle = "Décennie : {closest_state}",
    caption = "Source: Données BDAT et OCS_GE"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  ) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_fancy_orienteering()) +
  transition_states(decennie, transition_length = 2, state_length = 1) +
  ease_aes('linear')

# Étape 3 : exporter l’animation
animate(anim_plot, width = 800, height = 600, duration = 10, fps = 2, renderer = gifski_renderer("animation_igcs.gif"))


##############Statistique sur les propriétés des sol

#####Remobilisation du tableau avec les propriétés des sols

igcs_agri <- igcs %>%
  filter(id_prfl %in% igcs_agri_sf$id_prfl)

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

#triangle texturale BDAT

# Créer le df texture
texture_igcs <- igcs_agri %>%
  select(argile, limon, sand)

texture_igcs<-texture_igcs %>%
  st_drop_geometry() 

colnames(texture_igcs) <- c("CLAY", "SILT", "SAND")

texture_igcs <- texture_igcs %>% na.omit()

texture_igcs <- texture_igcs %>%
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

TT.plot(
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


#création de la base finale

##########BDAT

#selection des colonnes d'intéret BDAT

BDAT<-BDAT %>%
  select(bdatid,labo, no_campagne,annee.x,insee.x,INSEE_COM,canton,depart,region,CODE_US,x_commune,y_commune,
         pho, phk, argi, sabt, limt, corgox, cecmet,limf, limg, sabf,sabg,k2o,x.y, x.y, geometry)

#Renommer les colonnes
BDAT<-BDAT %>%
  rename(
    annee = annee.x,
    insee = insee.x,
    pH = pho,
    pHk = phk,
    CLAY = argi,
    SAND = sabt,
    SILT = limt,
    x=x.y,
    y=x.y,
    CEC = cecmet,
    C=corgox)

BDAT<-BDAT %>%
  select(-c(pHk,limf, limg, sabf,sabg))

summary(BDAT)
BDAT<-BDAT %>%
  mutate(
    annee=as.factor(annee),
    insee=as.factor(insee),
    k2o=as.numeric(k2o)
  )

#transformation en sf
BDAT<-st_sf(BDAT, geometry = BDAT$geometry, crs = 2154)

# Harmonisation IGCS ----------
#renommer les colonnes

igcs_agri<-igcs_agri %>%
  rename(
    pH = ph_eau,
    CEC = cec,
    C = carbone,
    CLAY = argile,
    SAND = sand,
    SILT = limon
    
  )

igcs_agri<-st_transform(igcs_agri, crs=2154)


cds<- st_coordinates(igcs_agri)
igcs_agri<-igcs_agri %>%
  mutate(x=cds[,1],
         y=cds[,2])

summary(igcs_agri)

#filtrer les horizons sans limites définies
igcs_agri<-igcs_agri %>%
  
  filter(!is.na(prof_nf),
         !is.na(prof_sp)
         )


#prédiction par des splines


#remplacement des NA

spl_df <- igcs_agri %>% 
  rename(id_profil = id_prfl,
         top       = prof_sp,
         bottom    = prof_nf
         ) %>%
  st_drop_geometry() %>%
  
# filtrage des profil composite de RMQS
  filter( is.na(typ_pr_)) %>%
  
  select(
    id_profil,
    top,
    bottom,
    no_hrzn,
    annee,
    x,
    y,
    pH

  ) %>%
  filter(!is.na(pH)) 

#filtrage des profils à horizons uniques
pf_uniq <- spl_df %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() == 1) %>%          
  ungroup()

#fitrage des profil en moyenne pondérée


#Filtrage des profils à slpiner
spl_dfs <- spl_df %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() > 3) %>%   #suppression des profils avec moins de 3 horizons       
  ungroup()



# Nettoyage (top & bottom) 
spl_dfs <- spl_dfs %>% 
  mutate(across(c(top, bottom), as.numeric)) %>% 
  arrange(id_profil, no_hrzn) %>% 
  group_by(id_profil) %>% 
  mutate(
   
    top = case_when(
      !is.na(top)                 ~ top,
      is.na(top) & no_hrzn == 1   ~ 0,
      TRUE                        ~ lag(bottom)
    ),
    
    bottom = coalesce(bottom, lead(top)),   
  ) %>% 
  ungroup() %>% 
  mutate(bottom = pmax(bottom, top))




spl_dat <- mpspline_tidy(obj = spl_dfs,
                         var_name = 'pH',
                         d = c(0, 30, 200),
                         lam = 0.05
                         )





tt =spl_dfs %>% filter(!is.na(pH    ))

table(tt$no_hrzn)



















#FUSION DES DEUX BASES
base_final <- bind_rows(igcs_agri, BDAT)
