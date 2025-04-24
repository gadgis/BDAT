
# Chargement de package ----
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mpspline2)

#Importation des couches ----
setwd("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees")
BDAT<-read.csv("bdat_53_x_y.csv", sep=",", header=TRUE)
igcs<-st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/DonneesIGCSStage.shp")
dpt <- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/dept_53.shp")
com <- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/COMMUNE.SHP")
ocsol<- st_read("Y:/BDAT/traitement_donnees/MameGadiaga/prétraitement/Analyse/Codes_Mame/Donnees/OCCUPATION_SOL.shp")
nomenclature<-read.csv("NomenclatureOCSGE.csv", sep=";", header=TRUE)

#Pretraitement sur la BDAT ----

##Correction de la base ----

###Correction reprojection et selection des points se trouvant en Mayenne----

####Selection des colonnes d'intéret
bdat<-BDAT %>%
  select(annee,insee,x,y,x_commune,y_commune,bdatid)

####Suppression dees lignes avec un X=NULL ou un Y=NULL
bdat<-bdat %>%
  filter(x!="NULL", y!="NULL")

####Identification des différents SCR 
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

####Transformation en sf et reprojection en EPSG:2154 
bdat_L93_sf<-st_as_sf(bdat_L93, coords = c("x","y"),crs=2154)

bdat_L2E_sf<-st_as_sf(bdat_L2E, coords = c("x","y"),crs=27572)%>%
  st_transform(crs=2154)

bdat_wgs84_sf<-st_as_sf(bdat_wgs84, coords = c("x","y"),crs=4326)%>%
  st_transform(crs=2154)

bdat_invrs_sf<-st_as_sf(bdat_invrs, coords = c("x","y"),crs=4326)%>%
  st_transform(crs=2154)

####Fusion des sous bases de SCR différents
bdat_sf<-rbind(bdat_wgs84_sf,bdat_L93_sf,bdat_invrs_sf,bdat_L2E_sf)

####Extraction des points se trouvant dans le département de la Mayenne
bdat_53<-st_filter(bdat_sf, dpt, .predicate=st_within)

###Verification de l'exactitude des codes INSEE----

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

# Nous allons avoir deux lots de points:
# un premier dont lINSEE est conforme que nous codifions 1
# un second dont l'INSEE n'est pas conforme codifié 2 

###Filrage des points avec un code INSEE conforme
bdat_INSEE_OK<-bdat_join %>%
  filter(verif==0)%>%
  mutate(codification =1)

###Filtrage des points avec un code INSEE non conforme
bdat_INSEE_ncf<-bdat_join %>%
  filter(verif!=0)%>%
  mutate(codification =2)

###Calculer la distance minimale entre chaque point de bdat_INSEE_ncf et les limites 
###de la commune qui lui est associée

bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  rowwise() %>%
  mutate(
    min_dist = min(as.numeric(
      st_distance(
        geometry,
        st_boundary(st_geometry(com[com$INSEE_COM == insee, ]))
      )
    )
  ))%>%

#conversion en km
  mutate(min_dist = min_dist / 1000)%>%

#conservation des points avec une distance inférieure à 10 km
  filter(min_dist <= 10)


###Fusion avec les point conforme
bdat_INSEE_ncf <- bdat_INSEE_ncf %>%
  select(-min_dist) 

bdat_ok<-rbind(bdat_INSEE_OK,bdat_INSEE_ncf)
bdat_ok<-bdat_ok %>%
  select(annee,insee,x,y,bdatid,INSEE_COM,codification, geometry)
  

##Identification des points en zone agricole et urbaine----

bdat_ocsol_sf<- st_join(bdat_ok, ocsol, join = st_intersects)

######points en zone agricole

bdat_agri_sf <- bdat_ocsol_sf %>%
  filter(CODE_US=="US1.1")

bdat_agri_sf<- bdat_agri_sf %>%
  select(annee, insee, x, y, bdatid, INSEE_COM, CODE_US,codification, geometry) 

######points en zone urbaine

bdat_urb_sf<- bdat_ocsol_sf %>%
  filter(CODE_US =="US5")

bdat_urb_sf<- bdat_urb_sf %>%
  select(annee, insee, x, y, bdatid, INSEE_COM, CODE_US, codification, geometry)



##Remobilisation du tableau avec les propriétés des sols----

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



# Pretraitement sur IGCS ----

##Identification des points en zone agricole et urbaine----

###selection des colonnes d'intéret

igcs$dt_cmp_ <- as.Date(igcs$dt_cmp_, format="%d/%m/%Y")

igcs_sf<-igcs %>%
  select(id_prfl, dt_cmp_)%>%
  group_by(id_prfl,dt_cmp_) %>%
  summarise()

igcs_sf<-igcs_sf %>%
  mutate(annee=as.factor(format(dt_cmp_, "%Y")))%>%
  st_transform(igcs_sf, crs=2154)

###Répartition des points suivant l'occupation du sol

igcs_ocsol_sf<- st_join(igcs_sf, ocsol, join = st_intersects)

######points en zone agricole

igcs_agri_sf <- igcs_ocsol_sf %>%
  filter(CODE_US=="US1.1")

igcs_agri_sf<- igcs_agri_sf %>%
  select(annee, id_prfl, CODE_US, geometry,dt_cmp_)


######points en zone urbaine

igcs_urb_sf<- igcs_ocsol_sf %>%
  filter(CODE_US =="US5")

igcs_urb_sf<- igcs_urb_sf %>%
  select(annee, id_prfl, CODE_US, geometry,dt_cmp_)



##Remobilisation du tableau avec les propriétés des sols ----

igcs_agri <- igcs %>%
  filter(id_prfl %in% igcs_agri_sf$id_prfl)%>%
  mutate(annee=as.factor(format(dt_cmp_, "%Y")))

##Harmonisation IGCS à 0-30 cm ----
#renommer les colonnes

igcs_agri<-igcs_agri %>%
  mutate(limon=1000-(argile+sand))%>%
  rename(
    pH = ph_eau,
    CEC = cec,
    C = carbone,
    CLAY = argile,
    SAND = sand,
    SILT = limon )%>%
  st_transform(igcs_agri, crs=2154)


cds<- st_coordinates(igcs_agri)
igcs_agri<-igcs_agri %>%
  mutate(x=cds[,1],
         y=cds[,2])

summary(igcs_agri)

###Filtrage des horizons sans limites définies----

igcs_agri<-igcs_agri %>%
  
  filter(!is.na(prof_nf),
         !is.na(prof_sp)
  )

###transformation en dataframe
igcs_agri_df <- igcs_agri %>% 
  rename(id_profil = id_prfl,
         top       = prof_sp,
         bottom    = prof_nf
  ) %>%
  st_drop_geometry()

###Filtrage des profils composite de RMQS---

#les profils RMQS sont ceux dont typ_pr = C ou F

igcs_agri_df<-igcs_agri_df%>%
  filter( is.na(typ_pr_)) 

###Harmonisation des données sur le pH----

igcs_agri_pH<-igcs_agri_df%>%
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

####selection des profils à horizons uniques
pf_uniq <- igcs_agri_pH %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() == 1) %>%          
  ungroup()

####selection des profiles à deux horizons dont la limite inférieur du premier
#horizon est comprise entre 30 et 35 cm
P2hrzn1_30_35 <- igcs_agri_pH %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() == 2) %>%          
  ungroup() %>%
  filter(no_hrzn==1, bottom >= 30 & bottom <= 35)


####selection des profils à deux horizons dont le premier horizon ne fait pas 30 cm 

pf_mp <- igcs_agri_pH %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() == 2) %>%          
  ungroup()%>%
  anti_join(                                   # garde ce qui n’est PAS joint
    P2hrzn1_30_35 %>% 
      distinct(id_profil),                    # id uniques à exclure
    by = "id_profil"
  )

p_mp <- pf_mp %>% 
  mutate(
    thick = bottom - top          # épaisseur cm
  ) %>% 
  group_by(id_profil,annee,x,y) %>% 
  summarise(
    pH = round(weighted.mean(pH, thick, na.rm = TRUE),1)   # moyenne pondérée
  ) %>% 
  ungroup()

####selection des profils dont les horizons sont supérieure ou égales à 3 et que le premier
#horizon est compris entre 30 et 35 cm
p3hrzn1_30_35 <- igcs_agri_pH %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() >= 3) %>%          
  ungroup() %>%
  filter(no_hrzn==1, bottom >= 30 & bottom <= 35)


####selection des profils à slpiner
spl_dfs <- igcs_agri_pH %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() >= 3) %>%   #suppression des profils avec moins de 3 horizons       
  ungroup()%>%
  anti_join(                                   # garde ce qui n’est PAS joint
    p3hrzn1_30_35 %>% 
      distinct(id_profil),                    # id uniques à exclure
    by = "id_profil"
  )



# Nettoyage (top & bottom) 


spl_layers <- mpspline(
  obj = spl_dfs,
  var_name = "pH",
  d = c(0, 30),
  lam = 0.1  
)



# Étape 1: Créer un dataframe à partir de la liste
igcs_spline <- imap_dfr(spl_layers, ~ {
  
  tibble(
    # Identifiant et métadonnées de base
    id_profil = .x$id_profil,
    
    # Horizon standard (est_dcm)
    pH_000_030_cm = round(.x$est_dcm[["000_030_cm"]],1),
    
    # Erreurs (est_err)
    RMSE = .x$est_err[["RMSE"]],
    RMSE_IQR = .x$est_err[["RMSE_IQR"]])
  
})

# Étape 2: Fusionner avec les métadonnées supplémentaires (

igcs_spline <- igcs_spline %>%
  left_join(
    spl_dfs %>% distinct(id_profil, annee, x, y),
    by = "id_profil"
  ) %>%
  rename(pH = pH_000_030_cm) %>%
  select(id_profil, annee, x, y, pH)

##haromisation des colonnes
pf_uniq <- pf_uniq %>%
  select(id_profil, annee, x, y, pH)

p3hrzn1_30_35 <- p3hrzn1_30_35 %>%
  select(id_profil, annee, x, y, pH)

P2hrzn1_30_35 <- P2hrzn1_30_35 %>%
  select(id_profil, annee, x, y, pH)

igcs_final_pH<-rbind(p_mp, pf_uniq, p3hrzn1_30_35, igcs_spline, P2hrzn1_30_35)

igcs_final_pH <- igcs_final_pH %>%
  mutate(
    source="IGCS"
  )
#Création de la base finale----

#BDAT

#selection des colonnes d'intéret BDAT_pH

BDAT_pH<-BDAT %>%
  select(bdatid,annee.x,insee.x,INSEE_COM,codification,
         pho, x.y, y.y)%>%
#Renommer les colonnes
  rename(
    annee = annee.x,
    insee = insee.x,
    pH = pho,
    x=x.y,
    y=y.y)%>%
  mutate(source="BDAT")%>%
  filter(pH!= "NULL")

BDAT_pH <- BDAT_pH %>%
  mutate(
    pH = as.numeric(pH),
    annee=as.factor(annee))

#FUSION DES DEUX BASES
base_final <- bind_rows(igcs_final_pH, BDAT_pH)
