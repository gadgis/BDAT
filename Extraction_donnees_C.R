
# Chargement de package ----
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(mpspline2)

#Importation des couches ----


igcs<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_agri.rds")
bdat<-readRDS("Y:/BDAT/traitement_donnees/MameGadiaga/resultats/BDAT_53.rds")

###Harmonisation des données sur le C----

igcs_agri_C<-igcs%>%
  select(
    id_profil,
    top,
    bottom,
    no_hrzn,
    annee,
    INSEE_COM,
    x,
    y,
    C,
    mat_org
  )%>% 
  mutate(C = coalesce(C, mat_org / 1.72))%>% 
  filter(!is.na(C))
# Harmonisation des horizons----


####selection des profils à horizons uniques
pf_uniq <- igcs_agri_C %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() == 1) %>%          
  ungroup()

####selection des profiles à deux horizons dont la limite inférieur du premier
#horizon est comprise entre 30 et 35 cm
P2hrzn1_30_35 <- igcs_agri_C %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() == 2) %>%          
  ungroup() %>%
  filter(no_hrzn==1, bottom >= 30 & bottom <= 35)


####selection des profils à deux horizons dont le premier horizon ne fait pas 30 cm 

pf_mp <- igcs_agri_C %>% 
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
  group_by(id_profil,annee,INSEE_COM,x,y) %>% 
  summarise(
    C = round(weighted.mean(C, thick, na.rm = TRUE),1)   # moyenne pondérée
  ) %>% 
  ungroup()

####selection des profils dont les horizons sont supérieure ou égales à 3 et que le premier
#horizon est compris entre 30 et 35 cm
p3hrzn1_30_35 <- igcs_agri_C %>% 
  filter(top < bottom) %>% 
  arrange(id_profil, top) %>% 
  group_by(id_profil) %>% 
  filter(n() >= 3) %>%          
  ungroup() %>%
  filter(no_hrzn==1, bottom >= 30 & bottom <= 35)


####selection des profils à slpiner
spl_dfs <- igcs_agri_C %>% 
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
  var_name = "C",
  d = c(0, 30),
  lam = 0.1  
)



# Étape 1: Créer un dataframe à partir de la liste
igcs_spline <- imap_dfr(spl_layers, ~ {
  
  tibble(
    # Identifiant et métadonnées de base
    id_profil = .x$id_profil,
    
    # Horizon standard (est_dcm)
    C_000_030_cm = round(.x$est_dcm[["000_030_cm"]],1),
    
    # Erreurs (est_err)
    RMSE = .x$est_err[["RMSE"]],
    RMSE_IQR = .x$est_err[["RMSE_IQR"]])
  
})

# Étape 2: Fusionner avec les métadonnées supplémentaires (

igcs_spline <- igcs_spline %>%
  left_join(
    spl_dfs %>% distinct(id_profil, INSEE_COM, annee, x, y),
    by = "id_profil"
  ) %>%
  rename(C = C_000_030_cm) %>%
  select(id_profil, annee, INSEE_COM, x, y, C)

##haromisation des colonnes
pf_uniq <- pf_uniq %>%
  select(id_profil, annee,INSEE_COM, x, y, C)
p_mp <- p_mp %>%
  select(id_profil, annee,INSEE_COM, x, y, C)

p3hrzn1_30_35 <- p3hrzn1_30_35 %>%
  select(id_profil, annee,INSEE_COM, x, y, C)

P2hrzn1_30_35 <- P2hrzn1_30_35 %>%
  select(id_profil, annee,INSEE_COM, x, y, C)

igcs_final_C<-rbind(p_mp, pf_uniq, p3hrzn1_30_35, igcs_spline, P2hrzn1_30_35)

igcs_final_C <- igcs_final_C %>%
  mutate(source="IGCS",
         INSEE_COM=as.numeric(INSEE_COM),
  )
#Création de la base finale----

#BDAT

#selection des colonnes d'intéret BDAT_C

BDAT_C<-bdat %>%
  select(bdatid,annee.x,insee.x,INSEE_COM,codification,
         corgox, x.y, y.y)%>%
  #Renommer les colonnes
  rename(
    annee = annee.x,
    insee = insee.x,
    C = corgox,
    x=x.y,
    y=y.y)%>%
  mutate(source="BDAT")%>%
  filter(C!= "NULL")

BDAT_C <- BDAT_C %>%
  mutate(
    C = as.numeric(C),
    annee=as.factor(annee))



#FUSION DES DEUX BASES

base_final <- bind_rows(igcs_final_C, BDAT_C)
base_final <- base_final %>%
  filter(!is.na(C))


saveRDS(base_final, "Y:/BDAT/traitement_donnees/MameGadiaga/resultats/igcs_bdat_C.rds")
