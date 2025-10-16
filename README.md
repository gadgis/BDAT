

# Quantification de l'apport du géoréférencement précis de données 

Ce repertoire est créeé pour héberger les scripts utilisés dans le cadre d'un mémoire de fin d'études sur le thème portant : Apport de la géolocalisation précise de données fournies par les agriculteurs pour l’analyse spatiale des propriétés des sols à l’échelle d’un département. 
Les données utilisées dans le cadre de ce stage sont celles de la BDAT (Base de Données Analyse de Terres) et celles de l'IGCS (Inventaire Gestion et Conservation des Sols). La BDAT est un programme du Gis Sol (Groupement d'intérêt scientifique sur les sols) qui recueille pour la France métropolitaine les résultats d’analyses de terre effectuées à la demande des agriculteurs auprès de laboratoires d’analyse de terre agréés par le ministère de l’agriculture. Les résultats récupérés sont le plus souvent géolocalisés à la commune d’origine du prélèvement et la restitution des statistiques sur les propriétés des sols se fait à l’échelle du canton ou de la petite région agricole

## Présentation du stage
Le stage vise à quantifier *l'apport du géoréférencement précis des données dans la spatialisation des propriétés des sols comparé à un géoréférencement à la commune comme dans la BDAT*.
En effet , le géoréférecemment comme la restitution des données sur les propriétés des sols de la BDAT est tributaire du support géographique qui n'est pas fixe car les communes peuvent se fusionnées comme se disloquées, ce qui contitue une limite dans l'explotation et la restitution des résultats. Ainsi pour monter l'importance d'avoir des données precises indépendantes du support géographique, ce présent travail se base sur l'hypothèse que les données géoréférencées apportent une précision plus grande que les données agrégées dans la cartographie des propriétés des sols. 
Pour répondre à cette hypothèse, des techniques de cartographie des sols par modélisation statistique (CSMS)sont utilisées dans des scénarios différents. 
Le travail s'articule autour de deux approches de CSMS:
 - La première approche est appelée « Données ponctuelles » où nous utilisons les données ponctuelles comme données d’entrée dans les modèles;
 - La deuxième consiste à agréger à la commune, par la moyenne, les valeurs des propriétés de sol et celles des covariables et de les affecter à leurs centroïdes. Cette méthode nous permet de simuler la BDAT en masquant la localisation de l’information spatiale. Les données ainsi obtenues et qui vont servir d’entrée dans les méthodes statistiques sont des moyennes communales comme dans la diffusion des données de la BDAT. Nous appelons cette méthode « Désagrégation » par le fait que la valeur au support communal est réaffectée au centroïde de la commune.
L'image suivante montre les des différents scénarios de spatialisation testés
<p align="center">
  <img src="carto/modeles_CSMS.png" width="500" alt="Modèles CSMS">
</p>
Dans ces différentes approches, les méthodes utilisées pour la cartographie sont : une régression suivant l’approche SCORPAN avec la forêt aléatoire (Random Forest), une méthode géostatistique (le krigeage) qui se fonde sur les corrélations spatiales et enfin une méthode hybride (le krigeage avec dérive externe) combinant les deux précédentes méthodes. Les méthodes géostatistique et hybride sont implémentées dans un contexte bayésien en utilisant l'approximation de Laplace imbriquée intégrée avec l'approche d'équation aux dérivées partielles stochastiques (INLA SPDE).
Une procédure de validation croisée (classique et spatiale) est également implémenté pour évaluer la capacité de généralisation des modèles prédictifs en calculant des indicateurs de qualité de prédiction comme le NSE(Nash-Sutcliffe Efficiency), l'EM (Erreur moyenne) le REQM(racine carrée de l’erreur quadratique moyenne ), le CCC (Coefficient de concordance).
Au cours de cette procédure de validation croisée, une expérience de dégradation progressive de l'information spatiale est menée pour évaluer la robustesse des modèles face aux variations de la densité d'échantillonnage. Cette expérience consiste à faire évoluer le jeu de calibration de façon progressive en faisant des sous-échantillonnages aléatoires, à ajuster un nouveau modèle à partir de ce jeu de données dégradé et à calculer les indicateurs de qualité de prédiction.
La figure ci-après montre la procédure de validation
<p align="center">
  <img src="dégradation/schema_degradation.png" width="500" alt="Schéma de dégradation">
</p>

Les blocs suivant décrivent les étapes scuccésives ayant perlis d'effectuer ce travail de l'extraction des données sols jusqu'à la visualisation des résultats.

## Model workflow (R scripts)

Cette partie....

### 1. Préparation des données sols




### 2. Extraction des covariables


### 3. Calibration des modèles et cartographie




### 4. Expérience de la dégradation de l'information spatiale




### 5. Expérience de Geomasking




## Réferences

