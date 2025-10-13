

# Quantification de l'apport du géoréférencement précis de données 

Ce repertoire est créeé pour héberger les scripts utilisés dans le cadre d'un mémoire de fin d'études sur le thème portant : Apport de la géolocalisation précise de données fournies par les agriculteurs pour l’analyse spatiale des propriétés des sols à l’échelle d’un département. 
Les données utilisées dans le cadre de ce stage sont celles de la BDAT (Base de Données Analyse de Terres) et celles de l'IGCS (Inventaire Gestion et Conservation des Sols). La BDAT est un programme du Gis Sol (Groupement d'intérêt scientifique sur les sols) qui recueille pour la France métropolitaine les résultats d’analyses de terre effectuées à la demande des agriculteurs auprès de laboratoires d’analyse de terre agréés par le ministère de l’agriculture. Les résultats récupérés sont le plus souvent géolocalisés à la commune d’origine du prélèvement et la restitution des statistiques sur les propriétés des sols se fait à l’échelle du canton ou de la petite région agricole

## Présentation du stage
Le stage...


## Model workflow (R scripts)

Cette partie....

### 1. Préparation des données sols

-   [11_soil_PFB_prep.R](11_soil_PFB_prep.R) - Connect to BIS 
-   [12_soil_BIS_master_tibble.R](12_soil_BIS_master_tibble.R) - Combine all BIS soil information 


### 2. Extraction des covariables

-   [20_cov_prep_gdal.R](20_cov_prep_gdal.R) - Assemble & prepare predictors (covariates as raster data):
    -   Designate coordinate system (projection)
    -   Resample covariates so they have the same origin, cell locations and extent
    -   Mask nodata areas (water and areas outside NL) of continuous covariates

### 3. Calibration des modèles et cartographie <img src="carto/modeles_CSMS.png" align="right" width="500"/>

-   [30_regression_matrix_BIS.R](30_regression_matrix_BIS.R) - Read in prepared covariate stack and soil point data with coordinates; overlay rasters and points and prepare regression matrix (extract covariate values from sampled/observed locations)

### 4. Expérience de dégradation de l'information spatiale <img src="dégradation/schema_degradation.png" align="right" width="500"/>

-   [40_train_RF_hyperparameter_tuning.R](40_train_RF_hyperparameter_tuning.R) - Tune random forest (RF) hyperparameters:
    -   Read in BIS regression matrix and select target (response) variable
    

### 5. Expérience de Geomasking




## Réferences

