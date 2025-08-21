dataINLA = train_aggr #calib_points 
data_test = test_data
name = "pH"
type = c("KO")
approach = c("Ponctuelle")
type_val = c("Classique")
NomsCoord = c("x", "y")
n = 4000
repets = 1
k = 30
rep = 1
fold_idx= 1
sample_sizes = c(600,1000)
kmax <- 30
ntree <- 350
NomsCoord <- c("x", "y")
types_validation <- c("Classique", "Spatiale")
drive = "/media/communs_infosol/" # ou "Y:/"
DistanceGeomasking = 0



name <- "pH"  #arg
sample_sizes <- 600 # c(500,1000 ) # c(600,800,1000,1200,1300,1400,1600,1800,2000,3000,4000,5000,6000,7000,7600)
repets <- 5

n=4000

n=400



nohup Rscript dégradation/script_principal_inla.R  arg 400 5 > test400.out&
nohup Rscript dégradation/script_principal_inla.R  arg 600 5 > test600.out&
nohup Rscript dégradation/script_principal_inla.R  arg 1000 5 > test2.out&
nohup Rscript dégradation/script_principal_inla.R  arg 2000 5 > test3.out&
nohup Rscript dégradation/script_principal_inla.R  arg 4000 5 > test4.out&
nohup Rscript dégradation/script_principal_inla.R  arg 5000 5 > test5000.out&
nohup Rscript dégradation/script_principal_inla.R  arg 7500 5 > test6.out&
  
  
nohup Rscript dégradation/script_principal_inla.R  pH 400 5 > test50a.out&
nohup Rscript dégradation/script_principal_inla.R  pH 600 5 > test500a.out&
nohup Rscript dégradation/script_principal_inla.R  pH 1000 5 > test2a.out&
  nohup Rscript dégradation/script_principal_inla.R  pH 2000 5 > test3a.out&
  nohup Rscript dégradation/script_principal_inla.R  pH 3000 5 > test3b.out&
  nohup Rscript dégradation/script_principal_inla.R  pH 4000 5 > test4a.out&
  nohup Rscript dégradation/script_principal_inla.R  pH 7500 5 > test6a.out&
  nohup Rscript dégradation/script_principal_inla.R  pH 10000 5 > test6a.out&
  
  