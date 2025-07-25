dataINLA = train_aggr #calib_points 
data_test = test_data
name = "arg"
type = c("KO")
approach = c("Ponctuelle")
type_val = c("Classique")
NomsCoord = c("x", "y")
n = 500
repets = 1
k = 2
rep = 1
fold_idx= 1
sample_sizes = 500

nohup Rscript dégradation/script_principal_inla.R  arg 600 15 > test500.out&
nohup Rscript dégradation/script_principal_inla.R  arg 1000 15 > test2.out&
nohup Rscript dégradation/script_principal_inla.R  arg 2000 15 > test3.out&
nohup Rscript dégradation/script_principal_inla.R  arg 4000 15 > test4.out&
nohup Rscript dégradation/script_principal_inla.R  arg 7500 15 > test6.out&
  
  
nohup Rscript dégradation/script_principal_inla.R  pH 600 15 > test500a.out&
nohup Rscript dégradation/script_principal_inla.R  pH 1000 15 > test2a.out&
nohup Rscript dégradation/script_principal_inla.R  pH 2000 15 > test3a.out&
nohup Rscript dégradation/script_principal_inla.R  pH 4000 15 > test4a.out&
nohup Rscript dégradation/script_principal_inla.R  pH 7500 15 > test6a.out&
  
  