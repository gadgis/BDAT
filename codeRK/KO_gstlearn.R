# krigeage ordinaire ----------

## Creation de la base de données -------
mydb <- Db_create()

for (var in c(NomsCoord,name) ) {
  mydb[var] <- datacov[,var]
}

# defining the coordinates and the active variable
err = mydb$setLocators(NomsCoord, 
                       locatorType = ELoc_X(), 
                       cleanSameLocator = TRUE)


err = mydb$setLocators(c(name),
                       locatorType = ELoc_Z(),
                       cleanSameLocator = TRUE
)


## Calcul la tranformation gaussienne --------

anam = AnamHermite(30)
err = anam$fitFromLocator(mydb)
err = anam$rawToGaussian(mydb, name)
anam


p1 = ggDefault() +
  plot.hist(mydb, name=name,  bins = 25, bg = "orange", color= "gray") +
  plot.decoration(xlab = "Raw variable", title="")
yy = sort(mydb[name])
zz = sort(mydb[paste0("Y.",name) ] )
p2 = ggDefault() +
  plot.XY(yy, zz, flagLine=TRUE, flagPoint=FALSE, color = "red") +
  plot.XY(yy, anam$transformToRawVec(yy), flagLine=TRUE, flagPoint=FALSE, color = "blue") +
  plot.decoration(xlab="Gaussian", ylab="Raw")
p3 = ggDefault() +
  plot.hist(mydb, name=paste0("Y.",name),  bins = 25, bg = "skyblue", color= "gray") +
  plot.decoration(xlab = "Gaussian variable", title="")

ggarrange(p1, p2, p3, nrow = 1, ncol = 3)



## Calculer le varigramme ------------

# myVarioParamOmni = VarioParam()
# mydir = DirParam_create(npas=30,
#                         dpas=50
# )
# myVarioParamOmni$addDir(mydir)


varioParamOmni = VarioParam_createOmniDirection(npas=20,
                                                dpas=75,
                                                toldis = 0.5
                                                ) 
varioexp = Vario(varioParamOmni)
err = varioexp$compute(mydb)


fitmod = Model()


types = ECov_fromKeys(c("NUGGET","SPHERICAL"))
err = fitmod$fit(varioexp, types=types)

# ggplot() + plot.varmod(varioexp, fitmod)
# 
# myvario = Vario(myVarioParamOmni)
# err = myvario$compute(mydb,ECalcVario_VARIOGRAM())
# mymodelM = Model_createFromDb(mydb)
# 
# err = mymodelM$fit(myvario,ECov_fromKeys(c("NUGGET" , "EXPONENTIAL" )))
# 


v <- ggplot() +
  plot.varmod(varioexp, fitmod)+
  plot.decoration(title="Variogramme du Krigeage Ordinaire")

print(v)

## krigeage ----------

# fr_gst <- sf_to_gstlearn(parc)
# 
# # creating a grid covering the metropolitan boundaries
# MyGrid = DbGrid_createCoveringDb(mydb,
#                                  dx = rep(5, 2)
#                                  )
# 
# 
# # selection of the cells inside the polygons
# err  = db_polygon(db = MyGrid,
#                   polygon = fr_gst)


MyGrid <- terra_to_gstlearn( parcR )

neighU = NeighUnique_create()

err = kriging(mydb, MyGrid, fitmod, neighU)


selectivity = Selectivity_createByKeys(c("Z"),
                                       flag_est=TRUE,
                                       flag_std=TRUE)
## Faire la transformation  inverse

err = ConditionalExpectation(MyGrid,
                             anam,
                             selectivity,
                             "K*.estim",
                             "K*.stdev",
                             nbsimu=50,
                             namconv=NamingConvention("CE",FALSE,TRUE,FALSE)
                             )


p = ggDefaultGeographic()
p = p + plot(MyGrid,"CE*estim")
p = p + plot(mydb,size = .1)
p = p + plot.decoration(title = "Conditional Expectation")
print(p)


## Validation croisée -----------


err = mydb$deleteColumns(c("CV.SK.Y.*", "CV.CE.Z.*"))
tab <- sample(1:k, size = mydb$getSampleNumber(), replace = TRUE)
mydb["code"] <- tab
err  = mydb$setLocators("code", locatorType = ELoc_C(), cleanSameLocator = TRUE)
# defining the variable to be estimated: Z
selectivity     = Selectivity_createByKeys(c("Z"), 
                                           flag_est=TRUE, flag_std=TRUE)

 neigh = NeighMoving_create(flag_xvalid = FALSE, nmaxi = 50, radius = 2500)

# simple kriging of the Gaussian variable
err = xvalid(mydb, 
             model = fitmod, 
             neigh = neigh, 
             flag_kfold = TRUE,
             flag_xvalid_est = -1,
             flag_xvalid_std = -1,
             namconv = NamingConvention("CV.SK")
)

err = ConditionalExpectation(db = mydb, 
                             anam = anam, 
                             selectivity = selectivity, 
                             name_est = paste("CV.SK", "Y", name, "estim", sep = "."),
                             name_std = paste("CV.SK", "Y", name,  "stdev", sep = "."),
                             namconv=NamingConvention("CE",FALSE,TRUE,FALSE))

err = mydb$setName("CE.Z-estim", paste("CV.CE", "Z", name, "estim", sep = "."))
err = mydb$setName("CE.Z-stdev", paste("CV.CE", "Z", name, "stdev", sep = "."))


resuXvalY <-  Myeval(mydb["CV.CE.Z.activ.estim"],   
                     mydb[name])


datacovPT$predKO <- mydb["CV.SK.Y.activ.estim"]


resuXvalKOgstlearn <-  Myeval(mydb["CV.SK.Y.activ.estim"],   
                    mydb["Y.activ"])





## Export vers terra ------

tt = MyGrid[]
x2 <- gstlearn_to_terra(MyGrid)
crs(x2) <- crs(uts)

x2 <- ifel(is.na(x2[[1]] ),NA, x2)

writeRaster(x2[[4]],file="output/predKOgstlearn.tif", overwrite = T)



p = tm_shape(x2) +
  tm_raster("CE.Z-estim") +
  tm_shape(uts) + tm_borders() +
  tm_shape(solsP) + tm_dots() +
tm_layout(legend.outside = T)


print(p)
