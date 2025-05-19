# INLA -------------

inla.setOption(
  num.threads = 15 ,
  inla.mode="experimental"
)
# il faut récupérer une licence sur Pardiso : inla.pardiso()
# inla.setOption(pardiso.license ="82F70E96BE7DA6A5956D4DF8F31E127ACCB33C981DE83F430BA469A2")





# Grille de prédictions
# transformation en sp depuis terra
pxl <- as.data.frame(parcR,xy=T)
colnames(pxl)[3] <- "qrf"
gridded(pxl) <- ~x+y




# Préparation du mesh pour la méthode INLA pour la discrétisation


# max.edge = diff(range(coords[,1]))/(3*5)
# bound.outer = diff(range(range(coords[,1])))/3
# 
# bndint <- inla.nonconvex.hull( points = dataINLA, convex=-.05)
# bndext <- inla.nonconvex.hull(points = dataINLA, convex=-.2)
# 
# # Use of inla.mesh.2d 
# prmesh = inla.mesh.2d(loc = coords,
#                       boundary = list(int = limit_zone,
#                                       out = bndext),
#                       max.edge = c(.5,2)*max.edge, 
#                       cutoff = cutoffValue,
#                       crs = crs(limit_zone)
# )
# 
# 
# px_mesh <- fm_mesh_2d_inla(
#   loc = dataINLA,
#   boundary = limit_zone,
#   max.edge = c(.5,2)*max.edge,
#   crs = st_crs(limit_zone)
# )


max.edge = diff(range(coords[,1]))/(3*5)

bound.outer = diff(range(coords[,1]) ) / 3

mesh3 = inla.mesh.2d(loc=coords,
                     max.edge = c(1,2)*max.edge,
                     offset=c(max.edge, bound.outer),
                     cutoff = max.edge/10)

ggplot() +
  # geom_sf(data=limit_zone,
  #         color='turquoise',fill='transparent')+  
  gg(mesh3) +
  geom_sf(data=as(dataINLA,"sf"),
          col='purple',
          size=1.7,alpha=0.5) 


# ggplot() +
#   geom_fm(data = px_mesh) +
#   geom_sf(
#     data = counts_df[counts_df$count > 0, ],
#     aes(color = count),
#     size = 1,
#     pch = 4
#   ) +
#   theme_minimal()
# # 
# ggplot() +
#   gg(mesh3) +
#   gg( dataINLA ) +
#   gg( dataINLA1,col=2 ) +
#   coord_equal()
# 
# # 
# mesh <- fm_mesh_2d_inla(loc=coords,boundary = limit_zone, max.edge = 50)
# ggplot() +
#   geom_fm(data = mesh)
# 

# Modelling inla
matern <-
  INLA::inla.spde2.pcmatern(mesh3,
                            alpha = 2,
                            prior.sigma = c(10, 0.01),# P(sigma > 1) = 0.5
                            prior.range = c(500, 0.01)  # P(range < 100000 m) = 0.9
  )

cmp <- activ ~ Intercept(1) +
  field(coordinates, model = matern)


fitKO <- inlabru::bru(components = cmp,
                      data = dataINLA,
                      family = "Gaussian",
                      domain = list(coordinates = mesh3),
                      options = list(
                        control.inla = list(int.strategy = "eb"),
                        verbose = FALSE)
)


summary(fitKO)


fi2plot = fitKO
int.plot <- plot(fi2plot, "Intercept")
spde.range <- spde.posterior(fi2plot, "field", what = "range")
spde.logvar <- spde.posterior(fi2plot, "field", what = "log.variance")
range.plot <- plot(spde.range)
var.plot <- plot(spde.logvar)

multiplot(range.plot, var.plot, int.plot)



predKO <- predict(
  fitKO,
  n.samples = nsim,
  pxl,
  ~Intercept + field   ,
  num.threads = 10
)



p = tm_shape(rast(predKO) ) +
  tm_raster(c("mean"),
            col.scale = tm_scale(style = "quantile",n = 10))


print(p)

writeRaster(rast(predKO)[["mean"]],file="output/predKOINLA.tif", overwrite = T)

# Validation croisée-------------
print("Validation croisée----------------")


datacov$predINLAKO = NA

resuXval <- 
  foreach(i = 1:k,
          .errorhandling='pass' ) %do% {
            
            print(i)
            
            # set to na to run a cross valid with inla
            # Mettre en NA les individus pour la validation crois?e
            
            dataINLA$elt <- dataINLA$activ
            dataINLA$elt[ fold[[i]] ]  <- NA
            
            
            
           cmp <- elt ~ Intercept(1) +
              field(coordinates,
                    model = matern)
            
            Myfit <- bru(cmp,
                         data = dataINLA,
                         family = "Gaussian",
                         options = list(
                           control.inla = list(int.strategy = "eb"),
                           verbose = FALSE)
            )
            
            # find the prediction in the output....
            fitted <- Myfit$summary.fitted.values$mean[1:length(dataINLA$activ)] 
            
            mask <- is.na(dataINLA$elt)
            
            datacov$predINLAKO[ fold[[i]] ] <-  fitted[mask] 
            
            fitted[mask] 
          }


resuXvalTKO <-  Myeval(datacov$predINLAKO,   datacov[,name] )






