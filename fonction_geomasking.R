#' r_calcul_quanti
#'
#' @param dataset (dataframe) tableau avec les coordonnées à flouter
#' @param DistFloutage (string) unité cartographique sur laquelle calculer les statistiques : région, département, petite région agricole (pra) ou canton
#' @param NomsCoord Vecteur avec les noms des colonnes contenant les coordonnées x et y
#'
#' @return un tableau avec les coordonnées floutées et le reste des collonne
#' @export NA crée une table avec les coodonnées floutées en remplaçant les valeurs des colonnes NomsCoord. La fonction travaille en m

geomasking <- function(
    
  dataset,
  DistGeom
){

  
  deplacer_aleatoire <- function(x, y, distance) {
    angle_rad <- runif(1, 0, 2 * pi)  # Angle aléatoire entre 0 et 2π
    pv <- runif( length(x) , .8, 1.2)  # une petite variation de distance
    nouveau_x <- x + pv*distance * cos(angle_rad)
    nouveau_y <- y + distance * sin(angle_rad)
    
    return(list(x = nouveau_x, y = nouveau_y))
  }
  
  
  res  <- deplacer_aleatoire(dataset[,NomsCoord][,1],
                     dataset[,NomsCoord][,2],
                     DistGeom)
  
  dataset[,NomsCoord][,1] <- res$x
  dataset[,NomsCoord][,2] <- res$y
  
  return(dataset)
  # 
  # x1 =dataset[,NomsCoord][,1]
  # y1 =dataset[,NomsCoord][,2]
  # 
  # 
  # x2 =res$x
  # y2 =res$y
  # 
  # summary(sqrt((x1-x2)^2 + (y1-y2)^2  ))
}

