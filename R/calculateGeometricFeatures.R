## usethis namespace: start
#' @useDynLib CloudGeometry, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @export
## usethis namespace: end
NULL

#' calcGF
#' @description Calculate Geometric Features
#'
#' @param pc3d matrix 3D point cloud in  matrix format with columns
#' with X Y Z coordinates, #' and rows of observations, so M rows
#' and 3 columns
#' @param  rk  **(1) in case of varRadius==FALSE**, rk is the radius of the sphere
#'  where  to find the nearest neighbours for each point.
#'
#'  **(2) in case of varRadius==TRUE** rk is the largest number of nearest
#'  neighbours to read from and that define the largest radius - from here 30
#'  smaller radii will be tested to find the one with the minimum eigen entropy
#'  value.
#' *the varRadius parameter is TRUE*, this is the k nearest neighbours that will
#' be evaluated for the largest starting radius.
#' @param varRadius boolean (default is false) - will use a different heuristic
#' and will add k-neighbours untile it minimizes the eigenEntropy and fixes the
#' radius.
#' @param progress a logical boolean (default is TRUE) - if true a progress bar
#' will be shown
#' @param verbose a logical boolean (default is FALSE) - if true a
#' lot of messages will be shown.
#' @param threads number of threads. If zero (default) is set, it will use t-2,
#' i.e. all threads except 2 if machine has more than 4 threads. Only works if
#' CloudGeometry was build with OpenMP support (see documentation)
#' @return matrix with M rows.
#'
#' Columns will have  the following geometric features:
#' \describe{
#'    \item{pointDensity}{Number of neighbours in 3d space of sphere of radius}
#'    \item{eigenValue1}{eigenValue 1 (largest value)}
#'    \item{eigenValue2}{eigenValue 2 (medium value)}
#'    \item{eigenValue3}{eigenValue 3 (smallest value)}
#'    \item{eigenEntropy}{Eigen Entropy: \eqn{- \sum_{i}^{3}\lambda_{i}*ln(\lambda_{i})}
#'    where eigen values are normalized as   \eqn{\lambda_{i}=\frac{e_{i}}{\sum_{i}^{3}e_{i}} } }
#'    \item{eigenSum}{Sum of Eigen values}
#'    \item{PCA1}{First principal component}
#'    \item{PCA2}{Second principal component}
#'    \item{linearity}{Linearity \eqn{\frac{\lambda_{1} - \lambda_{2} }{\lambda_{1}} }  }
#'    \item{planarity}{Planarity \eqn{\frac{\lambda_{2} - \lambda_{3} }{\lambda_{1}} }  }
#'    \item{sphericity}{Sphericity  \eqn{\frac{\lambda_{3}}{\lambda_{1}} } }
#'    \item{verticality}{ Angle of normal vector of best fit plane from the vertical normal vector. }
#'    \item{anisotropy}{  Anisotropy \eqn{\frac{\lambda_{1} - \lambda_{3} }{\lambda_{1}} }  }
#'    \item{omnivariance}{ Omnivariance \eqn{\sqrt[3]{\lambda_1\lambda_2\lambda_3 }}. }
#'    \item{change_of_curvature}{ Also called surface variation:
#'                          \eqn{\frac{e_{3}}{\sum_{i}^{3}e_{i}} }  . }
#'    \item{optRadius}{ Only available if "**varRadius=TRUE**" and shows the optimal radius value that minimizes the entropy of eigen values. }
#'
#'
#'
#' }
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @import RcppProgress
#' @export
#'
#' @examples #subset the first 100 rows of the lidar point cloud example to limit execution time.
#' #nn <- calcGF(lidar,5, FALSE, TRUE, TRUE)
#' # ## bind to original table
#' # data.table::fwrite( cbind(lidar, nn), "out.csv")
#'
calcGF <- function(pc3d, rk=1,  varRadius=FALSE, progress=T, verbose=F, threads=0){

  if(nrow(pc3d)<4){
    stop("Number of rows/points in matrix are too few: ",
         nrow(pc3d),  " points found.")
  }

  if(nrow(pc3d)<100 && varRadius){
    warning("Number of rows/points in matrix are too few for variable radius as at least 100 k-nearest neighbours are assessed as larger radius: ",
         nrow(pc3d),  " points found.")
    return(NA)
  }

  if( rk < 10 && varRadius){
    warning("rk < 10 does not make much sense with varRadius: remember that rk will
         mean k-nearest neighbours are used to define the larger radius:", rk ," points defined by user. Around 100 points are suggested.")
    return(NA)
  }

  if(ncol(pc3d)!=3){
    warning("There should be three  columns in matrix,
         with XYZ coordinates. Your matrix has ",
         ncol(pc3d), " columns.")
    return(NA)
  }

  if(verbose) message("Starting eigen calculations...")
  ne<-nnEigen(pc3d,  radius = rk, varRadius=varRadius,
              progress=progress, verbose=verbose, threads=threads)

  if(verbose) message("Almost done, doing final calculations...")

  if(nrow(ne)>0){
    gf <- data.frame(
      pointDensity = ne[,6],
      eigenValue1 = ne[,1],
      eigenValue2 = ne[,2],
      eigenValue3 = ne[,3],
      eigenEntropy = ne[,5],
      eigenSum =  ne[,4],
      PCA1 =  ne[,1] / ne[,4],
      PCA2 =  ne[,2] / ne[,4],
      linearity   = (ne[,1] - ne[,2]) / ne[,1],
      planarity   = (ne[,2] - ne[,3]) / ne[,1],
      sphericity  = ne[,3] / ne[,1],
      verticality = ne[,7],
      anisotropy   =  (ne[,1] - ne[,3]) / ne[,1],
      omnivariance =  (ne[,1] * ne[,2] *  ne[,3])^(1/3),
      change_of_curvature   =  ne[,3] / ne[,4]
    )
  } else {
    if(verbose)  message("No rows, maybe interrupted by user?")
    gf <- NA
  }

  if(varRadius){
    gf$optRadius <- ne[,8]
    gf$mcounts <- ne[,9]
    gf$fractionizer <- ne[,10]
  }

  gf

}
#
# h <- function(){
#   library(ggplot2)
#   a<-read.csv("A.csv")
#   a<-na.omit(a)
#   names(a)<-c("radius", "eigenEntropy", "pointID", "count")
#
#   a<-a[a$pointID!=0,]
#   a$pointID <- as.factor(a$pointID)
#
#   ggplot(a[as.integer(a$pointID)<10,], aes(y=eigenEntropy, x=radius,
#                                           group=pointID, color=pointID ) ) +
#     geom_line() +
#     xlab("Radius (m)") + theme_bw()
#
# }


