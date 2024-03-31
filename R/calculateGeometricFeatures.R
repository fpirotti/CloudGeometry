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
#' @param  radius maximum radius for which to find nearest neighbours. If
#' the varRadius parameter
#' @param varRadius boolean (default is false) - if to choose a radius for each
#' point depending on the minimization of eigenEntropy (see Weinmann et al. 2015)
#' @param progress a logical boolean (default is true) - if true a progress bar
#' will be shown
#' @param threads number of threads. If zero (default) is set, it will use t-2,
#' i.e. all threads except 2 if machine has more than 4 threads. Only works if
#' CloudGeometry was build with OpenMP support (see documentation)
#' @return matrix with M rows and columns with the following geometric features:
#' \describe{
#'    \item{nNeighbours}{number of neighbours in 3d space of sphere of radius}
#' }
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @import RcppProgress
#' @export
#'
#' @examples #subset the first 100 rows of the lidar point cloud example to limit execution time.
#' #nn <- calcGF(lidar,5, TRUE)
#' # ## bind
#' # data.table::fwrite( cbind(lidar, nn), "out.csv")
#'
calcGF <- function(pc3d, radius=1,  varRadius=FALSE, progress=T,  threads=0){

  if(nrow(pc3d)<4){
    stop("Number of rows in matrix are too few: ",
         nrow(pc3d),  " found.")
  }

  if(ncol(pc3d)!=3){
    stop("There should be three  columns in matrix,
         with XYZ coordinates. Your matrix has ",
         ncol(pc3d), " columns.")
  }

  ne<-nnEigen(pc3d,  radius = radius, varRadius=varRadius,
              progress=progress, threads=threads)

  message("Almost done, doing final calculations...")
  if(nrow(ne)>0){
    gf <- data.frame(
      nNeighbours = ne[,6],
      eigenValue1 = ne[,1],
      eigenValue2 = ne[,2],
      eigenValue3 = ne[,3],
      eigenEntropy = ne[,5],
      eigenSum =  ne[,4],
      PCA1 =  ne[,1] / ne[,4],
      PCA2 =  ne[,2] / ne[,4],
      linearity = (ne[,1] - ne[,2]) / ne[,1],
      planarity  = (ne[,2] - ne[,3]) / ne[,1],
      sphericity = ne[,3] / ne[,1],
      verticality = ne[,7],
      anisotropy   =  (ne[,1] - ne[,3]) / ne[,1],
      omnivariance =  (ne[,1] * ne[,2] *  ne[,3])^(1/3),
      change_of_curvature   =  ne[,3] / ne[,4]
    )
  } else {
    message("No rows, maybe interrupted by user?")
    gf <- NA
  }


  gf

}




