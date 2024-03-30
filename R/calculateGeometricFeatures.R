## usethis namespace: start
#' @useDynLib CloudGeometry, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @export
## usethis namespace: end
NULL


#' calcEigen
#' @description Calculate Eigen values
#'
#' @param pc3d matrix 3D point cloud in  matrix format with columns
#' with X Y Z coordinates, #' and rows of observations, so M rows
#' and 3 columns
#' @param  radius maximum radius for which to find nearest neighbours
#' @param progress a logical boolean - if true a progress bar
#' will be shown
#' @param threads number of threads. If zero (default) is set, it will use t-2,
#' i.e. all threads except 2 if machine has more than 4 threads. Only works if
#' CloudGeometry was build with OpenMP support (see documentation)
#' @return matrix with M rows and 3 columns with the eigenvalues of the 3D tensor
#' of the 3d coordinates of the nearest neighbours around each point.
#' (3 eigen values in descending order). If fewer than 3 neighbours are found,
#' NaN is returned for all 3 eigenvalues at that row.
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @import RcppProgress
#' @export
#'
#' @examples #subset the first 100 rows of the lidar point cloud example to limit execution time.
#' nn <- calcEigen(lidar[1:100,],5, TRUE)
calcEigen <- function(pc3d, radius=1,  progress=T, threads=0){
  if(nrow(pc3d)<4){
    stop("Number of rows in matrix are too few: ",
         nrow(pc3d),  " found.")
  }
  if(ncol(pc3d)!=3){
    stop("There should be three  columns in matrix, with XYZ coordinates.
         Your matrix has ", ncol(pc3d), " columns.")
  }

  pc3d[,1] <- pc3d[,1] - min(pc3d[,1])
  pc3d[,2] <- pc3d[,2] - min(pc3d[,2])
  nnEigen(pc3d,  radius = radius,  progbar=progress, threads=threads)

}

#' calcGF
#' @description Calculate Geometric Features
#'
#' @param pc3d matrix 3D point cloud in  matrix format with columns
#' with X Y Z coordinates, #' and rows of observations, so M rows
#' and 3 columns
#' @param  radius maximum radius for which to find nearest neighbours
#' @param progress a logical bolean - if true a progress bar
#' will be shown
#' @param threads number of threads. If zero (default) is set, it will use t-2,
#' i.e. all threads except 2 if machine has more than 4 threads. Only works if
#' CloudGeometry was build with OpenMP support (see documentation)
#' @return matrix with M rows and columns
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @import RcppProgress
#' @export
#'
#' @examples #subset the first 100 rows of the lidar point cloud example to limit execution time.
#' nn <- calcGF(lidar[1:100,],5, TRUE)
#'
calcGF <- function(pc3d, radius=1,  progress=T, threads=0){
  if(nrow(pc3d)<4){
    stop("Number of rows in matrix are too few: ",
         nrow(pc3d),  " found.")
  }
  if(ncol(pc3d)!=3){
    stop("There should be three  columns in matrix,
         with XYZ coordinates. Your matrix has ",
         ncol(pc3d), " columns.")
  }

  ne<-nnEigen(pc3d,  radius = radius,  progbar=progress, threads=threads)

  eigenSum = rowSums(ne)

  gf <- data.frame(
    eigenValue1 = ne[1,],
    eigenValue2 = ne[2,],
    eigenValue3 = ne[3,],
    eigenSum = eigenSum,
    linearity = (ne[1,] - ne[2,]) / ne[1,],
    planarity  = (ne[2,] - ne[3,]) / ne[1,],
    sphericity = ne[3,] / ne[1,],
    omnivariance =  (ne[1,] * ne[2,] *  ne[3,])^(1/3),
    anisotropy   =  (ne[1,] - ne[3,]) / ne[1,],
    change_of_curvature   =  ne[3,] / eigenSum
  )

  gf

}




