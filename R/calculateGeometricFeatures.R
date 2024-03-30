## usethis namespace: start
#' @useDynLib CloudGeometry, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
#' @export
## usethis namespace: end
NULL


#' calcGF
#' @description Calculate Geometric Features
#'
#' @param inMatrix a matrix with columns of X Y Z coordinates,
#' and rows of observations, so M rows and 3 columns
#' @param  radius maximum radius for which to find nearest neighbours
#' @param progress a logical bolean - if true a progress bar
#' will be shown
#' @param threads number of threads. If zero (default) is set, it will use t-2,
#' i.e. all threads except 2 if machine has more than 4 threads. Only works if
#' CloudGeometry was build with OpenMP support (see documentation)
#' @return matrix with M rows and 3 columns with the geometric
#' features (3 eigen values in ascending order)
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @import RcppProgress
#' @export
#'
#' @examples #subset the first 100 rows of the lidar point cloud example to limit execution time.
#' nn <- calcGF(lidar[1:100,],5, TRUE)
calcGF <- function(inMatrix, radius=1,  progress=T, threads=0){
  if(nrow(inMatrix)<4){
    stop("Number of rows in matrix are too few: ",
         nrow(inMatrix),  " found.")
  }
  if(ncol(inMatrix)!=3){
    stop("There should be three  columns in matrix, with XYZ coordinates. Your matrix has ", ncol(inMatrix), " columns.")
  }
  message("Starting eigan with radius=", radius, " and ", threads, " threads.")
  nnEigen(inMatrix,  radius = radius,  progbar=progress, threads=threads)


}


