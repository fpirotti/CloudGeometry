## usethis namespace: start
#' @useDynLib CloudGeometry, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
## usethis namespace: end
NULL

#' calculateGeometricFeatures
#'
#' @param inMatrix a matrix with columns of X Y Z coordinates,
#' and rows of observations
#' @param  radius maximum radius for which to find nearest neighbours
#' @param progress a logical boolean - if true a progress bar
#' will be shown
#' @return matrix with geometric features
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @importFrom RANN nn2
#' @export
#'
#' @examples calculateGeometricFeatures(lidar,5, TRUE)
calculateGeometricFeatures <- function(inMatrix, radius=1,  progress=T){
  rr<-RANN::nn2(inMatrix, k=100,  searchtype='radius',  radius = radius)
  rr$nn.idx[rr$nn.idx==0]<-NA
  geometricFeaturesCalculate(inMatrix, rr$nn.idx, progress)
}

