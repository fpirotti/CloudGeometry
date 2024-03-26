## usethis namespace: start
#' @useDynLib CloudGeometry, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
## usethis namespace: end
NULL

#' calculateGeometricFeatures
#'
#' @param inMatrix a matrix with columns of X Y Z coordinates,
#' and rows of observations
#'
#' @param progress a logical boolean - if true a progress bar
#' will be shown
#' @return matrix with geometric features
#' @importFrom Rcpp evalCpp
#' @export
#'
#' @examples calculateGeometricFeatures(lidar, TRUE)
calculateGeometricFeatures <- function(inMatrix,  progress=T){
  geometricFeaturesCalculate(inMatrix, progress)
}

