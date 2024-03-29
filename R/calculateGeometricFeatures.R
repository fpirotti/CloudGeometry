## usethis namespace: start
#' @useDynLib CloudGeometry, .registration = TRUE
#' @importFrom Rcpp evalCpp sourceCpp
## usethis namespace: end
NULL

#' calculateNeighboursIDX
#'
#' @param inMatrix a matrix with columns of X Y Z coordinates,
#' and rows of observations, so M rows and 3 columns
#' @param  radius maximum radius for which to find nearest neighbours
#' @param progress a logical boolean - if true a progress bar
#' will be shown
#' @return integer matrix with id values M x 100 (k=100) which is the
#' number of nearest neighbours
#' @importFrom RANN nn2
#' @export
#'
#' @examples #nf <- calculateGeometricFeatures(lidar[1:100,],5, TRUE)
calculateNeighboursIDX <- function(inMatrix, radius=1,  progress=T){
  rr<-RANN::nn2(inMatrix, k=100,  searchtype='radius',  radius = radius)
  rr$nn.idx[rr$nn.idx==0]<-NA
  rr$nn.idx
}


#' calculateGeometricFeatures
#'
#' @param inMatrix a matrix with columns of X Y Z coordinates,
#' and rows of observations, so M rows and 3 columns
#' @param  idx a matrix with nearest neighbour indices of the inMatrix cloud,
#' so M x N where N is the number of nearest neighbours (k value in RANN library)
#' @param progress a logical boolean - if true a progress bar
#' will be shown
#' @return matrix with geometric features
#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @export
#'
#' @examples
#' nn <- calculateNeighboursIDX(lidar[1:100,],5, TRUE)
#' nf <- calculateGeometricFeatures(lidar[1:100,],nn, TRUE)
calculateGeometricFeatures <- function(inMatrix, idx,  progress=T){
  if(nrow(inMatrix)!=nrow(idx)){
    stop("Number of rows in matrix with XYZ and in matrix
         with nearest neighbours ids do not match!  inMatrix=",
         nrow(inMatrix),  " vs idx=",
         nrow(idx) )
  }
  geometricFeaturesCalculate(inMatrix, idx, progress)
}
