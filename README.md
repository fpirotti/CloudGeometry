# CloudGeometry

<!-- badges: start -->

[![R-CMD-check](https://github.com/fpirotti/CloudGeometry/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fpirotti/CloudGeometry/actions/workflows/R-CMD-check.yaml) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

The goal of CloudGeometry is fast and parallel computation of geometric features for each point in a point cloud using an auto-adaptive or fixed radius.

Geometric features (GF) catch geometrical patterns around a certain position in 3D space based on information available from the neighboring points (local neighborhood). A feature-rich set of observations support training AI (machine learning and deep learning) models for segmentation and classification.

Geometric features can be characterized using the eigenvalues of the 3D structure tensor represented by the 3D covariance matrix derived from the 3D coordinates of all points within the local neighborhood ([Weinmann et al. 2014](#2), [Weinmann et al. 2015](#1)).

[![](images/clipboard-2941799751.png)](https://arxiv.org/abs/1107.0550)

Eigenvalues partitions for the local neighborhood PCA. Source: [Brodu&Lague 2012](#3)

Due to the the expected large volume of points in point clouds that are surveyed with laser scanners or created by photogrammetric procedures, this package uses a fast implementation using [Armadillo C++ library for linear algebra & scientific computing](https://arma.sourceforge.net/) to implement an optimized kdtree structure and implicit nearest neighbour search-by-radius. [OpenMP](https://www.openmp.org/ "https://www.openmp.org/") is used for parallelization to take advantage of multiple processors or cores, leading to 4x to 100x increased speed depending on the architecture of the machine and the characteristics of the 3D point cloud.

The following geometric features are currently implemented:

-   1st eigenvalue
-   2nd eigenvalue
-   3rd eigenvalue
-   Planarity
-   Linearity
-   Sphericity
-   Verticality
-   Sum of eigenvalues
-   Omnivariance
-   Eigenentropy
-   Anisotropy
-   PCA 1
-   PCA 2
-   Surface variation

## Installation

You can install the development version of CloudGeometry from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fpirotti/CloudGeometry")
```

## Example

This is a basic example which shows you how to process a lidar point cloud:

``` r

## basic example code

library(CloudGeometry)
gf <- calcGF(lidar, 2, TRUE)
```

![](images/clipboard-3962452338.png)

## Future developments

More feature descriptors such as Fast Point Feature Histograms (FPFH) as implemented in [point cloud library](https://pcl.readthedocs.io/) will be coming up.

## References

<a id="1">[1]</a> Weinmann, M., Jutzi, B., Hinz, S., Mallet, C., 2015. Semantic point cloud interpretation based on optimal neighborhoods, relevant features and efficient classifiers. ISPRS J. Photogramm. Remote Sens. 105, 286–304. <https://doi.org/10.1016/j.isprsjprs.2015.01.016>

<a id="2">[2]</a> Weinmann, M., Jutzi, B., Mallet, C., 2014. Semantic 3D scene interpretation: A framework combining optimal neighborhood size selection with relevant features. ISPRS Annals of the Photogrammetry, Remote Sensing and Spatial Information Sciences II–3, 181–188. <https://doi.org/10.5194/isprsannals-II-3-181-2014>

<a id="3">[3]</a>  Brodu, N., Lague, D., 2012. 3D terrestrial lidar data classification of complex natural scenes using a multi-scale dimensionality criterion: Applications in geomorphology. ISPRS Journal of Photogrammetry and Remote Sensing 68, 121–134. <https://doi.org/10.1016/j.isprsjprs.2012.01.006> 
