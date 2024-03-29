# CloudGeometry

<!-- badges: start -->

<!-- badges: end -->

The goal of CloudGeometry is to comput geometric features for each point in a point cloud based on a local neighbourhood around each point with user-defined radius.

Geometric features are used to describe a point neighbourhood and are widely used to increase the number of features that describe the 3D scene  [[Weinmann et al. 2015]](#1). This can be used as descriptors in training machine learning models for segmentation and classification.

The following

-   1st eigenvalue

-   2nd eigenvalue

-   3rd eigenvalue

-   Sum of eigenvalues

-   Omnivariance

-   Eigenentropy

-   Anisotropy

-   Planarity

-   Linearity

-   PCA 1

-   PCA 2

-   Surface variation

-   Sphericity

-   Verticality

## Installation

You can install the development version of CloudGeometry from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("fpirotti/CloudGeometry")
```

## Example

This is a basic example which shows you how to process a lidar point cloud:

``` r
library(CloudGeometry)
## basic example code
```

![](images/clipboard-3962452338.png)


<a id="1">[1]</a>   Weinmann, M., Jutzi, B., Hinz, S., Mallet, C., 2015. Semantic
point cloud interpretation based on optimal neighborhoods,
relevant features and efficient classifiers. ISPRS J. Photogramm.
Remote Sens. 105, 286–304.
https://doi.org/10.1016/j.isprsjprs.2015.01.016
