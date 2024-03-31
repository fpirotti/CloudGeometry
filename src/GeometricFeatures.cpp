#include "RcppArmadillo.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "kdtree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

using namespace Rcpp;
using namespace std;
using namespace arma;

arma::vec vertical = { 0, 0, 1 };

// [[Rcpp::depends(RcppProgress)]]
int  nn(int idr,
              double radius,
              arma::mat const &x,
              arma::mat &out,
              Kdtree::KdTree &tree,
              bool varRadius=false) {


  if ( Progress::check_abort() )
    return -1;

  std::vector<double> point(3);
  //Kdtree::KdNodeVector result;
  arma::mat result;
  point[0] =  x(idr, 0);
  point[1] =  x(idr, 1);
  point[2] =  x(idr, 2);

  arma::vec eigval;
  arma::mat eigvec;
  double eigenEntropy , minEE, finalRadius, eigenSum, verticality = 0;

  //////// if varRadius == true we have to calculate entropy and get
  // the lowest one... we start with largest and then
  // do some heuristics to limit the loop:  max 100
  if(varRadius){

    int mCounts = 1;
    // std::vector<float> distances;
    tree.arma_k_nearest_neighbors(point, 10, &result);
    eig_sym(eigval, eigvec,  arma::cov( result ));
    eigenSum  = arma::sum(eigval);
    eigval = eigval /  eigenSum;
    eigenEntropy  = ( -1.0*eigval.at(0)*log( (eigval.at(0)+0.00000000001) )
                        -eigval.at(1)*log( (eigval.at(1)+0.00000000001) )
                        -eigval.at(2)*log( (eigval.at(2)+0.00000000001) ) );

    minEE = eigenEntropy;
    finalRadius = radius;
    // double originalRadius = radius;
    double fractionizer = radius;
    int dir = 0;
    int changedSign = 0;

    while(mCounts< 50){
      mCounts++;
      if(result.n_rows < 5){
        if(dir < 0){
          changedSign++;
        }
        dir = 1;
      }

      if(result.n_rows > 500){
        if(dir > 0){
          changedSign++;
        }
        dir = -1;
      }

      if(changedSign >2){
        changedSign=0;
        fractionizer = fractionizer/2;
      }
      radius += (fractionizer*dir);
      tree.arma_range_nearest_neighbors(point, radius, &result);

      if(result.n_rows < 4){
        // Rprintf("Less than four results, breaking the search after %d searches\n", mCounts);
        continue;
      }
      eig_sym(eigval, eigvec,  arma::cov( result ));
      eigenSum  = arma::sum(eigval);
      eigval = eigval /  eigenSum;
      eigenEntropy  = ( -1.0*eigval.at(0)*log( (eigval.at(0)+0.00000000001) )
                            -eigval.at(1)*log( (eigval.at(1)+0.00000000001) )
                            -eigval.at(2)*log( (eigval.at(2)+0.00000000001) ) );

      if(minEE > eigenEntropy) {
        minEE = eigenEntropy;
        finalRadius = radius;
      }
    }
    radius = finalRadius;
    out(idr,8)=mCounts;
    out(idr,9)=fractionizer;
  }

  tree.arma_range_nearest_neighbors(point, radius, &result);

  if(result.n_rows >3){
    eig_sym(eigval, eigvec,  arma::cov( result ));
    eigenSum  = arma::sum(eigval);
    verticality = (dot( eigvec.col(0),vertical)/eigval.at(0) ) / arma::datum::pi * 180.0 ;
    // avoid log -inf adding a tiny amount
    out(idr,0) = eigval.at(2);
    out(idr,1) = eigval.at(1);
    out(idr,2) = eigval.at(0);
    out(idr,3) = eigenSum;

    eigval = eigval /  eigenSum;
    eigenEntropy  = ( -1.0*eigval.at(0)*log( (eigval.at(0)+0.00000000001) )
                          -eigval.at(1)*log( (eigval.at(1)+0.00000000001) )
                          -eigval.at(2)*log( (eigval.at(2)+0.00000000001) ) );

    out(idr,4) = eigenEntropy;
    out(idr,5) = result.n_rows;
    out(idr,6) = verticality;
    if(varRadius) {
      out(idr,7) = radius;
    }
  }

  return result.n_rows;

}




// [[Rcpp::depends(RcppProgress)]]
Kdtree::KdNodeVector cTree(arma::mat const &x,
                           bool progress = true) {

  int nRows = x.n_rows;

  // int counter = (int)(ceil(nRows / 100.0f + 1.0f));

  Kdtree::KdNodeVector nodes;
  nodes.resize(nRows);

  // if(progress){
    Progress p(nRows, progress);
    for (int idr=0; idr<nRows;idr++) {
      if ( p.increment() ) {
        std::vector<double> point(3);
        point[0] =  x(idr, 0);
        point[1] =  x(idr, 1);
        point[2] =  x(idr, 2);
        nodes[idr] = Kdtree::KdNode(point);
      }
    }
  return(nodes);
}

// [[Rcpp::export]]
arma::mat  nnEigen(arma::mat const &x,
                   double radius = 1.0,
                   bool varRadius=false,
                   bool progress = true,
                   bool verbose = false,
                   int threads=0) {


  int noNoNeighbours = 0;
  Kdtree::KdNodeVector result;
  // Kdtree::KdNodeVector nodes;
  int nRows = x.n_rows;

  arma::mat out(   nRows, 10  );
  out.fill(datum::nan);

  if(nRows < 4){
    out.clear();
    return(out);
  }

  // int counter = (int)(ceil(nRows / 100.0f + 1.0f));

  Kdtree::KdNodeVector a = cTree(x, progress);

  if(verbose) REprintf("Creating tree with %lu nodes\n", a.size() );
  Kdtree::KdTree tree(&a, 2, progress);


  if(tree.allnodes.size()==0){
    REprintf("Stopped by user\n" );
    out.clear();
    return(out);
  }

  if(verbose) REprintf("Finished tree with %lu  nodes\n", tree.allnodes.size() );
  if(verbose && varRadius) REprintf("Searching nearest points and defining ideal radius with max=%.3f\n", radius );
  if(verbose && !varRadius) REprintf("Searching nearest points around  radius of %.3f ... \n", radius );

  arma::vec eigval;
  Kdtree::KdNode node;
  int idr; //,  chunk;

#ifdef _OPENMP
  if ( threads > 0 ){
    omp_set_num_threads( threads );
    if(verbose) REprintf("\nNumber of threads = %i as requested.\n", threads);
  } else {
    int nt = omp_get_max_threads();
    int nto = nt;
    if(nt>4) nt = nt - 2 ;
    omp_set_num_threads( nt );
    if(verbose) REprintf(" with no specific number of threads requested, will use most of available threads (%i of %i)\n", nt, nto);
  }
#endif

  Progress p(nRows, progress);
  #pragma omp parallel for schedule(dynamic)
    for ( idr=0; idr<nRows;idr++) {
      if ( ! p.is_aborted() ){
        if ( p.increment() ) { // the only way to exit an OpenMP loop
          int nnc = nn(idr, radius, x, out, tree, varRadius );
          if(nnc<4) noNoNeighbours++;
        }
      }
    }


  p.cleanup();
  double sphereArea = (4/3* arma::datum::pi * pow(radius, 3.0) );
  out.col(5)  = out.col(5) / sphereArea;
  if(noNoNeighbours>0){
    if(verbose) REprintf(
"\n%d points out of %d have only 3 or less neighbours so NaN was assigned"
"\nto their eigenvalues... you might want to consider increasing the radius\n\n",
noNoNeighbours, out.n_rows );
   }

  return out;

}

