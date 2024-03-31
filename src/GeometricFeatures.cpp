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
  tree.arma_range_nearest_neighbors(point, radius, &result);
  arma::vec eigval;
  arma::mat eigvec;
  double eigenEntropy , eigenSum, verticality = 0;
  if(result.n_rows >3){
    eig_sym(eigval, eigvec,  arma::cov( result ));
    eigenSum  = arma::sum(eigval);
    verticality = dot( eigvec.col(2),vertical) / arma::datum::pi * 180.0 ;
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
  }

  return result.n_rows;

}




// [[Rcpp::depends(RcppProgress)]]
Kdtree::KdNodeVector cTree(arma::mat const &x,
                           bool progress = true) {

  int nRows = x.n_rows;

  int counter = (int)(ceil(nRows / 100.0f + 1.0f));

  Progress p(nRows, true);

  Kdtree::KdNodeVector nodes;
  nodes.resize(nRows);
  for (int idr=0; idr<nRows;idr++) {
    if ( p.increment() ) {
      std::vector<double> point(3);
      point[0] =  x(idr, 0);
      point[1] =  x(idr, 1);
      point[2] =  x(idr, 2);
      // nodes.push_back(Kdtree::KdNode(point));
      nodes[idr] = Kdtree::KdNode(point);
    }

  }

  // Kdtree::KdTree tree(&nodes);
  return(nodes);
}

// [[Rcpp::export]]
arma::mat  nnEigen(arma::mat const &x,
                   double radius = 1.0,
                   bool varRadius=false,
                   bool progress = true,
                   int threads=0) {

  int noNoNeighbours = 0;
  Kdtree::KdNodeVector result;
  // Kdtree::KdNodeVector nodes;
  int nRows = x.n_rows;

  // out matrix has
  //   3 eigenvalues
  //   eigenEntropy
  //
  arma::mat out(   nRows, 7  );
  out.fill(datum::nan);

  if(nRows < 4){
    out.clear();
    return(out);
  }

  int counter = (int)(ceil(nRows / 100.0f + 1.0f));

  Kdtree::KdNodeVector a = cTree(x, progress);

  REprintf("Creating tree with %lu nodes\n", a.size() );
  Kdtree::KdTree tree(&a);

  if(tree.allnodes.size()==0){
    REprintf("Stopped by user\n" );
    out.clear();
    return(out);
  }

  REprintf("Finished tree with %lu  nodes\n", tree.allnodes.size() );
  REprintf("Searching NN with radius of %.3f ... ", radius );

  arma::vec eigval;
  Kdtree::KdNode node;
  int idr,  chunk;

#ifdef _OPENMP
  if ( threads > 0 ){
    omp_set_num_threads( threads );
    REprintf(" with number of threads = %i as requested.\n", threads);
  } else {
    int nt = omp_get_max_threads();
    int nto = nt;
    if(nt>4) nt = nt - 2 ;
    omp_set_num_threads( nt );
    REprintf(" with no specific number of threads requested, will use most of available threads (%i of %i)\n", nt, nto);
  }

#endif


Progress p(nRows, true);

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
  if(noNoNeighbours>0){
    REprintf("%d points out of %d have only 3 or less neighbours so NaN was assigned to eigenvalues ... you might want to consider increasing your radius\n",
             noNoNeighbours, out.n_rows );
   }

  return out;

}

