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



// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppProgress)]]
int  nn(int idr,
              double radius,
              arma::mat const &x,
              arma::mat &out,
              Kdtree::KdTree &tree ) {


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
  if(result.n_rows >3){
    eig_sym(eigval, eigvec,  arma::cov(result));
    out(idr,0) = eigval.at(0);
    out(idr,1) = eigval.at(1);
    out(idr,2) = eigval.at(2);
  }

  return result.n_rows;

}




// [[Rcpp::depends(RcppProgress)]]
Kdtree::KdNodeVector cTree(arma::mat const &x,
                           bool progbar = true) {

  int nRows = x.n_rows;

  Rcout << "Starting the tree2 " << nRows << "\n";
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
                   bool progbar = true,
                   int threads=0) {

  int noNoNeighbours = 0;
  Kdtree::KdNodeVector result;
  // Kdtree::KdNodeVector nodes;
  int nRows = x.n_rows;

  arma::mat out(   nRows, 3  );
  out.fill(datum::nan);

  if(nRows < 4){
    Rcpp::stop(std::to_string(nRows).append(" points in table, not enough...") );
    return(out);
  }

  int counter = (int)(ceil(nRows / 100.0f + 1.0f));

  Kdtree::KdNodeVector a = cTree(x, progbar);
  Kdtree::KdTree tree(&a);

  Rcout << "Finished the tree with size=" << tree.allnodes.size() << " nodes\n";

  std::string mst("Searchin NN with radius ");
  mst = mst.append( std::to_string( radius ) );
  mst = mst.append("... [:bar] ETA: :eta");
  // RProgress::RProgress pb2(mst);

  Progress p(nRows, true);
  Rcout << mst << " \n";

  arma::vec eigval;
  Kdtree::KdNode node;
  int idr,  chunk;

#ifdef _OPENMP
  if ( threads > 0 ){
    omp_set_num_threads( threads );
    REprintf("\nNumber of threads requested=%i\n", threads);
  } else {
    int nt = omp_get_max_threads();
    if(nt>4) nt = nt - 2 ;
    omp_set_num_threads( nt );
    REprintf("\nNo specific number of threads requested, will use all available threads - 2 (%i of %i)\n", nt, omp_get_max_threads());
  }

#endif


#pragma omp parallel for schedule(dynamic)
  for ( idr=0; idr<nRows;idr++) {
    if ( ! p.is_aborted() ){
      if ( p.increment() ) { // the only way to exit an OpenMP loop
        int nnc = nn(idr, radius, x, out, tree  );
        if(nnc<4) noNoNeighbours++;
      }
    }
  }

  if(noNoNeighbours>0){
    Rcpp::warning( std::to_string(noNoNeighbours).append(" points have only  3 or less neighbours so NaN was assigned to eigenvalues ... you might want to consider increasing your radius") );
  }

  return out;
}

