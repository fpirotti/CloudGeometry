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
arma::mat varrad = arma::mat(100*26, 4);

// [[Rcpp::depends(RcppProgress)]]
int  nn(int idr,
              double radius,
              arma::mat const &x,
              arma::mat &out,
              Kdtree::KdTree &tree,
              bool varRadius=false) {


  if ( Progress::check_abort() )
    return -1;

  int k; // in case varRadius==true
  std::vector<double> point(3);
  //Kdtree::KdNodeVector result;
  arma::mat result;
  point[0] =  x(idr, 0);
  point[1] =  x(idr, 1);
  point[2] =  x(idr, 2);

  arma::vec eigval;
  arma::vec tempeigval;
  arma::mat eigvec;
  arma::mat tempeigvec;
  double eigenEntropy , minEE, finalRadius, eigenSum, verticality = 0;

  //////// if varRadius == true we have to calculate entropy and get
  // the lowest one... we start with largest and then
  // do some heuristics to limit the loop:  max 100
  if(varRadius){
    k = (int)radius;
    // std::vector<float> distances;
    //
    tree.arma_k_nearest_neighbors(point, k, &result);
    eig_sym(eigval, eigvec,  arma::cov( result.cols(0,2) ));
    tempeigval = eigval;
    tempeigvec = eigvec;
    eigenSum  = arma::sum(eigval);

    eigval = eigval /  eigenSum;
    eigenEntropy  = ( -1.0*eigval.at(0)*log( (eigval.at(0)+0.00000000001) )
                        -eigval.at(1)*log( (eigval.at(1)+0.00000000001) )
                        -eigval.at(2)*log( (eigval.at(2)+0.00000000001) ) );

    minEE = eigenEntropy;
    // first point is the most distant
    finalRadius = result(0,3);

    int mCounts = 0;
    int steps = 25;
    int step = (int)(floor( (k-10)/steps));

    // if(idr%((int)(x.n_rows/100))==0) {
    //   varrad( (int)(idr/((int)(x.n_rows/100)))*25, 0 ) = finalRadius;
    //   varrad( (int)(idr/((int)(x.n_rows/100)))*25, 1 ) = eigenEntropy;
    //   varrad( (int)(idr/((int)(x.n_rows/100)))*25, 2)  = idr;
    //   varrad( (int)(idr/((int)(x.n_rows/100)))*25, 3)  = 0;
    // }

    while(mCounts< steps){
      mCounts++;

      eig_sym(eigval, eigvec,  arma::cov( result.submat((step*mCounts - 1),0,(k-1),2) ) );
      eigenSum  = arma::sum(eigval);
      // return result.n_rows;


      tempeigval = eigval /  eigenSum;

      eigenEntropy  = ( -1.0*tempeigval.at(0)*log( (tempeigval.at(0)+0.00000000001) )
                          -tempeigval.at(1)*log( (tempeigval.at(1)+0.00000000001) )
                          -tempeigval.at(2)*log( (tempeigval.at(2)+0.00000000001) ) );


      // if(idr%((int)(x.n_rows/100))==0) {
      //   varrad( (int)(idr/((int)(x.n_rows/100)))*25+mCounts, 0 ) =  result((step*mCounts-1),3);
      //   varrad( (int)(idr/((int)(x.n_rows/100)))*25+mCounts, 1 ) = eigenEntropy;
      //   varrad( (int)(idr/((int)(x.n_rows/100)))*25+mCounts, 2)  = idr;
      //   varrad( (int)(idr/((int)(x.n_rows/100)))*25+mCounts, 3)  = mCounts;
      // }

      if(minEE > eigenEntropy) {
        minEE = eigenEntropy;
        finalRadius = result((step*mCounts-1),3);
        out(idr,5) = (k - step*mCounts);
        tempeigval = eigval;
        tempeigvec = eigvec;
      }
    }

    eigval = tempeigval;
    eigvec = tempeigvec;
    eigenSum  = arma::sum(eigval);
    out(idr,7)=finalRadius;


  } else{
    tree.arma_range_nearest_neighbors(point, radius, &result);
    if(result.n_rows >3){
      return result.n_rows;
    }
    eig_sym(eigval, eigvec,  arma::cov( result.cols(0,2) ));
    eigenSum  = arma::sum(eigval);
    out(idr,5) = result.n_rows;
  }

    // avoid log -inf adding a tiny amount
    out(idr,0) = eigval.at(2);
    out(idr,1) = eigval.at(1);
    out(idr,2) = eigval.at(0);
    out(idr,3) = eigenSum;

    eigval = eigval /  eigenSum;
    verticality = acos( 1 - dot( eigvec.col(0),vertical)/eigval.at(0) ) / arma::datum::pi * 180.0 ;

    eigenEntropy  = ( -1.0*eigval.at(0)*log( (eigval.at(0)+0.00000000001) )
                          -eigval.at(1)*log( (eigval.at(1)+0.00000000001) )
                          -eigval.at(2)*log( (eigval.at(2)+0.00000000001) ) );

    out(idr,4) = eigenEntropy;
    out(idr,6) = verticality;

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



  if(varRadius && radius<100){
    REprintf("At least 100 points are required for automatic radius definition \n" );
    REprintf("You have setup rk=%d and with varRadius==true this must be at least 100 \n", (int)radius );
    out.clear();
    return(out);
  }

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
    if(verbose) REprintf("\nNumber of threads = %i as requested.\n", omp_get_max_threads() );
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

  if(varRadius) {
    // varrad.save("A.csv", arma::csv_ascii );
    out.col(5)  = out.col(5) / out(idr,7);
  } else {
    out.col(5)  = out.col(5) / sphereArea;
  }


  if(noNoNeighbours>0){
    if(verbose) REprintf(
"\n%d points out of %d have only 3 or less neighbours so NaN was assigned"
"\nto their eigenvalues... you might want to consider increasing the radius\n\n",
noNoNeighbours, out.n_rows );
   }

  return out;

}

