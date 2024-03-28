#include "RcppArmadillo.h"
#include <RProgress.h>
using namespace Rcpp;
using namespace arma;

enum feature { PCA1 = 1,
               PCA2,
               PCA3,
               Linearity,
               Planarity,
               Spericity,
               Omnivariance,
               Anisotropy

};






// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec getEigen(arma::mat M) {
  return arma::eig_sym(M);
}


// [[Rcpp::export]]
arma::mat  geometricFeaturesCalculate(arma::mat&  x,
                                      arma::mat& ids,
                                      bool progbar) {

  int nRows = x.n_rows;
  int counter = (int)(ceil(nRows / 100.0f + 1.0f));
  RProgress::RProgress pb("Working... [:bar] ETA: :eta");
  if(progbar) {
    pb.set_total(nRows);
  }
//
//
//   // NumericVector out = no_init(nCols);

  Rcpp::NumericMatrix out( nRows, 6  );

//
//
//
arma::uvec indices;
arma::mat result;
arma::uvec res;


for (int idr=0; idr<ids.n_rows;idr++) {

  arma::rowvec rowids = ids.row(idr) ;
  arma::vec indices = rowids.elem( find(rowids > 0 ) );

  arma::mat subset(indices.size(), 3, arma::fill::zeros );
  int count=0;
  for(arma::vec::iterator i = indices.begin(); i != indices.end(); ++i){
    int id = *i;
    id--; //this is for 0-based vs 1-based

    subset.at(count, 0) =  x(id, 0);
    subset.at(count, 1) =  x(id, 1);
    subset.at(count, 2) =  x(id, 2);
    count++;
  }

  //res = getEigen(subset);
  // Rcpp::NumericMatrix::Row tmp = x(1, _);
  // Rcout << "  i=" <<  ids.row(idr) <<  " - " << std::endl;
  Rcout << "  i=" <<  subset <<  " - " << std::endl;
  // Rcout << "  i=" <<  subset <<  " - " << std::endl;
  break;
  // arma::mat mat;

   // for (int i=0; i<nRows;i++) {
   //   Rcpp::NumericMatrix::Row tmp = x(i,_);
   //   if(i%counter == 0){
   //   Rcpp::checkUserInterrupt();
   //   if(progbar) {
   //      // Rcout << "  i=" <<  i << " j=" << nRows <<  " - " << std::endl;
   //      pb.tick();
   //    }
   //  }
   // }

  }

  return result  ;
}
