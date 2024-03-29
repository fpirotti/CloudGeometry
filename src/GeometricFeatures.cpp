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








// [[Rcpp::export]]
arma::mat  geometricFeaturesCalculate(arma::mat const &x,
                                      arma::mat& ids,
                                      bool progbar) {

  int nRows = x.n_rows;
  int counter = (int)(ceil(nRows / 100.0f + 1.0f));

  RProgress::RProgress pb("Working... [:bar] ETA: :eta");
  if(progbar) {
    pb.set_total(100);
  }
//
//
//   // NumericVector out = no_init(nCols);

  arma::mat out(  3, nRows  );

//
//
//
arma::uvec indices;
arma::mat result;
arma::uvec res;
arma::vec eigval;


for (int idr=0; idr<ids.n_rows;idr++) {

  if(idr%counter==0){
    pb.tick();
  }
  arma::rowvec rowids = ids.row(idr) ;
  arma::vec indices = rowids.elem( find(rowids > 0 ) );
  // at least ten neighbours?
  if(indices.size()<10){
    continue;
  }

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


  eig_sym(eigval,  cov(subset));

  out.col(idr) = eigval; // (arma::arg(eigval));


  }

  return out  ;
}
