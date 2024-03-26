#include <Rcpp.h>
#include <RProgress.h>
using namespace Rcpp;

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
NumericMatrix geometricFeaturesCalculate(NumericMatrix  x,
                                         bool progbar) {

  int nRows = x.nrow();
  int counter = ceil(nRows / 1000);
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
   for (int i=0; i<nRows;i++) {
     NumericMatrix::Row tmp = x(i,_);
//     if(i%counter){
//       Rcpp::checkUserInterrupt();
      if(progbar) {
        Rcout << "  i=" <<  i << " j=" << nRows <<  " - " << std::endl;
        pb.tick();
      }
//     }
   }

  return out;
}
