#include <Rcpp.h>
#include <RProgress.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector geometricFeaturesCalculate(NumericVector a, bool progbar) {
  int na = a.size();
  NumericVector xab(na);

  RProgress::RProgress pb("Working... [:bar] ETA: :eta");
  if(progbar) {
    pb.set_total(na);
  }

  for (int i = 0; i < na; i++){
    Rcpp::checkUserInterrupt();
    if(progbar) {
      //Rcout << "  i=" <<  i << " j=" <<  " - " << std::endl;
      pb.tick();
    }
  }
  return xab;
}
