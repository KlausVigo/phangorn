#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
IntegerVector threshStateC(NumericVector x, NumericVector thresholds) {
  int n = x.size(), m = thresholds.size()-1L, j=0L;
  IntegerVector out(n);
  for (int i = 0; i < n; i++) {
    j=0L; 
    while(x[i]>thresholds[j] && j<m ) j++; 
    out[i] = j+1L;
  }
  return out;
}


