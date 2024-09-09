#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


#ifndef _PML_H_
#define _PML_H_



using namespace Rcpp;


arma::umat readPML(const List &xlist, int nSeq, int m);



class PML {
public:
  PML(RObject obj, int k, int nNodes) {
    weight = obj.attr("weight");
    nChar = (int) obj.attr("nr");
    nStates = (int) obj.attr("nc");
    NumericMatrix contr = obj.attr("contrast");
    Rcpp::List xlist(obj);
    nSeq = xlist.size();
    X = readPML(xlist, nSeq, nChar);
  }
  int getNR(void){ return nChar; }
  int getK(void){ return k; }
  arma::umat X;
  arma::cube Y;
  // arma::field Y;
  IntegerVector pscore_nodes;
  NumericVector weight; // Integer??
  int nChar;
  int nStates;
  int nSeq;
  NumericMatrix contrast;
  int k;
};



#endif // _PML_H_
