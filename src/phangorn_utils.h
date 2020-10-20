#ifndef PHANGORNUTILS_H
#define PHANGORNUTILS_H

#include <Rcpp.h>
using namespace Rcpp;

List allDescCPP(IntegerMatrix orig, int nTips);

List bipartCPP(IntegerMatrix orig, int nTips);

std::vector< std::vector<int> > bipCPP(IntegerMatrix orig, int nTips);

List allChildrenCPP(const IntegerMatrix orig);

List allSiblingsCPP(const IntegerMatrix & edge);

IntegerVector p2dna(NumericMatrix xx, double eps=0.999);

NumericVector node_height_cpp(IntegerVector edge1, IntegerVector edge2,
                              NumericVector edge_length);

NumericVector cophenetic_cpp(IntegerMatrix edge, NumericVector edge_length,
                             int nTips, int nNode);

IntegerVector threshStateC(NumericVector x, NumericVector thresholds);

int countCycles_cpp(IntegerMatrix M);

std::vector<int> getIndex(IntegerVector left, IntegerVector right, int n);

#endif
