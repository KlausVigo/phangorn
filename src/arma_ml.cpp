#include "pml.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp ;

const double ScaleEPS = 1.0/4294967296.0;
const double ScaleMAX = 4294967296.0;
const double LOG_SCALE_EPS = -22.18070977791824915926;

// to do optE (unrooted trees), optQuartet (unrooted trees),
// likelihoodQuartet (rooted trees), rnodes like function
// pmlE, for partitions, mixtures, codons

//' Test subsetting X[ind,] funktioniert!!!
//' @noRd
// [[Rcpp::export]]
arma::mat subset_X(const arma::mat& X, const arma::uvec& ind){
  arma::mat out = X.rows(ind-1);
  return out;
}

// (contr %*% P)[ind,]  matp
// [[Rcpp::export]]
arma::mat contr_P_ind(const arma::mat& contr, const arma::mat& P, const arma::uvec& ind){
  arma::mat tmp = (contr * P);
  arma::mat out = tmp.rows(ind-1);
  return out;
}


void scaleMatrix(arma::mat& X, arma::vec& sc){
  int i, j;
  int nc = X.n_cols, nr = X.n_rows;
  double tmp;
  for(i = 0; i < nr; i++) {
    tmp = 0.0;
    for(j = 0; j <  nc;j++){
      tmp += X(i, j);
    }
    while(tmp < ScaleEPS && tmp > 0.0){
      for(j = 0; j < nc; j++) X(i , j) *=ScaleMAX;
      sc[i] +=1L;
      tmp *= ScaleMAX;
    }
  }
}


// [[Rcpp::export]]
arma::mat getP_arma(const arma::vec& e_value, const arma::mat& e_vec, const arma::mat& e_vec_inv, double el, double g){
  arma::mat out = e_vec * arma::diagmat(exp(el*g*e_value)) * e_vec_inv;
  return out;
}


// [[Rcpp::export]]
arma::vec rowMinScale_arma(arma::mat& dat, int k){
 int n = dat.n_rows;
 double tmp;
 arma::vec res(n);
 for(int i = 0; i < n; i++){
   tmp = dat[i];
   for(int h = 1; h< k; h++) {if(dat(i , h) < tmp) tmp=dat(i, h);}
   if(tmp>0L){for(int h = 0; h< k; h++) dat(i , h) -= tmp;}
   res(i) = tmp;
   }
 return res;
}


// Workhorse log-likelihood computation
// [[Rcpp::export]]
arma::vec lll(const arma::umat& X, arma::cube& Y, arma::vec& SC,
              const arma::vec& e_value, const arma::mat& e_vec, const arma::mat& e_vec_inv,
              const arma::mat& contrast,
              const arma::vec& el, const arma::uvec& node, const arma::uvec& edge,
              double g, double w, int nTips, const arma::vec& bf){ //, int nr, int nc, int nnode
//  int n_nodes = Y.n_slices;
  int n_rows  = Y.n_rows;
  int n_cols  = Y.n_cols;
  int n_edge  = edge.n_elem;
  int mnodep1 = max(node) + 1L;
  int ni = mnodep1;
  arma::mat P = arma::mat(n_cols, n_cols);
  arma::mat tmp = arma::mat(n_rows, n_cols);
  arma::mat tmp_c = arma::mat(contrast.n_rows, n_cols);
  for(int i = 0; i < n_edge; i++) {
    P = getP_arma(e_value, e_vec, e_vec_inv, el[i], g);
    int ei = edge[i];
    if(ni != node[i]){
      if(ni < mnodep1){
        // scale here
        scaleMatrix(tmp, SC);
        Y.slice(ni) = tmp;
      }
      ni = node[i];
      if(ei < nTips){
        tmp = contr_P_ind(contrast, P, X.col(ei));
      }
      else{
        tmp = Y.slice(ei-nTips) * P;
      }
    }
    else{
      if(ei < nTips)
        tmp %= contr_P_ind(contrast, P, X.col(ei));
      else
        tmp %= Y.slice(ei-nTips) * P;
    }
  }
  // scale here
  scaleMatrix(tmp, SC);
  return tmp*bf;
}

// like lll with proper scaling allowing for edge optimisation
// [[Rcpp::export]]
arma::vec lll_E(const arma::umat& X, arma::cube& Y, arma::vec& SC,
               const arma::vec& e_value, const arma::mat& e_vec, const arma::mat& e_vec_inv,
               const arma::mat& contrast,
               const arma::vec& el, const arma::uvec& node, const arma::uvec& edge,
               double g, double w, int nTips, const arma::vec& bf,
               arma::mat& SCM){ //, int nr, int nc, int nnode
  //  int n_nodes = Y.n_slices;
  int n_rows  = Y.n_rows;
  int n_cols  = Y.n_cols;
  int n_edge  = edge.n_elem;
  int mnodep1 = max(node) + 1L;
  int ni = mnodep1;
  arma::mat P = arma::mat(n_cols, n_cols);
  arma::mat tmp = arma::mat(n_rows, n_cols);
  arma::mat tmp_c = arma::mat(contrast.n_rows, n_cols);
  SC *= 0;
  for(int i = 0; i < n_edge; i++) {
    P = getP_arma(e_value, e_vec, e_vec_inv, el[i], g);
    int ei = edge[i];
    if(ni != node[i]){
      if(ni < mnodep1){
        // scale here
        scaleMatrix(tmp, SC);
        SCM.col(ni) = SC;
        Y.slice(ni) = tmp;
        SC *= 0;
      }
      ni = node[i];
      if(ei < nTips){
        tmp = contr_P_ind(contrast, P, X.col(ei));
      }
      else{
        tmp = Y.slice(ei-nTips) * P;
        SC += SCM.col(ei-nTips);
      }
    }
    else{
      if(ei < nTips)
        tmp %= contr_P_ind(contrast, P, X.col(ei));
      else{
        tmp %= Y.slice(ei-nTips) * P;
        SC += SCM.col(ei-nTips);
      }
    }
  }
  // scale here
  scaleMatrix(tmp, SC);
  SCM.col(ni) = SC;
  Y.slice(ni) = tmp;

  return tmp*bf;
}


// lll with k rates (mixtures need always everything)
// could be building block for partitions
// needs a version for partitions and mixtures

// [[Rcpp::export]]
arma::vec PML_E(const arma::umat& X, arma::cube& Y, //arma::vec& SC,
                const arma::vec& e_value, const arma::mat& e_vec, const arma::mat& e_vec_inv,
                const arma::mat& contrast,
                const arma::vec& el, const arma::uvec& node, const arma::uvec& edge,
                const arma::vec& g, const arma::vec& w, int nTips, const arma::vec& bf,
                arma::cube& SCM, int nNodes){
  int n_rows  = Y.n_rows;
  int n_cols  = Y.n_cols;
  int k = g.n_elem;

  arma::mat SC = arma::mat(n_rows, k);
  arma::vec sc = arma::vec(n_rows);

  arma::mat tmp = arma::mat(n_rows, k);
  arma::vec res = arma::vec(n_rows);
  for(int i=0; i<k; i++){
    arma::vec tmp_sc = SC.col(i);
    tmp.col(i) = lll_E(X, Y,//Y.slices(i*nNodes, ((i+1L)*nNodes)-1L ),
                    tmp_sc, e_value, e_vec, e_vec_inv,
                    contrast, el, node, edge,
                    g[i], w[i], nTips, bf, SCM.slice(i));
  }
  sc = rowMinScale_arma(SC, k);
//  for(int i=0; i<n_rows; i++){
//    res[i]=0.0;
//    for(int j=0;j<k;j++) res[i] += w[j] * exp(LOG_SCALE_EPS * SC(i, j)) * tmp(i, j);
  for(int i=0; i<k; i++) res += w[i] * exp(LOG_SCALE_EPS * SC.col(i)) % tmp.col(i);
//  for(int i=0; i<n_rows; i++) res[i] = log(res[i]) + LOG_SCALE_EPS * sc[i];
  res = log(res) + LOG_SCALE_EPS * sc;
  return res;
}


// [[Rcpp::export]]
arma::vec NR_f_arma(const arma::vec& eva, double el, const arma::vec& w,
                const arma::vec& g, const arma::cube& X, int k){
  int n= X.n_rows;
  arma::vec res(n);
  for(int j=0; j<k; j++) res += w[j] * X.slice(j) * exp(eva * g[j] * el);
  return res.as_col();
}

// [[Rcpp::export]]
arma::vec NR_df_arma(const arma::vec& eva, double el, const arma::vec& w,
            const arma::vec& g, const arma::cube& X, int k, const arma::vec& f){
  int n= X.n_rows;
  arma::vec res(n);
  for(int j=0; j<k; j++){
    res += w[j] * X.slice(j) * ((eva * g[j] * el) % exp(eva * g[j] * el));
    // res += X.slice(j) * ((w[j] * eva * g[j] * el) % exp(eva * g[j] * el));
  }
  res/=f;
  return res.as_col();
}

// [[Rcpp::export]]
arma::vec NR_d2f_arma(const arma::vec& eva, double el, const arma::vec& w,
            const arma::vec& g, const arma::cube& X, int k, const arma::vec& f){
  int n= X.n_rows;
  arma::vec res(n);
  for(int j=0; j<k; j++){
    res += w[j] * X.slice(j) * ((eva * g[j]) % exp(eva * g[j] * el));
  }
  res/=f;
  return res.as_col();
}


// [[Rcpp::export]]
const arma::vec fs_3(const arma::vec& eva, double el, const arma::vec& w,
                     const arma::vec& g, const arma::cube& X, int ld, const arma::vec& weight,
                     const arma::vec& f0, double tau) //, int mkv)
{
  int mkv=0;
  double edle, ledle, newedle, eps=10;
  double ll=0.0, lll, delta=0.0, scalep = 1.0, l1=0.0, l0=0.0;
  double y, p0=0.0; //sum_wgt=0.0,
  int i, k=0;

  edle = el; //REAL(el)[0];
  int nr = X.n_rows;
  //  int nc = X.n_cols;

  arma::vec f = f0 + NR_f_arma(eva, edle, w, g, X, ld);
  if(mkv){
    p0 = sum(f);
    f /= (1-p0);
  }
  l0 = sum( weight % log(f));
  while ( (eps > 1e-05) &&  (k < 10) ) {
    if(scalep>0.6){
      arma::vec tmp = NR_df_arma(eva, edle, w, g, X, ld, f);
      ll=0.0;
      lll=0.0;
      for(i=0; i<nr ;i++){
        y = weight[i]*tmp[i];
        ll+=y;
        lll+=y*tmp[i];
      }
      delta = ((ll/lll) < 3) ? (ll/lll) : 3;
    } // end if
    ledle = log(edle) + scalep * delta;
    newedle = exp(ledle);
    if (newedle > 10.0) newedle = 10.0;
    if (newedle < tau) newedle = tau;
    f = f0 + NR_f_arma(eva, newedle, w, g, X, ld);
    if(mkv){
      p0 = sum(f);
      f /= (1-p0);
    }
    l1 = sum( weight % log(f));
    eps = l1 - l0;
    // some error handling
    if (eps < 0 || ISNAN(eps)) {
      if (ISNAN(eps))eps = 0;
      else {
        scalep = scalep/2.0;
        eps = 1.0;
      }
      newedle = edle;
      l1 = l0;
    }
    else scalep = 1.0;
    edle=newedle;
    l0 = l1;
    k ++;
  }
  // variance n
  arma::vec tmp2 = NR_d2f_arma(eva, edle, w, g, X, ld, f);
  lll = sum(weight % tmp2 % tmp2);
  arma::vec res(3);
  res[0] = edle;
  res[1] = 1.0 / lll;
  res[2] = l0; //l0
  return res;
}



// [[Rcpp::export]]
arma::umat readPML(const List &xlist, int nSeq, int m){
  arma::umat out(m, nSeq);
  unsigned int tmp;
  for(int i=0; i<nSeq; ++i) {
    Rcpp::IntegerVector y(xlist[i]);
    for(int j=0; j<m; ++j){
      tmp = unsigned(y[j]);
      out(j,i) = tmp;
      //      out[i].push_back(tmp);
    }
  }
  return out; // wrap(out)
}




RCPP_MODULE(PML_mod) {
  using namespace Rcpp;
  class_<PML>("PML")
    .constructor<RObject, int, int>("Default constructor")
    .property("get_nr", &PML::getNR)
  ;
}
