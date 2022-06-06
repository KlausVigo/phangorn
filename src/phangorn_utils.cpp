#include <Rcpp.h>
using namespace Rcpp;


#define DINDEX(i, j) n*(i - 1) - i * (i - 1)/2 + j - i - 1

int give_index3(int i, int j, int n)
{
  if (i > j) return(DINDEX(j, i));
  else return(DINDEX(i, j));
}


// [[Rcpp::export]]
NumericVector fhm_new(NumericVector v, int n){
  unsigned int l, i, j;
  unsigned int start, step, num_splits;
  unsigned int max_n = (unsigned int)n;
  double vi, vj;
  num_splits = (1 << n);
  step = 1;
  for(l = 0; l < max_n; l++){
    start = 0L;
    while(start < (num_splits-1L)){
      for(i = start; i < (start + step); i++){
        j = i + step;
        vi = v[i];
        vj = v[j];
        v[i] = vi + vj;
        v[j] = vi - vj;
      }
      start = start + 2*step;
    }
    step *= 2;
  }
  return v;
}



//library(tools)
//package_native_routine_registration_skeleton("package-root-directory")

// import: edge matrix, number of tips
// export: Descendants(x, 1:max(x$edge), "all")
// [[Rcpp::export]]
List allDescCPP(IntegerMatrix orig, int nTips) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent);
    // create list for results
    std::vector< std::vector<int> > out(m) ;
    for(int i = 0; i<nTips; i++){
        out[i].push_back(i + 1L);
    }
    std::vector<int> y;
    for(int i = 0; i<parent.size(); i++){
        out[parent[i]-1L].push_back(children[i]);
        if(children[i] > nTips){
            y = out[children[i] -1L];
            out[parent[i]-1L].insert( out[parent[i]-1L].end(), y.begin(), y.end() );
        }
    }
    // return the list
    return wrap(out);
}


//replaces void countCycle(int *M, int *l, int *m, int *res){
// [[Rcpp::export]]
int countCycle_cpp(IntegerMatrix M){
  int j, tmp;
  int l = M.nrow();
  int m = M.ncol();
  int res=0L;
  for (int i=0; i<l; i++) {
    tmp = 0;
    if(M[i] != M[i + (m -1) * l])tmp++;
    for (j=1; j<m; j++) {
      if(M[i + (j-1)* l] != M[i + j * l])tmp++;
    }
    if(tmp>2L) res += tmp;
  }
  return(res);
}

// [[Rcpp::export]]
IntegerVector countCycle2_cpp(IntegerMatrix M){
  int tmp;
  int l = M.nrow();
  int m = M.ncol();
  IntegerVector res(l);
  for (int i=0; i<l; i++) {
    tmp = 0L;
    if(M[i] != M[i + (m -1) * l])tmp=1L;
    for (int j=1; j<m; j++) {
      if(M[i + (j-1L)* l] != M[i + j * l])tmp++;
    }
    res[i]=tmp;
  }
  return(res);
}



// speed up some code for NJ
// [[Rcpp::export]]
IntegerVector out_cpp(NumericVector d, NumericVector r, int n){
  int i, j; //, k, l;
  double res, tmp;
//  k=1;
//  l=2;
  IntegerVector xx = IntegerVector::create(1, 2);
  res = d[1] - r[0] - r[1];
  for(i = 0; i < (n-1); i++){
    for(j = i+1; j < n; j++){
      tmp = d[i*n+j] - r[i] - r[j];
      if(tmp<res){
        xx[0]=i+1;
        xx[1]=j+1;
        res = tmp;
      }
    }
  }
  return(xx);
}



// [[Rcpp::export]]
std::vector<int> getIndex(IntegerVector left, IntegerVector right, int n){
  int k;
  std::vector<int> res;
  for (int i = 0; i < left.size(); i++){
    for (int j = 0; j < right.size(); j++){
      k = give_index3(left[i], right[j], n) + 1L;
      res.push_back(k);
      k++;
    }
  }
  return res;
}


// transfer bootstrap
// [[Rcpp::export]]
double Transfer_Index(const IntegerVector bp, const IntegerMatrix orig, int l) {
  IntegerVector parent = orig( _, 0);
  IntegerVector children = orig( _, 1);
  int m = max(parent), tmp, ei, ni;
  int p = bp.size();
  int lmp = l - p;
  int best = p - 1;
  double result;
  IntegerVector l0(m+1);
  IntegerVector l1(m+1);
  for(int i = 0; i<l; i++) l0[i] = 1;
  for(int i = 0; i<p; i++){
    l0[bp[i]] = 0;
    l1[bp[i]] = 1;
  }
  int node = parent[0];
  for(int i = 0; i<parent.size(); i++){
    ni = parent[i];
    ei = children[i];
    l0[ni] += l0[ei];
    l1[ni] += l1[ei];
    if(ni != node){
      tmp = std::min((p - l1[node]) + l0[node], (lmp - l0[node]) + l1[node]);
      best = std::min(best, tmp);
      if(best == 1){
        result = 1.0 - (best / (p-1.0));
        return(result);
      }
      node = ni;
    }
  }
  tmp = std::min((p - l1[node]) + l0[node], (lmp - l0[node]) + l1[node]);
  best = std::min(best, tmp);
  result = 1.0 - (best / (p-1.0));
  return(result);
}



// import: edge matrix, number of tips
// export: Descendants(x, 1:max(x$edge), "all")
// [[Rcpp::export]]
std::vector< std::vector<int> > bipartCPP(IntegerMatrix orig, int nTips) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent), j=0;
    int nnode = m - nTips;
    // create list for results
    std::vector< std::vector<int> > out(nnode) ;
    std::vector<int> y;
    for(int i = 0; i<parent.size(); i++){
        j = parent[i] - nTips - 1L;
        if(children[i] > nTips){
            y = out[children[i] - nTips -1L];
            out[j].insert( out[j].end(), y.begin(), y.end() );
        }
        else out[j].push_back(children[i]);
    }
    for(int i=0; i<nnode; ++i){
        sort(out[i].begin(), out[i].end());
    }
    return out;    // return the list
}


// replacement for bip maybe more error tolerant slightly slower
// import: edge matrix, number of tips
// export: Descendants(x, 1:max(x$edge), "all")
// [[Rcpp::export]]
std::vector< std::vector<int> > bipCPP(IntegerMatrix orig, int nTips) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent), j=0;
    // create list for results
    std::vector< std::vector<int> > out(m) ;
    std::vector<int> y;
    for(int i = 0; i<nTips; i++){
        out[i].push_back(i + 1L);
    }
    for(int i = 0; i<parent.size(); i++){
        j = parent[i] - 1L;
        if(children[i] > nTips){
            y = out[children[i] - 1L];
            out[j].insert( out[j].end(), y.begin(), y.end() );
        }
        else out[j].push_back(children[i]);
    }
    for(int i=0; i<m; ++i){
        sort(out[i].begin(), out[i].end());
    }
    return out;    // return the list
}


// [[Rcpp::export]]
int bip_shared(SEXP tree1, SEXP tree2, int nTips){
  List M1 = tree1;
  List M2 = tree2;
  IntegerMatrix E1 = M1["edge"];
  IntegerMatrix E2 = M2["edge"];

  std::vector< std::vector<int> > bp1 = bipartCPP(E1, nTips);
  std::vector< std::vector<int> > bp2 = bipartCPP(E2, nTips);

  std::sort(bp1.begin(), bp1.end());
  std::sort(bp2.begin(), bp2.end());

  int shared=0;
  for(auto i=0U, j=0U; i<bp1.size(), j<bp2.size(); ){
    if(bp1[i]==bp2[j]) {
      shared++;
      i++;
      j++;
    }
    else {
      if(bp1[i] < bp2[j]) i++;
      else j++;
    }
  }
  return shared;
}




// shorter and easier to understand replacement of C function
// import: edge matrix
// export: list of children
// [[Rcpp::export]]
List allChildrenCPP(const IntegerMatrix orig) {
    IntegerVector parent = orig( _, 0);
    IntegerVector children = orig( _, 1);
    int m = max(parent);
    // create list for results
    std::vector< std::vector<int> > out(m) ;
    for(int i = 0; i<parent.size(); i++){
        out[parent[i]-1L].push_back(children[i]);
    }
    return wrap(out);
}

//library(tools)
//package_native_routine_registration_skeleton("package-root-directory")
// [[Rcpp::export]]
List allSiblingsCPP(const IntegerMatrix & edge) {
  IntegerVector parent = edge( _, 0);
  int m = max(parent), l, left, right;
  int root = min(parent);
  List ch = allChildrenCPP(edge);
  std::vector< std::vector<int> > out(m); //max(edge)
  for(int h = root-1L; h<m; h++){
    IntegerVector tmp_ch = ch[h];
    l = tmp_ch.size();
    if(l>0){
      for(int j=0; j<l; j++){
        left = tmp_ch[j];
        for(int k=0; k<l; k++) {
          right = tmp_ch[k];
          if(left != right) out[left-1L].push_back(right);
        }
      }
    }
  }
  return wrap(out);
}



// [[Rcpp::export]]
IntegerVector p2dna(NumericMatrix xx, double eps=0.999){
    int nr = xx.nrow(); //xx.ncol(), nc = 4;
    double m=0.0;
    IntegerVector tmp = IntegerVector::create(1,2,4,8);
    IntegerVector res(nr);
    for(int i=0; i<nr; ++i){
        m=xx(i,0);
        for(int j=1; j<4; ++j){
            if(m<xx(i,j)) m=xx(i,j);
        }
        for(int j=0; j<4; ++j){
            if(xx(i,j) > (m * eps)) res(i)+=tmp[j];
        }
    }
    return res;
}


// [[Rcpp::export]]
NumericVector node_height_cpp(IntegerVector edge1, IntegerVector edge2,
                              NumericVector edge_length)
{
    NumericVector xx(max(edge2));
    for (int i = (edge2.size() - 1); i >= 0; i--)
        xx[edge2[i] - 1] = xx[edge1[i] - 1] + edge_length[i];
    return(max(xx) - xx);
}



/*
Fast cophenetic distance
*/


void copheneticHelpCpp(std::vector<int> left, std::vector<int> right, int h, NumericVector nh, int nTips, NumericVector dm){
    int ind;
    for(std::size_t i=0; i<left.size(); i++){
        for(std::size_t j=0; j<right.size(); j++){
            ind = give_index3(left[i], right[j], nTips);
            dm[ind] = 2.0*nh[h] - nh[left[i]-1L] - nh[right[j]-1L];
        }
    }
}


// [[Rcpp::export]]
NumericVector cophenetic_cpp(IntegerMatrix edge, NumericVector edge_length,
                             int nTips, int nNode){
    IntegerVector parents = edge( _, 0);
    IntegerVector children = edge( _, 1);
    NumericVector nh = node_height_cpp(parents, children, edge_length);
    List ch = allChildrenCPP(edge);
    std::vector< std::vector<int> > bip = bipCPP(edge, nTips);
    NumericVector dm( nTips * (nTips-1) /2 );
    int l=0, left, right;
    for(int h=nNode; h<(nNode + nTips); h++){
        // changed from NumericVector to IntegerVector
        IntegerVector tmp_ch = ch[h];
        l = tmp_ch.size();
        for(int j=0; j<(l-1L); j++){
            left = tmp_ch[j] - 1;
            for(int k=j+1L; k<l; k++) {
                right = tmp_ch[k] - 1;
                copheneticHelpCpp(bip[left], bip[right], h, nh, nTips, dm);
            }
        }
    }
    return dm;
}


// For phytools
//' @rdname phangorn-internal
//' @export
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

