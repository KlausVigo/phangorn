#include <Rcpp.h>
using namespace Rcpp;


#define DINDEX(i, j) n*(i - 1) - i * (i - 1)/2 + j - i - 1


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


// replacement for bipart maybe more error tolerant and slightly slower
// import: edge matrix, number of tips
// export: Descendants(x, 1:max(x$edge), "all")
// [[Rcpp::export]]
List bipartCPP(IntegerMatrix orig, int nTips) {
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
    return wrap(out);    // return the list
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
int give_index3(int i, int j, int n)
{
    if (i > j) return(DINDEX(j, i));
    else return(DINDEX(i, j));
}


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

