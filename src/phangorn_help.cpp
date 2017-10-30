#include <Rcpp.h>
using namespace Rcpp;

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
List bipCPP(IntegerMatrix orig, int nTips) {
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
    return wrap(out);    // return the list
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

