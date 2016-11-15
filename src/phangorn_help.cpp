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


// Descendants(x, 1:max(x$edge), "all")
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
    std::vector<int> x;
    std::vector<int> y;
    int tmp=parent[1];
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


// shorter and easier to understand replacement of C function 
// [[Rcpp::export]]
List allChildrenCPP(IntegerMatrix orig) {
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


