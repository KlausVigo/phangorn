#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
int timesTwo(int x) {
   return x * 2;
}


//
// rowMin2 
// sankoffQuartet
// sankoff3B
// pnodes
// mpr
//
//
//
//
//
//

void rowMin3(double *dat, int n, int k, double *res){
    int i, h;  
    double x;
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;               
        }        
    }


void sankoffInternal(double *dat, int n, double *cost, int k, double *result){
    int i, j, h; 
    double tmp[k], x;
    for(i = 0; i < n; i++){
        for(j = 0; j < k; j++){
            for(h = 0; h< k; h++){tmp[h] = dat[i + h*n] + cost[h + j*k];}
            x = tmp[0];
            for(h = 1; h< k; h++) {if(tmp[h]<x) {x=tmp[h];}}
            result[i+j*n] += x;
        }                   
    }        
}    

void sankoffTips2(int *x, double *tmp, int nr, int nc, int nrs, double *result){
    int i, j;
    for(i = 0; i < (nr); i++){ 
        for(j = 0; j < (nc); j++) result[i + j*(nr)] += tmp[x[i] - 1L + j*(nrs)];  
    }
}


/*
Rcpp::List sankoff0(Rcpp::List dlist, Rcpp::NumericVector scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips, SEXP contrast, SEXP nrs){
    int i, n = length(node); //, nt = length(tips);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], mn=INTEGER(mNodes)[0], nrc = INTEGER(nrs)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node), ntips=INTEGER(tips)[0];
    SEXP result, dlist2; //tmp, 
    double *res, *cost, *tmp; // *rtmp,
    tmp = (double *) R_alloc(ncx * nrc, sizeof(double));
    for(j=0; j<(ncx * nrc); j++) tmp[j] = 0.0;
    cost = REAL(scost);  
    sankoff4(REAL(contrast), nrc, cost, ncx, tmp); 

    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(dlist2 = allocVector(VECSXP, mn));
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    res = REAL(result);
    
    for(j=0; j<(nrx * ncx); j++) res[j] = 0.0; 
 
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        if(ni == nodes[i]){            
            if(ei < ntips) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, res);
            else sankoff4(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res);
            }
        else{          
            SET_VECTOR_ELT(dlist2, ni, result);
            UNPROTECT(1); 
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            for(j=0; j<(nrx * ncx); j++) res[j] = 0.0; 
            ni = nodes[i];
            if(ei < ntips) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, res);
            else sankoff4(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res); 
            }
    }
    SET_VECTOR_ELT(dlist2, ni, result);    
    UNPROTECT(2); 
    return(dlist2);
}


*/
