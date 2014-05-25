/* 
 * dist.c
 *
 * (c) 2008-2013  Klaus Schliep (klaus.schliep@gmail.com)
 * 
 * 
 * This code may be distributed under the GNU GPL
 *
 */


# define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include <Rinternals.h>



     
SEXP C_rowMin(SEXP sdat, SEXP sn, SEXP sk){
    int i, h, n=INTEGER(sn)[0], k=INTEGER(sk)[0];  
    double x, *res, *dat;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, n));
    res = REAL(result);
    PROTECT(sdat = coerceVector(sdat, REALSXP));
    dat = REAL(sdat);
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;               
        }
    UNPROTECT(2);
    return(result);        
}


void rowMin2(double *dat, int n,  int k, double *res){
    int i, h;  
    double x;
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;               
        }        
    }

   
void rowMinInt(int *dat, int n,  int k, double *res){
    int i, h;  
    int x;
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;               
        }        
    }


void sankoff4(double *dat, int n, double *cost, int k, double *result){
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


SEXP sankoffQuartet(SEXP dat, SEXP sn, SEXP scost, SEXP sk){
    int j, n=INTEGER(sn)[0], k = INTEGER(sk)[0];  
    double *cost, *res, *rtmp;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, n));
    rtmp = (double *) R_alloc(n*k, sizeof(double));
    res = (double *) R_alloc(n*k, sizeof(double));
    PROTECT(scost = coerceVector(scost, REALSXP));
    cost = REAL(scost);
    for(j=0; j<(n*k); j++) rtmp[j] = 0.0;
    for(j=0; j<(n*k); j++) res[j] = 0.0;   
    sankoff4(REAL(VECTOR_ELT(dat,0)), n, cost, k, rtmp);
    sankoff4(REAL(VECTOR_ELT(dat,1)), n, cost, k, rtmp);
    sankoff4(rtmp, n, cost, k, res);
    sankoff4(REAL(VECTOR_ELT(dat,2)), n, cost, k, res);
    sankoff4(REAL(VECTOR_ELT(dat,3)), n, cost, k, res);
    rowMin2(res, n, k, REAL(result));  //res, sn sk  
    UNPROTECT(2);    
    return(result);        
}    

/*
SEXP sankoffTwin(SEXP dat1, SEXP dat2, SEXP sn, SEXP scost, SEXP sk){
    int j, n=INTEGER(sn)[0], k = INTEGER(sk)[0];  
    double *cost, *tmp;
    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, n, k));
    tmp = REAL(result);
    PROTECT(scost = coerceVector(scost, REALSXP));
    cost = REAL(scost);
    for(j=0; j<(n*k); j++) tmp[j] = REAL(dat1)[j];
    sankoff4(REAL(dat2), n, cost, k, tmp);
    UNPROTECT(2);    
    return(result);        
} 
*/

SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips){
    R_len_t i, n = length(node), nt = length(tips);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], mn=INTEGER(mNodes)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP result, dlist2; //tmp, 
    double *res, *cost; // *rtmp,
    cost = REAL(scost);
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(dlist2 = allocVector(VECSXP, mn));
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    res = REAL(result);
    for(i = 0; i < nt; i++) SET_VECTOR_ELT(dlist2, INTEGER(tips)[i], VECTOR_ELT(dlist, INTEGER(tips)[i]));
    for(j=0; j<(nrx * ncx); j++) res[j] = 0.0; 
 
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        if(ni == nodes[i]){            
            sankoff4(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res);
            }
        else{          
            SET_VECTOR_ELT(dlist2, ni, result);
            UNPROTECT(1); 
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            for(j=0; j<(nrx * ncx); j++) res[j] = 0.0; 
            ni = nodes[i];
            sankoff4(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res); 
            }
    }
    SET_VECTOR_ELT(dlist2, ni, result);    
    UNPROTECT(2); 
    return(dlist2);
}


void sankoffTips(int *x, double *tmp, int nr, int nc, int nrs, double *result){
    int i, j;
    for(i = 0; i < (nr); i++){ 
        for(j = 0; j < (nc); j++) result[i + j*(nr)] += tmp[x[i] - 1L + j*(nrs)];  
    }
}

// faster und memory efficient, aehnlich wie logLik2
//SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips)
/*
 * 
library(phangorn) 
data(Laurasiatherian)
tree = NJ(dist.hamming(Laurasiatherian))
sankoffNew(tree, Laurasiatherian)
 
 */



SEXP sankoff3B(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips, SEXP contrast, SEXP nrs){
    R_len_t i, n = length(node); //, nt = length(tips);
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
// die naechte Zeile vielleicht raus
//    for(i = 0; i < nt; i++) SET_VECTOR_ELT(dlist2, INTEGER(tips)[i], VECTOR_ELT(dlist, INTEGER(tips)[i]));
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


    
SEXP pNodes(SEXP data, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge){
    R_len_t n = length(node); 
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0];
    int k, pj, i, j, start, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP result, dlist;  
    double *res, *tmp, *cost;
    cost = REAL(scost);
    pj = nodes[n-1L];
    start = n-1L;
    PROTECT(dlist = allocVector(VECSXP, length(data)));
    tmp = (double *) R_alloc(nrx*ncx, sizeof(double));    
    for(i=0; i<(nrx * ncx); i++) tmp[i] = 0.0;
    for(j=n-1L; j>=0; j--) {
        PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
        res = REAL(result);
        if (pj != nodes[j]) {
            for(i=0; i<(nrx * ncx); i++) tmp[i] = 0.0;
            sankoff4(REAL(VECTOR_ELT(dlist, nodes[j])), nrx, cost, ncx, tmp);
            for(i=0; i<(nrx * ncx); i++) res[i] = tmp[i] ;  
            pj = nodes[j];
            start = j;
        }
        else for(i=0; i<(nrx * ncx); i++) res[i] = tmp[i] ;
        k = start;
        while (k >= 0 && pj == nodes[k]) {
            if (k != j) 
                sankoff4(REAL(VECTOR_ELT(data, edges[k])), nrx, cost, ncx, res);                
            k--;
        }
        SET_VECTOR_ELT(dlist, edges[j], result);    
        UNPROTECT(1);
    }
    UNPROTECT(1);
    return(dlist);
}



/*

static double *data1;
static double *weight;

void sankoff_init(double *weights, int *nr, int *nc, int *nTips)
{
    int i;
    data1 = (double *) calloc(*nr * *nc * *nTips, sizeof(double));
    weight = (double *) calloc(*nr, sizeof(double));
    for(i=0; i<*nr; i++)weight[i] = weights[i]; 
}


void sankoff_free(){
    free(data1);
    free(weight);
}




Ziel weniger memory und schneller

static R_INLINE void sankoff4(double *dat, int n, double *cost, int k, double *result){
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


void sankoffTips(int *x, double *tmp, int nr, int nc, int nrs, double *result){
    int i, j;
    for(i = 0; i < (nr); i++){ 
        for(j = 0; j < (nc); j++) result[i + j*(nr)] += tmp[x[i] - 1L + j*(nrs)];  
    }
}


// mNodes raus
// SEXP fuer verschiedene IO mit R
// ohne zuviele SEXP definieren, Ausnahme dlist
void SANKOFF4(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips, SEXP contrast, SEXP nrs){
    R_len_t i, n = length(node), nt = length(tips);
    int nTips = (int)length(tips), nii;
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nrc = INTEGER(nrs)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    int rc = ncx*nrx;
    double *cost, *tmp; 
    tmp = (double *) R_alloc(ncx * nrc, sizeof(double)); 
    cost = REAL(scost);
    sankoff4(REAL(contrast), nrc, cost, ncx, tmp); 
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = 0L;     
    for(i = 0; i < n; i++) {
        ei = edges[i]; 
        nii = (nodes[i]-nTips) * rc; 
        if(ni == nodes[i]){            
            if(ei < nt) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, &data1[nii]);
            else sankoff4(&data1[(ei-1L)*rc], nrx, cost, ncx, &data1[nii]);
            }
        else{          
            ni = nodes[i];
            for(j=0; j<rc; j++) data1[nii+j] = 0.0; 
            if(ei < nt) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, &data1[nii]);
            else sankoff4(&data1[(ei-1L)*rc], nrx, cost, ncx, &data1[nii]); 
            }
    }
}

*/


