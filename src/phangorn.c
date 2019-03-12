/*
 * phangorn.c
 *
 * (c) 2008-2019  Klaus Schliep (klaus.schliep@gmail.com)
 *
 *
 * This code may be distributed under the GNU GPL
 *
 */

# define USE_RINTERNALS

#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>



void countCycle(int *M, int *l, int *m, int *res){
    int j, i, tmp;
    res[0]=0L;
    for (i=0; i<*l; i++) {
        tmp = 0;
        if(M[i] != M[i + (*m -1) * *l])tmp++;
        for (j=1; j<*m; j++) {
            if(M[i + (j-1)* *l] != M[i + j * *l])tmp++;
        }
        if(tmp>2L)res[0]+=tmp;
    }
}


void countCycle2(int *M, int *l, int *m, int *res){
    int j, i, tmp;
    for (i=0; i<*l; i++) {
        tmp = 0L;
        if(M[i] != M[i + (*m -1) * *l])tmp=1L;
        for (j=1; j<*m; j++) {
            if(M[i + (j-1L)* *l] != M[i + j * *l])tmp++;
        }
        res[i]=tmp;
    }
}

/*
void nodeH(int *edge, int *node, double *el, int *l,  double *res){
    int ei, i;
    for (i=*l-1L; i>=0; i--) {
        ei = edge[i] - 1L;
        res[ei] = res[node[i]-1L] + el[ei];
    }
}
*/

// C++
SEXP rowMax(SEXP sdat, SEXP sn, SEXP sk){
    int i, h, n=INTEGER(sn)[0], k=INTEGER(sk)[0];
    double x, *res, *dat;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, n));
    res = REAL(result);
    PROTECT(sdat = coerceVector(sdat, REALSXP));
    dat = REAL(sdat);
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] > x) x=dat[i + h*n];}
        res[i] = x;
        }
    UNPROTECT(2);
    return(result);
}

/*
static R_INLINE void getP00(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double tmp, res;
    for(i = 0; i < m; i++){
        tmp = exp(eva[i] * w * el);
        for(j=0; j<m; j++) evi[i + j*m] *= tmp;
    }
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++) res += ev[i + h*m] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }https://de.wikipedia.org/wiki/R%C3%BCdiger_Schulzki
}


// in getPM2
static R_INLINE void getPP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
//    double tmp[m];
    double *tmp;
    tmp = malloc(m * sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            result[i+j*m] = 0;
            for(h = 0; h < m; h++) result[i+j*m] += ev[i + h*m] * tmp[h] * evi[h + j*m];
        }
    }
    free(tmp);
}
*/

void getdP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double  res; //tmp[m],
    double *tmp;
    tmp = (double *) malloc(m * sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w  * el) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }
    free(tmp);
}


void getdP2(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double  res; //tmp[m],
    double *tmp;
    tmp = (double *) malloc(m * sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }
    free(tmp);
}


void getd2P(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double res; //tmp[m],
    double *tmp;
    tmp = (double *) malloc(m * sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w * el) * (eva[i] * w * el) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }
    free(tmp);
}


void getd2P2(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double  res; //tmp[m],
    double *tmp;
    tmp = (double *) malloc(m * sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = (eva[i] * w) * (eva[i] * w) * exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++)    res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }
    free(tmp);
}

/*
SEXP getPM2(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    if(!isNewList(eig)) error("'eig' must be a list");
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            getPP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1);
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
}
*/

SEXP getdPM(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getdP(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
}


SEXP getdPM2(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getdP2(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
}


SEXP getd2PM(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getd2P(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
}


SEXP getd2PM2(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw;
    int m=INTEGER(nc)[0], l=0;
    double *ws=REAL(w);
    double *edgelen=REAL(el);
    double *eva, *eve, *evei;
    SEXP P, RESULT;
    nel = length(el);
    nw = length(w);
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    PROTECT(RESULT = allocVector(VECSXP, nel*nw));
    double *p;
    if(!isNewList(eig)) error("'dlist' must be a list");
    for(j=0; j<nel; j++){
        for(i=0; i<nw; i++){
            PROTECT(P = allocMatrix(REALSXP, m, m));
            p = REAL(P);
            getd2P2(eva, eve, evei, m, edgelen[j], ws[i], p);
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1); //P
            l++;
        }
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
}


/*
static R_INLINE void emult(double *x, double *y, int n){
    for(int i=0; i<n; i++) x[i]*=y[i];
}



void tabulate(int *x, int *n, int *nbin, int *ans){
    int i, tmp;
    for (i=0; i < *nbin; i++) ans[i]=0L;
    for (i=0; i < *n; i++) {
        tmp = x[i];
        if( (tmp>0) & (tmp<(*nbin+1L)) )
        ans[tmp-1L] ++;
    }
}


void C_reorder(int *from, int *to, int *n, int *sumNode,  int *neworder, int *root){
    int i, j, sum=0, k, Nnode, ind, *ord, *csum, *tips, *stack, z=0;  // l,
    double *parent;
    int m=sumNode[0];
    parent = (double *) R_alloc((*n), sizeof(double));
    tips = (int *) R_alloc(m, sizeof(int));
    ord = (int *) R_alloc((*n), sizeof(int));
    csum = (int *) R_alloc( (m+1), sizeof(int));
    stack = (int *) R_alloc(m, sizeof(int));
    for(j=0;j<(*n);j++) parent[j] = (double)from[j];

    for(j=0;j<(*n);j++) ord[j] = j;
    for(j=0;j<m;j++) tips[j] = 0;

    rsort_with_index(parent, ord, *n);
    tabulate(from, n, sumNode, tips);
    csum[0]=0;
    for(i=0;i<(*sumNode);i++){
        sum+=tips[i];
        csum[i+1] = sum;
    }
    k = (*n)-1;
    Nnode = 0;
    stack[0] = *root;

    while(z > -1){
        j=stack[z];
        if(tips[j]>0){
            for(i=csum[j];i<csum[j+1];i++){
                ind = ord[i];
                neworder[k] = ind + 1;
                stack[z] = to[ind]-1;
                k -=1;
                z++;
            }
            Nnode += 1;
            }
        z--;
    }
    root[0]=Nnode;
}

*/


