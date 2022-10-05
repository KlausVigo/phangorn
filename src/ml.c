/*
 * ml.c
 *
 * (c) 2008-2021  Klaus Schliep (klaus.schliep@gmail.com)
 *
 *
 * This code may be distributed under the GNU GPL
 *
 */
#define USE_FC_LEN_T

#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif


#define LINDEX(i, k) (i - ntips - 1L) * (nr * nc) + k * ntips * (nr * nc)
// index for LL
#define LINDEX2(i, k) (i - *ntips - 1L) * (*nr* *nc) + k * *ntips * (*nr * *nc)
// index for scaling matrix SCM
#define LINDEX3(i, j) (i - *ntips - 1L) * *nr + j * *ntips * *nr  //nr statt *nr


double one = 1.0, zero = 0.0;
int ONE = 1L;
const double ScaleEPS = 1.0/4294967296.0;
const double ScaleMAX = 4294967296.0;
const double LOG_SCALE_EPS = -22.18070977791824915926;


static double *LL;  //, *WEIGHTS;
static int *SCM; //, *XXX;


SEXP ll_free2(void){
    free(LL);
    free(SCM);
    return R_NilValue;
}


SEXP ll_init2(SEXP nr, SEXP nTips, SEXP nc, SEXP k)
{
    int nrx = INTEGER(nr)[0], nTipsx = INTEGER(nTips)[0], ncx = INTEGER(nc)[0], kx = INTEGER(k)[0];
    int i;
    LL = (double *) calloc(nrx * ncx * kx * nTipsx, sizeof(double));
    SCM = (int *) calloc(nrx * kx * nTipsx, sizeof(int));  // * 2L
    for(i =0; i < (nrx * kx * nTipsx); i++) SCM[i] = 0L;
    return R_NilValue;
}

/*
LL likelihood for internal edges
SCM scaling coefficients
*/


void matm(int *x, double *contrast, int *nr, int *nc, int *nco, double *result){
    int i, j;
    for(i = 0; i < (*nr); i++){
        for(j = 0; j < (*nc); j++) result[i + j*(*nr)] *= contrast[x[i] - 1L + j*(*nco)];
    }
}


SEXP invSites(SEXP dlist, SEXP nr, SEXP nc, SEXP contrast, SEXP nco){
    R_len_t n = length(dlist);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], i, j;
    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    double *res;
    res = REAL(result);
    for(j=0; j < (nrx * ncx); j++) res[j] = 1.0;
    for(i=0; i < n; i++) matm(INTEGER(VECTOR_ELT(dlist, i)), REAL(contrast), INTEGER(nr), INTEGER(nc), INTEGER(nco), res);
    UNPROTECT(1); // result
    return(result);
}


void scaleMatrix(double *X, int *nr, int *nc, int *result){
    int i, j;
    double tmp;
    for(i = 0; i < *nr; i++) {
        tmp = 0.0;
        for(j = 0; j < *nc; j++) tmp += X[i + j* *nr];
        while(tmp < ScaleEPS){
           for(j = 0; j < *nc; j++) X[i + j* *nr] *=ScaleMAX;
           result[i] +=1L;
           tmp *= ScaleMAX;
       }
    }
}


// contrast to full dense matrix
void matp(int *x, double *contrast, double *P, int *nr, int *nc, int *nrs, double *result){
    int i, j;
    double *tmp;
    tmp = (double *) R_alloc((*nc) *(*nrs), sizeof(double));
    F77_CALL(dgemm)("N", "N", nrs, nc, nc, &one, contrast, nrs, P, nc, &zero, tmp, nrs FCONE FCONE);
    for(i = 0; i < (*nr); i++){
        for(j = 0; j < (*nc); j++) result[i + j*(*nr)] = tmp[x[i] - 1L + j*(*nrs)];
    }
}


void rowMinScale(int *dat, int n,  int k, int *res){
    int i, h;
    int tmp;
    for(i = 0; i < n; i++){
        tmp = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < tmp) tmp=dat[i + h*n];}
        if(tmp>0L){for(h = 0; h< k; h++) dat[i + h*n] -= tmp;}
        res[i] = tmp;
    }
}


void getP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double res;
    double *tmp;
    tmp = (double *) R_alloc(m, sizeof(double));
    for(i = 0; i < m; i++) tmp[i] = exp(eva[i] * w * el);
    for(i = 0; i < m; i++){
        for(j = 0; j < m; j++){
            res = 0.0;
            for(h = 0; h < m; h++) res += ev[i + h*m] * tmp[h] * evi[h + j*m];
            result[i+j*m] = res;
        }
    }
}


SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w){
    R_len_t i, j, nel, nw, k;
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
            if(edgelen[j]==0.0 || ws[i]==0.0){
                for(k=0; k<(m*m);k++)REAL(P)[k]=0.0;
                for(k=0; k<m; k++)REAL(P)[k+k*m]=1.0;
            }
            else getP(eva, eve, evei, m, edgelen[j], ws[i], REAL(P));
            SET_VECTOR_ELT(RESULT, l, P);
            UNPROTECT(1);
            l++;
        }
    }
    UNPROTECT(1);//RESULT
    return(RESULT);
}


void lll(SEXP dlist, double *eva, double *eve, double *evei, double *el, double g, int *nr, int *nc, int *node, int *edge, int nTips, double *contrast, int nco, int n, int *scaleTmp, double *bf, double *TMP, double *ans){
    int  ni, ei, j, i, rc; //    R_len_t i, n = length(node);
    double *rtmp, *P;
    ni = -1;
    rc = *nr * *nc;
    rtmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    for(j=0; j < *nr; j++) scaleTmp[j] = 0L;
    for(i = 0; i < n; i++) {
        getP(eva, eve, evei, *nc, el[i], g, P);
        ei = edge[i];
        if(ni != node[i]){
            if(ni>0)scaleMatrix(&ans[ni * rc], nr, nc, scaleTmp); // (ni-nTips)
            ni = node[i];
            if(ei < nTips)
                matp(INTEGER(VECTOR_ELT(dlist, ei)), contrast, P, nr, nc, &nco, &ans[ni * rc]);
            else
                F77_CALL(dgemm)("N", "N", nr, nc, nc, &one, &ans[(ei-nTips) * rc], nr, P, nc, &zero, &ans[ni * rc], nr FCONE FCONE);
        }
        else {
            if(ei < nTips)
                matp(INTEGER(VECTOR_ELT(dlist, ei)), contrast, P, nr, nc, &nco, rtmp);
            else
                F77_CALL(dgemm)("N", "N", nr, nc, nc, &one, &ans[(ei-nTips) * rc], nr, P, nc, &zero, rtmp, nr FCONE FCONE);
            for(j=0; j < rc; j++) ans[ni * rc + j] *= rtmp[j];
        }
    }
    scaleMatrix(&ans[ni * rc], nr, nc, scaleTmp);
    F77_CALL(dgemv)("N", nr, nc, &one, &ans[ni * rc], nr, bf, &ONE, &zero, TMP, &ONE FCONE);
}


// this seems to work perfectly
void lll3(SEXP dlist, double *eva, double *eve, double *evei, double *el, double g, int *nr, int *nc, int *node, int *edge,
    int nTips, double *contrast, int nco, int n, int *scaleTmp, double *bf, double *TMP, double *ans, int *SC){
    int  ni, ei, j, i, rc; //    R_len_t i, n = length(node);
    double *rtmp, *P;
    ni = -1L;
    rc = *nr * *nc;
    rtmp = (double *) R_alloc(*nr * *nc, sizeof(double));
    P = (double *) R_alloc(*nc * *nc, sizeof(double));
    for(j=0; j < *nr; j++) scaleTmp[j] = 0L;
    for(i = 0; i < n; i++) {
        getP(eva, eve, evei, *nc, el[i], g, P);
        ei = edge[i];
        if(ni != node[i]){
// test for node[i+1]
            if(ni>0)scaleMatrix(&ans[ni * rc], nr, nc, &SC[ni * *nr]); // (ni-nTips)
            ni = node[i];
            for(j=0; j < *nr; j++) SC[j + ni * *nr] = 0L;
            if(ei < nTips)
                matp(INTEGER(VECTOR_ELT(dlist, ei)), contrast, P, nr, nc, &nco, &ans[ni * rc]);
            else{
                F77_CALL(dgemm)("N", "N", nr, nc, nc, &one, &ans[(ei-nTips) * rc], nr, P, nc, &zero, &ans[ni * rc], nr FCONE FCONE);
                for(j=0; j < *nr; j++) SC[ni * *nr + j] = SC[(ei-nTips) * *nr + j];
            }
        }
        else {
            if(ei < nTips)
                matp(INTEGER(VECTOR_ELT(dlist, ei)), contrast, P, nr, nc, &nco, rtmp);
            else{
                F77_CALL(dgemm)("N", "N", nr, nc, nc, &one, &ans[(ei-nTips) * rc], nr, P, nc, &zero, rtmp, nr FCONE FCONE);
                for(j=0; j < *nr; j++) SC[ni * *nr + j] += SC[(ei-nTips) * *nr + j];
            }
            for(j=0; j < rc; j++) ans[ni * rc + j] *= rtmp[j];
        }
    }
    scaleMatrix(&ans[ni * rc], nr, nc, &SC[ni * *nr]);
    for(j=0; j < *nr; j++) scaleTmp[j] = SC[ni * *nr + j];

    F77_CALL(dgemv)("N", nr, nc, &one, &ans[ni * rc], nr, bf, &ONE, &zero, TMP, &ONE FCONE);
}


SEXP PML0(SEXP dlist, SEXP EL, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP nco, SEXP contrast, SEXP N){
    int nr=INTEGER(NR)[0], nc=INTEGER(NC)[0], k=INTEGER(K)[0], i, indLL;
    int nTips = INTEGER(NTips)[0], *SC;
    double *g=REAL(G), *tmp, logScaleEPS;
    SEXP TMP;
    double *eva, *eve, *evei;
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    SC = (int *) R_alloc(nr * k, sizeof(int));

    PROTECT(TMP = allocMatrix(REALSXP, nr, k)); // changed
    tmp=REAL(TMP);
    for(i=0; i<(k*nr); i++)tmp[i]=0.0;
    indLL = nr * nc * nTips;
    for(i=0; i<k; i++){
        lll(dlist, eva, eve, evei, REAL(EL), g[i], &nr, &nc, INTEGER(node), INTEGER(edge), nTips, REAL(contrast), INTEGER(nco)[0], INTEGER(N)[0], &SC[nr * i], REAL(bf), &tmp[i*nr], &LL[indLL *i]);
    }
    logScaleEPS = log(ScaleEPS);
    for(i=0; i<(k*nr); i++) tmp[i] = logScaleEPS * SC[i] + log(tmp[i]);
    UNPROTECT(1);
    return TMP;
}


SEXP PML4(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP nco, SEXP contrast, SEXP N){
    int nr=INTEGER(NR)[0], nc=INTEGER(NC)[0], k=INTEGER(K)[0], i, j, indLL;
    int nTips = INTEGER(NTips)[0], *SC, *sc;
    double *g=REAL(G), *w=REAL(W), *tmp, *res;
    SEXP TMP;
    double *eva, *eve, *evei;
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));
    SC = (int *) R_alloc(nr * k, sizeof(int));
    sc = (int *) R_alloc(nr, sizeof(int));
    tmp = (double *) R_alloc(nr * k, sizeof(double));
    PROTECT(TMP = allocVector(REALSXP, nr));

    res=REAL(TMP);
    for(i=0; i<(k*nr); i++)tmp[i]=0.0;
    indLL = nr * nc * nTips;
    for(i=0; i<k; i++){
        lll3(dlist, eva, eve, evei, REAL(EL), g[i], &nr, &nc, INTEGER(node), INTEGER(edge), nTips, REAL(contrast), INTEGER(nco)[0], INTEGER(N)[0],  &SC[nr * i], REAL(bf), &tmp[i*nr], &LL[indLL *i], &SCM[nr * nTips * i]);
    }
    rowMinScale(SC, nr, k, sc);
    for(i=0; i<nr; i++){
        res[i]=0.0;
        for(j=0;j<k;j++)res[i] += w[j] * exp(LOG_SCALE_EPS * SC[i+j*nr]) * tmp[i+j*nr];
    }
    for(i=0; i<nr; i++) res[i] = log(res[i]) + LOG_SCALE_EPS * sc[i];
    UNPROTECT(1);
    return TMP;
}


void moveLL5(double *LL, double *child, double *P, int *nr, int *nc, double *tmp){
    int j;
    F77_CALL(dgemm)("N", "N", nr, nc, nc, &one, child, nr, P, nc, &zero, tmp, nr FCONE FCONE);
    for(j=0; j<(*nc * *nr); j++) LL[j]/=tmp[j];
    F77_CALL(dgemm)("N", "N", nr, nc, nc, &one, LL, nr, P, nc, &zero, tmp, nr FCONE FCONE);
    for(j=0; j<(*nc * *nr); j++) child[j] *= tmp[j];
}


// dad / child * P
void helpDADI(double *dad, double *child, double *P, int nr, int nc, double *res){
    F77_CALL(dgemm)("N", "N", &nr, &nc, &nc, &one, child, &nr, P, &nc, &zero, res, &nr FCONE FCONE);
    for(int j=0; j<(nc * nr); j++) dad[j]/=res[j];
}


// braucht Addition skalierte Werte
void helpPrep(double *dad, double *child, double *eve, double *evi, int nr, int nc, double *tmp, double *res){
    F77_CALL(dgemm)("N", "N", &nr, &nc, &nc, &one, child, &nr, eve, &nc, &zero, res, &nr FCONE FCONE);
    F77_CALL(dgemm)("N", "N", &nr, &nc, &nc, &one, dad, &nr, evi, &nc, &zero, tmp, &nr FCONE FCONE);
    for(int j=0; j<(nc * nr); j++) res[j]*=tmp[j];
}


void helpDAD2(double *dad, int *child, double *contrast, double *P, int nr, int nc, int nco, double *res){
    matp(child, contrast, P, &nr, &nc, &nco, res);
    for(int j=0; j<(nc * nr); j++) res[j]=dad[j]/res[j];
}


void helpDAD5(double *dad, int *child, double *contrast, double *P, int nr, int nc, int nco, double *res){
    matp(child, contrast, P, &nr, &nc, &nco, res);
    for(int j=0; j<(nc * nr); j++) dad[j]/=res[j];
}


SEXP getDAD2(SEXP dad, SEXP child, SEXP contrast, SEXP P, SEXP nr, SEXP nc, SEXP nco){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], nrs=INTEGER(nco)[0]; //, j
    SEXP TMP, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpDAD2(REAL(VECTOR_ELT(dad, i)), INTEGER(child), REAL(contrast), REAL(VECTOR_ELT(P, i)), nrx, ncx, nrs, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT
    return(RESULT);
}


void helpPrep2(double *dad, int *child, double *contrast, double *evi, int nr, int nc, int nrs, double *res){
    int i, j;
    F77_CALL(dgemm)("N", "N", &nr, &nc, &nc, &one, dad, &nr, evi, &nc, &zero, res, &nr FCONE FCONE);
    for(i = 0; i < nr; i++){
        for(j = 0; j < nc; j++) res[i + j*nr] *= contrast[child[i] - 1L + j*nrs];
    }
}


SEXP getPrep2(SEXP dad, SEXP child, SEXP contrast, SEXP evi, SEXP nr, SEXP nc, SEXP nco){
    R_len_t i, n=length(dad);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], ncs=INTEGER(nco)[0];
    SEXP TMP, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpPrep2(REAL(VECTOR_ELT(dad, i)), INTEGER(child), REAL(contrast),  REAL(evi), nrx, ncx, ncs, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1);
        }
    UNPROTECT(1);
    return(RESULT);
}


// child *= (dad * P)
void goDown(double *dad, double *child, double *P, int nr, int nc, double *res){
    F77_CALL(dgemm)("N", "N", &nr, &nc, &nc, &one, dad, &nr, P, &nc, &zero, res, &nr FCONE FCONE);
    for(int j=0; j<(nc * nr); j++) child[j]*=res[j];
}

// dad *= (child * P)
void goUp(double *dad, int *child, double *contrast, double *P, int nr, int nc, int nco, double *res){
    matp(child, contrast, P, &nr, &nc, &nco, res);
    for(int j=0; j<(nc * nr); j++) dad[j]*=res[j];
}

// , double *w
void updateLLQ(SEXP dlist, int pa, int ch, double *eva, double *eve, double*evei,
               double el, double *g, int nr,
               int nc, int ntips, double *contrast, int nco, int k,
               double *tmp, double *P){
    int i;
    if(ch>ntips){
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, el, g[i], P);
            goDown(&LL[LINDEX(ch, i)], &LL[LINDEX(pa, i)], P, nr, nc, tmp);
        }
    }
    else{
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, el, g[i], P);
            goUp(&LL[LINDEX(pa, i)], INTEGER(VECTOR_ELT(dlist, ch-1L)), contrast, P, nr, nc, nco, tmp);
        }
    }
}


// double *w,
void updateLL2(SEXP dlist, int pa, int ch, double *eva, double *eve, double*evei,
    double el,  double *g, int nr,
    int nc, int ntips, double *contrast, int nco, int k,
    double *tmp, double *P){
    int i;

    if(ch>ntips){
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, el, g[i], P);
            goDown(&LL[LINDEX(pa, i)], &LL[LINDEX(ch, i)], P, nr, nc, tmp);
         }
    }
    else{
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, el, g[i], P);
            goUp(&LL[LINDEX(pa, i)], INTEGER(VECTOR_ELT(dlist, ch-1L)), contrast, P, nr, nc, nco, tmp);
        }
    }
}


void ExtractScale(int ch, int k, int *nr, int *ntips, double *res){
    int i;
    int j, blub, tmp;
    for(i = 0; i < k; i++){
        blub = LINDEX3(ch, i);
        for(j=0; j< *nr; j++) res[j +i * *nr] = SCM[blub+j];
    }
    for(i = 0; i< *nr; i++){
        tmp = res[i];
        for(j = 1; j<k; j++){
            if(res[i+j * *nr]<tmp)tmp = res[i+j * *nr];
        }
        for(j=0; j<k; j++) res[i+j * *nr] = pow(ScaleEPS, (res[i+j * *nr] - tmp));
    }
}


// in getDad
// dad / child * P
void helpDAD(double *dad, double *child, double *P, int nr, int nc, double *res){
    F77_CALL(dgemm)("N", "N", &nr, &nc, &nc, &one, child, &nr, P, &nc, &zero, res, &nr FCONE FCONE);
    for(int j=0; j<(nc * nr); j++) res[j]=dad[j]/res[j];
}

//
SEXP getDAD(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0]; //, j
    SEXP TMP, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpDAD(REAL(VECTOR_ELT(dad, i)), REAL(VECTOR_ELT(child, i)), REAL(VECTOR_ELT(P, i)), nrx, ncx, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT
    return(RESULT);
    }


SEXP getPrep(SEXP dad, SEXP child, SEXP eve, SEXP evi, SEXP nr, SEXP nc){
    R_len_t i, n=length(dad);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0]; //, j
    double *tmp;
    SEXP TMP, RESULT;
    tmp = (double *) R_alloc(nrx*ncx, sizeof(double));
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        helpPrep(REAL(VECTOR_ELT(dad, i)), REAL(VECTOR_ELT(child, i)), REAL(eve),  REAL(evi), nrx, ncx, tmp, REAL(TMP));
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT
    return(RESULT);
}


void NR_d2f(double *eva, int nc, double el, double *w, double *g, double *X, int ld, int nr, double *f, double *res){
    int i, j, k;
    double *tmp;
    tmp = (double *) R_alloc(nc, sizeof(double));

    for(k=0; k<nr; k++) res[k] = 0.0;

    for(j=0;j<ld;j++){
        for(i=0; i<nc ;i++) tmp[i] = (eva[i] * g[j]) * exp(eva[i] * g[j] * el);
        F77_CALL(dgemv)("N", &nr, &nc, &w[j], &X[j*(nr*nc)], &nr, tmp, &ONE, &one, res, &ONE FCONE);
        }
    for(i=0; i<nr ;i++) res[i]/=f[i];
}


void NR_df(double *eva, int nc, double el, double *w, double *g, double *X, int ld, int nr, double *f, double *res){
    int i, j, k;
    double *tmp;
    tmp = (double *) R_alloc(nc, sizeof(double));
    for(k=0; k<nr; k++) res[k] = 0.0;
    int nrc = (nc+1L) * nr;
    for(j=0;j<ld;j++){
        for(i=0; i<nc ;i++) tmp[i] = (eva[i] * g[j]  * el) * exp(eva[i] * g[j] * el);
        F77_CALL(dgemv)("N", &nr, &nc, &w[j], &X[j*nrc], &nr, tmp, &ONE, &one, res, &ONE FCONE);
        }
    for(i=0; i<nr ;i++) res[i]/=f[i];
}


// tmp ausserhalb?
void NR_f(double *eva, int nc, double el, double *w, double *g, double *X, int ld, int nr, double *res){
    int i, j;
    double *tmp; //*res,  *dF,
    tmp = (double *) R_alloc(nc, sizeof(double));
    for(j=0;j<ld;j++){
        for(i=0; i<nc ;i++) tmp[i] = exp(eva[i] * g[j] * el);
        // alpha = w[j]
        F77_CALL(dgemv)("N", &nr, &nc, &w[j], &X[j*nr*nc], &nr, tmp, &ONE, &one, res, &ONE FCONE);
    }
}


// in ancestral.pml and rnodes
SEXP LogLik2(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco){
    R_len_t i, n = length(node);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], nt=INTEGER(nTips)[0], mn=INTEGER(mNodes)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node);
    SEXP ans, result;
    double *res, *rtmp;
    if(!isNewList(dlist)) error("'dlist' must be a list");
    ni = nodes[0];
    PROTECT(ans = allocVector(VECSXP, mn));
    PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
    res = REAL(result);
    rtmp = (double *) R_alloc(nrx*ncx, sizeof(double));
    for(j=0; j < (nrx * ncx); j++) res[j] = 1.0;
    for(i = 0; i < n; i++) {
        ei = edges[i];
        if(ni != nodes[i]){
            SET_VECTOR_ELT(ans, ni, result);
            UNPROTECT(1); //result
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            ni = nodes[i];
            if(ei < nt)
               matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), res);
            else
               F77_CALL(dgemm)("N", "N", &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx,
                   REAL(VECTOR_ELT(P, i)), &ncx, &zero, res, &nrx FCONE FCONE);
            }
        else {
            if(ei < nt)
                matp(INTEGER(VECTOR_ELT(dlist, ei)), REAL(contrast), REAL(VECTOR_ELT(P, i)), INTEGER(nr), INTEGER(nc), INTEGER(nco), rtmp);
            else
                F77_CALL(dgemm)("N", "N", &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(ans, ei-nt)), &nrx,
                    REAL(VECTOR_ELT(P, i)), &ncx, &zero, rtmp, &nrx FCONE FCONE);
            for(j=0; j < (nrx*ncx); j++) res[j] *= rtmp[j];
        }
    }
    SET_VECTOR_ELT(ans, ni, result);
    UNPROTECT(2); // result ans
    return(ans);
}

// raus
//void matprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z)
//{
//    F77_CALL(dgemm)("N", "N", &nrx, &ncy, &ncx, &one, x, &nrx, y, &nry, &zero, z, &nrx FCONE FCONE);
//}


SEXP getM3(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc){
    R_len_t i, n=length(P);
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0], j;
    SEXP TMP, RESULT;
    double *tmp, *daddy;
    PROTECT(RESULT = allocVector(VECSXP, n));
    for(i=0; i<n; i++){
        PROTECT(TMP = allocMatrix(REALSXP, nrx, ncx));
        tmp = REAL(TMP);
//        matprod(REAL(VECTOR_ELT(child, i)), nrx, ncx, REAL(VECTOR_ELT(P, i)), ncx, ncx, tmp);
        F77_CALL(dgemm)("N", "N", &nrx, &ncx, &ncx, &one, REAL(VECTOR_ELT(child, i)), &nrx, REAL(VECTOR_ELT(P, i)), &ncx, &zero, tmp, &nrx FCONE FCONE);
        daddy = REAL(VECTOR_ELT(dad, i));
        for(j=0; j<(ncx * nrx); j++) tmp[j]*=daddy[j];
        SET_VECTOR_ELT(RESULT, i, TMP);
        UNPROTECT(1); // TMP
        }
    UNPROTECT(1); //RESULT
    return(RESULT);
    }



void fs3(double *eva, int nc, double el, double *w, double *g, double *X, int ld, int nr, double *weight,
         double *f0, double tau, double *res) //, int mkv)
{
    int mkv=0;
    double *tmp, *f, edle, ledle, newedle, eps=10;
    double ll=0.0, lll, delta=0.0, scalep = 1.0, l1=0.0, l0=0.0;
    double y, p0=0.0; //sum_wgt=0.0,
    int i, j, k=0;
    tmp = (double *) R_alloc(nr, sizeof(double));
    f = (double *) R_alloc(nr, sizeof(double));
    edle = el; //REAL(el)[0];

    for(i=0; i<nr; i++)f[i] = f0[i];
    NR_f(eva, nc, edle, w, g, X, ld, nr, f); // nc-1L !!!
    if(mkv){
      for(j=0; j<nc; j++)p0 += f[i];
      for(i=0; i<nr ;i++) f[i] /= (1-p0);
    }
    for(i=0; i<nr ;i++){
      l0 += weight[i] * log(f[i]);
      //sum_wgt += weight[i];
    }

    while ( (eps > 1e-05) &&  (k < 10) ) {
        if(scalep>0.6){
            NR_df(eva, nc-1L, edle, w, g, X, ld, nr, f, tmp);
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
        // some error handling avoid too big small edges & too big steps
        if (newedle > 10.0) newedle = 10.0;
        if (newedle < tau) newedle = tau; // 1e-8 phyML

        for(i=0; i<nr; i++)f[i] = f0[i];
        NR_f(eva, nc, newedle, w, g, X, ld, nr, f);
        if(mkv){
          for(j=0; j<nc; j++)p0 += f[i];
          for(i=0; i<nr ;i++) f[i] /= (1-p0);
        }
        l1 = 0.0;
        for(i=0; i<nr ;i++) l1 += weight[i] * log(f[i]);
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
    NR_d2f(eva, nc-1L, edle, w, g, X, ld, nr, f, tmp);
    lll=0.0;
    for(i=0; i<nr ;i++) lll+=weight[i]*tmp[i]*tmp[i];
    res[0] = edle;
    res[1] = 1.0 / lll;
    res[2] = l0; //l0
    // return res[3] = 1 or 0 for success
}


// in fs()
SEXP FS4(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP dad, SEXP child, SEXP ld, SEXP nr,
         SEXP weight, SEXP f0, SEXP tau, SEXP retA, SEXP retB)
{
    SEXP RESULT, EL, P;
    double *tmp, *wgt=REAL(weight), edle, *eva=REAL(VECTOR_ELT(eig,0));
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
    tmp = (double *) R_alloc(3L, sizeof(double));
    PROTECT(RESULT = allocVector(VECSXP, 4));
    edle = REAL(el)[0];
    fs3(eva, ncx, edle, REAL(w), REAL(g), REAL(X), INTEGER(ld)[0], nrx, wgt, REAL(f0), REAL(tau)[0], tmp);
    PROTECT(EL = ScalarReal(tmp[0]));
    PROTECT(P = getPM(eig, nc, EL, g));
    SET_VECTOR_ELT(RESULT, 0, EL);
    if(INTEGER(retA)[0]>0L)SET_VECTOR_ELT(RESULT, 1, getM3(child, dad, P, nr, nc));
    if(INTEGER(retB)[0]>0L)SET_VECTOR_ELT(RESULT, 2, getM3(dad, child, P, nr, nc));
// add variance ??
    SET_VECTOR_ELT(RESULT, 3, ScalarReal(tmp[2]));
    UNPROTECT(3);
    return (RESULT);
}

// in dist.ml
SEXP FS5(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP ld, SEXP nr, SEXP weight, SEXP f0, SEXP tau)
{
    SEXP RESULT; // EL, P;
    double *wgt=REAL(weight), *eva=REAL(VECTOR_ELT(eig,0));
    int ncx=INTEGER(nc)[0], nrx=INTEGER(nr)[0];
    PROTECT(RESULT = allocVector(REALSXP, 3));
    double edle = REAL(el)[0];
    fs3(eva, ncx, edle, REAL(w), REAL(g), REAL(X), INTEGER(ld)[0], nrx, wgt, REAL(f0), REAL(tau)[0], REAL(RESULT));
    UNPROTECT(1);
    return (RESULT);
}


SEXP optE(SEXP PARENT, SEXP CHILD, SEXP ANC, SEXP eig, SEXP EVI, SEXP EL,
                  SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP NTIPS, SEXP CONTRAST,
                  SEXP CONTRAST2, SEXP NCO,
                  SEXP dlist, SEXP WEIGHT, SEXP F0, SEXP TAU){
    int i, k=length(W), h, j, n=length(PARENT), m, lEL=length(EL);
    int nc=INTEGER(NC)[0], nr=INTEGER(NR)[0], ntips=INTEGER(NTIPS)[0];
    int *parent=INTEGER(PARENT), *child=INTEGER(CHILD), *anc=INTEGER(ANC);
    int loli, nco =INTEGER(NCO)[0];
    double *weight=REAL(WEIGHT), *f0=REAL(F0), *w=REAL(W), tau=REAL(TAU)[0];
    double *g=REAL(G), *evi=REAL(EVI), *contrast=REAL(CONTRAST), *contrast2=REAL(CONTRAST2);
    double *el; //=REAL(EL);
    double *eva, *eve, *evei, *tmp, *P;
    double  *X; // define it *blub=REAL(BLUB),
    double *blub = (double *) R_alloc(nr * k, sizeof(double));
    double oldel; //=el[ch-1L]
    int ancloli, pa, ch; //=anc[loli]
    double *res = (double *) R_alloc(3L, sizeof(double));
    tmp = (double *) R_alloc(nr * nc, sizeof(double));
    P = (double *) R_alloc(nc * nc, sizeof(double));
    X = (double *) R_alloc(k * nr * nc, sizeof(double));

    ExtractScale(parent[0], k, &nr, &ntips, blub);

    SEXP RESULT;
    PROTECT(RESULT = allocVector(REALSXP, lEL));
    el=REAL(RESULT);
    for(i = 0; i < lEL; i++) el[i] = REAL(EL)[i];
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));

    loli = parent[0];

    for(m = 0; m < n; m++){
        pa = parent[m];
        ch = child[m];
        oldel=el[ch-1L];

    while(loli != pa){
        ancloli=anc[loli];
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, el[loli-1L], g[i], P);
            moveLL5(&LL[LINDEX(loli, i)], &LL[LINDEX(ancloli, i)], P, &nr, &nc, tmp);
        }
        loli = ancloli;
    }
    // moveDad
    if(ch>ntips){
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, oldel, g[i], P);
            helpDADI(&LL[LINDEX(pa, i)], &LL[LINDEX(ch, i)], P, nr, nc, tmp);
            helpPrep(&LL[LINDEX(pa, i)], &LL[LINDEX(ch, i)], eve, evi, nr, nc, tmp, &X[i*nr*nc]);
            for(h = 0; h < nc; h++){
                for(j = 0; j < nr; j++){
                    X[j+h*nr + i*nr*nc] *= blub[j+i*nr];
                }
            }
        }
    }
    else{
        for(i = 0; i < k; i++){
            getP(eva, eve, evei, nc, oldel, g[i], P);
            helpDAD5(&LL[LINDEX(pa, i)], INTEGER(VECTOR_ELT(dlist, ch-1L)), contrast, P, nr, nc, nco, tmp);
            helpPrep2(&LL[LINDEX(pa, i)], INTEGER(VECTOR_ELT(dlist, ch-1L)), contrast2, evi, nr, nc, nco, &X[i*nr*nc]); //;
            for(h = 0; h < nc; h++){
                for(j = 0; j < nr; j++){
                    X[j+h*nr + i*nr*nc] *= blub[j+i*nr];
                }
            }
        }
    }
    fs3(eva, nc, oldel, w, g, X, k, nr, weight, f0, tau, res);
    updateLL2(dlist, pa, ch, eva, eve, evei, res[0], g, nr,
        nc, ntips, contrast, nco, k, tmp, P);
        el[ch-1L] = res[0];
        if (ch > ntips) loli  = ch;
        else loli = pa;
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
}


SEXP optQrtt(SEXP PARENT, SEXP CHILD, SEXP eig, SEXP EVI, SEXP EL,
          SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP NTIPS, SEXP CONTRAST,
          SEXP CONTRAST2, SEXP NCO,
          SEXP dlist, SEXP WEIGHT, SEXP F0, SEXP TAU){
    int i, k=length(W), h, j, m, lEL=length(EL);
    int nc=INTEGER(NC)[0], nr=INTEGER(NR)[0], ntips=INTEGER(NTIPS)[0];
    int *parent=INTEGER(PARENT), *child=INTEGER(CHILD), pa, ch;
    int nco =INTEGER(NCO)[0]; // loli,
    double *weight=REAL(WEIGHT), *f0=REAL(F0), *w=REAL(W), tau=REAL(TAU)[0];
    double *g=REAL(G), *evi=REAL(EVI), *contrast=REAL(CONTRAST), *contrast2=REAL(CONTRAST2);
    double *el, *X;
    double *eva, *eve, *evei, *tmp, *P;
    double *blub = (double *) R_alloc(nr * k, sizeof(double));
    double oldel;
    double *res = (double *) R_alloc(3L, sizeof(double));
    tmp = (double *) R_alloc(nr * nc, sizeof(double));
    P = (double *) R_alloc(nc * nc, sizeof(double));
    X = (double *) R_alloc(k * nr * nc, sizeof(double));

    ExtractScale(parent[0], k, &nr, &ntips, blub);

    SEXP RESULT;
    PROTECT(RESULT = allocVector(REALSXP, lEL));
    el=REAL(RESULT);
    for(i = 0; i < lEL; i++) el[i] = REAL(EL)[i];
    eva = REAL(VECTOR_ELT(eig, 0));
    eve = REAL(VECTOR_ELT(eig, 1));
    evei = REAL(VECTOR_ELT(eig, 2));

    for(m = 4L; m > -1L; m--){
        pa = parent[m];
        ch = child[m];
        oldel=el[m];
        // moveDad
        if(ch>ntips){
            for(i = 0; i < k; i++){
                getP(eva, eve, evei, nc, oldel, g[i], P);
                helpDADI(&LL[LINDEX(pa, i)], &LL[LINDEX(ch, i)], P, nr, nc, tmp);
                helpPrep(&LL[LINDEX(pa, i)], &LL[LINDEX(ch, i)], eve, evi, nr, nc, tmp, &X[i*nr*nc]);
                for(h = 0; h < nc; h++){
                    for(j = 0; j < nr; j++){
                        X[j+h*nr + i*nr*nc] *= blub[j+i*nr];
                    }
                }
            }
        }
        else{
            for(i = 0; i < k; i++){
                getP(eva, eve, evei, nc, oldel, g[i], P);
                helpDAD5(&LL[LINDEX(pa, i)], INTEGER(VECTOR_ELT(dlist, ch-1L)), contrast, P, nr, nc, nco, tmp);
                helpPrep2(&LL[LINDEX(pa, i)], INTEGER(VECTOR_ELT(dlist, ch-1L)), contrast2, evi, nr, nc, nco, &X[i*nr*nc]); //;
                for(h = 0; h < nc; h++){
                    for(j = 0; j < nr; j++){
                        X[j+h*nr + i*nr*nc] *= blub[j+i*nr];
                    }
                }
            }
        }
        fs3(eva, nc, oldel, w, g, X, k, nr, weight, f0, tau, res);
// go up
// if i=2 go down
        if(m==2)updateLLQ(dlist, ch, pa, eva, eve, evei, res[0], g, nr,
                  nc, ntips, contrast, nco, k, tmp, P);
        else updateLLQ(dlist, pa, ch, eva, eve, evei, res[0], g, nr,
                  nc, ntips, contrast, nco, k, tmp, P);
        el[m] = res[0];
    }
    UNPROTECT(1); //RESULT
    return(RESULT);
}
