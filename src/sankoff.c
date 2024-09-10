/*
 * dist.c
 *
 * (c) 2008-2024  Klaus Schliep (klaus.schliep@gmail.com)
 *
 *
 * This code may be distributed under the GNU GPL
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>



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

/*
void rowMin2(double *dat, int n,  int k, double *res){
    int i, h;
    double x;
    for(i = 0; i < n; i++){
        x = dat[i];
        for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
        res[i] = x;
    }
}
*/

double get_ps(double *dat, int n,  int k, double *weight){
  int i, h;
  double x, y=0.0;
  for(i = 0; i < n; i++){
    x = dat[i];
    for(h = 1; h< k; h++) {if(dat[i + h*n] < x) x=dat[i + h*n];}
    y += x*weight[i];
  }
  return(y);
}


void sankoffNode(double *dat, int n, double *cost, int k, double *result){
    int i, j, h;
    double x, tmp;
    for(i = 0; i < n; i++){
        for(j = 0; j < k; j++){
            x = dat[i] + cost[j*k];
            for(h = 1; h< k; h++){
                tmp = dat[i + h*n] + cost[h + j*k];
                if(tmp<x)x = tmp;
            }
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


// switch if else < nTip
double sankoffQuartet_new(SEXP dat, int n, double *cost, int k, double *weight,
                        int pos1, int pos2, int pos3, int pos4, int nTip,
                        double *contr, int nrc){  // SEXP dat2,
  double *res, *rtmp, erg=0.0;
  rtmp = (double *) R_alloc(n*k, sizeof(double));
  res = (double *) R_alloc(n*k, sizeof(double));
  for(int j=0; j<(n*k); j++) rtmp[j] = 0.0;
  for(int j=0; j<(n*k); j++) res[j] = 0.0;
//  if(pos1 < nTip) sankoffTips(INTEGER(VECTOR_ELT(dat, pos1)), contr, n, k, nrc, rtmp);
//  else
    sankoffNode(REAL(VECTOR_ELT(dat,pos1)), n, cost, k, rtmp);
//  if(pos2 < nTip) sankoffTips(INTEGER(VECTOR_ELT(dat,pos2)), contr, n, k, nrc, rtmp);
//  else
    sankoffNode(REAL(VECTOR_ELT(dat,pos2)), n, cost, k, rtmp);
  sankoffNode(rtmp, n, cost, k, res);
//  if(pos3 < nTip) sankoffTips(INTEGER(VECTOR_ELT(dat,pos3)), contr, n, k, nrc, res);
//  else
    sankoffNode(REAL(VECTOR_ELT(dat,pos3)), n, cost, k, res);
//  if(pos4 < nTip) sankoffTips(INTEGER(VECTOR_ELT(dat,pos4)), contr, n, k, nrc, res);
//  else
    sankoffNode(REAL(VECTOR_ELT(dat,pos4)), n, cost, k, res);
  erg = get_ps(res, n, k, weight);
  return(erg);
}


SEXP sankoff_nni_c(SEXP dat, SEXP sn, SEXP scost, SEXP sk, SEXP sweight,
                   SEXP POS, SEXP snpos, SEXP ntip, SEXP contrasts, SEXP snrc){ // SEXP dat2,
    int n=INTEGER(sn)[0], k=INTEGER(sk)[0], npos=INTEGER(snpos)[0];
    int nTip=INTEGER(ntip)[0], nrc=INTEGER(snrc)[0];
    double *cost, *res, *weight, *contr;
    int *pos;
    SEXP result;
    PROTECT(result = allocMatrix(REALSXP, npos, 2L));
    res = REAL(result);
    cost = REAL(scost);
    contr = REAL(contrasts);
    weight = REAL(sweight);
    pos = INTEGER(POS);
    for(int i=0; i<npos; ++i){
      res[i] = sankoffQuartet_new(dat, n, cost, k, weight, pos[i], pos[i+2*npos],
                                  pos[i+npos], pos[i+3*npos], nTip, contr, nrc);
      res[i+npos] = sankoffQuartet_new(dat, n, cost, k, weight, pos[i+npos],
                                  pos[i+2*npos], pos[i], pos[i+3*npos], nTip, contr, nrc);
    }
    UNPROTECT(1);
    return result;
}


// sankoffNew
SEXP sankoff_c(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips, SEXP contrast, SEXP nrs){
    R_len_t i, n = length(node); //, nt = length(tips);
    int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], mn=INTEGER(mNodes)[0], nrc = INTEGER(nrs)[0];
    int  ni, ei, j, *edges=INTEGER(edge), *nodes=INTEGER(node), ntips=INTEGER(tips)[0];
    SEXP result, dlist2; //tmp,
    double *res, *cost, *tmp; // *rtmp,
    tmp = (double *) R_alloc(ncx * nrc, sizeof(double));
    for(j=0; j<(ncx * nrc); j++) tmp[j] = 0.0;
    cost = REAL(scost);
    sankoffNode(REAL(contrast), nrc, cost, ncx, tmp);

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
            else sankoffNode(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res);
            }
        else{
            SET_VECTOR_ELT(dlist2, ni, result);
            UNPROTECT(1);
            PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
            res = REAL(result);
            for(j=0; j<(nrx * ncx); j++) res[j] = 0.0;
            ni = nodes[i];
            if(ei < ntips) sankoffTips(INTEGER(VECTOR_ELT(dlist,ei)), tmp, nrx, ncx, nrc, res);
            else sankoffNode(REAL(VECTOR_ELT(dlist2,ei)), nrx, cost, ncx, res);
            }
    }
    SET_VECTOR_ELT(dlist2, ni, result);
    UNPROTECT(2);
    return(dlist2);
}


SEXP sankoffMPR(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nnode){
  R_len_t i, n = length(node);
  int nrx=INTEGER(nr)[0], ncx=INTEGER(nc)[0], n0, nn=INTEGER(nnode)[0];
  int  ei, j, *nodes=INTEGER(node), *edges=INTEGER(edge);
  SEXP result, dlist2; //tmp,
  double *res, *cost; // *rtmp,
  cost = REAL(scost);
  n0 = nodes[n-1L];
  PROTECT(dlist2 = allocVector(VECSXP, n+1L));
  PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
  res = REAL(result);
  for(j=0;j<(nrx*ncx);j++)res[j]=0.0;
  for(j=n-1L; j>=0; j--) {
    if(nodes[j]!=n0){
      SET_VECTOR_ELT(dlist2, n0, result);
      UNPROTECT(1);
      n0 = nodes[j];
      PROTECT(result = allocMatrix(REALSXP, nrx, ncx));
      res = REAL(result);
      for(i=0; i<(nrx * ncx); i++) res[i] = 0.0;
      sankoffNode(REAL(VECTOR_ELT(dlist,nodes[j]+nn)), nrx, cost, ncx, res); //nodes[j] + nnode
    }
    ei = edges[j];
    sankoffNode(REAL(VECTOR_ELT(dlist,ei)), nrx, cost, ncx, res);
  }
  SET_VECTOR_ELT(dlist2, n0, result);
  UNPROTECT(2);
  return(dlist2);
}

