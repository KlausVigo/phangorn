#define USE_RINTERNALS

#include <Rmath.h>
#include <math.h>
#include <R.h> 
#include <Rinternals.h>


// use R_len_t stat int, e.g. nr
double huge = 1.0e300;
static int *data1, *data2;
static double *weight;


void fitch_free(){
    free(data1);
    free(data2);
    free(weight);
}

// type of fitch depending on nc e.g. int, long generic C++
void fitch_init(int *data, int *m, int *n, double *weights, int *nr)
{
    int i;
    data1 = (int *) calloc(*n, sizeof(int));
    data2 = (int *) calloc(*n, sizeof(int));  
    weight = (double *) calloc(*nr, sizeof(double));   
    for(i=0; i<*m; i++) data1[i] = data[i];  
    for(i=0; i<*nr; i++) weight[i] = weights[i];
}


SEXP getData(SEXP n, SEXP k){
    int i, m=INTEGER(n)[0], l=INTEGER(k)[0];  
    SEXP DAT, DAT2, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, 2L));
    PROTECT(DAT = allocMatrix(INTSXP, m, l));
    PROTECT(DAT2 = allocMatrix(INTSXP, m, l)); 
    for(i=0; i< m*l; i++) INTEGER(DAT)[i] = data1[i];
    for(i=0; i< m*l; i++) INTEGER(DAT2)[i] = data2[i];
    SET_VECTOR_ELT(RESULT, 0, DAT);
    SET_VECTOR_ELT(RESULT, 1, DAT2);
    UNPROTECT(3);
    return(RESULT); 
}


SEXP getWeight(SEXP n){
    int i, m=INTEGER(n)[0];  
    SEXP RESULT;
    PROTECT(RESULT = allocVector(REALSXP, m));
    for(i=0; i<m; i++) REAL(RESULT)[i] = weight[i];
    UNPROTECT(1);
    return(RESULT); 
}


int bitcount(int x){ 
    int count;
    for (count=0; x != 0; x>>=1)
       if ( x & 01)
           count++;
    return count;
}


void bitCount(int *x, int *count){
    count[0]=bitcount(x[0]);
} 


void addOne(int *edge, int *tip, int *ind, int *l, int *m, int *result){
    int add = 1L, j=0L, p, k, i, l2=*l+2L, ei;
    p = edge[*ind-1L];
    k = edge[*ind-1L + *l];
    for(i=0; i<*l; i++){
        ei = edge[i]; 
        if( (add==1L) && (ei==p) ){
            result[j] = *m;
            result[j+l2] = k;
            j++;  
            result[j] = *m;
            result[j+l2] = *tip;
            j++;
            add=0L;
        }
        if(i== (*ind-1L)) result[j+l2] = *m;
        else result[j+l2] = edge[i+ *l];
        result[j] = edge[i];       
        j++; 
    }
}


SEXP AddOne(SEXP edge, SEXP tip, SEXP ind, SEXP l, SEXP m){
    SEXP result;
    PROTECT(result = allocMatrix(INTSXP, INTEGER(l)[0]+2L, 2L));
    addOne(INTEGER(edge), INTEGER(tip), INTEGER(ind), INTEGER(l), INTEGER(m), INTEGER(result));
    UNPROTECT(1);
    return(result);
}


void fitch43(int *dat1, int *dat2, int *nr, int *pars, double *weight, double *w){
    int k, tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(!tmp){
            tmp = dat1[k] | dat2[k];
            (pars[k])++;
            (*w)+=weight[k];
        }
        dat1[k] = tmp;
    } 
}


void fitch44(int *res, int *dat1, int *dat2, int *nr, int *pars, double *weight, double *w){
    int k, tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(!tmp){
            tmp = dat1[k] | dat2[k];
            (pars[k])++;
            (*w)+=weight[k];
        }
        res[k] = tmp;
    } 
}


void fitch53(int *dat1, int *dat2, int *nr, double *weight, double *w){
    int k, tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(!tmp){
            tmp = dat1[k] | dat2[k];
            (*w)+=weight[k];
        }
        dat1[k] = tmp;
    } 
}


void fitch54(int *res, int *dat1, int *dat2, int *nr, double *weight, double *w){
    int k, tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(!tmp){
            tmp = dat1[k] | dat2[k];
            (*w)+=weight[k];
        }
        res[k] = tmp;
    } 
}


SEXP FITCHTRIP3(SEXP DAT3, SEXP nrx, SEXP edge, SEXP score, SEXP PS){ 
    R_len_t i, m = length(edge);  
    int nr=INTEGER(nrx)[0], k, tmp, ei, *edges=INTEGER(edge); 
    int d3=INTEGER(DAT3)[0] - 1;
    double *pvtmp;  
    double ps = REAL(PS)[0];
    SEXP pvec;
    PROTECT(pvec = allocVector(REALSXP, m));
    pvtmp = REAL(pvec);
    for(i=0; i<m; i++) pvtmp[i] = REAL(score)[i]; 
#ifdef SUPPORT_OPENMP    
#pragma omp parallel for private(i, ei, k, tmp) shared(edges, data1, data2, d3, nr, weight, ps, pvtmp)
#endif
    for(i=0; i<m; i++){
        ei = edges[i] - 1L;
//      pvtmp[i] = REAL(score)[ei];
        for(k = 0; k < nr; k++){
            tmp = data1[k + ei*nr] & data2[k + ei*nr];
            if(!tmp){
                tmp = data1[k + ei*nr] | data2[k + ei*nr];
                pvtmp[i]+=weight[k];
                
            }
            tmp = tmp & data1[k + d3*nr];
            if(!tmp){
               pvtmp[i]+=weight[k];                
            }
            if(pvtmp[i]>ps)break;
        }
//        if(pvtmp[i]<ps) ps = pvtmp[i] + 1.0e-8 random.addition order
    }
    UNPROTECT(1);
    return(pvec); 
}



void fitch8(int *dat, int *nr, int *pars, int *node, int *edge, int *nl, double *weight, double *pvec, double *pscore) 
{   
    int i, ni=0L, ri, le;
    i=0L;
    while(i<(*nl - 1L)){
        ni = node[i] - 1L; 
        le = edge[i] - 1L;
        ri = edge[i+1] - 1L; 
        pvec[ni] = pvec[le] + pvec[ri];
	fitch44(&dat[ni * (*nr)], &dat[le * (*nr)], &dat[ri * (*nr)], nr, pars, weight, &pvec[ni]); 
        i++;
        i++;                  
    }
    if(i == (*nl-1L)){
        le = edge[i] - 1L;
        pvec[ni] += pvec[le]; 
        fitch43(&dat[ni * (*nr)], &dat[le * (*nr)], nr, pars, weight, &pvec[ni]); 
    } 
    pscore[0]=pvec[ni];
}


void fitch9(int *dat, int *nr, int *node, int *edge, int *nl, double *weight, double *pvec, double *pscore) 
{   
    int i, ni=0L, ri, le;
    i=0L;
    while(i<(*nl - 1L)){
        ni = node[i] - 1L; 
        le = edge[i] - 1L;
        ri = edge[i+1] - 1L; 
        pvec[ni] = pvec[le] + pvec[ri];
	fitch54(&dat[ni * (*nr)], &dat[le * (*nr)], &dat[ri * (*nr)], nr, weight, &pvec[ni]); 
        i++;
        i++;                  
    }
    if(i == (*nl-1L)){
        le = edge[i] - 1L;
        pvec[ni] += pvec[le]; 
        fitch53(&dat[ni * (*nr)], &dat[le * (*nr)], nr, weight, &pvec[ni]); 
    } 
    pscore[0]=pvec[ni];
}


// in fitch
SEXP FITCH(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP weight, SEXP mx, SEXP q){   
    int *data, *nr=INTEGER(nrx), m=INTEGER(mx)[0], i, n=INTEGER(q)[0];   
    double *pvtmp;  
    SEXP DAT, pars, pvec, pscore, RESULT;
    PROTECT(RESULT = allocVector(VECSXP, 4L));
    PROTECT(pars = allocVector(INTSXP, *nr));
    PROTECT(pscore = allocVector(REALSXP, 1L));
    PROTECT(DAT = allocMatrix(INTSXP, nr[0], m));
    PROTECT(pvec = allocVector(REALSXP, m));
    pvtmp = REAL(pvec);
    data = INTEGER(DAT);
    for(i=0; i<m; i++) pvtmp[i] = 0.0;
    for(i=0; i<*nr; i++) INTEGER(pars)[i] = 0L;
    REAL(pscore)[0]=0.0;
    for(i=0; i<(*nr * n); i++)data[i] = INTEGER(dat)[i];
    
    fitch8(data, nr, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), REAL(weight), pvtmp, REAL(pscore));
    
    SET_VECTOR_ELT(RESULT, 0, pscore);
    SET_VECTOR_ELT(RESULT, 1, pars);
    SET_VECTOR_ELT(RESULT, 2, DAT);
    SET_VECTOR_ELT(RESULT, 3, pvec);
    UNPROTECT(5);
    return(RESULT); 
}


/*
ACCTRAN
*/
void fitchT(int *dat1, int *dat2, int *nr, double *pars, double *weight, double *w){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(tmp > 0L){
             dat1[k] = tmp;
             }
    } 
}


void fitchT3(int *dat1, int *dat2, int *nr, double *pars, double *weight, double *w){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
       tmp = dat1[k] & dat2[k];
       if(tmp==0L) {
             (*w)+=weight[k];
             pars[k] += 1;
             }
       if(tmp >0){
           if(tmp < dat2[k]){ 
              (*w)+= .5*weight[k];
              pars[k] += .5;
           }
       }

    } 
}


// return lower and upper bound for the number of changes 
// upper bound very conservative 
void countMPR(double *res, int *dat1, int *dat2, int *nr, double *weight, int *external){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];

        if(tmp==0){
            res[0] += weight[k];
            res[1] += weight[k];
        }
        else{ 
            if( external[0]==0L){ 
                 if( bitcount(dat1[k] | dat2[k])>1L ) res[1] += weight[k]; // dat1[k] != dat2[k]
            }   
            else{ 
                 if( tmp  < dat2[k] ) res[1] += weight[k];
            }
        }
    } 
}


void ACCTRAN2(int *dat, int *nr, double *pars, int *node, int *edge, int *nl, double *weight, double *pvec, int *nTips) 
{   
    int i;
    for (i=0; i< *nl; i++) {       
        if(edge[i]>nTips[0]) fitchT(&dat[(edge[i]-1L) * (*nr)], &dat[(node[i]-1) * (*nr)], nr, pars, weight, &pvec[i]); 
        }
}


void ACCTRAN3(int *dat, int *nr, double *pars, int *node, int *edge, int *nl, double *weight, double *pvec, int *nTips) 
{   
    int i;
    for (i=0; i< *nr; i++)pars[i]=0.0;
    for(i=0; i< *nl; i++)pvec[i] = 0.0;
    for (i=0; i< *nl; i++) {               
        fitchT3(&dat[(edge[i]-1L) * (*nr)], &dat[(node[i]-1) * (*nr)], nr, pars, weight, &pvec[i]); 
    }            
}


void fitchNNN(int d1, int d2){
    int tmp;
    tmp = d1 & d2;
    if(tmp) d1 = tmp;
    else d1 = d1 | d2;
}
// haeufig 0
void fitchTripletNew(int *res, int *dat1, int *dat2, int *dat3, int *nr) 
{   
    int k, v1, v2, v3;

    for(k = 0; k < (*nr); k++){
    v1 = dat1[k];
    fitchNNN(v1, dat2[k]);
    fitchNNN(v1, dat3[k]);

    v2 = dat1[k];
    fitchNNN(v2, dat3[k]);
    fitchNNN(v2, dat2[k]);

    v3 = dat2[k];
    fitchNNN(v3, dat3[k]);
    fitchNNN(v3, dat1[k]);

    res[k] = v1 & v2; // &v3[k];  
    res[k] = res[k] & v3; 
    }
}

void fitchN(int *dat1, int *dat2, int *nr){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
        tmp = dat1[k] & dat2[k];
        if(tmp) dat1[k] = tmp;
        else dat1[k] = dat1[k] | dat2[k];
    } 
}

// raus
void fitchN2(int *res, int *dat, int *node, int *edge, int *nr, int *nl) { 
    int i;
    for (i=0; i< *nl; i++) {
        fitchN(&res[(node[i]-1L) * (*nr)], &dat[(edge[i]-1L) * (*nr)], nr);              
    }
}

// MPR reconstruction nicht immer gleiches ergebnis
void fitchTriplet(int *res, int *dat1, int *dat2, int *dat3, int *nr) 
{   
    int k; // ni,
//    ni = 0;
    
    int *v1, *v2, *v3;
    v1 = (int *) R_alloc(*nr, sizeof(int));    
    v2 = (int *) R_alloc(*nr, sizeof(int));
    v3 = (int *) R_alloc(*nr, sizeof(int));

    for(k = 0; k < (*nr); k++) v1[k] = dat1[k];
    fitchN(v1, dat2, nr);
    fitchN(v1, dat3, nr);

    for(k = 0; k < (*nr); k++) v2[k] = dat1[k];
    fitchN(v2, dat3, nr);
    fitchN(v2, dat2, nr);

    for(k = 0; k < (*nr); k++) v3[k] = dat2[k];
    fitchN(v3, dat3, nr);
    fitchN(v3, dat1, nr);

    for(k = 0; k < (*nr); k++)res[k] = v1[k] & v2[k]; // &v3[k];  
    for(k = 0; k < (*nr); k++)res[k] = res[k] & v3[k];  
}


void prepRooted(int *res, int *nr, int *kids){ //int *data1, 
    fitchTriplet(res, &data1[*nr * (kids[0]-1L)], &data1[*nr * (kids[1]-1L)],  
        &data1[*nr * (kids[2]-1L)], nr);
}


void C_MPR(int *res, int *nr, int *parent, int *kids, int *nl) { 
    int p, k1, k2;
    int i = *nl -1;    
    while (i > 0L) {
        p = parent[i] - 1L;
        k1 = kids[i] - 1L;
        k2 = kids[i-1L] - 1L;
        fitchTriplet(&res[*nr * p], &data1[*nr* (k1)], &data1[*nr* (k2) ], &data2[*nr * p], nr);
        i -= 2L;
    }        
}


void fitchNACC2(int *root, int *dat, int *nr, double *pars, int *result, double *weight, double *pars1){
    int k;
    int tmp;
    for(k = 0; k < (*nr); k++){
//       result[k] = 0L;
       tmp = root[k] & dat[k];
       if(tmp==0L) {
             pars[0] += weight[k];
             pars1[k] += weight[k];
             }
       if(tmp >0){
           if(tmp < root[k]){ 
              pars[0] += .5*weight[k];
              pars1[k] += .5*weight[k];
              result[k] += 1L;
           }
       }
    }        
}


void fitchTripletACC4(int *root, int *dat1, int *dat2, int *dat3, int *nr, double *p1, double *p2, double *p3, double *weight, double *pars1, int *v1) 
{   
    int k;
       
    int tmp, a, b, c, t1, t2, t3;
    double d, f;
    for(k = 0; k < (*nr); k++){
        tmp = root[k];
        a = dat1[k] & dat2[k]; 
        b = dat1[k] & dat3[k];
        c = dat2[k] & dat3[k];
        if((a+b+c) == 0L){
           d = (2.0/3.0) * weight[k];
           p1[0] += d; 
           p2[0] += d;
           p3[0] += d;        
           pars1[k] += 2*weight[k]; 
           v1[k] = 2L; 
        }
        else{  
            f = 0.0;
            d = weight[k];
            t1 = 0.0;
            t2 = 0.0;
            t3 = 0.0;
            if( (dat1[k] & tmp)<tmp){ 
                t1 = d; 
                f+=1.0;
            }
            if( (dat2[k] & tmp)<tmp){ 
                t2 = d; 
                f+=1.0;
            }
            if( (dat3[k] & tmp)<tmp){ 
                t3 = d; 
                f+=1.0;
            }
            if(f>0.0){   
                pars1[k] += weight[k]; 
                p1[0] += t1/f; 
                p2[0] += t2/f;
                p3[0] += t3/f;
                v1[k] += 1L;
            }
        }
    }
}



SEXP FITCH345(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP ps){   
    int *nr=INTEGER(nrx), m=INTEGER(mx)[0], i;  
    double *pvtmp;  
    SEXP pars, pscore; 
    PROTECT(pars = allocVector(INTSXP, *nr));
    PROTECT(pscore = allocVector(REALSXP, 1L));
    pvtmp = (double *) R_alloc(m, sizeof(double)); 
    for(i=0; i<m; i++) pvtmp[i] = 0.0;
    for(i=0; i<*nr; i++) INTEGER(pars)[i] = 0L;
    REAL(pscore)[0]=0.0;
    fitch8(data1, nr, INTEGER(pars), INTEGER(node), INTEGER(edge), INTEGER(l), weight, pvtmp, REAL(pscore));
    
    UNPROTECT(2);
    if(INTEGER(ps)[0]==1)return(pscore);
    else return(pars); 
}

//, double *pvec




void FN4(int *dat, int *res, int *nr, int *node, int *edge, int *nl, int *pc, double *weight, double *tmpvec, double *pvec) { 
    int i=0L, ni, le, ri;
    while(i< *nl) {
        ni = node[i] - 1L;
        le = edge[i] - 1L;
        ri = edge[i+1L] - 1L;
        if(pc[i+1L]==0L){
	        pvec[ni] = tmpvec[le] + tmpvec[ri];
	        fitch54(&res[ni * (*nr)], &dat[(edge[i]-1L) * (*nr)], &dat[ri * (*nr)], nr, weight, &pvec[ni]);              
        }    
        else{ 
            pvec[ni] = tmpvec[le] + pvec[ri];
	        fitch54(&res[ni * (*nr)], &dat[le * (*nr)], &res[ri * (*nr)], nr, weight, &pvec[ni]);   
        }
        i++;
        i++;
    }
}


void sibs(int *node, int *n, int *start, int *end){
    int tmp, k, i;
    tmp=node[0]; 
    k=node[0];     
    start[k]=0L; 
    for (i = 0L; i < *n; i++) {
        tmp = node[i];
        if(tmp!=k){
            end[k] = i-1L;
            start[tmp] = i; 
            k=tmp;
        }   
    }
    end[tmp] = i-1L;
}


void fnindex(int *nodes, int* edges, int *nNodes,  int *start, int *end, int *root, int *res1, int *res2, int *pc){
    int i, j, p, k, ni, nj, m;
    k=0L;
    for(i=0; i<*nNodes; i++){
        m = *nNodes-(1L+i);
        p = nodes[m];
        ni = edges[m];
        for(j=start[p]; j<=end[p]; j++){
            nj = edges[j]; 
            if(ni!=nj){
                res1[k] = nj;
                res2[k] = ni;
                pc[k] = 0L;  
                k++;
            }
        } 
        if(p!=*root){
            res1[k] = p;
            res2[k] = ni;
            pc[k] = 1L;                
            k++;
        }     
    }
}


void fnhelp(int *node, int * edge, int *n, int *m, int *root, int *edge2, int *node2, int *pc){
    int *startv, *endv, i;
    startv = (int *) R_alloc(*m, sizeof(int));
    endv = (int *) R_alloc(*m, sizeof(int));
    for(i=0; i<*m; i++){
        startv[i] = 0L;
        endv[i] = 0L;
    }
    sibs(node, n, startv, endv);
    fnindex(node, edge, n, startv, endv, root, edge2, node2, pc);         
}


SEXP FNALL_NNI(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP my, SEXP root){   
    int *nr=INTEGER(nrx), m=INTEGER(mx)[0], i,  *n=INTEGER(l);  //*pars,
    double *pvtmp, *pvtmp2, pscore=0.0;  
    SEXP pvec1, pvec2, res; 
    int *pc, *edge2, *node2;
/* edge2, node2, pc ausserhalb definieren? */        
    edge2 = (int *) R_alloc(2L * *n, sizeof(int));
    node2 = (int *) R_alloc(2L * *n, sizeof(int));
    pc = (int *) R_alloc(2L * *n, sizeof(int));
    
//    pvtmp2 = (double *) R_alloc(m, sizeof(double));
    PROTECT(res = allocVector(VECSXP, 2L));
    PROTECT(pvec1 = allocVector(REALSXP, m));
    PROTECT(pvec2 = allocVector(REALSXP, m));
    
    pvtmp2 = REAL(pvec2);
    pvtmp = REAL(pvec1);
    for(i=0; i<m; i++){
        pvtmp[i] = 0.0;
        pvtmp2[i] = 0.0;
    }
    fnhelp(INTEGER(node), INTEGER(edge),  n, &m, INTEGER(root), edge2, node2, pc); 
    fitch9(data1, nr, INTEGER(node), INTEGER(edge), INTEGER(l), weight, pvtmp, &pscore); 
    FN4(data1, data2, nr, node2, edge2, INTEGER(my), pc, weight, pvtmp, pvtmp2); // pars,
    
    SET_VECTOR_ELT(res, 0, pvec1);
    SET_VECTOR_ELT(res, 1, pvec2);    
    UNPROTECT(3);
    return(res); 
}


SEXP FNALL5(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP my, SEXP root){   
    int *nr=INTEGER(nrx), m=INTEGER(mx)[0], i,  *n=INTEGER(l);  //*pars,
    double *pvtmp, *pvtmp2, pscore=0.0;  
    SEXP pvec; 
    // fnhelp
    int *pc, *edge2, *node2;
/* edge2, node2, pc ausserhalb definieren? */        
    edge2 = (int *) R_alloc(2L * *n, sizeof(int));
    node2 = (int *) R_alloc(2L * *n, sizeof(int));
    pc = (int *) R_alloc(2L * *n, sizeof(int));
    
//    pars = (int *) R_alloc(*nr, sizeof(int)); // raus     
    pvtmp2 = (double *) R_alloc(m, sizeof(double));
    
    PROTECT(pvec = allocVector(REALSXP, m));
    
    pvtmp = REAL(pvec);
    for(i=0; i<m; i++){
        pvtmp[i] = 0.0;
        pvtmp2[i] = 0.0;
    }
    fnhelp(INTEGER(node), INTEGER(edge),  n, &m, INTEGER(root), edge2, node2, pc);
//    fitch8(data1, nr, pars, INTEGER(node), INTEGER(edge), INTEGER(l), weight, pvtmp, &pscore);  
    fitch9(data1, nr, INTEGER(node), INTEGER(edge), INTEGER(l), weight, pvtmp, &pscore); 
//    FN3(data1, data2, nr, pars, node2, edge2, INTEGER(my), pc, weight, pvtmp, pvtmp2);
    FN4(data1, data2, nr, node2, edge2, INTEGER(my), pc, weight, pvtmp, pvtmp2); // pars,
    for(i=0; i<m; i++) pvtmp[i] += pvtmp2[i];
// return(pvtmp[edge])??    
    UNPROTECT(1);
    return(pvec); 
}

// inside optNNI Ziel 3* schneller  , double best
void fitchquartet(int *dat1, int *dat2, int *dat3, int *dat4, int *nr, double *weight, double *pars){   
    int k, tmp1, tmp2;  
    pars[0] = 0.0; 
    for(k = 0; k < *nr; k++){
        tmp1 = dat1[k] & dat2[k];
        tmp2 = dat3[k] & dat4[k];  
        if(!tmp1){
            tmp1 = dat1[k] | dat2[k];
            pars[0]+=weight[k];
        }
        if(!tmp2){
            tmp2 = dat3[k] | dat4[k];
            pars[0]+=weight[k];
        }
        tmp1 = tmp1 & tmp2;
        if(!tmp1){
            pars[0]+=weight[k];
        }
    }
}

// weight raus  double *weight,
void fitchQuartet(int *index, int *n, int *nr, double *psc1, double *psc2, double *weight, double *res){
    int i, e1, e2, e3, e4;
    for(i=0; i<*n; i++){ 
        e1=index[(i* 6L)] - 1L;
        e2=index[1L + (i* 6L)] - 1L;
        e3=index[2L + (i* 6L)] - 1L;
        e4=index[3L + (i* 6L)] - 1L;

        if(index[5L + (i* 6L)] == 1){
            fitchquartet(&data2[e1 * (*nr)], &data1[e2 * (*nr)], &data1[e3 * (*nr)], &data1[e4 * (*nr)], nr, weight, &res[i]);
            res[i] += psc2[e1] + psc1[e2] + psc1[e3] + psc1[e4]; // stimmt
        } 
        else{
            fitchquartet(&data1[e1 * (*nr)], &data1[e2 * (*nr)], &data1[e3 * (*nr)], &data1[e4 * (*nr)], nr, weight, &res[i]);
            res[i] += psc1[e1] + psc1[e2] + psc1[e3] + psc1[e4]; 
        }
    }
}



