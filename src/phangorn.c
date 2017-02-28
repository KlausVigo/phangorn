/* 
 * phangorn.c
 *
 * (c) 2008-2017  Klaus Schliep (klaus.schliep@gmail.com)
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



// off-diagonal
#define DINDEX(i, j) n*(i - 1) - i * (i - 1)/2 + j - i - 1
// with diagonal (+i), R index (+1)
#define DINDEX2(i, j) n*(i - 1) - i * (i - 1)/2 + j - 1

// index likelihood pml
// need to define nr, nc, nTips, nNodes k
#define LINDEX(i) (i-nTips) * (nr*nc) //+ k * nTips * (nr * nc)
#define LINDEX2(i, k) (i-nTips) * (nr*nc) + k * nTips * (nr * nc)
#define LINDEX3(i, k) (i-*nTips-1L) * (*nr* *nc) + k * *nTips * (*nr * *nc)

// index sankoff
#define SINDEX(i) i * (nr*nc) 


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



void nodeH(int *edge, int *node, double *el, int *l,  double *res){
    int ei, i;
    for (i=*l-1L; i>=0; i--) {
        ei = edge[i] - 1L;
        res[ei] = res[node[i]-1L] + el[ei];
    }
}


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
    }
}
*/

 
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


void getdP(double *eva, double *ev, double *evi, int m, double el, double w, double *result){
    int i, j, h;
    double  res; //tmp[m],
    double *tmp;
    tmp = malloc(m * sizeof(double));
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
    tmp = malloc(m * sizeof(double));
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
    tmp = malloc(m * sizeof(double));
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
    tmp = malloc(m * sizeof(double));
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
*/


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


/*
SEXP AllChildren(SEXP children, SEXP parent, SEXP M){
    int i, j, k, l=0L, m=INTEGER(M)[0], *tab, p;   
    R_len_t n=length(parent); 
    SEXP RESULT, TMP;
    tab = (int*)R_alloc(m, sizeof(int));
    for(i=0; i<m; i++)tab[i]=0L;
    j=0;    
    p = INTEGER(parent)[0];
    for(i=0; i<n; i++){
        if(INTEGER(parent)[i]!=p){
            p = INTEGER(parent)[i]; 
            j=j+1;
        } 
        tab[j] += 1L;
    }
//    for(i=0; i<n; i++) tab[INTEGER(parent)[i] - 1L] ++;  // 7 Zeilen weniger      
    PROTECT(RESULT = allocVector(VECSXP, m));

    i=0L;    
    while(l<n){    
        k=tab[i];        
        PROTECT(TMP = allocVector(INTSXP, k));  
        p = INTEGER(parent)[l]-1;
        for(j=0; j<k; j++){
            INTEGER(TMP)[j] = INTEGER(children)[l];
            l++;
        } 
        SET_VECTOR_ELT(RESULT, p, TMP);
        UNPROTECT(1);
        i++;
    }
    UNPROTECT(1);
    return(RESULT);
}
*/

void AllKids(int *children, int *parents, int *nTips, int *nNode, int *lp, int *kids, int *lkids, int *pkids){
    int i, k, m=nNode[0], p; // l=0L, *tab , j 
    int n=lp[0]; 
    for(i=0; i<m; i++){
        pkids[i]=0L;
        lkids[i]=0L;
    }
    for(i=0; i<lp[0]; i++)kids[i]=0L;
//    j=0;
    p = 0L;
    for(i=0; i<n; i++){
        p = parents[i] - 1L - nTips[0];
        pkids[p] += 1L;
    }
    for(i=0; i<*nNode; i++)lkids[i+1] = lkids[i] + pkids[i];
    
    i=0L;   
    k=0;
    p=0L;
    for(i=0; i<n; i++){
        if(parents[i]!=p){
            p=parents[i];
            k=lkids[p- nTips[0] -1L];
        }
        else k++;
        kids[k] = children[i];
    }
    
}


/*
library(phangorn)
tree =  rtree(10)

allDesc = function(x, node){
  x = reorder(x, "postorder")
  parent = x$edge[, 1]
  children = x$edge[, 2]
  .Call("AllDesc", as.integer(children), as.integer(parent), as.integer(max(parent)), as.integer(node)) 
}

allDesc(tree, 14)
 
 
SEXP AllDesc(SEXP child, SEXP parent, SEXP M, SEXP NODE){
    int i, m=INTEGER(M)[0]+1, *tab, *res, p, node=INTEGER(NODE)[0];   
    R_len_t n=length(parent); 
    SEXP RESULT;
    tab = (int*)R_alloc( m, sizeof(int));
    for(i=0; i<m; i++)tab[i]=0L;
    tab[node] = 1L;
    PROTECT(RESULT = allocVector(INTSXP, n));
        res = INTEGER(RESULT);
    for(i=0; i<n; i++)res[i]=0L;

   for(i=n-1L; i>=0L; i--){
        p = INTEGER(parent)[i];
        if(tab[p]==1L){
            res[i] = 1L;
            tab[INTEGER(child)[i]] = 1L; 
        } 
    }
    
    UNPROTECT(1);
    return(RESULT);
}
 */

// std::merge
void cisort(int *x, int *y, int a, int b, int *res){
   int i, j, k;    
   i=0;
   j=0;
   k=0;
   int xi=x[0];
   int yi=y[0];  
   while(k<((a)+(b))){
      if(i<(a)){
          if( (xi<yi) | (j==b) ){  //-1L
              res[k]=xi;      
              i++;
              if(i<(a))xi=x[i];   
              k++;     
          }
          else{
              j++;
              res[k]=yi;
              if(j<(b))yi=y[j];  
              k++;
          }
        }
        else{
              j++;
              res[k]=yi;
              if(j<(b))yi=y[j];  
              k++;
          }
    }
}  


// faster cophenetic 
void C_bipHelp(int *parents, int *children, int *ntips, int *mp, int *l, int *ltips, int *ptips){
   int p, k, i;
   for(i=0; i<*ntips; i++)ltips[i]=1L;
   for(i=*ntips; i<*mp; i++)ltips[i]=0L;
   for(i=0; i<*l; i++){
       p = parents[i]-1L;
       k = children[i]-1L;
       ltips[p]+=ltips[k];
   }
   for(i=0; i<(*mp+1); i++)ptips[i]=0L;
   for(i=0; i<*mp; i++)ptips[i+1]=ptips[i] + ltips[i];
}   


void C_bip2(int *parents, int *children, int *ntips, int *mp, int *l, int *ltips, int *ptips, int *tips){
    int eins=1L, i, j, p, pi, ci, ltmp; 
    int *tmp, *tmp2;  
    tmp = (int *) R_alloc(*mp, sizeof(int));
    tmp2 = (int *) R_alloc(*mp, sizeof(int));
    for(i=0; i<*ntips; i++)tips[i]=i+1L;
    for(i=*ntips; i<ptips[*mp]; i++)tips[i]=0L; 
    p=parents[0];

    tmp[0] = 0L; //children[0]; 
    ltmp=  0L; //1L;
    for(i=0; i<*l; i++){ 
        pi = parents[i]; 
        ci = children[i];
        if(pi==p){
             if(ci < (*ntips+1L)){
                 cisort(&ci, tmp, eins, ltmp, tmp2);            
                 ltmp += 1L;
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];
             }
             else{
                 cisort(&tips[ptips[ci-1L]], tmp, (ltips[ci-1L]), ltmp, tmp2);                       
                 ltmp += ltips[ci-1L]; //  lch[ci]; 
//               ltmp +=   lch[ci];
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];                                
             } 
//             kl[pi]=k; 
//             lch[pi] = ltmp;
        }  
        else{
            for(j=0; j<ltmp; j++) tips[ptips[p-1L]+j] = tmp2[j];//tmp2[j]
            if(ci < (*ntips+1)){ 
                 tmp[0]=ci;
                 ltmp=1L; 
            } 
            else{ 
                ltmp=ltips[ci-1L];
                for(j=0; j<ltmp; j++)tmp[j] = tips[ptips[ci-1L]+j]; // , ci-1L))[j];
            }
//            k += 1L; 
            p = pi;
        }
    }
    for(j=0; j<ltmp; j++) tips[ptips[p-1L]+j] = tmp2[j];
}   

// doppelt
int give_index3(int i, int j, int n)
{
    if (i > j) return(DINDEX(j, i));
    else return(DINDEX(i, j));
}

// faster and less memory consuming cophenetic
void copheneticHelp(int *left, int *right, int *ll, int *lr, int h, double *nh, int *nTips, double *dm){
    int i, j, ind;
    for(i=0; i<*ll; i++){
        for(j=0; j<*lr; j++){
            ind = give_index3(left[i], right[j], *nTips);
            dm[ind] = 2.0*nh[h] - nh[left[i]-1L] - nh[right[j]-1L]; 
        }   
    }
}     


void C_coph(int *tips, int *kids, int *ptips, int *pkids, int *ltips, int *lkids, int*Nnode, double *nh, int *nTips, double *dm){
    int h, j, k, lk, pk, lt, rt, leftk, rightk;
    for(h=0; h<*Nnode; h++){
        lk=lkids[h]; 
        pk=pkids[h];
        for(j=0; j<(lk-1L); j++){
            leftk=kids[pk+j] - 1L;
            lt=ptips[leftk];
            for(k=j+1L; k<lk; k++) {
                rightk=kids[pk+k] - 1L;
                rt = ptips[rightk];
                copheneticHelp(&tips[lt], &tips[rt], &ltips[leftk], &ltips[rightk], (*nTips+h), nh, nTips, dm);
            }
        }
    }
}


void C_cophenetic(int *children, int *parents, double *el, int *lp, int *m, int *nTips, int *nNode, double *res){
    double *nh, maxNH; 
    int i, lt; 
    int *kids, *lkids, *pkids;
    int *tips, *ltips, *ptips;
    nh = (double *) calloc(*m, sizeof(double)); 
    kids = (int *) R_alloc(*lp, sizeof(int));
    lkids = (int *) R_alloc(*nNode + 1L, sizeof(int));
    pkids = (int *) R_alloc(*nNode, sizeof(int));
    ltips = (int *) R_alloc(*m, sizeof(int));
    ptips = (int *) R_alloc(*m + 1L, sizeof(int));
    //nodeH(int *edge, int *node, double *el, int *l,  double *res)
    nodeH(children, parents, el, lp,  nh);
    maxNH=nh[0];
    for(i=1; i<*m; i++)if(maxNH<nh[i]) maxNH=nh[i];
    for(i=0; i<*m; i++)nh[i] = maxNH - nh[i]; 
//    tmp <- .C("AllKids", kids, parents, nTips, nNode, lp, integer(lp), integer(nNode+1L),
//              integer(nNode))
// void AllKids(int *children, int *parents, int *nTips, int *nNode, int *lp, int *kids, int *lkids, int *pkids){
    AllKids(children, parents, nTips, nNode, lp, kids, lkids, pkids);
//tmp2 = .C("C_bipHelp", parents, kids, nTips, m, lp, integer(m), integer(m+1L))
    C_bipHelp(parents, children, nTips, m, lp, ltips, ptips);
// tips <- .C("C_bip2", parents, kids, nTips, m, lp, ltips=tmp2[[6]], ptips=tmp2[[7]], integer(sum(tmp2[[6]])))[[8]]    
    lt = 0;
    for(i=0; i<*m; i++)lt += ltips[i];
    tips = (int *) R_alloc(lt, sizeof(int));
    C_bip2(parents, children, nTips, m, lp, ltips, ptips, tips);
    //void coph(int *tips, int *kids, int *ptips, int *pkids, int *ltips, int *lkids, int*Nnode, double *nh, int *nTips, double *dm)
    C_coph(tips, kids, ptips, lkids, ltips, pkids, nNode, nh, nTips, res);
}


// a bit faster 
SEXP C_bip(SEXP parent, SEXP child, SEXP nTips, SEXP maxP){ //, SEXP Nnode){
   int eins=1L, i, j, k, l=length(child), *tmp, *tmp2, *lch, *kl, pi, ci, p, nt=INTEGER(nTips)[0], mp=INTEGER(maxP)[0], ltmp; 
   SEXP ans, ktmp;
   tmp = (int *) R_alloc(mp, sizeof(int));
   tmp2 = (int *) R_alloc(mp, sizeof(int));
   lch = (int *) R_alloc(mp+1L, sizeof(int));
   kl = (int *) R_alloc(mp+1L, sizeof(int));
   PROTECT(ans = allocVector(VECSXP, mp)); //INTEGER(Nnode)[0])); 
   for(i=0; i<nt; i++) SET_VECTOR_ELT(ans, i, ScalarInteger(i+1L)); 
   p=INTEGER(parent)[0];
   pi = INTEGER(parent)[1];
   k=0L;
   kl[p]=0;
   lch[p]=1;
   tmp[0] = INTEGER(child)[0]; 
   ltmp=1L;
   for(i=1; i<l; i++){ 
        pi = INTEGER(parent)[i]; 
        ci = INTEGER(child)[i];
        if(pi==p){
             if(ci < (nt+1L)){
                 cisort(&ci, tmp, eins, ltmp, tmp2);            
                 ltmp += 1L;
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];
             }
             else{
                 cisort(INTEGER(VECTOR_ELT(ans, ci-1L)), tmp, (lch[ci]), ltmp, tmp2);                       
                 ltmp += lch[ci]; 
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];                                
             }
             kl[pi]=k; 
             lch[pi] = ltmp;
        }  
        else{
            PROTECT(ktmp = allocVector(INTSXP, ltmp));
            for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
// k???           
            SET_VECTOR_ELT(ans, p-1L, ktmp); 
            UNPROTECT(1); // ktmp

            if(ci < (nt+1)){ 
                 tmp[0]=ci;
                 ltmp=1L; 
            } 
            else{ 
                ltmp=lch[ci];
                for(j=0; j<ltmp; j++)tmp[j] = INTEGER(VECTOR_ELT(ans, ci-1L))[j];
            }
            k += 1L; 
            p = pi;
        }
   }
   PROTECT(ktmp = allocVector(INTSXP, ltmp));// mp+1L
   for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
   SET_VECTOR_ELT(ans, pi-1L, ktmp);
   UNPROTECT(2);
   return(ans);  
}


SEXP C_bipart(SEXP parent, SEXP child, SEXP nTips, SEXP maxP){ //, SEXP Nnode){
   int eins=1L, i, j, k, l=length(child), *tmp, *tmp2, *lch, *kl, pi, ci, p, nt=INTEGER(nTips)[0], mp=INTEGER(maxP)[0], ltmp; 
   SEXP ans, ktmp;
   int nnode=1L;
   for(i=1; i<l; i++){
       if(INTEGER(parent)[i-1L] != INTEGER(parent)[i])nnode+=1L;
   }
   tmp = (int *) R_alloc(mp, sizeof(int));
   tmp2 = (int *) R_alloc(mp, sizeof(int));
   lch = (int *) R_alloc(mp+1L, sizeof(int));
   kl = (int *) R_alloc(mp+1L, sizeof(int));
// Nnode  
   PROTECT(ans = allocVector(VECSXP, nnode)); //INTEGER(Nnode)[0]));  
   p=INTEGER(parent)[0];
   pi=INTEGER(parent)[1];
   k=0L;
   kl[p]=0;
   lch[p]=1;
   tmp[0] = INTEGER(child)[0]; 
   ltmp=1L;
   for(i=1; i<l; i++){ 
        pi = INTEGER(parent)[i]; 
        ci = INTEGER(child)[i];
        if(pi==p){
             if(ci < (nt+1L)){
                 cisort(&ci, tmp, eins, ltmp, tmp2);            
                 ltmp += 1L;
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];
             }
             else{
                 cisort(INTEGER(VECTOR_ELT(ans, kl[ci])), tmp, lch[ci], ltmp, tmp2);                       
                 ltmp += lch[ci]; 
                 for(j=0; j<ltmp; j++) tmp[j] = tmp2[j];                                
             }
             kl[pi]=k; 
             lch[pi] = ltmp;
        }  
        else{
            PROTECT(ktmp = allocVector(INTSXP, ltmp));
            for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
// k???           
            SET_VECTOR_ELT(ans, k, ktmp); 
            UNPROTECT(1); // ktmp

            if(ci < (nt+1)){ 
                 tmp[0]=ci;
                 ltmp=1L; 
            } 
            else{ 
                ltmp=lch[ci];
                for(j=0; j<ltmp; j++)tmp[j] = INTEGER(VECTOR_ELT(ans, kl[ci]))[j];
            }
            k += 1L; 
            p = pi;
        }
   }
// k ??   
   PROTECT(ktmp = allocVector(INTSXP, ltmp));// mp+1L
   for(j=0; j<ltmp; j++)INTEGER(ktmp)[j] = tmp2[j];
   SET_VECTOR_ELT(ans, k, ktmp);
   UNPROTECT(2);
   return(ans);  
}




