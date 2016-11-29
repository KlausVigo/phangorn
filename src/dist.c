/* 
 * dist.c
 *
 * (c) 2008-2016 Klaus Schliep (klaus.schliep@gmail.com)
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
// #include "dist.h" 


// off-diagonal
#define DINDEX(i, j) n*(i - 1) - i * (i - 1)/2 + j - i - 1
// with diagonal (+i), R index (+1)
#define DINDEX2(i, j) n*(i - 1) - i * (i - 1)/2 + j - 1

// #define threshold parameters


int give_index(int i, int j, int n)
{
    if (i > j) return(DINDEX(j, i));
    else return(DINDEX(i, j));
}
 

int give_index2(int i, int j, int n)
{
    if (i > j) return(DINDEX2(j, i));
    else return(DINDEX2(i, j));
}
 

 
void giveIndex(int *left, int* right, int *ll, int *lr, int *n, int *res){
    int i, j, k;
    k=0;
    for (i = 0; i < *ll; i++){
        for (j = 0; j < *lr; j++){
             res[k] = give_index(left[i], right[j], *n);
             k++;
             }    
        }
    }


void giveIndex2(int *left, int* right, int *ll, int *lr, int *n, int *res){
    int i, j, k;
    k=0;
    for (i = 0; i < *ll; i++){
        for (j = 0; j < *lr; j++){
             res[k] = give_index2(left[i], right[j], *n);
             k++;
             }    
        }
    }




void PD(int *x, int *y, int *n, int *weight){
    int i, k; //n =length(x) 
    for(i=0; i< *n; i++){
        k=give_index(x[i], y[i], *n);
        weight[k]++;
    }
}


void pwIndex(int *left, int* right, int *l, int *n, double *w, double *res){
    int i, k;
    k=0;
    for (i = 0; i < *l; i++){
        k = give_index2(left[i], right[i], *n);
        res[k] += w[i];        
        }
    }



SEXP PWI(SEXP LEFT, SEXP RIGHT, SEXP L, SEXP N, SEXP W, SEXP LI){
    int i, li=INTEGER(LI)[0];    
    SEXP res;  
    PROTECT(res = allocVector(REALSXP, li));
    for(i = 0; i < li; i++)REAL(res)[i] = 0.0;
    pwIndex(INTEGER(LEFT), INTEGER(RIGHT), INTEGER(L), INTEGER(N), REAL(W), REAL(res));    
    UNPROTECT(1);    
    return(res);
}



void C_fhm(double *v, int *n){
    unsigned int level, i, j; 
    unsigned int start, step, num_splits;
    unsigned int max_n = (unsigned int)*n;
    double vi, vj;
    num_splits = (1 << (*n));
    step = 1;
    for(level = 0; level < max_n; level++){
        start = 0;
        while(start < (num_splits-1)){
            for(i = start; i < (start + step); i++){
                j = i + step;
                vi = v[i];
                vj = v[j];
                v[i] = vi + vj;
                v[j] = vi - vj;
            }
            start = start + 2*step;    
        }
        step *= 2;        
    }
}


void distance_hadamard(double *d, int n) {
    unsigned int num_splits;
    unsigned int x, r, nr, p, b, e;
    unsigned int odd = 1;                // The inner while loop can only terminate with odd == 1 so we don't need to set it inside the for loop.
    double cost, best_cost;
        
    num_splits = (1 << (n - 1));
        
    for (x = 1; x < num_splits; ++x) {
        r = (x - 1) & x;                // r = x without LSB
        nr = (r - 1) & r;                // nr = r without LSB
            
        if (nr) {                        // If x contains 1 or 2 bits only, then it has already been computed as a pairwise distance.
            b = x - r;                    // b = LSB of x: the "fixed" taxon in the pair.
            best_cost = 1e20;
            e = 0;                        // e holds bits to the right of the current p.
                
            while (1) {
                p = r - nr;                // p = 2nd half of pair
                cost = d[nr + e] + d[p + b];
                if (cost < best_cost) best_cost = cost;
                    
                if (!nr && odd) break;    // Ensure we get the (LSB with reference taxon) pair when there are an odd number of taxa
                r = nr;
                e += p;
                nr = (r - 1) & r;        // nr = r without LSB
                odd ^= 1;
                }
                d[x] = best_cost;
            }
        }

        d[0] = 0.0;
    }
    
// int num_splits raus  
void pairwise_distances(double *dm, int n, double *d) {
    int k=0;
    unsigned int offset;
        for (int i = 0; i < (n-1); ++i) {
            for (int j = (i+1); j < n; ++j) {
// Calculate the offset within the array to put the next value
                offset = (1 << i);
                if (j < n - 1) {            // If i == n - 1 then this is a distance between the reference taxon and some other taxon.
                    offset += (1 << j);        // Note that "+" is safe since (1 << i) and (1 << j) are always different bits.
                }
                d[offset]=dm[k];
                k++;
            }
        }
    }


SEXP dist2spectra(SEXP dm, SEXP nx, SEXP ns) {   
    int n = INTEGER(nx)[0];
    int nsp = INTEGER(ns)[0];   
    double *res;
    SEXP result;
    PROTECT(result = allocVector(REALSXP, nsp));
    res = REAL(result);
    pairwise_distances(REAL(dm), n, res); //nsp, 
    distance_hadamard(res, n);
    UNPROTECT(1);
    return(result);
}


// speed up some code for NJ    
void out(double *d, double *r, int *n, int *k, int *l){
    int i, j; 
    double res, tmp;
    k[0]=1;
    l[0]=2;
    res = d[1] - r[0] - r[1];
    for(i = 0; i < (*n-1); i++){
        for(j = i+1; j < (*n); j++){
                tmp = d[i*(*n)+j] - r[i] - r[j];
                if(tmp<res){
                    k[0]=i+1;
                    l[0]=j+1;
                    res = tmp;
                    }
            }
                
        }        
    }


// hamming distance    
void distHamming(int *x, double *weight, int *nr, int *l, double *d){
    int i, j, k, m;
    k = 0L;
    for(i = 0; i< (*l-1L); i++){
        for(j = (i+1L); j < (*l); j++){
             for(m=0; m<(*nr); m++){
                 if(!(x[i*(*nr) + m] & x[j*(*nr) + m])) d[k] += weight[m]; 
                 } 
             k++;
        }
    }
}






