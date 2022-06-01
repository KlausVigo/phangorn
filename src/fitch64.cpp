#include <Rcpp.h>
#include "Fitch.h"
#include "phangorn_utils.h"

//#include <omp.h>
// parallel for
// simd

//using namespace Rcpp;

//Enable C++11, as we may want to have unsigned long long (uint64_t)
// [[Rcpp::plugins(cpp11)]]


/*
// Examples from Dirk to find minimum dist
double vecmin(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}

int vecminInd(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the position (+1 in R)
  return it - x.begin();
}
*/


std::vector< std::vector<uint64_t> > readFitch(const List &xlist, IntegerMatrix contr,
                                               int nSeq, int nChar, int nStates, int nBits, int m){
  int current_bit=0;
  std::vector< std::vector<uint64_t> > out(m * nSeq);
  std::vector<uint64_t> tmp; // tmp(nStates);
  for (int k = 0; k < nStates; ++k) tmp.push_back(0ull);
  for(int i=0; i<nSeq; ++i) {
    Rcpp::IntegerVector y(xlist[i]);
    current_bit=0;
    for(int j=0; j<nChar; ++j){
      for (int k = 0; k < nStates; ++k){
        if (contr(y[j], k) > 0){
          tmp[k] |= (1ull << current_bit);
        }
      }
      current_bit++;
      if (current_bit == BIT_SIZE){
        for (int k = 0; k < nStates; ++k){
          out[i].push_back(tmp[k]);
          tmp[k] = 0ull;
        }
        current_bit = 0;
      }
    }
    if (current_bit && (current_bit != BIT_SIZE)){
      for (; current_bit < BIT_SIZE; ++current_bit){
        for (int k = 0; k < nStates; ++k){
          tmp[k] |= (1ull << current_bit);
        }
      }
      for (int k = 0; k < nStates; ++k){
        out[i].push_back(tmp[k]);
        tmp[k] = 0ull;
      }
    }
    out[i].shrink_to_fit();
  }
  uint64_t tmp0=0ull;
  if(m>1){
    for(int i=nSeq; i<(m*nSeq); ++i) {
      for(int j=0; j<(nStates*nBits); ++j){
        out[i].push_back(tmp0);
      }
    out[i].shrink_to_fit();
    }
  }
  return out; // wrap(out)
}



IntegerMatrix getAnc(Fitch* obj, int i){
  int states = obj->nStates;
  int nBits = obj->nBits;
  uint64_t tmp;
  std::vector< std::vector<uint64_t> > vector = obj->X;
  uint64_t * seq;
  seq = vector[i - 1].data();
  IntegerMatrix res(BIT_SIZE*nBits, states);
  for (int i = 0; i < nBits; ++i){
    for (int j = 0; j < states; ++j){
      tmp = seq[j];
      for(int l=0; l<BIT_SIZE; ++l){
        if( (tmp >> l) & 1ull ) res(i*BIT_SIZE+l,j) = 1;
      }
    }
    seq += states;
  }
  return(res);
}

IntegerVector getAncAmb(Fitch* obj, int i){
  int states = obj->nStates;
  int nBits = obj->nBits;
  uint64_t tmp;
  std::vector< std::vector<uint64_t> > vector = obj->X;
  IntegerVector xx = IntegerVector::create(1, 2, 4, 8);
  uint64_t * seq;
  seq = vector[i - 1].data();
  IntegerVector res(BIT_SIZE*nBits);
  for (int i = 0; i < nBits; ++i){
    for (int j = 0; j < states; ++j){
      tmp = seq[j];
      for(int l=0; l<BIT_SIZE; ++l){
        if( (tmp >> l) & 1ull ) res(i*BIT_SIZE+l) += xx[j];
      }
    }
    seq += states;
  }
  return(res);
}



// Instead of looping through every single bit, you can instead loop through only the set bits, which can be faster if you expect bits to be sparsely set:
//  Assume the bit field is in (scalar integer) variable field.
/*
while (field){
  temp = field & -field;  //extract least significant bit on a 2s complement machine
  field ^= temp;  // toggle the bit off
  //now you could have a switch statement or bunch of conditionals to test temp
  //or get the index of the bit and index into a jump table, etc.
}
*/
// Works pretty well when the bit field is not limited to the size of a single data type, but could be of some arbitrary size. In that case, you can extract 32 (or whatever your register size is) bits at a time, test it against 0, and then move on to the next word.



void update_vector_generic(uint64_t * parent, const uint64_t * child1,
                           const uint64_t * child2, int nBits, int states){
  for (int i = 0; i < nBits; ++i){
    // OR the ANDs of all states
    uint64_t orvand = 0ull;
    // #pragma omp simd private(orvand)
    for (int j = 0; j < states; ++j) orvand |= (child1[j] & child2[j]);
    // store vectors at parent
    // #pragma omp simd
    for (int j = 0; j < states; ++j) {
      parent[j] = (child1[j] & child2[j]) | (~orvand & (child1[j] | child2[j]));
    }
    child1 += states;
    child2 += states;
    parent += states;
  }
}


void update_vector_4x4(uint64_t * parent, const uint64_t * child1,
                   const uint64_t * child2, int nBits, int states){
  uint64_t tmp0, tmp1, tmp2, tmp3;
  // #pragma omp simd
  for (int i = 0; i < nBits; ++i){
    uint64_t orvand = 0;
    tmp0 = (child1[0] & child2[0]);
    tmp1 = (child1[1] & child2[1]);
    tmp2 = (child1[2] & child2[2]);
    tmp3 = (child1[3] & child2[3]);
    orvand = tmp0 | tmp1 | tmp2 | tmp3;
    parent[0] = tmp0 | (~orvand & (child1[0] | child2[0]));
    parent[1] = tmp1 | (~orvand & (child1[1] | child2[1]));
    parent[2] = tmp2 | (~orvand & (child1[2] | child2[2]));
    parent[3] = tmp3 | (~orvand & (child1[3] | child2[3]));
    child1 += states;
    child2 += states;
    parent += states;
  }
}


void update_vector_2x2(uint64_t * parent, const uint64_t * child1,
                       const uint64_t * child2, int nBits, int states){
  uint64_t tmp0, tmp1;
  for (int i = 0; i < nBits; ++i){
    // OR the ANDs of all states
    uint64_t orvand = 0;
    tmp0 = (child1[0] & child2[0]);
    tmp1 = (child1[1] & child2[1]);
    orvand = tmp0 | tmp1;
    parent[0] = tmp0 | (~orvand & (child1[0] | child2[0]));
    parent[1] = tmp1 | (~orvand & (child1[1] | child2[1]));
    child1 += states;
    child2 += states;
    parent += states;
  }
}


void update_vector(uint64_t * parent, const uint64_t * child1,
                   const uint64_t * child2, int nBits, int states)
{
  if (states == 4)
    update_vector_4x4(parent, child1, child2, nBits, states);
  else if (states == 2)
    update_vector_2x2(parent, child1, child2, nBits, states);
  else
    update_vector_generic(parent, child1, child2, nBits, states);
}


void update_vector_single_4x4(uint64_t * parent, const uint64_t * child,
                          int nBits, int states){
  uint64_t tmp0, tmp1, tmp2, tmp3;
  for (int i = 0; i < nBits; ++i){
    uint64_t orvand = 0ull;
    tmp0 = (child[0] & parent[0]);
    tmp1 = (child[1] & parent[1]);
    tmp2 = (child[2] & parent[2]);
    tmp3 = (child[3] & parent[3]);
    orvand = tmp0 | tmp1 | tmp2 | tmp3;
    parent[0] = tmp0 | (~orvand & (child[0] | parent[0]));
    parent[1] = tmp1 | (~orvand & (child[1] | parent[1]));
    parent[2] = tmp2 | (~orvand & (child[2] | parent[2]));
    parent[3] = tmp3 | (~orvand & (child[3] | parent[3]));
    child += states;
    parent += states;
  }
}


void update_vector_single_2x2(uint64_t * parent, const uint64_t * child,
                              int nBits, int states){
  uint64_t tmp0, tmp1;
  for (int i = 0; i < nBits; ++i){
    uint64_t orvand = 0ull;
    tmp0 = (child[0] & parent[0]);
    tmp1 = (child[1] & parent[1]);
    orvand = tmp0 | tmp1;
    parent[0] = tmp0 | (~orvand & (child[0] | parent[0]));
    parent[1] = tmp1 | (~orvand & (child[1] | parent[1]));
    child += states;
    parent += states;
  }
}


void update_vector_single_generic(uint64_t * parent, const uint64_t * child,
                   int nBits, int states){
  for (int i = 0; i < nBits; ++i){
    uint64_t orvand = 0;
    for (int j = 0; j < states; ++j) orvand |= (child[j] & parent[j]);
    for (int j = 0; j < states; ++j) {
      parent[j] = (child[j] & parent[j]) | (~orvand & (child[j] | parent[j]));
    }
    child += states;
    parent += states;
  }
}


void update_vector_single(uint64_t * parent, const uint64_t * child,
                          int nBits, int states)
{
  if (states == 4)
    update_vector_single_4x4(parent, child, nBits, states);
  else if (states == 2)
    update_vector_single_2x2(parent, child, nBits, states);
  else
    update_vector_single_generic(parent, child, nBits, states);
}


void traverse(Fitch* obj, const IntegerMatrix & orig){
  int states = obj->nStates;
  int nBits = obj->nBits;

  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);

  int nl=desc.size();
  int unrooted = nl % 2;
  if(unrooted == 1) nl = nl-1;

  for(int k=0; k<nl; k+=2){
    update_vector(obj->X[anc[k] - 1].data(), obj->X[desc[k] - 1].data(),
                  obj->X[desc[k+1] - 1].data(), nBits, states);
  }
  if(unrooted){
    update_vector_single(obj->X[anc[nl] - 1].data(),
                         obj->X[desc[nl] - 1].data(), nBits, states);
  }
}


void traversetwice(Fitch* obj, const IntegerMatrix & orig, int nni){
  int states = obj->nStates;
  int nBits = obj->nBits;
  int nTips = obj->nSeq;
  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);

  if(nni > 0) nni = nTips - 1;
  else  nni = -1;

  int nl=desc.size();
  int unrooted = nl % 2;
//  int l = nl;
  if(unrooted == 1) nl = nl-1;

  for(int k=0; k<nl; k+=2){
    update_vector(obj->X[anc[k] - 1].data(), obj->X[desc[k] - 1].data(),
                  obj->X[desc[k+1] - 1].data(), nBits, states);
  }
  if(unrooted == 1){
    update_vector_single(obj->X[anc[nl] - 1].data(),
                         obj->X[desc[nl] - 1].data(), nBits, states);
    int a = desc[nl] -1;
    int b = desc[nl-1] -1;
    int c = desc[nl-2] -1;
    update_vector(obj->X[a + 2*nTips].data(), obj->X[b].data(),
                  obj->X[c].data(), nBits, states);
    update_vector(obj->X[b + 2*nTips].data(), obj->X[a].data(),
                  obj->X[c].data(), nBits, states);
    update_vector(obj->X[c + 2*nTips].data(), obj->X[a].data(),
                  obj->X[b].data(), nBits, states);
  }

  else{
    int a = desc[nl-1] -1;
    int b = desc[nl-2] -1;
    update_vector_single(obj->X[a + 2*nTips].data(),
                         obj->X[b].data(), nBits, states);
    update_vector_single(obj->X[b + 2*nTips].data(),
                         obj->X[a].data(), nBits, states);
  }
  nl -= 2;
  for(int i=nl; i>0; i-=2){
    int p = anc[i-1] -1;
    int c1 = desc[i-1] -1;
    int c2 = desc[i-2] -1;
    if(c1 > nni)update_vector(obj->X[c1 + 2*nTips].data(), obj->X[p + 2*nTips].data(),
       obj->X[c2].data(), nBits, states);
    if(c2 > nni)update_vector(obj->X[c2 + 2*nTips].data(), obj->X[p + 2*nTips].data(),
       obj->X[c1].data(), nBits, states);
  }
}


void acctran_help(uint64_t * child, const uint64_t * parent,
                                  int nBits, int states){
  for (int i = 0; i < nBits; ++i){
    uint64_t orvand = 0;
    for (int j = 0; j < states; ++j) orvand |= (child[j] & parent[j]);
    for (int j = 0; j < states; ++j) {
      child[j] = (child[j] & parent[j]) | (~orvand & child[j]);
    }
    child += states;
    parent += states;
  }
}


void acctran_traverse(Fitch* obj, const IntegerMatrix & orig){
  int states = obj->nStates;
  int nBits = obj->nBits;
  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);
  for(int i=0; i < anc.size(); ++i) {
    acctran_help(obj->X[desc[i] - 1].data(),
                 obj->X[anc[i] - 1].data(), nBits, states);
  }
}



void root_all_node(Fitch* obj, const IntegerMatrix orig)
{
  int states = obj->nStates;
  int nBits = obj->nBits;
  int nSeq = obj->nSeq;
//  std::vector< std::vector<uint64_t> > vector = obj->X;
  IntegerVector node = orig( _, 1);

  for(int i=0; i < node.size(); ++i) {
    int ni = node[i]-1;
// generic
    update_vector_single(obj->X[ni + 2*nSeq].data(),
                         obj->X[ni].data(), nBits, states);
  }
}


// needed for random.addition, SPR & TBR
void prep_spr(Fitch* obj, IntegerMatrix orig){
  traversetwice(obj, orig, 0L);
  root_all_node(obj, orig);
}


// generic, TODO: bitcount, 2x2, 4x4
double pscore_vector_generic(const uint64_t* x, const uint64_t* y, const NumericVector weight,
                     int nBits, int wBits, int states){
  double pscore = 0.0;
  uint64_t ones = ~0ull;
  uint64_t tmp = 0ull;
  for (int i = 0; i < wBits; ++i){
    uint64_t orvand = 0;
    for (int j = 0; j < states; ++j) orvand |= (x[j] & y[j]);
    tmp = ~orvand & ones;
    if(tmp>0ull){
      for(int l=0; l<64; ++l){
        if( (tmp >> l) & 1ull ) pscore += weight[i*BIT_SIZE + l];
      }
    }
    x += states;
    y += states;
  }
  for (int i = wBits; i < nBits; ++i){
    uint64_t orvand = 0;
    for (int j = 0; j < states; ++j) orvand |= (x[j] & y[j]);
    tmp = ~orvand & ones;
    pscore += popcnt64(tmp);
    x += states;
    y += states;
  }
  return(pscore);
}


double pscore_vector_4x4(const uint64_t* x, const uint64_t* y, const NumericVector weight,
                     int nBits, int wBits, int states){
  double pscore = 0.0;
  uint64_t ones = ~0ull;
  uint64_t tmp = 0ull;
  uint64_t orvand = 0;
  for (int i = 0; i < wBits; ++i){
    orvand = (x[0] & y[0]) | (x[1] & y[1]) | (x[2] & y[2]) | (x[3] & y[3]);
    tmp = ~orvand & ones;
    if(tmp>0ull){
      for(int l=0; l<64; ++l){
        if( (tmp >> l) & 1ull ) pscore += weight[i*BIT_SIZE + l];
      }
    }
    x += states;
    y += states;
  }
  for (int i = wBits; i < nBits; ++i){
    orvand = (x[0] & y[0]) | (x[1] & y[1]) | (x[2] & y[2]) | (x[3] & y[3]);
    tmp = ~orvand & ones;
    pscore += popcnt64(tmp);
    x += states;
    y += states;
  }
  return(pscore);
}


double pscore_vector_2x2(const uint64_t* x, const uint64_t* y, const NumericVector weight,
                         int nBits, int wBits, int states){
  double pscore = 0.0;
  uint64_t ones = ~0ull;
  uint64_t tmp = 0ull;
  uint64_t orvand = 0;
  for (int i = 0; i < wBits; ++i){
    orvand = (x[0] & y[0]) | (x[1] & y[1]);
    tmp = ~orvand & ones;
    if(tmp>0ull){
      for(int l=0; l<64; ++l){
        if( (tmp >> l) & 1ull ) pscore += weight[i*BIT_SIZE + l];
      }
    }
    x += states;
    y += states;
  }
  for (int i = wBits; i < nBits; ++i){
    orvand = (x[0] & y[0]) | (x[1] & y[1]);
    tmp = ~orvand & ones;
    pscore += popcnt64(tmp);
    x += states;
    y += states;
  }
  return(pscore);
}

double pscore_vector(const uint64_t* x, const uint64_t* y, const NumericVector weight,
                   int nBits, int wBits, int states)
{
  double res=0.0;
  if (states == 4)
    res=pscore_vector_4x4(x, y, weight, nBits, wBits, states);
  else if (states == 2)
    res=pscore_vector_2x2(x, y, weight, nBits, wBits, states);
  else
    res=pscore_vector_generic(x, y, weight, nBits, wBits, states);
  return(res);
}


int pscore_quartet(const uint64_t* a, const uint64_t* b, const uint64_t* c,
                   const uint64_t* d, const NumericVector weight,
                   int nBits, int wBits, int states)
{
  double pscore = 0.0;
  uint64_t ones = ~0ull;
  uint64_t tmp1, tmp2, tmp3;
  uint64_t e, f;
  for (int i = 0; i < wBits; ++i){
    uint64_t ou_ab = 0;
    uint64_t ou_cd = 0;
    uint64_t ou_ef = 0;
    for (int j = 0; j < states; ++j){
      ou_ab |= (a[j] & b[j]);
      ou_cd |= (c[j] & d[j]);
    }
    for (int j = 0; j < states; ++j) {
      e = (a[j] & b[j]) | (~ou_ab & (a[j] | b[j]));
      f = (c[j] & d[j]) | (~ou_cd & (c[j] | d[j]));
      ou_ef |= (e & f);
    }
    tmp1 = ~ou_ab & ones;
    tmp2 = ~ou_cd & ones;
    tmp3 = ~ou_ef & ones;
    if((tmp1 | tmp2 | tmp3) > 0ull){
      for(int l=0; l<BIT_SIZE; ++l){
        if( (tmp1 >> l) & 1ull ) pscore += weight[i*BIT_SIZE + l];
        if( (tmp2 >> l) & 1ull ) pscore += weight[i*BIT_SIZE + l];
        if( (tmp3 >> l) & 1ull ) pscore += weight[i*BIT_SIZE + l];
      }
    }
    a += states;
    b += states;
    c += states;
    d += states;
  }
  for (int i = wBits; i < nBits; ++i){
    uint64_t ou_ab = 0;
    uint64_t ou_cd = 0;
    uint64_t ou_ef = 0;
    for (int j = 0; j < states; ++j){
      ou_ab |= (a[j] & b[j]);
      ou_cd |= (c[j] & d[j]);
    }
    for (int j = 0; j < states; ++j) {
      e = (a[j] & b[j]) | (~ou_ab & (a[j] | b[j]));
      f = (c[j] & d[j]) | (~ou_cd & (c[j] | d[j]));
      ou_ef |= (e & f);
    }
    tmp1 = ~ou_ab & ones;
    tmp2 = ~ou_cd & ones;
    tmp3 = ~ou_ef & ones;
    pscore += popcnt64(tmp1) + popcnt64(tmp2) + popcnt64(tmp3);
    a += states;
    b += states;
    c += states;
    d += states;
  }
  return pscore;
}


IntegerMatrix pscore_nni(Fitch* obj, IntegerMatrix & M){
  int nr = M.nrow();
  IntegerMatrix res(nr, 3);
  std::vector< std::vector<uint64_t> > X = obj->X;
  int states = obj->nStates;
  int nBits = obj->nBits;
  int wBits = obj->wBits;
  NumericVector weight = obj->weight;
  int a=0, b=0, c=0, d=0;

  for (int i = 0; i < nr; i++) {
    a = M(i,0) - 1L;
    b = M(i,1) - 1L;
    c = M(i,2) - 1L;
    d = M(i,3) - 1L;
    res(i,0) = pscore_quartet(X[a].data(), X[b].data(), X[c].data(), X[d].data(),
                              weight, nBits, wBits, states);
    res(i,1) =  pscore_quartet(X[a].data(), X[c].data(), X[b].data(), X[d].data(),
                              weight, nBits, wBits, states);
    res(i,2) = pscore_quartet(X[b].data(), X[c].data(), X[a].data(), X[d].data(),
                              weight, nBits, wBits, states);
  }
  return(res);
}



/*
int get_quartet(Fitch* obj, IntegerVector & M){
  std::vector< std::vector<uint64_t> > X = obj->X;
  int states = obj->nStates;
  int nBits = obj->nBits;
  int wBits = obj->wBits;
  NumericVector weight = obj->weight;
  int res = pscore_quartet(X[M[0]].data(), X[M[1]].data(), X[M[2]].data(), X[M[3]].data(),
        weight, nBits, wBits, states);
  return(res);
}
*/

NumericVector pscore_vec(Fitch* obj, IntegerVector & edge_to, int node_from){
  // std::vector<double> res;
  int n = edge_to.size();
  NumericVector res(n);
  int states = obj->nStates;
  int nBits = obj->nBits;
  int wBits = obj->wBits;
  NumericVector weight = obj->weight;
  uint64_t * node_vec;
  node_vec = obj->X[node_from - 1L].data();
  //  #pragma omp parallel for num_threads(4)
  for(int i=0; i < edge_to.size(); ++i) {
    res[i] = pscore_vector(obj->X[edge_to[i]-1].data(),
                           node_vec, weight, nBits, wBits, states);
  }
  return(res);
}

// dist.hamming works for >31 states, TODO openMP
NumericVector hamming_dist(Fitch* obj){
  int i, j;
  size_t  ij;
  int states = obj->nStates;
  int nBits = obj->nBits;
  int wBits = obj->wBits;
  int nTips = obj->nSeq;
  R_xlen_t N;
  N = (R_xlen_t)nTips * (nTips-1)/2;
  std::vector< std::vector<uint64_t> > X = obj->X;
  NumericVector weight = obj->weight;
  NumericVector ans(N);
  ij = 0;
  i=0;
  j=1;
  for(j = 0 ; j < (nTips-1L) ; j++)
    for(i = j+1; i < nTips ; i++)
      ans[ij++] = pscore_vector(X[i].data(), X[j].data(), weight, nBits, wBits, states);
  return(ans);
}


IntegerVector sitewise_pscore(Fitch* obj, const IntegerMatrix & orig){
  int i,j;
  int states = obj->nStates;
  int nBits = obj->nBits;
  std::vector< std::vector<uint64_t> > vector = obj->X;
  IntegerVector pars(nBits * BIT_SIZE);

  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);

  int nl=desc.size();
  int unrooted = nl % 2;
  if(unrooted == 1) nl = nl-1;

  uint64_t * child1;
  uint64_t * child2;
  uint64_t * parent;

  // set all bits to one
  uint64_t ones = ~0ull;
  uint64_t tmp = 0ull;
  for(int k=0; k<nl; k+=2){
    child1 = vector[desc[k] - 1L].data();
    child2 = vector[desc[k+1] - 1L].data();
    parent = vector[anc[k] - 1L].data();
    for (i = 0; i < obj->nBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & child2[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & child2[j]) | (~orvand & (child1[j] | child2[j]));
      }
      child1 += states;
      child2 += states;
      parent += states;
      tmp = ~orvand & ones;

      for(int l=0; l<BIT_SIZE; ++l){
        pars[i*BIT_SIZE + l] += (int)((tmp >> l) & 1ull);
//        tmp >>= 1ull;
      }
    }
  }
  if(unrooted){
    child1 = vector[desc[nl] - 1].data();
    parent = vector[anc[nl] - 1].data();
    for (i = 0; i < obj->nBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & parent[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & parent[j]) | (~orvand & (child1[j] | parent[j]));
      }
      child1 += states;
      parent += states;
      uint64_t tmp = ~orvand & ones;
      for(int l=0; l<BIT_SIZE; ++l){
        pars[i*BIT_SIZE+l] += (int)((tmp >> l) & 1ull);
//        tmp >>= 1ullll;
      }
    }
  }
  return(pars);
}


double pscore(Fitch* obj, const IntegerMatrix & orig){
  int i,j;
  int states = obj->nStates;
  int nBits = obj->nBits;
  std::vector< std::vector<uint64_t> > vector = obj->X;
  double pars = 0;

  int p0 = obj->p0;

  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);

  int nl=desc.size();
  int unrooted = nl % 2;
  if(unrooted == 1) nl = nl-1;

  uint64_t * child1;
  uint64_t * child2;
  uint64_t * parent;

  // set all bits to one
  uint64_t ones = ~0ull;
  uint64_t tmp = 0ull;
  for(int k=0; k<nl; k+=2){
    child1 = vector[desc[k] - 1L].data();
    child2 = vector[desc[k+1] - 1L].data();
    parent = vector[anc[k] - 1L].data();
    for (int i = 0; i < obj->wBits; ++i){
//    for (i = 0; i < obj->nBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & child2[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & child2[j]) | (~orvand & (child1[j] | child2[j]));
      }
      child1 += states;
      child2 += states;
      parent += states;
      tmp = ~orvand & ones;
      for(int l=0; l<BIT_SIZE; ++l){
        if((tmp >> l) & 1ull) pars+= obj->weight[i*BIT_SIZE + l];
      }
    }
    for (int i = obj->wBits; i < nBits; ++i){
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & child2[j]);
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & child2[j]) | (~orvand & (child1[j] | child2[j]));
      }
      child1 += states;
      child2 += states;
      parent += states;
      tmp = ~orvand & ones;
      pars += popcnt64(tmp);
    }
  }
  if(unrooted){
    child1 = vector[desc[nl] - 1].data();
    parent = vector[anc[nl] - 1].data();
    for (i = 0; i < obj->wBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & parent[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & parent[j]) | (~orvand & (child1[j] | parent[j]));
      }
      child1 += states;
      parent += states;
      uint64_t tmp = ~orvand & ones;
      for(int l=0; l<BIT_SIZE; ++l){
        if((tmp >> l) & 1ull) pars+= obj->weight[i*BIT_SIZE + l];
      }
    }
    for (int i = obj->wBits; i < nBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & parent[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & parent[j]) | (~orvand & (child1[j] | parent[j]));
      }
      child1 += states;
      parent += states;
      uint64_t tmp = ~orvand & ones;
      pars += popcnt64(tmp);
    }
  }
  pars += p0;
  return(pars);
}


NumericVector pscore_node(Fitch* obj, const IntegerMatrix & orig){
  int i,j;
  int states = obj->nStates;
  int nBits = obj->nBits;
  std::vector< std::vector<uint64_t> > vector = obj->X;
  int nSeq = obj->nSeq;
  NumericVector pars(2 * nSeq);

  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);

  int nl=desc.size();
  int unrooted = nl % 2;
  if(unrooted == 1) nl = nl-1;

  uint64_t * child1;
  uint64_t * child2;
  uint64_t * parent;

  // set all bits to one
  uint64_t ones = ~0ull;
  uint64_t tmp = 0ull;
  for(int k=0; k<nl; k+=2){
    child1 = vector[desc[k] - 1L].data();
    child2 = vector[desc[k+1] - 1L].data();
    parent = vector[anc[k] - 1L].data();
    for (int i = 0; i < obj->wBits; ++i){
      //    for (i = 0; i < obj->nBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & child2[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & child2[j]) | (~orvand & (child1[j] | child2[j]));
      }
      child1 += states;
      child2 += states;
      parent += states;
      tmp = ~orvand & ones;
      for(int l=0; l<BIT_SIZE; ++l){
        if((tmp >> l) & 1ull) pars[anc[k] - 1L]+= obj->weight[i*BIT_SIZE + l];
      }
    }
    for (int i = obj->wBits; i < nBits; ++i){
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & child2[j]);
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & child2[j]) | (~orvand & (child1[j] | child2[j]));
      }
      child1 += states;
      child2 += states;
      parent += states;
      tmp = ~orvand & ones;
      pars[anc[k] - 1L] += popcnt64(tmp);
    }
  }
  if(unrooted){
    child1 = vector[desc[nl] - 1].data();
    parent = vector[anc[nl] - 1].data();
    for (i = 0; i < obj->wBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & parent[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & parent[j]) | (~orvand & (child1[j] | parent[j]));
      }
      child1 += states;
      parent += states;
      uint64_t tmp = ~orvand & ones;
      for(int l=0; l<BIT_SIZE; ++l){
        if((tmp >> l) & 1ull) pars[anc[nl] - 1L]+= obj->weight[i*BIT_SIZE + l];
      }
    }
    for (int i = obj->wBits; i < nBits; ++i){
      // OR the ANDs of all states
      uint64_t orvand = 0ull;
      for (j = 0; j < states; ++j) orvand |= (child1[j] & parent[j]);
      // store vectors at parent
      for (j = 0; j < states; ++j) {
        parent[j] = (child1[j] & parent[j]) | (~orvand & (child1[j] | parent[j]));
      }
      child1 += states;
      parent += states;
      uint64_t tmp = ~orvand & ones;
      pars [anc[nl] - 1L]+= popcnt64(tmp);
    }
  }
  return(pars);
}


NumericVector pscore_acctran(Fitch* obj, const IntegerMatrix & orig){
  int states = obj->nStates;
  int nBits = obj->nBits;
  int wBits = obj->wBits;
  NumericVector weight = obj->weight;
  int nSeq = obj->nSeq;
  NumericVector pars(2 * nSeq);

  IntegerVector anc = orig( _, 0);
  IntegerVector desc = orig( _, 1);

  for(int i=0; i <desc.size(); ++i) {
    pars[desc[i]-1] = pscore_vector(obj->X[anc[i]-1].data(),
                            obj->X[desc[i]-1].data(), weight, nBits, wBits, states);
  }
  return(pars);
}


RCPP_MODULE(Fitch_mod) {
    using namespace Rcpp;
    class_<Fitch>("Fitch")
        .constructor<RObject, int, int>("Default constructor")
//        .property("get_nr", &Fitch::getNR)
//        .property("get_nbits", &Fitch::getnBits)
//        .property("get_weight", &Fitch::getWeight)
//        .property("get_p0", &Fitch::getP0)
//        .method("get_quartet", &get_quartet)
//        .method("prep_nni", &prep_nni)
        .method("prep_spr", &prep_spr)
        .method("pscore_nni", &pscore_nni)
        .method("pscore", &pscore)
        .method("pscore_vec", &pscore_vec)
        .method("pscore_node", &pscore_node)
        .method("pscore_acctran", &pscore_acctran)
        .method("acctran_traverse", &acctran_traverse)
        .method("traverse", &traverse)
        .method("sitewise_pscore", &sitewise_pscore)
        .method("hamming_dist", &hamming_dist)
        .method("root_all_node", &root_all_node)
        .method("getAnc", &getAnc)
        .method("getAncAmb", &getAncAmb)
        .method("traversetwice", &traversetwice)
    ;
}

