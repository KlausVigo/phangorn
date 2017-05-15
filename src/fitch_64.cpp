#include <Rcpp.h>

using namespace Rcpp;

//Enable C++11, as we want to have unsigned long long (uint64_t)
//[[Rcpp::plugins(cpp11)]]

const int bits = 64; 


// Counting bits set, Brian Kernighan's wa
// from Bit Twiddling Hacks by Sean Eron Anderson seander@cs.stanford.edu
// count the number of bits set in v
uint64_t bitCount(uint64_t v){
  unsigned int c; // c accumulates the total bits set in v
  for (c = 0; v; c++)
  {
    v &= v - 1; // clear the least significant bit set
  }
  return c;
}


/* """This uses fewer arithmetic operations than any other known  
 implementation on machines with fast multiplication.
 It uses 12 arithmetic operations, one of which is a multiply.""" */

//inline ?
uint64_t bitcount64(uint64_t i) {
    i = i - ((i >> 1) & 0x5555555555555555);                                                                                                  
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);                                                                           
    return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;                                                                  
}

// set first state is ambiguous  
std::vector<uint64_t> contrast2uint(NumericMatrix xx) {
    int nr = xx.nrow(), nc = xx.ncol();
    double eps=1e-5;
    uint64_t tmp=0;
    // create vector for results
    std::vector<uint64_t> out;
    // add ambiguous state at 0 position
    for(int j=0; j<nc; ++j) tmp = tmp + (1 << j);
    out.push_back(tmp);
    for(int i=0; i<nr; ++i){
        tmp = 0L;
        for(int j=0; j<nc; ++j){
            if(xx(i,j) > eps){
                tmp = tmp + (1 << j);  
            }
        }
        out.push_back(tmp);
    }
    return out;
}

/*
 * needed as nibble can differ from 4 (nucleotides)
 * 0x77777777, 0x88888888 in White et al.
 */
uint64_t getEights(int nibble, int bits) {
    int l = bits / nibble;
    uint64_t tmp = (1 << (nibble - 1));
    uint64_t res = 0;
    for(int i=0; i<l; ++i){
        res += (tmp << (i * nibble));
    }
    return res;
}
 
 
uint64_t getSevens(int nibble, int bits) {
  int l = bits / nibble;
  uint64_t tmp = (1 << (nibble - 1)) - 1;
  uint64_t res = 0;
  for(int i=0; i<l; ++i){
    res += (tmp << (i * nibble));
  }
  return res;
}


uint64_t getOnes(int nibble, int bits) {
    int l = bits / nibble;
    uint64_t res = 0;
    for(int i=0; i<l; ++i){
        res += (1 << (i * nibble));
    }
    return res;
}


uint64_t getFifteen(int nibble, int bits) {
    uint64_t res = 0;
    for(int i=0; i<nibble; ++i){
        res += (1 << i);
    }
    return res;
}


std::vector< std::vector<uint64_t> > readFitch(RObject x, std::vector<uint64_t> contr, int nr, int bits, int nibble){
    // std::vector< std::vector<int> > statt List
    Rcpp::List xlist(x);
    uint64_t tmp; 
    int l=0L;
    int word = bits / nibble;
    int words = nr / word;
    int rest = nr % word;
    int all_words = words;
    if(rest>0)  all_words += 1;
    int n = xlist.size();
    std::vector< std::vector<uint64_t> > out(2 * n);
    for(int i=0; i<n; ++i) {
        Rcpp::NumericVector y(xlist[i]);
        l=0;
        for(int j=0; j<words; ++j){
            tmp = 0;
            for(int k=0; k<word; ++k){
                tmp += (contr[y[l]] << (k * nibble));
                l++;
            }
            out[i].push_back(tmp);
        }
        if(rest>0){
            tmp = 0;
            for(int k=0; k<rest; ++k){
                tmp += (contr[y[l]] << (k * nibble));
                l++;
            }
            // replace last free nibbles with ambiguous states 
            for(int k=rest; k<word; ++k){
                tmp += (contr[0] << (k * nibble));
                l++;
            }
            out[i].push_back(tmp);
        }
    }
    tmp=0;
    for(int i=n; i<(2*n); ++i) {
        for(int j=0; j<all_words; ++j){
            out[i].push_back(tmp);
        }
    }
    return out; // wrap(out)
}


class Fitch {
public:
    // initialize data 
    Fitch (RObject obj) {
        // std::vector< std::vector<int> > Y;
      weight = obj.attr("weight");
      nr = obj.attr("nr");
      NIBBLE = obj.attr("nc");
      NIBBLE_M1 = NIBBLE - 1;
      word = bits / NIBBLE; 
      rest = nr % word;
      n_words = (rest) ? ((nr / word) + 1) : (nr / word);
      EIGHTS = getEights(NIBBLE, bits);
      SEVENS = getSevens(NIBBLE, bits);
      ONES = getOnes(NIBBLE, bits);
      FIFTEEN = getFifteen(NIBBLE, bits);
      contr = contrast2uint(obj.attr("contrast"));
      X = readFitch(obj, contr, nr, bits, NIBBLE);
    }
    
    double pscore(IntegerMatrix orig)
    {   
      IntegerVector parent = orig( _, 0);
      IntegerVector child = orig( _, 1);
      uint64_t x, y, u; //, z;
      int j=0, pi=0, rc, lc, nl=parent.size();//max(parent);
      
      int unrooted = nl % 2;
      if(unrooted == 1){
        nl = nl-1;
      }
      //    std::vector<double> pvec(nl); 
      //  NumericVector pvec(nl);       
      double tmp=0.0, pars=0.0;
      for(int i=0; i<n_words; ++i){    
        j=0;
//        z = 0U;
        while(j < nl){ // (nl - 1L)){
          pi = parent[j] - 1; 
          lc = child[j] - 1;
          rc = child[j+1] - 1; 
          //      pvec[pi] = pvec[lc] + pvec[rc];
          x = X[lc][i];
          y = X[rc][i];
          u = ((((x & y & SEVENS) + SEVENS) | (x & y)) & EIGHTS) >> NIBBLE_M1;
//          z+= (u ^ ONES);
          X[pi][i] = (x & y) | ((x | y) & ((u + SEVENS) ^ EIGHTS));
//          if(z & EIGHTS){
            for(int k=0; k<word; ++k){
              pars += (u & 1U) ? 0.0 : weight[word * i + k]; //weight[(i << LENGTH_WORD) + j];
//              pars += (z & FIFTEEN) * weight[word * i + k];
              u >>= NIBBLE;
            }
//            z=0U;
//          }
          j+=2;
        }
        // unrooted case    
        if(unrooted){
          lc = child[nl] - 1;
          x = X[lc][i];
          y = X[pi][i];
          u = ((((x & y & SEVENS) + SEVENS) | (x & y)) & EIGHTS) >> NIBBLE_M1;
//          z+= (u ^ ONES);
          X[pi][i] = (x & y) | ((x | y) & ((u + SEVENS) ^ EIGHTS));
          for(int k=0; k<word; ++k){
            pars += (u & 1U) ? 0.0 : weight[word * i + k]; //weight[(i << LENGTH_WORD) + j];
            u >>= NIBBLE;
          }
          //      pvec[pi] += pvec[lc]; 
        }
      }  
      return pars;
    }

// minimal schneller  
// parallelize openMP intern or forks extern??
    double pscore2(IntegerMatrix orig)
    {   
        IntegerVector parent = orig( _, 0);
        IntegerVector child = orig( _, 1);
        uint64_t x, y, u;
        int j=0, pi=0, rc, lc, nl=parent.size();//max(parent);
        
        int unrooted = nl % 2;
        if(unrooted == 1){
            nl = nl-1;
        }
        //    std::vector<double> pvec(nl); 
        //  NumericVector pvec(nl);       
        double tmp=0.0, pars=0.0;
        j=0;
        while(j < nl){ // (nl - 1L)){
            pi = parent[j] - 1; 
            lc = child[j] - 1;
            rc = child[j+1] - 1; 
            for(int i=0; i<n_words; ++i){    
            
                //      pvec[pi] = pvec[lc] + pvec[rc];
                x = X[lc][i];
                y = X[rc][i];
                u = ((((x & y & SEVENS) + SEVENS) | (x & y)) & EIGHTS) >> NIBBLE_M1;
                X[pi][i] = (x & y) | ((x | y) & ((u + SEVENS) ^ EIGHTS));
                for(int k=0; k<word; ++k){
                    pars += (u & 1U) ? 0.0 : weight[word * i + k]; //weight[(i << LENGTH_WORD) + j];
                    u >>= NIBBLE;
                }
            }
            j+=2;
        }
        // unrooted case    
        if(unrooted){
            lc = child[nl] - 1;
            for(int i=0; i<n_words; ++i){
                x = X[lc][i];
                y = X[pi][i];
                u = ((((x & y & SEVENS) + SEVENS) | (x & y)) & EIGHTS) >> NIBBLE_M1;
                X[pi][i] = (x & y) | ((x | y) & ((u + SEVENS) ^ EIGHTS));
                for(int k=0; k<word; ++k){
                    pars += (u & 1U) ? 0.0 : weight[word * i + k]; //weight[(i << LENGTH_WORD) + j];
                    u >>= NIBBLE;
                }
            }    
                //      pvec[pi] += pvec[lc]; 
        }
        return pars;
    }
    

  int getNR(void){
    return nr;
  }

    NumericVector getWeight(void){
        return weight;
    }

    int getSeven(void){
        return SEVENS;
    }
    
    int getWord(void){
        return word;
    }

    int getWords(void){
        return n_words;
    }

    int getRest(void){
        return rest;
    }
    
    
private:
    std::vector< std::vector<uint64_t> > X;
    // std::vector< std::vector<int> > Y;
    NumericVector weight;
    int nr;
    int NIBBLE;
    int NIBBLE_M1;
    int word; 
    std::vector<uint64_t> contr;
    int n_words;
    int rest;
    uint64_t EIGHTS;
    uint64_t SEVENS;
    uint64_t ONES;
    uint64_t FIFTEEN;
};


RCPP_MODULE(Fitch_mod) {
    using namespace Rcpp;
    class_<Fitch>("Fitch")
        .constructor<RObject>("Default constructor")
//        .method("pscore", &Fitch::pscore)
        .property("nr", &Fitch::getNR)
        .property("weight", &Fitch::getWeight)
        .property("seven", &Fitch::getSeven)
        .property("word", &Fitch::getWord)
        .property("words", &Fitch::getWords)
//        .property("contr", &Fitch::getContr)
//        .property("X", &Fitch::getX)
        .method("pscore", &Fitch::pscore)
        .method("pscore2", &Fitch::pscore2)
    ;
}



/*

library(phangorn)
data(Laurasiatherian)
tree = NJ(dist.ml(Laurasiatherian)) 
f <- new(Fitch, Laurasiatherian)
f$pscore(tree$edge)
 
*/
