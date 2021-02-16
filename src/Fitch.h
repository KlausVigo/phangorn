#include <Rcpp.h>

#ifndef _FITCH_H_
#define _FITCH_H_



#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

#ifdef __GNUC__
#define GNUC_PREREQ(x, y) \
(__GNUC__ > x || (__GNUC__ == x && __GNUC_MINOR__ >= y))
#else
#define GNUC_PREREQ(x, y) 0
#endif

#if GNUC_PREREQ(4, 2) || \
    __has_builtin(__builtin_popcount)
#define HAVE_BUILTIN_POPCOUNT
#endif


using namespace Rcpp;


// Counting bits set, Brian Kernighan's way
// from Bit Twiddling Hacks by Sean Eron Anderson seander@cs.stanford.edu
// count the number of bits set in v
// fast if few bits set, what we assume
static inline int bitCount32(uint32_t v){
  int c;
  for (c = 0; v; c++){
    v &= v - 1; // clear the least significant bit set
  }
  return c;
}

static inline int bitCount64(uint64_t v){
  int c;
  for (c = 0; v; c++){
    v &= v - 1; // clear the least significant bit set
  }
  return c;
}


#if defined(HAVE_BUILTIN_POPCOUNT)

static inline int popcnt64(uint64_t x)
{
  return __builtin_popcountll(x);
}

static inline int popcnt32(uint32_t x)
{
  return __builtin_popcount(x);
}

#else

static inline int popcnt64(uint64_t x)
{
  return bitCount64(x);
}

static inline int popcnt32(uint32_t x)
{
  return bitCount32(x);
}

#endif



#define BIT_SIZE 64

std::vector< std::vector<uint64_t> > readFitch(const List &xlist,
             IntegerMatrix contr, int nSeq, int nChar, int nStates,
             int nBits, int m);

// IntegerMatrix preorder(const IntegerMatrix & edge, int nTips);


class Fitch {
public:
  // Fitch(Robject)
  Fitch (RObject obj, int w1, int m) {
    weight = obj.attr("weight");
    nChar = (int) obj.attr("nr");
    if( (nChar % BIT_SIZE) ){
      for(int i= (nChar % BIT_SIZE); i<BIT_SIZE; ++i) weight.push_back(0.0);
    }
    nStates = (int) obj.attr("nc");
    p0 = obj.attr("p0");
    wBits = (w1 / BIT_SIZE) + (w1 % BIT_SIZE != 0);
    nBits = (nChar / BIT_SIZE) + (nChar % BIT_SIZE != 0);
    IntegerMatrix contr = obj.attr("contrast");
    Rcpp::List xlist(obj);
    nSeq = xlist.size();
    X = readFitch(xlist, contr, nSeq, nChar, nStates, nBits, m);
  }

  int getNR(void){ return nChar; }
  NumericVector getWeight(void){ return weight; }
  int getP0(void){ return p0; }
  int getnSeq(void){ return nSeq; }
  int getnBits(void){ return nBits; }

  std::vector< std::vector<uint64_t> > X;
  IntegerVector pscore_nodes;
  NumericVector weight; // Integer??
  int nChar;
  int nSeq;
  int nStates;
  int nBits;
  int wBits;
  int m;
  int p0;
};



#endif // _FITCH_H_



