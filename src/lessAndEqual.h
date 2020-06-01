#ifndef R_R_H
#include <R.h>
#endif

#ifndef R_INTERNALS_H_
#include <Rinternals.h>
#endif

#include <cstring>

/* 
 NOTE: R_NaString is a different SEXP than mkChar("NA"), but holding the same string "NA". 
 We will treat R_NaString to be smaller than every usual string, including mkChar("NA"). 
 Real NaN becomes mkChar("NaN") by as.character();
 Real -Inf becomes mkChar("-Inf") by as.character();
 Real  Inf becomes mkChar("Inf") by as.character();
 */

class CharSEXP{
public:
    SEXP sexp;
    inline bool valid() {return( TYPEOF(sexp) == CHARSXP );}
    
    CharSEXP(SEXP x)
    {
        if (TYPEOF(x) == CHARSXP) sexp = x;
        else error("CharSEXP should be initialized with a CHARSXP type object");
    }
    
    CharSEXP()
    {
        sexp = R_NaString; 
    }
    
    friend inline bool operator< (const CharSEXP& lhs, const CharSEXP& rhs) 
    {
        if (lhs.sexp == R_NaString) return( rhs.sexp != R_NaString );
        if (rhs.sexp == R_NaString) return(false);
        return( 
            strcmp(const_cast<const char*>(CHAR(lhs.sexp)), const_cast<const char*>(CHAR(rhs.sexp)) )<0
        ); 
    }
    
    friend inline bool operator== (const CharSEXP& lhs, const CharSEXP& rhs) 
    {	
        return (lhs.sexp == rhs.sexp);  // R CHARSXP objects are cached (only one copy per string)
    }
    
};


/* for general T where operator< and operator== have been implemented; 
 * this is helpful for 
 *  (1) integers where NA's is just a special integer value;
 *  (2) unsigned char, where NA's is converted to 00;
 *  (3) CharSEXP, where NA's are properly handled by operator<.
 */
template <typename T>
 bool lessThan (const T& lhs, const T& rhs) { return lhs < rhs;}
template <typename T>
 bool equalTo (const T& lhs, const T& rhs) { return lhs == rhs;}


/* double Assumptions: NaN < NA_real_ < -Inf < Finite numbers < Inf */
template <>
inline bool lessThan<double>(const double& lhs, const double& rhs) 
{
    if (R_FINITE(lhs) && R_FINITE(rhs)) return lhs< rhs; // probably the most common case (both finite)
    
    bool rhsTest = R_IsNaN(rhs);        // rhs = NaN    
    if (R_IsNaN(lhs)) return !rhsTest; // lhs = NaN
    
    rhsTest = rhsTest || ISNA(rhs);     // rhs <= NA_real_
    if (ISNA(lhs)) return !rhsTest;     // lhs = NA
    
    rhsTest = rhsTest || (rhs == R_NegInf); // rhs <= -Inf
    if (lhs == R_NegInf) return !rhsTest;   // lhs = -Inf
    
    if(rhsTest) return false;       // lhs is finite or +Inf but rhs <= -Inf 
    return R_FINITE(lhs);           // lhs is finite or +Inf but rhs is +Inf
}

template <>
inline bool equalTo<double> (const double& lhs, const double& rhs) 
{return(
    (lhs == rhs) ||
    (ISNA(lhs) && ISNA(rhs)) ||
    (R_IsNaN(lhs) && R_IsNaN(rhs))
);}


