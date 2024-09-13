#include "rcSet.h"


vecMap<int>     		intVecMap;
vecMap<double> 			doubleVecMap;
vecMap<CharSEXP>		charsexpVecMap;
vecMap<unsigned char>	rawVecMap; 		// Rbyte is an alias of unsigned char

extern "C" {

SEXP grpDupAtomMat(SEXP x, SEXP MARGIN, SEXP fromLast)
{/* returns an integer vector of duplicated rows of numeric matrix x */
    SEXP out;
	int* dim;
    int nGrps=0;
	dim=INTEGER(Rf_getAttrib(x, R_DimSymbol));
	out = PROTECT(Rf_allocVector(INTSXP, dim[*INTEGER(MARGIN)-1]));

	switch (TYPEOF(x)) {
		case REALSXP:
			nGrps = doubleVecMap.grpDuplicatedMat	(REAL(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			break;
		case INTSXP:  // factor type is also covered here
			// if(!inherits(x, "factor"))
				nGrps = intVecMap.grpDuplicatedMat	(INTEGER(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			// else {;}
			break;
		case LGLSXP:
			nGrps = intVecMap.grpDuplicatedMat	(LOGICAL(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			break;
		case STRSXP: {
			CharSEXP* charSexpPtr = new CharSEXP [ dim[0]*dim[1] ];
			for(int i=dim[0]*dim[1]-1; i>=0; --i)
				charSexpPtr[i].sexp = STRING_ELT(x, i);

			nGrps = charsexpVecMap.grpDuplicatedMat	(charSexpPtr, dim, dim+1, INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );

			delete[] charSexpPtr;
			break;
		}
		case RAWSXP:
			nGrps = rawVecMap.grpDuplicatedMat	(RAW(x), dim, dim+1,  INTEGER(out), *INTEGER(MARGIN)==1, (bool)(*(LOGICAL(fromLast))) );
			break;
//		default:
//			error("C function 'grpDupAtomMat' only accepts REALSXP, LGLSXP, INTSXP and STRSXP");
	}

    SEXP nLevels;
    nLevels = PROTECT(Rf_allocVector(INTSXP, 1));
    INTEGER(nLevels)[0] = nGrps;
    Rf_setAttrib(out, Rf_install("nlevels"), nLevels);
    UNPROTECT(2);
	return out;
}

}

