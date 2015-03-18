#include <R.h>
#include <Rinternals.h>

// The initial code defining and initialising the translation table:
//
//"a" "r" "n" "d" "c" "q" "e" "g" "h" "i" "l" "k" "m" "f" "p" "s" "t" "w" "y" "v" "b" "z" "x"
//  "-" "?"
//
//    for (i = 0; i < 122; i++) tab_trans[i] = 0x00;
//
//	tab_trans[65] = 0x88; /* A */
//	tab_trans[71] = 0x48; /* G */
//	tab_trans[67] = 0x28; /* C */
// 	tab_trans[84] = 0x18; /* T */
// 	tab_trans[82] = 0xc0; /* R */
// 	tab_trans[77] = 0xa0; /* M */
// 	tab_trans[87] = 0x90; /* W */
// 	tab_trans[83] = 0x60; /* S */
// 	tab_trans[75] = 0x50; /* K */
// 	tab_trans[89] = 0x30; /* Y */
// 	tab_trans[86] = 0xe0; /* V */
// 	tab_trans[72] = 0xb0; /* H */
// 	tab_trans[68] = 0xd0; /* D */
//  	tab_trans[66] = 0x70; /* B */
// 	tab_trans[78] = 0xf0; /* N */
//
//	tab_trans[97] = 0x88; /* a */
//	tab_trans[103] = 0x48; /* g */
//	tab_trans[99] = 0x28; /* c */
// 	tab_trans[116] = 0x18; /* t */
// 	tab_trans[114] = 0xc0; /* r */
// 	tab_trans[109] = 0xa0; /* m */
// 	tab_trans[119] = 0x90; /* w */
// 	tab_trans[115] = 0x60; /* s */
// 	tab_trans[107] = 0x50; /* k */
// 	tab_trans[121] = 0x30; /* y */
// 	tab_trans[118] = 0xe0; /* v */
// 	tab_trans[104] = 0xb0; /* h */
// 	tab_trans[100] = 0xd0; /* d */
//  	tab_trans[98] = 0x70; /* b */
// 	tab_trans[110] = 0xf0; /* n */
//
//  	tab_trans[45] = 0x04; /* - */
//  	tab_trans[63] = 0x02; /* ? */


static const int tab_trans2[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 0-9 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 10-19 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 20-29 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 30-39 */
	0, 0, 0, 0, 0, 24, 0, 0, 0, 0, /* 40-49 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 50-59 */
	0, 0, 0, 25, 0, 1, 21, 5, 4, 7, /* 60-69 */
	14, 8, 9, 10, 0, 12, 11, 13, 3, 0, /* 70-79 */
	15, 6, 2, 16, 17, 0, 20, 18, 23, 19, /* 80-89 */
	22, 0, 0, 0, 0, 0, 0, 1, 21, 5, /* 90-99 */
	4, 7, 14, 8, 9, 10, 0, 12, 11, 13, /* 100-109 */
	3, 0, 15, 6, 2, 16, 17, 0, 20, 18, /* 110-119 */
	23, 19, 22, 0, 0, 0, 0, 0, 0, 0, /* 120-129 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 130-139 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /* 140-149 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  /* 150-159 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 160-169 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 170-179 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 180-189 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 190-199 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 200-209 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 210-219 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 220-229 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 230-239 */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, /* 240-249 */
	0, 0, 0, 0, 0, 0}; /* 250-255 */


static const unsigned char hook = 0x3e;
static const unsigned char lineFeed = 0x0a;
/* static const unsigned char space = 0x20; */


// needs buffer seq
// + buffer names
SEXP rawStream2phyDat(SEXP x)
{
	int N, i, j, k, n, startOfSeq;
	unsigned char *xr, *bufferNames;
    int  *rseq, *buffer, tmp;
	SEXP obj, nms, seq;

	PROTECT(x = coerceVector(x, RAWSXP));
	N = LENGTH(x);
	xr = RAW(x);

/* do a 1st pass to find the number of sequences

   this code should be robust to '>' present inside
   a label or in the header text before the sequences */

	n = j = 0; /* use j as a flag */
	if (xr[0] == hook) {
		j = 1;
		startOfSeq = 0;
	}
	i = 1;
	for (i = 1; i < N; i++) {
		if (j && xr[i] == lineFeed) {
			n++;
			j = 0;
		} else if (xr[i] == hook) {
			if (!n) startOfSeq = i;
			j = 1;
		}
	}

	PROTECT(obj = allocVector(VECSXP, n));
	PROTECT(nms = allocVector(STRSXP, n));

/* Refine the way the size of the buffer is set? */
	buffer = (int *)R_alloc(N, sizeof(int *));
    bufferNames = (unsigned char *)R_alloc(N, sizeof(unsigned char *));

	i = startOfSeq;
	j = 0; /* gives the index of the sequence */
	while (i < N) {
		/* 1st read the label... */
		i++;
		k = 0;
		while (xr[i] != lineFeed) bufferNames[k++] = xr[i++];
		bufferNames[k] = '\0';
		SET_STRING_ELT(nms, j, mkChar((char *)bufferNames));
		/* ... then read the sequence */
		n = 0;
		while (i < N && xr[i] != hook) {
			tmp = tab_trans2[xr[i++]];
/* If we are sure that the FASTA file is correct (ie, the sequence on
   a single line and only the IUAPC code (plus '-' and '?') is used,
   then the following check would not be needed; additionally the size
   of tab_trans could be restriced to 0-121. This check has the
   advantage that all invalid characters are simply ignored without
   causing error -- except if '>' occurs in the middle of a sequence. */
			if(tmp) buffer[n++] = tmp;
		}
		PROTECT(seq = allocVector(INTSXP, n));
		rseq = INTEGER(seq);
		for (k = 0; k < n; k++) rseq[k] = buffer[k];
		SET_VECTOR_ELT(obj, j, seq);
		UNPROTECT(1);
		j++;
	}
	setAttrib(obj, R_NamesSymbol, nms);
	UNPROTECT(3);
	return obj;
}
