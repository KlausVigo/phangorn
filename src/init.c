#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    phangorn_allDescCPP
    phangorn_bipartCPP
    phangorn_bipCPP
    phangorn_allChildrenCPP

  Most likely possible values need to be added below.
*/

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void ACCTRAN2(void *, void *, void *, void *, void *, void *);
extern void ACCTRAN3(void *, void *, void *, void *, void *, void *, void *, void *);
extern void C_fhm(void *, void *);
//extern void C_reorder(void *, void *, void *, void *, void *, void *);
extern void countCycle(void *, void *, void *, void *);
extern void countCycle2(void *, void *, void *, void *);
extern void distHamming(void *, void *, void *, void *, void *);
extern void fitch_free();
extern void fitch_init(void *, void *, void *, void *, void *);
extern void fitchQuartet(void *, void *, void *, void *, void *, void *, void *);
extern void fitchTriplet(void *, void *, void *, void *, void *);
extern void fitchTripletACC4(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void giveIndex(void *, void *, void *, void *, void *, void *);
extern void ll_free();
extern void ll_init(void *, void *, void *, void *);
//extern void nodeH(void *, void *, void *, void *, void *);
extern void out(void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP AddOnes(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_rowMin(SEXP, SEXP, SEXP);
extern SEXP C_sprdist(SEXP, SEXP, SEXP);
extern SEXP dist2spectra(SEXP, SEXP, SEXP);
extern SEXP FITCH(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FITCH345(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FITCHTRIP3(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FNALL_NNI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FNALL5(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FNALL6(SEXP, SEXP, SEXP, SEXP);
extern SEXP FS4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP FS5(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getd2PM(SEXP, SEXP, SEXP, SEXP);
extern SEXP getd2PM2(SEXP, SEXP, SEXP, SEXP);
extern SEXP getDAD(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getDAD2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getdPM(SEXP, SEXP, SEXP, SEXP);
extern SEXP getdPM2(SEXP, SEXP, SEXP, SEXP);
extern SEXP getPM(SEXP, SEXP, SEXP, SEXP);
extern SEXP getPrep(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP getPrep2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP invSites(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LogLik2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optE(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optQrtt(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _phangorn_bipartCPP(SEXP, SEXP);
extern SEXP _phangorn_bipCPP(SEXP, SEXP);
extern SEXP PML0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PML3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PML4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP pNodes(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP PWI(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// extern SEXP rawStream2phyDat(SEXP);
extern SEXP rowMax(SEXP, SEXP, SEXP);
extern SEXP sankoff3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sankoff3B(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sankoffMPR(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sankoffQuartet(SEXP, SEXP, SEXP, SEXP);
extern SEXP _phangorn_allDescCPP(SEXP, SEXP);
extern SEXP _phangorn_allChildrenCPP(SEXP);
//extern SEXP _phangorn_allSiblingsCPP(SEXP);
//extern SEXP _phangorn_preorder(SEXP);
extern SEXP _phangorn_p2dna(SEXP, SEXP);
extern SEXP _phangorn_threshStateC(SEXP, SEXP);
extern SEXP _phangorn_node_height_cpp(SEXP, SEXP, SEXP);
extern SEXP _phangorn_cophenetic_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"ACCTRAN2",         (DL_FUNC) &ACCTRAN2,          6},
    {"ACCTRAN3",         (DL_FUNC) &ACCTRAN3,          8},
    {"C_fhm",            (DL_FUNC) &C_fhm,             2},
    {"countCycle",       (DL_FUNC) &countCycle,        4},
    {"countCycle2",      (DL_FUNC) &countCycle2,       4},
    {"distHamming",      (DL_FUNC) &distHamming,       5},
    {"fitch_free",       (DL_FUNC) &fitch_free,        0},
    {"fitch_init",       (DL_FUNC) &fitch_init,        5},
    {"fitchQuartet",     (DL_FUNC) &fitchQuartet,      7},
    {"fitchTriplet",     (DL_FUNC) &fitchTriplet,      5},
    {"fitchTripletACC4", (DL_FUNC) &fitchTripletACC4, 11},
    {"giveIndex",        (DL_FUNC) &giveIndex,         6},
    {"ll_free",          (DL_FUNC) &ll_free,           0},
    {"ll_init",          (DL_FUNC) &ll_init,           4},
    {"out",              (DL_FUNC) &out,               5},
    {NULL, NULL, 0}
};
//    {"nodeH",            (DL_FUNC) &nodeH,             5},

static const R_CallMethodDef CallEntries[] = {
    {"AddOnes",            (DL_FUNC) &AddOnes,             5},
    {"C_rowMin",           (DL_FUNC) &C_rowMin,            3},
    {"C_sprdist",          (DL_FUNC) &C_sprdist,           3},
    {"dist2spectra",       (DL_FUNC) &dist2spectra,        3},
    {"FITCH",              (DL_FUNC) &FITCH,               8},
    {"FITCH345",           (DL_FUNC) &FITCH345,            6},
    {"FITCHTRIP3",         (DL_FUNC) &FITCHTRIP3,          5},
    {"FNALL_NNI",          (DL_FUNC) &FNALL_NNI,           7},
    {"FNALL5",             (DL_FUNC) &FNALL5,              7},
    {"FNALL6",             (DL_FUNC) &FNALL6,              4},
    {"FS4",                (DL_FUNC) &FS4,                14},
    {"FS5",                (DL_FUNC) &FS5,                10},
    {"getd2PM",            (DL_FUNC) &getd2PM,             4},
    {"getd2PM2",           (DL_FUNC) &getd2PM2,            4},
    {"getDAD",             (DL_FUNC) &getDAD,              5},
    {"getDAD2",            (DL_FUNC) &getDAD2,             7},
    {"getdPM",             (DL_FUNC) &getdPM,              4},
    {"getdPM2",            (DL_FUNC) &getdPM2,             4},
    {"getPM",              (DL_FUNC) &getPM,               4},
    {"getPrep",            (DL_FUNC) &getPrep,             6},
    {"getPrep2",           (DL_FUNC) &getPrep2,            7},
    {"invSites",           (DL_FUNC) &invSites,            5},
    {"LogLik2",            (DL_FUNC) &LogLik2,            10},
    {"optE",               (DL_FUNC) &optE,               17},
    {"optQrtt",            (DL_FUNC) &optQrtt,            16},
    {"_phangorn_bipartCPP", (DL_FUNC) &_phangorn_bipartCPP,  2},
    {"_phangorn_bipCPP",    (DL_FUNC) &_phangorn_bipCPP,     2},
    {"PML0",               (DL_FUNC) &PML0,               14},
    {"PML3",               (DL_FUNC) &PML3,               14},
    {"PML4",               (DL_FUNC) &PML4,               15},
    {"pNodes",             (DL_FUNC) &pNodes,              6},
    {"PWI",                (DL_FUNC) &PWI,                 6},
//    {"rawStream2phyDat",   (DL_FUNC) &rawStream2phyDat,    1},
    {"rowMax",             (DL_FUNC) &rowMax,              3},
    {"sankoff3",           (DL_FUNC) &sankoff3,            8},
    {"sankoff3B",          (DL_FUNC) &sankoff3B,          10},
    {"sankoffMPR",         (DL_FUNC) &sankoffMPR,          7},
    {"sankoffQuartet",     (DL_FUNC) &sankoffQuartet,      4},
    {"_phangorn_allDescCPP",       (DL_FUNC) &_phangorn_allDescCPP,        2},
    {"_phangorn_allChildrenCPP",   (DL_FUNC) &_phangorn_allChildrenCPP,    1},
//    {"_phangorn_allSiblingsCPP",   (DL_FUNC) &_phangorn_allSiblingsCPP,    1},
//    {"_phangorn_preorder",         (DL_FUNC) &_phangorn_preorder,          1},
    {"_phangorn_p2dna",     (DL_FUNC) &_phangorn_p2dna,   2},
    {"_phangorn_threshStateC",     (DL_FUNC) &_phangorn_threshStateC, 2},
    {"_phangorn_node_height_cpp",     (DL_FUNC) &_phangorn_node_height_cpp, 3},
    {"_phangorn_cophenetic_cpp",     (DL_FUNC) &_phangorn_cophenetic_cpp, 4},
    {NULL, NULL, 0}
};

void R_init_phangorn(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
