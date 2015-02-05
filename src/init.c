#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void countCycle(int *M, int *l, int *m, int *res);
void countCycle2(int *M, int *l, int *m, int *res);

/* called in fitch.R */
void fitch_free();
void fitch_init(int *data, int *m, int *n, double *weights, int *nr);
//void fitch8(int *dat, int *nr, int *pars, int *node, int *edge, int *nl, double *weight, double *pvec, double *pscore);
void fitchQuartet(int *index, int *n, int *nr, double *psc1, double *psc2, double *weight, double *res);
void fnhelp(int *node, int * edge, int *n, int *m, int *root, int *edge2, int *node2, int *pc);

/* called in parsimony.R */
void fitchTriplet(int *res, int *dat1, int *dat2, int *dat3, int *nr); 
void fitchTripletACC4(int *root, int *dat1, int *dat2, int *dat3, int *nr, double *p1, 
    double *p2, double *p3, double *weight, double *pars1, int *v1); 
void ACCTRAN2(int *dat, int *nr, double *pars, int *node, int *edge, int *nl, 
    double *weight, double *pvec, int *nTips);  
void ACCTRAN3(int *dat, int *nr, double *pars, int *node, int *edge, int *nl, 
    double *weight, double *pvec, int *nTips); 
void prepRooted(int *res, int *nr, int *kids);
void C_MPR(int *res, int *nr, int *parent, int *kids, int *nl); 
    

/* from distSeq.R */    
void distHamming(int *x, double *weight, int *nr, int *l, double *d);    
//void C_coph(SEXP children, SEXP tips, double *nh, int *nTips, int *lch, int *lkids, int *ltips, double *dm);    

/* from networx*/
//void neworder_cladewise(int *n, int *edge1, int *edge2, int *N, int *neworder);

/* from parsimony*/
void countMPR(double *res, int *dat1, int *dat2, int *nr, double *weight, int *external);
void reorder(int *from, int *to, int *n, int *sumNode,  int *neworder, int *root);

/* from phylo.R*/
void ll_free();
void ll_init(int *nr, int *nTips, int *nc, int *k);
void ll_free();
void ll_init2(int *data, double *weights, int *nr, int *nTips, int *nc, int *k);
//void cisort(int *x, int *y, int *a, int *b, int *res);
void moveLL(double *LL, double *child, double *P, int *nr, int *nc, double *tmp);

void out(double *d, double *r, int *n, int *k, int *l);
void fhm(double *v, int *n);
void giveIndex(int *left, int* right, int *ll, int *lr, int *n, int *res);

void nodeH(int *edge, int *node, double *el, int *l,  double *res);
void AllKids(int *children, int *parents, int *nTips, int *nNode, int *lp, int *kids, int *lkids, int *pkids);
//void C_bipHelp(int *parents, int *children, int *ntips, int *mp, int *l, int *ltips, int *ptips);
//void C_bip2(int *parents, int *children, int *ntips, int *mp, int *l, int *ltips, int *ptips, int *tips);
void C_cophenetic(int *children, int *parents, double *el, int *lp, int *m, int *nTips, int *nNode, double *res);


/* called in fitch.R */
SEXP FITCH(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP weight, SEXP mx, SEXP q);   
SEXP FITCH345(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP ps);
SEXP FITCHTRIP3(SEXP DAT3, SEXP nrx, SEXP edge, SEXP score, SEXP PS);
SEXP FNALL5(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP my, SEXP root);
// SEXP FNALL3(SEXP nrx, SEXP node, SEXP edge, SEXP node2, SEXP edge2, SEXP l, SEXP mx, SEXP my, SEXP q, SEXP pc);
SEXP FNALL_NNI(SEXP nrx, SEXP node, SEXP edge, SEXP l, SEXP mx, SEXP my, SEXP root);
SEXP AddOne(SEXP edge, SEXP tip, SEXP ind, SEXP l, SEXP m);
SEXP getData(SEXP n);

/* treemanipulation phangorn.c */
SEXP AllChildren(SEXP children, SEXP parent, SEXP M);
SEXP AllDesc(SEXP child, SEXP parent, SEXP M, SEXP NODE);

/* called in phylo.R */
SEXP C_bipart(SEXP parent, SEXP child, SEXP nTips, SEXP maxP); //, SEXP Nnode);
SEXP C_bip(SEXP parent, SEXP child, SEXP nTips, SEXP maxP);


SEXP sankoff3(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips);
SEXP sankoffMPR(SEXP dlist, SEXP plist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge);
SEXP C_rowMin(SEXP sdat, SEXP sn, SEXP sk);
SEXP pNodes(SEXP data, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge);
SEXP sankoffQuartet(SEXP dat, SEXP sn, SEXP scost, SEXP sk);

/* from distSeq.R */
SEXP PWI(SEXP LEFT, SEXP RIGHT, SEXP L, SEXP N, SEXP W, SEXP LI);

/* from networx*/
SEXP LogLik2(SEXP dlist, SEXP P, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP nTips, SEXP mNodes, SEXP contrast, SEXP nco);

/* from parsimony*/
//SEXP FNALL(SEXP dat, SEXP nrx, SEXP node, SEXP edge, SEXP node2, SEXP edge2, SEXP l, SEXP weight, SEXP mx, SEXP my, SEXP q, SEXP pc);

/* from phylo.R */
SEXP dist2spectra(SEXP dm, SEXP nx, SEXP ns);
SEXP getd2PM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getdPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getdPM2(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getd2PM2(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getPM2(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP getPM(SEXP eig, SEXP nc, SEXP el, SEXP w);
SEXP invSites(SEXP dlist, SEXP nr, SEXP nc, SEXP contrast, SEXP nco);
SEXP PML0(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N);
SEXP PML_NEW(SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N);
SEXP rowMax(SEXP sdat, SEXP sn, SEXP sk);
SEXP PML3(SEXP dlist, SEXP EL, SEXP W, SEXP G, SEXP NR, SEXP NC, SEXP K, SEXP eig, SEXP bf, SEXP node, SEXP edge, SEXP NTips, SEXP root, SEXP nco, SEXP contrast, SEXP N);
SEXP getLL(SEXP ax, SEXP bx, SEXP nrx, SEXP ncx, SEXP nTips);
SEXP getSCM(SEXP kk, SEXP nrx, SEXP nTips);
SEXP getM3(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc);
SEXP FS4(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP dad, SEXP child, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight, SEXP f0, SEXP retA, SEXP retB);
SEXP FS5(SEXP eig, SEXP nc, SEXP el, SEXP w, SEXP g, SEXP X, SEXP ld, SEXP nr, SEXP basefreq, SEXP weight, SEXP f0);
SEXP getDAD(SEXP dad, SEXP child, SEXP P, SEXP nr, SEXP nc);
SEXP getDAD2(SEXP dad, SEXP child, SEXP contrast, SEXP P, SEXP nr, SEXP nc, SEXP nco);
SEXP getPrep(SEXP dad, SEXP child, SEXP eve, SEXP evi, SEXP nr, SEXP nc);
SEXP getPrep2(SEXP dad, SEXP child, SEXP contrast, SEXP evi, SEXP nr, SEXP nc, SEXP nco);
SEXP getXX(SEXP nr, SEXP nTips);
/* from sankoff.R */
SEXP sankoff3B(SEXP dlist, SEXP scost, SEXP nr, SEXP nc, SEXP node, SEXP edge, SEXP mNodes, SEXP tips, SEXP contrast, SEXP nrs);

R_CallMethodDef callMethods[] = {
{"FITCH", (DL_FUNC) &FITCH, 8},
{"FITCH345", (DL_FUNC) &FITCH345, 6},
{"FITCHTRIP3", (DL_FUNC) &FITCHTRIP3, 5},
{"FNALL5", (DL_FUNC) &FNALL5, 7},    
{"FNALL_NNI", (DL_FUNC) &FNALL_NNI, 7}, 
{"getData", (DL_FUNC) &getData, 2}, 
{"AddOne", (DL_FUNC) &AddOne, 5},
{"AllChildren", (DL_FUNC) &AllChildren, 3},
{"C_bipart", (DL_FUNC) &C_bipart, 4},
{"C_bip", (DL_FUNC) &C_bip, 4},
{"sankoff3", (DL_FUNC) &sankoff3, 8},   
{"sankoffMPR", (DL_FUNC) &sankoffMPR, 7},  
{"C_rowMin", (DL_FUNC) &C_rowMin, 3},
{"pNodes", (DL_FUNC) &pNodes, 6},
{"sankoffQuartet", (DL_FUNC) &sankoffQuartet, 4},
{"PWI", (DL_FUNC) &PWI, 6},
{"LogLik2", (DL_FUNC) &LogLik2, 10},
{"dist2spectra", (DL_FUNC) &dist2spectra, 3},
{"getdPM", (DL_FUNC) &getdPM, 4},
{"getd2PM", (DL_FUNC) &getd2PM, 4},
{"getdPM2", (DL_FUNC) &getdPM2, 4},
{"getd2PM2", (DL_FUNC) &getd2PM2, 4},
{"getPM2", (DL_FUNC) &getPM2, 4},
{"getPM", (DL_FUNC) &getPM, 4},
{"invSites", (DL_FUNC) &invSites, 5},
{"PML0", (DL_FUNC) &PML0, 16},
{"PML3", (DL_FUNC) &PML3, 16},
{"PML_NEW", (DL_FUNC) &PML_NEW, 15},
{"rowMax", (DL_FUNC) &rowMax, 3},
{"getSCM", (DL_FUNC) &getSCM, 3},
{"getM3", (DL_FUNC) &getM3, 5},
{"FS4", (DL_FUNC) &FS4, 15},
{"FS5", (DL_FUNC) &FS5, 11},
{"getDAD", (DL_FUNC) &getDAD, 5},
{"getDAD2", (DL_FUNC) &getDAD2, 7},
{"getPrep", (DL_FUNC) &getPrep, 6},
{"getPrep2", (DL_FUNC) &getPrep2, 7},
{"sankoff3B", (DL_FUNC) &sankoff3B, 10},
{"AllDesc", (DL_FUNC) &AllDesc, 4},
{"getXX", (DL_FUNC) &getXX, 2},
{NULL, NULL, 0}
};
//{"FNALL3", (DL_FUNC) &FNALL5, 10},
//{"FNALL", (DL_FUNC) &FNALL, 12},
//{"coph", (DL_FUNC) &coph, 7},
     

R_CMethodDef cMethods[] = { 
{"countCycle", (DL_FUNC) &countCycle, 4},
{"countCycle2", (DL_FUNC) &countCycle2, 4},
{"fitch_free", (DL_FUNC) &fitch_free, 0},  
{"fitch_init", (DL_FUNC) &fitch_init, 5}, 
{"fnhelp", (DL_FUNC) &fnhelp, 8},
{"C_fhm", (DL_FUNC) &fhm, 2},
{"out", (DL_FUNC) &out, 5},
{"giveIndex", (DL_FUNC) &giveIndex, 6},
{"fitchQuartet", (DL_FUNC) &fitchQuartet, 7},
{"fitchTriplet", (DL_FUNC) &fitchTriplet, 5},
{"fitchTripletACC4", (DL_FUNC) &fitchTripletACC4, 11},
{"ACCTRAN2", (DL_FUNC) &ACCTRAN2, 9},
{"ACCTRAN3", (DL_FUNC) &ACCTRAN3, 9},
{"prepRooted", (DL_FUNC) &prepRooted, 3},
{"C_MPR", (DL_FUNC) &C_MPR, 5},
{"distHamming", (DL_FUNC) &distHamming, 5},
{"C_reorder", (DL_FUNC) &reorder, 6},
{"countMPR", (DL_FUNC) &countMPR, 6},
{"ll_free", (DL_FUNC) &ll_free, 0},
{"ll_init", (DL_FUNC) &ll_init, 4},
{"ll_free2", (DL_FUNC) &ll_free, 0},
{"ll_init2", (DL_FUNC) &ll_init2, 6},
{"moveLL", (DL_FUNC) &moveLL, 6},
{"nodeH", (DL_FUNC) &nodeH, 5},
{"AllKids", (DL_FUNC) &AllKids, 8},
{"C_cophenetic", (DL_FUNC) &C_cophenetic, 8},
{NULL, NULL, 0}
};
// fitchQuartet raus eins weniger
//{"fitch8", (DL_FUNC) &fitch8, 9},
//{"C_coph", (DL_FUNC) &C_coph, 10},
//{"C_bipHelp", (DL_FUNC) &C_bipHelp, 7},
//{"C_bip2", (DL_FUNC) &C_bip2, 8},
//{"C_cisort", (DL_FUNC) &cisort, 5},

void R_init_phangorn(DllInfo *info){
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
    
    

