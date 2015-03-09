#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>

#define MAXG 100
#define MAXW 10

#include "badpairs.h"
#include "admutils.h"
#include "regsubs.h"  
#include "egsubs.h"  


int loadindx(Indiv **xindlist, int *xindex, Indiv **indivmarkers, int numindivs) ;
int loadsnpx(SNP **xsnplist, SNP **snpmarkers, int numsnps, Indiv **indivmarkers) ;
void loadxdataind(double *xrow, SNP **snplist, int ind,  int ncols) ;           
void fixxrow(double *xrow, double *xmean, double *xfancy, int len)  ;
void dofancy(double *cc, int n, double *fancy) ;
int vadjust(double *rr, int n, double *pmean)  ;
void getcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
void getcolx(double *xcol, SNP *cupt, int *xindex, 
  int nrows, int col, double *xmean, double *xfancy)  ;          
void putcol(double *cc, double *xdata, int col, int nrows, int ncols)  ;
double dottest(char *sss, double *vec, char **eglist, int numeg, int *xtypes, int len) ;
double yll(double x1, double x2, double xlen) ;
void calcmean(double *wmean, double *vec, int len, int *xtypes, int numeg)  ;
void getrawcol(int *rawcol, SNP *cupt, int *xindex, int nrows) ;
void getrawcolx(int **ccc, SNP *cupt, int *xindex, int nrows, Indiv **indm) ;

void setmiss(SNP **snpm, int numsnps)  ;

void fixrho(double *a, int n) ;
void printdiag(double *a, int n) ;


double dofst(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, int mode) ;

double fstcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) ;

double divcol(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) ;

double fst(SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int type1, int type2, double *psd, int mode) ;

double dofstx(double *fstans, double *fstsd, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg) ;

void fstcolinb(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg) ;

void fstcolyy(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg) ;

double fstcoly(double *estn, double *estd, SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2) ;

double fstx(SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int type1, int type2, double *psd) ;

void
 setplimit(Indiv **indivmarkers, int numindivs, 
 char **eglist, int numeg, int plimit)  ;

void loadzdata(double **zdata,  SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg, int *ncolx, int *tagnums)  ;

void getpdata(int *rawcol, double *pm, double *pn, int *xtypes, int nrows, int numeg) ;

void getrscore(double *rscore, double *rho, double **zz, 
  int ncols, int a, int b, int c, int d, int numeg, int *blabels, int nblocks)   ;

double qcorr(double **zz, double *rho, 
  int ncols, int a, int b, int c, int d, int numeg, int *blabels, int nblocks)   ;
void xcopy(int rp[4], int a , int b, int c, int d)  ;
void settsc(int tpat[3][4], double tscore[3], int rpat[3][4], double rscore[3]) ;
void printsc(int pat[3][4], double tscore[3], char **eglist, double ymin)  ;
double dohzg(double *top, double *bot, SNP **xsnplist, int *xindex, int *xtypes, 
   int nrows, int ncols, int numeg)  ; 

void dohzgjack(double *fstest, double *fstsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int *bcols, int nblocks)  ;

void gethscore(double *hscore, double *scores, 
  int a, int b, int c, int d, int numeg) ;

double qhdiff(double *scores,  int a, int b, int c, int d, int numeg)  ;
void setblocks(int *block, int *bsize, int *nblock, SNP **snpm, int numsnps, double blocklen)  ;
int numblocks(SNP **snpm, int numsnps, double blocklen)  ;
void setmgpos(SNP **snpm, int numsnps, double *maxgdis)  ;
void setgfromp(SNP **snpm, int numsnps)   ;
void wjackest(double *est, double *sig, double mean, double *jmean, double *jwt, int n)  ;
void wjackvest(double *vest, double *var, int d, double *mean, double **jmean, double *jwt, int g)  ;
void corrwjack(double *xrho, double *xsig, double *z1, double *z2, int n, int *bcols, int nblocks);
double crho(double *stats)  ;

void ndfst5(double *zzest, double *zzsig, double **zn, double **zd, int ncols, int *bcols, int nblocks) ;
void regestit(double *ans, double *xn, double *xd) ;

void setwt(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, int nrows, 
  int *xindex, int *xtypes, char * outpop, char **eglist, int numeg) ;
void countg(int *rawcol, int **cc, int *xtypes, int n, int ntypes)   ;
void dohzgjack(double *hest, double *hsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int *bcols, int nblocks)   ;

void setbcols(SNP **xsnplist, int ncols, int *bcols)  ;
double
dofstnum(double *fst, double *fstnum, double *fstsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg,  int nblocks)  ;

double doinbreed(double *inb, double *inbest, double *inbsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, Indiv **indivmarkers)     ;

double
dofstnumx(double *fst, double *fstnum, double *fstsig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg,  int nblocks, Indiv **indm, int fstmode)  ;

void
dof3(double *f3, double *f3sig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale, int mode)  ;

void
dof4(double *f4, double *f4sig, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale, int mode)  ;

void f3y(double *estn,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3) ;

void f4y(double *estn,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3, int type4) ;

void f3sc(double *estn,  double *estd, SNP *cupt, Indiv **indm,  
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3) ;

void f2sc(double *estn,  double *estd, SNP *cupt, Indiv **indm,  
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3) ;

void f4yx(double *estn,  SNP *cupt, Indiv **indm,  
  int *xindex, int *xtypes, int nrows, int type1, int type2, int type3, int type4) ;

void f3yy(double *estmat,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg) ;  

int f3yyx(double *estmat,  SNP *cupt, 
  int *xindex, int *xtypes, int nrows, int numeg, Indiv **indm) ;  

void countpops(int ***counts, SNP **xsnplist, int *xindex, int *xtypes, int nrows, int ncols) ;

double doadmlin(double *jest, double *jsig, double *zlin, double *var, 
 SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale, Indiv **indm)  ;

double estmix(double *z, double *f3, int n) ;

void bump2(double *x, int a, int b, int n, double val)  ;
double dump2(double *x, int a, int b, int n)  ;
void bump3(double *x, int a, int b, int c, int n, double val)   ;
double dump3(double *x, int a, int b, int c, int n)  ;

void bump4(double *x, int a, int b, int c, int d, int n, double val)   ;
void bump4x(double *x, int a, int b, int c, int d, int n, double val)  ; // all 4 images
void set4x(double *x, int a, int b, int c, int d, int n, double val)   ;
void set4(double *x, int a, int b, int c, int d, int n, double val)  ; // all 4 images

double dump4(double *x, int a, int b, int c, int d, int n)  ;
double ff3val(double *ff3, int a, int b, int c, int n)  ;

// graph stuff
int loadgraph(char *readit, char ***peglist)   ;
void getgmix(double **vmix, int *lmix, int *nmix) ;
void putgmix(double **vmix) ; 
void getpwts(double *pwts, int *nrows, int *nedge) ;
void getenames(char **enames) ; 
void setsimp(double *ww, int n) ;
int edgenum(char *edgename) ;
void addvertex(char *vertname)  ;   
void supergetvnames(char **vnames, int *xvlist, int nxvlist)  ;
void superalloc(int **xv, int ***xe, int ***adv,  int **aedge) ;
void supersetup(int *xvlist, int *nxvlist, int **xelist, int *nxelist, int **admixv, int *admixedge, int *nxalist)  ;
void supereglist(int *eelist) ;
void supergetvar(double *svar, 
  int *xvlist, int nxvlist, int **xelist, int nxelist, int **admixv, int *admixedge, int nxalist)  ;
void superputvals(double **admixw, double *elen, 
  int *xvlist, int nxvlist, int **xelist, int nxelist, int **admixv, int *admixedge, int nxalist)  ;
void superest(double *xmean, double *xvar, double *svar, double *yobs,  int nxvlist, int *elist) ;
void supergetvals(double **admixw, double *elen, 
  int *xvlist, int nxvlist, int **xelist, int nxelist, int **admixv, int *admixedge, int nxalist)  ;
void superreest(double *s2, 
  int *xvlist, int nxvlist, int **xelist, int nxelist, int **admixv, int *admixedge, int nxalist)  ;
void setadmfix(char *fixname) ;
void setinbreed(int val) ;

