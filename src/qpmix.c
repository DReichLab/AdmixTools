#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include <nicklib.h>
#include <globals.h>
#include <mcmcpars.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "regsubs.h"  
#include "egsubs.h"  
#include "qpsubs.h" 
#include "eigsubs.h" 

/** 
 See /home/np29/biologyx/Admix5.0/src/qpmdir for examples on how to run 
 Can't find .tex file for theory  ; but pdf in ccpaper
*/

#define Y  0 
#define E  1
#define I1 2
#define I2 3
#define A  4
#define O  5

#define WVERSION   "950"
// was qpadmlin
// inbreed added
// mvname, weightname added
// Monte Carlo Integration to get std. errors on coeffs

#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

 double gslprecision = .0001;
 int gslmaxiter = 10000 ;
 int gsldetails = NO ;


char *parname = NULL ;
char *trashdir = "/var/tmp" ;

char *dumpname = NULL ;

//int verbose = NO ;
//int plotmode = NO ;
int debug = NO ; 
int vprint = NO ; 

int qtmode = NO ;

int sizeweight = NO ;
int inbreed = NO ;

Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int markerscore = NO ;
int seed = 0 ;
int chisqmode = NO ;  // approx p-value better to use F-stat
int missingmode = NO ;
int dotpopsmode = YES ;
int noxdata = YES ;  /* default as pop structure dubious if Males and females */
int pcorrmode = NO ;
int pcpopsonly = YES ;
int znval = -1 ;
int popsizelimit = -1 ;
// int fancynorm = YES ;
int gfromp = NO ;  // genetic distance from physical 

double plo = .001 ;
double phi = .999 ;
double pvhit = .001 ;
double pvjack = 1.0e-6 ;
double blgsize = 0.05 ;  // block size in Morgans */
double *chitot ;
int   *popsize ;
int  xmode = NO ;

int xchrom = -1 ;
int zchrom = -1 ;

char *genotypename = NULL ;
char  *snpname = NULL ;
char  *snpoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *outliername = NULL ;
char *outpop = NULL ; 
// list of outliers
// int outnum = -1 ;
double lambdascale = -1 ;  
int numboot = 1000 ; 
int maxboottry = -1 ; 
int mctrials = -1 ;


char *weightname = NULL ;
char **eglist ;
int  numeg = -1 ;
char *mvname = NULL; 
double fudgevar = -1 ; // null value

int deletefstats = NO ;
char *fstatslog = NULL ; 
char *fstatsname = NULL ; 
char *fstatsoutname = NULL ; 
char *basepop = NULL ; 
int keeptmp = NO ;

void readcommands(int argc, char **argv) ;

int calcmat(double *ccc, double **vccc, SNP **xsnplist, int *xindex, int *xtypes, 
  int nrows, int ncols, int numeg, int nblocks) ; 

void printmm(double *mmm, int dim) ;  // debug
int calccoeffs(double *coeffs, double *cvar, double *cmat, 
 double **vcmat, int dim, int nblocks)  ;
int writemv(char *mvname, long *kmean, long *kvar, int n) ; 
void dof2bug (double *f2score,  double *f2sig, SNP ** xsnplist, int *xindex,
	   int *xtypes, int nrows, int ncols, int numeg, int nblocks) ;
void calcx(double *zmean, double *zvar, double *cmat, double **vcmat, int dim, int nblocks)  ;
int mkcoeff(double *coeff, double *mat, int dim) ;
int estcoeffs(double *coeffs, double *cvar, double *zmean, double *zvar, int dim, int numboot, int maxtry) ;

void dumpz(double *z, double *zv, char *dumpname, int dim) ;  
double rank1opt(double *ff, double *z, double *zv, int dim)  ;
double setpars(double *yy, int dim)  ;
double chol2prod(double *full, double *yy, int dim)  ;
void solveco(double *coeff, double *mat, int dim) ;
void setr1mat(int **r1mat, int dim) ;
int setconstraint(double **x1mat, int **r1mat, double *wvec, int dim, int *pivarr)  ;
double lognormint(double *Q, double *L, double C, int n) ;
double phiint(double **x1mat, double *rhs, double *cmean, double *cvar, int dv, int numconstraint) ; 
double phival(double *ll, double *qq, double *pkon, int dim, double **constraint, double *rhs, int numc, int *pivot) ;
void mcest(double *qcoeffs, double *qvar,  
 double *coeffs, double *cvar, double *zmean, double *zvar, int dim, int mctrials) ;
double mcsamp(double *qcoeffs, double *qvar,  
 double *coeffs, double *cvar, double *zmean, double *zvar, int dim, int btrials, int printit) ;

void printcv(double *coeffs, double *var, int n) ;
double ess(double *wt, int n)  ;

int getgsldetails() ;
void setgsldetails(int val) ;


int dim2dv(int dim)  ;
int dv2dim(int dv)  ;
int mkfstats (char *parname) ;
void loadymvmix(double *ymean, double *yvar,  char *fstatsname, char **poplist,  int npops) ;

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  int  lsnplist, lindlist ;
  int i, j, k, k1, k2, k3, k4, kk, isbad; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double y1, y2, y, sig, tail, yy1, yy2 ;
  char ss[11] ;
  int *blstart, *blsize, nblocks ;
  int  xnblocks ; /* for xsnplist */
  int *bcols ;
  double maxgendis ;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n1, n2 ;

  int nindiv = 0, e, f, lag=1  ;
  double xc[9], xd[4], xc2[9] ;
  int nignore, numrisks = 1 ;
  double  *xrow, *xpt ; 
  SNP **xsnplist  ;
  int *tagnums  ;
  Indiv **xindlist ;
  int *xindex, *xtypes ;
  int nrows, ncols, m, nc ;
  int weightmode = NO ;
  double chisq, ynrows ;
  int *numhits, t ;  
  double *xmean, *xfancy ;
  double *divans, *divsd ; 
  double *hettop, *hetbot ; 
  int chrom,  numclear ;
  double gdis ;
  int outliter, *badlist, nbad ;
  double **zdata, *z1, *z2 ;
  int maxtag = -1 ;
  double **zz ; 
  double *ww  ;
  int *vvint ; 
  int  *rawcol ; ;
  int a, b, c, d, col, dv  ;
  double *qpscores ;
  double *hest, *hsig ;
  double mingenpos, maxgenpos ;
  int *qhit  ;  /* number of times pair is clade in quartet */
  int *qmiss ;  /* number of times pair migration event implied */
  int **qplist, numqp = 0, maxqp=10000 ;
  double *qpscore ;
  double jest, jsig, *wvar, *www ;
  double *coeffs, *cvar ; 
  long *jmean, *jvar ;
  double *zmean, *zvar ;
  double *rhs ;

  double ymin ;
  double *zlin ;
  double *f2, *f2sig, *fst ;

  int ng2, ng3, ng4 ;
  int cmsize ;
  double ymem ; 

  double *cmat, **vcmat ;
  int ret ; 
  double ychi, ytail, ysc, *ff ;
  double *wz, *wzv, *yy ;
  double *qcoeffs, *qvar ;
  int numconstraint ;
  double *wtsinit ; 
  double *evecs, *evals ;

  readcommands(argc, argv) ;
  printf("## qpmix version: %s\n", WVERSION) ;

  cputime(0) ;
  calcmem(0) ;

  if (parname == NULL) return 0 ;

  setinbreed(inbreed) ;
  setsizeweight(sizeweight) ;

  if (seed==0) seed = seednum() ; 
  printf("seed: %d\n", seed) ;
  SRAND(seed) ;

  if (fudgevar<0) { 
   fudgevar = .0001 ; 
  } 
  else { 
   printf("fudgevar: %12.6f\n", fudgevar) ;
  }

  if (xmode) xchrom = numchrom + 1 ;
  if (xchrom == (numchrom + 1)) noxdata = NO;
  if (basepop != NULL) printf("basepop set: %s\n", basepop) ;


   t = numlines(poplistname) ;
   ZALLOC(eglist, t+1, char *) ; 
   if (poplistname == NULL) fatalx("No poplistname!\n") ;
   numeg = loadlist(eglist, poplistname) ;

  d = numeg - 1 ; 
  cmsize = d*d+1 ;
// a little trickery last element is a count
  if (d <= 1) fatalx("bad d dimension\n") ; 

  dv = dim2dv(d) ;
  ZALLOC(zmean, d*d, double) ;
  ZALLOC(zvar, d*d*d*d, double) ; // should realy be dv (see calcx) 
  ZALLOC(wz, d*d, double) ;
  ZALLOC(wzv, d*d*d*d, double) ; // should realy be dv (see calcx) 
  ZALLOC(yy, d*d, double) ; // should really be dv 


 if (fstatsname == NULL) { 

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;

  ZALLOC(popsize, numeg, int) ;
  setpopsize(popsize, eglist, numeg, indivmarkers, numindivs) ;
  for (i=0; i<numeg; ++i) { 
   printf("pop: %20s  %5d\n", eglist[i], popsize[i]) ;  
  }
  ivlmaxmin(popsize, numeg, NULL, &t) ; 
  if (popsize[t] == 0) fatalx("pop %s has zero samples!\n", eglist[t]) ;
  if ((popsize[0] == 1) && (inbreed == YES)) fatalx("target pop %s has 1 sample and inbreed set!\n", eglist[0]) ;


  setgklist(indivmarkers, numindivs, eglist, numeg)   ;
  gkignore(indivmarkers, numindivs) ;

  k = getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    chrom = cupt->chrom;
    if ((xchrom > 0) && (chrom != xchrom))
      cupt->ignore = YES;
    if ((noxdata) && (chrom == (numchrom + 1)))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > (numchrom + 1))
      cupt->ignore = YES;

    if (chrom == zchrom)
      cupt->ignore = YES;
  }

  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;
  printf("numsnps: %d\n", numsnps) ;

  printf ("jackknife block size: %9.3f\n", blgsize);

 }

   if (numeg < 2) fatalx("Need at least 2 pops in %s\n", poplistname) ;

   if (numindivs>0) seteglist(indivmarkers, numindivs, poplistname);

 for (i=0; i<numeg; i++) {  
  printf("%3d %s\n",i, eglist[i]) ;
 }

 if (fstatsname == NULL) { 
  ZALLOC(xindex, numindivs, int) ;
  ZALLOC(xindlist, numindivs, Indiv *) ;
  ZALLOC(rawcol, numindivs, int) ;  
  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
  ZALLOC(xtypes, nrows, int) ;
  for (i=0; i<nrows; i++) {
   indx = xindlist[i] ;
   k = indxindex(eglist, numeg, indx -> egroup) ;
   xtypes[i] = k ;
// printf("zz %d %s :: %d\n", i, indx -> egroup, xtypes[i]) ; 
  }

  printf("before setwt numsnps: %d\n", numsnps) ;
  setwt(snpmarkers, numsnps, indivmarkers, nrows, xindex, xtypes, outpop, eglist, numeg) ;
  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;  //  rid ignorable snps
  printf("setwt numsnps: %d\n", numsnps) ;

  setindm(indivmarkers) ;

  if (weightname != NULL) {  
    printf("loading weights\n") ;
    getweights(weightname, snpmarkers, numsnps) ;
  }


  if (numsnps ==  0) fatalx("no valid snps\n") ;

  setmgpos(snpmarkers, numsnps, &maxgendis) ;  
  if ((maxgendis < .00001) || (gfromp == YES)) setgfromp(snpmarkers, numsnps) ;

  nblocks = numblocks(snpmarkers, numsnps, blgsize) ;
  ZALLOC(blstart, nblocks, int) ;
  ZALLOC(blsize, nblocks, int) ;
  printf("number of blocks for moving block jackknife: %d\n", nblocks) ;

  ZALLOC(xsnplist, numsnps, SNP *) ;
  ZALLOC(tagnums, numsnps, int) ;

  if (popsizelimit > 0) {  
   setplimit(indivmarkers, numindivs, eglist, numeg, popsizelimit) ; 
  }

  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;


  printf("snps: %d  indivs: %d\n", ncols, nrows) ;
  setblocks(blstart, blsize, &xnblocks, xsnplist, ncols, blgsize)  ;
  printf("number of blocks for moving block jackknife: %d\n", xnblocks) ;

  d = numeg - 1 ; 
  cmsize = d*d+1 ;
  ZALLOC(cmat, cmsize, double) ;
  vcmat = initarray_2Ddouble(nblocks, cmsize, 0) ; 
// a little trickery last element is a count

  dv = dim2dv(d) ;


  t = calcmat(cmat, vcmat, xsnplist, xindex, xtypes, nrows, ncols, numeg, nblocks) ; 
// estimate f-stat matrix vcmat is on single block
  if (t==0) fatalx("no f3 stats could be calculated!\n") ;
// big step 1
  calcx(zmean, zvar, cmat, vcmat, d, nblocks) ;
// estimate mean f-stat matrix (upper triangle) and  and covariance by jackknife

 }

 if (fstatsname != NULL) { 
  loadymvmix(zmean, zvar,  fstatsname, eglist,  numeg) ;
 }


  if (dumpname != NULL) dumpz(zmean, zvar, dumpname, d) ;
  if (maxboottry < 0) maxboottry = 10*numboot ;

   ZALLOC(ff, d*d, double) ; 
   vst(wz, zmean, 1000, dv) ; 
   vst(wzv, zvar, 1000*1000, dv*dv) ; 
   if (isposdef(wzv, dv) == NO) fatalx("zvar bug\n") ;
   mkfull(ff, wz, d) ; 
   printf("saturated matrix\n") ;
   printmatwl(ff, d, d, d) ;
   printnl() ;
// calculate max likelihood corank1 matrix 
   ychi = rank1opt(yy, wz, wzv, d) ;           
   chol2prod(ff, yy, d) ; 
   setpars(yy, d) ;
   ytail = rtlchsq(1, ychi) ; 
   printf("best corank1 matrix: chisq: %9.3f\n", ychi) ; 
   printmatwl(ff, d, d, d) ;
   printnl() ; 
   printf("zzyy:") ; printmatwl(yy, 1, 6, 6) ;  
   ZALLOC(evals, d, double) ;
   ZALLOC(evecs, d*d, double) ;

   eigvecs(ff, evals,  evecs, d) ;  
   printf("evals: ") ;
   printmatl(evals, 1, d) ;
   printnl() ;
   printmatl(evecs, d, d) ;
   printnl() ;


   ZALLOC(wtsinit, d, double) ; 
   xline(wtsinit, ff, d, d) ;
   printf("wtsinit: ") ; 
   printmat(wtsinit, 1, d) ;

  ZALLOC(coeffs, d, double) ;
  ZALLOC(cvar, d*d, double) ;

  estcoeffs(coeffs, cvar, wz, wzv, d, numboot, maxboottry) ;
/** 
 this gets a rough estimate of covariance of coeffs.  Not critical
 Big step 2
*/

/**
  copyarr(wtsinit, coeffs, d) ; 
  setidmat(cvar, d) ; 
  y = 0.05 * 0.05 ; // variance
  vst(cvar, cvar, y, d*d) ; // hack
*/

// here we bootstrap on coeff matrix and compute coeff mean and covariance.  
// Insist matrix is pos def                                

// now we sample coeffs using posterior 
  
  t = d + dv ;
  ZALLOC(ww, t*t, double) ;

  
  getdiag (ww, cvar, d);
  vsp (ww, ww, 1.0e-20, d);
  vsqrt (ww, ww, d);
  printnl() ;
  printnl() ;

  printf("old method:\n")  ;
  printcv(coeffs, cvar, d) ;

  getdiag(ww, cvar, d); 
  vsqrt(ww, ww, d) ;

//  printf("       coeffs: ") ; printmat(coeffs, 1, d) ;
//  printf("  std. errors: ");  printmat(ww, 1, d) ;
  printnl() ;

  printf("Target: %25s (old)\n", eglist[0]) ;
  for (i=1; i<numeg; ++i) { 
   printf("%25s ", eglist[i]) ;
   j = i-1; 
   printf("%9.3f  ", coeffs[j]) ;
   printf("%12.6f", ww[j]) ;
   printf(" (old)") ; 
   printnl() ;
  }

  printnl() ;
  ZALLOC(qcoeffs, d, double) ;
  ZALLOC(qvar, d*d, double) ;
  
  copyarr(coeffs, qcoeffs, d) ;
  copyarr(cvar, qvar, d*d) ;


  printnl() ;  
  if (mctrials < 0) { 
   mctrials = MIN(pow(2, d), 10000) ; 
   mctrials = MAX(mctrials, 1000) ;
  }

  
  printf("calling mcest; mctrials: %d\n", mctrials) ;  
  mcest(qcoeffs, qvar,  coeffs, cvar, zmean, zvar, d, mctrials) ;
// Big step 3 ;  Importance sample coeffs using iterative estimate of distribution


   printf("mcest: (coeff * 1000)\n") ;
   printcv(qcoeffs, qvar, d) ;

  getdiag(ww, qvar, d); 
  vsqrt(ww, ww, d) ;

  


  printf("Target: %s\n", eglist[0]) ;
  for (i=1; i<numeg; ++i) { 
   printf("%25s ", eglist[i]) ;
   j = i-1; 
   printf("%9.3f  ", qcoeffs[j]) ;
   printf("%12.6f", ww[j]) ; 
   printnl() ;
  }

  printnl() ;
  fflush(stdout) ; 

  ZALLOC(jmean, d, long) ;
  ZALLOC(jvar, d*d, long) ;

  vst(ww, qcoeffs, 1000, d) ; 
  fixitl(jmean, ww, d) ; 
  vst(ww, qvar, 1.0e6, d*d) ;   
  fixitl(jvar, ww, d*d) ;  
  writemv(mvname, jmean, jvar, d) ;

  printnl() ;
  printnl() ;

  printf("chisq: %12.6f  p-value: %12.6f\n", ychi, ytail) ;
  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qpmix: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

}

int calcmat(double *ccc, double **vccc, SNP **xsnplist, int *xindex, int *xtypes, 
	int nrows, int ncols, int numeg, int nblocks) 
{


  int t1, t2;
  int a, b, c, t;
  int isok, ret ;  
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  double ytop, ybot, wt, yy ;
  SNP *cupt;
  double *w1, *w2, *w3, *ww;
  int bnum;
  int dim = numeg - 1, cmsize ;
  static int ncall = 0 ; 
  int fsx[4] ;
  int numok = 0 ;

  ++ncall ;
  cmsize = dim*dim + 1 ; 
  ZALLOC (ww, cmsize, double);

//  dof2bug (&ya,  &yb,  xsnplist, xindex, xtypes,  nrows,  ncols,  numeg,  nblocks) ;

  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    verbose = NO ; 
    t = ranmod(1000) ; 
    if (t==-1) verbose = YES ;
    if (verbose) printf("zzrandom: %d\n", col) ;
    wt = cupt -> weight ;
    if (wt <= 0.0) continue ;
    if (cupt->ignore) continue;
    bnum = cupt->tagnumber;
    if (bnum < 0)
      continue;

    vzero(ww, dim*dim) ; 
    a = 0 ; 

   loadaa (cupt, xindex, xtypes, nrows, numeg);

   isok = YES ; 
    for (b=1; b < numeg; ++b) {  
     for (c=b; c < numeg; ++c) {  
      if (isok == NO) break ;
       fsx[0] = fsx[2] = 0 ; 
       fsx[1] = b ; 
       fsx[3] = c ;
       ytop = fstatx(fsx) ;
       ret = 1 ; 
       if (ytop < -99) ret = -9 ; 
       if (ytop > 100) verbose = YES ;  
 
       if (verbose) printf("zzfstatx: %d %d  %12.6f\n", b, c, ytop) ;
       if (ret<0) isok = NO ;
       k = (b-1)*dim + c - 1 ; ww[k] = ytop ;
       k = (c-1)*dim + b - 1 ; ww[k] = ytop ;
     }
    }
    if (verbose == YES) printf("isok: %d\n", isok) ;
    if (isok == NO) continue ;
    ++numok ;
    ww [dim*dim] = 1 ; 
    vst(ww, ww, wt, cmsize) ;

    vvp (vccc[bnum], vccc[bnum], ww, cmsize) ;
    vvp (ccc, ccc, ww, cmsize) ;
    if (verbose) { 
     printf("zzwr; %6.0f %6.0f\n", ww[dim*dim], wt) ; 
     printf("count: %6.0f\n", ccc[dim*dim]) ; 
     printmat(ccc, dim, dim) ; 
    }
  }

  free(ww) ;
  return numok ; 

} 
void calcx(double *cmean, double *cvar, double *cmat, double **vcmat, int dim, int nblocks) 
{


  double **vvv, *vv, y ; 
  double **vvt, *vt, *wjack ; 
  double *vmean, *vcovar ;
  int d2, k, dv ; 

  d2 = dim*dim ;
  dv = dim*(dim+1)/2 ;

  if (cmat[d2] <= 0.1) fatalx("(qqpmixpmix) no data!\n") ;

  vvv = initarray_2Ddouble(nblocks, d2+1, 0) ; 
  vvt = initarray_2Ddouble(nblocks, dv, 0) ; 
  ZALLOC(wjack, nblocks, double) ;
  ZALLOC(vv, d2, double) ;
  ZALLOC(vt, dv, double) ;
  y = cmat[d2] ; 
  vst(vv, cmat, 1.0/y, d2) ; 
  mktriang(vt, vv, dim) ;
  for (k=0; k<nblocks; ++k) { 
   vvm(vvv[k], cmat, vcmat[k], d2+1) ; 
   wjack[k]  = y = vvv[k][d2] ; 
   vst(vvv[k], vvv[k], 1.0/y, d2) ; 
   mktriang(vvt[k], vvv[k], dim) ;
  }

  wjackvest (cmean, cvar, dv, vt, vvt, wjack, nblocks);

/**
  printf("zzcalcx\n") ;
  printmat(vv, 2, 2) ;
  printmat(vt, 1, dv) ;
  printmat(cmean, 1, dv) ;
    printnl() ;
  printmatl(cvar, dv, dv) ;
  printnl() ;
*/
  


  free(wjack) ;
  free(vv) ;
  free(vt) ;
  free2D(&vvv, nblocks) ;
  free2D(&vvt, nblocks) ;

}

void matnorm(double *xmat, double *mat, int n)  
// xmat can equal mat
{
   double ymul, ytrace ;

   ytrace = trace(mat, n) ;
   if (ytrace == 0.0) fatalx("(matnorm) zero trace!\n") ;
   ymul = (double) n / ytrace ;
   vst(xmat, mat, ymul, n * n ) ;

}

int mkcoeff(double *coeff, double *mat, int dim) 
{
 double *ll, y ; 

 if (isposdef(mat, dim) == NO) return -1 ;
 
 ZALLOC(ll, dim, double) ;
 vclear(ll, 1.0, dim) ;

 solvit(mat, ll, dim, coeff) ;  
 y = asum(coeff, dim) ; 
 vst(coeff, coeff, 1.0/y, dim) ;  // coeff can be  negative.  Don't use bal1

 return 1 ;

 free (ll) ;

}

int
calccoeffs(double *coeffs, double *cvar, double *cmat, double **vcmat, int dim, int nblocks) 
{
   double *mean, **jmean, *jwt, *tans, *tmat, *ll ; 
   double wt, y ; 
   int k, ret, len ;

   ZALLOC(mean, dim, double) ;
   ZALLOC(tans, dim, double) ;
   ZALLOC(ll, dim, double) ;
   ZALLOC(tmat, dim*dim, double) ;
   ZALLOC(jwt, nblocks, double) ;

   jmean = initarray_2Ddouble(nblocks, dim, 0) ; 

   vzero(coeffs, dim) ; 
   setidmat(cvar, dim) ;
   vclear(ll, 1.0, dim) ;
   matnorm(tmat, cmat, dim) ; // not needed
   ret = solvit(tmat, ll, dim, tans) ;

   if (ret<0) { 
    printf("*** warning: coefficient matrix not positive definite ***\n") ;
    printmat(tmat, dim, dim) ; 
    return ret ;
   }
   bal1(tans, dim) ; 
   copyarr(tans, coeffs, dim) ;
   copyarr(tans, mean, dim) ;

   len = dim*dim ; 

   for (k=0;  k<nblocks; ++k) { 
    wt = vcmat[k][len] ;
    vvm(tmat, cmat, vcmat[k], dim*dim) ;
    matnorm(tmat, tmat, dim) ;
    ret = solvit(tmat, ll, dim, tans) ;

    if (ret < 0) { 
     copyarr(mean, tans, dim) ;
     wt = 0.0 ;  
     printf("*** block %d not pos. def.\n", k ) ; 
     printmat(tmat, dim, dim) ;
    }
    bal1(tans, dim) ; 
    copyarr(tans, jmean[k], dim) ;
    jwt[k] = wt ;
   }
   wjackvest(coeffs, cvar, dim, mean, jmean, jwt, nblocks) ;
/**
   printf("zzans:\n") ; 
  
   printmat(mean, 1, dim) ;
   printmat(coeffs, 1, dim) ;
   printmat(cvar, dim, dim) ;
*/
    
   free(mean) ;
   free(tans) ;
   free(tmat) ;
   free(jwt) ;

   free2D(&jmean, nblocks) ;

   return 1 ; 

}
void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   if (parname==NULL) { 
    fprintf(stderr, "no parameters\n") ;
	    return ;
   }

   pcheck(parname,'p') ;
   printf ("%s: parameter file: %s\n", basename(argv[0]), parname);
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "outpop:", &outpop) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getdbl(ph, "blgsize:", &blgsize) ;
   getdbl(ph, "lambdascale:", &lambdascale) ;
   getint(ph, "xmode:", &xmode) ;
   getint (ph, "seed:", &seed) ; 
   getint(ph, "noxdata:", &noxdata) ; 
   getint(ph, "inbreed:", &inbreed) ; 
   getint (ph, "numboot:", &numboot) ; 
   getint (ph, "maxboottry:", &maxboottry) ; 
   getint (ph, "mctrials:", &mctrials) ; 
   getint (ph, "gsldetails:", &gsldetails) ; 
   getint (ph, "gslmaxiter:", &gslmaxiter) ; 
   getdbl (ph, "gslprecision:", &gslprecision) ; 
   getdbl (ph, "fudgevar:", &fudgevar) ; 
   getdbl (ph, "fudgev:", &fudgevar) ; 
   getstring (ph, "basepop:", &basepop);

   getstring (ph, "fstatsname:", &fstatsname);
   getstring (ph, "fstatsoutname:", &fstatsoutname);
   getstring (ph, "fstatslog:", &fstatslog);
   getint (ph, "keeptmp:", &keeptmp);

   getint(ph, "popsizelimit:", &popsizelimit) ; 
   getint(ph, "gfromp:", &gfromp) ;  // gen dis from phys
   getint(ph, "chrom:", &xchrom) ;  
   getstring(ph, "weightname:", &weightname) ;
   getstring(ph, "mvname:", &mvname) ;
   getstring(ph, "dumpname:", &dumpname) ;
   getint(ph, "vprint:", &vprint) ;

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}

void printmm(double *mmm, int dim) 
// debug routine
{
  int len = dim*dim ; 
  double *ww, *ll, *cc  ; 
  double ysize, y, tr ;     
  int ret ; 

  ZALLOC(ww, len, double) ;
  ZALLOC(ll, dim, double) ;
  ZALLOC(cc, dim, double) ;

  ysize = mmm[len] ; 
  printf("mmm: %d %9.3f\n", dim, ysize) ; 
  if (ysize < 1.0e-6) fatalx("nothing here\n") ;
  vst(ww, mmm, 1.0/ysize, len) ; 
  y = trace(ww, dim) ; 
  y = (double) dim / y ; 
  vst(ww, ww, y, len) ; 

  printmat(ww, dim, dim) ; 
  vclear(ll, dim, 1.0) ;
  ret = solvit(ww, ll, dim, cc) ;
  if (ret < 0) printf("cmat not pos. def.\n") ;


  y = asum(cc, dim) ; 
  vst(cc, cc, 1.0/y, dim) ; 
  printf("coeffs: ") ; printmat(cc, 1, dim) ;


 free(ww) ;
 free(ll) ;
 free(cc) ;

}
int
writemv(char *mvname, long *kmean, long *kvar, int n)
{
 FILE *fff ;
 int trap ; 

 if (mvname == NULL) return NO ; 
 trap = openit_trap(mvname, &fff, "w") ;

 if (trap == NO) return trap ;

 printlmatwfile(kmean, 1, n, n, fff) ;
 printlmatwfile(kvar, n, n, n, fff) ;
 fclose(fff) ;
 
 return trap ;
}

void
dof2bug (double *f2score,  double *f2scoresig, SNP ** xsnplist, int *xindex,
	   int *xtypes, int nrows, int ncols, int numeg, int nblocks)
{

  int t1, t2;
  int a, b, c, d;
  int c1[2], c2[2], *cc;
  int *rawcol, *popall, *pop0, *pop1;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  SNP *cupt;
  double top, bot, *djack, *wjack, gtop, gbot, *wbot, *wtop;
  double *btop, *bbot, wt;
  double ytop, ybot;
  double y1, y2, yscal;
  double *w1, *w2, *ww, m1, m2;
  int bnum, totnum;

  if (nrows == 0)
    fatalx ("badbug\n");

  ZALLOC (w1, nblocks, double);
  ZALLOC (w2, nblocks, double);
  ZALLOC (ww, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  ZALLOC (djack, nblocks, double);
  ZALLOC (btop, nblocks, double);
  ZALLOC (bbot, nblocks, double);
  ZALLOC (wtop, nblocks, double);
  ZALLOC (wbot, nblocks, double);


  totnum = 0;
  gtop = gbot = 0 ; 
  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    if (cupt->ignore)
      continue;

    bnum = cupt->tagnumber;
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");

    f2sc (&ytop, &ybot, cupt, indivmarkers, xindex, xtypes, nrows, 2, 0, 1);
    if (isnan (ytop))
      fatalx ("zznan\n");
    if (ybot < -0.5)
      continue;

    ybot = 1.0 ;

    btop[bnum] += ytop;
    bbot[bnum] += ybot;

    gtop += ytop ;
    gbot += ybot ;
    ++wjack[bnum];
    ++totnum;

  }

//  printf("zzzf2 %12.6f %12.6f %12.6f\n", gtop, gbot, gtop/gbot) ;

  gtop = asum (btop, nblocks);
  gbot = asum (bbot, nblocks);

  *f2score = mean = gtop / gbot;

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    wtop[k] = gtop - top;
    wbot[k] = gbot - bot;
    wbot[k] += 1.0e-10;
    djack[k] = wtop[k] / wbot[k];	// delete-block estimate
  }

  wjackest (&jest, &jsig, mean, djack, wjack, nblocks);

  *f2scoresig = jsig;

  free (w1);
  free (w2);
  free (ww);

  free (wtop);
  free (wbot);
  free (djack);
  free (wjack);

  free (btop);
  free (bbot);

}

int 
estcoeffs(double *coeffs, double *cvar, double *zmean, double *zvar, int dim, int numboot, int maxtry) 
// bootstrap over observed f-stats and cmpute max likelihood rank 1 and ceeffs 
{

 int dv ; 
 int nboot=0, ntry = 0, t, k ; 
 double *rvec, *fmat, *yy ;
 double **bootcoeffs, *bco ;
 double y ; 
 int gsld, gslt ;  
 
 printf("zzgsl: %d\n", gslmaxiter) ;

 gsld = getgsldetails() ;  // controls output; not important
 dv = dim*(dim+1)/2 ;

 ZALLOC(rvec, dv, double) ;
 ZALLOC(yy, dv, double) ;
 ZALLOC(fmat, dim*dim, double) ;
 bootcoeffs = initarray_2Ddouble(numboot, dim, 0) ;

 for (k=0; k<numboot; ++k) {  
  genmultgauss(rvec, 1, dv, zvar) ; 
  vvp(rvec, rvec, zmean, dv) ;
  gslt = gsld ; 
   
//  if (k<2) gslt = YES ;
  setgsldetails(gslt) ; 

  rank1opt(yy, rvec, zvar, dim) ;           
  chol2prod(fmat, yy, dim) ;  // corank 1
  xline(bootcoeffs[k], fmat, dim, dim) ;
 }
   vzero(coeffs, dim) ; 
   vzero(cvar, dim*dim) ; 

   for (k=0; k<numboot; ++k) { 
    vvp(coeffs, coeffs, bootcoeffs[k], dim) ; 
   }
   y = 1.0/(double) numboot ; 
   vst(coeffs, coeffs, y, dim) ;  // mean 

   for (k=0; k<numboot; ++k) { 
    bco = bootcoeffs[k] ;
    vvm(bco, bco, coeffs, dim) ; 
// subtract mean
    addouter(cvar, bco, dim) ;
   } 
   vst(cvar, cvar, y, dim*dim) ;


   free(rvec) ;
   free(yy) ;
   free(fmat) ;

   free2D(&bootcoeffs, numboot) ;

   return numboot ;
}

void dumpz(double *z, double *zv, char *dumpname, int dim) 
{

  FILE *fff ; 
  double *ww ; 
  int dv ; 
  
 if (dumpname == NULL) return ;  
 dv = dim*(dim+1)/2 ;

 ZALLOC(ww, dv*dv+10, double) ; 
 openit(dumpname, &fff, "w") ;
 fprintf(fff, "##dumpz: %d\n", dim) ;

 vst(ww, z, 1000, dv) ;
 printmatwfile(ww, 1, dv, dv, fff) ;

 vst(ww, zv, 1000*1000, dv*dv) ;
 printmatwfile(ww, dv, dv, dv, fff) ;

 free(ww) ; 
 fclose(fff) ; 

}
void q2v(double *lin, double *q, double *m, double *v, int dim) 
{
  int t ; 
 
  t = isposdef(q, dim) ; 
  if (t==NO) fatalx("(q2c) v not posdef\n") ;
  pdinv(v, q, dim) ; 
  mulmat(m, v, lin, dim, dim, 1) ;
}

void v2q(double *lin, double *q, double *m, double *v, int dim) 
{
  int t ; 
 
  t = isposdef(v, dim) ; 
  if (t==NO) fatalx("(v2q) v not posdef\n") ;
  pdinv(q, v, dim) ; 
  mulmat(lin, q, m, dim, dim, 1) ;
}

void
weightsamp(double **sampco, double *logweight, double *cmean, double *cvar, double *zmean, double *zvar, int numsamp, int dim) 
{
/** we generate 
numsamp samples with weights logweight[k] 
1) generate v mvg(cmman, cvar) log likelihood L(v) 
2) W(v) = log int (X; zmean, zvar) satisfying constraints  
3 weught = W(v)-L(v) 
*/  

 return  ;
 

}
double lognormint(double *Q, double *L, double C, int n) 
// log int exp -2 \sum_Q_{ij} x_i x_j -2 \sum L_i x_i + C 
{

  double y, yc, yl, yn, yldet ;
  double *mean ; 
  static double logtwopi = -1  ; 

  if (logtwopi < 0) { 
   y = pi() ; 
//  printf("pi: %15.10f", y) ;
   logtwopi = log(2.0*y) ;
  }

  ZALLOC(mean, n, double) ;
  solvit(Q, L, n, mean) ; 
  yc = C - vdot(mean, L, n) ; 
  yl = -0.5*yc ; 
  yn = (double) n ; 

  yl += (yn/2)*logtwopi ;  
  yldet = pdinv(NULL, Q, n) ; 
  yl -= 0.5 * yldet ; 
  
  free(mean) ;
  return yl ; 

}

int dv2dim(int dv) 
{
  int t ; 
 
  t = nnint(sqrt(8*dv+1))  ; 
  return (t-1)/2 ;

} 

int dim2dv(int dim) 
{

 return dim*(dim+1)/2 ;

}
void setzf(double *zfmean, double *zfvar, double *zmean, double *zvar, int dim) 
{
  int *tab, d2 = dim*dim ; 
  int a, b, c, d, x, u, v, dv, i, j ;
  double y ; 

  mkfull(zfmean, zmean, dim) ; 

  ZALLOC(tab, d2, int) ;

  x = 0 ;  
  for (a=0; a<dim; ++a) {
   for (b=a; b<dim; ++b) { 
     tab[a*dim+b] = tab[b*dim+a] = x ; 
     ++x ; 
  }}
  vclear(zfvar, -99, d2) ;
  dv = dim2dv(dim) ;
  for (a=0; a<dim; ++a) {
   for (b=0; b<dim; ++b) { 
    i = a*dim+b ;
    u = tab[i] ; 
    for (c=0; c<dim; ++c) {
     for (d=0; d<dim; ++d) { 
      j = c*dim+d ;
      v = tab[j] ; 
      y = zvar[u*dv+v] ; 
      zfvar[i*d2 + j] = y ;
   }}}}
   vst(zfmean, zfmean, 1000, d2) ;
   vst(zfvar, zfvar, 1000*1000, d2*d2) ;
   printf("setzf:\n") ;  
   printmatw(zfmean, dim, dim, dim) ;
   printnl() ;
   printmatw(zfvar, d2, d2, d2) ;
   printnl() ;

 free(tab) ;

}
void mcest(double *qmean, double *qvar,  
 double *coeffs, double *cvar, double *zmean, double *zvar, 
 int dim, int mctrials) 
{
    double *mz, *mzvar ; 
    double y, ylike, yy ;
    int dv, numc, t, k ;
    double yeff ;
    double *zfmean, *zfvar ; 
    int d2 = dim*dim ; 

    double *mcmean, *mcvar, *xmean, *xvar ;
    int ntrials = 10000 ;

    dv = dim2dv(dim) ;

    ZALLOC(mcmean, dim, double) ;
    ZALLOC(mcvar, dim*dim, double) ;

    ZALLOC(xmean, dim, double) ;
    ZALLOC(xvar, dim*dim, double) ;

    ZALLOC(zfmean, d2, double) ;
    ZALLOC(zfvar, d2*d2, double) ;

    setzf(zfmean, zfvar, zmean, zvar, dim) ;

    copyarr(coeffs, xmean, dim) ;
    copyarr(cvar, xvar, dim*dim) ;
/** these are oldstyle estimates */

    for (;;) { 
     for (k=0; k<dim; ++k) { 
      y = xvar[k*dim+k] ; 
      y = MAX(y, fudgevar) ; 
      xvar[k*dim+k] = y ;  // fudge to avoid excessive concentration of prior 
     }
     t = NO ;
     if (ntrials > mctrials) t = vprint ; 
//   if (ntrials == 1000) t = YES ;
     yeff = mcsamp(mcmean, mcvar,  xmean, xvar, zfmean, zfvar, dim, ntrials, t) ;
     printf("big iter: trials: %d ; efficiency %12.6f\n", ntrials, yeff) ;
     printcv(mcmean, mcvar, dim) ;
     copyarr(mcmean, xmean, dim) ; 
     copyarr(mcvar, xvar, dim*dim) ;
     if (ntrials > mctrials) break ; 
     ntrials *= 2 ; 
    }
    copyarr(mcmean, qmean, dim) ; 
    copyarr(mcvar, qvar, dim*dim) ;

    free(mcmean) ;
    free(mcvar) ;
    free(xmean) ;
    free(xvar) ;
    free(zfmean) ;
    free(zfvar) ;

    return ; 
}

void printcv(double *coeffs, double *var, int n) 
{
   double *w1, *w2 ;

   ZALLOC(w1,  n, double) ;
   ZALLOC(w2,  n*n, double) ;

   vst(w1, coeffs, 1000, n) ;
   vst(w2, var, 1000*1000, n*n) ;
   printmatw(w1, 1, n, n) ;
   printnl() ;
   printmatw(w2, n, n, n) ;
   printnl() ;

   free(w1) ;
   free(w2) ;

}

int mkfstats(char *parname)  
// modified from qpAdm
{
 char sss[512] ; 
 char fsx[256] ; 
 char ppp[256] ;  
 char pops[256] ;  
 char tpar[256] ; 
 char fslog[256] ; 
 char **poplist, **popx ; 
 int numpop ; 

 int pid ; 
 int k ; 

 int retkode, trap ;

 FILE *fff ; 

 ZALLOC(poplist, numeg+1, char *) ; 

 popx = poplist ; 
 numpop = numeg ;

 if (basepop != NULL) {  
   poplist[0] = strdup(basepop) ;
   ++popx ; 
   ++numpop ;
 }
 
 copystrings(eglist, popx, numeg) ;

 pid = getpid() ; 
 strcpy(fsx, mytemp("fsx")) ; 
 strcpy(ppp, mytemp("ppp")) ; 
 strcpy(pops, mytemp("pops")) ; 
 strcpy(tpar, mytemp("tpar")) ; 
 strcpy(fslog, mytemp("fslog")) ; 
 if (fstatslog != NULL) strcpy(fslog, fstatslog) ;

/**
 strcpy(pops, "tjunk") ; 
 strcpy(ppp, "tppp") ; 
*/

 printf("tmp: %s\n", fsx) ;
 fstatsname = strdup(fsx) ; 
 deletefstats = YES ;  

 openitntry(ppp, &fff, "w", 10) ; 
 sprintf(sss, "cat %s", ppp) ;
// system(sss) ;
 sprintf(sss, "fgrep -v fstatsoutname: %s | fgrep -v poplistname: > %s", parname, tpar) ; 
 system(sss) ; 
// printf("zz2 %s\n", sss) ; 
 fflush(stdout) ; 
 copyfs(tpar, fff) ; 
 fprintf(fff, "fstatsoutname: %s\n", fsx) ; 
 writestrings(pops, poplist, numeg) ;
 
 
// printf("cmd1: %s\n", sss) ; 
 fprintf(fff, "poplistname: %s\n", pops) ; 
 fclose(fff) ; 

 sprintf(sss, "qpfstats -p %s > %s", ppp, fslog) ; 
 retkode =  system(sss) ; 
 trap = openit_trap(fsx, &fff, "r") ;
 fclose(fff) ; 
 if (trap == NO) retkode = -88 ;
//  printf("zz: %d\n", trap) ; 

 
//  printf("exiting mkfstats\n") ; 
 fflush(stdout) ;
 if ((fstatsoutname != NULL) && (retkode >= 0)) { 
  sprintf(sss, "cp %s %s", fsx, fstatsoutname) ; 
  system(sss) ;
  printf("fstats file: %s made\n", fstatsoutname) ;
 } 
 if (retkode < 0) return -1 ;
 freeup(poplist, numeg + 1) ;

 if (keeptmp == YES) return 1  ; 

 remove(ppp) ;
 remove(pops) ;
 remove(tpar) ; 

 if (fstatslog == NULL) remove(fslog) ;
 return 1 ; 

}

void load4(int *x, int a, int b, int c, int d)
{
 x[0] = a ;
 x[1] = b ;
 x[2] = c ;
 x[3] = d ;
}



void loadymvmix(double *ymean, double *yvar,  char *fstatsname, char **poplist,  int npops) 
{

  char **eglist ; 
  int numeg, np, numfs ; 
  double *ff3, *ff3var, *vest, *vvar ; 
  int a, b, c, d, i, j, t, tt, k, dim ; 
  int nl, nr ; 
  int **fsindex, *fs ; 
  double y, ysig, yv ; 

  

  ZALLOC(eglist, MAXPOPS, char *) ;
  np = fstats2popl(fstatsname, eglist) ; 

  printf("fstats pops: %s\n", fstatsname) ;
  printstrings(eglist, np) ;
  printnl() ;

  numeg = npops ;

  t = np*np ;
  ZALLOC(ff3, t, double) ; 
  ZALLOC(ff3var, t*t, double) ; 
  loadfstats(fstatsname, ff3, ff3var, eglist, np) ; 


  ZALLOC(vest, t, double) ; 
  ZALLOC(vvar, t*t, double) ;

  setvv(vest, vvar, ff3, ff3var, NULL, np) ; 

  t = npops * npops ;  numfs = 0 ; 
 
  fsindex = initarray_2Dint(t, 4, -1) ; 
  a = c = indxstring(eglist, np, poplist[0]) ; 
  for (i=1; i<npops; ++i) { 
   b = indxstring(eglist, np, poplist[i]) ; 
   if (b<0) {
     fatalx("%s not found\n", poplist[i]) ;
   }
   for (j=i; j<npops; ++j) { 
    d = indxstring(eglist, np, poplist[j]) ; 
   if (d<0) {
     fatalx("%s not found\n", poplist[j]) ;
   }
    fs = fsindex[numfs] ; 
    load4(fs, a, b, c, d) ; 
    ivmaxmin(fs, 4, NULL, &tt) ; 
    if (tt<0) fatalx("pop not found in fstats: %s %s %s %s\n", 
       poplist[0], poplist[i], poplist[0], poplist[j]) ; 
    ++numfs ; 

  }} 

  printf("fsindex: %d\n", numfs) ; 
  printimat2D(fsindex, numfs, 4) ;
  printnl() ;
  fflush(stdout) ;
  if (numfs==0) fatalx("numfs bug\n") ;

  vv2ww(ymean, yvar, vest, vvar, np, fsindex, numfs) ;

  printf("ymean stats %12.6f %12.6f\n", asum(ymean, numfs), asum2(ymean, numfs) ) ;
  printf("zza %15.9f %15.9f\n", vest[0], vvar[0]) ; 
  printf("zzb %15.9f %15.9f\n", vest[1], vvar[1]) ; 
  printf("zzc %15.9f %15.9f\n", ymean[0], yvar[0]) ; 


  free(eglist) ; 
  free(ff3) ; 
  free(ff3var) ; 
  free(vest) ; 
  free(vvar) ; 

  free2Dint(&fsindex, t) ; 

}

void setzmat(double *zmat, int n) 
{

  int a, b, x ;

  vzero(zmat, n*n*n) ;  

  for (a=0; a<n; ++a) { 
   for (b=0; b<n; ++b) { 
    x = a*n+b ; 
    zmat[a*n*n+x] = 1 ;
   }
  }
}

void mulmat3(double *outmat, double *mata, double *matb, double *matc, 
 int na, int nb, int nc, int nd)  
{

  double *w1 ;  

  ZALLOC(w1, na*nb, double) ;
 
  mulmat(w1, mata, matb, na, na, nb) ; 
  mulmat(outmat, w1, matc, na, nc, nd) ; 

  free(w1) ;

}
double phi2(double *wts, double *mean, double *var, int n) 
{
  double *wmat, *qmat ;
  double *wvec, y, yldet ;  
  int d2 = n*n ;
  int a, b, i, j, s, t ; 
  static long ncall = 0 ;

  ++ncall ;
   
  ZALLOC(wvec, n, double) ;
  ZALLOC(wmat, d2, double) ;
  ZALLOC(qmat, d2, double) ;

  mulmat(wvec, wts, mean, 1, n, n) ;

  for (a=0; a<n; ++a) { 
   for (b=0; b<n; ++b) { 

  for (i=0; i<n; ++i) { 
   for (j=0; j<n; ++j) { 

   s = a*n + i ; 
   t = b*n + j ;

   y = var[s*d2+t] ; 
   wmat[a*n+b] += wts[i]* wts[j] * y ; 

 }}}} 

  if (ncall == 1) {
   printf("zzphi2\n") ;
   printmat(wts, 1, n) ;
   printmat(wvec, 1, n) ;
   printnl() ;
   printmatw(mean, n, n, n) ;
   printnl() ;
   printmatw(var, n*n, n*n, n*n) ;
   printnl() ;
   printmatw(wmat, n, n, n) ;
  }


   yldet = pdinv(qmat, wmat, n) ;
   y = scx(qmat, NULL, wvec, n) ;

  free(wmat) ;
  free(qmat) ;
  free(wvec) ;

   return -0.5*(y + yldet) ; 
}

double mcsamp (double *qmean, double *qvar,  
 double *coeffs, double *cvar, double *zfmean, double *zfvar, 
 int dim, int ntrials, int printit) 
{
    double *mz, *mzvar, *mzinv ; 
    double *mzrvec, *ww ;
    double y, ylike, yy ;
    int d2, dv, numc, t ; 
    double *yaprob, *yphi, *yweight ;
    int iter ;  
    int mctrials ; 
    double yeff ;

    dv = dim2dv(dim) ;
    d2 = dim*dim ; 

    printf("zfmean + zfvar\n") ;
    printmatw(zfmean, dim, dim, dim) ;
    printnl() ;
    printmatw(zfvar, d2, d2, d2) ;
    printnl() ;

    printf("cmean + cvar trials: %d\n", ntrials) ;
    printcv(coeffs, cvar, dim) ;
    printnl() ;

    mctrials = ntrials ; 

    ZALLOC(yaprob, mctrials, double) ;
    ZALLOC(yphi, mctrials, double) ;
    ZALLOC(yweight, mctrials, double) ;


    ZALLOC(mzrvec, mctrials*dim, double) ;
    ZALLOC(mz, dim, double) ;
    ZALLOC(mzvar, d2, double) ;
    ZALLOC(mzinv, d2, double) ;

    copyarr(coeffs, mz, dim-1) ; 
    squish(mzvar, cvar, dim-1, dim, dim-1) ;
//  printmat(mz, 1, dim-1) ; 
//  printmatl(mzvar, dim-1, dim-1) ;
    pdinv(mzinv, mzvar, dim-1) ;

   for (iter = 0; iter < mctrials; ++iter) { 
    verbose = NO ;
    ww = mzrvec + iter * dim ; 
    genmultgauss(ww, 1, dim-1, mzvar) ; 
    yaprob[iter] = ylike = -0.5*scx(mzinv, NULL, ww, dim-1) ;
    vvp(ww, ww, coeffs, dim-1) ;  // add mean 
    y = 1.0 - asum(ww, dim-1) ; 
    ww[dim-1] = y ; 
    yphi[iter]  = phi2(ww, zfmean, zfvar, dim) ;         
  }
   vvm(yweight, yphi, yaprob, mctrials) ; 
   vmaxmin(yweight, mctrials, &y, NULL) ; 
   vsp(yweight, yweight, -y, mctrials) ; 

// printf("mcweights:\n") ;
// printmatl(yweight, 1, mctrials) ;
// now must exponentiate and balance;  
   vexp(yweight, yweight, mctrials) ; 
   bal1(yweight, mctrials) ; 
   yeff = ess(yweight, mctrials) ;
// printf("zzess: %d %9.3f\n", mctrials, yeff) ;
   yeff /= (double) mctrials ;
// printf ("sampling efficiency:  %12.6f\n", y) ; 
// if (y<.001) printf("efficiency very low.  Bad fit?\n") ;
   vzero(qmean, dim) ;
   vzero(qvar, dim*dim) ;
   for (iter = 0; iter < mctrials; ++iter) { 
    ww = mzrvec + iter * dim ;
    y = yweight[iter] ;
    vst(mz, ww, y, dim) ; 
    vvp(qmean, qmean, mz, dim) ;
   }
// qmean is now mean
   vmaxmin(yaprob, mctrials, &y, NULL) ;
   vsp(yaprob, yaprob, -y, mctrials) ;
   vmaxmin(yphi, mctrials , &y, NULL) ;
   vsp(yphi, yphi, -y, mctrials) ;
   for (iter = 0; iter < mctrials; ++iter) { 
    ww = mzrvec + iter * dim ;
    y = yweight[iter] ;
    vvm(mz, ww, qmean, dim) ;
    vzero(mzvar, dim*dim) ;
    addouter(mzvar, mz, dim) ; 
    vst(mzvar, mzvar, y, dim*dim) ; 
    vvp(qvar, qvar, mzvar, dim*dim) ;
    if (printit) { 
     printf("verboseprint: %6d ", iter) ; 
     printmatx(ww, 1, dim) ; 
     printf(" : ") ; 
     printf("%12.6f ", yaprob[iter]*1000) ;
     printf("%9.3f ",  yphi[iter]*1000) ;
     printf("%12.6f ", yweight[iter] * 1000*1000 ) ;
     printnl() ;
    }
   }

   free(mz) ; 
   free(mzvar) ;
   free(mzinv) ;
   free(mzrvec) ;
   free(yaprob) ;
   free(yphi) ;
   free(yweight) ;

    printf("qmean + qvar trials: %d\n", ntrials) ;
    printcv(qmean, qvar, dim) ;

   return yeff ; 

}


