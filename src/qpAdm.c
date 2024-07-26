#include <stdio.h> // mean s.deeed
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include <nicklib.h>
#include <getpars.h>

#include "globals.h"
#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "regsubs.h"
#include "egsubs.h"
#include "qpsubs.h"
#include "f4rank.h"
#include "eigsubs.h" 


#define WVERSION   "2050" 
// best analysis added
// hires added;  including on summ line
// chrom: 23 added
// allsnps option added
// more hires output added
// nochrom added 
// detail = YES default
// summ: line added
// gendstat added.  Generalized D-stat
// bug in doq4vecb fixed  denominator sometimes wrong... 
// coverage stats addede
// check for nested model.  
// best logic changed if infeasibles 
// better treatment of boundary case (nl=1) 
// instem added
// numchrom added
// format change for large number of sources
// calcadmfix added;
// hiprec_covar added
// fancyf4 added  -- remove bias on f4 calculation
// code cleeanup (Mac) 10/2/19 
// fixed up boundary case (m=n) in ranktest
// mkfstats fstatsname added; made more robust
// keeptmp added (default NO) 
// Cokie Parker speedup + loadsnpx speedup
// doratio added 
// finddworstb added 
// basepop added
// print popsizes 
// mbname added
// numiter -> 50 (ranktest) 
// one call to addscaldiag eliminated (Robert Maier) 
// Monte Carlo importance sampling (phiint) etc
// diagvarplus added inedependent of yscale 
// oldmode/newmode  added ; defaults for yscale changed.  Old default works though

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

double yscale = -1, diagvarplus = -1 ; 

char *parname = NULL;
char *trashdir = "/var/tmp";
int colcaac = YES;
int hires = NO;
int maxrank, rank;
int nsnpused = 0;
int deletefstats = NO ;
char *fstatslog = NULL ; 
int keeptmp = NO ;
int oldmode = YES ;
int newmode = -99 ;
int positive_coeffs = NO ;


Indiv **indivmarkers = NULL;
SNP **snpmarkers;
int numsnps, numindivs;
int seed = 0;
int qtmode = NO;
int chisqmode = NO;		// approx p-value better to use F-stat
int dotpopsmode = YES;
int noxdata = YES;		/* default as pop structure dubious if Males and females */
int popsizelimit = -1;
int gfromp = NO;		// genetic distance from physical 
int xchrom = -1;
int zchrom = -1;
// if bankermode  bankers MUST be in quartet  at most one type 1 in quartet
int lopos, hipos ; 

int allsnps = -99 ;
int oldallsnpsmode = NO ;

int hiprec_covar = NO ;

int ratcoords[2] ;  // ratio 0/1 in bootstrap 
int doratio = NO ;

double blgsize = 0.05;		// block size in Morgans */
double *chitot;
char *popleft, *popright;
char **popllist, **poprlist;
int nleft, nright;

char *basepop = NULL ; 

char *instem = NULL ; 
char *indivname = NULL;
char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *badsnpname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
char *mvname = NULL ;  
int inbreed = -99;
double lambdascale;
double *jmean;
long *wkprint;
int *zprint ; 
int numboot = 1000 ; 
int mctrials = -1 ; 
double *coeffs, *cvar ; 


char *fstatsname = NULL ; 
char *fstatsoutname = NULL ; 
char *tmpdir = NULL ;

int fancyf4 = YES;

char *outputname = NULL;
char *weightname = NULL;
char *blockname = NULL;
FILE *ofile;

double **btop, **bbot, *gtop, *gbot;
int bnblocks, bdim = 0;
int details = YES;

int *ktable = NULL;

void readcommands (int argc, char **argv);

double
doq4vecb (double *ymean, double *yvar, int ***counts, int *bcols,
	  int nrows, int ncols,
	  int lbase, int *llist, int nl, int rbase, int *rlist, int nr,
	  int nblocks);

void
calcevar (double *mean, double *var, double *yvar, int nl, int nr,
	  double **btop, double **bbot, double *gtop, double *gbot,
	  int nblocks, int dim);

void 
calcevarboot(double *mean, double *var, double *ymean, double *yvar, int nl, int nr, int dim, int numboot)  ;

int getf4 (int **xx, int *indx, double *ans);
int getf4old (int **xx, int *indx, double *ans);
void calcadmfix (double *ans, double *A, int n, int *vf) ;
void calcadm (double *ans, double *A, int n);
void printss (double *mean, double *wt, int nl, int nr);
double scorel(double *d, double *V, double *lam, int n)  ;
double gendstat(double *ca, int b1, int b2, int nr, int nl, double *mean, double *var, int *ktable) ;
int isnested(int a, int b) ; 
int usage (char *prog, int exval); 
void loadymv(double *ymean, double *yvar,  char *fstatsname, char **popllist, char **poprlist, int nleft, int nright) ;
void  setktable(int *ktable, int nl,  int nr) ; 
int  mkfstats(char *parname)   ; 
double f4mean(int a, int b, int nl, int nr, double *mean, double *var, int *ktable)   ;
double f4var(int a1, int a2, int b1, int b2, int nl, int nr, double *mean, double *var, int *ktable)  ;
double worstb(double *bb, double *aa, int nl, int nr, double *mean, double *var, int *ktable)   ;
void setpopsize(int *popsize, char **poplist, int numeg, Indiv **indivmarkers, int numindivs)  ;
void writemv(char *mvname, long *kmean, long *kvar, int n) ;
double phiint(double **x1mat, double *rhs, double *cmean, double *cvar, int fsize, int numconstraint)  ;
double tttphi(double *zmean, double *var, int aleft, int aright)  ;
void mcest(double *qcoeffs, double *qvar,  double *coeffs, double *cvar, double *zmean, double *zvar,  int nl, int nr, int mctrials) ;
void printcv(double *coeffs, double *var, int n) ;  

int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int **snppos;
  int *snpindx;
  char **eglist;
  int *nsamppops;
  int numeg;
  int *popsize ; 
  int i, j, k, k1, k2, k3, k4, kk;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;
  double zscore, y1, y2, y, yfeas, ysig, tail, ttail, yy1, yy2, yy;
  int *blstart, *blsize, nblocks;
  int xnblocks;			/* for xsnplist */
  int *bcols;
  double maxgendis;
  int **xtop;
  int npops = 0;
  int nr, nl, minnm, d, dd;
  double *ymean, *yvar, *yvarinv, *ww, *lambda;
  int nindiv = 0, e, f, lag = 1;
  int nignore, numrisks = 1;
  SNP **xsnplist;
  int *tagnums;
  Indiv **xindlist;
  int *xindex, *xtypes;
  int nrows, ncols, m, nc;
  int weightmode = NO;
  int chrom;
  int maxtag = -1;
  int *rawcol;;
  int a, b, x, s, t, col, wt;
  int *xind;
  int xmax, dof;

  int ***counts;
  int jjj, j1, j2;
  int *rlist, lbase, lpop, rbase;
  int *vfix;
  double *yscbest, *ychi, *yinfeasbest;
  int *sbest, *infeassbest, *nfeas, isfeas;
  double *var;
  long *jvar, *jmp ;
  double *w0, *w1, *w2;
  double *wbest;
  int dim, b1, b2;
  double pval ;
  double lfudge ; 
  int numchromp ;  
  int pos ;
  double ymem ; 
  int retkode ; 
  int np ; 
  double *ff3, *ff3var ; 
  double *bworst ;
  double *qcoeffs, *qvar ;
  double bigserr ;


  F4INFO **f4info, *f4pt, *f4pt2, **g4info, *ggpt;


  lopos = -1 ; hipos = BIGINT - 1 ;  
  ivclear(ratcoords, -1, 2) ;

  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;

  printcmdline(argc, argv) ;
  printf ("## qpAdm version: %s\n", WVERSION);
  if (seed==0) seed = seednum() ; 
  printf("seed: %d\n", seed) ;
  SRAND(seed) ;

  if ((yscale < 0) && (diagvarplus < 0)) yscale = 0.0001 ;  // old default
  if (yscale < 0)  yscale = 0.0 ;  // diagvarplus set 
  if (diagvarplus < 0)  diagvarplus = 0.0 ;  

  if (doratio) { 
   ratcoords[0] = 0 ;
   ratcoords[1] = 1 ;
  }
  
  if (newmode == NO) oldmode = YES ;
  if (newmode == YES) oldmode = NO ;
  
  if (tmpdir != NULL) setenv("STMP", tmpdir, 1) ;
  if (keeptmp) deletefstats = NO ;

  if (parname == NULL)
    return 0;

  if (mvname != NULL) printf("mvname: %s\n", mvname) ;

  numchromp = numchrom+ 1 ;

  setfancyf4(fancyf4) ; 

  if (basepop != NULL) printf("basepop set: %s\n", basepop) ;

  if (allsnps == -99) { 
   allsnps = NO ; 
   if (fstatsname == NULL) { 
    printf("allsnps set NO.  It is recommended that allsnps be set explicitly\n") ;
   }
  }

  if (inbreed == -99) { 
   inbreed = allsnps ; 
   if (fstatsname == NULL) {
    printf(" *** recommended that inbreed be explicitly set ***\n") ;
   }
  }

  setinbreed(inbreed) ;

  nleft = numlines (popleft);
  ZALLOC (popllist, nleft, char *);
  nright = numlines (popright);
  ZALLOC (poprlist, nright, char *);
  nleft = loadlist (popllist, popleft);
  nright = loadlist (poprlist, popright);

  if (nleft==0) fatalx("no pops in left!\n") ;
  if (nright==0) fatalx("no pops in right!\n") ;

  if ((allsnps == YES) && (oldallsnpsmode == NO) && (fstatsname == NULL)) {
   if (mkfstats(parname) < 0) { 
    printf("qpfstats failure ... terminating\n") ;
    return -1 ; 
   }; 
  }

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 

  if ((allsnps == YES) && (oldallsnpsmode == YES)) {
   printf("oldallsnpsmode deprecated!\n") ;
  }

  if (outputname != NULL)
    openit (outputname, &ofile, "w");

 if (fstatsname == NULL) { 
  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);
  setindm (indivmarkers);

  numeg = npops = nleft + nright;
  ZALLOC (eglist, npops, char *);
  ZALLOC(popsize, npops, int) ;  
  copystrings (popllist, eglist, nleft);
  copystrings (poprlist, eglist + nleft, nright);

  ckdup (eglist, numeg);

  for (k=0 ; k<numindivs; ++k) {
   indx = indivmarkers[k] ;
   t = indxstring(eglist, numeg, indx -> egroup) ;
   indx -> gkode = t ; // correct when t < 0
  }

  gkignore(indivmarkers, numindivs) ;

  k = getgenos (genotypename, snpmarkers, indivmarkers,
		numsnps, numindivs, nignore);

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    chrom = cupt->chrom;
    pos = nnint(cupt -> physpos) ;

    if ((xchrom > 0) && (chrom != xchrom))
      cupt->ignore = YES;

    if ((xchrom > 0) && (pos < lopos)) cupt -> ignore = YES ; 
    if ((xchrom > 0) && (pos > hipos)) cupt -> ignore = YES ; 

    if ((zchrom > 0) && (chrom == zchrom))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > numchromp)
      cupt->ignore = YES;
    if ((chrom == numchromp) && (xchrom != numchromp))
      cupt->ignore = YES;
  }

 } 

 if (npops == 0) {
  numeg = npops = nleft + nright;
  ZALLOC (eglist, npops, char *);
  ZALLOC(popsize, npops, int) ;  
  copystrings (popllist, eglist, nleft);
  copystrings (poprlist, eglist + nleft, nright);
 }



  printnl ();


  if ((indivmarkers == NULL) && (indivname != NULL))  numindivs = getindivs (indivname, &indivmarkers);
// needed for popsize 

  setpopsize(popsize, eglist, numeg, indivmarkers, numindivs) ;
  printf ("left pops:\n");
  for (k = 0; k < nleft; ++k) {
    printf ("%20s", popllist[k]);
    t = popsize[k] ;
    if (t>=0) printf (" %6d", t) ; 
    printnl() ;
  }

  printnl ();
  printf ("right pops:\n");
  for (k = 0; k < nright; ++k) {
    printf ("%20s", poprlist[k]);
    t = popsize[nleft+k] ;
    if (t>=0) printf (" %6d", t) ; 
    printnl() ;
  }

  printnl ();

  ivlmaxmin(popsize, numeg, NULL, &t) ; 
  if (popsize[t] == 0) fatalx("pop %s has zero samples!\n", eglist[t]) ;


  if (nleft > nright)
    fatalx ("nleft not less than or equal to nright\n");

  for (k = 0; k < nright; ++k) {
    j = indxindex (popllist, nleft, poprlist[k]);
    if (j >= 0)
      fatalx ("population in both left and right: %s\n", poprlist[k]);
  }



  ZALLOC (nsamppops, numeg, int);

  nr = nright - 1;
  nl = nleft - 1;

  ZALLOC (ktable, nleft * nright + 100, int);
  setktable(ktable, nl,  nr) ; 

  d = nl * nr;
  dd = d * d;
  t = dd + d + 10 ; 
  ZALLOC (ymean, d + 10, double);
  ZALLOC (ww, t, double);
  ZALLOC (lambda,  nl, double);
  ZALLOC (yvar, dd, double);


 if (fstatsname == NULL) { 
  for (i = 0; i < numindivs; i++) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    k = indxindex (eglist, numeg, indx->egroup);
    if (k < 0) {
      indx->ignore = YES;
      continue;
    }
    indx->affstatus = YES;
    ++nsamppops[k];
  }

  for (i = 0; i < numeg; i++) {
    t = nsamppops[i];
    if (t == 0)
      fatalx ("No samples: %s\n", eglist[i]);
//    printf ("%3d %20s %4d\n", i, eglist[i], t);
  }

  printf ("jackknife block size: %9.3f\n", blgsize);


  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);
  ZALLOC (rawcol, numindivs, int);
  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);
  if (nrows == 0)
    fatalx ("badbug: no data\n");

  ZALLOC (xtypes, nrows, int);
  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (eglist, numeg, indx->egroup);
    xtypes[i] = k;
  }


  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    cupt->weight = 1.0;
  }
  numsnps = rmsnps (snpmarkers, numsnps, NULL);	//  rid ignorable snps
  if (numsnps == 0)
    fatalx ("no valid snps\n");

  setmgpos (snpmarkers, numsnps, &maxgendis);
  if ((maxgendis < .00001) || (gfromp == YES))
    setgfromp (snpmarkers, numsnps);

  ZALLOC (xsnplist, numsnps, SNP *);
  ZALLOC (tagnums, numsnps, int);

  if (popsizelimit > 0) {
    setplimit (indivmarkers, numindivs, eglist, numeg, popsizelimit);
  }

  ncols = loadsnpxx (xsnplist, snpmarkers, numsnps, xindlist, nrows);

  printf ("snps: %d  indivs: %d\n", ncols, nrows);

 nblocks = setblocksz (&blstart, &blsize, xsnplist, ncols, blgsize, blockname) ;

// loads tagnumber
  printf ("number of blocks for block jackknife: %d\n", nblocks);
  xnblocks = nblocks += 10 ; 

  ZALLOC (counts, ncols, int **);
  for (k = 0; k < ncols; ++k) {
    counts[k] = initarray_2Dint (numeg, 2, 0);
  }

  countpops (counts, xsnplist, xindex, xtypes, nrows, ncols);

  printf("## ncols: %d\n", ncols) ;
  for (i=0; i<numeg; ++i) { 
   t = 0 ; 
   for (a=0; a<ncols; ++a) { 
    s = counts[a][i][0] + counts[a][i][1]  ;
    if (s>0) ++t ;
   }
   printf("coverage: %20s %6d\n", eglist[i] , t) ;
  }

  ZALLOC (bcols, ncols, int);	// blocks for columns -1 => unused
  ivclear (bcols, -1, ncols);
  for (k = 0; k < ncols; k++) {
    cupt = xsnplist[k];
    t = bcols[k] = cupt->tagnumber;
    if (t>=xnblocks) fatalx("bad logic bug: %d %d %d\n", k, t, xnblocks) ;
  }
 }

  lbase = 0;
  rbase = nleft;
  ZALLOC (rlist, nright - 1, int);
  for (a = 1; a < nright; ++a) {
    rlist[a - 1] = nleft + a;
  }

  ZALLOC (xind, nl, int);
  for (jjj = 1; jjj < nleft; ++jjj) {
    xind[jjj - 1] = jjj;
  }

/**
 mnmin = MIN(nl, nr) ;
 t = MIN(mnmin, 4) ;
*/


  minnm = MIN (nl, nr);

  maxrank = nl;
  rank = nl - 1;

  ZALLOC (f4info, maxrank + 1, F4INFO *);
  for (x = 0; x <= maxrank; ++x) {
    ZALLOC (f4info[x], 1, F4INFO);
    f4info_init (f4info[x], nl, nr, popllist, poprlist, x);
  }


 if (fstatsname == NULL) { 
  y = doq4vecb (ymean, yvar, counts, bcols,
		nrows, ncols, lbase, xind, nl, rbase, rlist, nr, xnblocks);
// y is jackknife dof 
  printf ("Effective number of blocks: %9.3f\n", y);
  printf ("numsnps used: %d\n", nsnpused);

 }

 else {
  loadymv(ymean, yvar,  fstatsname, popllist, poprlist, nleft,  nright) ; 
//  printf("loadymv exited\n") ;
 }

 retkode = checkmv(ymean, yvar, nl, nr) ; 

 if (retkode == -1) printf("f4 stats all zero.  Rank 0!.  Aborting run\n") ;
 if (retkode == -2) printf("f4 variance absurdly small.  Aborting run\n") ;

 if (retkode < 0) return -1 ;
 if (diagvarplus < 0) diagvarplus = yscale ;
 addscaldiag(yvar, diagvarplus, nl*nr) ;    

  ZALLOC (vfix, nl, int);

  for (x = rank; x <= maxrank; ++x) {
    f4pt = f4info[x];
    f4pt->dofjack = y;
    doranktest (ymean, yvar, nl, nr, x, f4pt);
    if (x > 0) {
      f4pt2 = f4info[x - 1];
      f4pt->dofdiff = f4pt2->dof - f4pt->dof;
      f4pt->chisqdiff = f4pt2->chisq - f4pt->chisq;
    }
  }


  ZALLOC (jmean, nl, double);
  ZALLOC (wkprint, nl * nl, long);
  ZALLOC (zprint, nl * nl, int);
  ZALLOC (coeffs, nl, double);
  ZALLOC (cvar, nl*nl, double);

  for (x = rank; x <= maxrank; ++x) {
    if (x == rank)
      printf ("codimension 1\n");
    if (x == maxrank)
      printf ("full rank\n");
    f4pt = f4info[x];
    printf4info (f4pt);
    pval = rtlchsq (f4pt->dofdiff, f4pt->chisqdiff) ; 
    printnl ();
/**
    if (x==maxrank) continue ;
    printf("++++++++ %d\n", maxrank-x)  ; 
    scorel(f4pt -> resid, yvar, &lfudge, nl*nr) ; // eigenanalysis  
*/
  }

  f4pt = f4info[rank];
  calcadm (ww, f4pt->A, nl);
  printf ("best coefficients: ");
  if (hires)
    printmatl (ww, 1, nl);
  else
    printmatw (ww, 1, nl, nl);

  ZALLOC (wbest, nl, double);
  ZALLOC (w0, nl * nr, double);
  ZALLOC (w1, nl * nr, double);
  ZALLOC (w2, nl * nr, double);

  copyarr (ww, wbest, nl);
  ZALLOC (var, nl * nl, double);
  ZALLOC (jvar, nl * nl, long);
  ZALLOC (jmp, nl, long);

//  printss (ymean, ww, nl, nr);

  

  if (fstatsname == NULL) { 
   calcevar (jmean, var, yvar, nl, nr, btop, bbot, gtop, gbot, bnblocks, bdim);
  }
  else { 
 //printf("zzrandom %d\n", ranmod(1000*1000*1000)) ;
   calcevarboot(jmean, var, ymean, yvar, nl, nr, nl*nr, numboot) ; 
   printf("bootstrap saimpling of F-coeffs: %d\n", numboot) ;
   printf("zzjmean ") ;  printmatw(jmean, 1, nl, nl) ;
  }
  copyarr(jmean, coeffs, nl) ;
  copyarr(var, cvar, nl*nl) ;
  getdiag (ww, var, nl);
  vsp (ww, ww, 1.0e-20, nl);
  vsqrt (ww, ww, nl);
  vmaxmin(ww, nl, &bigserr, NULL) ;
  printf ("      std. errors: ");
  
  if (hires)
    printmatwl (ww, 1, nl, nl);
  else
    printmatw (ww, 1, nl, nl);


  printnl ();

  vst(w0, wbest, 1000, nl) ; 
  fixitl(jmp, w0, nl) ; 
  vst(w0, var, 1.0e6, nl*nl) ;  
  fixitl(jvar, w0, nl*nl) ;  
  writemv(mvname, jmp, jvar, nl) ;

  if (hiprec_covar == NO) {

   vst (ww, var, 1.0e6, nl * nl);

/**
   printf("zzz\n") ; 
   printmatl(ww, nl, nl) ;
   printnl() ;
   printmatl(var, nl, nl) ;
   printnl() ;
*/

   fixitl (jvar, ww, nl * nl);
   printf ("error covariance (* 1,000,000)\n");
   printlmat (jvar, nl, nl);
  }

  else {
   vst (ww, var, 1.0e9, nl * nl);
   fixitl (jvar, ww, nl * nl);
   printf ("error covariance (* 1,000,000,000)\n");
   printlmat (jvar, nl, nl);

   eigvecs (var, lambda, ww, nl);
   printf("eigenanalysis:\n") ;
   vst(lambda, lambda, 1.0e9, nl) ;
   printf("evals: (*10e9)\n") ;  printmat(lambda, 1, nl) ;
   printnl() ;
   printf("evecs:\n") ;  printmat(ww, nl, nl) ;
   printnl() ;
  }

  if (hires) {
    printf ("hires: ");
    vst (ww, jmean, 1.0e6, nl);
    fixit (zprint, ww, nl);
    printimatx(zprint, 1, nl);
    printf ("  ");
    t = mktriang (ww, var, nl);
    vst (ww, ww, 1.0e8, t);
    fixitl (wkprint, ww, t);
    printlmat (wkprint, 1, t);
    printf ("[mean *1.0e6, var *1.0e8]\n");
  }

  printnl ();
  printnl ();

  if ((oldmode == NO) && (bigserr > 1.0)) { 
   printf("large std. error on coeffs.  Reverting to oldmode\n") ;
   oldmode = YES ;
  }
         

  if (oldmode == NO) printf("summ (old):") ;  // we will write summ: later 
  if (oldmode == YES) printf("summ:") ;
  printf(" %s ", popllist[0]) ; 
  printf(" %3d ", nl) ; 
  printf(" %12.6f ", pval) ; 
  if (hires) {
   printmatwlx(jmean, 1, nl, nl) ; 
  }
  else {
   printmatwx(jmean, 1, nl, nl) ; 
  }
  t = mktriang (ww, var, nl);
  vst (ww, ww, 1.0e6, t);
  fixitl (wkprint, ww, t);
  printlmat (wkprint, 1, t);
  printnl ();
  
  if (nl == 1) { 
   printf("single source. terminating\n") ; 
   printf("run qpDstat if you want Z scores for f4 stats\n") ;
   printf ("## end of run\n");
   return 0 ;
  }

//  we now loop over fix patterns.  
  xmax = pow (2, nl) - 1;
  ZALLOC (g4info, xmax + 2, F4INFO *);
  for (x = 0; x <= xmax; ++x) {
    ZALLOC (g4info[x], 1, F4INFO);
    f4info_init (g4info[x], nl, nr, popllist, poprlist, nl - 1);
  }

  ZALLOC (sbest, nl, int);
  ZALLOC (nfeas, nl, int);
  ZALLOC (yscbest, nl, double);
  ZALLOC (ychi, nl, double);
  vclear (yscbest, -1.0e5, nl);	// p-value large is good
  printf (" %12s  ", "fixed pat");
  printf ("wt  %3s", "dof");
  printf (" %9s", "chisq");
  printf (" %15s", "tail prob");
  printnl ();
  for (wt = 0; wt < nl; ++wt) {
    for (k = 0; k <= xmax; ++k) {
      f4pt = g4info[k];
      dekodeitb (vfix, k, nl, 2);
      t = intsum (vfix, nl);
      if (t != wt)
	continue;
      if (wt == 0) {
	f4pt = g4info[0] = f4info[nl - 1];
      }
      else {
        doranktestfix (ymean, yvar, nl, nr, nl - 1, f4pt, vfix);
      }
      if (wt==0) calcadm (ww, f4pt->A, nl);
      else calcadmfix(ww, f4pt -> A, nl, vfix) ;  
      vmaxmin (ww, nl, NULL, &y);
      yfeas = y ; 
      dof = nnint (f4pt->dof);
      if (dof > 0) {
	ttail = tail = rtlchsq (dof, f4pt->chisq);
      }
      else {
	ttail = tail = 0;
      }

      isfeas = NO ; 
      if (yfeas >= -0.001) { 
       isfeas = YES ; 
       ++nfeas[wt] ;
      }
     
      if (isfeas == NO)  { 
	tail = 0;		// not feasible
      }
      if (ttail < 1.0e-30)
	ttail =  0;

      if (nfeas[wt] == 1) { 
	yscbest[wt] = tail;
	sbest[wt] = k;
	ychi[wt] = f4pt->chisq;
      }

     
      if ((nfeas[wt] == 0) && (ttail>yscbest[wt]))  { 
	yscbest[wt] = ttail;
	sbest[wt] = k;
	ychi[wt] = f4pt->chisq;
      }

      if (tail > yscbest[wt])  {
	yscbest[wt] = tail;
	sbest[wt] = k;
	ychi[wt] = f4pt->chisq;
      }

      printf (" %12s  %1d   %3d %9.3f %15.6g ", binary_string (k, nl), wt,
	      dof, f4pt->chisq, ttail);
      printmatwx (ww, 1, nl, nl);
      if (isfeas == NO) 
	printf (" infeasible");
      printnl ();
      if (verbose) {
	printf4info (f4pt);
	printnl ();
      }
/**
   if (y>yscbest[wt]) { 
    sbest[wt] = k ;
    yscbest[wt] = y ;
   }
*/
    }
  }

  for (wt = 0; wt < nl; ++wt) {
    k = sbest[wt];
    printf ("best pat: %12s ", binary_string (k, nl));
    printf (" %15.6g ", yscbest[wt]);
    if (wt == 0) {
      printf ("             -  -\n");
      continue;
    }
    y = ychi[wt] - ychi[wt - 1];
    if (isnested(sbest[wt-1], sbest[wt])) {  
     tail = rtlchsq (1, y);
     printf (" chi(nested): %9.3f p-value for nested model: %15.6g", y,
	    tail);
    }
    else { 
     printf("not nested") ;
    }
    if (nfeas[wt] == 0) printf(" infeasible") ; 
    printnl() ;
  }
  printnl ();
  fflush (stdout);
  if (details) {
    printf ("coeffs: ");
    printmatw (wbest, 1, nl, nl);
    printnl() ;
    printf("## dscore:: f_4(Base, Fit, Rbase, right2)\n") ; 
    printf("## genstat:: f_4(Base, Fit, right1, right2)\n") ; 
    printnl() ;
    dim = nl * nr;
    for (b = 0; b < nr; ++b) {
      vzero (w0, nl * nr);
      vzero (w1, nl * nr);
      for (a = 0; a < nl; ++a) {
	k = ktable[a * nr + b];
	printf ("details: ");
	printf ("%20s %20s ", popllist[a + 1], poprlist[b + 1]);
	y1 = w0[a] = ymean[k];
	w1[k] = wbest[a];
	printf ("%12.6f", y1);
	yy = yvar[k * dim + k];
	ysig = sqrt (yy);
	printf ("%12.6f", y1 / ysig);	// Z-score
	printnl ();
      }
      printf ("dscore: %20s ", poprlist[b + 1]);
      y1 = vdot (w0, wbest, nl);
// w0 is empirical f4;  wbest coeffs;  ideally should be zero. 
      mulmat (w2, yvar, w1, dim, dim, 1);
      y2 = vdot (w1, w2, nl * nr);	// variance
      ysig = y1 / sqrt (y2);
      printf ("f4: %12.6f Z: %12.6f\n", y1, ysig);
      printnl ();
    }
    for (b1 = -1; b1 < nr; ++b1) { 
     for (b2=b1+1; b2 < nr; ++b2) { 
      y1 =  gendstat(wbest, b1, b2, nl, nr, ymean, yvar, ktable)  ;  
      printf("gendstat: ") ; 
      printf("%20s ", poprlist[b1+1]) ;
      printf("%20s ", poprlist[b2+1]) ;
      printf("%9.3f", y1) ;
      printnl() ; 
     }
    }
    printnl() ;
    
   ZALLOC(bworst, nr, double) ;
   y = worstb(bworst, wbest,  nl,  nr, ymean, yvar, ktable)   ;
   printf("worst Z-score with right hand mix\n") ;
   printf("f4(Target, Fit, Base, mix of Right pops;  Z: %9.3f sum: %9.3f\n", y, asum(bworst, nr))  ; 
   for (k=0; k<nr; ++k) {  
    printf("%30s ", poprlist[k+1]) ; 
    printf("%9.3f", bworst[k]) ;
    printnl() ;
   }
   printnl() ;
  }


/**
  for (a=0; a<nl; ++a) { 
   for (b=0; b<nr; ++b) { 
   y1 = f4mean(a, b, nl, nr, ymean, yvar, ktable) ;
   y2 = f4var(a, a, b, b, nl, nr, ymean, yvar, ktable) ;
   printf("zzcheck %s %s ", popllist[a+1], poprlist[b+1]) ;
   ysig = y1 / sqrt (y2+1.0e-20);
   printf ("f4: %12.6f Z: %12.6f", y1, ysig);
   printnl ();
 }}
*/

  if (deletefstats) { 
// fstatsname has been created  
   printf("removing %s\n", fstatsname) ;
   remove(fstatsname) ; 
  }

  if (oldmode == YES) {
   printf("oldmode set: terminating\n") ;
   ymem = calcmem(1)/1.0e6 ;
   printf("##end of qpAdm: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
   return 0;
  }
  
 


  ZALLOC(qcoeffs, nl, double) ;
  ZALLOC(qvar, nl*nl, double) ;
  
  copyarr(jmean, qcoeffs, nl) ;
  copyarr(var, qvar, nl*nl) ;

 printf("oldmode: \n") ;
   printcv(qcoeffs, qvar, nl) ;

  printnl() ;  
  if (mctrials < 0) { 
   mctrials = MIN(pow(2, nl), 10000) ; 
   mctrials = MAX(mctrials, 1000) ;
  }

  
  mcest(qcoeffs, qvar,  coeffs, cvar, ymean, yvar, nl, nr, mctrials) ;
  printf("mctrials: %d\n", mctrials) ;  


  printf("mcest (Bayesian):\n") ;
    printcv(qcoeffs, qvar, nl) ;

   getdiag(ww, qvar, nl); 
   vsqrt(ww, ww, nl) ;


  printf("Target: %s\n", eglist[0]) ;
  for (i=1; i<nleft; ++i) { 
   printf("%25s ", eglist[i]) ;
   j = i-1; 
   printf("%9.3f  ", qcoeffs[j]) ;
   printf("%12.6f", ww[j]) ; 
   printnl() ;
  }

  printf("summ:") ;
  printf(" %s ", popllist[0]) ; 
  printf(" %3d ", nl) ; 
  printf(" %12.6f ", pval) ; 
  if (hires) {
   printmatwlx(qcoeffs, 1, nl, nl) ; 
  }
  else {
   printmatwx(qcoeffs, 1, nl, nl) ; 
  }
  t = mktriang (ww, qvar, nl);
  vst (ww, ww, 1.0e6, t);
  fixitl (wkprint, ww, t);
  printlmat (wkprint, 1, t);
  printnl ();
  
  fflush(stdout) ;
//tttphi(ymean, yvar, nl, nr)  ;

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qpAdm: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;
}

void addco(double *co, double yadd, int n) 
{
// special purpose. Last diagonal element not changed 
 int i, k ; 

 for (i=0; i<(n-1) ; ++i) {  
  k = i*n + i ; 
  co[k] += yadd ; 
 }
}


void
calcadmfix (double *ans, double *A, int n, int *vf)
{
  double *coeff, *rhs,  y, ytrace, *ctran, *prod, *fvals, *rr ;
  int t, baditer = 0, *vfix, nfix, k ;

  if (A == NULL) {
    vzero (ans, n);
    return;
  }

  if (n==1) { 
   ans[0] = 1 ; 
   return ; 
  }
  ZALLOC (coeff, n * n, double);
  ZALLOC (ctran, n * n, double);
  ZALLOC (prod, n * n, double);
  ZALLOC (rhs, n, double);
  ZALLOC (rr, n, double);
  ZALLOC (fvals, n, double);
  ZALLOC (vfix, n, int);

  transpose (coeff, A, n, n - 1);
  ytrace = trace(coeff, n) ;
  vclear (coeff + n * (n - 1), 1.0, n);
  rhs[n - 1] = 1;
  vzero(ans, n) ;
  transpose(ctran, coeff, n, n) ;

  mulmat(prod, ctran, coeff, n, n, n) ; 
  mulmat(rr, ctran, rhs, n, n, 1) ; 

  nfix = 0 ; 
  for (k=0; k<n; ++k) { 
   if (vf[k] == 1) {  
    vfix[nfix] = k ;  
    ++nfix ; 
   }
  }
  

  for (;;) { 
   if (baditer>=10) break ;
   t = solvitfix (prod, rr, n, ans, vfix, fvals, nfix);
   y = asum(ans, n) ;  
   if (!isfinite(y)) t = 1 ; 
   if (t==0) break ; 
    ++baditer ; 
    if (verbose) printf("singular matrix: %3d %12.6f\n", baditer, ytrace) ;
    addco(coeff, ytrace*.001, n) ;
  }
  if (asum(ans, n) > 0) bal1(ans, n) ;

  free (coeff);
  free (rhs);
  free(ctran) ; 
  free(prod) ; 
  free(rr) ;  
  free(fvals) ; 
  free(vfix) ; 

}


void
calcadm (double *ans, double *A, int n)
{
  double *coeff, *rhs,  y, ytrace ;
  int t, baditer = 0 ;

  if (A == NULL) {
    vzero (ans, n);
    return;
  }

  if (n==1) { 
   ans[0] = 1 ; 
   return ; 
  }
  ZALLOC (coeff, n * n, double);
  ZALLOC (rhs, n, double);

  transpose (coeff, A, n, n - 1);
  ytrace = trace(coeff, n) ;
  vclear (coeff + n * (n - 1), 1.0, n);
  rhs[n - 1] = 1;
  vzero(ans, n) ;

  

  for (;;) { 
   if (baditer>=10) break ;
   t = linsolv (n, coeff, rhs, ans);
   y = asum(ans, n) ;  
   if (!isfinite(y)) t = 1 ; 
   if (t==0) break ; 
    ++baditer ; 
    if (verbose) printf("singular matrix: %3d %12.6f\n", baditer, ytrace) ;
    addco(coeff, ytrace*.001, n) ;
  }
  if (asum(ans, n) > 0) bal1(ans, n) ;

  free (coeff);
  free (rhs);

}

void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n;

  if (argc == 1) { usage(basename(argv[0]), 1); }
  while ((i = getopt (argc, argv, "p:vVdh")) != -1) {

    switch (i) {

    case 'h':
      usage(basename(argv[0]), 0);
      break ; 

    case 'p':
      parname = strdup (optarg);
      break;

    case 'v':
      printf ("version: %s\n", WVERSION);
      break;

    case 'V':
      verbose = YES;
      break;

    case 'd':
      details = YES;
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }


  if (parname == NULL) {
    fprintf (stderr, "no parameters\n");
    return;
  }

  pcheck (parname, 'p');
  printf ("%s: parameter file: %s\n", argv[0], parname);
  ph = openpars (parname);
  dostrsub (ph);

  getstring (ph, "instem:", &instem);
  getstring (ph, "genotypename:", &genotypename);
  getstring (ph, "snpname:", &snpname);
  getstring (ph, "indivname:", &indivname);
  getstring (ph, "popleft:", &popleft);
  getstring (ph, "popright:", &popright);
  getstring (ph, "basepop:", &basepop);

  getstring (ph, "output:", &outputname);
  getstring (ph, "outputname:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "blockname:", &blockname);
  getstring (ph, "fstatsname:", &fstatsname);
  getstring (ph, "fstatsoutname:", &fstatsoutname);
  getstring (ph, "fstatslog:", &fstatslog);
  getstring (ph, "tmpdir:", &tmpdir);
  getstring (ph, "mvname:", &mvname);

  getdbl (ph, "blgsize:", &blgsize);

  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "chrom:", &xchrom);
  getint (ph, "nochrom:", &zchrom);
  getint (ph, "numchrom:", &numchrom);
  getint (ph, "lopos:", &lopos);
  getint (ph, "hipos:", &hipos);
  getint (ph, "maxrank:", &maxrank);
  getint (ph, "hires:", &hires);
  getint (ph, "allsnps:", &allsnps);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "fancyf4:", &fancyf4);
  getint (ph, "useallsnps:", &allsnps);
  getint (ph, "keeptmp:", &keeptmp);
  getint (ph, "details:", &details);
  getint (ph, "seed:", &seed) ; 
  getint (ph, "hiprec_covar:", &hiprec_covar) ; 
  getint (ph, "numboot:", &numboot) ; 
  getint (ph, "doratio:", &doratio) ; 
  getint (ph, "mctrials:", &mctrials) ; 
  getint (ph, "oldmode:", &oldmode) ; 
  getint (ph, "newmode:", &newmode) ; 
  getint (ph, "positive_coeffs:", &positive_coeffs) ; 
  getdbl (ph, "diagplus:", &yscale);
  getdbl (ph, "diagvarplus:", &diagvarplus);



printf ("### THE INPUT PARAMETERS\n");
printf ("##PARAMETER NAME: VALUE\n");
writepars (ph);

}

void
setktable(int *ktable,  int nl,  int nr) 
{
int a, b, k ; 
k = 0;
for (a = 0; a < nl; ++a) {
 for (b = 0; b < nr; ++b) {
  ktable[a * nr + b] = k;
  ++k;
 }
}

}

double
doq4vecb (double *ymean, double *yvar, int ***counts, int *bcols,
  int nrows, int ncols, int lbase, int *llist, int nl, int rbase,
  int *rlist, int nr, int nblocks)
// return dof estimate
{

int a, b, c, d;
int k, g, i, col, j;
double ya, yb, y  ;
double *top, *bot, *wjack, *wbot, *wtop, **bjtop;
double *bwt;
double ytop, ybot;
double y1, y2, yscal;
double *w1, *w2, *ww, m1, m2;
double *mean;
static int nbad = 0 ; 

int bnum, totnum;
int ret, ret2;
double *f4, *f4bot;
int dim = nl * nr;
int isok;
int **xtop, *xt ;  

if (nrows == 0)
fatalx ("badbug\n");


ZALLOC (f4, dim, double);
ZALLOC (f4bot, dim, double);
ZALLOC (mean, dim, double);
ZALLOC (wjack, nblocks, double);
ZALLOC (wtop, dim, double);
ZALLOC (wbot, dim, double);

ZALLOC (gtop, dim, double);
ZALLOC (gbot, dim, double);

btop = initarray_2Ddouble (nblocks, dim, 0);
bbot = initarray_2Ddouble (nblocks, dim, 0);

bjtop = initarray_2Ddouble (nblocks, dim, 0);
xtop = initarray_2Dint (dim, 4, 0);

k = 0;
for (a = 0; a < nl; ++a) {
for (b = 0; b < nr; ++b) {
xt = xtop[k];
xt[0] = lbase;
xt[1] = llist[a];
xt[2] = rbase;
xt[3] = rlist[b];
++k;
// printf("xtk: %3d ", k) ; printimat(xt, 1, 4) ;
}
}

totnum = 0;
for (col = 0; col < ncols; ++col) {
bnum = bcols[col];
if (bnum < 0)
continue;
if (bnum >= nblocks)
 fatalx ("logic bug:: %d %d %d\n", col, bnum, nblocks);

vzero (f4, dim);
vzero (f4bot, dim);

isok = YES;
for (a = 0; a < dim; ++a) {
ret = getf4 (counts[col], xtop[a], &y);

/**
ret2 = getf4old (counts[col], xtop[a], &y2);
if ((ret != ret2) || (fabs(y-y2) > .001)) { 
 ++nbad ; 
 if (nbad < 100) printf("zzbad %d %d  :: $%d %d :: %9.3f %9.3f\n", col, a, ret, ret2, y, y2) ; 
}
*/

if (ret==2) {
 printf("bad quad: "); printimat(xtop[a], 1, 4) ;
 fatalx("bad quad\n") ;
}

f4[a] = y ;  
f4bot[a] = 1 ; 
if (ret<0) {
 f4[a] = f4bot[a] =  0 ;
}
 if ((allsnps == NO) && (ret < 0))  isok = NO;
}
if (allsnps == YES) isok = YES;
if (isok == NO) continue;
vvp (btop[bnum], btop[bnum], f4, dim);
vvp (bbot[bnum], bbot[bnum], f4bot, dim);
++wjack[bnum];
++totnum;
}

sum2D (gtop, btop, nblocks, dim);
sum2D (gbot, bbot, nblocks, dim);

vsp (gbot, gbot, 1.0e-20, dim);
vvd (mean, gtop, gbot, dim);

for (k = 0; k < nblocks; k++) {
top = btop[k];
bot = bbot[k];
vvm (wtop, gtop, top, dim);
vvm (wbot, gbot, bot, dim);
vvd (bjtop[k], wtop, wbot, dim);
}


free (f4);
free (f4bot);
free (wtop);
free (wbot);

for (k = 0; k < nblocks; ++k) {
if (allsnps == NO) break;
wjack[k] = asum (bbot[k], dim);
}

wjackvest (ymean, yvar, dim, mean, bjtop, wjack, nblocks);
if (verbose)  { 
 printf("wjackvest! %d\n", dim) ;
 printmatl (ymean, 1, dim) ; 
 printnl() ; 
 printmatl (yvar, dim, dim) ; 
 for (k=0; k<dim; ++k) { 
  printf("diag: %3d %15.9f\n", k, yvar[k*dim+k]) ;
 }
}

// addscaldiag(yvar, yscale, dim) ;    

bnblocks = nblocks;
bdim = dim;

y1 = asum (wjack, nblocks);
y2 = asum2 (wjack, nblocks);

free (wjack);
free (mean);
free2D (&bjtop, nblocks);
free2Dint (&xtop, dim);

nsnpused = totnum;
return (y1 * y1) / y2;	// a natural estimate for degrees of freedom in F-test. 

}


int
getf4old (int **xx, int *indx, double *ans)
{

int a, i;
double ya, yb,  y0, y1, ytot, ff[4];
double h0, h1;
int isok, f4mode ;
// f4mode == NO => f2, or f3

*ans = 0.0;
if (indx == NULL) {
*ans = 1.0;
return 2;
}

isok = f4mode = YES ; 

 if (indx[0] == indx[2]) f4mode = NO ;
 if (indx[0] == indx[3]) f4mode = NO ;
 if (indx[1] == indx[2]) f4mode = NO ;
 if (indx[1] == indx[3]) f4mode = NO ;

for (i = 0; i < 4; ++i) {
 ff[i] = -10*(i+1) ;  // silly value ;
 a = indx[i];
 if (a < 0) {
  *ans = 1.0;
  return 2;
 }


 // *ans - 0 ;
 y0 = (double) xx[a][0];
 y1 = (double) xx[a][1];
 ytot = y0 + y1;
 if (ytot <= 0.0) {
  isok = NO ;
  continue ;
 }
 ff[i] = y0 / ytot;
}

if ((isok == NO) && (f4mode == NO)) return -1 ;

ya = fabs(ff[0]-ff[1])  ;  
yb = fabs(ff[2]-ff[3])  ;  

if (f4mode && (MIN(ya, yb) < .00001)) { 
 return 1 ; 
}

if (isok == NO) return -1 ;

*ans = (ff[0] - ff[1]) * (ff[2] - ff[3]);
if (f4mode == YES) return 1 ;

a = indx[0];
h0 = hfix (xx[a]);
a = indx[1];
h1 = hfix (xx[a]);
if (indx[0] == indx[2])
*ans -= h0;
if (indx[0] == indx[3])
*ans += h0;
if (indx[1] == indx[3])
*ans -= h1;
if (indx[1] == indx[2])
*ans += h1;
return 1;

}
void mv2D(double *mm, double *vv, double **data, int n, int dim) 
{
  double y ; 
  double *ww ;  
  int k ; 

  if (n<=2) fatalx("(mv2D) number of rows too small: %d\n", n) ; 
  y = (double) n ; 

  sum2D(mm, data, n, dim) ; 
  vst(mm, mm, 1.0/y, dim) ; 
 
  ZALLOC(ww, dim, double) ; 
  vzero(vv, dim*dim) ; 
  for (k=0; k<n; ++k) { 
   vvm(ww, data[k], mm, dim) ; 
   addouter(vv, ww, dim) ;
  }
  
  y = (double) (n-1) ; 
  vst(vv, vv, 1.0/y, dim*dim) ; 

 free(ww) ; 
}
void 
calcevarboot(double *bootmean, double *bootvar, double *ymean, double *yvar, int nl, int nr, int dim, int numboot)  
{

 double *rvec, *rp, *svec,  *ww, *mean, **tmean, *totmean  ; 
 int k, ka, kb ; 
 F4INFO *f4wk, f4tt;
 double yn ; 
 double *wwcovar ;
 double *rathist ; 
 double ya, yb, rat, yratmean, yratsig ;

 ZALLOC(rvec, 2*dim*numboot, double) ; 
 ZALLOC(ww, dim, double) ; 
 ZALLOC(rathist, 2*numboot, double) ; 
 ZALLOC(wwcovar, dim*dim, double) ; 
 ZALLOC(totmean, nl, double) ; 
 tmean = initarray_2Ddouble(2*numboot, nl, 0) ; 

 genmultgauss(rvec, numboot, dim, yvar) ; 
 svec = rvec + numboot*dim  ; 
 for (k=0; k<numboot; ++k) { 
  rp = rvec + k*dim ; 
  copyarr(rp, ww, dim) ; 
  vvp(rp, rp, ymean, dim) ; 
  vvm(ww, ymean , ww, dim) ; 
  copyarr(ww, svec + k*dim, dim) ; 
 }
 f4wk = &f4tt;
 f4info_init (f4wk, nl, nr, popllist, poprlist, nl - 1);

 doranktest (ymean, yvar, nl, nr, nl - 1, f4wk);
 calcadm(totmean, f4wk -> A, nl) ;
 printf("totmean:  ") ; 

 if (hires == NO) printmat(totmean, 1, nl) ; 
 else printmatl(totmean, 1, nl) ; 

 for (k=0; k<numboot; ++k) { 
  mean = rvec + k*dim ; 
  doranktest (mean, yvar, nl, nr, nl - 1, f4wk);
  calcadm (tmean[2*k], f4wk->A, nl);
  if (k==0) { 
   printf("zzevarboot\n") ; 
   printmatl(tmean[2*k], 1, nl) ;
   printmatl(f4wk->A, nl, nl-1) ;
  }
  mean = svec + k*dim ; 
  doranktest (mean, yvar, nl, nr, nl - 1, f4wk);
  calcadm (tmean[2*k+1], f4wk->A, nl);
  if (k==0) { 
   printf("zzevarboot2\n") ; 
   printmatl(tmean[2*k+1], 1, nl) ;
   printmatl(f4wk->A, nl, nl-1) ;
  }
 }

 ka = ratcoords[0] ;
 kb = ratcoords[1] ;
 if (kb>=0)  {        
  for (k=0; k< 2*numboot; ++k) { 
   ya = tmean[k][0] ;
   yb = tmean[k][1] ;
   yb = MAX(yb, 1.0e-10) ; 
   rathist[k] = ya/yb ;
   printf("zz %d ", k) ; 
   printf("%9.3f ", ya) ;
   printf("%9.3f ", yb) ;
   printf("%9.3f ", rathist[k]) ;
   printnl() ; 
  }
  calcms(rathist, 2*numboot, &yratmean, &yratsig) ; // mean s.dev
  printf("coeff ratio: %d %d ", ka, kb) ; 
  printf(" %9.3f %9.3f\n", yratmean, yratsig) ;
 }

 yn = (double) (2*numboot) ; 
 mv2D(ww, wwcovar, tmean, 2*numboot, nl) ; 
 printf("boot mean: ") ; 

 if (hires == NO) printmat(ww, 1, nl) ; 
 else printmatl(ww, 1, nl) ; 

 copyarr(ww, bootmean, nl) ;
 copyarr(wwcovar, bootvar, nl*nl) ; 


 free(ww) ; 
 free(rathist) ; 
 free(wwcovar) ; 
 free(rvec) ; 
 free(totmean) ; 
 free2D(&tmean, 2*numboot) ; 
}

void 
estrat(double *mean, double *var, int dim, double *prat,  double *psig, int numboot) 
{
 
 double *rvec, *rp, *svec,  *ww ; 
 int k, ka, kb ; 
 double *rathist ; 
 double ya, yb, rat, yratmean, yratsig ;

 ZALLOC(rvec, 2*dim*numboot, double) ; 
 ZALLOC(ww, dim, double) ; 
 ZALLOC(rathist, 2*numboot, double) ; 

 genmultgauss(rvec, numboot, dim, var) ; 
 svec = rvec + numboot*dim  ; 


 svec = rvec + numboot*dim ; 
 for (k=0; k<numboot; ++k) { 
  rp = rvec + k*dim ; 
  copyarr(rp, ww, dim) ; 
  vvp(rp, rp, mean, dim) ; 
  vvm(ww, mean , ww, dim) ; 
  copyarr(ww, svec + k*dim, dim) ; 
 }

 ka = ratcoords[0] ;
 kb = ratcoords[1] ;
 if (kb>=0)  {        
  for (k=0; k< numboot; ++k) { 
   rp = rvec + k*dim ; 
   ya = rp[0] ;
   yb = rp[1] ;
   yb = MAX(yb, 1.0e-10) ; 
   rathist[2*k] = ya/yb ;
   rp = svec + k*dim ; 
   ya = rp[0] ;
   yb = rp[1] ;
   yb = MAX(yb, 1.0e-10) ; 
   rathist[2*k+1] = ya/yb ;
  }
  calcms(rathist, 2*numboot, prat, psig) ; // mean s.dev
  printf("coeff ratio: %d %d ", ka, kb) ; 
  printf(" %9.3f %9.3f\n", *prat, *psig) ;
 }

 free(rvec) ; 
 free(ww) ; 
 free(rathist) ; 


}


void
calcevar (double  *jmean, double *jvar, double *yvar, int nl, int nr, double **btop,
  double **bbot, double *gtop, double *gbot, int nblocks, int dim)
{

double *ymean, *mean, *wtop, *wbot, *totmean, **tmean, *wjack;
F4INFO *f4wk, f4tt;
int k;
double y, yratmean, yratsig ; 

if (dim <= 0)
fatalx ("(calcevar) dimension unset\n");

ZALLOC (mean, dim, double);
ZALLOC (ymean, dim, double);
ZALLOC (wtop, dim, double);
ZALLOC (wbot, dim, double);

ZALLOC (totmean, nl, double);
tmean = initarray_2Ddouble (nblocks, nl, 0);

// addscaldiag(yvar, yscale, dim) ;

f4wk = &f4tt;
f4info_init (f4wk, nl, nr, popllist, poprlist, nl - 1);

vvd (mean, gtop, gbot, dim);
doranktest (mean, yvar, nl, nr, nl - 1, f4wk);
calcadm (totmean, f4wk->A, nl);


for (k = 0; k < nblocks; ++k) {
vvm (wtop, gtop, btop[k], dim);
vvm (wbot, gbot, bbot[k], dim);
vvd (mean, wtop, wbot, dim);
doranktest (mean, yvar, nl, nr, nl - 1, f4wk);
calcadm (tmean[k], f4wk->A, nl);
// printf("zzm: %d %9.3f\n", k, asum(tmean[k], nl)) ;
}

ZALLOC (wjack, nblocks, double);
for (k = 0; k < nblocks; ++k) {
wjack[k] = asum (bbot[k], dim);
 y = asum(tmean[k], nl) ; 
 if (isnan(y)) { 
  copyarr(totmean, tmean[k], nl) ; 
  wjack[k] = 0 ;
 }

}

wjackvest (jmean, jvar, nl, totmean, tmean, wjack, nblocks);
printf ("Jackknife mean:  ");
if (nl==1) jmean[0] = 1 ; 
printmatwl (jmean, 1, nl, nl);


if (doratio) {   
 estrat(jmean, jvar, nl, &yratmean, &yratsig, numboot) ;
}

free (wjack);
free (mean);
free (ymean);
free (wtop);
free (wbot);
free (totmean);
free2D (&tmean, nblocks);

return;


}

void
printss (double *mean, double *wt, int nl, int nr)
// experimental 
{
double *wa, *wb, *pt;
int a;
double y2;

ZALLOC (wa, nr, double);
ZALLOC (wb, nr, double);

for (a = 0; a < nl; ++a) {
pt = mean + a * nr;
vst (wb, pt, wt[a], nr);
vvp (wa, wa, wb, nr);
}

y2 = asum2 (wa, nr) / (double) nr;
vst (wb, wa, 1.0 / sqrt (y2), nr);

printf ("ssres:\n");
printmatwl (wa, 1, nr, nr);
printmatwl (wb, 1, nr, nr);
printf ("\n");

free (wa);
free (wb);

}

double sderiv1(double u, double theta, double *ff)  
// 2 * log likeighood
{
double t2, t3 ; 

t2 = theta*theta ; 
t3 = t2*theta ;
ff[0] = -log(theta) - u/theta ; 
ff[1] = -(1/theta) + (u/t2) ; 
ff[2] = (1/t2) - (2*u)/t3 ; 

return asum(ff, 3) ;

}



double sderiv(double *u, double *theta, int n, double *ff)  
{
double ww[3] ;
int k ; 

vzero(ff, 3) ; 

for (k=0; k<n; ++k) { 
sderiv1(u[k], theta[k], ww) ;   
vvp(ff, ff, ww, 3) ;
}


return asum(ff, 3) ;

}

double scorel(double *d, double *V, double *lam, int n) 
{

double *evecs, *evals, *uvals, *u1, *w1 ; 
int iter, k, a, b ;  
double *theta, base, step, *ff ; 
double  ymin, ymax, yinc, y, tmax, tmin, ynext ;  
double y1, y2, ylast ;

ZALLOC(evecs, n*n, double) ; 
ZALLOC(evals, n, double) ; 
ZALLOC(uvals, n, double) ;
ZALLOC(u1, n, double) ;
ZALLOC(theta, n, double) ;
ZALLOC(ff, n, double) ;

eigvecs(V, evals, evecs, n) ; 
mulmat(u1, evecs, d, n, n, 1) ; 

vst(u1, u1, 1000, n) ;
vst(evals, evals, 1000*1000, n) ;  
copyarr(evals, theta, n) ; 

//printf("unorms: %12.6f %12.6f\n", asum2(u1, n), 1000*1000*asum2(d, n)) ;
vvt(uvals, u1, u1, n) ;   
vmaxmin(theta, n, &tmax, &tmin) ;  
vmaxmin(evals, n, &ymax, &ymin) ;  
ymin =  -tmin ;
yinc = 0.001 ; 

printf("input: \n") ; 
printmatwl(uvals, 1, n, n) ;  
printmatwl(theta, 1, n, n) ;  
printnl() ;

y2 = 0.0 ; 
for (k=0; k<n; ++k) { 
printf("eig %3d ", k) ; 
printf("%12.6f ", u1[k]) ; 
y1 = sqrt(evals[k]) ; 
y = u1[k]/y1 ; 
printf("%12.6f %9.3f", y1, y) ; 
y2 += y*y ; 
printnl() ;
}
printf("chisq: %9.3f\n", y2) ;

base = 0 ; 

ylast = -1.0e20 ;
for (iter = 1; iter <= 100; ++iter) {  
vsp(theta, evals, base, n) ; 
sderiv(uvals, theta, n, ff) ;
if (ff[1] > 0) {  
ymin = base ; 
}
else  {  
ymax = base ; 
}

y = 0.5*(ymin+ymax) ; 
step = -ff[1]/ff[2] ; 
ynext = base + step ;  
printf("newton: %12.6f %12.6f %12.6f", ynext, y, step) ;
printf(" %12.6f %12.6f", ymin, ymax) ;
if (ynext<ymin)  ynext = y  ; 
if (ynext>ymax)  ynext = y ; 
if (ff[2]>0) ynext = y ; 
printf(" %12.6f ", ynext) ;
printnl() ; 

printf("iter: %3d base: %15.9f val: %15.9f", iter, base, ff[0]) ;
printf(" minmax: %12.6f %12.6f", ymin, ymax) ;
printmat(ff, 1, 3) ;
base = ynext ;
y = ff[0] - ylast ; 
if (fabs(y) < 1.e-8) break ;  
ylast = ff[0] ;
}
vsp(theta, evals, base, n) ; 
sderiv(uvals, theta, n, ff) ;
*lam = base ;

y = ff[0] ; 



free(evecs) ; 
free(evals) ; 
free(uvals) ;
free(u1) ;
free(theta) ;
free(ff) ;

return y ;

}

double gendstat(double *ca, int b1, int b2, int nl, int nr, double *mean, double *var, int *ktable) 
// compute Z-score for generalized D-stat 
// sum_i ca[i] f_4(l0, i; b1, b2) 
// convention  b1 = -1 => r0 
{

 double *va, *vb, *vv, *ww ;
 int a, b, k, dim ; 
 double ym, yvar, ysig, y ;

 dim = nl*nr ;
 ZALLOC(va, nl, double) ;
 ZALLOC(vb, nr, double) ;

 ZALLOC(vv, dim, double) ; 

 copyarr(ca, va, nl) ; 

 if (b2<0) fatalx("bad qqsig") ; 
 vb[b2] = 1 ; 
 if (b1>=0) vb[b1] = -1 ; // b1 = -1 => base
 // now make tensor
 for (a=0; a<nl; ++a) { 
  for (b=0; b<nr; ++b) { 
   k = a*nr + b ;
   if (ktable != NULL) k = ktable[a*nr+b] ; 
   if (k<0)  fatalx("bad1\n") ;
   if (k>=dim) fatalx("bad2\n") ;
   vv[k] = va[a]*vb[b] ;  
/**
   if (vv[k] !=0)  {  
    printf("zzgd %d %d %9.3f ", a,  b, vv[k]) ; 
    printf(" %9.3f ", mean[k]) ;  
    printnl() ;
   }
*/
  }
 }
 ym = vdot(vv, mean, dim) ;
 yvar = scx(var, NULL, vv, dim) ;
 ysig = sqrt(yvar) ;
 y = ym/ysig ;  
 if (isnan(y)) { 
  printf("zzm %12.6f %12.6f %12.6f\n", ym, ysig, y) ;
  fatalx("(gendtat)\n") ; 
 }
 return y ; 

 free(va) ;
 free(vb) ;
 free(vv) ;  

}
int isnested(int a, int b) 
{
// YES if 1 in A => 1 in b 
 int x ; 
 x = ~b ; 
 x &= a ; 
 if (x==0) return YES ; 
 return NO ; 

}
int usage (char *prog, int exval)
{

  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");
  (void)fprintf(stderr, "   -V          ... toggle details mode ON.\n");

  exit(exval);
};

void load4(int *x, int a, int b, int c, int d)
{
 x[0] = a ;
 x[1] = b ;
 x[2] = c ;
 x[3] = d ;
}


void loadymv(double *ymean, double *yvar,  char *fstatsname, char **popllist, char **poprlist, int nleft, int nright) 
{

  char **eglist ; 
  int numeg, np, numfs, nh2 ; 
  double *ff3, *ff3var, *vest, *vvar ; 
  int a, b, c, d, i, j, t, tt, k, dim ; 
  int nl, nr ; 
  int **fsindex, *fs ; 
  double y, ysig, yv ; 

  ZALLOC(eglist, MAXPOPS, char *) ; 
  numeg = np = fstats2popl(fstatsname, eglist) ; 
  ZALLOC(ff3, np*np, double) ; 
  ZALLOC(ff3var, np*np*np*np, double) ; 
  loadfstats(fstatsname, ff3, ff3var, eglist, numeg) ; 

  nl = nleft -1 ; 
  nr = nright -1 ; 

  nh2 = numeg*(numeg-1) ; nh2 /= 2 ; 
  ZALLOC(vest, nh2, double) ; 
  ZALLOC(vvar, nh2*nh2, double) ;

  setvv(vest, vvar, ff3, ff3var, NULL, numeg) ; 

  t = nleft*nright ;  numfs = 0 ; 
 
  fsindex = initarray_2Dint(t, 4, -1) ; 
  a = indxstring(eglist, numeg, popllist[0]) ; 
  c = indxstring(eglist, numeg, poprlist[0]) ; 
  for (i=1; i<nleft; ++i) { 
   b = indxstring(eglist, numeg, popllist[i]) ; 
   if (b<0) fatalx("%s not found\n", popllist[i]) ;
   for (j=1; j<nright; ++j) { 
    d = indxstring(eglist, numeg, poprlist[j]) ; 
   if (d<0) fatalx("%s not found\n", poprlist[j]) ;
    fs = fsindex[numfs] ; 
    load4(fs, a, b, c, d) ; 
    ivmaxmin(fs, 4, NULL, &tt) ; 
    if (tt<0) fatalx("pop not found in fstats: %s %s %s %s\n", 
       popllist[0], popllist[i], poprlist[0], poprlist[j]) ; 
    ++numfs ; 

  }} 

  vv2ww(ymean, yvar, vest, vvar, numeg, fsindex, numfs) ;
/**
  printf("ymean stats %12.6f %12.6f\n", asum(ymean, numfs), asum2(ymean, numfs) ) ;
  printf("zza %15.9f %15.9f\n", vest[0], vvar[0]) ; 
  printf("zzb %15.9f %15.9f\n", vest[1], vvar[1]) ; 
  printf("zzc %15.9f %15.9f\n", ymean[0], yvar[0]) ; 
  printimat2D(fsindex, numfs, 4) ; 

  dim = nl*nr ; 
  for (i=0; i<numfs; ++i) { 
   y = ymean[i] ; 
   yv = yvar[i*numfs+i] ; 
   ysig = sqrt(yv) ; 
   printf("%12.6f ", y) ; 
   printf("%12.6f ", ysig) ; 
   printf("%9.3f ", y/ysig) ; 
   printnl() ; 
  } 

*/

  free(eglist) ; 
  free(ff3) ; 
  free(ff3var) ; 
  free(vest) ; 
  free(vvar) ; 

  free2Dint(&fsindex, t) ; 

}

int mkfstats(char *parname)  
{
 char sss[256] ; 
 char fsx[256] ; 
 char ppp[256] ;  
 char pops[256] ;  
 char tpar[256] ; 
 char fslog[256] ; 
 char **poplist, **popx ; 
 int numeg ; 

 int pid ; 
 int k ; 

 int retkode, trap ;

 FILE *fff ; 

 numeg = nleft + nright ; 
 ZALLOC(poplist, numeg+1, char *) ; 

 popx = poplist ; 

 if (basepop != NULL) {  
   poplist[0] = strdup(basepop) ;
   ++popx ; 
   ++numeg ;
 }
 
 copystrings(poprlist, popx, nright) ;
 popx += nright ; 
 copystrings(popllist, popx, nleft) ;

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
 freeup(poplist, nleft + nright + 1) ;

 if (keeptmp == YES) return 1  ; 

 remove(ppp) ;
 remove(pops) ;
 remove(tpar) ; 

 if (fstatslog == NULL) remove(fslog) ;
 return 1 ; 

}

double f4mean(int a, int b, int nl, int nr, double *mean, double *var, int *ktable)  

{

  int k ; 
  
   k = a*nr + b ;
   if (ktable != NULL) k = ktable[a*nr+b] ; 

   return mean[k] ;


}

double f4var(int a1, int a2, int b1, int b2, int nl, int nr, double *mean, double *var, int *ktable) 
{

 int k1, k2, dim = nl*nr ;

 k1 = a1*nr + b1 ;
 if (ktable != NULL) k1 = ktable[a1*nr+b1] ; 
 k2 = a2*nr + b2 ;
 if (ktable != NULL) k2 = ktable[a2*nr+b2] ; 
 return var[k1*dim + k2] ;

}

double worstb(double *bb, double *aa, int nl, int nr, double *mean, double *var, int *ktable)  
{

  int a1, a2, b1, b2 ; 
  double zl, zq, y, zans, ztop, zbot ;
  double *ll, *qq, *w1 ;
  int ret ; 

  ZALLOC(ll, nr, double) ;
  ZALLOC(qq, nr*nr, double) ;
  ZALLOC(w1, nr, double) ;

  for (b1=0; b1<nr; ++b1) { 
   zl = 0; 
   for (a1=0; a1<nl; ++a1) { 
     zl += aa[a1]*f4mean(a1, b1, nl, nr, mean, var, ktable) ; 
   }
   ll[b1] = zl ;  
   for (b2=0; b2<nr; ++b2) { 
    zq  = 0 ;
    for (a1=0; a1<nl; ++a1) { 
     for (a2=0; a2<nl; ++a2) { 
      zq += aa[a1]*aa[a2]*f4var(a1, a2, b1, b2, nl, nr, mean, var, ktable) ; 
    }}
    qq[b1*nr+b2] = zq ; 
  }} 
// defensive programming  
  y = trace(qq, nr) ; 
  vclear(w1, y*1.0e-20, nr) ; 
  adddiag(qq, w1, nr) ; 

// now solve equation 
  ret = solvit(qq, ll, nr, bb) ; 
  y = asum(bb, nr) ; 
  y = unclip(y, -1.0e-6, 1.0e-6) ; 
  vst(bb, bb, 1.0/y, nr) ; 
  ztop = vdot(bb, ll, nr) ; 
  zbot = scx(qq, NULL, bb, nr) ; 
  zans = ztop/sqrt(zbot) ;

/** 
  printf("zzworst\n") ;
  vst(qq, qq,  1000*1000, nr*nr) ; 
  vst(ll, ll,  1000*1000, nr) ; 
  printmatw(qq, nr, nr, nr) ; 
  printnl() ; 
  printmatw(ll, 1, nr, nr) ; 
  printnl() ; 

*/

  free(ll) ; 
  free(qq) ; 
  free(w1) ;

  return zans ; 

}

void
writemv(char *mvname, long *kmean, long *kvar, int n)
{
 FILE *fff ;

 if (mvname == NULL) return ; 
 openit(mvname, &fff, "w") ;

 printlmatwfile(kmean, 1, n, n, fff) ;
 printlmatwfile(kvar, n, n, n, fff) ;
 fclose(fff) ;
}

void setr1mat(int **r1mat, int al, int ar) 
{
  int a, b, x = 0, z ;
  int **xtab ;
  

  x = 0 ; 
  for (a = 0 ; a < ar; ++a) { 
   for (b = 0 ; b < al; ++b) { 
   x = b*ar + a ; 
   r1mat[a][x] = 1 ; 
   if (b==0) ++r1mat[a][x] ; 
  }} 


}

int 
setconstraint(double **x1mat, int **r1mat, double *wvec,  int al, int ar) 
{

  int fsize, t ; 
  int a, b, x = 0 ;
  int *piv ;
  static int ncall = 0 ; 
  
  ++ncall ;
  fsize = al*ar ;
  for (a = 0 ; a < ar; ++a) { 
   vzero(x1mat[a], fsize) ;
   b  = 0 ; 
   for (x = 0 ; x<fsize; ++x) { 
    t = r1mat[a][x] ; 
    if (t<=0) continue ; 
    x1mat[a][x] = wvec[b] ; 
    ++b ; 
   }
  }

  if (ncall == -1) { 
   printf("setc first call\n") ;
   printimat2D(r1mat, ar, fsize) ; 
   printmat(wvec, 1, al) ;
   printnl() ;
   printmat2D(x1mat, al, fsize) ; 
  }

  return 1 ; 




}

double tttphi(double *zmean, double *zvar, int aleft, int aright) 
{

 int d = aleft, fsize = aleft*aright ;
 int **r1mat ;
 double **x1mat ; 
 double *yy, *rhs, y, wz[2] ; 
 int k, numconstraint ; 


// now we sample coeffs using posterior 
  r1mat = initarray_2Dint(aright, fsize, -1) ;
  x1mat = initarray_2Ddouble(aright, fsize, -1) ;

  setr1mat(r1mat, aleft, aright) ;
  ZALLOC(yy, 101, double) ; 
  vclear(yy, -1.0e10, 101) ;

  ZALLOC(rhs, fsize, double) ; 
  for (k=1; k <= 100; ++k) { 
   if (d>2) break ;
   y = (double) k / 100.0 ; 
   wz[0] = y ; wz[1] = 1.0 - y ; 
   numconstraint = setconstraint(x1mat, r1mat, wz,  aleft, aright) ;
/**
   printimat2D(r1mat, d, fsize) ;
   printmat2D(x1mat, d, fsize) ;
*/
   vzero(rhs, fsize) ;
   yy[k]  =  phiint(x1mat, rhs, zmean, zvar, fsize, numconstraint) ; 
  }

  vmaxmin(yy, 100, &y, NULL) ;  
  vsp(yy, yy, -y, 100) ; 
  for (k=1; k<100; ++k) { 
   if (d>2) break ;
   y = (double) k / 100.0 ; 
   printf("zzphi: %9.3f %9.3f\n", y, yy[k]) ; 
  }

  return 0.0 ; 
  
}

// move to statsubs.c 
double phiint(double **x1mat, double *rhs, double *cmean, double *cvar, int fsize, int numconstraint)  
{
   double *ww, *qq, *ll, y, yldet ;  
   int pivot, a, k, l, t ;
   static int ncall = 0 ;  
   double *cmat, *cmatt, *zvar, *zmean ;    ;  
   int numc = numconstraint, n = fsize ;
   ++ncall ; 
   
   t = numc + n ; 
   ZALLOC(cmat, numc * fsize, double) ; 
   ZALLOC(cmatt, numc * fsize, double) ; 
   ZALLOC(zvar, numc * numc, double) ; 
   ZALLOC(zmean, numc, double) ; 
   ZALLOC(ww, t * t, double) ; 
   ZALLOC(qq, t * t, double) ; 
  
   for (k=0; k<numc; ++k) { 
    copyarr(x1mat[k], cmat + k*n, n) ;
   }
   transpose(cmatt, cmat, numc, n) ; 
   mulmat(ww, cmat, cvar, numc, n, n) ; 
   mulmat(zvar, ww, cmatt, numc, n, numc) ; 
   mulmat(zmean, cmat, cmean, numc, n, 1) ; 

   yldet = pdinv(qq, zvar, numc) ;
   y = scx(qq, zmean, rhs, numc) ;
   if (verbose || (ncall <= -1)) {
    printf("phiint: %9.3f %9.3f\n", yldet, y)  ;
    vst(zmean, zmean, 1000, numc) ;
    vst(zvar, zvar, 1000*1000, numc*numc) ;
    printmatw(cmean, 1, n, n) ; 
    printmatw(cvar, n, n, n)  ;
    printnl() ;
    printmatw(cmat, numc, n, n) ;
    printnl() ;
    printmat(zmean, 1, numc) ;
    printmat(zvar, numc, numc) ;
    fflush(stdout) ;
   }

   free(qq) ; 
   free(ww) ; 
   free(cmat) ; 
   free(cmatt) ; 
   free(zmean) ; 
   free(zvar) ; 

   return -0.5*(y + yldet) ; 

}

void mcest(double *qcoeffs, double *qvar,  
 double *coeffs, double *cvar, double *zmean, double *zvar, 
 int nl, int nr, int mctrials) 
{
    double *mz, *mzvar ; 
    double *mzrvec, *rhs, *ww ;
    double y, ylike, yy ;
    int **r1mat ; 
    int *pivarr ; 
    double **x1mat ;
    int fsize, numc, t ; 
    double *yaprob, *yphi, *yweight, *yyweight ;
    int iter ;  
    int numneg = 0, numtrial = 0 ; 


    fsize = nl*nr ;        

    r1mat = initarray_2Dint(nr, fsize, -1) ;
    x1mat = initarray_2Ddouble(nr, fsize, -1) ;
    setr1mat(r1mat, nl, nr) ;
 
    ZALLOC(rhs, fsize, double) ; 
    ZALLOC(mz, nl, double) ; 
    ZALLOC(mzvar, nl*nl, double) ; 
    ZALLOC(mzrvec, mctrials*nl, double) ; 

    ZALLOC(yaprob, mctrials, double) ;
    ZALLOC(yphi, mctrials, double) ;
    ZALLOC(yweight, mctrials, double) ;
    ZALLOC(yyweight, mctrials, double) ;



//  printf("zzmz:\n") ;
    copyarr(coeffs, mz, nl-1) ; 
    squish(mzvar, cvar, nl-1, nl, nl-1) ;
  
    if (verbose) {
     printf("zzmcest printcv:\n") ;
     printcv(coeffs, cvar, nl) ;
     printnl() ;  
     printcv(mz, mzvar, nl-1) ;
    }
    

   for (iter = 0; iter < mctrials; ++iter) { 
    verbose = NO ;
    ww = mzrvec + iter * nl ; 
    for (;;) { 
     ++numtrial ; 
     genmultgauss(ww, 1, nl-1, mzvar) ; 
     vvp(ww, ww, coeffs, nl-1) ;  // add mean 
     yy = 1.0 - asum(ww, nl-1) ; 
     ww[nl-1] = yy ; 
     vmaxmin(ww, nl,  NULL, &y) ; 
     if ((positive_coeffs == NO) || (y >= 0)) break ;
     ++numneg ;
    } 
    yaprob[iter] = ylike = -0.5*scx(mzvar, NULL, ww, nl-1) ;
    
   numc = setconstraint(x1mat, r1mat, ww,  nl, nr) ;
   vzero(rhs, fsize) ;
   yphi[iter] = y  = phiint(x1mat, rhs, zmean, zvar, fsize, numc) ;
  }
   vvm(yweight, yphi, yaprob, mctrials) ; 
   vmaxmin(yweight, mctrials, &y, NULL) ; 
   vsp(yweight, yweight, -y, mctrials) ; 
// printf("mcweights:\n") ;
// printmatl(yweight, 1, mctrials) ;
// now must exponentiate and balance;  
   copyarr(yweight, yyweight, mctrials) ;
   vexp(yweight, yweight, mctrials) ; 
   bal1(yweight, mctrials) ; 
   y = ess(yweight, mctrials) ;
   y /= (double) mctrials ;
   printf ("sampling efficiency:  %12.6f\n", y) ; 
   if (y<.001) printf("efficiency very low.  Bad fit?\n") ;
   vzero(qcoeffs, nl) ;
   vzero(qvar, nl*nl) ;
   for (iter = 0; iter < mctrials; ++iter) { 
    ww = mzrvec + iter * nl ;
    y = yweight[iter] ;
    vst(mz, ww, y, nl) ; 
    vvp(qcoeffs, qcoeffs, mz, nl) ;
   }
// qcoeffs is now mean
   for (iter = 0; iter < mctrials; ++iter) { 
    ww = mzrvec + iter * nl ;
    y = yweight[iter] ;
    vvm(mz, ww, qcoeffs, nl) ;
    vzero(mzvar, nl*nl) ;
    addouter(mzvar, mz, nl) ; 
    vst(mzvar, mzvar, y, nl*nl) ; 
    vvp(qvar, qvar, mzvar, nl*nl) ;
    if (verbose) {
     printf("zzw %d %9.3f ", iter, yyweight[iter]) ;
     printmat(ww, 1, nl) ;
    } 
 
   }

   if (positive_coeffs) { 
    printf("(mcest): number of negative coeff trials: %d\n", numneg) ;
   } 

   free(mz) ; 
   free(mzvar) ;
   free(mzrvec) ;
   free(rhs) ;
   free2Dint(&r1mat, nr) ;
   free2D(&x1mat, nr) ;
   free(yaprob) ;
   free(yphi) ;
   free(yweight) ;
   free(yyweight) ;

}

void printcv(double *coeffs, double *var, int n) 
{
   double *w1, *w2 ;

   ZALLOC(w1,  n, double) ;
   ZALLOC(w2,  n*n, double) ;

   vst(w1, coeffs, 1000, n) ;
   vst(w2, var, 1000*1000, n*n) ;
   printf("mean *1000:\n") ;
   printmatw(w1, 1, n, n) ;
   printnl() ;
   printf("covariance *1000*1000:\n") ; 
   printmatw(w2, n, n, n) ;
   printnl() ;
   getdiag(w1, w2, n) ; 
   vsqrt(w1, w1, n) ;
   printf("std. errors *1000:\n") ;  
   printmat(w1, 1, n) ;

   free(w1) ;
   free(w2) ;

}

