#include <stdio.h>
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


#define WVERSION   "1000" 
// best analysis added
// hires added
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

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

double yscale = .0001 ;

char *parname = NULL;
char *trashdir = "/var/tmp";
int colcaac = YES;
int hires = NO;
int maxrank, rank;
int nsnpused = 0;

Indiv **indivmarkers;
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

int allsnps = NO;
int hiprec_covar = NO ;

double blgsize = 0.05;		// block size in Morgans */
double *chitot;
char *popleft, *popright;
char **popllist, **poprlist;

char *instem = NULL ; 
char *indivname = NULL;
char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *badsnpname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
int inbreed = NO;
double lambdascale;
double *jmean;
long *wkprint;
int *zprint ; 

int ldregress = 0;
double ldlimit = 9999.0;	/* default is infinity */
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
calcevar (double *var, double *yvar, int nl, int nr,
	  double **btop, double **bbot, double *gtop, double *gbot,
	  int nblocks, int dim);

int getf4 (int **xx, int *indx, double *ans);
int getf4old (int **xx, int *indx, double *ans);
void calcadmfix (double *ans, double *A, int n, int *vf) ;
void calcadm (double *ans, double *A, int n);
void printss (double *mean, double *wt, int nl, int nr);
double scorel(double *d, double *V, double *lam, int n)  ;
double gendstat(double *ca, int b1, int b2, int nr, int nl, double *mean, double *var, int *ktable) ;
void addscaldiag(double *mat, double scal, int n) ; 
int isnested(int a, int b) ; 
int usage (char *prog, int exval); 

int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int **snppos;
  int *snpindx;
  char **eglist;
  int *nsamppops;
  int numeg;
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
  int nleft, nright;
  int *rlist, lbase, lpop, rbase;
  int *vfix;
  double *yscbest, *ychi, *yinfeasbest;
  int *sbest, *infeassbest, *nfeas, isfeas;
  double *var;
  long *jvar;
  double *w0, *w1, *w2;
  double *wbest;
  int dim, b1, b2;
  double pval ;
  double lfudge ; 
  int numchromp ;  
  int pos ;
  double ymem ; 


  F4INFO **f4info, *f4pt, *f4pt2, **g4info, *ggpt;


  lopos = -1 ; hipos = BIGINT - 1 ;  

  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;

  printf ("## qpAdm version: %s\n", WVERSION);
  if (seed==0) seed = seednum() ; 
  printf("seed: %d\n", seed) ;
  if (parname == NULL)
    return 0;

  numchromp = numchrom+ 1 ;

  setinbreed (inbreed);
  setfancyf4(fancyf4) ; 

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 

  if (outputname != NULL)
    openit (outputname, &ofile, "w");

  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);
  setindm (indivmarkers);

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

  nleft = numlines (popleft);
  ZALLOC (popllist, nleft, char *);
  nright = numlines (popright);
  ZALLOC (poprlist, nright, char *);
  nleft = loadlist (popllist, popleft);
  nright = loadlist (poprlist, popright);

  printnl ();
  printf ("left pops:\n");
  for (k = 0; k < nleft; ++k) {
    printf ("%s\n", popllist[k]);
  }

  printnl ();
  printf ("right pops:\n");
  for (k = 0; k < nright; ++k) {
    printf ("%s\n", poprlist[k]);
  }

  printnl ();


  if (nleft > nright)
    fatalx ("nleft not less than or equal to nright\n");

  for (k = 0; k < nright; ++k) {
    j = indxindex (popllist, nleft, poprlist[k]);
    if (j >= 0)
      fatalx ("population in both left and right: %s\n", poprlist[k]);
  }

  numeg = npops = nleft + nright;
  ZALLOC (eglist, npops, char *);
  copystrings (popllist, eglist, nleft);
  copystrings (poprlist, eglist + nleft, nright);

  ckdup (eglist, numeg);

  ZALLOC (nsamppops, numeg, int);

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
    printf ("%3d %20s %4d\n", i, eglist[i], t);
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

  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);

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

  lbase = 0;
  rbase = nleft;
  ZALLOC (rlist, nright - 1, int);
  for (a = 1; a < nright; ++a) {
    rlist[a - 1] = nleft + a;
  }
  nr = nright - 1;
  nl = nleft - 1;

  ZALLOC (ktable, nleft * nright + 100, int);

  d = nl * nr;
  dd = d * d;
  ZALLOC (ymean, d + 10, double);
  ZALLOC (ww, d + 10 + nl * nl, double);
  ZALLOC (lambda,  nl, double);
  ZALLOC (yvar, dd, double);

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


  y = doq4vecb (ymean, yvar, counts, bcols,
		nrows, ncols, lbase, xind, nl, rbase, rlist, nr, xnblocks);
// y is jackknife dof 
  printf ("dof (jackknife): %9.3f\n", y);
  printf ("numsnps used: %d\n", nsnpused);

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

//  printss (ymean, ww, nl, nr);

  calcevar (var, yvar, nl, nr, btop, bbot, gtop, gbot, bnblocks, bdim);
  getdiag (ww, var, nl);
  vsp (ww, ww, 1.0e-20, nl);
  vsqrt (ww, ww, nl);
  printf ("      std. errors: ");
  
  if (hires)
    printmatl (ww, 1, nl);
  else
    printmat (ww, 1, nl);



  printnl ();

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

  printf("summ: %s ", popllist[0]) ; 
  printf(" %3d ", nl) ; 
  printf(" %12.6f ", pval) ; 
  printmatwx(jmean, 1, nl, nl) ; 
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
  }

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

  getstring (ph, "output:", &outputname);
  getstring (ph, "outputname:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "blockname:", &blockname);

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
  getint (ph, "fancyf4:", &fancyf4);
  getint (ph, "useallsnps:", &allsnps);
  getint (ph, "details:", &details);
  getint (ph, "seed:", &seed) ; 
  getint (ph, "hiprec_covar:", &hiprec_covar) ; 
  getdbl (ph, "diagplus:", &yscale);



printf ("### THE INPUT PARAMETERS\n");
printf ("##PARAMETER NAME: VALUE\n");
writepars (ph);

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
ktable[a * nr + b] = k;
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

 addscaldiag(yvar, yscale, dim) ;

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

void
calcevar (double *var, double *yvar, int nl, int nr, double **btop,
  double **bbot, double *gtop, double *gbot, int nblocks, int dim)
{
double *ymean, *mean, *wtop, *wbot, *totmean, **tmean, *wjack;
F4INFO *f4wk, f4tt;
int k;
double y ; 

if (dim <= 0)
fatalx ("(calcevar) dimension unset\n");

ZALLOC (mean, dim, double);
ZALLOC (ymean, dim, double);
ZALLOC (wtop, dim, double);
ZALLOC (wbot, dim, double);

ZALLOC (totmean, nl, double);
tmean = initarray_2Ddouble (nblocks, nl, 0);

addscaldiag(yvar, yscale, dim) ;

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

wjackvest (jmean, var, nl, totmean, tmean, wjack, nblocks);
printf ("Jackknife mean:  ");
if (nl==1) jmean[0] = 1 ; 
printmatwl (jmean, 1, nl, nl);


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
 double ym, yvar, ysig ;

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
// printf("zzm %12.6f %12.6f %12.6f\n", ym, ysig, ym/ysig) ;
 return ym/ysig ; 

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
