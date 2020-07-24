#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "regsubs.h"
#include "egsubs.h"
#include "qpsubs.h"
#include "eigsubs.h"
#include "globals.h"

#define Y  0
#define E  1
#define A  2

//  (YRI, CEU, Papua, .... )               


#define WVERSION   "200"

// useweight added  
// allsnps added
// doscale NO added
// small bug (error check in dofstats fixed 

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *rootname = NULL;
char *trashdir = "/var/tmp";
int details = NO;
int hires = NO;			/* %10.4f */
int qtmode = NO;
int inbreed = NO;
int allsnpsmode = NO;
char *f3name = NULL;
char *fstatsoutname = NULL;

Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;
int seed = 0;
int missingmode = NO;
int noxdata = YES;		/* default as pop structure dubious if Males and females */
int doanalysis = YES;		/* if no just print stats */
int nostatslim = 10;
int znval = -1;
int popsizelimit = -1;
int gfromp = NO;		// genetic distance from physical 

int forcezmode = NO;
double blgsize = 0.05;		// block size in Morgans */ double *chitot ;
double diag = 0.0;
int fstdmode = NO;		// YES denominators done as in qp3test
int xchrom = -1;
int *xpopsize;

int isinit = NO;
int lsqmode = NO;
double f2weight = 1.0;		// lsqmode only

char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *poplistname = NULL;
char *outliername = NULL;

char *dumpname = NULL;
char *dumpf3block = NULL;

char *outpop = NULL;
// list of outliers
char *basepop = NULL;
int basenum = -1;
double baseval = 0.0;

// outnum used for weights 
// basenum for f3 status  


int doscale = YES ; 
double lambdascale = -1.0;
int *f2ind, *ind2f;
double *vest, *vvar, *vvinv, *xvvar, *xvvinv;
double **vmix;
int *lmix, nmix;
int nh2, numeg;
int *ezero = NULL;
double wtmin = .0001;
double minvar = 0.0;		// minvalue for variance term
int quartet = NO;
int xnumeg;


char *outputname = NULL;
char *weightname = NULL;
FILE *ofile;
char **eglist;
char **egshort;
char **enames;
double zthresh = 3.0;
double f2diag = 0.0;
int useweights = YES;
double ymem ; 

void readcommands (int argc, char **argv);
void indiaestit (double *f2, double *f3, double *f4, int n);
void sol2 (double *co, double *rhs, double *ans);
void mkww (double *f3, double *f2, int k, double *ww, int n);
char *getshort (char *ss, int n);
int doff3 (double *ff3, double *ff3var, SNP ** xsnplist, int *xindex,
	   int *xtypes, int nrows, int ncols, int numeg, int nblocks,
	   double scale);
void map4x (double *aa, double *bb, int n2, int *indx);
void map4y (double *aa, double *bb, int n2, int *indx);
void getmv (int a, int b, int c, int d, double *mean, double *var,
	    double *xest, double *xvar);
void bumpm (double *y, int a, int c, double val, double *xest);
void bumpw (double *w, int a, int c, double val);
void loadppwts (double *ppwts, double *pwts, int n);
double calcxx (double *xxans, double *qmat, double *ppwts, double *rhs,
	       int nrow, int ncol);
void setwww (double **tmix, double *www, int n);
void getwww (double **tmix, double *www, int n);
double scorit (double *www, int n, double *pfix, double *ans);
void printvals (double **tmix, double *edgelen, int nedge);
void printfit (double *ww);
int iscanon (int a, int b, int c, int d);

void setrand (double *ww, int n);
void setsimp (double *ww, int n);
void forcez (double *cc, double *rr, int *ez, int n);
void calcfit (double *ff3fit, double *ww, int numeg);
void dump1 (FILE * dumpfile, double *ww, int n);
void loadpars (char *loadname, double *www, int nwts, double *xxans,
	       int nedge);
void read1 (FILE * loadfile, double *ww, int n);
void print4 (double *ff3, int a, int b, int c, int d, int numeg);
void printf3 (char *sss, FILE * fff, double *ww, char **eglist, int n);
int listsubset (int **x, int n, int k);
void balw (double **ww, int **vv, int n, int *nw);
void estff3 (double *fv, double *v, int numv, int *elist, int n);
void rcsquish (double *xmat, double *mat, int *cols, int oldn, int newn);
double ff4val (double *ff3, int a, int b, int c, int d, int numeg);
void dumpit (char *dumpname, double *ff3, double *ff3var, char **eglist,
	     int numeg);
void checkpd (double *a, int n);
void dumpf3 (char *dumpf3name, double **btop, double **bbot, int nblock);
int usage (char *prog, int exval);
void load4(int *x, int a, int b, int c, int d)  ;
void loadco(double *co, int *fs, int *fsb, int np) ;

int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int **snppos;
  int *snpindx;
  char **snpnamelist, **indnamelist;
  int i, j, k, k1, k2, k3, k4, kk ;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;
  double y1, y2, y, sig, tail, yy1, yy2;
  int *blstart, *blsize, nblocks;
  int xnblocks;			/* for xsnplist */
  int *bcols;
  double maxgendis;

  int t, num, n, n1, n2, n3, n4;

  int nindiv = 0, e, f, lag = 1;
  int nignore, numrisks = 1;
  double *xrow, *xpt;
  SNP **xsnplist;
  Indiv **xindlist;
  int *xindex, *xtypes;
  int nrows, ncols, nedge, m, nc;
  int chrom, numclear;
  int a, b, c, d, col, u, v, x;
  int *popsizes;
  double *qpscore;

  int popx[4];
  double *f2, *f2sig, *fst;
  double *ff3, *gg3, *ff3var, *ff3fit, *ff3sig;

  double *f3, *f3sig, *www, *ww2;
  int ng2, ng3, ng4;
  double *pwts;
  double *ppwts;
  double ytail;
  int *xpopsize;
  FILE *f3file = NULL ;
  double scale;

  int **findex, np, *fsnum, nfstats ; 
  int *tt, *jwork ; 
  int **basis, *basisnum,  *basisfn,  nbasis ; 
  double **fcoeffs ; 
  double *fsmean, *fssig ; 
  double *fbmean, *fbcovar ; 
  int kret, bad ;
  double *hrate, *hvalid ; 


  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;
  

  printf ("## qpfstats version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;
  if (xchrom == 23)
    noxdata = NO;

/**
  if (seed == 0)
    seed = seednum ();
  SRAND (seed);
  h:0
rintf ("seed: %d\n", seed);
*/

  setinbreed(inbreed) ;
  setallsnpsmode (allsnpsmode);

  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);
  k = getgenos (genotypename, snpmarkers, indivmarkers,
		numsnps, numindivs, nignore);

  if (poplistname != NULL) {
    ZALLOC (eglist, numindivs, char *);
    numeg = loadlist (eglist, poplistname);
    seteglist (indivmarkers, numindivs, poplistname);
  }

  else {
    setstatus (indivmarkers, numindivs, NULL);
    ZALLOC (eglist, MAXPOPS, char *);
    numeg = makeeglist (eglist, MAXPOPS, indivmarkers, numindivs);
  }
  setindm (indivmarkers);

  np = numeg ; 
  n = 0;  n2  = n3 = n4 = 0 ; 
  t = np*np*np*np ;  
  ZALLOC(findex, t, int *) ; 
  ZALLOC(fsnum, t, int) ; 
  ZALLOC(jwork, t, int) ; 
  ZALLOC(basisfn, np*np, int) ; 
  for (k=0; k<t; ++k) { 
   ZALLOC(findex[k], 4, int) ; 
   ivclear(findex[k], -1, 4) ;
  }
  for (a=0; a<np; ++a) { 
   for (b=a+1; b<np; ++b) { 
     load4(findex[n], a, b, a, b) ; 
     ++n ; ++n2 ; 
  }} 
  for (c=0; c<np; ++c) { 
   for (a=0; a<np; ++a) { 
     for (b=a+1; b<np; ++b) { 
      if (a==c) continue ; 
      if (b==c) continue ; 
      load4(findex[n], c, a, c, b) ; 
      ++n ; ++n3 ; 
   }}} 
   for (a=0; a<np; ++a) { 
     for (b=a+1; b<np; ++b) { 
      for (c=a+1; c<np; ++c) { 
       for (d=c+1; d<np; ++d) { 
       if (c==b) continue ; 
       if (d==b) continue ; 
        load4(findex[n], a, b, c, d) ; 
        ++n ; ++n4 ; 
   }}}} 
   for (t=0; t<n; ++t) { 
     fsnum[t] = kodeitb(findex[t], 4, np) ; 
   }

   nfstats = n ; 

// set up basis 
   t = np*np ; 
   x = 0 ; nbasis = 0 ; 
  basis = initarray_2Dint(t, 2, -1) ; 

  ZALLOC(basisnum, t, int) ;  
  ivclear(basisnum, -1, t) ;
  ZALLOC(ind2f, t, int) ;
  ivclear(ind2f, -1, t) ;
  basenum = 0 ; 
  
  
  for (a=1; a<np; ++a) { 
   for (b=a; b<np; ++b) { 
    basis[nbasis][0] = a ; 
    basis[nbasis][1] = b ; 
    basisnum[a*np+b] = nbasis ; 
    basisnum[b*np+a] = nbasis ; 
    ind2f[nbasis] = a*np + b ; 
    ivzero(jwork, 4) ; 
    jwork[1] = a ; 
    jwork[3] = b ;  
    kk = kodeitb(jwork, 4, np) ; 
    t = findfirst(fsnum, n, kk) ; 
    if (t<0) fatalx("logic bug\n") ;
    basisfn[nbasis] = t ; 
    ++nbasis ; 
  }} 
  printf("np: %d n: %d nbasis: %d coefsize: %d\n", np, n, nbasis, n*nbasis) ; 
  fcoeffs = initarray_2Ddouble(n, nbasis, 0) ;  

/**
  printimatw(basisnum, 1, np*np, np) ; 
  fflush(stdout) ; 
*/
                          
  for (t=0; t<n; ++t) { 
    loadco(fcoeffs[t], findex[t], basisnum, np) ;
  }
  

  ZALLOC (egshort, numeg, char *);
  for (i = 0; i < numeg; i++) {
    egshort[i] = strdup (getshort (eglist[i], 5));
    printf ("%3d %s\n", i, eglist[i]);
  }

  outnum = 0;
  basenum = 0;

//   outpop = strdup (eglist[0]);

  if (outpop == NULL) { 
   outpop = strdup("NULL") ; 
  }

  printf ("outpop: %s \n", outpop) ; 


  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    chrom = cupt->chrom;
    if ((xchrom > 0) && (chrom != xchrom))
      cupt->ignore = YES;
    if ((noxdata) && (chrom == 23))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > 23)
      cupt->ignore = YES;
  }

  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);
  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);
  ZALLOC (xtypes, nrows, int);

  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (eglist, numeg, indx->egroup);
    xtypes[i] = k;
    t = strcmp (indx->egroup, outpop);
    if (t == 0)
      xtypes[i] = outnum;
  }

  ZALLOC (xpopsize, numeg, int);
  for (i = 0; i < nrows; i++) {
    k = xtypes[i];
    ++xpopsize[k];
  }

  for (i = 0; i < numeg; i++) {
    printf ("population: %3d %20s %4d\n", i, eglist[i], xpopsize[i]);
  }
  for (i = 0; i < numeg; i++) {
    if (xpopsize[i] == 0)
      fatalx ("zero popsize\n");
  }


  printf ("before setwt numsnps: %d\n", numsnps);
  setwt (snpmarkers, numsnps, indivmarkers, nrows, xindex, xtypes, outpop,
	 eglist, numeg);

  numsnps = rmsnps (snpmarkers, numsnps, NULL);	//  rid ignorable snps
  printf ("setwt numsnps: %d\n", numsnps);
  if (numsnps == 0)
    fatalx ("no valid snps\n");

/**
  for (k = 0; k < numsnps; ++k) {
    if (useweights)
      break;
    cupt = snpmarkers[k];
    cupt->weight = 1;
  }
*/

  setmgpos (snpmarkers, numsnps, &maxgendis);
  if ((maxgendis < .00001) || (gfromp == YES))
    setgfromp (snpmarkers, numsnps);

  nblocks = numblocks (snpmarkers, numsnps, blgsize);
  xnblocks = nblocks + 10 ; 

  ZALLOC (blstart, nblocks, int);
  ZALLOC (blsize, nblocks, int);
  printf ("number of blocks for moving block jackknife: %d\n", nblocks);

  ZALLOC (xsnplist, numsnps, SNP *);

  if (popsizelimit > 0) {
    setplimit (indivmarkers, numindivs, eglist, numeg, popsizelimit);
  }

  ZALLOC (popsizes, numeg, int);
  cntpops (popsizes, indivmarkers, numindivs, eglist, numeg);
  setpopsizes (popsizes, eglist, numeg);

  for (k=0; k<numeg; ++k) { 
   t = popsizes[k] ; 
   if (t==0) fatalx("pop: %s has sample size 0\n", eglist[k]) ; 
   if (inbreed && (t==1)) fatalx("pop: %s has sample size 1 and inbreed set\n", eglist[k]) ; 
  }

  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);


  printf ("snps: %d  indivs: %d\n", ncols, nrows);
  setblocks (blstart, blsize, &xnblocks, xsnplist, ncols, blgsize);

  ng2 = numeg * numeg;
  ng3 = numeg * ng2;
  ng4 = numeg * ng3;

  ZALLOC (f2, ng2, double);
  ZALLOC (f2sig, ng2, double);
  ZALLOC (fst, ng2, double);
  ZALLOC (ff3, ng2, double);
  ZALLOC (ff3var, ng4, double);

  ZALLOC(hrate, numeg, double) ; 
  ZALLOC(hvalid, numeg, double) ; 
 
   calchet (hrate, hvalid, xsnplist, xindex, xtypes, nrows, ncols, numeg) ; 

   for (k=0; k<numeg; ++k) { 
    printf("pop: %15s ", eglist[k]) ;  
    printf(" hetrate: %9.3f", hrate[k]) ; 
    printf(" valid snps: %9.0f", hvalid[k]) ; 
    printf(" samples: %4d", popsizes[k]) ; 
    if (hvalid[k] < .01) printf(" *** ") ;
    printnl() ; 
   }


  scale =
    dofstnumx (fst, f2, f2sig, xsnplist, xindex, xtypes, nrows, ncols, numeg,
	       xnblocks, indivmarkers, NO);

  if (lambdascale <= 0.0)
    lambdascale = scale;

  if (doscale == NO) lambdascale = 1.0 ; 

  if (fabs(lambdascale - scale) > 1.0e-5) { 
    vst (f2, f2, lambdascale / scale, ng2);
    vst (f2sig, f2sig, lambdascale / scale, ng2);
  }

  if (isnan(lambdascale)) fatalx("scaling bug!\n") ;
  printf ("lambdascale: %9.3f\n", lambdascale);

  printf("statistics multiplied by 1000\n") ;

  printf ("fst:");
  printnl ();
  printmatz (fst, egshort, numeg);
  printnl ();
  printf ("f2:");
  printnl ();
  printmatz (f2, egshort, numeg);
  printnl ();
  printf3 ("f2:", f3file, f2, egshort, numeg);


   ZALLOC(fsmean, nfstats, double) ; 
   ZALLOC(fssig, nfstats, double) ; 

   ZALLOC(fbmean, nbasis, double) ; 
   ZALLOC(fbcovar, nbasis*nbasis, double) ; 

   if (basenum < 0) basenum = 0 ; 

   cputimes(0, 2) ; 
   kret = dofstats(fbmean, fbcovar, fcoeffs, nbasis, 
      fsmean, fssig, findex, nfstats, 
      xsnplist, xindex, xtypes, nrows, ncols, numeg, nblocks, lambdascale) ; 

   if (kret<0) { 
    bad = kret+1000*1000 ; 
    printf("*** warning ***\n") ; 
    printf("fstat with no data.  Unable to compute heterozygosity?: ") ; 
    printimat(findex[bad], 1, 4) ; 
    printf("bad quadruple: ") ;
    for (k=-0; k<4; ++k) { 
      a = findex[bad][k] ;  
      printf(" %s ", eglist[a]) ;
    }
    printnl() ;
   }

  y = cputimes(1, 2) ; 
  printf("time in dofstats: %9.3f\n", y) ; 


     for (i=0; i<nbasis; ++i) { 
      a = basis[i][0] ; 
      b = basis[i][1] ; 
      ff3[a*numeg+b] = ff3[b*numeg+a] = fbmean[i] ;  
       for (j=0; j<nbasis; ++j) { 
        c = basis[j][0] ; 
        d = basis[j][1] ; 
        y = fbcovar[i*nbasis+j] ; 
        set4x(ff3var, a, b, c, d, numeg, y) ; 
       }
     } 

     if (hires) dumpfstatshr(fstatsoutname, ff3, ff3var, eglist, numeg, ind2f, basenum) ; 
     else dumpfstats(fstatsoutname, ff3, ff3var, eglist, numeg, ind2f, basenum) ; 

   if (details) printf("raw fstat table (before fit): scale:  %9.3f\n", lambdascale) ;
   for (k=0; k<nfstats; ++k) { 

    if (details == NO) break ; 

    tt = findex[k] ; 

    a = tt[0] ; 
    b = tt[1] ; 
    c = tt[2] ; 
    d = tt[3] ; 

   printf("fstat: ") ; 
/**
   printf(" %s",  eglist[a]) ;
   printf(" %s",  eglist[b]) ;
   printf(" %s",  eglist[c]) ;
   printf(" %s",  eglist[d]) ;
*/
  
   printimatx(tt, 1, 4) ; 

   printf(" %12.6f", fsmean[k]) ;
   printf(" %12.6f", fssig[k]) ;
   printf(" %9.3f", fsmean[k]/fssig[k]) ;
   printnl() ; 

  }
  printf ("##end of qpfstats\n");


  for (i=0; i<nbasis; ++i) {
   a = basis[i][0] ; 
   b = basis[i][1] ; 
   t = basisfn[i] ; 
   printf("basis: %3d %3d ", a, b) ; 
   y = fbcovar[i*nbasis+i] ; y = sqrt(y) ; 
   printf("%12.6f %12.6f", fbmean[i], y) ; 
   printf(" :: ") ; 
   printf("%12.6f %12.6f", fsmean[t], fssig[t]) ; 
   printnl() ;
    
  }


  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qpfstats: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;

}


void
estff3 (double *fv, double *v, int numv, int *elist, int n)
// extract covarance for pops
{
  int a, b, xa, xb;

  for (a = 0; a < n; a++) {
    for (b = 0; b < n; b++) {
      xa = elist[a];
      xb = elist[b];
      fv[a * n + b] = v[xa * numv + xb];
    }
  }
}


void
balw (double **ww, int **vv, int n, int *nw)
// balance admixture weights
{
  int j, k, t;

  for (k = 0; k < n; ++k) {
    t = 0;
    for (j = 0; j < MAXW; ++j) {
      if (vv[k][j] < 0)
	break;
      ++t;
    }
    nw[k] = t;
    for (j = t; j < MAXW; ++j) {
      ww[k][j] = 0.0;
    }
    bal1 (ww[k], t);
  }
}

int
iscanon (int a, int b, int c, int d)
// test quartet for is it canonical?  
{
  if (a >= b)
    return NO;
  if (c >= d)
    return NO;
  if (a > c)
    return NO;
  if ((a == c) && (b > d))
    return NO;

  return YES;

}

void
print4 (double *ff3, int a, int b, int c, int d, int numeg)
{
  double y;

  y = ff4val (ff3, a, b, c, d, numeg);
  printf ("%4s ", get3 (egshort[a]));
  printf ("%4s ", get3 (egshort[b]));
  printf ("   ");
  printf ("%4s ", get3 (egshort[c]));
  printf ("%4s ", get3 (egshort[d]));
  printf ("%6d", nnint (1000.0 * y));
  printnl ();
}

double
ff4val (double *ff3, int a, int b, int c, int d, int numeg)
{
  double y1, y2, y3, y4;

  y1 = dump2 (ff3, a, c, numeg);
  y2 = dump2 (ff3, a, d, numeg);
  y3 = dump2 (ff3, b, c, numeg);
  y4 = dump2 (ff3, b, d, numeg);

  return y1 + y4 - (y2 + y3);

}

void
printvals (double **tmix, double *edgelen, int nedge)
{
  char sss[MAXSTR];
  int k;

  for (k = 0; k < nmix; ++k) {
    getmixstr (k, sss);
    printf ("%40s ", sss);
    printmat (tmix[k], 1, lmix[k]);
  }
  printnl ();

  for (k = 0; k < nedge; ++k) {
    printf ("%9s ", enames[k]);
  }
  printnl ();
  printmatw (edgelen, 1, nedge, nedge);


}

void
setwww (double **tmix, double *www, int n)
// copy tmix to vector 
{
  int k, l;
  double *ww;

  if (n != intsum (lmix, nmix))
    fatalx ("dimension bug\n");

  ww = www;
  for (k = 0; k < nmix; ++k) {
    l = lmix[k];
    copyarr (tmix[k], ww, l);
    bal1 (ww, l);
    ww += l;
  }

}

void
getwww (double **tmix, double *www, int n)
// copy vector to tmix  
{
  int k, l;
  double *ww;

  if (n != intsum (lmix, nmix))
    fatalx ("dimension bug\n");

  ww = www;
  for (k = 0; k < nmix; ++k) {
    l = lmix[k];
    copyarr (ww, tmix[k], l);
    ww += l;
  }
}

void
normvec (double *www, int n)
{
  double **tmix;

  tmix = initarray_2Ddouble (nmix, MAXW, 0.0);

  getwww (tmix, www, n);
  setwww (tmix, www, n);
  free2D (&tmix, nmix);
}

void
printf3 (char *sss, FILE * fff, double *ww, char **eglist, int n)
{
  int i, j, x;

  if (fff == NULL)
    return;
  fprintf (fff, "%s\n", sss);
  fprintf (fff, " %4s", "   ");
  for (i = 0; i < n; i++) {
    fprintf (fff, " %4s", get3 (eglist[i]));
  }
  fprintf (fff, "\n");
  for (i = 0; i < n; i++) {
    fprintf (fff, "%4s", get3 (eglist[i]));
    for (j = 0; j < n; j++) {
      x = nnint (1000 * ww[i * n + j]);
      fprintf (fff, " %4d", x);
    }
    fprintf (fff, "\n");
  }
  fprintf (fff, "\n");
}

void
indiaestit (double *f2, double *f3, double *f4, int n)
{

  double z1, z2, z3, z4, z5, z6;
  double u1a0, u2a0, v1m, v2m;
  double y1, xa02, xb, xmd;
  double u1, u2, v1, v2;
  double xu1, xu2, xa0, xa2, xm, xd;
  double xa1, xt;
  double r, s, x1, x2;
  double co[4], rr[2], ans[2];

  int a2, b, d;
  int c, c1, c2;

  a2 = E;
  b = Y;
  d = A;

  z1 = dump2 (f2, a2, b, n);
  z2 = dump2 (f2, a2, d, n);
  z3 = dump2 (f2, b, d, n);

  y1 = z2 - z3;			// a0 + a2 -  b  
  xa02 = 0.5 * (z1 + y1);
  xb = z1 - xa02;
  xmd = z3 - xb;

  printf ("a0 + a2: %9.3f\n", xa02);
  printf ("      b: %9.3f\n", xb);
  printf ("  m + d: %9.3f\n", xmd);

  for (c = d + 1; c < n; ++c) {
    z4 = dump3 (f3, a2, b, c, n);
    u1a0 = xa02 - z4;
    z5 = dump3 (f3, d, b, c, n);
    v1m = xmd - z5;
    printf ("%10s %10.4f %10.4f\n", eglist[c], u1a0, v1m);
  }

}

void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n, t;

  while ((i = getopt (argc, argv, "p:g:s:o:l:vVxh")) != -1) {

    switch (i) {

    case 'h':
      usage(basename(argv[0]), 0);
      break ; 

    case 'p':
      parname = strdup (optarg);
      break;

    case 's':
      seed = atoi (optarg);
      break;

    case 'l':
      lambdascale = atof (optarg);
      break;

    case 'v':
      printf ("version: %s\n", WVERSION);
      break;

    case 'V':
      verbose = YES;
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
  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);
  dostrsub (ph);

  getstring (ph, "genotypename:", &genotypename);
  getstring (ph, "snpname:", &snpname);
  getstring (ph, "indivname:", &indivname);
  getstring (ph, "outpop:", &outpop);
  getstring (ph, "output:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "poplistname:", &poplistname);
  getstring (ph, "fstatsoutname:", &fstatsoutname);

  getstring (ph, "f3log:", &f3name);
  getstring (ph, "root:", &rootname);
  getdbl (ph, "blgsize:", &blgsize);
  getdbl (ph, "diag:", &diag);
  getdbl (ph, "minvar:", &minvar);
  getdbl (ph, "lambdascale:", &lambdascale);
  getint (ph, "fstdmode:", &fstdmode);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "useallsnps:", &allsnpsmode);
  getint (ph, "allsnps:", &allsnpsmode);
  getint (ph, "hires:", &hires);

  getint (ph, "noxdata:", &noxdata);
  t = -1;
  getint (ph, "xdata:", &t);
  if (t >= 0)
    noxdata = 1 - t;
  getint (ph, "chrom:", &xchrom);
  getint (ph, "doanalysis:", &doanalysis);
  getint (ph, "quartet:", &quartet);

  getint (ph, "nostatslim:", &nostatslim);
  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "seed:", &seed);
  getint (ph, "details:", &details);
  getint (ph, "useweights:", &useweights);
  getint (ph, "scale:", &doscale);
  getint (ph, "doscale:", &doscale);
  

  printf ("### THE INPUT PARAMETERS\n");
  printf ("##PARAMETER NAME: VALUE\n");
  writepars (ph);

}

void
sol2 (double *co, double *rhs, double *ans)
// solve 2x2 system
{

  double c00, c01, c10, c11, r0, r1, ra, xa;

  c00 = co[0];
  c01 = co[1];
  c10 = co[2];
  c11 = co[3];
  r0 = rhs[0];
  r1 = rhs[1];

  ra = r0 * c11 - r1 * c01;
  xa = c00 * c11 - c10 * c01;
  ans[0] = ra / xa;
  if (c11 != 0.0) {
    ra = r1 - c10 * ans[0];
    xa = c11;
  }
  else {
    ra = r1 - c00 * ans[0];
    xa = c01;
  }
  ans[1] = ra / xa;

}


int
doff3 (double *ff3, double *ff3var, SNP ** xsnplist, int *xindex, int *xtypes,
       int nrows, int ncols, int numeg, int nblocks, double scale)
{

  int t1, t2;
  int a, b, c;
  int ng2, ng3;
  int c1[2], c2[2], *cc;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  SNP *cupt;
  double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop;
  double **btop, **bbot, wt;
  double *w1, *w2, *w3;
  double ytop, ybot;
  double y1, y2, yscal;
  int bnum;
  int numegm = numeg - 1;
  int u, v, x;
  double *estmat;
  int ngood = 0, kret;

  ng2 = numeg * numeg;
  ng3 = numeg * numeg * numeg;

  ZALLOC (w1, ng3, double);
  ZALLOC (w2, ng3, double);
  ZALLOC (w3, ng3 * numeg, double);
  ZALLOC (gtop, ng3, double);
  ZALLOC (gbot, ng3, double);
  ZALLOC (wtop, ng3, double);
  ZALLOC (wbot, ng3, double);
  ZALLOC (estmat, ng3, double);
  ZALLOC (djack, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  btop = initarray_2Ddouble (nblocks, ng3, 0.0);
  bbot = initarray_2Ddouble (nblocks, ng3, 0.0);

  ZALLOC (vest, nh2, double);
  ZALLOC (vvar, nh2 * nh2, double);

// printf("zz ") ;  printimat(ind2f, 1, nh2) ;

  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    if (cupt->ignore)
      continue;
    wt = cupt->weight;
    if (wt <= 0.0)
      continue;
    bnum = cupt->tagnumber;
    if (bnum < 0)
      continue;
    ++wjack[bnum];
    top = btop[bnum];
    bot = bbot[bnum];

    // oldf3yyx(estmat,  cupt, xindex, xtypes, nrows, numeg, indivmarkers) ;
    kret = f3yyx (estmat, cupt, xindex, xtypes, nrows, numeg, indivmarkers);
    if (kret < 0)
      continue;
    ++ngood;

    a = basenum;
    for (u = 0; u < nh2; u++) {
      x = ind2f[u];
      b = x / numeg;
      c = x % numeg;
      ytop = dump3 (estmat, a, b, c, numeg);
      if (fstdmode == NO) {


/**
        y1 = wt*ytop ; 
        if ((u==0) & (y1>1.01)) {  
         printf("zzu: %6d %9.3f %9.3f %9.3f\n", col, y1, wt, ytop) ;
        }
*/

	top[u] += wt * ytop;
	bot[u] += 1.0;
      }
      else {
	top[u] += ytop;
	bot[u] += 1.0 / wt;
      }
    }
  }
  dumpf3 (dumpf3block, btop, bbot, nblocks);

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    vvp (gtop, gtop, top, nh2);
    vvp (gbot, gbot, bot, nh2);
  }

  vsp (w2, gbot, 1.0e-10, nh2);
  vvd (w3, gtop, w2, nh2);

  vzero (ff3, numeg * numeg);
  for (u = 0; u < nh2; u++) {
    x = ind2f[u];
    b = x / numeg;
    c = x % numeg;
    y1 = w3[u];
    ff3[b * numeg + c] = ff3[c * numeg + b] = y1;
  }

  printf ("ff3 (unscaled):\n");
  printmatz (ff3, eglist, numeg);


  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    vvm (wtop, gtop, top, nh2);
    vvm (wbot, gbot, bot, nh2);
    vsp (wbot, wbot, 1.0e-10, nh2);
    vvd (top, wtop, wbot, nh2);	// delete-block estimate
  }

  vsp (gbot, gbot, 1.0e-10, nh2);
  vvd (gtop, gtop, gbot, nh2);


  wjackvest (vest, vvar, nh2, gtop, btop, wjack, nblocks);
  vst (vest, vest, scale, nh2);
  vst (vvar, vvar, scale * scale, nh2 * nh2);

/**
    vzero(w3, numeg*numeg) ;
    for (u=0; u<nh2 ; u++) { 
       x = ind2f[u] ; 
       b = x / numeg ;
       c = x % numeg ;
       y1 = vest[u] ;                      
       bump2(w3, b, c, numeg,y1)  ;
       if (b<c) bump2(w3, c, b, numeg,y1)  ;
    }
    printf("vest (unscaled):\n") ;
    printmatz(w3, eglist, numeg) ;

    copyarr(vvar, w3, nh2*nh2) ;

    choldc(w3, nh2, w2) ;
    printmat(w2, 1, nh2) ;
    map4x(vvar, w3, numeg, ind2f) ;
    for (b=0; b<numeg; ++b)  {  
     set4(w3, 0, b, 0, b, numeg, 1.0) ;
    }
    choldc(w3, ng2, w2) ;
    printmat(w2, 1, ng2) ;

    printnl() ;
    printnl() ;
*/

  vst (ff3, ff3, scale, numeg * numeg);	// correct ff3 not scaled 
  map4x (vvar, ff3var, numeg, ind2f);

//  vst(ff3var, ff3var, scale*scale, numeg*numeg*numeg*numeg) ;


  free (w1);
  free (w2);
  free (w3);

  free (gbot);
  free (wtop);
  free (wbot);
  free (estmat);
  free (djack);
  free (wjack);

  free2D (&btop, nblocks);
  free2D (&bbot, nblocks);

  return ngood;

}

void
map4y (double *aa, double *bb, int n2, int *indx)
// map 4d array (n1 x n1 x n1 x n1  -> b  n2 x n2 x n2 x n2 
{
  int u, v, a, b, c, d;
  int x;
  double y1;
  int nh2;

  nh2 = n2 * (n2 - 1);
  nh2 /= 2;
  vzero (aa, nh2 * nh2);
  for (u = 0; u < nh2; ++u) {
    for (v = u; v < nh2; ++v) {

      x = indx[u];
      a = x / n2;
      b = x % n2;
      x = indx[v];
      c = x / n2;
      d = x % n2;

      y1 = dump4 (bb, a, b, c, d, n2);
      aa[u * nh2 + v] = y1;
      aa[v * nh2 + u] = y1;
    }
  }
}

void
bumpm (double *y, int a, int c, double val, double *xest)
{

  int k;

  if (a == basenum)
    return;
  if (c == basenum)
    return;
  k = f2ind[a * numeg + c];
  *y += val * xest[k];

}

void
bumpw (double *w, int a, int c, double val)
{

  int k;

  if (a == basenum)
    return;
  if (c == basenum)
    return;
  k = f2ind[a * numeg + c];
  w[k] += val;

}

void
dumpf3 (char *dumpf3name, double **btop, double **bbot, int nblock)
{
#define SQ  40
  int fdes, ret;
  int k, ng2, ng4;
  char sss[SQ];

  if (dumpf3name == NULL)
    return;

  ridfile (dumpf3name);
  fdes = open (dumpf3name, O_CREAT | O_TRUNC | O_RDWR, 006);
  cclear ((unsigned char *) sss, CNULL, SQ);
  sprintf (sss, "numeg: %d nblock: %d nh2: %d", numeg, nblock, nh2);
  ret = write (fdes, sss, SQ * sizeof (char));
  if (ret < 0) {
    perror ("write failure");
    fatalx ("bad write:  %s", sss);
  }
  for (k = 0; k < numeg; ++k) {
    cclear ((unsigned char *) sss, CNULL, SQ);
    strncpy (sss, eglist[k], SQ);
    ret = write (fdes, sss, SQ * sizeof (char));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  %s", sss);
    }
  }
  for (k = 0; k < nblock; ++k) {
    ret = write (fdes, btop[k], nh2 * sizeof (double));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  block %d", k);
    }
    ret = write (fdes, bbot[k], nh2 * sizeof (double));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  block %d", k);
    }
  }
  close (fdes);
  chmod (dumpf3name, 0644);
}

void
dumpit (char *dumpname, double *ff3, double *ff3var, char **eglist, int numeg)
{
#define SQ  40
  int fdes, ret;
  int k, ng2, ng4;
  char sss[SQ];
  if (dumpname == NULL)
    return;

  ridfile (dumpname);
  fdes = open (dumpname, O_CREAT | O_TRUNC | O_RDWR, 006);
  cclear ((unsigned char *) sss, CNULL, SQ);
  sprintf (sss, "numeg: %d", numeg);
  ret = write (fdes, sss, SQ * sizeof (char));
  if (ret < 0) {
    perror ("write failure");
    fatalx ("bad write:  %s", sss);
  }
  for (k = 0; k < numeg; ++k) {
    cclear ((unsigned char *) sss, CNULL, SQ);
    strncpy (sss, eglist[k], SQ);
    ret = write (fdes, sss, SQ * sizeof (char));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  %s", sss);
    }
  }
  ng2 = numeg * numeg;
  ng4 = ng2 * ng2;
  ret = write (fdes, ff3, ng2 * sizeof (double));
  if (ret < 0) {
    perror ("write failure");
  }
  ret = write (fdes, ff3var, ng4 * sizeof (double));
  close (fdes);
  chmod (dumpname, 0644);

}

void
checkpd (double *a, int n)
{
  double *b, *d;
  ZALLOC (b, n * n, double);
  ZALLOC (d, n, double);
  vsp (b, a, 1.0e-20, n * n);
  vst (b, b, 1.0e10, n * n);
  getdiag (d, b, n);
  printmat (d, 1, n);
  choldc (b, n, b);
  printnl ();
  printnl ();
  getdiag (d, b, n);
  printmat (d, 1, n);
  free (b);
  free (d);
}

int usage (char *prog, int exval)
{

  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -b <val>    ... use <va> as base value.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -g <>   ... .\n");
  (void)fprintf(stderr, "   -s <val>   ... use <val> as seed.\n");
  (void)fprintf(stderr, "   -o <>   ... .\n");
  (void)fprintf(stderr, "   -l <val>    ... use <val> as lambda scale.\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");
  (void)fprintf(stderr, "   -x          ... toggle doAnalysis ON.\n");

  exit(exval);
};

void load4(int *x, int a, int b, int c, int d)  
{
 x[0] = a ; 
 x[1] = b ; 
 x[2] = c ; 
 x[3] = d ; 
}

void loadco(double *co, int *fs, int *fsb, int np) 

{
  int a, b, c, d, t  ;  

  a = fs[0] ; 
  b = fs[1] ; 
  c = fs[2] ; 
  d = fs[3] ; 
 
  t = fsb[a*np+c] ;  
  if (t>=0) co[t] += 1 ; 
  t = fsb[b*np+d] ;  
  if (t>=0) co[t] += 1 ; 
  t = fsb[a*np+d] ;  
  if (t>=0) co[t] -= 1 ; 
  t = fsb[b*np+c] ;  
  if (t>=0) co[t] -= 1 ; 


}
