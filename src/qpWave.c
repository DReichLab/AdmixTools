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

/** 
 // like qp4wave but all left pops analyzed together
 // chrom: 23 now supported
 // bugfix in doq4vecb (same as qpAdm)
 // instem:   
 // cleaned up boundary case (m=n in ranktest) 
 // fstats properly supported
*/


#define WVERSION   "1520" 

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

double yscale = 0.0001 ; 

char *parname = NULL;
char *trashdir = "/var/tmp";
int qtmode = NO;
int colcaac = YES;
int hires = NO;
int maxrank = 99;
int nsnpused = 0;

Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;
int seed = 0;
int chisqmode = NO;		// approx p-value better to use F-stat
int dotpopsmode = YES;
int noxdata = YES;		/* default as pop structure dubious if Males and females */
int popsizelimit = -1;
int gfromp = NO;		// genetic distance from physical 
int xchrom = -1;
// if bankermode  bankers MUST be in quartet  at most one type 1 in quartet

int allsnps = -99 ;
int oldallsnpsmode = NO ;
int deletefstats = NO ; 
char *fstatslog = NULL ; 
char *tmpdir = NULL ;
int keeptmp = NO ;

int fancyf4 = YES;

double blgsize = 0.05;		// block size in Morgans */
double *chitot;
char *popleft, *popright;
char **popllist, **poprlist;
int nleft, nright;

char *instem = NULL ; 
char *indivname = NULL;
char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *badsnpname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
char *blockname = NULL;
int inbreed = -99;
double lambdascale;

int numboot = 1000 ; 

char *fstatsname = NULL ; 
char *fstatsoutname = NULL ; 

char *outputname = NULL;
char *weightname = NULL;
FILE *ofile;

double **btop, **bbot, *gtop, *gbot;
int bnblocks, bdim = 0;


void readcommands (int argc, char **argv);

double
doq4vecb (double *ymean, double *yvar, int ***counts, int *bcols,
	  int nrows, int ncols, int lbase, int *llist, int nl, int rbase,
	  int *rlist, int nr, int nblocks);

int getf4 (int **xx, int *indx, double *ans);
int usage (char *prog, int exval);

void loadymv(double *ymean, double *yvar,  char *fstatsname, char **popllist, char **poprlist, int nleft, int nright) ;
int  mkfstats(char *parname)   ; 

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
  double zscore, y1, y2, y, ysig, tail, yy1, yy2, yy;
  int *blstart, *blsize, nblocks;
  int xnblocks;			/* for xsnplist */
  int *bcols;
  double maxgendis;
  int **xtop;
  int npops = 0;
  int nr, nl, minnm, d, dd;
  double *ymean, *yvar, *yvarinv, *ww;
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
  int a, b, x, t, col;
  double T2, dof;
  int *xind;
  int numchromp ; 
  int retkode ;


  int ***counts, **ccc;
  int jjj, j1, j2;
  int *rlist, lbase, lpop, rbase;

  F4INFO **f4info, *f4pt, *f4pt2;
 
  double ymem ; 


  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;

  printf ("## qpWave version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;

  if (tmpdir != NULL) setenv("STMP", tmpdir, 1) ;

  numchromp = numchrom + 1 ; 

  setinbreed (inbreed);

  nleft = numlines (popleft);
  ZALLOC (popllist, nleft, char *);
  nright = numlines (popright);
  ZALLOC (poprlist, nright, char *);
  nleft = loadlist (popllist, popleft);
  nright = loadlist (poprlist, popright);

  if (nleft==0) fatalx("no pops in left!\n") ;
  if (nright==0) fatalx("no pops in right!\n") ;

  if (allsnps == -99) { 
   allsnps = NO ; 
   if (fstatsname == NULL) { 
    printf("allsnps set NO.  It is recommended that allsnps be set explicitly\n") ;
   }
  }
  if (inbreed == -99) { 
   inbreed = allsnps ; 
   printf(" *** recommended that inbreed be explicitly set ***\n") ;
  }

  setinbreed(inbreed) ;


  if ((allsnps == YES) && (oldallsnpsmode == NO) && (fstatsname == NULL)) {
   if (mkfstats(parname) < 0) { 
    printf("qpfstats failure ... terminating\n") ;
    return -1 ; 
   }; 
  }


  if ((allsnps == YES) && (oldallsnpsmode == YES)) {
   printf("oldallsnpsmode deprecated!\n") ;
  }

  if (outputname != NULL)
    openit (outputname, &ofile, "w");

  if (fstatsname == NULL) { 

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 


  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);
  setindm (indivmarkers);
  setfancyf4(fancyf4) ; 

  k = getgenos (genotypename, snpmarkers, indivmarkers,
		numsnps, numindivs, nignore);

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    chrom = cupt->chrom;
    if ((xchrom > 0) && (chrom != xchrom))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > numchromp)
      cupt->ignore = YES;
    if ((chrom == numchromp) && (xchrom != numchromp))
      cupt->ignore = YES;
  }
 }

  if (popleft == NULL) fatalx("no popleft\n") ;
  if (popright == NULL) fatalx("no popright\n") ;


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

  printf ("number of blocks for block jackknife: %d\n", nblocks);
  xnblocks = nblocks += 10 ; 

  ZALLOC (counts, ncols, int **);
  for (k = 0; k < ncols; ++k) {
    counts[k] = initarray_2Dint (numeg, 2, 0);
  }

  countpops (counts, xsnplist, xindex, xtypes, nrows, ncols);


  ZALLOC (bcols, ncols, int);	// blocks for columns -1 => unused
  ivclear (bcols, -1, ncols);
  for (k = 0; k < ncols; k++) {
    cupt = xsnplist[k];
    bcols[k] = cupt->tagnumber;
  }

 }

  lbase = 0;
  rbase = nleft;
  ZALLOC (rlist, nright - 1, int);
  for (a = 1; a < nright; ++a) {
    rlist[a - 1] = nleft + a;
  }
  nr = nright - 1;
  nl = nleft - 1;

  d = nl * nr;
  dd = d * d;
  ZALLOC (ymean, d + 10, double);
  ZALLOC (ww, d + 10, double);
  ZALLOC (yvar, dd, double);

  ZALLOC (xind, nl, int);
  for (jjj = 1; jjj < nleft; ++jjj) {
    xind[jjj - 1] = jjj;
  }

/**
 mnmin = MIN(nl, nr) ;
 t = MIN(mnmin, 4) ;
*/

  if (maxrank < 0)
    maxrank = 4;

  minnm = MIN (nl, nr);
  maxrank = MIN (minnm, maxrank);

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
 if (retkode == -2) printf("f4 variance absursly small.  Aborting run\n") ;

 if (retkode < 0) return -1 ;


  for (x = 0; x <= maxrank; ++x) {
    f4pt = f4info[x];
    f4pt->dofjack = y;
    doranktest (ymean, yvar, nl, nr, x, f4pt);
    if (x > 0) {
      f4pt2 = f4info[x - 1];
      f4pt->dofdiff = f4pt2->dof - f4pt->dof;
      f4pt->chisqdiff = f4pt2->chisq - f4pt->chisq;
    }
  }


  for (x = 0; x <= maxrank; ++x) {
    f4pt = f4info[x];
    printf4info (f4pt);
  }

  if (deletefstats) { 
// fstatsname has been created  
   printf("removing %s\n", fstatsname) ;
   remove(fstatsname) ; 
  }

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qpWave: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;
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

  while ((i = getopt (argc, argv, "p:vVh")) != -1) {

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
  getdbl (ph, "blgsize:", &blgsize);
  getint (ph, "numchrom:", &numchrom);

  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "chrom:", &xchrom);
  getint (ph, "maxrank:", &maxrank);
  getint (ph, "allsnps:", &allsnps);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "fancyf4:", &fancyf4);
  getint (ph, "keeptmp:", &keeptmp);
  getdbl (ph, "diagplus:", &yscale);
  getstring (ph, "blockname:", &blockname);
  getstring (ph, "fstatsname:", &fstatsname);
  getstring (ph, "fstatsoutname:", &fstatsoutname);
  getstring (ph, "fstatslog:", &fstatslog);
  getstring (ph, "tmpdir:", &tmpdir);

  getint (ph, "numboot:", &numboot) ; 


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
double ya, yb, y;
double *top, *bot, *wjack, *wbot, *wtop, **bjtop;
double *bwt;
double ytop, ybot;
double y1, y2, yscal;
double *w1, *w2, *ww, m1, m2;
double *mean;

int bnum, totnum;
int ret;
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
fatalx ("lobotgic bug\n");

vzero (f4, dim);
vzero (f4bot, dim);

isok = YES;
for (a = 0; a < dim; ++a) {
ret = getf4 (counts[col], xtop[a], &y);
if (ret==2) {
 printf("bad quad: "); printimat(xtop[1], 1, 4) ;
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
int usage (char *prog, int exval)
{

  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");

  exit(exval);
};

// should find a way to share code here with qpAdm
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
 char **poplist ; 
 int numeg ; 

 int pid ; 
 int k ; 

 int retkode, trap ;

 FILE *fff ; 

 numeg = nleft + nright ; 
 ZALLOC(poplist, numeg, char *) ; 
 
 copystrings(poprlist, poplist, nright) ;
 copystrings(popllist, poplist+nright, nleft) ;

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
 sprintf(sss, "fgrep -v fstatsoutname: %s | fgrep -v poplistname: > %s", parname, tpar) ; 
 system(sss) ; 
 copyfs(tpar, fff) ; 
 fprintf(fff, "fstatsoutname: %s\n", fsx) ; 
 writestrings(pops, poplist, numeg) ;
 
 
// printf("cmd1: %s\n", sss) ; 
 fprintf(fff, "poplistname: %s\n", pops) ; 
 fclose(fff) ; 
/**
 sprintf(sss, "cat %s\n", ppp) ; 
 system(sss) ; 
*/
 sprintf(sss, "qpfstats -p %s > %s", ppp, fslog) ; 
 retkode =  system(sss) ; 
 trap = openit_trap(fsx, &fff, "r") ;
 fclose(fff) ; 
 if (trap == NO) retkode = -88 ;
 
 // printf("exiting mkfstats\n") ; 
 if ((fstatsoutname != NULL) && (retkode >= 0)) { 
  sprintf(sss, "cp %s %s", fsx, fstatsoutname) ; 
  system(sss) ;
  printf("fstats file: %s made\n", fstatsoutname) ;
 } 
 if (retkode < 0) return -1 ;
 freeup(poplist, numeg) ;

 if (keeptmp == YES) return 1  ; 

 remove(ppp) ;
 remove(pops) ;
 remove(tpar) ; 

 if (fstatslog == NULL) remove(fslog) ;
 return 1 ; 

}
