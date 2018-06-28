#include <stdio.h>

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

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
*/


#define WVERSION   "410"

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

double yscale = 0.0 ; 

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

int allsnps = NO;

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

int ldregress = 0;
double ldlimit = 9999.0;	/* default is infinity */

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


  int ***counts, **ccc;
  int jjj, j1, j2;
  int nleft, nright;
  int *rlist, lbase, lpop, rbase;

  F4INFO **f4info, *f4pt, *f4pt2;


  readcommands (argc, argv);
  printf ("## qpWave version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;

  setinbreed (inbreed);

  if (outputname != NULL)
    openit (outputname, &ofile, "w");

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 

  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);
  setindm (indivmarkers);

  k = getgenos (genotypename, snpmarkers, indivmarkers,
		numsnps, numindivs, nignore);

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    chrom = cupt->chrom;
    if ((xchrom > 0) && (chrom != xchrom))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > 23)
      cupt->ignore = YES;
    if ((chrom == 23) && (xchrom != 23))
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

  nblocks = numblocks (snpmarkers, numsnps, blgsize);
  nblocks += 10;
  ZALLOC (blstart, nblocks, int);
  ZALLOC (blsize, nblocks, int);

  ZALLOC (xsnplist, numsnps, SNP *);
  ZALLOC (tagnums, numsnps, int);

  if (popsizelimit > 0) {
    setplimit (indivmarkers, numindivs, eglist, numeg, popsizelimit);
  }

  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);

  printf ("snps: %d  indivs: %d\n", ncols, nrows);
  setblocks (blstart, blsize, &xnblocks, xsnplist, ncols, blgsize);
// loads tagnumber
  printf ("number of blocks for block jackknife: %d\n", xnblocks);

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


  y = doq4vecb (ymean, yvar, counts, bcols,
		nrows, ncols, lbase, xind, nl, rbase, rlist, nr, nblocks);
// y is jackknife dof 
  printf ("dof (jackknife): %9.3f\n", y);
  printf ("numsnps used: %d\n", nsnpused);

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

  printf ("## end of run\n");
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

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i) {

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

  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "chrom:", &xchrom);
  getint (ph, "maxrank:", &maxrank);
  getint (ph, "allsnps:", &allsnps);


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


double
hfix (int *x)
{
// correction factor counts in x
  double ya, yb, yt, h;
  ya = (double) x[0];
  yb = (double) x[1];
  yt = ya + yb;
  if (yt <= 1.5)
    fatalx ("(hfix)\n");
  h = ya * yb / (yt * (yt - 1.0));
  return h / yt;
}


int
getf4 (int **xx, int *indx, double *ans)
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


 *ans - 0 ;
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

