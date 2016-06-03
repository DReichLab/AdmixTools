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


#define WVERSION   "352"
// best analysis added
// hires added
// chrom: 23 added
// allsnps option added
// more hires output added

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

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
// if bankermode  bankers MUST be in quartet  at most one type 1 in quartet

int allsnps = NO;

double blgsize = 0.05;		// block size in Morgans */
double *chitot;
char *popleft, *popright;
char **popllist, **poprlist;

char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
int inbreed = NO;
double lambdascale;
double *jmean;
int *wkprint;

int ldregress = 0;
double ldlimit = 9999.0;	/* default is infinity */

char *outputname = NULL;
char *weightname = NULL;
FILE *ofile;

double **btop, **bbot, *gtop, *gbot;
int bnblocks, bdim = 0;
int details = NO;
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
void calcadm (double *ans, double *A, int n);
void printss (double *mean, double *wt, int nl, int nr);

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
  int a, b, x, t, col, wt;
  int *xind;
  int xmax, dof;

  int ***counts;
  int jjj, j1, j2;
  int nleft, nright;
  int *rlist, lbase, lpop, rbase;
  int *vfix;
  double *yscbest, *ychi;
  int *sbest;
  double *var;
  int *jvar;
  double *w0, *w1, *w2;
  double *wbest;
  int dim;


  F4INFO **f4info, *f4pt, *f4pt2, **g4info, *ggpt;


  readcommands (argc, argv);
  printf ("## qpAdm version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;

  setinbreed (inbreed);

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

  ZALLOC (ktable, nleft * nright + 100, int);

  d = nl * nr;
  dd = d * d;
  ZALLOC (ymean, d + 10, double);
  ZALLOC (ww, d + 10 + nl * nl, double);
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
		nrows, ncols, lbase, xind, nl, rbase, rlist, nr, nblocks);
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
  ZALLOC (wkprint, nl * nl, int);

  for (x = rank; x <= maxrank; ++x) {
    if (x == rank)
      printf ("codimension 1\n");
    if (x == maxrank)
      printf ("full rank\n");
    f4pt = f4info[x];
    printf4info (f4pt);
    printnl ();
  }

  f4pt = f4info[rank];
  calcadm (ww, f4pt->A, nl);
  printf ("best coefficients: ");
  if (hires)
    printmatl (ww, 1, nl);
  else
    printmat (ww, 1, nl);
  ZALLOC (wbest, nl, double);
  ZALLOC (w0, nl * nr, double);
  ZALLOC (w1, nl * nr, double);
  ZALLOC (w2, nl * nr, double);

  copyarr (ww, wbest, nl);
  ZALLOC (var, nl * nl, double);
  ZALLOC (jvar, nl * nl, int);

  printss (ymean, ww, nl, nr);

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
  vst (ww, var, 1.0e6, nl * nl);
  fixit (jvar, ww, nl * nl);
  printf ("error covariance (* 1000000)\n");
  printimatl (jvar, nl, nl);

  if (hires) {
    printf ("hires: ");
    vst (ww, jmean, 1.0e6, nl);
    fixit (wkprint, ww, nl);
    printimatx (wkprint, 1, nl);
    printf ("  ");
    t = mktriang (ww, var, nl);
    vst (ww, ww, 1.0e8, t);
    fixit (wkprint, ww, t);
    printimatx (wkprint, 1, t);
    printnl ();
    printf ("[mean *1.0e6, var *1.0e8\n");
  }

  printnl ();
  printnl ();

//  we now loop over fix patterns.  
  xmax = pow (2, nl) - 1;
  ZALLOC (g4info, xmax + 2, F4INFO *);
  for (x = 0; x <= xmax; ++x) {
    ZALLOC (g4info[x], 1, F4INFO);
    f4info_init (g4info[x], nl, nr, popllist, poprlist, nl - 1);
  }

  ZALLOC (sbest, nl, int);
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
      calcadm (ww, f4pt->A, nl);
      vmaxmin (ww, nl, NULL, &y);
      dof = nnint (f4pt->dof);
      if (dof > 0)
	tail = rtlchsq (dof, f4pt->chisq);
      else
	tail = 0;
      if (y < -0.001)
	tail = 0;		// not feasible
      if (tail < 1.0e-30)
	tail = 0;
      if (tail > yscbest[wt]) {
	yscbest[wt] = tail;
	sbest[wt] = k;
	ychi[wt] = f4pt->chisq;
      }
      printf (" %12s  %1d   %3d %9.3f %15.6g ", binary_string (k, nl), wt,
	      dof, f4pt->chisq, tail);
      printmatx (ww, 1, nl);
      if (y < -.001)
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
    tail = rtlchsq (1, y);
    printf (" chi(nested): %9.3f p-value for nested model: %15.6g\n", y,
	    tail);
  }
  printnl ();
  fflush (stdout);
  if (details) {

    printf ("coeffs: ");
    printmat (wbest, 1, nl);
    dim = nl * nr;
    for (b = 0; b < nr; ++b) {
      vzero (w0, nl * nr);
      vzero (w1, nl * nr);
      for (a = 0; a < nl; ++a) {
	k = ktable[a * nr + b];
	printf ("details: ");
	printf ("%15s %15s ", popllist[a + 1], poprlist[b + 1]);
	y1 = w0[a] = ymean[k];
	w1[k] = wbest[a];
	printf ("%12.6f", y1);
	yy = yvar[k * dim + k];
	ysig = sqrt (yy);
	printf ("%12.6f", y1 / ysig);	// Z-score
	printnl ();
      }
      printf ("dscore: %15s ", poprlist[b + 1]);
// printf(" ") ; printmatl(w0, 1, nl) ;
      y1 = vdot (w0, wbest, nl);
// w0 is empirical f4;  wbest coeffs;  ideally should eb zero. 
      mulmat (w2, yvar, w1, dim, dim, 1);
      y2 = vdot (w1, w2, nl * nr);	// variance
      ysig = y1 / sqrt (y2);
      printf ("f4: %12.6f Z: %12.6f\n", y1, ysig);
      printnl ();
    }
  }


  printf ("## end of run\n");
  return 0;
}

void
calcadm (double *ans, double *A, int n)
{
  double *coeff, *rhs;

  if (A == NULL) {
    vzero (ans, n);
    return;
  }
  ZALLOC (coeff, n * n, double);
  ZALLOC (rhs, n, double);

  transpose (coeff, A, n, n - 1);
  vclear (coeff + n * (n - 1), 1.0, n);
  rhs[n - 1] = 1;

  linsolv (n, coeff, rhs, ans);



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

  while ((i = getopt (argc, argv, "p:vVd")) != -1) {

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
  getint (ph, "hires:", &hires);
  getint (ph, "allsnps:", &allsnps);
  getint (ph, "useallsnps:", &allsnps);
  getint (ph, "details:", &details);


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
  double *f4;
  int **xtop, *xt;
  int dim = nl * nr;
  int isok;

  if (nrows == 0)
    fatalx ("badbug\n");


  ZALLOC (f4, dim, double);
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
    }
  }

  totnum = 0;
  for (col = 0; col < ncols; ++col) {
    bnum = bcols[col];
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");
    vzero (f4, dim);
    isok = YES;
    if (allsnps)
      isok = NO;
    for (a = 0; a < dim; ++a) {
      ret = getf4 (counts[col], xtop[a], &y);
      if ((allsnps == NO) && (ret < 0)) {
	isok = NO;
	break;
      }
      if ((allsnps == YES) && (ret < 0))
	continue;
      if (allsnps == YES)
	isok = YES;
      bbot[bnum][a] += 1;
      if (ret == 2)
	fatalx ("bad pop in numerator\n");
      f4[a] = y;
    }
    if (isok == NO)
      continue;
    vvp (btop[bnum], btop[bnum], f4, dim);
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
  free2Dint (&xtop, dim);
  free (wtop);
  free (wbot);

  for (k = 0; k < nblocks; ++k) {
    if (allsnps == NO)
      break;
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
  double y0, y1, ytot, ff[4];
  double h0, h1;

  *ans = 0.0;
  if (indx == NULL) {
    *ans = 1.0;
    return 2;
  }

  for (i = 0; i < 4; ++i) {
    a = indx[i];
    if (a < 0) {
      *ans = 1.0;
      return 2;
    }
    y0 = (double) xx[a][0];
    y1 = (double) xx[a][1];
    ytot = y0 + y1;
    if (ytot <= 0.0)
      return -1;
    ff[i] = y0 / ytot;
  }
  *ans = (ff[0] - ff[1]) * (ff[2] - ff[3]);
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

  if (dim <= 0)
    fatalx ("(calcevar) dimension unset\n");

  ZALLOC (mean, dim, double);
  ZALLOC (ymean, dim, double);
  ZALLOC (wtop, dim, double);
  ZALLOC (wbot, dim, double);

  ZALLOC (totmean, nl, double);
  tmean = initarray_2Ddouble (nblocks, nl, 0);

  f4wk = &f4tt;
  f4info_init (f4wk, nl, nr, popllist, poprlist, nl - 1);

  vvd (mean, gtop, gbot, dim);
  doranktest (mean, yvar, nl, nr, nl - 1, f4wk);
  calcadm (totmean, f4wk->A, nl);

//  printf("zzmean:  ") ;  printmatl(totmean, 1, nl) ;

  for (k = 0; k < nblocks; ++k) {
    vvm (wtop, gtop, btop[k], dim);
    vvm (wbot, gbot, bbot[k], dim);
    vvd (mean, wtop, wbot, dim);
    doranktest (mean, yvar, nl, nr, nl - 1, f4wk);
    calcadm (tmean[k], f4wk->A, nl);
  }

  ZALLOC (wjack, nblocks, double);
  for (k = 0; k < nblocks; ++k) {
    wjack[k] = asum (bbot[k], dim);
  }
  wjackvest (jmean, var, nl, totmean, tmean, wjack, nblocks);
  printf ("Jackknife mean:  ");
  printmatl (jmean, 1, nl);


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
