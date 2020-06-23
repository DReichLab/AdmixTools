#include <stdio.h>

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "regsubs.h"
#include "egsubs.h"
#include "qpsubs.h"


#define WVERSION   "400"

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *trashdir = "/var/tmp";
//int verbose = NO ;
int qtmode = NO;
int colcalc = YES;
// int fstdetails = NO ;
int hires = NO;
//int fancynorm = NO ;
//int outnum = 0 ;

char *bwname = NULL;

Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;
int isinit = NO;
int markerscore = NO;
int seed = 0;
int chisqmode = NO;		// approx p-value better to use F-stat
int missingmode = NO;
//int plotmode = NO ;
int dotpopsmode = YES;
int noxdata = YES;		/* default as pop structure dubious if Males and females */
int pcorrmode = NO;
int pcpopsonly = YES;
int nostatslim = 10;
int znval = -1;
int popsizelimit = -1;
int gfromp = NO;		// genetic distance from physical 
int msjack = NO;
char *msjackname = NULL;
int msjackweight = YES;		// weighted jackknife
int bankermode = NO;
int forceclade = NO;
int numbanker = 0;
int xchrom = -1;
int xmode = NO;
// if bankermode  bankers MUST be in quartet  at most one type 1 in quartet

int jackweight = YES;
double jackquart = -1.0;

double plo = .001;
double phi = .999;
double pvhit = .001;
double pvjack = 1.0e-6;
double blgsize = 0.05;		// block size in Morgans */
double *chitot;
int *xpopsize;

char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
int inbreed = NO;
double bfudge = 1.0;
int *f2ind, *ind2f, ng2, nh2;
int optmode = YES;
double *zfmat, *zsampsize, zffval;
double *aweights = NULL;

int clinetest = NO ; 

char *outputname = NULL;
char *weightname = NULL;
FILE *ofile;

void readcommands (int argc, char **argv);
int readpopx (char *pname, char ***plists, int npops);
void doq4 (double *q4stat, double *q4statsig, int ***counts, int *bcols,
	   int nrows, int ncols, int *xtop, int numeg, int nblocks);
void doq4w (double *q4stat, double *q4statsig, int ***counts, int *bcols,
	    int nrows, int ncols, int *xtop, int numeg, int nblocks,
	    double *zfmat, double *zsampsize, double *aweights);

void doq4diff (double *q4stat, double *q4statsig, int ***counts, int *bcols,
	       int nrows, int ncols, int *xtop, int *xbot, int *xtop2,
	       int *xbot2, int numeg, int nblocks);

void
doq4var (double *q4stat, double *q4var, int ***counts, int *bcols,
	  int nrows, int ncols, int *xtop,  int *xtop2,
	  int numeg, int nblocks) ; 

int getf4 (int **xx, int *indx, double *ans);
int getf4ssize (int **xx, int *indx, double *pabba, double *pbaba,
		double *ssize);
void mkfmat (double *zmat, double *zampsize, int ***counts, int *bcols,
	     int nrows, int ncols, int *xtop, int numeg, int nblocks);
void update_b1 (double *zfmat, double *zsampsize, double zffval,
		double *aweight, double *bweight1, int nblocks);

int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int **snppos;
  int *snpindx;
  char **snpnamelist, **indnamelist;
  char **eglist;
  int *nsamppops;
  int lsnplist, lindlist, numeg;
  int i, j, k, k1, k2, k3, k4, kk;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;
  double y1, y2, y, ysig, tail, yy1, yy2;
  char ss[11];
  int *blstart, *blsize, nblocks;
  int xnblocks;			/* for xsnplist */
  int *bcols;
  double maxgendis;
  int xind[4];
  int xtop[4], xbot[4];
  int xtop2[4], xbot2[4];

  int ch1, ch2;
  int fmnum, lmnum;
  int num, n1, n2;

  int nindiv = 0, e, f, lag = 1;
  double xc[9], xd[4], xc2[9];
  int nignore, numrisks = 1;
  double *xrow, *xpt;
  SNP **xsnplist;
  int *tagnums;
  Indiv **xindlist;
  int *xindex, *xtypes;
  int nrows, ncols, m, nc;
  double zn, zvar;
  int weightmode = NO;
  double chisq, ynrows;
  int *numhits, t;
  int chrom, numclear;
  double gdis;
  double *pmean, *pnum, rscore[3], dstat[3], hscore[3], rrr[3], ww[4];
  int a, b, c, d, col;
  double *qpscores;
  double *hest, *hsig;
  double mingenpos, maxgenpos;
  int *qhit;			/* number of times pair is clade in quartet */
  int *qmiss;			/* number of times pair migration event implied */
  int **qplist, numqp = 0, maxqp = 10000;
  double *qpscore;


  double **dsctop, **dscbot;
  double **abx, **bax;
  int popx[4];
  double tn[4 * 5], td[4 * 4];
  double zzsig[5], zzest[5], zsc[5];
  double ymin;
  double *f2, *f2sig, *fst;

  double *f3, *f4, *f3sig, *f4sig;
  int t1, t2;
  int ***counts, **ccc;
  char ***plists;
  char *px;
  int nplist = 0, trun, nplisth;
  double *bsig, **xx;
  int naweights = 0;

  double **fmat, **sampsize;
  double *aweight, *bweight, *ffval;
  double q4stat[2], q4var[4] , w1[4], w2[4] ; 

  readcommands (argc, argv);
  printf ("## qpdslow version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;

  if (bwname != NULL) {
    t = numlines (bwname);
    ZALLOC (aweights, t, double);
    ZALLOC (bsig, t, double);
    xx = initarray_2Ddouble (2, t, 0);
    naweights = getxx (xx, t, 2, bwname);
    for (k = 0; k < naweights; ++k) {
      t = nnint (xx[0][k]);
      if (t != k)
	fatalx ("bad bw: %d %d\n", t, k);
      bsig[k] = xx[1][k];
    }

    copyarr (bsig, aweights, naweights);
  }

  nostatslim = MAX (nostatslim, 3);
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
    if (chrom >= 23)
      cupt->ignore = YES;
  }
// autosomes only

  nplist = numlines (popfilename);
  ZALLOC (plists, nplist, char **);
  num = readpopx (popfilename, plists, 4);
  nplist = num;
  printf ("nplist: %d\n", nplist);
  if (nplist == 0)
    return 1;


  ZALLOC (eglist, nplist * 4, char *);
  numeg = 0;

  for (trun = 0; trun < nplist; ++trun) {
    for (k = 0; k < 4; ++k) {
      px = plists[trun][k];
      t = strcmp (px, "NULL");
      if (t == 0)
	continue;
      t = indxindex (eglist, numeg, px);
      if (t < 0) {
	eglist[numeg] = strdup (px);
	++numeg;
      }
    }
  }


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
  ZALLOC (blstart, nblocks, int);
  ZALLOC (blsize, nblocks, int);
//  printf("number of blocks for moving block jackknife: %d\n", nblocks) ;

  ZALLOC (xsnplist, numsnps, SNP *);
  ZALLOC (tagnums, numsnps, int);

  if (popsizelimit > 0) {
    setplimit (indivmarkers, numindivs, eglist, numeg, popsizelimit);
  }

  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);

  printf ("snps: %d  indivs: %d\n", ncols, nrows);
  setblocks (blstart, blsize, &xnblocks, xsnplist, ncols, blgsize);
// loads tagnumber
  printf ("number of blocks for moving block jackknife: %d\n", xnblocks);

  ZALLOC (counts, ncols, int **);
  for (k = 0; k < ncols; ++k) {
    counts[k] = initarray_2Dint (numeg, 2, 0);
  }

  if (xmode && (xchrom == (numchrom + 1))) {
    countpopsx (counts, xsnplist, xindlist, xindex, xtypes, nrows, ncols);
    printf ("countpopsx called\n");
  }
  else {
    countpops (counts, xsnplist, xindex, xtypes, nrows, ncols);
    printf ("countpops called\n");
  }

  fmat = initarray_2Ddouble (nplist, nblocks, 0);
  sampsize = initarray_2Ddouble (nplist, nblocks, 0);

  ZALLOC (bcols, ncols, int);	// blocks for columns -1 => unused
  ivclear (bcols, -1, ncols);
  for (k = 0; k < ncols; k++) {
    cupt = xsnplist[k];
    bcols[k] = cupt->tagnumber;
  }

  for (a = 0; a < nplist; ++a) {
    for (i = 0; i < 4; ++i) {
      xtop[i] = indxindex (eglist, numeg, plists[a][i]);
    }

    mkfmat (fmat[a], sampsize[a], counts, bcols,
	    nrows, ncols, xtop, numeg, nblocks);
  }

  for (a = 0; a < nplist; ++a) {
    for (i = 0; i < 4; ++i) {
      xtop[i] = indxindex (eglist, numeg, plists[a][i]);
    }
    zfmat = fmat[a];
    zsampsize = sampsize[a];

    y = (double) nblocks;
//  bal1(aweights, nblocks) ;  vst(aweights, aweights, y, nblocks) ;
//  bal1(zsampsize, nblocks) ;  vst(zsampsize, zsampsize, y, nblocks) ;

/**
    for (k=0; k< naweights; ++k) { 
     printf("aweights: %d %9.3f %12.6f\n", k, aweights[k], zsampsize[k]) ;
    }
bwname:           aw.out
optmode:          YES
*/

    y1 = vdot (zfmat, zsampsize, nblocks);
    y2 = asum (zsampsize, nblocks);
    zffval = y1 / y2;
    printf ("D: %12.6f\n", zffval);

    if (aweights != NULL) {
      if (nblocks != naweights)
	fatalx ("block mismatch %s %d\n", nblocks, naweights);
      doq4w (&y, &ysig, counts, bcols,
	     nrows, ncols, xtop, numeg, nblocks, zfmat, zsampsize, aweights);
    }
    else {
      doq4 (&y, &ysig, counts, bcols, nrows, ncols, xtop, numeg, nblocks);
    }

    printf ("result: ");
    for (t = 0; t < 4; ++t) {
      printf ("%10s ", plists[a][t]);
    }
    printf ("%15.9f %12.6f  %9.3f", y, ysig, y / ysig);
    if (aweights != NULL)
      printf ("   fudge: %9.3f", bfudge);
    printnl ();
  }

  for (a=0; a< nplist ; ++a) { 
    if (clinetest == NO) break ; 
     for (b=a+1; b< nplist ; ++b) { 
      printf("clinetest:: ") ;  
      printstringsx(plists[a], 4) ;  
      printf(" :: ") ;
      printstringsx(plists[b], 4) ;  

     for (i = 0; i < 4; ++i) {
      xtop[i] = indxindex (eglist, numeg, plists[0][i]);
      xtop2[i] = indxindex (eglist, numeg, plists[1][i]);
     }
     doq4var (q4stat, q4var, counts, bcols,
	  nrows, ncols, xtop,  xtop2,  numeg,  nblocks) ; 

    vst(w1, q4stat, 1000, 2) ; 
    vst(w2, q4var, 1000*1000, 4) ; 
    printf("mean*1000:\n") ; printmatl(w1, 1, 2) ; 
    printf("covar*1M: ") ; printmatl(w2, 2, 2) ; 

    for (k=0; k<=100; ++k) { 
     y = ww[0] = (double) k / 100.0 ; 
     ww[1] = 1.0-y ; 
     y1 = vdot(q4stat, ww, 2) ; 
     y2 = scx(q4var, NULL, ww, 2) ;  y2 = sqrt(y2) ; 

     printf("mixtable: %9.3f ", y) ; 
     printf(" %12.6f %12.6f %9.3f\n", y1, y2, y1/y2) ;

    }
  }}  


  printf ("## end of qpdslow\n");
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

  while ((i = getopt (argc, argv, "f:p:vV")) != -1) {

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

    case 'f':
      bfudge = atof (optarg);
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
  getstring (ph, "popfilename:", &popfilename);
  getstring (ph, "output:", &outputname);
  getstring (ph, "outputname:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "bwname:", &bwname);
  getdbl (ph, "blgsize:", &blgsize);
  getdbl (ph, "bfudge:", &bfudge);

  getint (ph, "optmode:", &optmode);
  getint (ph, "nostatslim:", &nostatslim);
  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "chrom:", &xchrom);
  getint (ph, "xmode:", &xmode);
  getint (ph, "clinetest:", &clinetest) ; 

  printf ("### THE INPUT PARAMETERS\n");
  printf ("##PARAMETER NAME: VALUE\n");
  writepars (ph);

}

double
getrs (double *f4, double *f4sig, double *f2, double *f2sig,
       int a, int b, int c, int d, int numeg, double *rho)
{

  double y1, y2, ya, yb;

  y1 = dump4 (f4, a, b, c, d, numeg);
  y2 = dump4 (f4sig, a, b, c, d, numeg);
  ya = f2[a * numeg + b];
  yb = f2[c * numeg + d];

  *rho = y1 / sqrt (ya * yb);

  return y1 / y2;

}

int
islegal (int *xind, int *types, int mode, int n)
{
// at least 1 type 1 
// check all bankers set.  OR at least 2 bankers 
  int cc[3];
  int zz[MAXPOPS];
  int i, x, k, t;


  if (mode == NO)
    return YES;
  ivzero (zz, n);
  ivzero (cc, 3);
  for (i = 0; i < 4; i++) {
    t = xind[i];
    x = types[t];
    ++cc[x];
    ++zz[t];
  }

/**
  printf("zzq1 ") ;  printimat(xind, 1, 4) ;
  printimat(cc, 1, 4) ;
  printimat(zz, 1, n) ;
  printimat(types, 1, n) ;
*/

  if (cc[1] == 0)
    return NO;
  if (cc[2] >= 2)
    return YES;
//  if (cc[1] == 2) return NO ; 
  for (t = 0; t < n; ++t) {
    if (types[t] != 2)
      continue;
    if (zz[t] == 0)
      return NO;
  }
  return YES;
}

int
isclade (int *rr, int *zz)
// is a clade banker only? 
{

  int a, b;

  a = rr[0];
  b = rr[1];
  if ((zz[a] == 2) && (zz[b] == 2))
    return YES;

  a = rr[2];
  b = rr[3];
  if ((zz[a] == 2) && (zz[b] == 2))
    return YES;

  return NO;

}

void
setabx (double **abx, double **bax, int ***counts, int ncols, int numeg)
{


  int i, j, k, t1, t2, a, h;
  double y1, y2;
  int **ccc;

  clear2D (&abx, nh2, ncols, 0.0);
  clear2D (&bax, nh2, ncols, 0.0);

  for (k = 0; k < ncols; ++k) {

    ccc = counts[k];

    for (i = 0; i < numeg; i++) {
      for (j = i + 1; j < numeg; j++) {

	t1 = intsum (ccc[i], 2);
	t2 = intsum (ccc[j], 2);
	if (t1 == 0)
	  continue;
	if (t2 == 0)
	  continue;

	a = ccc[i][0];
	y1 = (double) a / (double) t1;
	a = ccc[j][0];
	y2 = (double) a / (double) t2;
	h = f2ind[i * numeg + j];

	abx[h][k] = y1 * (1 - y2);
	bax[h][k] = y2 * (1 - y1);

      }
    }
  }
}



int
readpopx (char *pname, char ***plists, int npops)
// reads lists of n pops on a line
{
  FILE *fff;
  char line[MAXSTR + 1], c;
  char *spt[MAXFF], *sx;
  char **pp;
  int nsplit, t, num = 0;

  openit (pname, &fff, "r");
  line[MAXSTR] = '\0';
  while (fgets (line, MAXSTR, fff) != NULL) {
    subcolon (line);
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < npops)
      fatalx ("length mismatch %s\n", line);
    ZALLOC (plists[num], npops + 1, char *);
    plists[num][npops] == NULL;
    pp = plists[num];
    for (t = 0; t < npops; ++t) {
      pp[t] = strdup (spt[t]);
    }
    ++num;
    freeup (spt, nsplit);
  }
  fclose (fff);
  return num;
}

void
q4jest (double *q4stat, double *q4statsig,
	double *fmat, double *blockvar, int nblocks)
{

  int a, b, c, d;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  double top, bot, *djack, *wjack, gtop, gbot, *wbot, *wtop;
  double *btop, *bbot, wt;
  double ytop, ybot, ytot, zabba, zbaba;
  double y1, y2, yscal;
  double *w1, *w2, *ww, m1, m2;
  int bnum, totnum;
  int ret;
  double t1 = 0;
  double t2 = 0, t;

  ZALLOC (w1, nblocks, double);
  ZALLOC (w2, nblocks, double);
  ZALLOC (ww, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  ZALLOC (djack, nblocks, double);
  ZALLOC (btop, nblocks, double);
  ZALLOC (bbot, nblocks, double);
  ZALLOC (wtop, nblocks, double);
  ZALLOC (wbot, nblocks, double);

  for (a = 0; a < nblocks; ++a) {
    t = blockvar[a];
    if (t <= 0)
      continue;
    btop[a] = fmat[a] / t;
    bbot[a] = 1 / t;
    wjack[a] = 1 / t;
  }

  gtop = asum (btop, nblocks);
  gbot = asum (bbot, nblocks);

  mean = gtop / gbot;
  printf ("mean: %12.6f\n", mean);

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    wtop[k] = gtop - top;
    wbot[k] = gbot - bot;
    wbot[k] += 1.0e-10;
    djack[k] = wtop[k] / wbot[k];	// delete-block estimate
  }

  weightjack (&jest, &jsig, mean, djack, wjack, nblocks);

  *q4stat = jest;
  *q4statsig = jsig;

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

void
doq4w (double *q4stat, double *q4statsig, int ***counts, int *bcols,
       int nrows, int ncols, int *xtop, int numeg, int nblocks, double *fmat,
       double *sampsize, double *wts)
{
  double y1;
  double *w1;
  int a;

  if (wts == NULL) {
    doq4 (q4stat, q4statsig, counts, bcols,
	  nrows, ncols, xtop, numeg, nblocks);
    return;
  }
  if (nrows == 0)
    fatalx ("badbug\n");

  ZALLOC (w1, nblocks, double);

  if (optmode) {
    y1 = bfudge;
    update_b1 (fmat, sampsize, zffval, aweights, &y1, nblocks);
    bfudge = y1;
    printf ("estimated bfudge:  %9.3f\n", bfudge);
  }

  for (a = 0; a < nblocks; ++a) {
    if (sampsize[a] == 0)
      continue;
    w1[a] = wts[a] + bfudge / sampsize[a];
  }
  q4jest (q4stat, q4statsig, zfmat, w1, nblocks);

  free (w1);
}

void
doq4 (double *q4stat, double *q4statsig, int ***counts, int *bcols,
      int nrows, int ncols, int *xtop, int numeg, int nblocks)
{

  int a, b, c, d;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  double top, bot, *djack, *wjack, gtop, gbot, *wbot, *wtop;
  double *btop, *bbot, wt;
  double ytop, ybot;
  double y1, y2, yscal;
  double *w1, *w2, *ww, m1, m2;
  int bnum, totnum;
  int ret;

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
  for (col = 0; col < ncols; ++col) {
    ret = getf4 (counts[col], xtop, &ytop);
    if (ret < 0)
      continue;
    if (ret == 2)
      fatalx ("bad pop\n");

    bnum = bcols[col];
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");

    btop[bnum] += ytop;
    bbot[bnum] += 1.0;
    ++wjack[bnum];
    ++totnum;

  }

  gtop = asum (btop, nblocks);
  gbot = asum (bbot, nblocks);

  mean = gtop / gbot;

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    wtop[k] = gtop - top;
    wbot[k] = gbot - bot;
    wbot[k] += 1.0e-10;
    djack[k] = wtop[k] / wbot[k];	// delete-block estimate
  }

  weightjack (&jest, &jsig, mean, djack, wjack, nblocks);

  *q4stat = jest;
  *q4statsig = jsig;

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

void
doq4var (double *q4stat, double *q4var, int ***counts, int *bcols,
	  int nrows, int ncols, int *xtop,  int *xtop2,
	  int numeg, int nblocks)
// calculates mean and covar of m1, m2  
{

  int a, b, c, d;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  double top, bot, *djack, *wjack, gtop, gbot, gtop2, gbot2, *wbot, *wtop;
  double *btop, *bbot, wt;
  double *btop2, *bbot2;
  double ytop, ybot, ytop2, ybot2;
  double y1, y2, yscal;
  double *w1, *w2, *ww, m1, m2, *ww1, *ww2, var[4], mm[2],  **jmean ;  
  int bnum, totnum, ret;

  if (nrows == 0)
    fatalx ("badbug\n");

  ZALLOC (w1, nblocks, double);
  ZALLOC (w2, nblocks, double);
  ZALLOC (ww, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  ZALLOC (djack, nblocks, double);
  ZALLOC (btop, nblocks, double);
  ZALLOC (bbot, nblocks, double);
  ZALLOC (btop2, nblocks, double);
  ZALLOC (bbot2, nblocks, double);
  ZALLOC (wtop, nblocks, double);
  ZALLOC (wbot, nblocks, double);
  ZALLOC (ww1, nblocks, double);
  ZALLOC (ww2, nblocks, double);
  jmean = initarray_2Ddouble(nblocks, 2, 0) ; 


  totnum = 0;
  for (col = 0; col < ncols; ++col) {
    ret = getf4 (counts[col], xtop, &ytop);
    if (ret < 0)
      continue;
    if (ret == 2)
      fatalx ("bad pop in numerator\n");

    ret = getf4 (counts[col], xtop2, &ytop2);
    if (ret < 0)
      continue;

    if (ret == 2)
      fatalx ("bad pop in numerator\n");

    bnum = bcols[col];
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");

    ybot = ybot2 = 1.0 ; 

    btop[bnum] += ytop;
    bbot[bnum] += 1.0 ; 
    btop2[bnum] += ytop2 ; 
    bbot2[bnum] += 1.0 ;
    ++wjack[bnum];
    ++totnum;

  }

  gtop = asum (btop, nblocks);
  gbot = asum (bbot, nblocks);

  gtop2 = asum (btop2, nblocks);
  gbot2 = asum (bbot2, nblocks);

  mm[0]  = gtop / gbot;
  mm[1]  = gtop2 / gbot2;

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    wtop[k] = gtop - top;
    wbot[k] = gbot - bot;
    wbot[k] += 1.0e-10;
    jmean[k][0] = wtop[k] / wbot[k];
    top = btop2[k];
    bot = bbot2[k];
    wtop[k] = gtop2 - top;
    wbot[k] = gbot2 - bot;
    wbot[k] += 1.0e-10;
    jmean[k][1] = wtop[k] / wbot[k];
  }

  wjackvest(q4stat, q4var, 2, mm, jmean, wjack, nblocks) ; 

  free (w1);
  free (w2);
  free (ww);

  free (wtop);
  free (wbot);
  free (djack);
  free (wjack);

  free (btop);
  free (bbot);
  free (btop2);
  free (bbot2);

  free(ww1) ; 
  free(ww2) ;    
  free2D(&jmean, nblocks) ; 

}


void
doq4diff (double *q4stat, double *q4statsig, int ***counts, int *bcols,
	  int nrows, int ncols, int *xtop, int *xbot, int *xtop2, int *xbot2,
	  int numeg, int nblocks)
// calculates mean+Z of (m1 +m2 -1)
{

  int a, b, c, d;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  double top, bot, *djack, *wjack, gtop, gbot, gtop2, gbot2, *wbot, *wtop;
  double *btop, *bbot, wt;
  double *btop2, *bbot2;
  double ytop, ybot, ytop2, ybot2;
  double y1, y2, yscal;
  double *w1, *w2, *ww, m1, m2;
  int bnum, totnum;
  int ret;

  if (nrows == 0)
    fatalx ("badbug\n");

  ZALLOC (w1, nblocks, double);
  ZALLOC (w2, nblocks, double);
  ZALLOC (ww, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  ZALLOC (djack, nblocks, double);
  ZALLOC (btop, nblocks, double);
  ZALLOC (bbot, nblocks, double);
  ZALLOC (btop2, nblocks, double);
  ZALLOC (bbot2, nblocks, double);
  ZALLOC (wtop, nblocks, double);
  ZALLOC (wbot, nblocks, double);

  totnum = 0;
  for (col = 0; col < ncols; ++col) {
    ret = getf4 (counts[col], xtop, &ytop);
    if (ret < 0)
      continue;
    if (ret == 2)
      fatalx ("bad pop in numerator\n");
    ret = getf4 (counts[col], xbot, &ybot);
    if (ret < 0)
      continue;

    ret = getf4 (counts[col], xtop2, &ytop2);
    if (ret < 0)
      continue;
    if (ret == 2)
      fatalx ("bad pop in numerator\n");
    ret = getf4 (counts[col], xbot2, &ybot2);
    if (ret < 0)
      continue;

    bnum = bcols[col];
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");

    btop[bnum] += ytop;
    bbot[bnum] += ybot;
    btop2[bnum] += ytop2;
    bbot2[bnum] += ybot2;
    ++wjack[bnum];
    ++totnum;

  }

  gtop = asum (btop, nblocks);
  gbot = asum (bbot, nblocks);

  gtop2 = asum (btop2, nblocks);
  gbot2 = asum (bbot2, nblocks);

  mean = gtop / gbot;
  mean += gtop2 / gbot2;
  mean -= 1.0;

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    wtop[k] = gtop - top;
    wbot[k] = gbot - bot;
    wbot[k] += 1.0e-10;
    djack[k] = wtop[k] / wbot[k];
    top = btop2[k];
    bot = bbot2[k];
    wtop[k] = gtop2 - top;
    wbot[k] = gbot2 - bot;
    wbot[k] += 1.0e-10;
    djack[k] += wtop[k] / wbot[k];	// delete-block estimate
    djack[k] -= 1.0;
  }

  weightjack (&jest, &jsig, mean, djack, wjack, nblocks);

  *q4stat = jest;
  *q4statsig = jsig;

  free (w1);
  free (w2);
  free (ww);

  free (wtop);
  free (wbot);
  free (djack);
  free (wjack);

  free (btop);
  free (bbot);
  free (btop2);
  free (bbot2);

}

int
getf4ssize (int **xx, int *indx, double *pabba, double *pbaba, double *ssize)
{

  int a, i;
  double y0, y1, ytot, ff[4];
  double h0, h1;
  double xc[4][2];
  double yabba, ybaba;

  *pabba = *pbaba = *ssize = 0;
  for (i = 0; i < 4; ++i) {
    a = indx[i];
    if (a < 0) {
      return 2;
    }
    y0 = (double) xx[a][0];
    y1 = (double) xx[a][1];
    ytot = y0 + y1;
    if (ytot <= 0.0)
      return -1;
    xc[i][0] = y0 / ytot;
    xc[i][1] = y1 / ytot;
    ff[i] = y0 / ytot;
  }
// *f4 = y1 = (ff[0]-ff[1]) * (ff[2]-ff[3]) ;

  ybaba = xc[0][0] * xc[1][1] * xc[2][0] * xc[3][1];
  ybaba += xc[0][1] * xc[1][0] * xc[2][1] * xc[3][0];

  yabba = xc[0][0] * xc[1][1] * xc[2][1] * xc[3][0];
  yabba += xc[0][1] * xc[1][0] * xc[2][0] * xc[3][1];
  *ssize = yabba + ybaba;
  *pbaba = ybaba;
  *pabba = yabba;
  return 1;
}

void
mkfmat (double *zmat, double *zampsize, int ***counts, int *bcols,
	int nrows, int ncols, int *xtop, int numeg, int nblocks)
{

  int a, b, c, d, s;
  int k, g, i, col, j;
  double ya, yb, y, jest, jsig, mean;
  double top, bot, *djack, *wjack, gtop, gbot, *wbot, *wtop;
  double *btop, *bbot, wt;
  double ytop, ybot;
  double y1, y2, yscal;
  double *w1, *w2, *ww, m1, m2, ff;
  int bnum, totnum;
  int ret;
  double yabba, ybaba, ytot, zabba, zbaba;

  if (nrows == 0)
    fatalx ("badbug\n");

  ZALLOC (w1, nblocks, double);
  ZALLOC (w2, nblocks, double);

  vzero (zmat, nblocks);
  vzero (zampsize, nblocks);
// initialize to zero 
  totnum = 0;
  for (col = 0; col < ncols; ++col) {
    ret = getf4ssize (counts[col], xtop, &zabba, &zbaba, &ytot);
//  if (col<20) printf("zzret: %d %d %9.3f %9.3f\n", col, ret, zabba, zbaba) ;
    if (ret < 0)
      continue;
    if (ret == 2)
      fatalx ("bad pop\n");

//  if (ytot<0.1) continue ;  // ytot is integer #ABBA + # BABA ; otherwise not used
    bnum = bcols[col];
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");

    w1[bnum] += (zbaba - zabba);
    y1 = zabba + zbaba;
    w2[bnum] += y1;
    zampsize[bnum] += y1;
/**
    s = bnum ; 
    if (fabs(y1) > 0) { 
       if (s<=2) printf("zznzd %d %d %9.3f %9.3f %9.3f %9.3f\n", col, s, zbaba, zabba, w1[s], w2[s]) ;
    }
*/
  }

  vsp (w2, w2, 1.0e-10, nblocks);	// avoid divide by zerp
  vvd (zmat, w1, w2, nblocks);	// fstat by blocks;  This is D-statistic

  free (w1);
  free (w2);


}

void
pmoment (double x, double a, double b, double *pm1, double *pvar)
{
// x = y + n   y has prior mean 0, variance a;  n has variance b
  double c1, c2;
  c1 = -2 * x / b;
  c2 = (1 / a) + (1 / b);

  *pm1 = -c1 / (2 * c2);
  *pvar = 1 / c2;

}

double
dd1 (double a, double *b, double *data, int n)
{
  double *v;
  double y, y1, y2;
  int j;

  ZALLOC (v, n, double);
  vsp (v, b, a, n);
  y = 0.0;

  for (j = 0; j < n; ++j) {
    y1 = -0.5 / v[j];;
    y2 = 0.5 * data[j] * data[j] / (v[j] * v[j]);
    y += (y1 + y2);
  }
  free (v);
  return y;
}

double
dd2 (double a, double *b, double *data, int n)
{
  double *v;
  double y, y1, y2;
  int j;

  ZALLOC (v, n, double);
  vsp (v, b, a, n);
  y = 0.0;

  for (j = 0; j < n; ++j) {
    y1 = 0.5 / (v[j] * v[j]);;
    y2 = -data[j] * data[j] / (v[j] * v[j] * v[j]);
    y += (y1 + y2);
  }
  free (v);
  return y;
}

double
scorit (double a, double *b, double *data, int n)
{
  double *v;
  double y, y1, y2;
  int j;

  ZALLOC (v, n, double);
  vsp (v, b, a, n);
  y = 0.0;

  for (j = 0; j < n; ++j) {
    y1 = -0.5 * log (v[j]);
    y2 = -0.5 * data[j] * data[j] / v[j];
    y += (y1 + y2);
  }
  free (v);
  return y;
}

double
optit (double ainit, double *b, double *data, double n)
{
  double aa, aold, anew;
  double ybase, ylast, y, ynewt;
  double d1, d2, y2, mm, vv;
  int iter, j, bad;
  int noclimb = 0;

  aa = aold = ainit;
  ybase = ylast = scorit (aa, b, data, n);
  if (isnan (ybase))
    fatalx ("trouble\n");
  ylast -= ybase;

  for (iter = 1; iter <= 100; ++iter) {
    d1 = dd1 (aold, b, data, n);
    d2 = dd2 (aold, b, data, n);
    anew = aold - (d1 / d2);
    bad = NO;
    ynewt = -100;
    if (anew <= 0.0)
      bad = YES;
    else {
      ynewt = y = scorit (anew, b, data, n) - ybase;
      if (isnan (y))
	fatalx ("trouble\n");
      if (y < ylast)
	bad = YES;
    }
    if (bad == YES) {
// EM iteration
      y2 = 0.0;
      for (j = 0; j < n; ++j) {
	pmoment (data[j], aold, b[j], &mm, &vv);
	y2 += (mm * mm);
	y2 += vv;
	if (isnan (y2))
	  fatalx ("trouble\n");
      }
      anew = y2 / (double) n;
      y = scorit (anew, b, data, n) - ybase;
    }
    printf ("optit: %4d %d %9.3f %15.9f %15.9f\n", iter, bad, anew, y, ynewt);
    if (fabs (ylast - y) < .001)
      ++noclimb;
    aold = anew;
    ylast = y;
    if (noclimb > 3)
      break;
  }

  return aold;

}

void
update_b1 (double *zfmat, double *zsampsize, double zffval, double *aweight,
	   double *bweight1, int nblocks)
{

  int k, i, kk;
  double t1, t2, tb, d1, d2, ssize, bb, yscal;
  double *zd, *v2;

  ZALLOC (zd, nblocks, double);
  ZALLOC (v2, nblocks, double);

  for (k = 0, kk = 0; k < nblocks; ++k) {
    ssize = zsampsize[k];
    if (ssize <= 0.0)
      continue;
    yscal = sqrt (ssize);
    zd[kk] = yscal * (zfmat[k] - zffval);
    v2[kk] = aweight[k] * yscal;
    ++kk;
  }
  t2 = *bweight1;
  tb = optit (t2, v2, zd, kk);

  *bweight1 = tb;

  free (zd);
  free (v2);
}
