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
#include <globals.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "regsubs.h"
#include "egsubs.h"
#include "qpsubs.h"

/** 
 incorporates a jackknife to estimate branch lengths of a tree in f_st 
 space together with standard errors.  This is a test for whether a simple 
 phylogenetic split is supported by the data 
*/


#define WVERSION   "980"
// clade hits and misses (migrations?)
// forcclade added
// outpop NONE forced 
// print number of samples / pop
// popfilename added 
// bbestmode  
// No abc 
// f4mode = YES (no normalization) 
// printsd
// nochrom added
// BAD bug fix (sign of dstat)
// repeated pops trapped in popfilename quads
// id2pops added
// nzdata;  require 10 blocks with abba/baba counts
// xmode ;  call countpopsx which deals with gender on X
// syntactic sugar (strip :) in popfilename  
// numchrom added
// formatting bug (overflow of abba counts) fixed
// instem added 
// cputime etc printed
// error message if nzdata no large enough

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *trashdir = "/var/tmp";
int colcalc = YES;
// int fstdetails = NO ;
int qtmode = NO;
int hires = NO;
int printsd = NO;
int fourier = NO ; 

Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;
int isinit = NO;
int markerscore = NO;
int seed = 0;
int chisqmode = NO;             // approx p-value better to use F-stat
int missingmode = NO;
int dotpopsmode = YES;
int noxdata = YES;              /* default as pop structure dubious if Males and females */
int xdata = NO;              /* default as pop structure dubious if Males and females */
int pcorrmode = NO;
int pcpopsonly = YES;
int nostatslim = 10;
int znval = -1;
int popsizelimit = -1;
int gfromp = NO;                // genetic distance from physical 
int msjack = NO;
char *msjackname = NULL;
int msjackweight = YES;         // weighted jackknife
int bankermode = NO;
int forceclade = NO;
int bbestmode = YES;
int numbanker = 0;
int xchrom = -1;
int zchrom = -1;
int f4mode = NO;
int xmode = NO;
int printssize = YES;
int locount = -1, hicount = 9999999;
// if bankermode  bankers MUST be in quartet ;  at least one type 1 in quartet

int jackweight = YES;
double jackquart = -1.0;

double plo = .001;
double phi = .999;
double pvhit = .001;
double pvjack = 1.0e-6;
double blgsize = 0.05;          // block size in Morgans */
double *chitot;
int *xpopsize;

char *instem = NULL ; 
char *genotypename = NULL;
char *snpname = NULL;
char *indivname = NULL;
char *snpoutfilename = NULL;
char *badsnpname = NULL;
char *poplistname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
char *blockname = NULL;
char *outpop = NULL;
// list of outliers
int fullmode = NO;
int inbreed = NO;
double lambdascale;
int *f2ind, *ind2f, ng2, nh2;

int dotree = NO;

char *outputname = NULL;
char *weightname = NULL;
char *id2pops = NULL;
FILE *ofile;

void readcommands (int argc, char **argv);
double getrs (double *f4, double *f4sig, double *f2, double *f2sig,
              int a, int b, int c, int d, int numeg, double *rho);
int islegal (int *xind, int *types, int mode, int n);
int isclade (int *rr, int *zz);
void setabx (double **abx, double **bax, int ***countcols, int ncols,
             int numeg);
void setf2 (double **f2, int ***countcols, int ncols, int numeg);
void getabc (double *dzscx, double *dzsc, double **abx, double **bax,
             int ncols, int a, int b, int c, int numeg, int *bcols,
             int nblocks);
void getdsc (double *dzscx, double *dzsc, double **abx, double **bax,
             int ncols, int a, int b, int c, int d, int numeg, int *bcols,
             int nblocks, int *sz);
void gettreelen (double *tlenz, double *tlen, double **f2, double **abx,
                 double **bax, int ncols, int *ttt, int numeg, int *bcols,
                 int nblocks);
void regesttree (double *ans, double *xn, double xd);
void ckdups (char **list, int n);
void printbadql (char ***qlist, int tt, char **eglist, int numeg);

int usage (char *prog, int exval);


int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int **snppos;
  int *snpindx;
  char **snpnamelist, **indnamelist;
  char **eglist;
  int *nsamppops;
  int *ztypes;
  int lsnplist, lindlist, numeg;
  int i, j, k, k1, k2, k3, k4, kk;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;
  double y1, y2, y, sig, tail, yy1, yy2;
  char ss[11];
  int *blstart, *blsize, nblocks;
  int xnblocks;                 /* for xsnplist */
  int *bcols;
  double maxgendis;
  int xind[4];

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
  double *xmean, *xfancy;
  double *divans, *divsd;
  double *hettop, *hetbot;
  int chrom, numclear;
  double gdis;
  int outliter, *badlist, nbad;
  double **zdata, *z1, *z2;
  int maxtag = -1;
  double **zz;
  double *pmean, *pnum, rscore[3], dstat[3], rrr[3], ww[3],
    serr[3];
  int ssize[3][3], *sz;
  int tpat[3][4], rpat[3][4], *rrtmp, *rp;
  int *rawcol;;
  int a, b, c, d, col;
  int aa, bb, cc, dd;
  double *qpscores;
  double *hest, *hsig;
  double mingenpos, maxgenpos;
  int *qhit;                    /* number of times pair is clade in quartet */
  int *qmiss;                   /* number of times pair migration event implied */
  int **qplist, numqp = 0, maxqp = 10000;
  double *qpscore;
  char ***qlist, *sx;
  int nqlist = 0;
  int bbest[3];
  double absscore[3];
  double ascore[4], astat[4];


  double **dsctop, **dscbot;
  double **abx, **bax, **f2;
  int popx[4];
  double tn[4 * 5], td[4 * 4];
  double zzsig[5], zzest[5], zsc[5];
  double ymin;

  double *f3, *f4, *f3sig, *f4sig;
  int t1, t2, tt;
  int ***counts, **ccc;

  double tlenz[5], tlen[5];
  int lenz[5];
  double ymem ; 


  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;
  
  printf ("## qpDstat version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;
  if ((poplistname == NULL) && (popfilename == NULL))
    fatalx ("poplistname, popfilename both null\n");

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 

  if (xdata) noxdata = NO ; 
  if (xchrom == (numchrom + 1)) noxdata = NO;

  if (!bankermode)
    forceclade = NO;
//if (fancynorm) printf("fancynorm used\n") ;
//else printf("no fancynorm used\n") ;
  setjquart (NO, jackweight, jackquart);

  if (noxdata == NO) printf("X used!\n") ;
  nostatslim = MAX (nostatslim, 3);

  setinbreed (inbreed);

  if (outputname != NULL)
    openit (outputname, &ofile, "w");

  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  numindivs = getindivs (indivname, &indivmarkers);
  if (id2pops != NULL) {
    setid2pops (id2pops, indivmarkers, numindivs);
  }
  setindm (indivmarkers);

  k = getgenos (genotypename, snpmarkers, indivmarkers,
                numsnps, numindivs, nignore);

  ZALLOC (eglist, numindivs, char *);
  ZALLOC (ztypes, numindivs, int);
  if (popfilename == NULL) {
    if (bankermode == NO)
      numeg = loadlist (eglist, poplistname);
    else {
      numeg = loadlist_type (eglist, poplistname, ztypes, 0);
      numbanker = 0;
      for (k = 0; k < numeg; ++k) {
        if (ztypes[k] == 2)
          ++numbanker;
      }
      printf ("bankermode: ");
      printimat (ztypes, 1, numeg);
    }
  }
  if (popfilename != NULL) {
    bbestmode = NO;
    nqlist = numlines (popfilename);
    ZALLOC (qlist, 4, char **);
    for (k = 0; k < 4; ++k) {
      ZALLOC (qlist[k], nqlist, char *);
    }
    nqlist =
      getnamesstripcolon (&qlist, nqlist, 4, popfilename, locount, hicount);
    numeg = 0;
    printf ("number of quadruples %d\n", nqlist);

    if (nqlist==0) { 
     printf("no quads!\n") ; 
     return 0 ;
    }

    fflush (stdout);
    for (k = 0; k < 4; ++k) {
      for (j = 0; j < nqlist; ++j) {
        sx = qlist[k][j];
        t1 = indxstring (eglist, numeg, sx);
        if (t1 >= 0)
          continue;
        eglist[numeg] = strdup (sx);
        ++numeg;
        setstatus (indivmarkers, numindivs, sx);
      }
    }
  }

  ckdups (eglist, numeg);

  ZALLOC (nsamppops, numeg, int);

  for (i = 0; i < numindivs; i++) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    k = indxindex (eglist, numeg, indx->egroup);
    if (k < 0)
      continue;
    ++nsamppops[k];
  }

  if (popfilename == NULL)
    seteglist (indivmarkers, numindivs, poplistname);   // null if poplistname == NULL

  for (i = 0; i < numeg; i++) {
    t = nsamppops[i];
    if (t == 0)
      fatalx ("No samples: %s\n", eglist[i]);
    printf ("%3d %20s %4d\n", i, eglist[i], t);
  }

  printf ("jackknife block size: %9.3f\n", blgsize);

  outpop = strdup ("NONE");
  outnum = 999;

// copied from qp3Pop.c 
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


  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);
  ZALLOC (rawcol, numindivs, int);
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


  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    cupt->weight = 1.0;
  }
  numsnps = rmsnps (snpmarkers, numsnps, NULL); //  rid ignorable snps
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

/**
  tt = xsnplist[1000] -> tagnumber ; 
  if (tt>nblocks) fatalx("yuck!\n") ; 
*/

  printf ("nrows, ncols: %d %d\n", nrows, ncols);
  ZALLOC (counts, ncols, int **);
  for (k = 0; k < ncols; ++k) {
    counts[k] = initarray_2Dint (numeg, 2, 0);
  }


  if (xmode && (xchrom == (numchrom + 1))) {
    countpopsx (counts, xsnplist, xindlist, xindex, xtypes, nrows, ncols);
  }
  else {
    countpops (counts, xsnplist, xindex, xtypes, nrows, ncols);
  }

  if (verbose) {
    for (k = 0; k < MIN (ncols, 100); ++k) {
      printf ("zzcount: %3d ", k);
      for (j = 0; j < numeg; ++j) {
        printf ("%3d %3d  ", counts[k][j][0], counts[k][j][1]);
      }
      printnl ();
    }
  }

  ng2 = numeg * numeg;
  nh2 = numeg * (numeg - 1);
  nh2 /= 2;

  ZALLOC (f2ind, numeg * numeg, int);
  ZALLOC (ind2f, nh2, int);

  k = 0;
  for (a = 0; a < numeg; ++a) {
    for (b = a + 1; b < numeg; ++b) {
      ind2f[k] = a * numeg + b;
      f2ind[a * numeg + b] = k;
      f2ind[b * numeg + a] = k;
      ++k;
    }
  }

  abx = initarray_2Ddouble (nh2, ncols, 0.0);
  bax = initarray_2Ddouble (nh2, ncols, 0.0);
  if (dotree)
    f2 = initarray_2Ddouble (nh2, ncols, 0.0);

  ZALLOC (bcols, ncols, int);   // blocks for columns -1 => unused
  ivclear (bcols, -1, ncols);
  for (k = 0; k < ncols; k++) {
    cupt = xsnplist[k];
    tt = bcols[k] = cupt->tagnumber;
    if (tt>nblocks) fatalx("bad tagnumber!\n") ; 
  }

  vzero(rscore,3) ; 
  vzero(dstat,3) ; 
   vzero(rrr,3) ; 
   vzero(ww,3) ; 



  setabx (abx, bax, counts, ncols, numeg);
  if (dotree)
    setf2 (f2, counts, ncols, numeg);

  ivzero (bbest, 3);
  for (a = 0; a < numeg; a++) {
    if (popfilename != NULL)
      break;
    for (b = a + 1; b < numeg; b++) {
      for (c = b + 1; c < numeg; c++) {
        for (d = c + 1; d < numeg; d++) {

          xind[0] = a;
          xind[1] = b;
          xind[2] = c;
          xind[3] = d;

          t = islegal (xind, ztypes, bankermode, numeg);
          if (t == NO)
            continue;


          getdsc (&rscore[0], &dstat[0], abx, bax, ncols,
                  a, b, c, d, numeg, bcols, nblocks, ssize[0]);

          getdsc (&rscore[1], &dstat[1], abx, bax, ncols,
                  a, c, b, d, numeg, bcols, nblocks, ssize[1]);

          getdsc (&rscore[2], &dstat[2], abx, bax, ncols,
                  a, d, b, c, numeg, bcols, nblocks, ssize[2]);

          for (k = 0; k < 3; ++k) {
            serr[k] = 1.0;
            if (rscore[k] != 0)
              serr[k] = dstat[k] / rscore[k];
          }

          xcopy (rpat[0], a, b, c, d);
          xcopy (rpat[1], a, c, b, d);
          xcopy (rpat[2], a, d, b, c);

          ivzero (bbest, 3);

          if (bbestmode) {
            vabs (absscore, rscore, 3);
            vlmaxmin (absscore, 3, NULL, &k);
            bbest[k] = 1;
            rp = rpat[k];
            aa = rp[0];
            bb = rp[1];
            cc = rp[2];
            dd = rp[3];
            getabc (&astat[0], &ascore[0], abx, bax, ncols,
                    aa, bb, cc, numeg, bcols, nblocks);
            getabc (&astat[1], &ascore[1], abx, bax, ncols,
                    aa, bb, dd, numeg, bcols, nblocks);
            getabc (&astat[2], &ascore[2], abx, bax, ncols,
                    cc, dd, aa, numeg, bcols, nblocks);
            getabc (&astat[3], &ascore[3], abx, bax, ncols,
                    cc, dd, bb, numeg, bcols, nblocks);

          }

          vclip (rscore, rscore, -100.0, 100.0, 3);

          for (k = 0; k < 3; ++k) {
            copyiarr (rpat[k], popx, 4);
            printf ("result: ");
            sz = ssize[k];

            if (bankermode) {
              for (i = 0; i < 4; i++) {
                t = popx[i];
                printf ("%1d", ztypes[t]);
              }
              printf ("  ");
            }
            for (i = 0; i < 4; i++) {
              t = popx[i];
              printf ("%10s ", eglist[t]);
            }

            if (f4mode == NO)
              printf (" %10.4f", dstat[k]);
            if (f4mode == YES)
              printf (" %12.6f", dstat[k]);

            if (printsd) {
              printf (" %12.6f ", serr[k]);
            }

            printf (" %9.3f ", rscore[k]);
            if (dotree) {
              gettreelen (tlenz, tlen, f2, abx, bax, ncols, rpat[k],
                          numeg, bcols, nblocks);
              ivzero (lenz, 5);
              for (i = 0; i < 5; i++) {
                if (tlenz[i] > 0)
                  continue;
                lenz[i] = nnint (-10.0 * tlenz[i]);
                lenz[i] = MIN (lenz[i], 9);
              }
              for (i = 0; i < 5; i++) {
                printf ("% 1d", lenz[i]);
              }
            }
            if (bbest[k] == 1) {
              printf ("%6s", " best ");
            }
            else {
              printf ("%6s", " - ");
            }
            if (printssize) {
              printf (" %6d ", sz[0]);
              printf ("%6d ", sz[1]);
              printf ("%6d ", sz[2]);
            }
            printnl ();
          }

        }
      }
    }
  }


  for (tt = 0; tt < nqlist; ++tt) {


    a = indxindex (eglist, numeg, qlist[0][tt]);
    b = indxindex (eglist, numeg, qlist[1][tt]);
    c = indxindex (eglist, numeg, qlist[2][tt]);
    d = indxindex (eglist, numeg, qlist[3][tt]);

    if ((a == b) || (a == c) || (a == d) || (b == c) || (b == d) || (c == d)) {
      printbadql (qlist, tt, eglist, numeg);
      continue;
    }


    vzero(rscore, 3) ; 
    getdsc (&rscore[0], &dstat[0], abx, bax, ncols,
            a, b, c, d, numeg, bcols, nblocks, ssize[0]);

    serr[0] = 1.0;
    if (rscore[0] != 0.0)
      serr[0] = dstat[0] / rscore[0];

    xcopy (rpat[0], a, b, c, d);


    vclip (rscore, rscore, -100.0, 100.0, 3);

    k = 0;
    sz = ssize[k];

//       printf("zz1 %d ", ncols) ;  printimat(sz, 1, 3) ;

    copyiarr (rpat[k], popx, 4);
    printf ("result: ");

    for (i = 0; i < 4; i++) {
      t = popx[i];
      printf ("%10s ", eglist[t]);
    }

    if (f4mode == NO)
      printf (" %10.4f ", dstat[k]);
    if (f4mode == YES)
      printf (" %12.6f ", dstat[k]);

    if (printsd) {
      printf (" %12.6f", serr[k]);
    }
    printf (" %9.3f ", rscore[k]);

    if (dotree) {
      gettreelen (tlenz, tlen, f2, abx, bax, ncols, rpat[k],
                  numeg, bcols, nblocks);
      ivzero (lenz, 5);
      for (i = 0; i < 5; i++) {
        if (tlenz[i] > 0)
          continue;
        lenz[i] = nnint (-10.0 * tlenz[i]);
        lenz[i] = MIN (lenz[i], 9);
      }
      for (i = 0; i < 5; i++) {
        printf ("% 1d", lenz[i]);
      }
    }

    if (printssize) {
      printf (" %6d ", sz[0]);
      printf ("%6d ", sz[1]);
      printf ("%6d ", sz[2]);
    }
    printnl ();

  }


  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qpDstat: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;
}

void
printbadql (char ***qlist, int tt, char **eglist, int numeg)
{
  int a, b, c, d;

  a = indxindex (eglist, numeg, qlist[0][tt]);
  b = indxindex (eglist, numeg, qlist[1][tt]);
  c = indxindex (eglist, numeg, qlist[2][tt]);
  d = indxindex (eglist, numeg, qlist[3][tt]);

  printf ("***warning: repeated population: ");
  printf ("%s ", eglist[a]);
  printf ("%s :  ", eglist[b]);
  printf ("%s ", eglist[c]);
  printf ("%s ", eglist[d]);
  printnl ();
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
  while ((i = getopt (argc, argv, "L:H:p:hvV")) != -1) {

    switch (i) {

    case 'h':
      usage(basename(argv[0]), 0);
      break ;  

    case 'L':
      locount = atoi (optarg);
      break;

    case 'H':
      hicount = atoi (optarg);
      break;

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
  getstring (ph, "poplistname:", &poplistname);
  getstring (ph, "popfilename:", &popfilename);
  //  getstring(ph, "outpop:", &outpop) ;
  getstring (ph, "output:", &outputname);
  getstring (ph, "outputname:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "msjackname:", &msjackname);
  getstring (ph, "id2pops:", &id2pops);
  getint (ph, "msjackweight:", &msjackweight);
  getdbl (ph, "blgsize:", &blgsize);
  getint (ph, "dotree:", &dotree);
  getint (ph, "nochrom:", &zchrom);
  getint (ph, "locount:", &locount);
  getint (ph, "hicount:", &hicount);

  getint (ph, "numchrom:", &numchrom);
  getint (ph, "noxdata:", &noxdata);
  getint (ph, "xdata:", &xdata);
  getint (ph, "colcalc:", &colcalc);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "printssize:", &printssize);

  getint (ph, "nostatslim:", &nostatslim);
  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);      // gen dis from phys
  getint (ph, "fullmode:", &fullmode);
  getint (ph, "bankermode:", &bankermode);
  getint (ph, "forceclade:", &forceclade);      // in bankermode one clade must be all banker = 2
  getint (ph, "bbestmode:", &bbestmode);        // print "best" for triples
// getint(ph, "fancynorm:", &fancynorm) ; 
// getint(ph, "usenorm:", &fancynorm) ;  /* changed 11/02/06 */
  getdbl (ph, "jackquart:", &jackquart);
  getint (ph, "jackweight:", &jackweight);
  getint (ph, "chrom:", &xchrom);
  getint (ph, "xmode:", &xmode);
  getint (ph, "f4mode:", &f4mode);
  getint (ph, "printsd:", &printsd);
  getstring (ph, "blockname:", &blockname);
  getint (ph, "fourierjack:", &fourier);


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

  clear2D (&abx, nh2, ncols, -1.0);
  clear2D (&bax, nh2, ncols, -1.0);

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

void
setf2 (double **f2, int ***counts, int ncols, int numeg)
{


  int i, j, k, t1, t2, a, h;
  double y1, y2, x0, x1;
  int **ccc;
  double hest, htest, h1, yt1, h2, yt2, yy;

  clear2D (&f2, nh2, ncols, -10.0);

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

        yy = (y1 - y2) * (y1 - y2);
        x0 = ccc[i][0];
        x1 = ccc[i][1];
        yt1 = t1;
        hest = h1 = x0 * x1 / (yt1 * (yt1 - 1));
        htest = h1 / yt1;
        yy -= htest;

        x0 = ccc[j][0];
        x1 = ccc[j][1];
        yt1 = t2;
        hest = h2 = x0 * x1 / (yt1 * (yt1 - 1));
        htest = h2 / yt1;
        yy -= htest;

        h = f2ind[i * numeg + j];
        f2[h][k] = yy;

      }
    }
  }
}

void
getabc (double *dzscx, double *dzsc, double **abx, double **bax, int ncols,
        int a, int b, int c, int numeg, int *bcols, int nblocks)
{
  double *ztop, *zbot, *znum, *jmean, *jwt;
  int k, s, h1, h2;
  double y1, y2, yabba, ybaba, yt, yb, ytop, ybot, ymatchd;
  double ymean, est, sig;

  *dzscx = *dzsc = 0.0;

  fflush (stdout);

  ZALLOC (ztop, nblocks, double);
  ZALLOC (zbot, nblocks, double);
  ZALLOC (znum, nblocks, double);
  ZALLOC (jmean, nblocks, double);
  ZALLOC (jwt, nblocks, double);

  h1 = f2ind[a * numeg + c];
  h2 = f2ind[b * numeg + c];
  for (k = 0; k < ncols; ++k) {
    s = bcols[k];
    if (s < 0)
      continue;
    y1 = abx[h1][k] + bax[h1][k];       // mismatch (a,c)
    y2 = abx[h2][k] + bax[h2][k];       // mismatch (b,c)

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    ymatchd = y2 - y1;          // match(a,c) -match(b,c) 

    ztop[s] += ymatchd;
    zbot[s] += 1.0;
    znum[s] += 1.0;

  }
  ytop = asum (ztop, nblocks);
  ybot = asum (zbot, nblocks);
  ybot += 1.0e-20;
  ymean = ytop / ybot;
// printf("zz %d %d %d %12.6f\n", a, b, c, ymean) ;
  for (s = 0; s < nblocks; ++s) {
    yt = ytop - ztop[s];
    yb = ybot - zbot[s];
    jmean[s] = yt / yb;
    jwt[s] = znum[s];
  }

  weightjack (&est, &sig, ymean, jmean, jwt, nblocks);
  *dzscx = est / (sig + 1.0e-20);       // Z score
  *dzsc = est;

  free (ztop);
  free (zbot);
  free (znum);
  free (jmean);
  free (jwt);
}



void
getdsc (double *dzscx, double *dzsc, double **abx, double **bax, int ncols,
        int a, int b, int c, int d, int numeg, int *bcols, int nblocks,
        int *sz)
{
  double *ztop, *zbot, *znum, *jmean, *jwt;
  int k, s, h1, h2;
  double y, y1, y2, yabba, ybaba, yt, yb, ytop, ybot;
  double ymean, est, sig;
  double sgn;
  int ntot, nbaba, nabba, nzdata = 0;
  int *nzd;
// nzdata is # blocks with  nonzero abba or baba counts 
  double zbaba, zabba;
  double rho, mblocks ;

  *dzscx = *dzsc = 0.0;
  sgn = 1.0;
  ntot = nbaba = nabba = 0;
  zbaba = zabba = 0;
  if (sz != NULL)
    ivzero (sz, 3);

  if (a > b)
    sgn *= -1.0;
  if (c > d)
    sgn *= -1.0;

  fflush (stdout);

  ZALLOC (nzd, nblocks, int);
  ZALLOC (ztop, nblocks, double);
  ZALLOC (zbot, nblocks, double);

  ZALLOC (znum, nblocks, double);
  ZALLOC (jmean, nblocks, double);
  ZALLOC (jwt, nblocks, double);

  h1 = f2ind[a * numeg + b];
  h2 = f2ind[c * numeg + d];
  for (k = 0; k < ncols; ++k) {
    s = bcols[k];
    if (s < 0)
      continue;
    y1 = abx[h1][k];
    y2 = bax[h2][k];

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    yabba = y1 * y2;

    y1 = bax[h1][k];
    y2 = abx[h2][k];

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    yabba += y1 * y2;

    y1 = abx[h1][k];
    y2 = abx[h2][k];

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    ybaba = y1 * y2;

    y1 = bax[h1][k];
    y2 = bax[h2][k];

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    ybaba += y1 * y2;

    zbaba += ybaba;
    zabba += yabba;

    ztop[s] += (ybaba - yabba) * sgn;
    if (f4mode == NO)
      zbot[s] += ybaba + yabba;
    if (f4mode == YES)
      zbot[s] += 1.0;

    znum[s] += 1.0;
    ++ntot;

    y1 = fabs (yabba + ybaba);
    if (y1 > 0.0)
      nzd[s] = 1;

/**
     if ((y1>0.0)) {
       printf("zznzd %d %d %9.3f %9.3f %9.3f %9.3f %d %d\n", k, s, ybaba, yabba, ztop[s], zbot[s]) ; 
     }
*/
  }
  nzdata = intsum (nzd, nblocks);
// number of blocks with any data at all 
  free (nzd);
  if (verbose)
    printf ("nzdata:  %d\n", nzdata);
  if (nzdata < 5) {
    free (ztop);
    free (zbot);
    free (znum);
    free (jmean);
    free (jwt);
    *dzscx = 0;                 // Z score
    *dzsc = 0;
    enuf("insufficient number of blocks -- you need 5 blocks with non zero ABBA/BABA program terminating\n") ;
    return;
  }
  ytop = asum (ztop, nblocks);
  ybot = asum (zbot, nblocks);
  ybot += 1.0e-20;
  ymean = ytop / ybot;
  for (s = 0; s < nblocks; ++s) {
    yt = ytop - ztop[s];
    yb = ybot - zbot[s];
    jmean[s] = yt / yb;
    jwt[s] = znum[s];
  }
// printf("## zzmean: %d %12.6f\n", nblocks, ymean) ; 

/**
   for (s=0; s<nblocks; ++s) { 
    y1 = ztop[s] ; 
    y2 = zbot[s] ; 
    printf("zzbug: %d %9.6f %9.6f %9.6f\n", s, y1, y2, y1/(y2+1.0e-20)) ; 
   }
*/
// abort() ;
  mblocks = nblocks ; 
  for (k=nblocks-1; k>=0; --k) { 
   if (jwt[k] > 0) break ;  
   --mblocks ; 
  }

  if (fourier) {
   weightjackfourier (&est, &sig, ymean, jmean, jwt, mblocks, &rho);
   printf("fourier called.  rho: %9.3f\n", rho) ; 
   verbose = NO ; 
  }

  else {
   weightjack (&est, &sig, ymean, jmean, jwt, mblocks);
  }

  *dzscx = est / (sig + 1.0e-20);       // Z score
  *dzsc = est;

  free (ztop);
  free (zbot);
  free (znum);
  free (jmean);
  free (jwt);

  if (sz == NULL)
    return;

  if (sgn < 0) {
    yt = zbaba;
    zbaba = zabba;
    zabba = yt;
  }
  sz[0] = nnint (zbaba);
  sz[1] = nnint (zabba);
  sz[2] = ntot;


//   printf("zz %12.6f ",  est) ;  
//   printimat(sz, 1, 3) ;
}

void
gettreelen (double *tlenz, double *tlen, double **f2, double **abx,
            double **bax, int ncols, int *ttt, int numeg, int *bcols,
            int nblocks)
{
#define NPAR 5
#define ND   6
  double *jmean, *jwt;
  int k, s, h1, h2, hh[6], i;
  double y, y1, y2, yabba, ybaba, yt, yb, ytop, ybot, yden, gden;
  double ymean, jest, jsig;
  double *gn, gd, **xn, *xx, *xden, *qqest, **xqest;
  double *wjack, *djack;
  int a, b, c, d, h;

  a = ttt[0];
  b = ttt[1];
  c = ttt[2];
  d = ttt[3];

  ZALLOC (gn, ND, double);
  ZALLOC (qqest, NPAR, double);
  gd = 0.0;

  xn = initarray_2Ddouble (nblocks, ND, 0.0);
  ZALLOC (xden, nblocks, double);
  ZALLOC (wjack, nblocks, double);
  ZALLOC (djack, nblocks, double);

  xqest = initarray_2Ddouble (nblocks, NPAR, 0.0);

  h1 = hh[0] = f2ind[a * numeg + b];
  hh[1] = f2ind[a * numeg + c];
  hh[2] = f2ind[a * numeg + d];
  hh[3] = f2ind[b * numeg + c];
  hh[4] = f2ind[b * numeg + d];
  h2 = hh[5] = f2ind[c * numeg + d];
  gden = 0;

  for (k = 0; k < ncols; ++k) {

    s = bcols[k];
    if (s < 0)
      continue;

    y1 = abx[h1][k];
    y2 = bax[h2][k];

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    yabba = y1 * y2;
    y1 = bax[h1][k];
    y2 = abx[h2][k];

    if (y1 < 0.0)
      continue;
    if (y2 < 0.0)
      continue;

    yabba += y1 * y2;
    y1 = abx[h1][k];
    y2 = abx[h2][k];
    ybaba = y1 * y2;
    y1 = bax[h1][k];
    y2 = bax[h2][k];
    ybaba += y1 * y2;
    yden = yabba + ybaba;
    xden[s] += yden;
    gden += yden;
    wjack[s] += 1.0;

    for (i = 0; i < ND; i++) {
      h = hh[i];
      y = f2[h][k];
      xn[s][h] += y;
      gn[h] += y;
    }
  }
  regesttree (qqest, gn, gden);
  copyarr (qqest, tlen, NPAR);
  printf ("qqest: ");
  printmatw (qqest, 1, 5, 5);
  verbose = NO;

  for (s = 0; s < nblocks; s++) {
    xx = xn[s];
    vvm (xx, gn, xx, ND);
    regesttree (xqest[s], xx, gden - xden[s]);
  }
  for (a = 0; a < NPAR; a++) {
    for (k = 0; k < nblocks; ++k) {
      djack[k] = xqest[k][a];
    }
    weightjack (&jest, &jsig, qqest[a], djack, wjack, nblocks);
    tlenz[a] = jest / jsig;
    tlen[a] = jest;
  }

  free2D (&xqest, nblocks);
  free2D (&xn, nblocks);
  free (xden);
  free (djack);
  free (wjack);
  free (gn);
  free (qqest);
}

void
regesttree (double *ans, double *xn, double xd)
{
  int a, b, c, k;
  double *co, *rr;
  double f;

  ZALLOC (co, 6 * 5, double);
  ZALLOC (rr, 6, double);

  k = 0;
  a = 0;
  b = 1;
  f = xn[k] / xd;
  co[k * 5 + 0] = co[k * 5 + 1] = 1;
  rr[k] = f;

  k = 1;
  a = 0;
  b = 2;
  f = xn[k] / xd;
  co[k * 5 + 0] = co[k * 5 + 2] = co[k * 5 + 3] = 1;
  rr[k] = f;

  k = 2;
  a = 0;
  b = 3;
  f = xn[k] / xd;
  co[k * 5 + 0] = co[k * 5 + 2] = co[k * 5 + 4] = 1;
  rr[k] = f;

  k = 3;
  a = 1;
  b = 2;
  f = xn[k] / xd;
  co[k * 5 + 1] = co[k * 5 + 2] = co[k * 5 + 3] = 1;
  rr[k] = f;

  k = 4;
  a = 1;
  b = 3;
  f = xn[k] / xd;
  co[k * 5 + 1] = co[k * 5 + 2] = co[k * 5 + 4] = 1;
  rr[k] = f;

  k = 5;
  a = 2;
  b = 3;
  f = xn[k] / xd;
  co[k * 5 + 3] = co[k * 5 + 4] = 1;
  rr[k] = f;

  regressit (ans, co, rr, 6, 5);

  free (co);
  free (rr);

}

void
ckdups (char **list, int n)
{
  int t, i, j;

  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      t = strcmp (list[i], list[j]);
      if (t == 0)
        fatalx ("dup in list: %s\n", list[i]);
    }
  }
}
int usage (char *prog, int exval)
{ 

  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -L <val>    ... use <val> as low value.\n");
  (void)fprintf(stderr, "   -H <val>    ... use <val> as high value.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");

  exit(exval);
};


