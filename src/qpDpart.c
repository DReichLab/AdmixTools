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


#define WVERSION   "100"
// outpop NONE forced 
// print number of samples / pop
// popfilename added 
// bbestmode  
// No abc 
// fmode = YES (no normalization) 
// printsd
// nochrom added
// nzdata;  require 10 blocks with abba/baba counts
// xmode ;  call countpopsx which deals with gender on X
// syntactic sugar (strip :) in popfilename  
// numchrm added

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
char *pattern[2], *rawpattern = NULL ;  
int npattern = 0 ;


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
int fmode = NO;
int xmode = NO;
int printssize = YES;
int locount = -1, hicount = 9999999;
double fscale = 1.0e6 ; 

int jackweight = YES;
double jackquart = -1.0;

double plo = .001;
double phi = .999;
double pvhit = .001;
double pvjack = 1.0e-6;
double blgsize = 0.05;    
double *chitot;
int *xpopsize;

char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *poplistname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
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

void ckdups (char **list, int n);
int setpattern(char **pattern, char *rawpattern)  ; 
void setcolprobs(double **colprobs, int ***counts, int ncols, int numeg) ;
void  patternscore (double *pysc, double *pysig,   char **pattern, double **colprobs, int *egindex, 
 int ncols, int numeg, int *bcols, int nblocks) ;


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
  double *pmean, *pnum, rscore[3], dstat[3], hscore[3], rrr[3], ww[3],
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


  double tn[4 * 5], td[4 * 4];
  double zzsig[5], zzest[5], zsc[5];
  double ymin;

  double *f3, *f4, *f3sig, *f4sig;
  int t1, t2, tt;
  int ***counts, **ccc;

  double tlenz[5], tlen[5];
  int lenz[5];
  double **colprobs ;
  int *egindex ;  
  double ysc, ysig ; 



  readcommands (argc, argv);
  printf ("## qpDpart version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;
  if ((poplistname == NULL) && (popfilename == NULL))
    fatalx ("poplistname, popfilename both null\n");

  if (xchrom == (numchrom + 1)) noxdata = NO;

  nostatslim = MAX (nostatslim, 3);
  npattern = setpattern(pattern, rawpattern) ;

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
  if (popfilename == NULL) fatalx("no popfilename!\n") ;
    nqlist = numlines (popfilename);
    ZALLOC (qlist, npattern, char **);
    for (k = 0; k < npattern; ++k) {
      ZALLOC (qlist[k], nqlist, char *);
    }
    nqlist =
      getnamesstripcolon (&qlist, nqlist, npattern, popfilename, locount, hicount);
    numeg = 0;
    printf ("number of population sets %d\n", nqlist);
    fflush (stdout);
    for (k = 0; k < npattern; ++k) {
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

  ckdups (eglist, numeg); // bugcheck

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

  nblocks = numblocks (snpmarkers, numsnps, blgsize);
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
  printf ("number of blocks for jackknife: %d\n", xnblocks);
  nblocks = xnblocks;

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

  colprobs = initarray_2Ddouble(ncols, numeg, -1) ;   
  setcolprobs(colprobs, counts, ncols, numeg) ;


  ZALLOC (bcols, ncols, int);   // blocks for columns -1 => unused
  ivclear (bcols, -1, ncols);
  for (k = 0; k < ncols; k++) {
    cupt = xsnplist[k];
    bcols[k] = cupt->tagnumber;
  }


  ZALLOC(egindex, npattern, int) ;  

  if (fmode == YES) printf("fstats multiplied by %12.0f\n", fscale) ; 
  for (tt = 0; tt < nqlist; ++tt) {

    for (a=0; a<npattern; ++a) {  
      egindex[a] = indxindex(eglist, numeg, qlist[a][tt]) ; 
    }

    patternscore (&ysc, &ysig,  pattern, colprobs, egindex, ncols, numeg, bcols, nblocks);

    printf ("result: ");

    for (i = 0; i < npattern; i++) {
      t = egindex[i];
      printf ("%10s ", eglist[t]);
    }

    if (fmode) { 
     printf (" %12.0f", ysc*fscale) ; 
    }

    else { 
     printf(" %12.3f", ysc) ; 
    }

    if (printsd) {
      printf (" %12.6f", ysig) ;
    }
    printf (" %9.3f ", ysc/ysig);

    printnl() ; 


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

  while ((i = getopt (argc, argv, "l:h:p:vV")) != -1) {

    switch (i) {

    case 'l':
      locount = atoi (optarg);
      break;

    case 'h':
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
  getint (ph, "printsd:", &printsd);
  getint (ph, "nochrom:", &zchrom);
  getint (ph, "locount:", &locount);
  getint (ph, "hicount:", &hicount);
  getint (ph, "numchrom:", &numchrom);

  getint (ph, "noxdata:", &noxdata);
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
  getint (ph, "fmode:", &fmode);
  getlongstring (ph, "pattern:", &rawpattern);
  getdbl (ph, "fscale:", &fscale);


  printf ("### THE INPUT PARAMETERS\n");
  printf ("##PARAMETER NAME: VALUE\n");
  writepars (ph);

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

int setpattern(char **pattern, char *rawpattern) 
{
  char *s1 ; 
  char *spt[MAXFF] ; 
  int nsplit, l0, l1  ; 

  if (rawpattern == NULL) fatalx("no pattern specified!\n") ;
  s1 = strdup(rawpattern) ;  
  substring(&s1, " ", "") ; 
  substring(&s1, ":", " ") ; 
  nsplit = splitup(s1, spt, MAXFF) ; 
  if (nsplit != 2) fatalx("bad pattern %s\n", rawpattern) ;
  l0 = strlen(spt[0]) ;
  l1 = strlen(spt[1]) ;
  if (l0 != l1)  fatalx("bad pattern(length mismatch) %s\n", rawpattern) ;
  pattern[0] = strdup(spt[0]) ; 
  pattern[1] = strdup(spt[1]) ; 
  freeup(spt, nsplit) ; 
  return l1 ;

}
void setcolprobs(double **colprobs, int ***counts, int ncols, int numeg) 
{
// set colprobs;  polarity of allele is irrelevant as patternscore will try both ways
  int k, j, c1, c2, t ;

  clear2D(&colprobs, ncols, numeg, -99) ; 
  for (k=0; k<ncols; ++k) { 
   for (j=0; j<numeg; ++j) {
    c1 = counts[k][j][0] ;
    c2 = counts[k][j][1] ;
    t = c1 + c2 ; 
    if (t<=0) continue ; 
    colprobs[k][j] = (double) c1 / (double) t ; 
   }
  }
}

void spat(char *cpat, int *ipat)  
{
  int k, n, t ; 
  n = strlen(cpat) ;  

  for (k=0; k<n; ++k) { 
   t = ipat[k] = (int) cpat[k] - (int) 'A' ;
   if (t<0) fatalx("bad pattern %s\n", cpat) ; 
   if (t>1) fatalx("bad pattern %s\n", cpat) ; 
  }


  return ; 

}

double abbaprob(double *cp, int *pat, int lenpat)
{
  double y1, y2, p, q ; 
  int k, s ; 
     
  y1 = y2 = 1 ;
  for (k=0; k<lenpat;++k) { 
   p = cp[k] ; 
   q = 1-p ; 
   s = pat[k] ;  
   if (s==0) { 
    y1 *= p ; 
    y2 *= q ; 
   } 
   if (s==1) { 
    y1 *= q ; 
    y2 *= p ; 
   } 
  }
  return y1 + y2 ;
}

void  patternscore (double *pysc, double *pysig,  char **pattern, double **colprobs, int *egindex, 
 int ncols, int numeg, int *bcols, int nblocks) 

{

 int npattern ; 
 double *cp, *wp ; 
 int *pat1, *pat2 ; 
 int k, j, s, t ; 
 double *ztop, *zbot, *znum, *jmean, *jwt ; 
 double y, y1, y2, yabba, ybaba, yt, yb, ytop, ybot;
 double ymean, est, sig;
 double z1, z2, zsum=0, ymin ; 

 npattern = strlen(pattern[0]) ;  
 ZALLOC(pat1, npattern, int) ;
 ZALLOC(pat2, npattern, int) ;
 ZALLOC(wp, npattern, double) ;

 spat(pattern[0], pat1) ;
 spat(pattern[1], pat2) ;

  ZALLOC (ztop, nblocks, double);
  ZALLOC (zbot, nblocks, double);

  ZALLOC (znum, nblocks, double);
  ZALLOC (jmean, nblocks, double);
  ZALLOC (jwt, nblocks, double);

 for (k=0; k<ncols; ++k) { 
  s = bcols[k] ;  
  if (s<0) continue ; 
  cp = colprobs[k] ; 
  for (j=0; j<npattern; ++j) { 
   t = egindex[j] ;  
   wp[j] = cp[t] ;  
  }
  vmaxmin(wp, 4, NULL, &ymin) ; 
  if (ymin < -1) continue ; 
  y1 = abbaprob(wp, pat1, npattern) ; 
  y2 = abbaprob(wp, pat2, npattern) ; 
  ztop[s] += (y1 - y2) ;   
  z1 = (wp[0]-wp[1])*(wp[2]-wp[3]) ;  
  zsum += z1 ; 
  t = ranmod(1000) ; 
  if (t==-1) printf("zzcheck %d %d %12.6f   %12.6f\n", k, s, y1-y2, z1) ;
  znum[s] += 1 ; 
  if (fmode) zbot[s] += 1 ; 
  else zbot[s] += (y1 + y2) ; 
 }

  z2 = asum(znum, nblocks) ; 
//  printf("zzsum: %9.3f\n", z2) ; 

  ytop = asum (ztop, nblocks);
  ybot = asum (zbot, nblocks);
  ybot += 1.0e-20;
  ymean = ytop / ybot;
//   printf("zzmean: %12.6f %12.6f\n", ymean, zsum/z2) ; 
  for (s = 0; s < nblocks; ++s) {
    yt = ytop - ztop[s];
    yb = ybot - zbot[s];
    jmean[s] = yt / yb;
    jwt[s] = znum[s];
  }

  weightjack (pysc, pysig, ymean, jmean, jwt, nblocks);

  
  free (ztop);
  free (zbot);
  free (znum);
  free (jmean);
  free (jwt);
  free (wp);



 free(pat1) ; 
 free(pat2) ; 
 free(cp) ;  
 

 return ;


}
