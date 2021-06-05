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
#include "eigsubs.h"

#define Y  0
#define E  1
#define A  2

//  (YRI, CEU, Papua, .... )               


#define WVERSION   "651"
// popsizelimit
// dzeromode.  But this is a bad idea.  Must include monomorphic snps if we are to get unbiasedness
// snpdetailsname added
// main program to do f3 test
// nochrom: added
// count of number of non-mono snps added
// missing pop now just error message
// inbreed bug in f3 score snps dropped unnecessarily
// code cleanup (graph stuff deleted) and numchrom
// outgroupmode  -- no denominator
// count of non mono snps bugfix
// loadaa f3scz rewritten  

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *trashdir = "/var/tmp";
int details = NO;
int qtmode = NO;
char *f3name = NULL;
double jackquart = -1.0;

Indiv **indivmarkers;
SNP **snpmarkers;
int numsnps, numindivs;
int seed = 0;
int missingmode = NO;
int noxdata = YES;		/* default as pop structure dubious if Males and females */
int nostatslim = 10;
int znval = -1;
int popsizelimit = -1;
int gfromp = NO;		// genetic distance from physical 
int locount = -1, hicount = 9999999;
int jackweight = YES;
int pubjack = NO;
int outgroupmode = NO;

double blgsize = 0.05;		// block size in Morgans */ double *chitot ;
int xchrom = -1;
int zchrom = -1;
int *xpopsize;
int dzeromode = YES;

int isinit = NO;
double f2weight = 1.0;		// lsqmode only

char *instem = NULL ; 
char *indivname = NULL;
char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *badsnpname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
char *snpdetailsname = NULL;
char *blockname = NULL;
FILE *fff;


char *dumpname = NULL;
char *loadname = NULL;

char *outpop = NULL;
// list of outliers
char *basepop = NULL;
int basenum = -1;
double baseval = 0.0;

// outnum used for weights 
// basenum for f3 status  


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
int inbreed = NO;
int xnumeg;


char *outputname = NULL;
char *weightname = NULL;
FILE *ofile;
char **eglist;
char **egshort;
char **enames;
double zthresh = 3.0;
double f2diag = 0.0;

void readcommands (int argc, char **argv);
char *getshort (char *ss, int n);

void dopop3 (char **eglist, SNP ** xsnplist, int ncols, int nblocks);
int readpopx (char *pname, char ***plists, int n);
int dof3score (double *f3score, double *f3scoresig, SNP ** xsnplist,
	       int *xindex, int *xtypes, int nrows, int ncols, int numeg,
	       int nblocks);
void dof2score (double *f3score, double *f3scoresig, SNP ** xsnplist,
		int *xindex, int *xtypes, int nrows, int ncols, int numeg,
		int nblocks);
void estjackq (double *pjest, double *pjsig, double *btop, double *bbot,
	       double *wjack, int nblocks);

int usage (char *prog, int exval); 
int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int **snppos;
  int *snpindx;
  char **snpnamelist, **indnamelist;
  int lsnplist, lindlist;
  int i, j, k, k1, k2, k3, k4, kk;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;
  double y1, y2, y, sig, tail, yy1, yy2;
  char ss[11];
  int *blstart, *blsize, nblocks;
  int xnblocks;			/* for xsnplist */
  int *bcols;
  int **subsets;
  double maxgendis;
  char **eglist;
  int numeg;

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
  int ncols, nedge, m, nc;
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

  char ***qlist, ***plists, *sx;
  int nplist=0, trun, nqlist=0;
  double ymem; 

  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;

  printf ("## qp3Pop version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;
  if (xchrom == (numchrom + 1)) noxdata = NO;
  setjquart (pubjack, jackweight, jackquart);
  setinbreed (inbreed);

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
    if ((noxdata) && (chrom == (numchrom + 1)))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > (numchrom + 1))
      cupt->ignore = YES;
    if (chrom == zchrom)
      cupt->ignore = YES;
  }

  if (popfilename == NULL) fatalx("no popfilename!\n ") ;
    nplist= nqlist = numlines (popfilename);
    ZALLOC (qlist, 3, char **);
    for (k = 0; k < 3; ++k) {
      ZALLOC (qlist[k], nqlist, char *);
    }
    nqlist =
      getnamesstripcolon (&qlist, nqlist, 3, popfilename, locount, hicount);
    numeg = 0;
    printf ("number of triples %d\n", nqlist);

    if (nqlist==0) { 
     printf("no triples!\n") ; 
     return 0 ;
    }

    fflush (stdout);
    nplist = nqlist ; 
    ZALLOC (plists, nplist, char **);
     for (j = 0; j < nplist; ++j) {
      ZALLOC(plists[j], 3, char *) ; 
     }
    for (k = 0; k < 3; ++k) {
      for (j = 0; j < nplist; ++j) {
        sx = qlist[k][j];
        plists[j][k] = strdup(sx) ;
      }
    }

    for (k = 0; k < 3; ++k) {
      freeup (qlist[k], nqlist);
    }

  printf ("nplist: %d\n", nplist);
  if (nplist == 0)
    return 1;


  ZALLOC (eglist, nplist * 3, char *);
  numeg = 0;
  for (trun = 0; trun < nplist; ++trun) {
    for (k = 0; k < 3; ++k) {
      t = indxindex (eglist, numeg, plists[trun][k]);
      if (t < 0) {
	eglist[numeg] = strdup (plists[trun][k]);
	++numeg;
      }
    }
  }

  if (popsizelimit > 0) {
    setplimit (indivmarkers, numindivs, eglist, numeg, popsizelimit);
  }

  ZALLOC (xsnplist, numsnps, SNP *);

  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);

  printf ("snps: %d\n", ncols);
  nblocks = setblocksz (&blstart, &blsize, xsnplist, ncols, blgsize, blockname) ;

// loads tagnumber
  printf ("number of blocks for block jackknife: %d\n", nblocks);
  xnblocks = nblocks += 10 ; 

  for (trun = 0; trun < nplist; ++trun) {
    dopop3 (plists[trun], xsnplist, ncols, nblocks);
  }

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qp3Pop: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;

}

void
dopop3 (char **eglist, SNP ** xsnplist, int ncols, int nblocks)
{
  Indiv **xindlist;
  Indiv *indx;
  int *xindex, *xtypes;
  int nrows;
  int t, k, i, trun;
  int numeg = 3;
  double f3score, f3scoresig;
  double f2score, f2scoresig, y, y1, y2, p, q;
  int nsnp = 0;
  static int ncall = 0;


  ++ncall;
  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);

  setstatusv (indivmarkers, numindivs, NULL, NO);
  setstatuslist (indivmarkers, numindivs, eglist, numeg);

  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);

  if (nrows == 0) {
    for (i = 0; i < numeg; ++i) {
      printf ("pop:  %s\n", eglist[i]);
    }
    fatalx ("fatal error  (probably missing pop)\n");
  }

  ZALLOC (xtypes, nrows, int);

  if (ncall == 1) {
    printf ("%9s ", "");
    printf ("%20s ", "Source 1");
    printf ("%20s ", "Source 2");
    printf ("%20s ", "Target");

    printf (" %12s", "f_3");
    printf ("   %12s", "std. err");
    printf ("   %9s", "Z");
    printf (" %7s", "SNPs");
    printnl ();
  }


  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (eglist, numeg, indx->egroup);
    xtypes[i] = k;
  }

  nsnp = dof3score (&f3score, &f3scoresig, xsnplist, xindex, xtypes,
		    nrows, ncols, numeg, nblocks);

  if (nsnp < 0)
    printf ("%9s ", "no data");
  if (nsnp >= 0)
    printf ("%9s ", "result: ");

  for (t = 0; t < numeg; ++t) {
    printf ("%20s ", eglist[t]);
  }

  printf (" %12.6f", f3score);
  printf ("   %12.6f", f3scoresig);
  printf ("   %9.3f", f3score / f3scoresig);

/**
  if (f3score < 0.0) { 
   dof2score(&f2score, &f2scoresig, xsnplist, xindex, xtypes, 
    nrows, ncols, numeg, nblocks)  ;
// estimate f_2 using type 2 as outgroup 
   printf("   %9.3f", f2score) ;
   f2score = MAX(1.0e-6, f2score) ;
   
   y = -f3score/f2score ;
   y1 = 1.0-4.0*y ; 
   if (y1 <= 0.0)  printf(" ??") ;
   else {  
    y2 = sqrt(y1) ;
    p  = 0.5*(1.0+y2) ;
    q  = 0.5*(1.0-y2) ;
    printf("   %9.3f", p) ;
    printf(" %9.3f", q) ;
   }
  }
  if (f3score >= 0.0) printf("  %9s  %9s %9s", "-", "-", "-") ;
*/
  printf (" %7d", nsnp);
  printnl ();
  fflush (stdout);

  free (xtypes);
  free (xindex);
  free (xindlist);

  return;

}

void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n, t;

  if (argc == 1) { usage(basename(argv[0]), 1); }
  while ((i = getopt (argc, argv, "f:b:p:s:L:H:vVh")) != -1) {

    switch (i) {

    case 'h':
     usage(basename(argv[0]), 0);
      break;

    case 'p':
      parname = strdup (optarg);
      break;

    case 'L':
      locount = atoi (optarg);
      break;

    case 'H':
      hicount = atoi (optarg);
      break;

    case 'f':
      snpdetailsname = strdup (optarg);
      break;

    case 's':
      seed = atoi (optarg);
      break;

    case 'b':
      baseval = atof (optarg);
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
  getstring (ph, "snpdetailsname:", &snpdetailsname);
  getstring (ph, "outpop:", &outpop);
  getstring (ph, "output:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "popfilename:", &popfilename);
  getstring (ph, "f3log:", &f3name);
  getdbl (ph, "blgsize:", &blgsize);
  getint (ph, "numchrom:", &numchrom);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "outgroupmode:", &outgroupmode);

  getint (ph, "noxdata:", &noxdata);
  t = -1;
  getint (ph, "xdata:", &t);
  if (t >= 0)
    noxdata = 1 - t;
  getint (ph, "chrom:", &xchrom);
  getint (ph, "nochrom:", &zchrom);

  getint (ph, "nostatslim:", &nostatslim);
  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "seed:", &seed);
  getint (ph, "details:", &details);
  getdbl (ph, "baseval:", &baseval);
  getint (ph, "jackweight:", &jackweight);
  getint (ph, "pubjack:", &pubjack);
  getstring (ph, "dumpname:", &dumpname);
  getstring (ph, "loadname:", &loadname);
  getstring (ph, "blockname:", &blockname);
  getdbl (ph, "jackquart:", &jackquart);

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
    plists[num][npops] = NULL;
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
dof2score (double *f2score, double *f2scoresig, SNP ** xsnplist, int *xindex,
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


  setjquart (NO, YES, -1.0);

  totnum = 0;
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

    btop[bnum] += ytop;
    bbot[bnum] += ybot;
    ++wjack[bnum];
    ++totnum;

  }

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
dof3score (double *f3score, double *f3scoresig, SNP ** xsnplist, int *xindex,
	   int *xtypes, int nrows, int ncols, int numeg, int nblocks)
{
// returns number of snps

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
  FILE *fff;
  double xn[3], xmean[3], xh[3];
  int ret ;  

  *f3score = 0;
  *f3scoresig = -1;

  if (snpdetailsname != NULL)
    openit (snpdetailsname, &fff, "w");

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


  setjquart (pubjack, jackweight, jackquart);

  totnum = 0;
  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    if (cupt->ignore)
      continue;

    bnum = cupt->tagnumber;
    if (bnum < 0)
      continue;
    if (bnum >= nblocks)
      fatalx ("logic bug\n");

    ret = f3sc (&ytop, &ybot, cupt, indivmarkers, xindex, xtypes, nrows, 2, 0, 1);
// -1 => no good, 0 => mon 
    if (isnan (ytop))
      fatalx ("zznan\n");
    if (ybot < -0.5)
      continue;
    if (ret<0) continue ;
    if (outgroupmode) ybot = 0.001 ;   
//  if ((dzeromode == NO) && (ybot <= 0.001)) continue;
    if (snpdetailsname != NULL) {

      finfo (xn + 0, xmean + 0, xh + 0, 0);
      finfo (xn + 1, xmean + 1, xh + 1, 1);
      finfo (xn + 2, xmean + 2, xh + 2, 2);	// Sardin Karit CEU

      fprintf (fff, "%20s ", cupt->ID);
      for (a = 0; a <= 2; ++a) {
	fprintf (fff, "  %6.0f", xn[a]);
	fprintf (fff, " %9.3f", xmean[a]);
	fprintf (fff, " %9.3f", xh[a]);
      }
      fprintf (fff, "  %9.3f ", ytop);
      //   fprintf(fff, "%9.3f ", ybot/2.0) ;         
      fprintf (fff, " %c ", cupt->alleles[0]);
      fprintf (fff, " %c", cupt->alleles[1]);
      fprintf (fff, " %d", ret) ;
      fprintf (fff, "\n");


    }

    if (outgroupmode) ybot = 0.001 ;   

    btop[bnum] += ytop;
    bbot[bnum] += ybot;
    ++wjack[bnum];
    if (ret>0) ++totnum;			// monomorphic snps not counted

  }

  if (totnum <= 1)
    return -1;
/**
   printf("totnum: %d ", totnum) ;
   printmat(wjack, 1, 10) ;
*/



  gtop = asum (btop, nblocks);
  gbot = asum (bbot, nblocks);

  *f3score = mean = gtop / gbot;

  for (k = 0; k < nblocks; k++) {
    top = btop[k];
    bot = bbot[k];
    wtop[k] = gtop - top;
    wbot[k] = gbot - bot;
    wbot[k] += 1.0e-10;
    djack[k] = wtop[k] / wbot[k];	// delete-block estimate
  }

  estjackq (&jest, &jsig, btop, bbot, wjack, nblocks);
  *f3score = jest;
  *f3scoresig = jsig;

  free (w1);
  free (w2);
  free (ww);

  free (wtop);
  free (wbot);
  free (djack);
  free (wjack);

  free (btop);
  free (bbot);
  if (snpdetailsname != NULL)
    fclose (fff);
  return totnum;

}


int usage (char *prog, int exval)
{

  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -f <nam>    ... use <nam> as snp details name.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -L <index>  ... locount -n popfilename.\n");
  (void)fprintf(stderr, "   -H <index>  ... hicount -n popfilename.\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");

  exit(exval);
};

