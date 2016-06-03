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
#include "eigsubs.h"

#define Y  0
#define E  1
#define A  2

//  (YRI, CEU, Papua, .... )               


#define WVERSION   "120"
// popsizelimit
// dzeromode.  But this is a bad idea.  Must include monomorphic snps if we are to get unbiasedness
// snpdetailsname added
// main program to do f3 test
// nochrom: added

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *rootname = NULL;
char *trashdir = "/var/tmp";
int details = NO;
int qtmode = NO;
char *f3name = NULL;

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
int jackweight = YES;
int pubjack = NO;
double jackquart = -1.0;

int forcezmode = NO;
double blgsize = 0.05;		// block size in Morgans */ double *chitot ;
double diag = 0.0;
int bigiter = 100;
int startiter = 50;
int xchrom = -1;
int zchrom = -1;
int *xpopsize;
int dzeromode = YES;

int isinit = NO;
double f2weight = 1.0;		// lsqmode only

char *genotypename = NULL;
char *snpname = NULL;
char *snpoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *graphname = NULL;
char *graphoutname = NULL;
char *popfilename = NULL;
char *outliername = NULL;
char *snpdetailsname = NULL;
FILE *fff;

char *dumpname = NULL;
char *loadname = NULL;

char *outpop = NULL;
// list of outliers
char *basepop = NULL;
int basenum = -1;
double baseval = 0.0;
char **lines;

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

void dopop3out (char **fglist, SNP ** xsnplist, int ncols, char *line,
		char *outpop);
int readpopx (char *pname, char ***plists, int n);
void dof3score (double *f3score, double *f3scoresig, SNP ** xsnplist,
		int *xindex, int *xtypes, int nrows, int ncols, int numeg,
		int nblocks);
void dof2score (double *f3score, double *f3scoresig, SNP ** xsnplist,
		int *xindex, int *xtypes, int nrows, int ncols, int numeg,
		int nblocks);
void estjackq (double *pjest, double *pjsig, double *btop, double *bbot,
	       double *wjack, int nblocks);

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
  int nedge, m, nc;
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
  int ***counts;

  char ***plists;
  int nplist, trun;
  int nrows, ncols;

  readcommands (argc, argv);
  printf ("## qpBound version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;
  if (xchrom == 23)
    noxdata = NO;
  if (outpop == NULL)
    fatalx ("no outpop\n");
  setinbreed (inbreed);

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
    if ((noxdata) && (chrom == 23))
      cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > 23)
      cupt->ignore = YES;
    if (chrom == zchrom)
      cupt->ignore = YES;
  }

  nplist = numlines (popfilename);
  ZALLOC (plists, nplist, char **);
  ZALLOC (lines, nplist, char *);
  num = readpopx (popfilename, plists, 3);
  nplist = num;
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

  if (outputname != NULL)
    openit (outputname, &ofile, "w");
  outnum = 0;
  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);
  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);
  ZALLOC (xtypes, nrows, int);
  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (eglist, numeg, indx->egroup);
    xtypes[i] = k + 1;		// dangerous bend
    t = strcmp (indx->egroup, outpop);
    if (t == 0)
      xtypes[i] = outnum;
    else
      fatalx ("outpop bug\n");
  }
  ZALLOC (xsnplist, numsnps, SNP *);
  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);



/**
  ZALLOC(counts, ncols, int **) ;
  for (k=0; k<ncols; ++k) {
   counts[k] = initarray_2Dint( numeg, 2, 0) ;
  }
  countpops(counts, xsnplist, xindex, xtypes, nrows, ncols) ;
*/

  for (trun = 0; trun < nplist; ++trun) {
    dopop3out (plists[trun], xsnplist, ncols, lines[trun], outpop);
  }

  if (outputname != NULL)
    fclose (ofile);

  printf ("##end of qpBound\n");
  return 0;

}

void
dopop3out (char **fglist, SNP ** xsnplist, int ncols, char *line,
	   char *outpop)
{
  Indiv **xindlist;
  Indiv *indx;
  int *xindex, *xtypes;
  int nrows;
  int t, k, i, trun;
  double f3score, f3scoresig;
  double f2score, f2scoresig, y, y1, y2, p, q;
  char *eglist[4];
  int numeg = 4;
  double ytop, ybot, yxbot;
  double ztop, zbot;
  int col;
  SNP *cupt;
  double zztop[6], yytop[6];
  double u, s1, s2, atop, btop, alphabot, betabot, alphatop;
  double ya, yb, za, zb, yt;
  char obuff[1024], *sx;
  int nsnp = 0;


  copystrings (fglist, eglist, 3);
  eglist[3] = strdup (outpop);

  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);


  setstatusv (indivmarkers, numindivs, NULL, NO);
  setstatuslist (indivmarkers, numindivs, eglist, numeg);

  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);

  if (nrows == 0) {
    for (i = 0; i < numeg; ++i) {
      printf ("zz %s\n", eglist[i]);
    }
    fatalx ("fatal error (probably missing pop)\n");
  }

  ZALLOC (xtypes, nrows, int);


  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (eglist, numeg, indx->egroup);
    xtypes[i] = k;
  }

  ztop = zbot = 0.0;
  vzero (zztop, 6);
  for (col = 0; col < ncols; ++col) {
    cupt = xsnplist[col];
    if (cupt->ignore)
      continue;
    loadaa (cupt, xindex, xtypes, nrows, numeg);

    f3scz (&ytop, &ybot, cupt, indivmarkers, xindex, xtypes, nrows, 2, 0, 1);
    if (isnan (ytop))
      fatalx ("zznan\n");
    if (ybot < -0.5)
      continue;
    f3scz (&yytop[0], &yxbot, cupt, indivmarkers, xindex, xtypes, nrows, 3, 0,
	   1);
    if (yxbot < -0.5)
      continue;
    f3scz (&yytop[1], &yxbot, cupt, indivmarkers, xindex, xtypes, nrows, 3, 0,
	   2);
    if (yxbot < -0.5)
      continue;
    f3scz (&yytop[2], &yxbot, cupt, indivmarkers, xindex, xtypes, nrows, 3, 1,
	   2);
    if (yxbot < -0.5)
      continue;
    f2scz (&yytop[3], &yxbot, cupt, indivmarkers, xindex, xtypes, nrows, 3, 0,
	   3);
    if (yxbot < -0.5)
      continue;
    f2scz (&yytop[4], &yxbot, cupt, indivmarkers, xindex, xtypes, nrows, 3, 1,
	   3);
    if (yxbot < -0.5)
      continue;
    f2scz (&yytop[5], &yxbot, cupt, indivmarkers, xindex, xtypes, nrows, 3, 2,
	   3);
    if (yxbot < -0.5)
      continue;
    ztop += ytop;
    zbot += ybot;
    if ((ytop > 0) || (ybot > 0))
      ++nsnp;			// monomorphic snps not counted
    vvp (zztop, zztop, yytop, 6);
  }
//verbose = YES ; 
  ztop /= zbot;
  vst (zztop, zztop, 1.0 / zbot, 6);
  u = zztop[0];
  vsp (yytop, zztop, -u, 6);
  s1 = yytop[1];		/* alpha a */
  s2 = yytop[2];
  atop = yytop[3];
  btop = yytop[4];
  alphabot = s1 / atop;
  betabot = s2 / btop;
  alphatop = 1.0 - betabot;

  y1 = -ztop - s1;
  if (s2 > s1) {
    alphabot = MAX (alphabot, y1 / (s2 - s1));
  }
  if (s2 < s1) {
    alphatop = MIN (alphatop, y1 / (s2 - s1));
  }


  sx = obuff;
  sx += sprintf (sx, "%s", line);
//printf(" %12.6f", ztop) ;
  sx += sprintf (sx, " %9.3f", alphabot);
  sx += sprintf (sx, " %9.3f", alphatop);
/**
// next code is computing bounds on h (drift -> C) 
  za = alphatop; zb = 1.0-za ;
  ya = s1/za; yb = s2/zb; yt = -za*zb*(ya+yb) ; y1 = ztop - yt ; 
  za = alphabot; zb = 1.0-za ;
  ya = s1/za; yb = s2/zb; yt = -za*zb*(ya+yb) ; y2 = ztop - yt ; 
  sx += sprintf(sx, "     %9.3f %9.3f", y1, y2) ;
  sx += sprintf(sx, " %7d", nsnp) ;
*/
  printf ("%s", obuff);
  printnl ();
  if (verbose)
    printmatwl (yytop, 1, 6, 6);
  if (outputname != NULL) {
    fprintf (ofile, "%s\n", obuff);
    fflush (ofile);
  }

  free (xtypes);
  free (xindex);
  free (xindlist);
  freeup (eglist, 4);
  destroyaa ();

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

  while ((i = getopt (argc, argv, "f:b:p:g:s:o:vVx")) != -1) {

    switch (i) {

    case 'p':
      parname = strdup (optarg);
      break;

    case 'f':
      snpdetailsname = strdup (optarg);
      break;

    case 'g':
      graphname = strdup (optarg);
      break;

    case 'o':
      graphoutname = strdup (optarg);
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

    case 'x':
      doanalysis = NO;
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
  getstring (ph, "graphname:", &graphname);
  getstring (ph, "graphoutname:", &graphoutname);
  int numeg = 4;
  getstring (ph, "snpdetailsname:", &snpdetailsname);
  getstring (ph, "outpop:", &outpop);
  getstring (ph, "output:", &outputname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "popfilename:", &popfilename);
  getstring (ph, "f3log:", &f3name);
  getstring (ph, "root:", &rootname);
  getdbl (ph, "blgsize:", &blgsize);
  getdbl (ph, "diag:", &diag);
  getdbl (ph, "f2diag:", &f2diag);
  getdbl (ph, "minvar:", &minvar);
  getint (ph, "bigiter:", &bigiter);
  getint (ph, "inbreed:", &inbreed);
  getint (ph, "startiter:", &startiter);
  getint (ph, "fancynorm:", &fancynorm);

  getint (ph, "noxdata:", &noxdata);
  t = -1;
  getint (ph, "xdata:", &t);
  if (t >= 0)
    noxdata = 1 - t;
  getint (ph, "chrom:", &xchrom);
  getint (ph, "nochrom:", &zchrom);
  getint (ph, "doanalysis:", &doanalysis);
  getint (ph, "quartet:", &quartet);

  getint (ph, "nostatslim:", &nostatslim);
  getint (ph, "popsizelimit:", &popsizelimit);
  getint (ph, "gfromp:", &gfromp);	// gen dis from phys
  getint (ph, "seed:", &seed);
  getint (ph, "details:", &details);
  getint (ph, "forcezmode:", &forcezmode);
  getint (ph, "dzeromode:", &dzeromode);
  getdbl (ph, "baseval:", &baseval);
  getdbl (ph, "jackquart:", &jackquart);
  getint (ph, "jackweight:", &jackweight);
  getint (ph, "pubjack:", &pubjack);
  getstring (ph, "dumpname:", &dumpname);
  getstring (ph, "loadname:", &loadname);

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
// reads lists of npops pops on a line
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
    striptrail (line, '\n');
    lines[num] = strdup (line);
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

void
dof3score (double *f3score, double *f3scoresig, SNP ** xsnplist, int *xindex,
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
  FILE *fff;
  double xn[3], xmean[3], xh[3];

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

    f3sc (&ytop, &ybot, cupt, indivmarkers, xindex, xtypes, nrows, 2, 0, 1);
    if (isnan (ytop))
      fatalx ("zznan\n");
    if (ybot < -0.5)
      continue;
    if ((dzeromode == NO) && (ybot <= 0.001))
      continue;
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
      fprintf (fff, "%c", cupt->alleles[1]);
      fprintf (fff, "\n");


    }

    btop[bnum] += ytop;
    bbot[bnum] += ybot;
    ++wjack[bnum];
    ++totnum;

  }

  if (totnum <= 1)
    fatalx ("no data...\n");
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

}
