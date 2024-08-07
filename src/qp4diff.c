
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

 Repeated pops (f2, f3, now supported
*/


#define WVERSION   "465" 

// overlap NO added
// allsnps: YES and instem: added
// firstf4mult added (Mark request) 
// print number of snps used 
// weightjack trapped (too few blocks) 

#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL ;
char *trashdir = "/var/tmp" ;
int xverbose ; 
int qtmode = NO ;
int colcalc = YES ;
int hires = NO ;
int fancyf4 = NO;

Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int isinit = NO ;
int overlap = YES ; // f4 must work for both halves
int allsnps = -99  ; // allsnps YES => overlap NO   
int markerscore = NO ;
int seed = 0 ;
int chisqmode = NO ;  // approx p-value better to use F-stat
int missingmode = NO ;
int dotpopsmode = YES ;
int noxdata = YES ;  /* default as pop structure dubious if Males and females */
int pcorrmode = NO ;
int pcpopsonly = YES ;
int nostatslim = 10 ;
int znval = -1 ;
int popsizelimit = -1 ;
int gfromp = NO ;  // genetic distance from physical 
int msjack = NO ;
char *msjackname = NULL ;
int msjackweight = YES ;  // weighted jackknife
int bankermode = NO ;  
int forceclade = NO ;
int numbanker = 0 ;   
int xchrom = -1 ;
// if bankermode  bankers MUST be in quartet  at most one type 1 in quartet

int jackweight = YES ; 
double jackquart = -1.0 ;

double plo = .001 ;
double phi = .999 ;
double pvhit = .001 ;
double pvjack = 1.0e-6 ;
double blgsize = 0.05 ;  // block size in Morgans */
double *chitot ;
int    *xpopsize ;

char *instem = NULL ; 
char  *indivname = NULL ;
char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *snpoutfilename = NULL ;
char *badsnpname = NULL ;
char *popfilename = NULL ;
char *outliername = NULL ;
char *blockname = NULL;
int inbreed  = NO ; 
double lambdascale ;  
int *f2ind, *ind2f, ng2, nh2 ;

int ldregress = 0 ;
double ldlimit = 9999.0 ;  /* default is infinity */
double firstf4mult = 1.0 ; 
/* we only consider markers as in possible LD if gdis <= ldlimit */

char *outputname = NULL ;
char *weightname = NULL ;
FILE *ofile ;

void readcommands(int argc, char **argv) ;
int readpopx(char *pname, char ***plists, int npops)  ;
int doq4diff(double *q4rat, double *q4ratsig, int ***counts, int *bcols, 
 int nrows, int ncols, int *xtop, int *xbot, int numeg, int nblocks) ;


int usage (char *prog, int exval); 

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  char **eglist ;
  int *nsamppops ;
  int  lsnplist, lindlist, numeg ;
  int i, j, k, k1, k2, k3, k4, kk; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double y1, y2, y, ysig, tail, yy1, yy2 ;
  char ss[11] ;
  int *blstart, *blsize, nblocks ;
  int  xnblocks ; /* for xsnplist */
  int *bcols ;
  double maxgendis ;
  int xind[4] ;
  int xtop[4], xbot[4] ;
  int xtop2[4], xbot2[4] ;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n1, n2 ;

  int nindiv = 0, e, f, lag=1  ;
  double xc[9], xd[4], xc2[9] ;
  int nignore, numrisks = 1 ;
  double  *xrow, *xpt ; 
  SNP **xsnplist  ;
  int *tagnums  ;
  Indiv **xindlist ;
  int *xindex, *xtypes ;
  int nrows, ncols, m, nc ;
  double zn, zvar ;
  int weightmode = NO ;
  double chisq, ynrows ;
  int *numhits, t ;  
  double *divans, *divsd ; 
  double *hettop, *hetbot ; 
  int chrom,  numclear ;
  double gdis ;
  int outliter, *badlist, nbad ;
  double **zdata, *z1, *z2 ;
  int maxtag = -1 ;
  double **zz ; 
  double *pmean, *pnum, rscore[3], dstat[3], hscore[3], rrr[3], ww[3] ;
  int tpat[3][4] , rpat[3][4], *rrtmp ;
  int  *rawcol ; ;
  int a, b, c, d, col  ;
  double *qpscores ;
  double *hest, *hsig ;
  double mingenpos, maxgenpos ;
  int *qhit  ;  /* number of times pair is clade in quartet */
  int *qmiss ;  /* number of times pair migration event implied */
  int **qplist, numqp = 0, maxqp=10000 ;
  double *qpscore ;
  int nsnps ; 


  double **dsctop, **dscbot ;
  double **abx, **bax ;
  int popx[4] ;
  double tn[4*5], td[4*4] ;
  double zzsig[5], zzest[5], zsc[5] ;
  double ymin ;
  double *f2, *f2sig, *fst ;

  double *f3, *f4, *f3sig, *f4sig ;
  int t1, t2 ;
  int ***counts, **ccc ; 
  char ***plists ; 
  char *px ;
  int nplist = 0, trun, nplisth ;
  double ymem ;
  int numchromp ; 


  readcommands(argc, argv) ;

  cputime(0) ;
  calcmem(0) ;

  printf("## qp4diff version: %s\n", WVERSION) ;

  numchromp = numchrom+1 ;
  if (parname == NULL) return 0 ;
  xverbose = verbose ;
  if (xverbose) printf("verbose set\n") ;
  if (allsnps == YES) overlap = NO ; 
  if (allsnps == NO) overlap = YES ; 
  if (fancyf4) printf("fancyf4 on\n") ;  
  if (fancyf4==NO) printf("fancyf4 off\n") ;  

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 

  y = fabs(firstf4mult-1.0) ; 
  if (y>.0001) printf("*** firstf4mult: %9.3f ***\n", firstf4mult) ;
  

  nostatslim = MAX(nostatslim, 3) ;
  setinbreed(inbreed) ;
  setfancyf4(fancyf4) ; 

  if (outputname != NULL)  openit(outputname, &ofile, "w") ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  setindm(indivmarkers) ;

  k = getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  for (i=0; i<numsnps; i++)  
  {  
    cupt = snpmarkers[i] ; 
    chrom = cupt -> chrom ;
    if ((xchrom>0) && (chrom != xchrom)) cupt -> ignore = YES ;
    if (chrom == 0) cupt -> ignore = YES ;
    if (chrom >= numchromp) cupt -> ignore = YES ;
  }

  nplist = numlines(popfilename) ;
  ZALLOC(plists, nplist, char **) ;
  num = readpopx(popfilename, plists, 8) ;
  nplist = num ;
  printf("nplist: %d\n", nplist) ;
  if (nplist == 0) return 0 ;


  ZALLOC(eglist, nplist*8, char *)  ;  
  numeg = 0 ;

  for (trun=0; trun<nplist; ++trun) {  
   for (k=0; k<8; ++k) { 
    px = plists[trun][k] ;
    t = strcmp(px, "NULL") ;
    if (t==0) continue ;
    t = indxindex(eglist,  numeg, px) ;
    if (t<0) {  
     eglist[numeg] = strdup(px) ;
     ++numeg ;
    }
   }
  }


   ZALLOC(nsamppops, numeg, int) ;

   for (i=0; i<numindivs; i++) {
    indx = indivmarkers[i] ;
    if (indx -> ignore) continue ;
    k = indxindex(eglist, numeg, indx -> egroup) ;
    if (k<0) { 
     indx -> ignore = YES ;
     continue ;
    }
    indx -> affstatus = YES ;
    ++nsamppops[k] ;
   }

 for (i=0; i<numeg; i++) {  
  t = nsamppops[i] ;
  if (t==0) fatalx("No samples: %s\n", eglist[i]) ;
  printf("%3d %20s %4d\n",i, eglist[i], t) ;
 }

 printf("jackknife block size: %9.3f\n", blgsize) ;


  ZALLOC(xindex, numindivs, int) ;
  ZALLOC(xindlist, numindivs, Indiv *) ;
  ZALLOC(rawcol, numindivs, int) ;  
  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
  if (nrows == 0) fatalx("badbug: no data\n") ;  

  ZALLOC(xtypes, nrows, int) ;
  for (i=0; i<nrows; i++) {
   indx = xindlist[i] ;
   k = indxindex(eglist, numeg, indx -> egroup) ;
   xtypes[i] = k ;
  }


   for (i=0; i<numsnps; i++)  {  
    cupt = snpmarkers[i] ; 
    cupt -> weight = 1.0 ;
   }
  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;  //  rid ignorable snps
  if (numsnps ==  0) fatalx("no valid snps\n") ;

  setmgpos(snpmarkers, numsnps, &maxgendis) ;  
  if ((maxgendis < .00001) || (gfromp == YES)) setgfromp(snpmarkers, numsnps) ;


  ZALLOC(xsnplist, numsnps, SNP *) ;
  ZALLOC(tagnums, numsnps, int) ;

  if (popsizelimit > 0) {  
   setplimit(indivmarkers, numindivs, eglist, numeg, popsizelimit) ; 
  }

  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;

  printf("snps: %d  indivs: %d\n", ncols, nrows) ;
  nblocks = setblocksz (&blstart, &blsize, xsnplist, ncols, blgsize, blockname) ;

// loads tagnumber
  printf ("number of blocks for block jackknife: %d\n", nblocks);
  xnblocks = nblocks += 10 ; 

  ZALLOC(counts, ncols, int **) ;          
  for (k=0; k<ncols; ++k) { 
   counts[k] = initarray_2Dint( numeg, 2, 0) ;
  }

  countpops(counts, xsnplist, xindex, xtypes, nrows, ncols) ;

  
  ZALLOC(bcols, ncols, int) ;  // blocks for columns -1 => unused
  ivclear(bcols, -1, ncols) ;
  for (k=0; k<ncols; k++)  { 
   cupt = xsnplist[k] ;
   bcols[k] = cupt -> tagnumber ;
  }
   for (a=0; a<nplist; ++a) { 
    for (i=0; i<4; ++i)  { 
     xtop[i] = indxindex(eglist, numeg, plists[a][i]) ;
     xbot[i] = indxindex(eglist, numeg, plists[a][i+4]) ;
    }
    nsnps = doq4diff(&y, &ysig, counts, bcols, 
     nrows, ncols, xtop, xbot, numeg, nblocks) ; 
    if (nsnps < -100) {
     y=0.0; ysig = 1.0 ; printf("badoctet: ") ;
    }
    else printf("result: ") ;
    for (t=0; t<8 ; ++t) { 
     printf("%10s ", plists[a][t]) ;
     if (t==3) printf(" : ") ;
    }
    printf("%12.6f %12.6f  %9.3f", y, ysig, y/ysig) ;
    printf(" %7d", nsnps) ; 
    printnl() ;
   }



  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qp4diff: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

  return 0 ;
}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  if (argc == 1) { usage(basename(argv[0]), 1); }
  while ((i = getopt (argc, argv, "p:vVh")) != -1) {

    switch (i)
      {

     case 'h':
      usage(basename(argv[0]), 0);
      break;

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   if (parname==NULL) { 
    fprintf(stderr, "no parameters\n") ;
    return ;
   }

   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring (ph, "instem:", &instem);
   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "popfilename:", &popfilename) ;
   getstring(ph, "output:", &outputname) ;
   getstring(ph, "outputname:", &outputname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getdbl(ph, "blgsize:", &blgsize) ;

   getint(ph, "overlap:", &overlap) ; 
   getint(ph, "allsnps:", &allsnps) ; 
   getint(ph, "nostatslim:", &nostatslim) ; 
   getint(ph, "popsizelimit:", &popsizelimit) ; 
   getint(ph, "gfromp:", &gfromp) ;  // gen dis from phys
   getint(ph, "fancyf4:", &fancyf4) ;
   getint(ph, "numchrom:", &numchrom) ;
   getdbl(ph, "firstf4mult:", &firstf4mult) ;

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}

double
getrs(double *f4, double *f4sig, double *f2, double *f2sig, 
  int a, int b, int c, int d, int numeg, double *rho) 
{

        double y1, y2, ya, yb ;

        y1 = dump4(f4, a, b, c, d, numeg) ;
        y2 = dump4(f4sig, a, b, c, d, numeg) ;
        ya = f2[a*numeg+b] ;                 
        yb = f2[c*numeg+d] ;                 

        *rho = y1/sqrt(ya*yb) ;

        return y1/y2 ;

}
int islegal(int *xind, int *types,  int mode, int n) 
{
// at least 1 type 1 
// check all bankers set.  OR at least 2 bankers 
  int cc[3] ;
  int zz[MAXPOPS] ;
  int i, x, k, t ;

  
  if (mode == NO) return YES ;
  ivzero(zz, n) ;
  ivzero(cc, 3) ;
  for (i=0; i<4; i++) {  
   t = xind[i] ;
   x = types[t]  ;
   ++cc[x] ;
   ++zz[t] ;
  }

/**
  printf("zzq1 ") ;  printimat(xind, 1, 4) ;
  printimat(cc, 1, 4) ;
  printimat(zz, 1, n) ;
  printimat(types, 1, n) ;
*/

  if (cc[1] == 0) return NO ; 
  if (cc[2] >= 2) return YES ;
//  if (cc[1] == 2) return NO ; 
  for (t=0; t<n; ++t) { 
   if (types[t] != 2) continue ;
   if (zz[t] == 0) return NO ;
  }
  return YES ;
}
int isclade(int *rr, int *zz) 
// is a clade banker only? 
{

 int a, b ;

 a = rr[0] ; b = rr[1] ;
 if ((zz[a] == 2) && (zz[b] == 2)) return  YES ;

 a = rr[2] ; b = rr[3] ;
 if ((zz[a] == 2) && (zz[b] == 2)) return YES ;

 return NO ;

}
void setabx(double **abx, double **bax, int ***counts, int ncols, int numeg) 
{

   
   int i, j, k, t1, t2, a, h ;
   double y1, y2 ; 
   int **ccc ;
   
   clear2D(&abx, nh2, ncols, 0.0) ;
   clear2D(&bax, nh2, ncols, 0.0) ;    

  for (k=0; k<ncols; ++k) {

   ccc = counts[k] ;

   for (i=0; i< numeg; i++) {  
    for (j=i+1; j< numeg; j++) {  

      t1 = intsum(ccc[i], 2) ;
      t2 = intsum(ccc[j], 2) ;
      if (t1==0) continue ; 
      if (t2==0) continue ; 

      a = ccc[i][0]  ;
      y1 = (double) a / (double) t1 ; 
      a = ccc[j][0]  ;
      y2 = (double) a / (double) t2 ; 
      h = f2ind[i*numeg+j] ;
      
      abx[h][k] = y1 * (1-y2) ;           
      bax[h][k] = y2 * (1-y1) ;           

    }
   }
  }
}



int readpopx(char *pname, char ***plists, int npops)  
// reads lists of n pops on a line
{
 FILE *fff ; 
 char line[MAXSTR+1], c ;
 char *spt[MAXFF], *sx  ;
 char **pp ;
 int nsplit,  t, num = 0  ;

 openit(pname, &fff, "r") ;
 line[MAXSTR] = '\0' ;
 while (fgets(line, MAXSTR, fff) != NULL)  {
   subcolon(line) ; 
   nsplit = splitup(line, spt, MAXFF) ;
   if (nsplit == 0) continue ;
   sx = spt[0] ; 
   if (sx[0] == '#')  {  
    freeup(spt, nsplit) ;
    continue ;
   }
   if (nsplit < npops) fatalx("length mismatch %s\n", line) ;
   ZALLOC(plists[num], npops+1, char *) ;
   plists[num][npops] = NULL ;
   pp = plists[num] ;
   for (t=0; t<npops; ++t) {  
    pp[t] = strdup(spt[t]) ;
   }
   ++num ; 
   freeup(spt, nsplit) ;
  }
  fclose(fff) ;
  return num ;
}
int
doq4diffx(double *q4diff, double *q4diffsig, int ***counts, int *bcols, 
 int nrows, int ncols, int *xtop1, int *xtop2, int numeg, int nblocks) 
// don't demand overlap
{

   int a, b, c, d ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   double top, bot, *djack, *wjack, gtop, gbot, *wbot, *wtop ;
   double gtop1, gtop2, gbot1, gbot2 ;
   double *btop1, *bbot1, wt ;
   double *btop2, *bbot2 ;
   double top1, bot1 ;
   double top2, bot2 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   double *w1, *w2, *ww, m1, m2 ;  
   int bnum, totnum  ;
   int ret, ret1, ret2, cntblocks=0 ;
   int num1, num2, num12 ;
   static int ncall = 0 ;
   
   if (nrows==0) fatalx("badbug\n") ;

   ZALLOC(wjack, nblocks, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(btop1, nblocks, double) ;
   ZALLOC(bbot1, nblocks, double) ;
   ZALLOC(btop2, nblocks, double) ;
   ZALLOC(bbot2, nblocks, double) ;
   ZALLOC(wtop, nblocks, double) ;
   ZALLOC(wbot, nblocks, double) ;

   
   ++ncall ;
   if (ncall == 1) printf("overlap NO mode (same as allsnps: YES)\n") ;
   num1 = num2 = num12 = 0 ; 

   totnum = 0 ;
   for (col=0; col<ncols;  ++col)  {
    bnum = bcols[col] ;  
    if (bnum<0) continue ;
    if (bnum>=nblocks) fatalx("logic bug\n") ;
    ret1 = ret = getf4(counts[col], xtop1, &y1) ;
    y1 *= firstf4mult ;
    if (ret == 2) fatalx("bad pop\n") ;
    if (ret>=0) { 
     btop1[bnum] += y1 ;
     bbot1[bnum] += 1 ; 
     ++wjack[bnum] ;
     ++num1 ; 
    }
    ret2 = ret = getf4(counts[col], xtop2, &y2) ;
    if (ret == 2) fatalx("bad pop\n") ;
    if (ret>=0) { 
     btop2[bnum] += y2 ;
     bbot2[bnum] += 1 ; 
     ++wjack[bnum] ;
     ++num2 ; 
     if ((ret1>=0) && (ret2>=0)) ++num12 ; 
    }
   }

    gtop1 = asum(btop1, nblocks) ;
    gbot1 = asum(bbot1, nblocks) ;
    gtop2 = asum(btop2, nblocks) ;
    gbot2 = asum(bbot2, nblocks) ;
    mean = (gtop1/gbot1) - (gtop2/gbot2) ;

    if (xverbose) { 
      printf("mean1: %12.6f\n", gtop1/gbot1) ;
      printf("mean2: %12.6f\n", gtop2/gbot2) ;
    }

    for (k=0; k<nblocks; k++) {  
     if (wjack[k] > 0.5) ++ cntblocks ; 
     top = btop1[k] ; 
     bot = bbot1[k] ;
     wtop[k] = gtop1-top ; 
     wbot[k] = gbot1-bot ;
     wbot[k] += 1.0e-10 ;
     djack[k] = wtop[k]/wbot[k] ;  // delete-block estimate
     top = btop2[k] ; 
     bot = bbot2[k] ;
     wtop[k] = gtop2-top ; 
     wbot[k] = gbot2-bot ;
     wbot[k] += 1.0e-10 ;
     djack[k] -= wtop[k]/wbot[k] ;  // delete-block estimate
    }
    if (cntblocks <=1) return -999 ;
      
    weightjack(&jest, &jsig, mean, djack,  wjack, nblocks) ;

    *q4diff = jest ;
    *q4diffsig = jsig ;

    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;

    free(btop1) ;
    free(bbot1) ;
    free(btop2) ;
    free(bbot2) ;

//  printf("zznum %d %d %d\n", num1, num2, num12) ;
    return num1 + num2 - num12 ;  // simple inclusion exclusion 
}

int
doq4diff(double *q4diff, double *q4diffsig, int ***counts, int *bcols, 
 int nrows, int ncols, int *xtop, int *xbot, int numeg, int nblocks) 
{

   int a, b, c, d ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   double top, bot, *djack, *wjack, gtop, gbot, *wbot, *wtop ;
   double *btop, *bbot, wt ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   double *w1, *w2, *ww, m1, m2 ;  
   int bnum, totnum  ;
   int ret, cntblocks=0 ;
   
   if (nrows==0) fatalx("badbug\n") ;
   if (overlap == NO) {
     return doq4diffx(q4diff, q4diffsig, counts, bcols, 
      nrows,  ncols, xtop, xbot,  numeg,  nblocks)  ;
   }

   ZALLOC(w1, nblocks, double) ;
   ZALLOC(w2, nblocks, double) ;
   ZALLOC(ww, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(btop, nblocks, double) ;
   ZALLOC(bbot, nblocks, double) ;
   ZALLOC(wtop, nblocks, double) ;
   ZALLOC(wbot, nblocks, double) ;

   totnum = 0 ;
   for (col=0; col<ncols;  ++col)  {
    ret = getf4(counts[col], xtop, &y1) ;

    y1 *= firstf4mult ; 
    if (ret < 0) continue ;
    if (ret == 2) fatalx("bad pop\n") ;

    ret = getf4(counts[col], xbot, &y2) ;
    if (ret < 0) continue ;
    if (ret == 2) fatalx("bad pop\n") ;

    bnum = bcols[col] ;  
    if (bnum<0) continue ;
    if (bnum>=nblocks) fatalx("logic bug\n") ;

      btop[bnum] += (y1-y2) ;
      bbot[bnum] += 1.0  ;
      ++wjack[bnum] ;
      ++totnum  ;
      
    }

    gtop = asum(btop, nblocks) ;
    gbot = asum(bbot, nblocks) ;

    mean = gtop/gbot ;

    for (k=0; k<nblocks; k++) {  
     if (wjack[k] >= 0.5) ++cntblocks ; 
     top = btop[k] ; 
     bot = bbot[k] ;
     wtop[k] = gtop-top ; 
     wbot[k] = gbot-bot ;
     wbot[k] += 1.0e-10 ;
     djack[k] = wtop[k]/wbot[k] ;  // delete-block estimate
    }

   if (cntblocks <= 1) return -999 ; 
      
    weightjack(&jest, &jsig, mean, djack,  wjack, nblocks) ;

    *q4diff = jest ;
    *q4diffsig = jsig ;

    free(w1) ; 
    free(w2) ; 
    free(ww) ; 

    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;

    free(btop) ;
    free(bbot) ;

    return totnum ;

}

void
doq4ratdiff(double *q4rat, double *q4ratsig, int ***counts, int *bcols, 
 int nrows, int ncols, int *xtop, int *xbot, int *xtop2, int *xbot2, int numeg, int nblocks) 
// calculates mean+Z of (m1 +m2 -1)
{

   int a, b, c, d ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   double top, bot, *djack, *wjack, gtop, gbot, gtop2, gbot2, *wbot, *wtop; 
   double *btop, *bbot, wt ;
   double *btop2, *bbot2 ;
   double ytop, ybot, ytop2, ybot2 ;
   double y1, y2, yscal ;
   double *w1, *w2, *ww, m1, m2 ;  
   int bnum, totnum  ;
   int ret ;
   
   if (nrows==0) fatalx("badbug\n") ;

   ZALLOC(w1, nblocks, double) ;
   ZALLOC(w2, nblocks, double) ;
   ZALLOC(ww, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(btop, nblocks, double) ;
   ZALLOC(bbot, nblocks, double) ;
   ZALLOC(btop2, nblocks, double) ;
   ZALLOC(bbot2, nblocks, double) ;
   ZALLOC(wtop, nblocks, double) ;
   ZALLOC(wbot, nblocks, double) ;

   totnum = 0 ;
   for (col=0; col<ncols;  ++col)  {
    ret = getf4(counts[col], xtop, &ytop) ;
    if (ret < 0) continue ;
    if (ret == 2) fatalx("bad pop in numerator\n") ;
    ret = getf4(counts[col], xbot, &ybot) ;
    if (ret < 0) continue ;

    ret = getf4(counts[col], xtop2, &ytop2) ;
    if (ret < 0) continue ;
    if (ret == 2) fatalx("bad pop in numerator\n") ;
    ret = getf4(counts[col], xbot2, &ybot2) ;
    if (ret < 0) continue ;

    bnum = bcols[col] ;  
    if (bnum<0) continue ;
    if (bnum>=nblocks) fatalx("logic bug\n") ;

      btop[bnum] += ytop ;
      bbot[bnum] += ybot ;
      btop2[bnum] += ytop2 ;
      bbot2[bnum] += ybot2 ;
      ++wjack[bnum] ;
      ++totnum  ;
      
    }

    gtop = asum(btop, nblocks) ;
    gbot = asum(bbot, nblocks) ;

    gtop2 = asum(btop2, nblocks) ;
    gbot2 = asum(bbot2, nblocks) ;

    mean = gtop/gbot ;
    mean += gtop2/gbot2 ;
    mean -= 1.0 ;

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     wtop[k] = gtop-top ; 
     wbot[k] = gbot-bot ;
     wbot[k] += 1.0e-10 ;
     djack[k] = wtop[k]/wbot[k] ;  
     top = btop2[k] ; 
     bot = bbot2[k] ;
     wtop[k] = gtop2-top ; 
     wbot[k] = gbot2-bot ;
     wbot[k] += 1.0e-10 ;
     djack[k] += wtop[k]/wbot[k] ;  // delete-block estimate
     djack[k] -= 1.0 ;
    }
      
    weightjack(&jest, &jsig, mean, djack,  wjack, nblocks) ;

    *q4rat = jest ;
    *q4ratsig = jsig ;

    free(w1) ; 
    free(w2) ; 
    free(ww) ; 

    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;

    free(btop) ;
    free(bbot) ;
    free(btop2) ;
    free(bbot2) ;

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

