#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <globals.h>
#include <mcmcpars.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "regsubs.h"  
#include "egsubs.h"  
#include "qpsubs.h" 

/** 
 See /home/np29/biologyx/Admix5.0/src/qpmdir for examples on how to run 
 Can't find .tex file for theory  
*/

#define Y  0 
#define E  1
#define I1 2
#define I2 3
#define A  4
#define O  5

#define WVERSION   "310" 
// was qpadmlin
// inbreed added

#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL ;
char *trashdir = "/var/tmp" ;

//int verbose = NO ;
//int plotmode = NO ;
int qtmode = NO ;

int colcalc = YES ;
int inbreed = NO ;

Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
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
// int fancynorm = YES ;
int gfromp = NO ;  // genetic distance from physical 
int msjack = NO ;
char *msjackname = NULL ;
int msjackweight = YES ;  // weighted jackknife

double plo = .001 ;
double phi = .999 ;
double pvhit = .001 ;
double pvjack = 1.0e-6 ;
double blgsize = 0.05 ;  // block size in Morgans */
double *chitot ;
int    *xpopsize ;
int  xmode = NO ;

int xchrom = -1 ;
int zchrom = -1 ;

char *genotypename = NULL ;
char  *snpname = NULL ;
char  *snpoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *outliername = NULL ;
char *outpop = NULL ; 
// list of outliers
// int outnum = -1 ;
double lambdascale = -1 ;  

int ldregress = 0 ;
double ldlimit = 9999.0 ;  /* default is infinity */
/* we only consider markers as in possible LD if gdis <= ldlimit */

char *outputname = NULL ;
char *weightname = NULL ;
FILE *ofile ;
char **eglist ;

void readcommands(int argc, char **argv) ;
void indiaestit(double *f2, double *f3, double *f4, int n)  ;
double estmix(double *z,  double *f3, int n) ;
double halfsamp(SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale)  ;

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  int  lsnplist, lindlist, numeg ;
  int i, j, k, k1, k2, k3, k4, kk; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double y1, y2, y, sig, tail, yy1, yy2 ;
  char ss[11] ;
  int *blstart, *blsize, nblocks ;
  int  xnblocks ; /* for xsnplist */
  int *bcols ;
  double maxgendis ;

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
  double *xmean, *xfancy ;
  double *divans, *divsd ; 
  double *hettop, *hetbot ; 
  int chrom,  numclear ;
  double gdis ;
  int outliter, *badlist, nbad ;
  double **zdata, *z1, *z2 ;
  int maxtag = -1 ;
  double **zz ; 
  double *pmean, *pnum, rscore[3], tscore[3], hscore[3], rrr[3], ww[3] ;
  int tpat[3][4] , rpat[3][4] ;
  int  *rawcol ; ;
  int a, b, c, d, col  ;
  double *qpscores ;
  double *hest, *hsig ;
  double mingenpos, maxgenpos ;
  int *qhit  ;  /* number of times pair is clade in quartet */
  int *qmiss ;  /* number of times pair migration event implied */
  int **qplist, numqp = 0, maxqp=10000 ;
  double *qpscore ;
  double jest, jsig, *wvar, *www ;


  double **zzn, **zzd ;
  int popx[4] ;
  double tn[4*5], td[4*4] ;
  double zzsig[5], zzest[5], zsc[5] ;
  double ymin ;
  double *zlin ;
  double *f2, *f2sig, *fst ;

  double *f3, *f4, *f3sig, *f4sig ;
  int ng2, ng3, ng4 ;


  readcommands(argc, argv) ;
  printf("## qpmix version: %s\n", WVERSION) ;
  if (parname == NULL) return 0 ;

  if (outputname != NULL)  openit(outputname, &ofile, "w") ;
  setinbreed(inbreed) ;

  if (xmode) xchrom = numchrom + 1 ;
  if (xchrom == (numchrom + 1)) noxdata = NO;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  k = getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

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

  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;
  printf("numsnps: %d\n", numsnps) ;

   ZALLOC(eglist, numindivs, char *) ; 
   numeg = loadlist(eglist, poplistname) ;
   seteglist(indivmarkers, numindivs, poplistname);

 for (i=0; i<numeg; i++) {  
  printf("%3d %s\n",i, eglist[i]) ;
 }

 printf("jackknife block size: %9.3f\n", blgsize) ;

  if (outpop == NULL) { 
   printf("setting outpop as target!\n") ;
   outpop = strdup(eglist[0]) ;
  }

  for (k=0; k<numeg; ++k)  {  
    t = strcmp(eglist[k],outpop) ;
    if (t==0) outnum = k ;
  }

  if (outnum < 0) { 
   eglist[numeg] = strdup(outpop) ;  
   outnum = numeg ;
  }


/* we also want outpop */

 nindiv=0 ;
 for (i=0; i<numindivs; i++) {
  indx = indivmarkers[i] ;
  t = strcmp(indx -> egroup, outpop) ; 
  if (t==0) indx -> affstatus = YES ;
  if(indx -> affstatus == YES) ++nindiv  ;
 }

   for (i=0; i<numsnps; i++)  {  
    cupt = snpmarkers[i] ; 
    if (cupt -> ignore) continue ;
    if (numvalidgtx(indivmarkers, cupt, YES) <= 1) { 
      printf("nodata: %20s\n", cupt -> ID) ;
      cupt -> ignore = YES ;
    }
   }

  ZALLOC(xindex, numindivs, int) ;
  ZALLOC(xindlist, numindivs, Indiv *) ;
  ZALLOC(rawcol, numindivs, int) ;  
  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;
  ZALLOC(xtypes, nrows, int) ;
  for (i=0; i<nrows; i++) {
   indx = xindlist[i] ;
   k = indxindex(eglist, numeg, indx -> egroup) ;
   xtypes[i] = k ;
   t = strcmp(indx -> egroup, outpop) ;
   if (t==0) xtypes[i] = outnum ;
  }


  printf("before setwt numsnps: %d\n", numsnps) ;
  setwt(snpmarkers, numsnps, indivmarkers, nrows, xindex, xtypes, outpop, eglist, numeg) ;
  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;  //  rid ignorable snps
  printf("setwt numsnps: %d\n", numsnps) ;
  if (numsnps ==  0) fatalx("no valid snps\n") ;

  setmgpos(snpmarkers, numsnps, &maxgendis) ;  
  if ((maxgendis < .00001) || (gfromp == YES)) setgfromp(snpmarkers, numsnps) ;

  nblocks = numblocks(snpmarkers, numsnps, blgsize) ;
  ZALLOC(blstart, nblocks, int) ;
  ZALLOC(blsize, nblocks, int) ;
  printf("number of blocks for moving block jackknife: %d\n", nblocks) ;

  ZALLOC(xsnplist, numsnps, SNP *) ;
  ZALLOC(tagnums, numsnps, int) ;

  if (popsizelimit > 0) {  
   setplimit(indivmarkers, numindivs, eglist, numeg, popsizelimit) ; 
  }

  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;


  printf("snps: %d  indivs: %d\n", ncols, nrows) ;
  setblocks(blstart, blsize, &xnblocks, xsnplist, ncols, blgsize)  ;
  printf("number of blocks for moving block jackknife: %d\n", xnblocks) ;

  zz = initarray_2Ddouble(numeg*numeg, ncols, 0.0) ;
  qplist = initarray_2Dint(numeg*numeg*numeg*numeg, 4, 0) ;
  ZALLOC (qpscore,numeg*numeg*numeg*numeg, double) ;

  ZALLOC(pmean, numeg, double) ; 
  ZALLOC(pnum, numeg, double) ;  
  ZALLOC(rawcol, nrows, int) ;  

  ng2 = numeg*numeg ;  
  ng3 = numeg*ng2   ;  
  ng4 = numeg*ng3   ;  
  ZALLOC(f2, ng2, double) ;
  ZALLOC(f2sig, ng2, double) ;
  ZALLOC(fst, ng2, double) ;
  ZALLOC(wvar, ng2, double) ;
  ZALLOC(www, ng2, double) ;

  y = dofstnum(fst, f2, f2sig, xsnplist, xindex, xtypes, nrows, ncols, numeg, nblocks) ;     
  if (lambdascale < 0.0) lambdascale = y ;
  printf("lambdascale: %9.3f\n", lambdascale) ;

  printf("f2:") ; 
  printnl() ;
  printmatw(f2, numeg, numeg, numeg) ;
  printnl() ;
  
  ZALLOC(f3, ng3, double) ;
  ZALLOC(f3sig, ng3, double) ;

  ZALLOC(zlin, numeg, double) ;
//  verbose = YES ;
  y = doadmlin(&jest, &jsig, zlin, wvar, 
    xsnplist, xindex, xtypes, 
    nrows, ncols, numeg, nblocks, lambdascale, indivmarkers) ;     

  printf("f2 error estimate: %12.6f s.err: %12.6f jest: %12.6f\n", y, jsig, jest) ;
  y = halfsamp(xsnplist, xindex, xtypes, nrows, ncols, numeg, nblocks, lambdascale) ;     
  printf("halfsamp sigmage: %9.3f\n", y) ;

  d = numeg-1 ;
  printf("coeffs: \n") ;
  printmatw(zlin, 1, d, d) ;
  printnl() ;

  vst(www, wvar, 1.0e6, ng2) ;
  printf("variance * 1000000: \n") ;
  printmatw(www, d, d, d) ;

  printnl() ;
  printnl() ;

  for (i=1; i<numeg; ++i) { 
   printf("%15s ", eglist[i]) ;
   j = i - 1; 
   printf("%9.3f  ", zlin[j]) ;
   y = sqrt(wvar[j*d+j]) ;
   printf("%12.6f", y) ;
   printnl() ;
  }


  return 0 ;
}
void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

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

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "outpop:", &outpop) ;
   getstring(ph, "output:", &outputname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "msjackname:", &msjackname) ;
   getint(ph, "msjackweight:", &msjackweight) ;
   getdbl(ph, "blgsize:", &blgsize) ;
   getdbl(ph, "lambdascale:", &lambdascale) ;
   getint(ph, "xmode:", &xmode) ;

   getint(ph, "noxdata:", &noxdata) ; 
   getint(ph, "colcalc:", &colcalc) ; 
   getint(ph, "inbreed:", &inbreed) ; 

   getint(ph, "nostatslim:", &nostatslim) ; 
   getint(ph, "popsizelimit:", &popsizelimit) ; 
   getint(ph, "gfromp:", &gfromp) ;  // gen dis from phys
// getint(ph, "fancynorm:", &fancynorm) ; 
// getint(ph, "usenorm:", &fancynorm) ;  /* changed 11/02/06 */

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}

double
halfsamp(SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale)     

{

   int t1, t2 ;
   int a, b, c ;
   int ng3, ng2 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j, d, t ; 
   double ya, yb, y, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot, wt ;
   double *w1, *w2, *w3 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   double xest, xsig, ynominal ;
   int bnum  ;

   double *f3, *f3sig  ;
   double *estmat, *zl, *rhs, errest ;
   double *vmean, **vjmean ;
   int *btype ;
   double *zl1, *zl2, *hest, ysig ; 
   double *xtop, *xbot ;
   
   
   ng2 = numeg*numeg ;
   ng3 = numeg*numeg*numeg ;

   ZALLOC(btype, nblocks, int)  ;  

// half sample  
   for (k=0; k<nblocks; ++k)  {  
    btype[k] = ranmod(2) + 1 ;
   }

   ZALLOC(zl1, numeg, double) ;
   ZALLOC(zl2, numeg, double) ;
   ZALLOC(hest, nblocks, double) ;

   ZALLOC(f3, ng3, double) ;
   ZALLOC(f3sig, ng3, double) ;
   ZALLOC(w1, ng3, double) ;
   ZALLOC(w2, ng3, double) ;
   ZALLOC(estmat, ng3, double) ;
   ZALLOC(w3, ng3, double) ;
   ZALLOC(gtop, ng3, double) ;
   ZALLOC(gbot, ng3, double) ;
   ZALLOC(wtop, ng3, double) ;
   ZALLOC(wbot, ng3, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   ZALLOC(xtop, nblocks, double) ;
   ZALLOC(xbot, nblocks, double) ;
   ZALLOC(rhs, numeg, double) ;
   btop  = initarray_2Ddouble(nblocks, ng3,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, ng3,  0.0) ;

   d = numeg - 1 ;

   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    if (cupt -> ignore) continue ;
    wt = cupt -> weight ;  
    if (wt <= 0.0) continue ;
    bnum = cupt -> tagnumber ; 
    if (bnum<0) continue ;
    ++wjack[bnum] ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;

    f3yyx(estmat,  cupt, xindex, xtypes, nrows, numeg, indivmarkers) ;
    vst(estmat, estmat, wt*scale, ng3) ;
    vvp(top, top, estmat, ng3) ;
    vsp(bot, bot, 1.0, ng3) ;

   }

    for (k=0; k<nblocks; k++) {  
     if (btype[k] == 2) continue ;
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, ng3) ;
     vvp(gbot, gbot, bot, ng3) ;
    }

   vsp(w2, gbot, 1.0e-10, ng3) ;
   vvd(f3, gtop, w2, ng3) ;



  for (i=0; i<numeg; ++i) { 
   printf("f3: base number %d:\n", i) ;
   printmatw(f3+i*numeg*numeg, numeg, numeg, numeg) ;
  }

  y1 = estmix(zl1+1, f3, numeg) ;
  printf("nominal error 1: %12.6f\n", y1) ;

  vzero(w2, ng3) ;
  for (a=0; a < d ; ++a) {
   for (b=0; b < d ; ++b) {
    for (c=0; c < d ; ++c) {
     y = dump3(f3, a, b, c, numeg) ;
     bump3(w2, a, b, c, d, y) ;
    }
   }
  }

  y2 = estmix(zl2+1, w2, d) ;
  printf("nominal error 2: %12.6f\n", y2) ;

  ytop = ybot = 0.0 ;
  t = 0 ;
  for (k=0; k<nblocks; k++) {  
   if (btype[k] == 1) {      
     wjack[k] = 0.0 ;
     djack[k] = 0.0 ; 
     continue ;
   }
   top = btop[k] ; 
   bot = bbot[k] ; 
   mulmat(rhs, top, zl1, numeg, numeg, 1) ;
   y1 = vdot(zl1, rhs, numeg)  ;
   mulmat(rhs, top, zl2, numeg, numeg, 1) ;
   y2 = vdot(zl2, rhs, numeg)  ;
   hest[t] = djack[k] = (y2 - y1)/bot[0] ;
   xtop[t] = y2-y1 ;
   xbot[t] = bot[0] ;
   wjack[t] = wjack[k] ;
   ++t ;
   ytop += (y2-y1) ; 
   ybot += bot[0] ;
  }

  mean = ytop/ybot ;

  for (j=0; j<t; ++j) { 
   y1 = ytop - xtop[j] ;
   y2 = ybot - xbot[j] ;
   djack[j] = y1/y2 ;
  }

  printf("half data size: %d\n", t) ;
  wjackest(&xest, &xsig, mean, djack, wjack, t) ;
  printf("mean: %12.6f  jmean:  %12.6f\n", mean, xest) ;
  ysig = xest/xsig ;              

    free(hest) ; 
    free(zl1) ; 
    free(zl2) ; 
    free(btype) ; 

    free(xtop) ; 
    free(xbot) ; 

    free(w1) ; 
    free(w2) ; 
    free(w3) ; 

    free(estmat) ;
    free(rhs) ;

    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(djack) ;
    free(wjack) ;
    free(f3) ;
    free(f3sig) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);

    return ysig ;

}

