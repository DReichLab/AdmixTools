#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "regsubs.h"  
#include "egsubs.h"  
#include "qpsubs.h" 

#define Y  0 
#define E  1
#define A  2

//  (YRI, CEU, Papua, .... )               


#define WVERSION   "2051 " 
// lsqmode 
// ff3fit added
// reroot added
// deleted wtmin check
// fstdmode added (denominator a la fst)
// lambdascale can be a parameter
// format change (10.4f on graph
// inbreed
// phylipname added  
// graphdotname added
// map4x bug fixed
// zthresh param
// isfixed supported
// outliername added


#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL ;
char *rootname = NULL ;
char *trashdir = "/var/tmp" ;
int verbose = NO ;
int details = NO ;
int fstdetails = NO ;
int qtmode = NO ;
int fstdmode = NO ;
int plotmode = NO ;
int fancynorm = YES ;
int inbreed = NO ;
char *f3name = NULL ;

Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int seed = 0 ;
int missingmode = NO ;
int noxdata = YES ;  /* default as pop structure dubious if Males and females */
int doanalysis = YES ;  /* if no just print stats */
int nostatslim = 10 ;
int znval = -1 ;
int popsizelimit = -1 ;
int gfromp = NO ;  // genetic distance from physical 
int hires = NO ;

int forcezmode = NO ;
double blgsize = 0.05 ;  // block size in Morgans */ double *chitot ;
double diag = 0.0 ;
int  bigiter = 1 ;
int  xchrom  = -1 ;
int  xnochrom  = -1 ;
int    *xpopsize ;

int isinit = NO ;
int lsqmode = NO ;
double f2weight = 1.0 ;  // lsqmode only

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *snpoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *graphname = NULL ;
char *graphoutname = NULL ;
char *graphdotname = NULL ;
char *poplistname = NULL ;
char *outliername = NULL ;
char *phylipname = NULL ;

char *dumpname = NULL ;
char *loadname = NULL ;

char *outpop = NULL ; 
// list of outliers
int outnum = -1 ;
char *basepop = NULL ; 
int basenum = -1 ;  

// outnum used for weights 
// basenum for f3 status  


double lambdascale = -1.0 ;  
int *f2ind, *ind2f ; 
double *vest, *vvar, *vvinv, *xvvar, *xvvinv ;
double **vmix ; 
int *lmix, nmix ;
int nh2, numeg ;
int *ezero = NULL ; 
double wtmin = .0001 ;
double minvar = 0.0 ;  // minvalue for variance term
int quartet = NO ;
int xnumeg ;


char *outputname = NULL ;
char *weightname = NULL ;
FILE *ofile ;
char **eglist ;
char **egshort ;
char **enames ;
double zthresh = 3.0 ;
double f2diag  = 0.0 ;

void readcommands(int argc, char **argv) ;
void indiaestit(double *f2, double *f3, double *f4, int n)  ;
void sol2(double *co, double *rhs, double *ans) ;  
void mkww(double *f3, double *f2, int k, double *ww, int n)  ;
char *getshort(char *ss, int n) ;
void doff3(double *ff3, double *ff3var, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale)  ;
void map4x(double *aa, double *bb, int n2, int *indx) ;
void map4y(double *aa, double *bb, int n2, int *indx) ;
void getmv (int a, int b, int c, int d, double *mean, double *var, double *xest, double *xvar)  ; 
void bumpm(double *y, int a, int c, double val, double *xest) ;
void bumpw(double *w, int a, int c, double val) ;
void loadppwts(double *ppwts, double *pwts, int n) ;
double calcxx(double *xxans, double *qmat, double *ppwts, double *rhs, int nrow, int ncol)  ;
void setwww(double **tmix, double *www, int n)   ;
void getwww(double **tmix, double *www, int n)   ;
double scorit(double *www, int n, double *pfix, double *ans) ;
void nagopt(double *params, int n) ;
void printvals(double **tmix, double *edgelen, int nedge)  ;
void printfit(double *ww) ;
void initvmix(double *wwinit, int nwts, int numiter) ;
int iscanon(int a, int b, int c, int d) ;

void setrand(double *ww, int n)  ;
void setsimp(double *ww, int n)  ;
void forcez(double *cc, double *rr, int *ez, int n) ;
void calcfit(double *ff3fit, double *ww, int numeg)  ;
void dumppars(char *dumpname, double *www, int nwts, double *xxans, int nedge) ;
void dump1(FILE *dumpfile, double *ww, int n) ;
void loadpars(char *loadname, double *www, int nwts, double *xxans, int nedge) ;
void read1(FILE *loadfile, double *ww, int n) ;
double ff4val(double *ff3, int a, int b, int c, int d, int numeg) ;  
void print4(double *ff3, int a, int b, int c, int d, int numeg) ;
void printf3(char *sss, FILE *fff, double *ww, char **eglist, int n)  ;
int listsubset(int **x, int n, int k)  ;

int main(int argc, char **argv)
{

  char sss[MAXSTR] ;
  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  int  lsnplist, lindlist ;
  int i, j, k, k1, k2, k3, k4, kk ; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double y1, y2, y, sig, tail, yy1, yy2 ;
  char ss[11] ;
  int *blstart, *blsize, nblocks ;
  int  xnblocks ; /* for xsnplist */
  int *bcols ;
  int **subsets ;
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
  int nrows, ncols, nedge, m, nc ;
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
  int a, b, c, d, col, u, v, x  ;
  double *qpscores ;
  double *hest, *hsig ;
  double mingenpos, maxgenpos ;
  int *qhit  ;  /* number of times pair is clade in quartet */
  int *qmiss ;  /* number of times pair migration event implied */
  int **qplist, numqp = 0, maxqp=10000 ;
  int *popsizes ; 
  double *qpscore ;
  double scale ;


  double **zzn, **zzd ;
  int popx[4] ;
  double tn[4*5], td[4*4] ;
  double zzsig[5], zzest[5], zsc[5] ;
  double ymin ;
  double *f2, *f2sig, *fst ;
  double *ff3, *ff3var, *ff3fit, *ff3sig ;

  double *f3, *f4, *f3sig, *f4sig, *www, *ww2 ;
  int ng2, ng3, ng4 ;
  double *pwts ;
  double *ppwts ;
  char **cnames ;
  char **vnames ;
  double yscore ; 
  double *xxans ;
  double *wwtemp ;
  int nwts, nforce ;
  int iter, dof ;
  double ytail ;
  int *xpopsize ;
  FILE *f3file = NULL ;
  FILE *phylipfile = NULL ;
  int dmode = NO ;


  readcommands(argc, argv) ;
  printf("## qpgraph version: %s\n", WVERSION) ;
  if (parname == NULL) return 0 ;
  if (outpop != NULL) printf("outpop:  %s\n", outpop) ;

   if (seed==0) seed = seednum() ;
   SRAND(seed) ;
   printf("seed: %d\n", seed) ;

   setinbreed(inbreed) ;

  if (xchrom == 23) noxdata = NO ;
  if (lsqmode) printf("simple lsqmode\n") ;
  if (outputname != NULL)  openit(outputname, &ofile, "w") ;
  if (f3name != NULL) openit (f3name, &f3file, "w") ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  setindm(indivmarkers) ;
  k = getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  for (i=0; i<numsnps; i++)  {  
     cupt = snpmarkers[i] ;
     chrom = cupt -> chrom ;
     if ((noxdata == YES) && (cupt -> chrom == 23)) cupt -> ignore = YES ;
     if ((noxdata == NO) && (cupt -> chrom != 23)) cupt -> ignore = YES ;
     if ((xchrom>0) && (xchrom != chrom)) cupt -> ignore = YES ;
     if ((xnochrom>0) && (xnochrom == chrom)) cupt -> ignore = YES ;
  }


  printnl()  ;

  
  if ((doanalysis == NO)  && (poplistname != NULL)) { 
   ZALLOC(eglist, MAXPOPS, char *) ;
   numeg = loadlist(eglist, poplistname) ;
  }
  else {
   if (graphname == NULL) fatalx("no graph given\n") ;
   printf("graph: %s\n", graphname) ;
/**
   sprintf(sss, "cat %s", graphname) ;
   system(sss) ;  
   fflush(stdout) ;
*/
   printnl()  ;
   printnl()  ;
   numeg = loadgraph(graphname, &eglist) ;
   if (rootname != NULL) reroot(rootname) ;
   nedge = getnumedge() ;
   ZALLOC(enames, nedge, char *) ;
   getenames(enames) ;
  }
  xnumeg = numeg ;
  for (i=0; i<numeg; i++) {
   setstatus(indivmarkers, numindivs, eglist[i]);
  }

 ZALLOC(egshort, numeg, char *) ;
 for (i=0; i<numeg; i++) {  
  egshort[i] = strdup(getshort(eglist[i], 5)) ;
  printf("%3d %s\n",i, eglist[i]) ;
 }

 if (!doanalysis) printf("jackknife block size: %9.3f\n", blgsize) ;

  outnum = -1 ; 
  basenum = 0 ; 

  getrootlabel(sss) ;
  printf("root label: %s\n", sss) ;
  if (outpop == NULL) outpop = strdup(sss) ;

  t = strcmp(outpop, "NULL") ; 
  if (t==0) outnum = -99 ;

  for (k=0; k<numeg; ++k)  {  
    t = strcmp(eglist[k],outpop) ;
    if (t==0) outnum = k ;
    t = strcmp(eglist[k],sss) ;
    if (t==0) basenum = k ;
  }
  if (outnum == -1) { 
   eglist[numeg] = strdup(outpop) ;  
   outnum = numeg ;
   ++numeg ;
  }
  x = (1<<numeg) ;
  subsets = initarray_2Dint(x, numeg, 0) ;

  if (outnum >= 0) outpop = eglist[outnum] ;
  basepop = eglist[basenum] ;

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

  ZALLOC(xpopsize, numeg, int) ;
  for (i=0; i<nrows; i++) {
    k = xtypes[i] ;
    ++xpopsize[k] ;
  }

  for (i=0; i<numeg; i++) {
    printf("population: %3d %20s %4d\n",i, eglist[i], xpopsize[i]) ;
  }
  for (i=0; i<numeg; i++) {
   if (xpopsize[i] == 0) fatalx("zero popsize\n") ;
  }


  printf("before setwt numsnps: %d  outpop: %s\n", numsnps, outpop) ;
  setwt(snpmarkers, numsnps, indivmarkers, nrows, xindex, xtypes, outpop, eglist, numeg) ;
  numsnps = rmsnps(snpmarkers, numsnps) ;  //  rid ignorable snps
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

  ZALLOC(popsizes, numeg, int) ;
  cntpops(popsizes, indivmarkers, numindivs, eglist, numeg) ;
  setpopsizes(popsizes, eglist, numeg) ;

  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;


  printf("snps: %d  indivs: %d\n", ncols, nrows) ;
  setblocks(blstart, blsize, &xnblocks, xsnplist, ncols, blgsize)  ;
 
  nblocks = xnblocks ;

  zz = initarray_2Ddouble(numeg*numeg, ncols, 0.0) ;
  qplist = initarray_2Dint(numeg*numeg*numeg*numeg, 4, 0) ;
  ZALLOC (qpscore,numeg*numeg*numeg*numeg, double) ;

  ZALLOC(pmean, numeg, double) ; 
  ZALLOC(pnum, numeg, double) ;  
  ZALLOC(rawcol, nrows, int) ;  

  vmix = initarray_2Ddouble(MAXG, MAXW, 0.0) ;
  ZALLOC(lmix, MAXG, int) ; 
  nmix = 0 ;

  ng2 = numeg*numeg ;  
  ng3 = numeg*ng2   ;  
  ng4 = numeg*ng3   ;  

  ZALLOC(hest, ng2, double) ;  
  ZALLOC(hsig, ng2, double) ;  

  ZALLOC(f2, ng2, double) ;
  ZALLOC(f2sig, ng2, double) ;
  ZALLOC(fst, ng2, double) ;
  ZALLOC(www, ng2, double) ;
  ZALLOC(ww2, MAXG, double) ;

  ZALLOC(ff3, ng2, double) ;
  ZALLOC(ff3fit, ng2, double) ;
  ZALLOC(ff3sig, ng2, double) ;
  ZALLOC(ff3var, ng4, double) ;
// differences from basepoint 
  dmode = NO ; 
  if (fstdmode) dmode = 2 ;

  scale = dofstnumx(fst, f2, f2sig, xsnplist, xindex, xtypes, nrows, ncols, numeg, xnblocks, indivmarkers, dmode) ;     
  if (lambdascale <= 0.0) lambdascale = scale ; 
  else { 
   vst(f2, f2, lambdascale/scale, ng2) ;
   vst(f2sig, f2sig, lambdascale/scale, ng2) ;
  }
  printf("lambdascale: %9.3f\n", lambdascale) ;


  printf("fst:") ;
  printnl() ;
  printmatz(fst, egshort, numeg) ;
  printnl() ;
  printf("f2:") ; 
  printnl() ;
  printmatz(f2, egshort, numeg) ;
  printnl() ;
  printf3("f2:", f3file, f2, egshort, numeg) ;

  if (phylipname != NULL) { 
   openit(phylipname, &phylipfile, "w") ;
   fprintf(phylipfile, "%6d\n",numeg) ; 
   sss[10] = CNULL ;
    for (k1=0; k1<numeg; ++k1) { 
     strncpy(sss, eglist[k1], 10) ;  
     fprintf(phylipfile, "%10s", sss) ;
     for (k2=0; k2<numeg; ++k2) {  
      y1 = f2[k1*numeg+k2] ; 
      y2 = f2[k2*numeg+k1] ; 
      fprintf(phylipfile, "%6.3f", (0.5*(y1+y2))) ;
     }
     fprintf(phylipfile, "\n") ;
    }
   fclose(phylipfile) ;
  }

  
  ZALLOC(f3, ng3, double) ;
  ZALLOC(f4, ng4, double) ;
  ZALLOC(f3sig, ng3, double) ;
  ZALLOC(f4sig, ng4, double) ;

  nh2 = numeg*(numeg-1) ; 
  nh2 /= 2 ;
  ZALLOC(ind2f, nh2, int) ;
  ZALLOC(f2ind, ng2, int) ;

  ZALLOC(vest, nh2, double) ;
  ZALLOC(vvar, nh2*nh2, double) ;
  ZALLOC(xvvar, nh2*nh2, double) ;  // adjusted
  ZALLOC(vvinv, nh2*nh2, double) ;
  ZALLOC(xvvinv, nh2*nh2, double) ;

   k = 0 ;
   for (a=0; a<numeg; ++a) {  
    if (a==basenum) continue ; 
     for (b=a; b<numeg; ++b) {  
       if (b==basenum) continue ; 
       ind2f[k] = a*numeg + b ;
       f2ind[a*numeg+b] = k ;
       f2ind[b*numeg+a] = k ;
       ++k ;
     }
   }

  doff3(ff3, ff3var, xsnplist, xindex, xtypes, nrows, ncols, numeg, nblocks, lambdascale) ;     
// side effect vvar made


  for (k=0; k<MIN(4, numeg-3); ++k) { 
   x = listsubset(subsets, numeg-1, k) ;
  }

  for (u=0; u<nh2; ++u) {  
   y = vvar[u*nh2+u] ;  
   vvar[u*nh2+u] = MAX(y, minvar) ;
  }

  copyarr(vvar, xvvar, nh2*nh2) ;

  for (u=0; u<nh2; ++u) {  
   y = vvar[u*nh2+u] ;  
   xvvar[u*nh2+u] += y*f2diag ; 
  }
  copyarr(vvar, xvvinv, nh2*nh2) ;
// wipe out diagonal
  for (u=0; u<nh2; ++u) {  
   for (v=0; v<nh2; ++v) {  
     if (u==v) continue ;
     xvvinv[u*nh2+v] = 0.0 ;
   }
  }


// print some values
  printf("ff3:") ;
  printnl() ;
  printmatz(ff3, egshort, numeg) ;
  printnl() ;
  printf3("ff3:", f3file, ff3, egshort, numeg) ;

  for (a=0; a<numeg; ++a)  { 
   for (b=0; b<numeg; ++b)  { 
     y= dump4(ff3var, a, b, a, b, numeg) ;
     y += 1.0e-12 ;
     ff3sig[a*numeg+b] = 10*sqrt(y) ;
   }
  }

  printf("ff3sig*10:") ;
  printnl() ;
  printmatz(ff3sig, egshort, numeg) ;
  printnl() ;
  printf3("ff3sig:", f3file, ff3sig, egshort, numeg) ;

  for (a=0; a<numeg; ++a)  { 
   for (b=0; b<numeg; ++b)  { 
    if (verbose == NO) break ; 
    printf("%4d %4d %9.3f %15.9f  ", a, b, 
    dump2(f2, a, b, numeg)  ,
    dump2(f2sig, a, b, numeg)   ) ;
    getmv( a, b, a, b, &y1, &y2, vest, vvar) ;

    y2 = sqrt(y2) ;
    printf("%9.3f", y1) ;
    printf(" %15.9f", y2) ;
    printnl() ;
   }
  }

 if (quartet)   {   
  for (a=0; a<numeg; ++a)  { 
   for (b=a+1; b<numeg; ++b)  { 
    for (c=b+1; c<numeg; ++c)  { 
     for (d=c+1; d<numeg; ++d)  { 
      print4(ff3, a, b, c, d, numeg) ;
      print4(ff3, a, c, b, d, numeg) ;
      print4(ff3, a, d, b, c, numeg) ;
      printnl() ;
     }
    }
   }
  }
 }

  if (doanalysis == NO) { 
   printf("no analysis\n") ;
   printf("end of run\n") ;
   return 0 ;
  }
  printf("starting analysis\n") ;
  pdinv(vvinv, xvvar, nh2) ;
  pdinv(xvvinv, xvvinv, nh2) ;
  for (a=1; a<numeg; a++) { 
   u = f2ind[a*numeg+a] ;
   xvvinv[u*nh2+u]  *= f2weight ;
  }
// xvvinv used in lsq mode

/**
  printmatw(vvar, nh2, nh2, nh2) ;
  printnl() ;
  printmatw(vvinv, nh2, nh2, nh2) ;
  printnl() ;
*/

  ZALLOC(pwts, numeg * MAXG, double) ;
// initialization of weights done in loadgraph
  getgmix(vmix, lmix, &nmix) ;

  nwts = intsum(lmix, nmix) ;
  getpwts(pwts, &nrows, &nedge) ;
  ZALLOC(xxans, nedge, double) ;
  ZALLOC(wwtemp, nwts+1, double) ;

  if (loadname == NULL) {
   initvmix(wwtemp, nwts, 100*nmix) ; // side effect sets vmix

   for (k=0; k<nmix; ++k)  {  
    printmat(vmix[k], 1, lmix[k]) ;
   }
   printnl() ;
   printnl() ;
  }

  if (loadname != NULL)  { 
   loadpars(loadname, wwtemp, nwts, xxans, nedge) ;
/**
   printf("zz2") ;
   printmat(wwtemp, 1, nwts) ;
   printmat(xxans, 1, nedge) ;
   printnl() ;
*/
   getwww(vmix, wwtemp, nwts) ;
   putgmix(vmix) ;
   calcfit(ff3fit,  xxans, numeg) ;
// now print answers  
   printf("initial\n") ;
   printf("ff3fit:\n") ;
   printmatz(ff3fit, egshort, numeg) ;
   printnl() ;
   printf("ff3diff:\n") ;
   vvm(ff3fit, ff3fit, ff3, ng2) ;
   printmatz(ff3fit, egshort, numeg) ;
   printnl() ;
  printf3("ff3diff:", f3file, ff3fit, egshort, numeg) ;
   printvals(vmix, xxans, nedge) ;
   printfit(xxans) ;

  }



  ZALLOC(ppwts, nh2*nedge, double) ;  // coeffs of x_{ab} = (p_0 - p_a) (p_0-p_b)
  loadppwts(ppwts, pwts, nedge) ;  

/**
  for (k=0; k<nedge; ++k) { 
   printf("%9s ", enames[k]) ; 
  }
  printnl() ;
  printmatw(pwts, nrows, nedge, nedge) ;
  printnl() ;
  printmatw(ppwts, nh2, nedge, nedge) ;
*/

  yscore = calcxx(xxans, vvinv, ppwts, vest, nh2, nedge) ;

//  printf ("score: %9.3f ", yscore) ;  printmatw(xxans, 1, nedge, nedge) ;

  for (iter=1; iter<=bigiter; ++iter) {

  nagopt(wwtemp, nwts) ;
  getwww(vmix, wwtemp, nwts) ;
  putgmix(vmix) ;
  if (details) verbose = YES ; 
  y = scorit(wwtemp, nwts, &y1, ww2) ;

  if (verbose) printnl() ;
  dof = nh2 - nedge + intsum(ezero, nedge) ;
  dof -= (nwts-nmix) ;
  ytail = rtlchsq(dof, y) ;
  if (verbose || (iter == bigiter)) 
   printf(" score: %12.3f  dof: %d nominal p-val %12.6f constraint: %12.6f\n", y, dof, ytail, y1)  ;

  if (forcezmode)  {  
   nforce = 0 ;
   for (k=0; k<nedge; ++k) { 
    if (ww2[k] <= 0.0) { 
     ezero[k] = 1 ; 
     if (verbose) 
     printf("forcing %s to zero\n", enames[k]) ;
     ++nforce ;
    }
   }
   if (nforce > 0) {  
     verbose = NO ;
     nagopt(wwtemp, nwts) ;
     getwww(vmix, wwtemp, nwts) ;
     putgmix(vmix) ;
     if (details) verbose = YES ; 
     y = scorit(wwtemp, nwts, &y1, ww2) ;
     if (verbose) printnl() ;
     dof = nh2 - nedge + intsum(ezero, nedge) ;
     dof -= (nwts-nmix) ;
     ytail = rtlchsq(dof, y) ;
     if (verbose) 
     printf(" score: %12.3f  dof: %d nominal p-val %12.6f constraint: %12.6f\n", y, dof, ytail, y1)  ;
   }
  }
 }
 putgmix(vmix) ;
// verbose = YES ;
 y = scorit(wwtemp, nwts, &y1, ww2) ;
 calcfit(ff3fit,  ww2, numeg) ;
// now print answers  
	  printf("ff3fit:\n") ;
	  printmatz(ff3fit, egshort, numeg) ;
	  printnl() ;
          printf3("ff3fit:", f3file, ff3fit, egshort, numeg) ;
	  printf("ff3diff:\n") ;
	  vvm(ff3fit, ff3fit, ff3, ng2) ;
	  printmatz(ff3fit, egshort, numeg) ;
	  printnl() ;
          printf3("ff3diff:", f3file, ff3fit, egshort, numeg) ;
	  printvals(vmix, ww2, nedge) ;
	  printfit(ww2) ;
          putewts(ww2) ;  // load graph

	  dumppars(dumpname, wwtemp, nwts, ww2, nedge) ;
	  dumpgraph(graphoutname) ;
          dumpdotgraph(graphdotname) ;

	  printf("## end of run\n") ;
	  return 0 ;
}
void
calcfit(double *ff3fit, double *ww, int numeg) 
{
  double *pwts, *ppwts, *bestf ;
  double *xwts ; 
  int u, x ;
  int a, b, c, d ; 
  int nrows, nedge, ng2 ;
  double y1, y2, x1, x2, diff, sig, z ;
  double *zz ;

  ZALLOC(pwts, numeg * MAXG, double) ;
  getpwts(pwts, &nrows, &nedge) ;
  ZALLOC(ppwts, nh2*nedge, double) ;  // coeffs of x_{ab} = (p_0 - p_a) (p_0-p_b)
  loadppwts(ppwts, pwts, nedge) ;  
  ZALLOC(bestf, nh2, double) ;

  xwts = ppwts ;
  ng2 = numeg+numeg ;
  vzero(ff3fit, ng2) ;
  for (u=0; u<nh2; ++u) {  
   y1 = bestf[u] = vdot(xwts, ww, nedge) ;
   x = ind2f[u] ; 
   b = x / numeg ;
   c = x % numeg ;
   ff3fit[b*numeg+c] = y1 ;
   ff3fit[c*numeg+b] = y1 ;
   xwts += nedge ;
  }


  free(pwts) ; 
  free(ppwts) ; 
  free(bestf) ;

}
void printfit(double *ww) 
{
  double *pwts, *ppwts, *bestf ;
  double *xwts ; 
  int u ;
  int a, b, c, d, t ;
  int nrows, nedge ;
  double y1, y2, x1, x2, diff, sig, z ;
  double *zz ;
  FILE *outff ;

  ZALLOC(pwts, numeg * MAXG, double) ;
  getpwts(pwts, &nrows, &nedge) ;
  ZALLOC(ppwts, nh2*nedge, double) ;  // coeffs of x_{ab} = (p_0 - p_a) (p_0-p_b)
  loadppwts(ppwts, pwts, nedge) ;  
  ZALLOC(bestf, nh2, double) ;

  xwts = ppwts ;
  for (u=0; u<nh2; ++u) {  
   bestf[u] = vdot(xwts, ww, nedge) ;
   xwts += nedge ;
  }
  if (outliername != NULL) openit(outliername, &outff, "w") ;

  printnl() ;
  printnl() ;

  printf("%10s %10s ", "", "") ;              
  printf("fst: %12s %12s ", "fitted", "estim") ;
  printf("%12s ", "diff") ;
  printf("%12s ", "std. err") ; 
  printf("%9s ", "Z score") ;
  printnl() ; 

  for (a=0; a<numeg; a++)  { 
   for (b=a+1; b<numeg; b++)  { 
    getmv( a, b, a, b, &y1, &y2, vest, vvar) ;
    getmv( a, b, a, b, &x1, &x2, bestf, vvar) ;
    diff = y1 - x1 ;
    sig = sqrt(y2)  ;
    z = diff/sig ;
    printf("%10s %10s ", egshort[a], egshort[b]) ;
    printf("f2: %12.6f %12.6f ", x1, y1) ;
    printf("%12.6f ", diff) ;
    printf("%12.6f ", sig) ;
    printf("%9.3f ", z) ;
    printnl() ;
   }
  }
  printnl() ;
  printnl() ;
  for (a=1; a<numeg; a++) { 
   for (b=1; b<numeg; b++) { 
     if (details == NO) break ; 
     u = f2ind[a*numeg+b] ;  
     printf("%10s %10s ff3fit: ", egshort[a], egshort[b]) ;
     printf("%12.6f ", bestf[u]) ;
     printf("%12.6f ", vest[u]) ;
     printnl() ;
   }
  }

  printf("outliers:") ;
  printnl() ;

  for (a=0; a<numeg; a++)  { 
   for (b=0 ; b<numeg; b++)  { 
    for (c=0;  c<numeg; c++)  { 
     for (d=0 ; d<numeg; d++)  { 
      if (iscanon(a,b,c,d) == NO) continue ;
      getmv( a, b, c, d, &y1, &y2, vest, vvar) ;
      getmv( a, b, c, d, &x1, &x2, bestf, vvar) ;
      diff = y1 - x1 ;
      sig = sqrt(y2)  ;
      z = diff/sig ;
      if (fabs(z) < zthresh) continue ;
      printf("%10s %10s ", egshort[a], egshort[b]) ;
      printf("%10s %10s ", egshort[c], egshort[d]) ;
      printf(" %12.6f %12.6f ", x1, y1) ;
      printf("%12.6f ", diff) ;
      printf("%12.6f ", sig) ;
      printf("%9.3f ", z) ;
      printnl() ;
      if (outliername != NULL) {
       fprintf(outff,"%10s %10s  ", eglist[a], eglist[b]) ;
       fprintf(outff,"%10s %10s ", eglist[c], eglist[d]) ;
       fprintf(outff," %12.6f %12.6f ", x1, y1) ;
       fprintf(outff,"%9.3f ", z) ;
       t = numqq(a, b, c, d) ; 
       fprintf(outff,"%d ", t) ;
       fprintf(outff,"\n") ;
       
       printnl() ;
      }
     }
    }
   }
  }

  printnl() ;
  free(pwts) ;
  free(ppwts) ;
  free(bestf) ;
  if (outliername != NULL) fclose(outff) ;
}
int numqq(a, b, c, d) 
{
 int xx[100] ;
 ivzero(xx, 100) ;
 xx[a] = xx[b] = xx[c] = xx[d] = 1 ;
 return intsum(xx, 100) ;
}
int iscanon(int a, int b, int c, int d) 
// test quartet for is it canonical?  
{
 if (a>=b) return NO ;
 if (c>=d) return NO ;
 if (a>c) return NO ; 
 if ((a==c) && (b>d)) return NO ;

 return YES ;

}
void
print4(double *ff3, int a, int b, int c, int d, int numeg) 
{
   double y ;

   y = ff4val(ff3, a, b, c, d, numeg) ;
   printf("%4s ", get3(egshort[a])) ;
   printf("%4s ", get3(egshort[b])) ;
   printf("   ")  ;  
   printf("%4s ", get3(egshort[c])) ;
   printf("%4s ", get3(egshort[d])) ;
   printf("%6d", nnint(1000.0*y)) ;
   printnl() ;
}
double ff4val(double *ff3, int a, int b, int c, int d, int numeg) 
{
    double y1, y2, y3, y4 ;

    y1 = dump2(ff3, a, c, numeg) ;
    y2 = dump2(ff3, a, d, numeg) ;
    y3 = dump2(ff3, b, c, numeg) ;
    y4 = dump2(ff3, b, d, numeg) ;

    return y1 +y4 - (y2 + y3) ;

}

void printvals(double **tmix, double *edgelen, int nedge) 
{
 char sss[MAXSTR]  ;
 int k, n ;

 for (k=0; k<nmix; ++k) { 
  getmixstr(k, sss) ;
  printf("%40s ", sss) ; 
  printmat(tmix[k], 1, lmix[k]) ;
 }
  printnl() ;

  for (k=0; k<nedge; ++k) { 
   printf("%9s ", enames[k]) ; 
  }
  printnl() ; 
  printmatwf(edgelen, 1, nedge, nedge, "%9.4f ") ;

}

void setwww(double **tmix, double *www, int n)  
// copy tmix to vector 
{
  int k, l ;
  double *ww ; 

  if (n != intsum(lmix, nmix)) fatalx("dimension bug\n") ;

  ww = www ;
  for (k=0; k<nmix; ++k) {  
   l = lmix[k] ;
   copyarr(tmix[k], ww, l) ;
   bal1(ww, l) ;
   ww += l ;
  }

}

void getwww(double **tmix, double *www, int n)  
// copy vector to tmix  
{
  int k, l  ;
  double *ww ; 

  if (n != intsum(lmix, nmix)) fatalx("dimension bug\n") ;

  ww = www ;
  for (k=0; k<nmix; ++k) {  
   l = lmix[k] ;
   copyarr(ww, tmix[k], l ) ;
   ww += l ;
  }
}
void dumppars(char *dumpname, double *www, int nwts, double *xxans, int nedge) 
{
 FILE *dumpfile ;

 if (dumpname == NULL) return ;
 openit (dumpname, &dumpfile, "w") ;
 dump1(dumpfile, www, nwts) ;
 dump1(dumpfile, xxans, nedge ) ;
 
 fclose(dumpfile) ;
}
void
dump1(FILE *dumpfile, double *ww, int n) 
{
  int k ; 
  for (k=0; k<n; ++k) fprintf(dumpfile, "%12.6f\n", ww[k]) ;

}
void loadpars(char *loadname, double *www, int nwts, double *xxans, int nedge) 
{
 FILE *loadfile ;

 if (loadname == NULL) return ;
 printf("zzloadpars %d %d\n", nwts, nedge) ;
 openit (loadname, &loadfile, "r") ;
 read1(loadfile, www, nwts) ;
 read1(loadfile, xxans, nedge ) ;
 printmat(www, 1, nwts) ;  

 printnl() ; 
 printnl() ; 
 printmat(xxans, 1, nedge) ;  
 
 fclose(loadfile) ;

}

void
read1(FILE *loadfile, double *ww, int n) 
{

  int k ; 
  vclear(ww, -1.0, n) ;
  for (k=0; k<n; ++k) { 
   fscanf(loadfile, "%lf\n", &ww[k]) ;
   printf("zzfs %d %9.3f\n", k, ww[k]) ;
  }

}

void normvec(double *www, int n) 
{
  double **tmix ;

  tmix = initarray_2Ddouble(nmix, MAXW, 0.0) ;

  getwww(tmix, www, n) ;
  setwww(tmix, www, n) ;
  free2D(&tmix, nmix) ;
}


double scorit(double *www, int n, double *pfix, double *ans) 
// pfix and ans can be NULL
{

  double **tmix ;
  double *pwts, *ppwts, *xxans ;
  double yfix, yscore, y, yy ;

  int k, l, j ;
  double *ww ; 
  double *ppinv ;
  int nrows, nedge ;

/**
  vmaxmin(www, n, NULL, &y) ;  
  if ((n>0) && (y < wtmin)) {
    printf("zzwww\n") ;
    printmat(www, 1, n) ;
    return 99999.0 ;  
  }
*/
  if (n != intsum(lmix, nmix)) fatalx("dimension bug\n") ;
  tmix = initarray_2Ddouble(MAXG, MAXW, 0.0) ;
  getwww(tmix, www, n) ;
  yfix = 0.0 ;
  for (k=0; k<nmix; ++k) {  
   l = lmix[k] ;
   y = asum(tmix[k], l) ;
   for (j=0; j<l; ++j) { 
    yy = tmix[k][j] ;  
    if (yy<0.0) { 
     tmix[k][j] = fabs(yy) ;
     yfix += 10000*fabs(yy) ;
    }
   }
   vabs(tmix[k], tmix[k], l) ;
   
   bal1(tmix[k], l) ;
   yfix += (y-1.0)*(y-1.0) ;
   ww += l ;
  }

  putgmix(tmix) ;
  ZALLOC(pwts, numeg * MAXG, double) ;
  getpwts(pwts, &nrows, &nedge) ;
  ZALLOC(xxans, nedge, double) ;
  ZALLOC(ppwts, nh2*nedge, double) ;  // coeffs of x_{ab} = (p_0 - p_a) (p_0-p_b)
  loadppwts(ppwts, pwts, nedge) ;  

  if (verbose) { 
   printf("zzfix: %12.3f\n", yfix) ;
   printf("path and path products\n") ;
   for (k=0; k<nedge; ++k) { 
    printf("%9s ", enames[k]) ; 
   }
   printnl() ;
   printmatw(pwts, nrows, nedge, nedge) ;
   printnl() ;
   printmatw(ppwts, nh2, nedge, nedge) ;
   printnl() ;
  }
  ppinv = vvinv ; 
  if (lsqmode) ppinv = xvvinv ;
  yscore = calcxx(xxans, ppinv, ppwts, vest, nh2, nedge) ;

  if (pfix != NULL) *pfix = yfix ;
  if (ans != NULL)  copyarr(xxans, ans, nedge) ;
  free2D(&tmix, MAXG) ;
  free(pwts) ;
  free(ppwts) ;

  return yscore ;

}
double 
calcxx(double *xxans, double *qmat, double *ppwts, double *rhs, 
 int nrow, int ncol)  

/** 
 Set L(a) = P(a, *) * X(*) 
 minimize Q(a,b) (L(a) - R(a))(L(b)-R(b))
 returns min value 
*/

{            

 double *cc, *rr, *ll, *w1, *w2 ;
 double *pa, *pb ; 
 int a,  b, i, j, kret, u, x, c ;
 double y ;

 if (ezero == NULL) { 
   ZALLOC(ezero, ncol, int) ;
   getezero(ezero) ;
 }

 ZALLOC(cc, ncol*ncol, double) ;
 ZALLOC(rr, ncol, double) ;
 ZALLOC(ll, nrow, double) ;
 ZALLOC(w1, nrow, double) ;
 ZALLOC(w2, nrow, double) ;

 for (a=0; a<nrow; a++) { 
  for (b=0; b<nrow; b++) { 
   y = qmat[a*nrow+b] ;
   pa = ppwts+a*ncol ; 
   pb = ppwts+b*ncol ;
   for (i=0; i<ncol ; i++) { 
    rr[i] += y*pa[i]*rhs[b] ;
    if (verbose && (a==0) && (b==0)) { 
     printf("zzbad %d %d %12.6f %12.6f %12.6f\n", a, b, y, pa[i], rhs[b]) ;
    }
    for (j=0; j<ncol ; j++) { 
     cc[i*ncol+j] += y*pa[i]*pb[j] ;
    }
   }
  }
 }

 y = (double) ncol / trace(cc, ncol) ; 
 if (diag>0) { 
   for (i=0; i<ncol ; i++) { 
    cc[i*ncol+i] += diag/y ;
   }
 }
 forcez(cc, rr, ezero, ncol) ;
 kret = solvit(cc, rr, ncol, xxans) ;
 if (verbose) { 
  printf("verbose details:\n") ;
  printmatw(cc, ncol, ncol, ncol) ;
  printnl() ;
  printmatw(rr, 1, ncol, ncol) ;
  printnl() ;
 }

 if (kret<0) { 
  free(cc) ; 
  free(rr) ; 
  free(ll) ; 
  free(w1) ; 
  free(w2) ;
  return 10000.0 ;  
  printf("zz %d %d\n", nrow, ncol) ;
  printmatw(cc, ncol, ncol, ncol) ;
  printnl() ; 
  printmatw(rr, 1, ncol, ncol) ;
  printnl() ; 
//fatalx("problem not well posed\n") ;
 }
 mulmat(ll, ppwts, xxans, nrow, ncol, 1) ;
 vvm(w2, ll, rhs, nrow) ;
 if (verbose) { 
  for (u=0; u<nrow; u++) {
   x = ind2f[u] ; 
   b = x / numeg ;
   c = x % numeg ;
   printf("%10s %10s ", egshort[b], egshort[c]) ;
   printf(" %9.3f ", rhs[u]) ;
   printf(" %9.3f ", ll[u]) ;
   printf(" %9.3f ", w2[u]) ;
   printnl() ;
  }
 }
 mulmat(w1, qmat, w2, nrow, nrow, 1) ;
 y = vdot(w1, w2, nrow) ;

 free(cc) ; 
 free(rr) ; 
 free(ll) ; 
 free(w1) ; 
 free(w2) ;

 return y ;
}
void forcez(double *cc, double *rr, int *ez, int n) 
{
// force variable k zero if set
  int j, k ;

  for (k=0; k<n; ++k) { 
   if (ez[k] == 0) continue ;
   rr[k] = 0 ;
   for (j=0; j<n; j++) { 
    cc[k*n+j] = 0 ;
    cc[j*n+k] = 0 ; 
   }
   cc[k*n+k] = 1.0 ; 
  }
}
void loadppwts(double *ppwts, double *pwts, int n) 
{
    int u, x, b, c ;
    for (u=0; u<nh2 ; u++) { 
       x = ind2f[u] ; 
       b = x / numeg ;
       c = x % numeg ;
       --b ; 
       --c ;
       vvt(ppwts+u*n, pwts+b*n, pwts+c*n, n) ;
    }
}
void
getmv (int a, int b, int c, int d, double *mean, double *var, double *xest, double *xvar) 
{
    double *w1, *w2 ;
    int ng2 ; 
    double y1, y2 ;
    int s, t, k  ;

// E (a-b)(c-d) using f_3 stats
 
    ZALLOC(w1, nh2, double) ;
    ZALLOC(w2, nh2, double) ;

    
    y1 = 0.0 ;
    bumpm(&y1, a, c, 1.0, xest) ;
    bumpm(&y1, b, d, 1.0, xest) ;
    bumpm(&y1, a, d, -1.0, xest) ;
    bumpm(&y1, b, c, -1.0, xest) ;
    *mean = y1 ;

    bumpw(w1, a, c, 1.0) ;
    bumpw(w1, b, d, 1.0) ;
    bumpw(w1, a, d, -1.0) ;
    bumpw(w1, b, c, -1.0) ;

    mulmat(w2, xvar, w1, nh2, nh2, 1) ;
    *var = vdot(w1, w2, nh2) ;

    free(w1) ;
    free(w2) ;

}
void mkww(double *f3, double *f2, int k, double *ww, int n) 
{
  int i, j ;
  for (i=0; i<n; i++) { 
   for (j=0; j<n; j++) { 
    ww[i*n+j] = dump3(f3, k, i, j, n) ;
   }
   ww[i*n+i] = dump2(f2, k, i, n) ;
  }
}
void
printf3(char *sss, FILE *fff, double *ww, char **eglist, int n) 
{
  int i, j , x ;

  if (fff==NULL) return ;
  fprintf(fff, "%s\n", sss)  ;   
  fprintf(fff, " %4s", "   ") ;
  for (i=0; i<n; i++)  {  
   fprintf(fff, " %4s", get3(eglist[i])) ; 
  }
  fprintf(fff, "\n") ;
  for (i=0; i<n; i++)  {  
    fprintf(fff, "%4s", get3(eglist[i])) ;
    for (j=0; j<n; j++) { 
     x = nnint(1000*ww[i*n+j]) ; 
     fprintf(fff, " %4d", x)  ;              
    }
    fprintf(fff, "\n") ;
  }
  fprintf(fff, "\n") ;
}
void indiaestit(double *f2, double *f3, double *f4, int n) 
{

  double z1, z2, z3, z4, z5, z6 ;
  double u1a0, u2a0, v1m, v2m ;
  double y1, xa02, xb, xmd ;
  double u1, u2, v1, v2 ;
  double xu1, xu2, xa0, xa2, xm, xd ;
  double xa1, xt ;
  double r, s, x1, x2 ;
  double co[4], rr[2], ans[2] ;

  int a2, b, d ; 
  int c, c1, c2 ;

  a2 = E ; 
  b =  Y ; 
  d =  A ;

  z1 = dump2(f2, a2, b, n) ;
  z2 = dump2(f2, a2, d, n) ; 
  z3 = dump2(f2, b, d, n) ;

  y1 = z2-z3 ;  // a0 + a2 -  b  
  xa02 = 0.5 * (z1 + y1) ;  
  xb   = z1 - xa02 ;
  xmd  = z3 - xb ; 

  printf("a0 + a2: %9.3f\n", xa02) ;
  printf("      b: %9.3f\n", xb) ;
  printf("  m + d: %9.3f\n", xmd) ;

  for (c=d+1; c<n; ++c)  {
  z4 = dump3(f3, a2, b, c, n) ;
   u1a0 = xa02-z4 ;
   z5 = dump3(f3, d, b, c, n) ;
   v1m  = xmd - z5 ;
   printf("%10s %10.4f %10.4f\n", eglist[c], u1a0, v1m) ;
  }

}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n, t ;

  while ((i = getopt (argc, argv, "z:p:g:s:o:d:x:l:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'g':
	graphname = strdup(optarg) ;
	break;

      case 'o':
	graphoutname = strdup(optarg) ;
	break;

      case 'd':
	graphdotname = strdup(optarg) ;
	break;

      case 'x':
	outliername = strdup(optarg) ;
	break;

      case 'l':
	lambdascale = atof(optarg) ;
	break;

      case 'z':
	zthresh = atof(optarg) ;
	break;

      case 's':
	seed = atoi(optarg) ;
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
   getstring(ph, "graphname:", &graphname) ;
   getstring(ph, "graphoutname:", &graphoutname) ;
   getstring(ph, "graphdotname:", &graphdotname) ;
   getstring(ph, "phylipname:", &phylipname) ;
   getstring(ph, "outpop:", &outpop) ;
   getstring(ph, "output:", &outputname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "f3log:", &f3name) ;
   getstring(ph, "root:", &rootname) ;
   getdbl(ph, "blgsize:", &blgsize) ;
   getdbl(ph, "diag:", &diag) ;
   getdbl(ph, "f2diag:", &f2diag) ;
   getdbl(ph, "zthresh:", &zthresh) ;
   getdbl(ph, "minvar:", &minvar) ;
   getdbl(ph, "lambdascale:", &lambdascale) ;
   getint(ph, "bigiter:", &bigiter) ;
   getint(ph, "inbreed:", &inbreed) ;

   getint(ph, "noxdata:", &noxdata) ; 
   t = -1 ;
   getint(ph, "xdata:", &t) ; if (t>=0) noxdata = 1-t ;
   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "nochrom:", &xnochrom) ;
   getint(ph, "doanalysis:", &doanalysis) ; 
   getint(ph, "quartet:", &quartet) ; 

   getint(ph, "nostatslim:", &nostatslim) ; 
   getint(ph, "popsizelimit:", &popsizelimit) ; 
   getint(ph, "gfromp:", &gfromp) ;  // gen dis from phys
   getint(ph, "seed:", &seed) ;  
   getint(ph, "details:", &details) ;  
   getint(ph, "forcezmode:", &forcezmode) ;  
   getint(ph, "lsqmode:", &lsqmode) ;  
   getint(ph, "fstdmode:", &fstdmode) ;  
   getstring(ph, "dumpname:", &dumpname) ;
   getstring(ph, "loadname:", &loadname) ;
   getint(ph, "hires:", &hires) ;

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}
void sol2(double *co, double *rhs, double *ans) 
// solve 2x2 system
{

  double c00, c01, c10, c11, r0, r1, ra, xa ;

  c00=co[0] ; 
  c01=co[1] ;
  c10=co[2] ; 
  c11=co[3] ;
  r0=rhs[0] ;
  r1=rhs[1] ;

  ra = r0*c11-r1*c01 ;
  xa = c00*c11-c10*c01 ;  
  ans[0] = ra/xa ;
  if (c11 != 0.0) { 
   ra = r1-c10*ans[0] ;
   xa = c11 ;
  }
  else {             
   ra = r1-c00*ans[0] ;
   xa = c01 ;
  }
   ans[1] = ra/xa ;

}
void
doff3(double *ff3, double *ff3var, SNP **xsnplist, int *xindex, int *xtypes, 
 int nrows, int ncols, int numeg, int nblocks, double scale)     
{

   int t1, t2 ;
   int a, b, c ;
   int ng2, ng3 ;
   int c1[2], c2[2], *cc ;
   int *rawcol, *popall, *pop0, *pop1 ;
   int k, g, i, col, j ; 
   double ya, yb, y, jest, jsig, mean ;
   SNP *cupt ;
   double *top, *bot, *djack, *wjack, *gtop, *gbot, *wbot, *wtop ;
   double **btop, **bbot, wt ;
   double *w1, *w2, *w3 ;
   double ytop, ybot ;
   double y1, y2, yscal ;
   int bnum  ;
   int numegm = numeg -1 ; 
   int  u, v, x, kret ; 
   double *estmat ;
   
   ng2 = numeg*numeg ;
   ng3 = numeg*numeg*numeg ;
   
   ZALLOC(w1, ng3, double) ;
   ZALLOC(w2, ng3, double) ;
   ZALLOC(w3, ng3*numeg, double) ;
   ZALLOC(gtop, ng3, double) ;
   ZALLOC(gbot, ng3, double) ;
   ZALLOC(wtop, ng3, double) ;
   ZALLOC(wbot, ng3, double) ;
   ZALLOC(estmat, ng3, double) ;
   ZALLOC(djack, nblocks, double) ;
   ZALLOC(wjack, nblocks, double) ;
   btop  = initarray_2Ddouble(nblocks, ng3,  0.0) ;
   bbot  = initarray_2Ddouble(nblocks, ng3,  0.0) ;

   ZALLOC(vest, nh2, double) ;  
   ZALLOC(vvar, nh2*nh2, double) ;  

// printf("zz ") ;  printimat(ind2f, 1, nh2) ;

   for (col=0; col<ncols;  ++col)  {
    cupt = xsnplist[col] ;
    if (cupt -> ignore) continue ;
    wt = cupt -> weight ;  
    if (wt <= 0.0) continue ;
    bnum = cupt -> tagnumber ; 
    if (bnum<0) continue ;
    top = btop[bnum] ; 
    bot = bbot[bnum] ;

    kret = f3yyx(estmat,  cupt, xindex, xtypes, nrows, numeg, indivmarkers) ;
    if (kret<0) continue ;

    ++wjack[bnum] ;

    a = basenum ; 
    for (u=0; u<nh2 ; u++) { 
       x = ind2f[u] ; 
       b = x / numeg ;
       c = x % numeg ;
       ytop = dump3(estmat, a, b, c, numeg) ;
       if (fstdmode == NO) { 
        top[u] += wt*ytop ; 
        bot[u] += 1.0 ;
       }
       else { 
        top[u] += ytop ; 
        bot[u] += 1.0/wt ;
       }
     }
    }

    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvp(gtop, gtop, top, nh2) ;
     vvp(gbot, gbot, bot, nh2) ;
    }

   vsp(w2, gbot, 1.0e-10, nh2) ;
   vvd(w3, gtop, w2, nh2) ;
   vzero(ff3, numeg*numeg) ;
    for (u=0; u<nh2 ; u++) { 
       x = ind2f[u] ; 
       b = x / numeg ;
       c = x % numeg ;
       y1 = w3[u] ;                      
       bump2(ff3, b, c, numeg,y1)  ;
       if (b<c) bump2(ff3, c, b, numeg,y1)  ;
    }


    for (k=0; k<nblocks; k++) {  
     top = btop[k] ; 
     bot = bbot[k] ;
     vvm(wtop, gtop, top, nh2) ;
     vvm(wbot, gbot, bot, nh2) ;
     vsp(wbot, wbot, 1.0e-10, nh2) ;
     vvd(top, wtop, wbot, nh2) ;  // delete-block estimate
    }
    vsp(gbot, gbot, 1.0e-10, nh2) ;
    vvd(gtop, gtop, gbot, nh2) ;
    
    
    wjackvest(vest, vvar, nh2, gtop, btop, wjack, nblocks) ;
    vst(vest, vest, scale, nh2) ;
    vst(vvar, vvar, scale*scale, nh2*nh2) ;

/**
    vzero(w3, numeg*numeg) ;
    for (u=0; u<nh2 ; u++) { 
       x = ind2f[u] ; 
       b = x / numeg ;
       c = x % numeg ;
       y1 = vest[u] ;                      
       bump2(w3, b, c, numeg,y1)  ;
       if (b<c) bump2(w3, c, b, numeg,y1)  ;
    }
    printf("vest:\n") ;
    printmatz(w3, eglist, numeg) ;
    copyarr(vvar, w3, nh2*nh2) ;

    choldc(w3, nh2, w2) ;
    printmat(w2, 1, nh2) ;
    map4x(vvar, w3, numeg, ind2f) ;
    for (b=0; b<numeg; ++b)  {  
     set4(w3, 0, b, 0, b, numeg, 1.0) ;
    }
    choldc(w3, ng2, w2) ;
    printmat(w2, 1, ng2) ;

    printnl() ;
    printnl() ;
*/

    vst(ff3, ff3, scale, numeg*numeg) ;
    map4x(vvar, ff3var, numeg, ind2f) ;
//  vst(ff3var, ff3var, scale*scale, numeg*numeg*numeg*numeg) ;
    

    free(w1) ; 
    free(w2) ; 
    free(w3) ; 

    free(gbot) ; 
    free(wtop) ; 
    free(wbot) ; 
    free(estmat) ; 
    free(djack) ;
    free(wjack) ;

    free2D(&btop, nblocks);
    free2D(&bbot, nblocks);

}
void
map4y(double *aa, double *bb, int n2, int *indx) 
// map 4d array (n1 x n1 x n1 x n1  -> b  n2 x n2 x n2 x n2 
{
   int u, v, a, b, c, d ;
   int x ;
   double y1 ;
   int nh2 ;

   nh2 = n2*(n2-1) ; 
   nh2 /= 2 ;
   vzero(aa, nh2*nh2) ;
   for (u = 0 ; u <nh2; ++u) {  
    for (v = u ; v <nh2; ++v) {  

     x = indx[u] ; 
     a = x / n2 ;
     b = x % n2 ;
     x = indx[v] ; 
     c = x / n2 ;
     d = x % n2 ;
    
     y1 = dump4(bb, a, b, c, d, n2) ;
     aa[u*nh2+v] = y1 ;
     aa[v*nh2+u] = y1 ;
   }
  }
}
void bumpm(double *y, int a, int c, double val, double *xest)
{

  int k ;

  if (a==basenum) return ;
  if (c==basenum) return ;
  k = f2ind[a*numeg+c] ;
  *y  += val*xest[k] ;

}

void bumpw(double *w, int a, int c, double val)
{

  int k ;

  if (a==basenum) return ;
  if (c==basenum) return ;
  k = f2ind[a*numeg+c] ;
  w[k] += val ;

}
void
initvmix(double *wwinit, int nwts, int numiter) 
{
  double *ww, yold, ybest ;
  double *ww2 ;
  int iter, k ;
  int xniter ;
  double y ;
  int nedge, ng2 ;
  double *ff3fit ;

  if (nwts == 0) return ;
  nedge = getnumedge() ;
  ng2 = numeg*numeg ;  

  ZALLOC(ff3fit, ng2, double) ;

  ZALLOC(ww, nwts, double)  ;  
  ZALLOC(ww2, nedge, double)  ;  

  getgmix(vmix, lmix, &nmix) ;
  setwww(vmix, wwinit, nwts) ;
  copyarr(wwinit, ww, nwts) ;

//verbose = YES ;
  ybest = scorit(ww, nwts, NULL, ww2) ;
  if (verbose) 
  printf("initv: %4d %9.3f\n", -1, ybest) ;
  
  calcfit(ff3fit,  ww2, numeg) ;
// now print answers  
  if (verbose) {
	  printf("ff3fit:\n") ;
	  printmatz(ff3fit, egshort, numeg) ;
	  printnl() ;
	  printvals(vmix, ww2, nedge) ;
	  printfit(ww2) ;
   }


  for (iter=1; iter <= numiter; ++iter) {
    setrand(ww, nwts) ;
   y = scorit(ww, nwts, NULL, NULL) ;
// printf("inititer: %3d %12.3f\n", iter, y) ;
   if (y<ybest) {
    ybest = y ; 
    copyarr(ww, wwinit, nwts) ;
   }
  }
  getwww(vmix, wwinit, nwts) ;      
  putgmix(vmix)  ;
  getgmix(vmix, lmix, &nmix) ;
  setwww(vmix, wwinit, nwts) ;

  xniter = numiter/nwts ;
  if (verbose) 
  printf("initv: %4d %9.3f\n", 0, ybest) ;
  for (iter = 1; iter <= xniter ; ++iter) {  
   for (k=0; k<nmix; ++k) { 
    getwww(vmix, wwinit, nwts) ;      
    setsimp(vmix[k], lmix[k]) ;
    setwww(vmix, ww, nwts) ;
    y = scorit(ww, nwts, NULL, NULL) ;
    if (y<ybest) {
     ybest = y ;
     copyarr(ww, wwinit, nwts) ; 
    }
   }
  
   if (verbose) 
   printf("initv: %4d %9.3f\n", iter, ybest) ;
  }
  getwww(vmix, wwinit, nwts) ;      
  putgmix(vmix)  ;
  getgmix(vmix, lmix, &nmix) ;
  setwww(vmix, wwinit, nwts) ;
  free(ww) ;
  free(ww2) ;
  free(ff3fit) ;
}
int listsubset(int **x, int n, int k) 
// inefficient ;
{
   int *w ; 
   int mx, i, j, t, tt ;

   ZALLOC(w, n, int) ;
   mx = (1<<n) - 1 ;

   t = 0 ;
   for (i=0; i <= mx ; ++i) {  
    tt = i ;
    for (j=0; j<n; ++j) {  
     w[j] = tt & 1 ; 
     tt = tt >> 1 ;
    }
    if (intsum(w,n) != k) continue ;
    copyiarr(w, x[t], n) ;
    ++t ;
   }

  free(w) ;
  return t ;

}
