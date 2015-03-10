#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include  <math.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  
#include "ldsubs.h"  

#define WVERSION   "1002" 
/** an attempt to 
 1) threshold for LD in parentals 
 2) Use a weighted Z-score
 Also compute ld rolloff (runmode 2)
 binning change as in simpsim2 
 wts forced to have  mean 0
 minparentcount  added
 nochrom added
 zdipmode added
 mincount default now 5
 Jackknife output for 22 autosomes 
 Bug in clearcorr fixed
 Bug in dumpit callfixed
 xnochrom added
 checkmap added (default) 
 First admixtools release
 Count of outer + inner loop made
*/


#define MAXFL  50   
#define MAXSTR  512
#define MAXAPOPS  100

#define BINSIZE  .0005
#define MAXDIS   0.05  

extern int packmode ;
int zdipcorrmode = NO ; 

// just to keep qpsubs quiet;  not relevant for rolloff

char *trashdir = "/var/tmp" ;
int qtmode = NO ;
int debug = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
SNP **snpm ;
int numsnps, numindivs ; 
int checkmap = YES ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *admixpop = NULL ;
char *admixlist = NULL ;
char **admixpoplist = NULL  ;
int nadmixpop = 0; 
char *dumpname = NULL ;
char *loadname = NULL ;
char *oname = NULL ;
char *weightname = NULL ; 
int crossmode = NO ;
int runmode = 0 ;
int ldmode = NO ;
int ransample = -1 ;
int flatmode = NO  ;
double chithresh = -1.0; // default 6.0 
int jackknife = NO ;

FILE *ofile ;
int  *fchrom, *lchrom ;
int  numeg = 0 ;
int nums = 2,   ns2, ns4 ;
int numbigiter = 2 ;
int xchrom = -1 ;
int xnochrom = -1 ;
int seed = 0 ;

int noxdata = YES ;
int samelambda = NO ;
int fixlambda = NO ;
int flatinit = NO ;  // different default
int popestimate = NO ; 
int nxlim = 100000 ;
int mincount =  5 ;
int minparentcount = 10 ;
int norun = NO ;

double maxdis = MAXDIS   ;
double binsize = BINSIZE ; 

// if YES same across pops



char  unknowngender = 'U' ;

int readcommands(int argc, char **argv) ;
void setfc(SNP **snpm, int numsnps) ;
void sethmm(double *theta, double lambda, int ind) ;
int addr4(int a, int b, int c, int d) ;
double doalpha(int ind) ;
double dobeta(int ind) ;

void cleargcounts()   ;
void reest()   ;
void calcgcounts(int ind)  ;
void printpars()   ;
void printgams()   ;
void calcgtl(int ind) ;
double doscore(double **thet, double *lam) ;
void dumppars(char *dumpname,  double **wmix, double **theta, int nums, int numeg) ;
void dorc(double *ans, double *res, double binsize, double maxdis) ;
void dorc2(double *ans, double *res, double binsize, double maxdis) ;
void dorcx(double *ans, double *rr, double *dd, double binsize, double maxdis) ;
void dorcxx(double *ans, double *rr, double *dd, double binsize, double maxdis) ;
void dumpit(char *oname, double *ww, CORR **corrbins, int len, double bsize, char *hdr) ;
double cntit1(double *xc, SNP *cupt, Indiv **indm, int numindivs, int t) ; 
double cntit2(double *xc, SNP *cupt, SNP *cupt2, Indiv **indm, int numindivs, int t) ; 
void clipwt(SNP **snpmarkers, int numsnps) ;
double estyran(SNP **snpmarkers, int numsnps, int nxlim) ;
void cx1(double *xc, double *y1, double *y2) ;
double cminor(double *xx) ;
int cmap(SNP **snppmarkers, int numsnps) ;
int calcldscore(SNP *cupt, SNP *cupt2, double *xscore) ;

int main(int argc, char **argv)
{

  char s1[MAXSTR], s2[MAXSTR] ;
  int x, x2, i, j, k, g, s, t, tt, iter ; 
  int a, b, c ;
  SNP *cupt, *cupt2 ;
  Indiv *indx ;

  int numvind, nignore, numrisks = 1 ;
  char **eglist ;
  double y, y1, y2, yy, z1, z2, znorm, yn ;
  double *yy2 ;
  int ind, chrom, num0=0 ;
  double xc[9], xww[3] ;
  char sss[100] ;

  double *ttt ;
  double yscore, ybase ;
  double lambdabase ;
  double *scarray, *wacc, ywynn ;
  int ncount = 0 ;
  double lastscore = 0.0 ;
  double *tvecs, *pp ;
  double *w1, *w2, *res, alpha, *ww, *ww1, *ww2, *dd  ;
  double *xnum, *xden ;
  int numbins  ;  
  double co1, co2, dis, wt, yran ;     
  double wtmean, wtnum ; 
  CORR **corrbins ;
  CORR *corrpt  ;
  CORR ***corrjbins ;  //jackknife
  int nloop = 0 ;
  int nncount = 0 ;
  int nsnp1 = 0, nsnp2 = 0 ;
  int numchrom = 22 ;
  double **jwt ;
  long snpdo1, snpdo2 ;
 

  snpdo1 = snpdo2 = 0 ;
  ofile = stdout; 
  packmode = YES ;
  

  printnl() ;
  printf("## rolloff.  Version %s\n", WVERSION) ;
  printnl() ;

  ZALLOC(admixpoplist, MAXAPOPS, char *) ;  
  if (readcommands(argc, argv) < 0) return -1 ;

  if (weightname != NULL) {
   poplistname = NULL ;
   printf("weights set on input\n") ;
  }
  if (runmode == 2) ldmode = YES ;
  if (zdipcorrmode) setzdipmode(1) ;

  if (seed == 0) seed = seednum() ;
  SRAND(seed) ;
  printf("seed: %d\n", seed) ;

  if (oname == NULL) jackknife = NO ;

  

  if (chithresh < -0.5) {
   chithresh = 6.0 ;
  }

  if (ldmode) { 
   weightname = poplistname = NULL ; 
   runmode = 2 ; 
   chithresh = 0.0 ;
  }

  printf("chithresh: %9.3f\n", chithresh) ;

/**
  if (poplistname == NULL) fatalx("no poplistname\n") ; 
*/
  if (admixlist != NULL)  {  
   nadmixpop = loadlist(admixpoplist, admixlist) ;
  }
  if (nadmixpop == 0) {  
   if (admixpop == NULL) fatalx("no admixpop\n") ; 
    nadmixpop = 1 ;
    admixpoplist[0] = strdup(admixpop) ;
  }

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;

  printf("genotypename:  %s\n", genotypename) ;
  fflush(stdout) ;

   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

   if (checkmap && (cmap(snpmarkers, numsnps) == NO)) { 
    printf("running rolloff without a real map.  If you mean to do this set checkmap: NO\n") ;
    fatalx("no real map\n") ;
   }

   ZALLOC(eglist, numindivs, char *) ; 

   numeg = 0 ;
   if ((weightname == NULL)  && (runmode != 2) && (poplistname != NULL))
     numeg = loadlist(eglist, poplistname) ;

   if (weightname != NULL) { 
    printf("loading weights...\n") ; 
    getweights(weightname, snpmarkers, numsnps) ;
   }

   t = 0 ;
   for (i=0; i<numsnps; ++i) { 
    cupt = snpmarkers[i] ;
    cupt -> tagnumber = i ;
    chrom = cupt -> chrom ; 
    
    if (chrom<1) cupt -> ignore = YES ;
    if (chrom>22) cupt -> ignore = YES ;
    if ((xchrom>0) && (chrom != xchrom)) cupt -> ignore = YES ;
    if (chrom == xnochrom) cupt -> ignore = YES ;
    if (cupt -> ignore == NO) ++t ;
   }
   printf("valid snps: %d\n", t) ;

   numsnps = rmsnps(snpmarkers, numsnps, NULL) ;   

   setfc(snpmarkers, numsnps) ; // does allocation
   for (k=0; k<numeg; ++k) { 
    printf("group %d %s\n", k, eglist[k]) ;
   }

   ncount = 0 ;
   for (i=0; i<numindivs; ++i) { 
    indx = indivmarkers[i] ;
    indx -> affstatus = -1 ;
    t = indxindex(eglist, numeg, indx -> egroup) ;
    if (t>=0) { 
     indx -> affstatus = t ;
    }
    tt = indxindex(admixpoplist, nadmixpop, indx -> egroup) ;
    if (tt>=0) { 
     if (t>0) fatalx(" %s in admixing and admixed pops\n", indx -> egroup) ;
     indx -> affstatus = 99 ;
     ++ncount ;   
    }
    if (indx -> affstatus  <  0) indx -> ignore = YES ;
    if (indx -> ignore) continue ;

   }
   printf("number admixed: %d number of references: %d\n", ncount, numeg) ;
   mincount = MIN(mincount, ncount) ;  

   numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;

   printf("numsnps: %d  numindivs: %d\n", numsnps, numindivs) ;

   ZALLOC(w1, numsnps, double) ;
   ZALLOC(w2, numsnps, double) ;
   ZALLOC(res, numsnps, double) ;
   ZALLOC(dd, numsnps, double) ;

   numbins = nnint(maxdis/binsize) + 5 ;
   ZALLOC(ww, numsnps +numbins, double) ;

   ZALLOC(xnum, numbins, double) ;
   ZALLOC(xden, numbins, double) ;
   ZALLOC(yy2, numbins, double) ;

   vclear(ww, 1.0, numsnps) ;

   ZALLOC(corrbins, numbins, CORR *) ;

   for (x=0; x<numbins; ++x) { 
     ZALLOC(corrbins[x], 1, CORR) ;
     corrpt = corrbins[x] ; 
     clearcorr(corrpt) ;
   }

   ZALLOC(corrjbins, numchrom+1, CORR **) ;
   jwt = initarray_2Ddouble(numchrom+1, numbins, 0.0) ; 
   for (k=1; k<=numchrom; ++k) {
    ZALLOC(corrjbins[k], numbins, CORR *) ;

    for (x=0; x<numbins; ++x) { 
     ZALLOC(corrjbins[k][x], 1, CORR) ;
     corrpt = corrjbins[k][x] ; 
     clearcorr(corrpt) ;
    }
   }

   yn = (double) numsnps ;

/** 
 w1 is (normalized) parental allele freq difference 
 w2 is (normalized) (geno - expected)   where expected is freq average
*/

   for (x=0; x<numsnps; ++x)  { 

   if (weightname != NULL) { 
    break ;
   }
     
    cupt = snpmarkers[x] ;
    if (cupt -> ignore) continue ; 
    if (numeg < 2) break ;

    cntit1(xc, cupt, indivmarkers, numindivs, 0) ;  
     y1 = 2*xc[2] + xc[1] ;  
     y2 = 2*asum(xc,3) ;    
     if (y2<=0.01) { 
      cupt -> ignore = YES ;
      continue ;
     }
     cupt -> cauc_freq = y1/y2 ;

     cntit1(xc, cupt, indivmarkers, numindivs, 1) ;  
     y1 = 2*xc[2] + xc[1] ;  
     y2 = 2*asum(xc,3) ;    
     if (y2<=0.01) { 
      cupt -> ignore = YES ;
      continue ;
     }
     ++nsnp1 ; 
     cupt -> af_freq = y1/y2 ;
     y1 = cupt -> cauc_freq - cupt -> af_freq ;
     y = fabs(y1) ;    
     if (y<0.001) { 
      cupt -> ignore = YES ;
      continue ;
     }
     
     y = cupt -> cauc_freq + cupt -> af_freq ; 
     y *= 0.5 ;  
     y2 = y*(1.0-y) ; 
     if (y2<.001) {  
      cupt -> ignore = YES ;
      continue ;
     }
     ++nsnp2 ; 
     znorm = sqrt(y) ; // ?? involved in weight ?
     switch (runmode) { 
       case 0:  cupt -> weight = y1 / znorm ;
        break ;
       case 1:  cupt -> weight = y1 ; 
        break ; 
       case 2:  cupt -> weight = 1.0 ;
        break ;
       default: fatalx("bad runmode\n") ;
     }
     if (isnan(cupt -> weight)) fatalx("bad weight\n") ;
   }

   numsnps = rmsnps(snpmarkers, numsnps, NULL) ;   

   wtmean = wtnum = 0.0 ;
   for (x=0; x<numsnps; ++x) {  
     cupt = snpmarkers[x] ;
     if (cupt -> ignore) continue ; 
     if (isnan(cupt -> weight)) fatalx("bad weight\n") ;
     wtmean += cupt -> weight ;
     wtnum  += 1.0 ;
   }

   wtmean /= wtnum ;

// printf("raw weight mean:  %9.3f num: %9.0f\n", wtmean, wtnum) ;
   clipwt(snpmarkers, numsnps) ;

   if (weightname != NULL) {
    ncount = 0 ;
    wtmean = wtnum = 0.0 ;
    for (x=0; x<numsnps; ++x) {  
     cupt = snpmarkers[x] ;
     if (cupt -> ignore) continue ; 
     wtmean += cupt -> weight ;
     wtnum  += 1.0 ;
    }
    wtmean /= wtnum ;
    printf("weight mean:  %9.3f num: %9.0f\n", wtmean, wtnum) ;

  }


  if ((weightname == NULL) && (!ldmode)) {
   ncount = 0 ;
   for (x=0; x<numsnps; ++x) {  
     cupt = snpmarkers[x] ;

     cntit1(xc, cupt, indivmarkers, numindivs, 99) ;  
     y1 = 2*xc[2] + xc[1] ;  
     y2 = 2*asum(xc,3) ;    
     t = nnint(y1) ;
     if (t<mincount)  cupt -> ignore = YES ;
     t = nnint(y2-y1) ;
     if (t<mincount) cupt -> ignore = YES ;
//   if (x<1000) printf("%s %d %9.0f %9.0f %d\n", cupt -> ID, cupt -> ignore, y1, y2-y1, mincount) ;

     if (!cupt -> ignore) ++ncount ;
   }
  }

   ncount = 0 ;
   for (x=0; x<numsnps; ++x) {  
     cupt = snpmarkers[x] ;
     if (cupt -> weight == 0.0) cupt -> ignore = YES ;
     if (!cupt -> ignore) ++ncount ;
   }

   //printf("number of snps used:  %d\n", ncount) ;
   if (ncount < 100) fatalx("too few snps; numsnps: %d\n", numsnps) ;


  ncount = 0 ;
  if (runmode == 2) {
     yran = estyran(snpmarkers, numsnps, nxlim) ;
     printf("mean chi-sq: %9.3f random sample size: %d\n", yran, nxlim) ;
  }
  if (norun) { 
   printf("norun!\n") ;
   return 0 ;
  }
  printf("starting main loop. numsnps: %d\n", numsnps) ;
  for (x=0; x<numsnps; ++x)  { 
   cupt = snpmarkers[x] ;
   if (cupt -> ignore) continue ;
   chrom = cupt -> chrom ;
/**
   t = ranmod(1000) ;
   if (t==0) printf("processing snp %6d: %20s topheap: %x \n", x, cupt -> ID, topheap()) ;
*/
   fflush(stdout) ;
   ++snpdo1 ;
   for (x2=x+1; x2<numsnps; ++x2)  { 
    cupt2 = snpmarkers[x2] ;
    if (cupt2 -> ignore) continue ;
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      s = (int) (dis/binsize) ; // as in simpsim2 

      if (ransample > 0) {
       t = ranmod(ransample) ;
       if (t != 0) continue ;
      }

      ++snpdo2 ;

      if (runmode == 2)  {
       t = calcldscore(cupt, cupt2, &yy) ;
       if (t==0) continue ;
       yy /= yran ;
       wt = 1.0 ;
       xnum[s] += wt * yy ;  
       xden[s] += wt * wt ;
       yy2[s] += yy * yy ;
       corrpt = corrbins[s] ; 
       addcorr(corrpt, wt, yy) ;
       continue ;
      }
      wt = cupt -> weight * cupt2 -> weight ;
      if (wt == 0.0) continue ;
      ++ncount ;
//    if (ncount > 100) break ; 
      if ( (weightname == NULL) && (chithresh > 0.001)) {
//  set chithresh 0 if weightname NULL ?? 
 //    printf("qq1\n") ;
       y = cntit2(xc, cupt, cupt2, indivmarkers, numindivs, 0) ;  
       if (y<minparentcount) continue ;
       y1 = lddip(xc) ;  if (y1>=chithresh) continue ;
       y = cntit2(xc, cupt, cupt2, indivmarkers, numindivs, 1) ;  
       if (y<minparentcount) continue ;
       y2 = lddip(xc) ;  if (y2>=chithresh) continue ;
      }

      y = cntit2(xc, cupt, cupt2, indivmarkers, numindivs, 99) ;  
//    printf("qq2 %9.3f\n", y) ;
      if (y<4) continue ;
      cx1(xc, &y1, &y2) ;
/**
 bug ?? 
       if (ncount < -1) {
       printf("zz2 %9.3f %9.3f\n", y1, y2) ;
       printmat(xc, 3, 3) ;
      }
*/

      if (nnint(y1) < mincount) continue ;
      if (nnint(y2) < mincount) continue ;

      xww[0] = xww[2] = 0.5  ;
      xww[1] = 1.0  ;
      if (flatmode) addouter(xc, xww, 3) ;

//    t = ranmod(100) ; 
//    if (t==0) verbose = YES ;

      t = 999 ;

      if (runmode == 2) { 
       if (nnint(y1) < mincount) continue ;
       if (nnint(y2) < mincount) continue ;
       if (flatmode) vsp(xc, xc, 0.5, 9) ;
       yy = lddip(xc) ; 
       yy /= yran ;
      }

      else { 
       ++nloop ;
       yy = zdip(xc) ; 
      }

      ++nncount ;

//    if (nloop == -1) return 0 ;

      if (t==0) {       
       y1 = lddip(xc) ;
       printf("zzz %3d %12.6f %12.6f %12.6f\n", s, wt, yy, y1) ;
       printnl() ; 
       printmat(xc, 3, 3) ;
      }

      xnum[s] += wt * yy ;  
      xden[s] += wt * wt ;
      yy2[s] += yy * yy ;
      corrpt = corrbins[s] ; 
      addcorr(corrpt, wt, yy) ;
      ++jwt[chrom][s] ;
      for (k=1; k<= numchrom; ++k) {
        if (k==chrom) continue ;
        if (jackknife == NO) break ;
        corrpt = corrjbins[k][s] ;
        addcorr(corrpt, wt, yy) ;
      }
   }
  }

/**
 xden = sqrt(\sum wt^2)) 
 yy2  = sqrt(\sum z^2) 

 T = \sum(wt*z)  
 B = xden*yy2

 So ww = T/B is correlation <ww, z> 

*/

  vsp(xden, xden, 1.0e-10, numbins) ; 
  vsp(yy2, yy2, 1.0e-10, numbins) ; 

   if (runmode != 2) vsqrt(xden, xden, numbins) ;
   vsqrt(yy2, yy2, numbins) ;

   if (runmode != 2) {
    vvt(yy2, xden, yy2, numbins) ;
    vvd(ww, xnum, yy2, numbins) ;
    if (verbose) {
     printf("correlation (weight and yy)\n") ;    
     printmat(ww, 1, numbins-5) ;
    }
   }


// printf("num0: %d\n", num0) ;
   vvd(ww, xnum, xden, numbins) ;
   
   sprintf(sss, " ##Z-score and correlation:: %s  binsize: %12.6f", poplistname, binsize) ;
   dumpit(oname, ww, corrbins, numbins-5, binsize, sss) ;
   for (k=1; k<=numchrom; ++k) { 
    if (jackknife == NO) break ;
    sprintf(s1, "%s:%d", oname, k) ;
    sprintf(s2, "## Jackknife output: chrom %d", k) ;
    for (s=0; s<numbins; ++s) { 
     ww[s] = jwt[k][s] ;
    }
    dumpit(s1, ww, corrjbins[k], numbins-5, binsize, s2) ;
   }

  printf("## snps processes (outer + inner loop) %ld %ld\n", snpdo1, snpdo2) ;
  printf("##end of run\n") ;
  return 0 ;
}

void dumpit(char *oname, double *ww, CORR **corrbins, int len, double bsize, char *hdr) 
{
   double y, yreg ; 
   int k, ret ;
   FILE *fff ;
   CORR *corrpt ;

   if (oname == NULL) return ;
   openit(oname, &fff, "w") ;

   if (hdr != NULL) fprintf(fff, "%s\n", hdr) ;
   for (k=0; k<len; ++k) { 

    y = (k+1) * bsize ;
    corrpt = corrbins[k] ;
    fprintf(fff, "%9.3f ", 100.0*y) ;
    fprintf(fff, "%12.6f ", ww[k]) ;   
//  fprintf(fff, "%12.6f ", corrpt -> Z) ;
    if (runmode != 2) {
     calccorr(corrpt, 1, NO) ;
     yreg = corrpt -> S12 / (corrpt -> S11 + 1.0e-20) ;
     fprintf(fff, "%12.6f ", yreg) ; 
     fprintf(fff, "%12.6f ", corrpt -> corr) ;   
     fprintf(fff, "%12.0f ", corrpt -> S0) ;
    }
    fprintf(fff, "\n") ;

   }


   fclose(fff) ;

}
void
dorcx(double *ans, double *rr, double *dd, double binsize, double maxdis) 
// accumulate dot product into bins
{
   int x, x2, z, t, s, k ; 
   SNP *cupt, *cupt2 ;
   double dis ;
   double *xnum, *xd1, *xd2, y1, y2 ;
   double *rrr ;

   z =nnint(maxdis/binsize) ;  
   ZALLOC(xnum, z+1, double) ;
   ZALLOC(xd1, z+1, double) ;
   ZALLOC(xd2, z+1, double) ;

   vclear(xd1, 1.0e-10, z+1) ;
   vclear(xd2, 1.0e-10, z+1) ;

   for (x=0; x<numsnps; ++x)  {
    cupt = snpmarkers[x] ; 
    for (x2=x+1; x2<numsnps; ++x2)  {
      cupt2 = snpmarkers[x2] ; 
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      s = nnint(dis/binsize) ;
      y2 = dd[x]*dd[x2] ;
      for (k=0; k<numindivs; ++k) { 
       rrr = rr + k*numsnps ;

       switch (runmode) { 
        case 0:  
         y2 = dd[x]*dd[x2] ;
         y1 = rrr[x]*rrr[x2] ;
         break ;
      
        case 1:  
         y2 = rrr[x2]*dd[x2] ;
         y1 = rrr[x]*dd[x]  ;
         break ;

/**
        case 2:  
         y2  = 1.0 ;
         y1 = rrr[x]*rrr[x2]  ;
         break ;

        case 3:  
         y2  = 1.0 ;
         y1 = dd[x]*dd[x2] ;
         break ;
*/

        default: 
         fatalx("bad runmode\n") ; 
       }

       if (isnan(y1)) fatalx("bad y1\n") ;
       if (isnan(y2)) fatalx("bad y2\n") ;
       xnum[s] += y1*y2 ;
       xd1[s] += y1*y1 ;
       xd2[s] += y2*y2 ;
      }
    }
   }
  
   
   for (s=0; s <= z; ++s) {
    ans[s] = xnum[s]/sqrt(xd1[s]*xd2[s]) ;
   }

   free(xnum) ;
   free(xd1) ;
   free(xd2) ;
}

void
dorc2(double *ans, double *res, double binsize, double maxdis) 
// accumulate square of dot product into bins
{
   int x, x2, z, t ; 
   SNP *cupt, *cupt2 ;
   double dis, y ;

   z =nnint(maxdis/binsize) ;  
   vzero(ans, z+1) ;
   for (x=0; x<numsnps; ++x)  {
    cupt = snpmarkers[x] ; 
    for (x2=x+1; x2<numsnps; ++x2)  {
      cupt2 = snpmarkers[x2] ; 
      t = cupt2 -> chrom - cupt -> chrom ;
      if (t != 0) break ;
      dis = cupt2 -> genpos - cupt -> genpos ;
      if (dis >= maxdis) break ;
      z = nnint(dis/binsize) ;
      y  = res[x]*res[x2] ;
      ans[z] += y*y ;             
    }
   }
}

int readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
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

         
   if (parname == NULL) { 
    printf("no parameter file\n") ;
    return -1 ;
   }
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "weightname:", &weightname) ;
   getstring(ph, "admixpop:", &admixpop) ;
   getstring(ph, "admixlist:", &admixlist) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "dumpname:", &dumpname) ;
   getstring(ph, "loadname:", &loadname) ;
   getint(ph, "popestimate:", &popestimate) ;
   getint(ph, "ransample:", &ransample) ;
   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "nochrom:", &xnochrom) ;
   getstring(ph, "output:",   &oname) ;
   getint(ph, "runmode:", &runmode) ;
   getint(ph, "ldmode:", &ldmode) ;
   getint(ph, "seed:", &seed) ;
   getint(ph, "nxlim:", &nxlim) ;
   getint(ph, "mincount:", &mincount) ;
   getint(ph, "minparentcount:", &minparentcount) ;
   getint(ph, "norun:", &norun) ;
   getint(ph, "flatmode:", &flatmode) ;
   getint(ph, "zdipcorrmode:", &zdipcorrmode) ;
   getint(ph, "jackknife:", &jackknife) ;

   getdbl(ph, "maxdis:", &maxdis) ;
   getdbl(ph, "binsize:", &binsize) ;
   getdbl(ph, "chithresh:", &chithresh) ;
   getint(ph, "checkmap:", &checkmap) ;

   writepars(ph) ;
   closepars(ph) ;

   return 1 ;

}
void setfc(SNP **snpm, int numsnps) 
// also allocates fchrom , lchrom
{
  int i, j, chrom ;
  SNP *cupt ;
  double dis, totdis, pos1, pos2 ;

  ZALLOC(fchrom, 25, int) ;
  ZALLOC(lchrom, 25, int) ;

  ivclear(fchrom, 999999, 25) ;
  ivclear(lchrom, -999999, 25) ;
  
/* initialize real marker array */
  for (i=0; i<numsnps; i++) {
    cupt = snpm[i] ;
//  if (cupt -> isfake) continue ;
    if (cupt -> ignore) continue ;
    chrom = cupt -> chrom ;
//  if ((noxdata) && (chrom==23)) continue  ;
    fchrom[chrom] = MIN(fchrom[chrom],i) ;
    lchrom[chrom] = MAX(lchrom[chrom],i) ;
  }
}
double wynn(double *v, int n, double *acc, int *nacc)   
{
 double *x0, *x1, *xn ; 
 double t, amax, amin, tlast ;
 int iter = 0, j, nn  ;

 vmaxmin(v, n, &amax, &amin) ;  
 if (amax<=amin) {  
  vclear(acc, amax, n/2) ;
  *nacc = n/2 ;
  return amax ;
 }

 ZALLOC(x0, n, double) ; 
 ZALLOC(x1, n, double) ;
 ZALLOC(xn, n, double) ;
 copyarr(v, x1, n) ;  
 nn = n ;  
 for (;;) {  
  for (j=0; (j+1) < nn ; ++j) {  
   t = x0[j+1] + 1.0/(x1[j+1]-x1[j]) ;
   xn[j] = t ;
  }
  --nn ; 
  if (nn<2) break ;  
  copyarr(x1, x0, n) ;
  copyarr(xn, x1, n) ;

  for (j=0; (j+1) < nn ; ++j) {  
   t = x0[j+1] + 1.0/(x1[j+1]-x1[j]) ;
   xn[j] = t ;
  }
  --nn ; 
  if (nn<2) break ;  
  copyarr(x1, x0, n) ;
  copyarr(xn, x1, n) ;
  tlast = acc[iter] = t ; 
  ++iter ;
 }
 free(x0) ; free(x1) ; free(xn) ;
 *nacc = iter ;
 return tlast ;
}
double cntit1(double *xc, SNP *cupt, Indiv **indm, int numindivs, int t) 
{
  Indiv *indx ;
  int k, g ;

  vzero(xc, 3) ;
  if (cupt -> ignore) return 0;
  for (k=0; k<numindivs; ++k) { 
   indx = indm [k] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus  != t) continue ;
   g = getgtypes(cupt, k) ;
   if (g<0) continue ;
   ++xc[g] ;
  }
  return asum(xc,3) ;
}


double cntit2(double *xc, SNP *cupt, SNP *cupt2, Indiv **indm, int numindivs, int t) 
{
  Indiv *indx ;
  int k, e, f ;

  vzero(xc, 9) ;
  if (cupt -> ignore) return 0;
  if (cupt2 -> ignore) return 0;
  for (k=0; k<numindivs; ++k) { 
   indx = indm [k] ;
   if (indx -> ignore) continue ;
   if (indx -> affstatus  != t) continue ;
   e = getgtypes(cupt, k) ;
   if (e<0) continue ;
   f = getgtypes(cupt2, k) ;
   if (f<0) continue ;
   ++xc[3*e+f] ;
  }
  return asum(xc,9) ;

}
void clipwt(SNP **snpmarkers, int numsnps) 
{
  int k ;
  SNP *cupt ;
  double y ;

  for (k=0; k<numsnps; ++k) { 
   cupt = snpmarkers[k] ;
   y = cupt -> weight ;
   if (isnan(y)) fatalx("bad weight\n") ;
// printf("%12s %d %d %9.3f\n", cupt -> ID, k, cupt -> ignore, y) ;
   if (y< -1.0e10) cupt -> weight = 0.0 ;
   y = fabs(y) ;
   if (y< -1.0e-10) cupt -> weight = 0.0 ;
   if (k < -1) 
   printf("qq1 %12s %d %d %9.3f %9.3f\n", cupt -> ID, k, cupt -> ignore, y, cupt -> weight) ;
  }

}
double estyran(SNP **snpmarkers, int numsnps, int nxlim) 
{
   int nx, x, x2, t ; 
   double ychi, yran, y ;  
   double y1, y2 ;
   SNP *cupt, *cupt2 ;
   double xc[9] ; 
   int nloop = 0 ;

// random choice
   nx = 0; 
   yran = 0; 

   for (;;) {
    if (xchrom > 0) break ;
    if (nxlim == 0) break ;
    ++nloop ; 
    if ((nloop > 1000) && (nx == 0)) fatalx("looping\n") ;
    x = ranmod(numsnps) ;
    x2 = ranmod(numsnps) ;
    cupt = snpmarkers[x] ; 
    cupt2 = snpmarkers[x2] ; 
    if (cupt -> ignore) continue ;
    if (cupt2 -> ignore) continue ;
    t = cupt2 -> chrom - cupt -> chrom ;
    if (t == 0) continue ;
    t = calcldscore(cupt, cupt2, &ychi) ;
    if (t==0) continue ;

    yran += ychi  ;
    ++nx ;
    if (nx >= nxlim) break ;
   }
   if (nx <= 0) return 1.0 ;
   return yran / (double) nx ;
 
}
int calcldscore(SNP *cupt, SNP *cupt2, double *xscore)
{
  double xc[9] ;
  double y1, y2 ;

  if (cupt -> ignore) return 0  ;
  if (cupt2 -> ignore) return 0  ;
    cntit2(xc, cupt, cupt2, indivmarkers, numindivs, 99) ;  
    cx1(xc, &y1, &y2) ;
    if (nnint(y1) < mincount) return 0 ;
    if (nnint(y2) < mincount) return 0 ;
    vsp(xc, xc, 0.5, 9) ;
    *xscore = lddip(xc) ;
    return 1 ;

}
void cx1(double *xc, double *y1, double *y2) 
{
   double x1[3], x2[3] ;
   int e, f ;

   vzero(x1,3) ;
   vzero(x2,3) ;
   for (e=0; e<3; ++e) {  
    for (f=0; f<3; ++f) {  
     x1[e] += xc[3*e+f] ;
     x2[f] += xc[3*e+f] ;
    }
   }
  
   *y1 = cminor(x1) ;
   *y2 = cminor(x2) ;
    
}

double cminor(double *xx) 
{
  double y1, y2 ;
  y1 = xx[1]+2*xx[2] ; 
  y2 = 2*asum(xx, 3) - y1 ;
  return MIN(y1, y2) ;
}
int cmap(SNP **snppmarkers, int numsnps) 
{
   int t, k ; 
   double y1, y2 ; 
   SNP *cupt ;  
   for (k=1 ; k<=10; ++k) { 
    t = ranmod(numsnps) ;
    cupt = snpmarkers[t] ; 
    y1 = cupt -> genpos ; 
    y2 = cupt -> physpos / 1.0e8  ; 
    if (fabs(y1-y2) > .001) return YES ;
   }
   return NO ;
}

