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

#define WVERSION   "570" 

/** 
 reads data, basic stats
 good place for oddball counts
 hetcount added
 malexhet set yes
 snpmissing indmissing snpmf
 pubgeno added 
 treemix output
 countmisspop added
*/


#define MAXFL  50   
#define MAXSTR  512

extern int packmode ;

double tiny = .01 ;
double small = .1 ;
int lopos = -1 ; 
int hipos = 1000*1000*1000 ;  

char *treename = NULL ;
int treemode = NO ;
char *trashdir = "/var/tmp" ;
//int verbose = NO ;
int qtmode = NO ;
int printvalids = NO ;
int snpcounts = NO ;
int downsamp = -1 ;
int seed = 0 ;

int snpmftest = YES ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int killmono = YES ;
int pubgeno = NO ; 

int xhets = NO ;
int checkonly = NO ;
int doatrend = NO ;
int maxchrom = 999 ; 
int xchrom = -1 ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *goodsnpname = NULL ;
char *badpairsname = NULL ;
char *poplistname = NULL ;

char *outputname = NULL ;
FILE *ofile ;

double thresh = 16.0 ;
double epithresh = 6.0 ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
double lddip(double *xxx) ;
void doxhets(Indiv **indm, int numind) ;
void doychrom(Indiv **indm, int numind) ;
double atrend(SNP *cupt, Indiv ** indm, int numind) ;
void cntitcasev(double *xc, SNP *c1, SNP *c2, Indiv **indm, int numind,  int val) ;
void cntit1(double *xc, SNP *c1) ;
void countmissingsnp(SNP *cupt, int *pahit, int *pamiss) ;
void countmissingind(Indiv  *indx, int *pahit, int *pamiss) ; 
void countmissingpop(char *pop, int *pahit, int *pamiss)  ;
void countmissingpopall(char **poplist, int npops, int *pahit, int *pamiss) ; 
void allelecnts(SNP *cupt, int *pn0, int *pn1) ;
double cntitmf(SNP *cupt) ;
double z2x2(double *a) ;
void writetreemix(char *treename, SNP **snpmarkers, int numsnps, 
 Indiv **indivmarkers, int numindivs, char **poplist, int npops) ;

char **eglist ;
int numeg = 0 ;

int main(int argc, char **argv)
{

  int x, i,j,n,nv,e,f,t, n0, n1; 
  SNP *cupt ;
  SNP *cupt1, *cupt2 ;
  Indiv *indx ;
  int numvind  ;

  int nindiv = 0 ;
  int nignore, numrisks = 1 ;
  double y, cc[3], xxx[9], ychi, xx1[9], xx2[9] ; 
  int xcc[3]  ;
  double p1, p2, pdiff ;
  double *eq, **rhs ; 
  double *w1, *ee ; 
  double ans1[9], ans2[9], y1, y2, yl1, yl2 ;
  int numstates, ignore ;
  int g, k, c0[2], c1[2] ;
  int amin, amax, xmin, xmax ;
  int ahit, amiss, xhit, xmiss ;
  int chrom ;
  int nvalsnpa, nvalsnpx, nvalind ;
  int malexhet = YES ;  // don't change data
  int ncase ; 
  int nmono, npoly ;
  double ymem ;

  cputime(0) ;
  calcmem(0) ;

  ofile = stdout; 
  packmode = NO ;
  readcommands(argc, argv) ;

  printf ("## easystats version: %s\n", WVERSION);
  if (downsamp) {    
   snpcounts = YES ;
   printf("downsample: %d\n", downsamp) ;
   if (seed==0) seed = seednum() ; 
   printf("seed: %d\n", seed) ;
   SRAND(seed) ;
  }
  if (treename != NULL) { 
   treemode = YES ;
  }
  if (treemode && (treename == NULL)) treename = strdup("Treedata") ;
  printvers(argv[0], WVERSION) ;
  if (outputname != NULL) openit(outputname, &ofile, "w") ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  numeg = 0 ; 
  if (poplistname != NULL) 
  { 
    ZALLOC(eglist, numindivs, char *) ; 
    numeg = loadlist(eglist, poplistname) ;
    seteglist(indivmarkers, numindivs, poplistname);
    for (i=0; i<numindivs; ++i)  {     
     indx = indivmarkers[i] ; 
     if (indx -> affstatus == NO) indx -> ignore = YES ;
    }
  }
  else  {
   ncase = setstatus(indivmarkers, numindivs, "Case") ;
   if (ncase == 0) doatrend = 0 ;
  }

  setgenotypename(&genotypename, indivname) ;

  printf("genotypename:  %s\n", genotypename) ;

  if (genotypename != NULL)  {
   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  }

 if (treemode) {
 writetreemix(treename, snpmarkers, numsnps, 
  indivmarkers, numindivs, eglist, numeg) ;

   printf("##end of run\n") ; 
   return 0 ;
 }

  numvind = numvalidind(indivmarkers, numindivs) ;

/**
  for (i=0; i<numsnps; i++) { 
   cupt = snpmarkers[i] ;  
   if (cupt -> ignore) continue ;
   for (k=0; k<numindivs; ++k) { 
    indx = indivmarkers[k] ;
    if (indx -> ignore) putgtypes(cupt, k, -1) ;
   }
  }
*/
  nmono = npoly = nignore = 0 ;
  for (i=0; i<numsnps; i++) { 
   cupt = snpmarkers[i] ;  
   t = nnint(cupt -> physpos) ;
   if (t<lopos) cupt -> ignore = YES ;
   if (t>hipos) cupt -> ignore = YES ;
   if (cupt -> chrom > maxchrom) cupt -> ignore = YES ;
   if ((xchrom > 0) && (cupt -> chrom != xchrom) ) cupt -> ignore = YES ;
   if (cupt -> ignore) continue ; 
// if (cupt -> chrom > 23) cupt -> ignore = YES ; 
   allelecnts(cupt, &n0, &n1) ;
   if (MIN(n0, n1) == 0) {
     if (killmono) cupt -> ignore = YES ;
     if (MAX(n0, n1) > 0) ++nmono ;
   }
   else ++npoly ; 
   if (cupt -> ignore) ++nignore ;
  }
  printf("\n\n") ;
  y1 = (double) nmono / (double) (npoly + nmono) ;
  printf("numindivs: %d valid: %d numsnps: %d nignore: %d monomorphic: %d poly: %d monorate: %10.4f\n" ,
    numindivs, numvind, numsnps, nignore, nmono, npoly, y1) ; 

  for (i=0; i<numsnps; i++) { 
   if (xhets) break ;
   cupt = snpmarkers[i] ;  
   if (cupt -> ignore) continue ;
   if (cupt -> chrom > 22) continue ;
   y = hwcheckx(cupt, indivmarkers, cc)  ;  
   fixit(xcc, cc, 3) ;

   printf("hw %20s %3d %9.3f %12.0f %4d %4d %4d   %9.3f\n",  
    cupt -> ID, cupt -> chrom, cupt -> genpos, cupt -> physpos,  
    xcc[0], xcc[1], xcc[2], y) ;

  }
  for (i=0; i<numindivs; i++) { 
   ivzero(c1, 2) ; 
   ivzero(c0, 2) ;   
   ivzero(xcc, 3) ;
   indx = indivmarkers[i] ;
   if (indx -> ignore) continue ;
   for (j=0; j<numsnps; j++) {  
    cupt = snpmarkers[j] ; 
//  gettln(cupt, indx, NULL,  NULL, &numstates, &ignore) ; 
    if (cupt -> ignore) continue ;
    if (indx -> ignore) continue ;
//  if (ignore) continue ; 
//  if (numstates != 3) continue ;  
    g = getgtypes(cupt, i) ;  
    if (g<0) continue ;
    ++xcc[g] ;
    k = 0 ;
    if (cupt -> chrom == 23) k=1 ;
    if (g==1) ++c1[k] ;  
    if (g!=1) ++c0[k] ;     
   }
   printf("hetcount %20s ", indx -> ID) ;  
   printf(" %c ", indx -> gender) ;
   for (k=0; k<=1; ++k)  {  
    p1 = c1[k] ; 
    p2 = p1 + c0[k] + 1.0e-8 ;
    printf("  %5d %5d", c1[k], nnint(p2)) ;
    printf("  %9.3f", p1/p2) ;
   }
   printf(" :: ") ;
   printimat(xcc, 1, 3) ;
  }
  fflush(stdout) ;

  //fastdupcheck(snpmarkers, indivmarkers, numsnps, numindivs) ;  
  if (xhets) {
     doxhets(indivmarkers, numindivs) ;
     doychrom(indivmarkers, numindivs) ;
   }

  amax = xmax = -1 ; 
  amin = xmin = numsnps ;
  nvalsnpa = nvalsnpx = nvalind = 0 ;
  if (printvalids)  { 
    printf("#snpmissing:%20s %3s %12s" , "ID", "chrom", "physpos") ;             
    printf(" %6s %6s", "hit", "miss") ;
    printf(" %12s", "miss frac") ; 
    printnl () ;
   for (i=0; i<numsnps; i++) { 
    cupt = snpmarkers[i] ;  
    if (cupt -> ignore) continue ;
    chrom = cupt -> chrom ;  
//  if (chrom > 23) cupt -> ignore = YES ;
    if (chrom<23) ++nvalsnpa ;          
    if (chrom==23) ++nvalsnpx ;          
//  k = numvalidgtypes(cupt) ;
    countmissingsnp(cupt, &ahit, &amiss) ;
    printf("snpmissing: %20s %3d %12.0f" , cupt ->ID, cupt -> chrom, cupt -> physpos) ;             
    printf(" %6d %6d", ahit, amiss) ;
    y1 = amiss ; y2 = ahit + amiss ;
    printf(" %12.6f", y1/(y2+1.0e-8)) ;
    cntit1(cc, cupt) ;
    for (j=0; j<3; ++j) { 
     t = nnint(cc[j]) ; 
     printf(" %6d", t) ;
    }
    printnl() ;
   }
   printnl() ;
   printnl() ;
   printf("snpmftest\n") ;
   for (i=0; i<numsnps; i++) { 
    if (!killmono) break ;
    cupt = snpmarkers[i] ;  
    if (cupt -> ignore) continue ;
    y = cntitmf(cupt) ;
    printf("snpmf: %20s %3d %12.0f" , cupt ->ID, cupt -> chrom, cupt -> physpos) ;             
    printf(" %12.3f", y) ;
    printnl() ;
   }
   printnl() ;
   printnl() ;
    printf("#indmissing:%20s"  , "ID") ; 
    printf(" %6s %6s", "hit", "miss") ;
    printf(" %12s", "miss frac") ; 
    printnl () ;
   for (i=0; i<numindivs; i++) { 
    indx = indivmarkers[i] ;
    if (indx -> ignore) continue ;
    ++nvalind ;
    countmissingind(indx, &ahit, &amiss) ;
    printf("indmissing: %20s", indx ->ID) ;
    printf(" %6d %6d", ahit, amiss) ;
    y1 = amiss ; y2 = ahit + amiss ;
    printf(" %12.6f", y1/(y2+1.0e-8)) ;
    printf("  %20s", indx -> egroup) ;
    printnl() ;
   }
   printf("nvalsnpa: %d  nvalsnpx: %d  nvalidind:  %d\n", nvalsnpa, nvalsnpx, nvalind) ;
   printnl() ;
   printnl() ;
  }
   printnl() ;
   printnl() ;
  if (snpcounts) {

   for (i=0; i<numsnps; i++) { 
    cupt = snpmarkers[i] ;  
    if (cupt -> ignore) continue ;
    allelecnts(cupt, &n0, &n1) ;
    if (downsamp > 0) { 
     t = ranhprob(n0+n1, n1, downsamp) ; 
     if (t < 0) continue ; 
     n0 = downsamp - t ; 
     n1 = t ; 
    }
    printf("snpcounts: %20s %3d %12.0f" , cupt ->ID, cupt -> chrom, cupt -> physpos) ;             
    printf(" %c", cupt -> alleles[0]) ;
    printf(" %c", cupt -> alleles[1]) ;
    printf(" %5d %5d", n0, n1) ;
    printnl() ;
   }

  }
   printnl() ;
   printnl() ;
  if (snpmftest) { 
   printf("snpmftest\n") ;
   for (i=0; i<numsnps; i++) { 
    if (!killmono) break ;
    cupt = snpmarkers[i] ;  
    if (cupt -> ignore) continue ;
    y = cntitmf(cupt) ;
    printf("snpmf: %20s %3d %12.0f" , cupt ->ID, cupt -> chrom, cupt -> physpos) ;             
    printf(" %12.3f", y) ;
    printnl() ;
   }
  }


  fflush(stdout) ;
  t = 9999999 ;
  for (k=0; k<numeg; ++k) { 
   countmissingpop(eglist[k], &ahit, &amiss) ;
   printf("popmisssing: %20s ", eglist[k]) ; 
    printf(" %6d %6d", ahit, amiss) ;
    t = MIN(t, ahit) ;
    y1 = amiss ; y2 = ahit + amiss ;
    printf(" %12.6f", y1/(y2+1.0e-8)) ;
    printnl() ;
  }
  printnl() ;

   countmissingpopall(eglist, numeg, &ahit, &amiss) ;
   printf("popmisssingall: %20s ", "All" ) ; 
   printf(" %6d %6d", ahit, amiss) ;
   y1 = amiss ; y2 = ahit + amiss ;

  printf(" %12.6f", y1/(y2+1.0e-8)) ;
  if (t==0) printf(" *** some pop is all missing ***\n") ;
  printnl() ;
  fflush(stdout) ;

  if (pubgeno) { 
   for (i=0; i<numindivs; ++i) { 
    indx = indivmarkers[i] ;
    printf("pubg: %s %c %s ", indx -> ID,indx -> gender, indx -> egroup) ;
    for (j = 0; j< numsnps; ++j) { 
      cupt = snpmarkers[j] ; 
      if (cupt -> ignore) continue  ;
      g = getgtypes(cupt, i) ;  
      if (g<0) g = 9 ; 
      printf("%d", g) ;
    }
    printnl() ;
   }
  }


 if (treemode) { 
  writetreemix(treename, snpmarkers, numsnps, 
   indivmarkers, numindivs, eglist, numeg) ;
 }

  if (checkonly) {
   ymem = calcmem(1)/1.0e6 ;
   printf("##end of easystats: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
   return 0 ;
  }

  for (i=0; i<numsnps; i++) { 
   if (xhets) break ;  
   for (j=i+1; j<numsnps; j++) { 
    cupt1 = snpmarkers[i] ;  
    cupt2 = snpmarkers[j] ;  

    if (cupt1 -> chrom != cupt2 -> chrom) break ;

    if (cupt1 -> ignore) continue ;
    if (cupt1 -> chrom > 22) continue ;
    if (cupt2 -> ignore) continue ;
    if (cupt2 -> chrom > 22) continue ;

    cntit(xxx, cupt1, cupt2) ;
    ychi = lddip(xxx) ;
    if (fabs(ychi)<thresh) continue ;  
    p1 = cupt1 -> physpos ;
    p2 = cupt2 -> physpos ;
    pdiff = p2 - p1 ;
    printf("chi %20s %20s %9.3f %10.0f %10.0f  %10.0f",  cupt1->ID, cupt2 -> ID, ychi, p1, p2, pdiff)  ;
    if (pdiff > 300000) printf ("  +++") ;
    printf("\n") ;
    printmat(xxx, 3, 3) ;
   }
  }
  for (i=0; i<numsnps; i++) { 
   if (doatrend == NO) break ;
   if (!killmono) break ;
   cupt = snpmarkers[i] ;  
   if (cupt -> ignore) continue  ;
   ychi = atrend(cupt, indivmarkers, numindivs) ;
   cupt -> score = ychi ; 
   p1 = cupt -> physpos ;
   printf("atrend:  %20s %10.0f   %9.3f\n", cupt -> ID, p1, ychi) ;
  }


  ZALLOC(eq, 9*numindivs, double) ;
  rhs =  initarray_2Ddouble(9, 2, 0.0) ;
  

  for (i=0; i<numsnps; i++) { 
   if (xhets) break ;  
   for (j=i+1; j<numsnps; j++) { 
    cupt1 = snpmarkers[i] ;  
    cupt2 = snpmarkers[j] ;  

    if (cupt1 -> ignore) continue ;
    if (cupt1 -> chrom > 22) continue ;
    if (cupt2 -> ignore) continue ;
    if (cupt2 -> chrom > 22) continue ;

    cntitcasev(xx1, cupt1, cupt2, indivmarkers, numindivs, YES) ;
    cntitcasev(xx2, cupt1, cupt2, indivmarkers, numindivs, NO) ;
    vsp(xx1, xx1, small,  9) ;
    vsp(xx2, xx2, small,  9) ;
/**
    printf("\n\n") ;
    printmat(xx1,3, 3) ;
    printf("\n\n") ;
    printmat(xx2,3, 3) ;
    printf("\n\n") ;
*/
    nv = 5 ;  
    n = 0 ;
    ZALLOC(w1, 9, double) ;
    for (e=0; e<3; e++)  {  
     for (f=0; f<3; f++)  {   // sum right = sum left
      y1 = xx1[3*e+f] ;
      y2 = xx2[3*e+f] ;
      ee = eq + n*nv ; 
      rhs[n][0] = y1 ;
      rhs[n][1] = y2 ;
      if (e<2) {
       ee[e] = 1.0 ;  
      }
      else {ee[0] = ee[1] = -1 ;}
      if (f<2) {
       ee[2+f] = 1.0 ;
      }
      else {ee[2] = ee[3] = -1 ; }   
      ee[4] = 1.0 ;
      ++n ;
/**
      gaussa(w1, nv) ;
      vst(w1, w1, tiny, nv) ;
      vvp(ee, ee, w1, nv) ;
*/
     }
    }
    yl1 = logregressit(ans1, eq, rhs, n, nv)  ;
    nv = 9 ;
    n = 0 ;
    vzero(eq, 9*numindivs) ;
    for (e=0; e<3; e++)  {  
     for (f=0; f<3; f++)  {   // sum right = sum left
      y1 = xx1[3*e+f] ;
      y2 = xx2[3*e+f] ;
      ee = eq + n*nv ; 
      rhs[n][0] = y1 ;
      rhs[n][1] = y2 ;
      ee[3*e+f] = 1.0 ;  
      ++n ;
/**
      gaussa(w1, nv) ;
      vst(w1, w1, tiny, nv) ;
      vvp(ee, ee, w1, nv) ;
*/
     }
    }
    yl2 = logregressit(ans2, eq, rhs, n, nv)  ;
    y = 2*(yl2-yl1) ;
    if (y<epithresh) continue ;
/**
    printf("yl1: %9.3f\n", yl1) ;
    printmat(ans1, 1, 5) ;
    printf("yl2: %9.3f\n", yl2) ;
    printmat(ans2, 1, 9) ;
*/
    y1 = cupt1 -> score ;
    y2 = cupt2 -> score ;
    printf("epi %20s %20s %9.3f %9.3f    %9.3f\n", cupt1 -> ID, cupt2 -> ID, y1, y2, y) ;  
    printmat(xx1, 3, 3) ;
    printf("\n") ;
    printmat(xx2, 3, 3) ;
    printf("\n") ;

   }
  }
   ymem = calcmem(1)/1.0e6 ;
   printf("##end of easystats: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

  printf("##end of run\n") ;
  return 0 ;
}
void writetreemix(char *treename, SNP **snpmarkers, int numsnps, Indiv **indivmarkers, int numindivs, char **poplist, int npops) 
{
  FILE *fff ;
  int ***counts , **cp;
  int i, k, nrows, ncols, a, b, tmin ;
  SNP **xsnplist ; 
  int *xindex, *xtypes ;
  Indiv **xindlist, *indx ;
  char ss[1024], *sx ;

  if (treename == NULL) return ; 
  
  openit(treename, &fff, "w") ;
 

  ZALLOC (xindex, numindivs, int);
  ZALLOC (xindlist, numindivs, Indiv *);

  nrows = loadindx (xindlist, xindex, indivmarkers, numindivs);
  ZALLOC (xtypes, nrows, int);
  for (i = 0; i < nrows; i++) {
    indx = xindlist[i];
    k = indxindex (poplist, npops, indx->egroup);
    xtypes[i] = k;
  }

  ZALLOC (xsnplist, numsnps, SNP *);
  ncols = loadsnpx (xsnplist, snpmarkers, numsnps, indivmarkers);

  printf ("snps: %d  indivs: %d\n", ncols, nrows);

  ZALLOC (counts, ncols, int **);
  for (k = 0; k < ncols; ++k) {
    counts[k] = initarray_2Dint (npops, 2, 0);
  }

  countpops (counts, xsnplist, xindex, xtypes, nrows, ncols);
  for (k=0; k<npops; ++k) { 
    fprintf(fff, "%s ", poplist[k]) ;
  }
  fprintf(fff,"\n") ; 
  for (i=0; i<ncols; i++) { 
   cp = counts[i] ;
   tmin = 999999 ;
   sx = ss ;
   for (k=0; k<npops; ++k) {   
    a=cp[k][0] ; 
    b=cp[k][1] ; 
    tmin = MIN(tmin, a+b) ; 
    sx += sprintf(sx, "%d,%d", a,b) ;  
    if (k<(npops-1)) sx += sprintf(sx, " ") ;     
   } 
   if (tmin>0) fprintf(fff, "%s\n", ss) ;


  }

 fclose(fff) ;

}
double atrend(SNP *cupt, Indiv ** indm, int numind) 
{

 double *xa, *xb ;
 Indiv *indx ;
 int cc[4], g, gc, n, i, xmin ; 
 double rho, ychi, yn ;

 if (cupt -> ignore) return 0.0 ;

 ZALLOC(xa, numind, double) ;
 ZALLOC(xb, numind, double) ;

 n = 0 ;
 ivzero(cc, 4) ;
 for (i=0; i<numind; i++) { 
  g = getgtypes(cupt, i) ;
  if (g<0) continue ;
  gc = 2-g ; 
  cc[0] += g ; 
  cc[1] += gc ;
  indx = indm[i] ; 
  if ((indx -> gender == 'M') && (cupt -> chrom == 23)) g /= 2 ;
  if (indx-> ignore) continue ; 
  xb[n] = indx -> affstatus ;
  cc[2] += xb[n] ; 
  cc[3] += 1-xb[n] ;
  xa[n] = g ;  
  ++n ;  
 }
/**
 printf("zz %20s ", cupt -> ID) ;
 printimat(cc, 1, 4) ;
*/

 yn = (double) n ;
 ivmaxmin(cc, 4, NULL, &xmin) ; 
 ychi = 0.0 ;
 if (xmin>0) {  
  rho = corr(xa, xb, n) ; 
  ychi = yn*rho*rho ; 
 }

 free(xa) ; 
 free(xb) ;
 return ychi ;

}
void doychrom(Indiv **indm, int numind) 
{
  int fc = numsnps+100, lc = -1 ;
  int i, k, g ;
  SNP *cupt ;
  Indiv *indx ;
  double hsum , hnum ;
  double msum , mnum ;

  for (k=0; k<numsnps; k++) { 
   cupt = snpmarkers[k] ;
   if (cupt -> chrom != 24) continue ;
   fc = MIN(fc, k) ; 
   lc = MAX(lc, k) ;
  }
  if (lc<0) return ;
   for (i=0; i<numind; i++) { 
    indx = indm[i] ;
    hsum = hnum = 0 ;
    msum = mnum = 0 ;
    for (k=fc; k<=lc; ++k) {  
     cupt = snpmarkers[k] ;
     g = getgtypes(cupt, i) ;
     ++mnum ; 
     if (g<0) { 
       ++msum ;
       continue ;  
     }
     ++hnum ; 
     if (g==1) { 
       ++hsum ;
     }
    }
    mnum += .001 ;
    hnum += .001 ;
    printf("ychrom %20s %c %4d %4d  %9.3f", indx->ID, indx -> gender, 
       nnint(hsum), nnint(hnum), hsum/hnum) ;
    printf("   %4d %4d  %9.3f", 
       nnint(msum), nnint(mnum), msum/mnum) ;
    printf(" %20s", indx -> egroup) ;
    printnl() ;
   }
}

void doxhets(Indiv **indm, int numind) 
{
  int fc = numsnps+100, lc = -1 ;
  int i, k, g ;
  SNP *cupt ;
  Indiv *indx ;
  double hsum , hnum ;
  int ismale ;
  int *badss, **ggg, **fff ;

  ZALLOC(badss, numsnps, int) ;
  ggg = initarray_2Dint(numsnps, 3, 0) ;
  fff = initarray_2Dint(numsnps, 3, 0) ;
  for (k=0; k<numsnps; k++) { 
   cupt = snpmarkers[k] ;
   if (cupt -> chrom != 23) continue ;
   fc = MIN(fc, k) ; 
   lc = MAX(lc, k) ;
  }
   for (i=0; i<numind; i++) { 
    indx = indm[i] ;
    if (indx -> ignore) continue ;
    ismale = NO ;
    if (indx -> gender == 'M') ismale = YES ;
    hsum = hnum = 0 ;
    for (k=fc; k<=lc; ++k) {  
     cupt = snpmarkers[k] ;
     g = getgtypes(cupt, i) ;
     if (g<0) continue ;  
     ++hnum ; 
     if (ismale) ++ ggg[k][g] ;
     if (!ismale) ++ fff[k][g] ;
     if (g==1) { 
       ++hsum ;
       if (ismale) {
        ++badss[k] ;
       }
     }
    }
    if (hnum<0.1) continue ;
    printf("hetx %20s %c %4d %4d  %9.3f", indx->ID, indx -> gender, 
       nnint(hsum), nnint(hnum), hsum/hnum) ;
    printf(" %20s", indx -> egroup) ;
    printnl() ;
   }
    for (k=fc; k<=lc; ++k) {  
     cupt = snpmarkers[k] ;
     printf("badhetsnpcnt: %20s %12.0f %6d   ", cupt -> ID, cupt -> physpos, badss[k]) ;
     for (g=0; g<3; ++g) { 
      printf(" %6d", ggg[k][g]) ;
     }
     printf("   ") ; 
     for (g=0; g<3; ++g) { 
      printf(" %6d", fff[k][g]) ;
     }
     printnl() ;
   }
   free(badss) ;
   free2Dint(&ggg, numsnps) ;
   free2Dint(&fff, numsnps) ;
}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vVt")) != -1) {

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

      case 't':
	treemode = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("%s: parameter file: %s\n", argv[0], parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

/**
DIR2:  /fg/nfiles/admixdata/ms2
SSSS:  DIR2/outfiles 
genotypename: DIR2/autos_ccshad_fakes
eglistname:    DIR2/eurlist  
output:        eurout
*/
   getint(ph, "packmode:", &packmode) ; // controls internals 

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "genooutfilename:", &genooutfilename) ;
   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "output:", &outputname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "goodsnpname:", &goodsnpname) ;
   getstring(ph, "badpairsname:", &badpairsname) ;                 
   getstring(ph, "poplistname:", &poplistname) ;
   getint(ph, "xhets:", &xhets) ;                 
   getint(ph, "checkonly:", &checkonly) ;                 
   getint(ph, "printvalids:", &printvalids) ;                 
   getint(ph, "snpcounts:", &snpcounts) ;                 
   getint(ph, "downsamp:", &downsamp) ;                 
   getint(ph, "downsample:", &downsamp) ;                 
   getint(ph, "seed:", &seed) ;                 
   getint(ph, "snpmftest:", &snpmftest) ;                 
   getint(ph, "doatrend:", &doatrend) ;                 
   getint(ph, "killmono:", &killmono) ;                 
   getint(ph, "pubgeno:", &pubgeno) ;                 
   getint(ph, "lopos:", &lopos) ;                 
   getint(ph, "hipos:", &hipos) ;                 
   getint(ph, "maxchrom:", &maxchrom) ;                 
   getint(ph, "chrom:", &xchrom) ;                 
   getstring(ph, "treemix:", &treename) ;
   writepars(ph) ;
   closepars(ph) ;

}

void cntitcasev(double *xc, SNP *c1, SNP *c2, Indiv **indm, int numind,  int val) 
{
 int n, i, e, f ;
 Indiv *indx ;
 n = MIN(c1->ngtypes, c2->ngtypes) ;
 vzero(xc, 9) ;
 for (i=0; i<numind; i++) {  
  indx = indm[i] ;
  if (indx->affstatus != val) continue ;
  if (indx-> ignore) continue ;
  e = getgtypes(c1, i) ;
  f = getgtypes(c2, i) ;
  if (e<0) continue ;
  if (f<0) continue ;
  ++xc[3*e+f] ;
 }
}
void allelecnts(SNP *cupt, int *pn0, int *pn1)
{
  double xc[3] ;

  *pn0 = *pn1 = 0 ;
  if (cupt -> ignore) return ;
  cntit1(xc, cupt) ;
  *pn0 = nnint(2*xc[0] + xc[1]) ;
  *pn1 = nnint(2*xc[2] + xc[1]) ;
}
void cntit1(double *xc, SNP *cupt) 
{
 int n, i, e ;
 Indiv *indx ;

 vzero(xc, 3) ;
 if (cupt -> ignore) return ;
 for (i=0; i<numindivs; i++) {  
  indx = indivmarkers[i] ;  
  if (indx -> ignore) continue ;
  e = getgtypes(cupt, i) ;
  if (e<0) continue ;
  ++xc[e] ;
 }
}
double cntitmf(SNP *cupt) 
{
 int g, i, e ;
 double y, xm[3], xf[3], cc[4] ; 
 Indiv *indx ;
 static int ncall = 0 ;

 if (cupt -> ignore) return 0.0 ;
 vzero(xm, 3) ;
 vzero(xf, 3) ;
 for (i=0; i<numindivs; i++) {  
  indx = indivmarkers[i] ;  
  if (indx -> ignore) continue ;
  g = getgtypes(cupt, i) ;
  if (g<0) continue ;
  if (indx -> gender == 'M') ++xm[g] ;
  else ++xf[g] ;
 }
 if (cupt -> chrom == 23) vst(xm, xm, 0.5, 3) ;
 cc[1] = (xm[1] + 2*xm[2]) ;
 cc[0] = (xm[1] + 2*xm[0]) ;
 cc[3] = (xf[1] + 2*xf[2]) ;
 cc[2] = (xf[1] + 2*xf[0]) ;
 ++ncall ;
 vsp(cc, cc, 1.0e-8, 4) ;
 y =  z2x2(cc) ;      // positive if females have higher freq than males 
 if (ncall < -1) { 
   printmat(cc, 2, 2) ;
   printf("%12.3f ", y) ;
 }
 return y ;
}
void
countmissingsnp(SNP *cupt, int *pahit, int *pamiss) 
{
   int k, ahit, amiss, g ;
   Indiv *indx ;

   *pahit = *pamiss = 0 ;
   if (cupt -> ignore) return ;
   ahit = amiss = 0 ;
   for (k=0; k<numindivs; ++k) { 
    indx = indivmarkers[k] ; 
    if (indx -> ignore) continue ;
     g = getgtypes(cupt, k) ;
     if (g<0) ++amiss ;
     else ++ahit ;
   }
   *pahit = ahit ;
   *pamiss = amiss ;
   return ;

}
void
countmissingind(Indiv  *indx, int *pahit, int *pamiss) 
{
   int k, ahit, amiss, g, n ;
   SNP *cupt ;  

   *pahit = *pamiss = 0 ;
   if (indx -> ignore) return ;
   ahit = amiss = 0 ;
   n = indx -> idnum ;
   for (k=0; k<numsnps; ++k) { 
    cupt = snpmarkers[k] ; 
    if (cupt -> ignore) continue ;
     g = getgtypes(cupt, n) ;
     if (g<0) ++amiss ;
     else ++ahit ;
   }
   *pahit = ahit ;
   *pamiss = amiss ;
   return ;

}
void
countmissingpopall(char **poplist, int npops, int *pahit, int *pamiss) 
{
   int j, k, ahit, amiss, g, miss1, n, t, ishit, ismiss ;
   Indiv *indx ; 
   char *pop ; 
   SNP *cupt ;  
   int *mm, *hh ;

   *pahit = *pamiss = 0 ;
   ahit = amiss = 0 ;
   for (k=0; k<numsnps; ++k) { 
    cupt = snpmarkers[k] ; 
    if (cupt -> ignore) continue ;
    miss1 = NO ; 
    for (j=0;  j< npops; ++j) { 
    pop = poplist[j] ;
    ishit = NO ; 
    ismiss = YES ; 
    for (n=0; n<numindivs; ++n) {
     indx = indivmarkers[n] ; 
     t = strcmp(indx -> egroup, pop) ;
     if (t != 0) continue ;
     g = getgtypes(cupt, n) ;
     if (g>=0) ismiss = NO ;
    }
    if (ismiss == YES) {
     miss1 = YES ;
     break ; 
    }
   }
   if (miss1 == YES) ++amiss ;
   else ++ ahit ;
  }

   *pahit = ahit ;
   *pamiss = amiss ;
   return ;

}
 
void
countmissingpop(char *pop, int *pahit, int *pamiss) 
{
   int k, ahit, amiss, g, n, t, ishit ;
   Indiv *indx ; 
   SNP *cupt ;  

   *pahit = *pamiss = 0 ;
   ahit = amiss = 0 ;
   for (k=0; k<numsnps; ++k) { 
    cupt = snpmarkers[k] ; 
    if (cupt -> ignore) continue ;
    ishit = NO ; 
    for (n=0; n<numindivs; ++n) {
     indx = indivmarkers[n] ; 
     t = strcmp(indx -> egroup, pop) ;
     if (t != 0) continue ;
     g = getgtypes(cupt, n) ;
     if (g>=0) ishit = YES ; ;
    }
    if (ishit) ++ahit ; 
    else ++amiss ;
   }

   *pahit = ahit ;
   *pamiss = amiss ;
   return ;

}
