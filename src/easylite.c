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

#include "admutils.h"
#include "mcio.h"  
#include "regsubs.h" 
#include "egsubs.h" 
#include "qpsubs.h" 

#define WVERSION   "550" 

/** 
 reads data, basic stats simplified easystats
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
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int killmono = NO ;

int maxchrom = 999 ; 
int xchrom = -1 ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;

FILE *ofile ;

double thresh = 16.0 ;
double epithresh = 6.0 ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
double lddip(double *xxx) ;
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

  int i,j,n,nv,e,f,t, n0, n1; 
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
  int nhit, nmiss ;

  ofile = stdout; 
  packmode = NO ;
  readcommands(argc, argv) ;

//printvers(argv[0], WVERSION) ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  numeg = 0 ; 
   ncase = setstatus(indivmarkers, numindivs, "Case") ;

  setgenotypename(&genotypename, indivname) ;

  printf("genotypename:  %s\n", genotypename) ;

   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  numvind = numvalidind(indivmarkers, numindivs) ;

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

  for (i=0; i<numindivs; i++) { 
   ivzero(c1, 2) ; 
   ivzero(c0, 2) ;   
   ivzero(xcc, 3) ;
   nmiss = nhit = 0 ; 
   indx = indivmarkers[i] ;
   if (indx -> ignore) continue ;
   for (j=0; j<numsnps; j++) {  
    cupt = snpmarkers[j] ; 
//  gettln(cupt, indx, NULL,  NULL, &numstates, &ignore) ; 
    if (cupt -> ignore) continue ;
//  if (numstates != 3) continue ;  
    g = getgtypes(cupt, i) ;  
    if (g<0) { 
     ++nmiss ;
     continue ;
    }
    ++xcc[g] ;
    ++nhit ;
    k = 0 ;
    if (cupt -> chrom == 23) k=1 ;
    if (g==1) ++c1[k] ;  
    if (g!=1) ++c0[k] ;     
   }
   printf("indcount: %20s ", indx -> ID) ;  
   printf(" %c ", indx -> gender) ;
   printimatx(xcc, 1, 3) ;
   printf(" missing: %6d\n", nmiss) ;
  }
  printf("\n\n") ; 

   for (j=0; j<numsnps; j++) {  
    ivzero(c1, 2) ; 
    ivzero(c0, 2) ;   
    ivzero(xcc, 3) ;
    nmiss = nhit = 0 ; 
    cupt = snpmarkers[j] ; 
    if (cupt -> ignore) continue ;
      for (i=0; i<numindivs; i++) { 
      indx = indivmarkers[i] ;
      if (indx -> ignore) continue ;
      g = getgtypes(cupt, i) ;  
      if (g<0) { 
       ++nmiss ;
       continue ;
      }
     ++xcc[g] ;
     ++nhit ;
    }
   printf("snpcount: %20s ", cupt -> ID) ;  
   printf(" %2d ", cupt -> chrom) ;
   printf(" %9.0f :: ", cupt -> physpos) ;
   printimatx(xcc, 1, 3) ;
   printf(" missing: %6d\n", nmiss) ;
  }
  fflush(stdout) ;


  printf("##end of easylite\n") ;
  return 0 ;
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

   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getint(ph, "killmono:", &killmono) ;                 
   getint(ph, "lopos:", &lopos) ;                 
   getint(ph, "hipos:", &hipos) ;                 
   getint(ph, "maxchrom:", &maxchrom) ;                 
   getint(ph, "chrom:", &xchrom) ;                 
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
