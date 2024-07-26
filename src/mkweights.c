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
 Can't find .tex file for theory  ; but pdf in ccpaper
*/

#define WVERSION   "100"
// was qpadmlin
// inbreed added

#define MAXFL  50   
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL ;
char *trashdir = "/var/tmp" ;

int qtmode = NO ;

int colcalc = YES ;
int inbreed = NO ;

Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int markerscore = NO ;
int seed = 0 ;
int noxdata = NO ;  // default as why not output X weights

double plo = .001 ;
double phi = .999 ;
double pvhit = .001 ;
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


char *outputname = NULL ;
char *weightname = NULL ;
FILE *wfile ;
char **eglist ;

void readcommands(int argc, char **argv) ;
void calcwts(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char *outpop) ;

int main(int argc, char **argv)
{

  int i, j, k, k1, k2, k3, k4, kk; 
  SNP *cupt ; 
  Indiv *indx ;
  double y1, y2, y, sig, tail, yy1, yy2 ;
  int chrom ; 
  double ymem ; 


  int nignore, numrisks = 1 ;

  readcommands(argc, argv) ;
  printf("## mkweights version: %s\n", WVERSION) ;

  cputime(0) ;
  calcmem(0) ;

  SRAND(77) ;  // dealing with hets and inbreed YES 

  if (parname == NULL) return 0 ;

  setinbreed(inbreed) ;

  printf("## weightoutname: %s\n", weightname) ;

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
    if ((noxdata) && (chrom == (numchrom + 1))) cupt->ignore = YES;
    if (chrom == 0)
      cupt->ignore = YES;
    if (chrom > (numchrom + 1))
      cupt->ignore = YES;
    if (chrom == zchrom)
      cupt->ignore = YES;
  }

   for (i = 0; i < numsnps; i++) {
     cupt = snpmarkers[i];
     cupt -> weight = 1.0 ; 
   }

   calcwts(snpmarkers, indivmarkers, numsnps, numindivs, outpop) ;
   if (weightname != NULL) openit(weightname, &wfile, "w") ; 
   else wfile = stdout ; 
   for (i = 0; i < numsnps; i++) {
     cupt = snpmarkers[i];
     if (cupt -> ignore) continue ;
     fprintf(wfile, "%20s %9.3f\n", cupt -> ID, cupt -> weight) ;
   }
  ymem = calcmem(1)/1.0e6 ;
  printf("##end of mkweight: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

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
   getstring(ph, "popname:", &outpop) ;
   getstring(ph, "pop:", &outpop) ;
   getstring(ph, "weightoutname:", &weightname) ;
   getstring(ph, "weightname:", &weightname) ; // deprecated 
   getint(ph, "inbreed:", &inbreed) ; 

   getint(ph, "chrom:", &xchrom) ;  

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);

}

int isxm(SNP *cupt, Indiv *indx) 

{
 if (cupt -> chrom != (numchrom+1)) return NO ; 
 if (indx -> gender != 'M') return NO ; 

 return YES ; 

}

void setab(int *pa, int *pb, int g, SNP *cupt, Indiv *indx)  
{
 int a, b ; 
 *pa = *pb = 0 ;
 
 if (g<0) return ; 
 a = g/2 ; b = 2-a ; 

 if (inbreed || isxm(cupt, indx)) { 
  if (g==1) {  
   a = ranmod(2) ;
   b = 1 - a ; 
  }
  else { 
    a = g/2 ; 
    b = 1-a ;
  }
 }
  *pa = a ; 
  *pb = b ; 
  return ;
}


void
calcwts(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char *outpop) 
{
 int *indarr ; 
 SNP *cupt ; 
 Indiv *indx ;
 int i, j, k,  nsamp, t ;
 int tt, bb ; 
 int a, b, na, nb, g ; 
 double ysum, ynum, ymul, ymean ;

 if (outpop == NULL) return ; 
 ZALLOC(indarr, numindivs, int) ;

 nsamp = 0 ; 
 for (i=0 ; i< numindivs; ++i) { 
  indx = indivmarkers[i] ;
  t = strcmp(indx -> egroup, outpop) ;
  if (t!=0) continue ; 
  indarr[nsamp] = i ; 
  ++nsamp ;
 }
 if (nsamp == 0) fatalx("no samples with label %s!\n", outpop) ;
 ysum = ynum = 0 ;
 for (k=0; k<numsnps; ++k) { 
  cupt = snpmarkers[k] ; 
  cupt -> weight = 1.0 ; 
  if (cupt -> ignore) continue ; 
  na = nb = 0 ; 
  for (i=0; i< nsamp; ++i) { 
   j = indarr[i] ; 
   indx = indivmarkers[j] ;
   g = getgtypes(cupt, j) ; 
   setab(&a, &b, g, cupt, indx) ;
   na += a ; 
   nb += b ; 
  }
     
  bb = (na+1) * (nb+1) ; 
  tt = (na+nb+2) * (na + nb + 1) ; 

  cupt -> weight = (double) tt / (double) bb ; // exp of 1/(p*(1-p)) with prior Beta 2, 2 
  ysum += cupt -> weight ;
  ++ynum ; 
//  printf("zzdebug %d %d %d %d %d %d\n", k, cupt -> chrom, a, b, na, nb) ;
 }

 ymean = ysum/ynum ; 
 for (k=0; k<numsnps; ++k) { 
  cupt = snpmarkers[k] ; 
  if (cupt -> ignore) continue ; 
  cupt -> weight /= ymean ;  // mean is now 1 
 }

 free(indarr) ; 
 return ; 
}
