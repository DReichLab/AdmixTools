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

#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  

#define WVERSION   "2501" 
#define MAXFL  50   
#define MAXSTR  512

// phasedmode added
// bugfix for pedcols
// dup indivs allowed if 1 is ignore
// O2 version

char *trashdir = "/var/tmp" ;
extern int verbose  ;
int qtmode = NO ;
Indiv **indm1, **indm2 ;  
SNP **snpm1, **snpm2 ;
SNP **outsnps ; 
Indiv **indivmarkers ;
int numsnps, numindivs ;
int nums1, nums2 ; 
int numi1, numi2 ; 

char  *snp1 = NULL ;
char  *ind1 = NULL ;
char  *geno1 = NULL ;

char  *snp2 = NULL ;
char  *ind2 = NULL ;
char  *geno2 = NULL ;

char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;
char *badsnpname = NULL ;
char *deletesnpoutname = NULL ;

int packout = -1 ;
int tersem  = YES ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
char *omode = "packedancestrymap" ;
extern int packmode ;
int ogmode = NO ;
int docheck = YES ;
int strandcheck = YES ;
int phasedmode = NO ;

int xchrom = -1 ;
int lopos = -999999999 ; 
int hipos = 999999999 ;
int minchrom = 1 ; 
int maxchrom = 97 ;
extern long rlen, packlen ; 

/** 
 docheck  YES  allele flipping check
 hashcheck  YES  ... NO => snp names, Indiv names changed MUST retain order and number
 strandcheck (default YES) if NO then alleles are assumed on same strand
*/

char  unknowngender = 'U' ;

void setomode(enum outputmodetype *outmode, char *omode)  ;
void readcommands(int argc, char **argv) ;
int checkmatch(SNP *cupt1, SNP *cupt2) ;

int 
mergeind(SNP **outsnps, Indiv **indivmarkers,
  int numsnps, int numindivs, int *x1, int *x2)   ;


int main(int argc, char **argv)
{
  unsigned char *packg1, *packg2 ;

  int **snppos ;
  int *snpindx ;
  int  lsnplist, lindlist, numeg ;
  int i,j; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n1, n2 ;
  int nkill = 0 ;
  int t, k, x ;
  int *x1, *x2 ; 

  int nignore, numrisks = 1 ;

  char **genolist ;
  int numgenolist ;
  int maxmiss ; 
  double ymem ; 

  tersem = YES ;     // no snp counts

  readcommands(argc, argv) ;

  setomode(&outputmode, omode) ;
  packmode = YES ;
  settersemode(tersem) ;

  cputime(0) ;
  calcmem(0) ;

  nums1 = 
    getsnps(snp1, &snpm1, 0.0, NULL, &nignore, numrisks) ;

  putped(1) ;
  freeped() ;

  nums2 = 
    getsnps(snp2, &snpm2, 0.0, NULL, &nignore, numrisks) ;

  putped(2) ;
  freeped() ;

  numi1 = getindivs(ind1, &indm1) ;
  numi2 = getindivs(ind2, &indm2) ;
  t = MIN(numi1, numi2) ;
  ZALLOC(x1, t, int) ;
  ZALLOC(x2, t, int) ;
  ivclear(x1, -1, t) ;
  ivclear(x2, -1, t) ;

   numindivs = 0 ;
   ZALLOC(indivmarkers, numi1 + numi2, Indiv *) ;
   for (x=0; x<numi1; ++x)  {  
// intersection
    indx = indm1[x] ;  
    k = indindex(indm2, numi2, indx -> ID) ;
    if (k<0) indx -> ignore = YES ;
    if (indx -> ignore != YES) {  
     indivmarkers[numindivs] = indx ;
     x1[numindivs] = x ;
     x2[numindivs] = k ;
     ++numindivs ;
    }
   }
// x1 index in ind1, x2 in ind2 

 t = nums1 + nums2 ; 
 ZALLOC(outsnps, t, SNP *) ;

  for (i=0; i<t; i++) {
    ZALLOC(outsnps[i], 1, SNP) ;
    cupt = outsnps[i] ;
    clearsnp(cupt) ;
    ZALLOC(cupt -> modelscores, 1, double) ;
    ZALLOC(cupt -> totmodelscores, 1, double) ;
  }
  numsnps = 0 ;

  for (x=0; x<nums1; ++x)  {  
   cupt1 = snpm1[x] ;
   if (cupt1 -> ignore) continue ;
   outsnps[numsnps] = cupt1 ;  
   cupt1 -> isrfake = 1 ;
   ++numsnps ;
  }

  for (x=0; x<nums2; ++x)  {  
   cupt2 = snpm2[x] ;
   if (cupt2 -> ignore) continue ;
   k = snpindex(snpm1, nums1, cupt2 -> ID) ;  
   if (k>=0) continue ;
   outsnps[numsnps] = cupt2 ;  
   cupt2 -> isrfake = 2 ;
   ++numsnps ;
  }


  freesnpindex() ;
  setgenotypename(&geno1, ind1) ;
  getped(1) ;
  getgenos(geno1, snpm1, indm1, 
     nums1, numi1, nignore) ;

  packg1 = (unsigned char *) getpackgenos() ;
  clearpackgenos() ;

  setgenotypename(&geno2, ind2) ;
  getped(2) ;
  getgenos(geno2, snpm2, indm2, 
     nums2, numi2, nignore) ;

 mergeind(outsnps, indivmarkers, numsnps, numindivs, x1, x2)    ;
 printf("## numsnps: %d    %d %d\n", numsnps, nums1, nums2) ;

 sortsnps(outsnps, outsnps, numsnps) ;

  outfiles (snpoutfilename, indoutfilename, genooutfilename,
    outsnps, indivmarkers, numsnps, numindivs, packout, ogmode);


  ymem = calcmem(1)/1.0e6 ;
  printf("##end of snpunion: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0 ;
}
int checkmatch(SNP *cupt1, SNP *cupt2) 
{

  char a1, a2, b1 , b2 ;

  if (docheck == NO) return 1 ;

  a1 = cupt1 -> alleles[0] ;
  a2 = cupt1 -> alleles[1] ;

  b1 = cupt2 -> alleles[0] ;
  b2 = cupt2 -> alleles[1] ;

  a1 = toupper(a1) ;
  a2 = toupper(a2) ;

  b1 = toupper(b1) ;
  b2 = toupper(b2) ;

  if ((a1 == 'X') && (a2 == 'X') && (strandcheck)) return -1 ;  // flipcheck impossible

  if (strandcheck) {

   if ((a1 == 'A') && (a2 == 'T')) return -1 ;
   if ((a1 == 'T') && (a2 == 'A')) return -1 ;
   if ((a1 == 'C') && (a2 == 'G')) return -1 ;
   if ((a1 == 'G') && (a2 == 'C')) return -1 ;

   if ((b1 == 'A') && (b2 == 'T')) return -1 ;
   if ((b1 == 'T') && (b2 == 'A')) return -1 ;
   if ((b1 == 'C') && (b2 == 'G')) return -1 ;
   if ((b1 == 'G') && (b2 == 'C')) return -1 ;
  }

  if ((a1 == b1) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b2 ;
   return 1 ;
  }

  if ((a1 == b2) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b1 ;
   return 2 ;
  }


  if ((a1 == b1) && (a2 == b2)) return 1 ;
  if ((a1 == b2) && (a2 == b1)) return 2 ;

  if ((a1 == b1) && (b2 =='X')) return 1 ;
  if ((a2 == b1) && (b2 =='X')) return 2 ;

  if (strandcheck == NO) return 0 ;

  b1 = compbase(b1) ;
  b2 = compbase(b2) ;

  if ((a1 == b1) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b2 ;
   return 1 ;
  }

  if ((a1 == b2) && (a2 == 'X')) {           
   cupt1 -> alleles[1] = b1 ;
   return 2 ;
  }

  if ((a1 == b1) && (a2 == b2)) return 1 ;
  if ((a1 == b2) && (a2 == b1)) return 2 ;

  if ((a1 == b1) && (b2 =='X')) return 1 ;
  if ((a2 == b1) && (b2 =='X')) return 2 ;

  return 0 ;

}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:vVf")) != -1) {

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

      case 'f':
	phasedmode = YES ;                      
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring(ph, "geno1:", &geno1) ;
   getstring(ph, "snp1:", &snp1) ;
   getstring(ph, "ind1:", &ind1) ;

   getstring(ph, "geno2:", &geno2) ;
   getstring(ph, "snp2:", &snp2) ;
   getstring(ph, "ind2:", &ind2) ;

   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "genooutfilename:", &genooutfilename) ; 
   getstring(ph, "outputformat:", &omode) ;
   getstring(ph, "deletesnpoutname:", &deletesnpoutname) ;

   getint(ph, "docheck:", &docheck) ;
   getint(ph, "hashcheck:", &hashcheck) ;
   getint(ph, "strandcheck:", &strandcheck) ;
   getint(ph, "phasedmode:", &phasedmode) ;

   writepars(ph) ;
   closepars(ph) ;

}
int 
mergeind(SNP **outsnps, Indiv **indivmarkers,
  int numsnps, int numindivs, int *x1, int *x2)   
{
   SNP *cupt1, *cupt2, *cupt ;
   int k, x, g, a, i, t  ;
   double y ;
   static unsigned char *packg, *buff ;
   int *bb ;

   ZALLOC(bb, numindivs, int) ;
   y = (double) (numindivs * 2) / (8 * (double) sizeof (char)) ;
   rlen = nnint(ceil(y)) ;
   rlen = MAX(rlen, 48)  ;
   packlen = numsnps*rlen ;
   ZALLOC(packg, packlen, unsigned char) ;
   clearepath((char *) packg) ;
  
    buff = packg ;
    for (x=0; x<numsnps; ++x) {  
     cupt = outsnps[x] ;
     cupt -> ignore = NO ;
     a = cupt -> isrfake ;
     for (i=0; i<numindivs; ++i) {  
      if (a==1) g = getgtypes(cupt, x1[i]) ; 
      if (a==2) g = getgtypes(cupt, x2[i]) ; 
      bb[i] = g ;
     }
     cupt -> pbuff = (char *) buff ;
     cupt -> isrfake = NO ;
     cupt -> ngtypes = numindivs ;
     for (i=0; i<numindivs; ++i) {  
      g = bb[i] ;
      putgtypes(cupt, i, g) ;
     }
     buff += rlen ;
    }
   free(bb) ;
   return numindivs ;
}
