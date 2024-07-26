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

#define WVERSION   "165" 

/** 
 reads data, basic checks ; hashnode hashes genotypes for each sample 
 genohash option
*/


#define MAXFL  50   
#define MAXSTR  512


typedef struct {
  char ID[IDSIZE];
  char *egroup ;
  char gender;     /* 'M' or 'F' */
  int idnum ;
  int ignore ;        /* YES => do not use */
  int droplo ; 
  int drophi ; 
  int hashets ;
  int isdg ;
  int issg ;
  int **gcount ; 
  int *dropx ;
} Indz; 

int dropdetails = NO ;

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
Indiv **indivmarkers;
Indz  **indzarr, *indz ;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int killmono = NO ;
int pubgeno = NO ; 

int xhets = NO ;
int checkonly = NO ;
int doatrend = NO ;
int maxchrom = 22 ;  // autosomes only
int xchrom = -1 ;

int hashmode = YES ;
int genohash = YES ;
int misshash = YES ;

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *genooutfilename = NULL ;
char  *indoutfilename = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;

char *outputname = NULL ;
FILE *ofile ;

double thresh = 16.0 ;
double epithresh = 6.0 ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;

void allelecnts(SNP *cupt, int *pn0, int *pn1) ;
void printshort(char *ss, int len) ;
int qisdg(Indz *indz) ;
int ckdrop(Indz *indz) ;
void printdropdetails(Indz *indz) ;
int hashlongarr(long *a, int len)  ;
int hashintarr(int *a, int len)  ;



char **eglist ;;
int numeg = 0 ;
void checkduppos(SNP **snpm, int nsnps) ;
void splitlong(unsigned long x, int *pa, int *pb) ;

int main(int argc, char **argv)
{

  int i,j,n,nv,e,f,t, tt, n0, n1; 
  int u, v ; 
  SNP *cupt ;
  SNP *cupt1, *cupt2 ;
  Indiv *indx ;
  int numvind  ;
  int *hasharr, hashlen ,  zhash, zlen, maxhashlen  ;
  long *lhasharr, zlhash ; 


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
  int ahit, amiss, xhit ;
  int chrom ;
  int nvalsnpa, nvalsnpx, nvalind ;
  int malexhet = YES ;  // don't change data
  int ncase ; 
  int nmono, npoly ;
  int nhit, nmiss ;
  double *xmiss, ymiss, yhit, ymthresh ;
  double ymem ; 
  int lo, hi ;

  ofile = stdout; 
  packmode = NO ;
  readcommands(argc, argv) ;

  cputime(0) ;
  calcmem(0) ;

  if ((genohash == NO) && (misshash==NO)) hashmode = YES  ;



//printvers(argv[0], WVERSION) ;
  setfancyhash(NO) ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;

  numindivs = getindivs(indivname, &indivmarkers) ;
  numeg = 0 ; 
   ncase = setstatus(indivmarkers, numindivs, "Case") ;
   if (ncase == 0) doatrend = 0 ;

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
  checkduppos(snpmarkers, numsnps) ;
  printf("\n\n") ;
  y1 = (double) nmono / (double) (npoly + nmono) ;
  printf("numindivs: %d valid: %d numsnps: %d nignore: %d monomorphic: %d poly: %d monorate: %10.4f\n" ,
    numindivs, numvind, numsnps, nignore, nmono, npoly, y1) ; 

  ZALLOC(xmiss, numsnps, double) ;
  ZALLOC(indzarr, numindivs, Indz *) ;
  t = numsnps/32 ; t += 32 ; 
  maxhashlen = t ; 
  ZALLOC(hasharr, t, int) ;
  ZALLOC(lhasharr, t, long) ;
  for (i=0; i<numindivs; i++) { 
   ZALLOC(indzarr[i], 1, Indz) ;
   ivzero(c1, 2) ; 
   ivzero(c0, 2) ;   
   ivzero(xcc, 3) ;
   nmiss = nhit = 0 ; 
   indx = indivmarkers[i] ;

   indz = indzarr[i] ;   
   strcpy(indz -> ID, indx -> ID) ;
   indz -> ignore = indx -> ignore ;
   indz -> egroup = strdup(indx -> egroup) ;
   indz -> ignore = indx -> ignore ;
   indz -> gender = indx -> gender ;
   indz -> idnum = indx -> idnum ;
   indz -> droplo = indz -> drophi = -1 ;
   indz -> hashets = NO ;
   qisdg(indz) ;
   indz -> gcount = initarray_2Dint(maxchrom+1, 2, 0) ;  
   ZALLOC(indz -> dropx, maxchrom+1, int) ;

   if (indx -> ignore) continue ;
   hashlen = zlen = zhash = 0  ;  
   zlhash = 0 ; 
   for (j=0; j<numsnps; j++) {  
    cupt = snpmarkers[j] ; 
//  gettln(cupt, indx, NULL,  NULL, &numstates, &ignore) ; ;
    if (cupt -> ignore) continue ;
//  if (numstates != 3) continue ;  
    chrom = cupt -> chrom ; 
    tt = g = getgtypes(cupt, i) ;  
    t = 0 ; 
    if (g<0) t = 1 ;
    if (tt<0) tt = 3;
    zhash = (zhash << 1) + t ; 
    zlhash = (zlhash << 2) + tt ; 
    zlen += 1 ;  
    if (zlen == 32) { 
     zlen = 0 ; 
     hasharr[hashlen] = zhash ; 
     lhasharr[hashlen] = zlhash ; 
     ++hashlen ; 
     if (hashlen >= maxhashlen) fatalx("hasharr overflow\n") ; 
     zhash = 0 ;
     zlhash = 0 ;
    }
    if (g<0) { 
     ++nmiss ;
     ++xmiss[j] ;  
     ++indz -> gcount[chrom][1] ; 
     continue ;
     
    }
    ++xcc[g] ;
    ++nhit ;
    ++indz -> gcount[chrom][0] ; 
    k = 0 ;
    if (cupt -> chrom == 23) k=1 ;
    if (g==1) {
     ++c1[k] ;  
     indz -> hashets = YES ;
    }
    if (g!=1) ++c0[k] ;     
   }
   printf("indcount: %20s ", indx -> ID) ;  
   printf(" %c ", indx -> gender) ;
    if (misshash == YES) { 
     zhash = hashintarr(hasharr, hashlen) ; 
     printf(" %08X ", zhash) ;
    } 
    if (genohash == YES) { 
     zhash = hashlongarr(lhasharr, hashlen) ; 
     printf(" G%08X ", zhash) ;
    } 
   printimatx(xcc, 1, 3) ;
   printf(" missing: %6d", nmiss) ;
   printshort(indz -> egroup, 20) ;
   if ((indz -> isdg) && (indz -> hashets == NO)) {
    printf(" *** baddg") ; 
   }
   if ((indz -> issg) && (indz -> hashets == YES)) {
    printf(" *** badsg") ; 
   }
    ckdrop(indz) ;
    if (indz -> droplo > 0) {
      printf(" *** chrom. dropout ") ;
      printdropdetails(indz) ;
    } 
    printnl() ;
  }
  printf("\n\n") ; 

   ymiss = asum(xmiss, numsnps) / (double) numsnps ; 
   yhit = 1.0 - ymiss ;
   ymthresh =  yhit / 5.0 ;
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
   printf(" missing: %6d", nmiss) ;
   if ( (double) nhit < ymthresh) printf(" *** himiss") ;
   printnl() ;
  }
  fflush(stdout) ;

   printnl() ;
  ymem = calcmem(1)/1.0e6 ;
  printf("##end of easycheck: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

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
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getint(ph, "lopos:", &lopos) ;                 
   getint(ph, "hipos:", &hipos) ;                 
   getint(ph, "maxchrom:", &maxchrom) ;                 
   getint(ph, "chrom:", &xchrom) ;                 
   getint(ph, "dropdetails:", &dropdetails) ;                 
   getint(ph, "genohash:", &genohash) ;                 
   getint(ph, "misshash:", &misshash) ;                 
   writepars(ph) ;
   closepars(ph) ;

}

int qisdg(Indz *indz) 
{

 int xt = 0 ;  

if (strstr(indz -> ID, ".DG")) xt = 1 ;    
if (strstr(indz -> ID, ".HO")) xt = 1 ;    
if (strstr(indz -> egroup, ".DG"))  xt = 1 ;  

if (strstr(indz -> ID, ".SG")) xt = 2 ;    
if (strstr(indz -> ID, ".BY")) xt = 2 ;    
if (strstr(indz -> ID, ".WGC")) xt = 2 ;    
if (strstr(indz -> ID, ".AA")) xt = 2 ;    
if (strstr(indz -> ID, ".EC")) xt = 2 ;    
if (strstr(indz -> ID, ".REF")) xt = 2 ;    
if (strstr(indz -> egroup, ".SG"))  xt = 2 ;  

 indz -> isdg = indz -> issg = NO ;
 if (xt == 1) indz -> isdg = YES ;
 if (xt == 2) indz -> issg = YES ;

 return xt ;

}
void printshort(char *ss, int len) 
{
 char sss[MAXSTR] ;  

 strncpy(sss, ss, len) ; 
 sss[len] = CNULL ; 
  
 printf (" %s", sss) ;

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


void allelecnts(SNP *cupt, int *pn0, int *pn1)
{
  double xc[3] ;

  *pn0 = *pn1 = 0 ;
  if (cupt -> ignore) return ;
  cntit1(xc, cupt) ;
  *pn0 = nnint(2*xc[0] + xc[1]) ;
  *pn1 = nnint(2*xc[2] + xc[1]) ;
}

int ckdrop(Indz *indz) 
{
   int lo=999, hi=0 ;
   int gc[2], k, a0 ; 
   double yx, y0, y1 ;
   int nthresh ;

   ivzero(gc, 2) ;
   indz -> droplo = indz -> drophi = 0 ;
   for (k=1; k <= maxchrom; ++k) { 
    ivvp(gc, gc, indz -> gcount[k], 2) ;
   }
   y0 = gc[0] ; 
   y1 = gc[1] ; 
   for (k=1; k <= maxchrom; ++k) { 
    yx = intsum(indz -> gcount[k], 2) ; 
    yx = (yx*y0)/(y0+y1) ;
    nthresh = nnint(yx/5.0) ; 
    if (nthresh <= 10) continue ;
    a0 = indz -> gcount[k][0] ; 
    if (a0 <= nthresh) {  
     lo = MIN(lo, k) ;
     hi = MAX(hi, k) ;
     indz -> dropx[k] = YES ;
    }

   }

   if (hi > 0) { 
    indz -> droplo = lo ;
    indz -> drophi = hi ;
    return 1 ;
   }

   return 0 ;

}

void printdropdetails(Indz *indz) 
{
  int k ; 
 if (dropdetails == NO) { 
  printf(" %d-%d ", indz -> droplo, indz -> drophi) ; 
  return ;
 }

 for (k=1; k<= maxchrom; ++k) { 
  if (indz -> dropx[k] == NO) continue ;
  printf (" %d %d ; ", k, indz -> gcount[k][0]) ;
 }
 return ;

} 


int hashlongarr(long *a, int len) 
{
 int thash, k ;  
 int *xx, u, v ; 

 ZALLOC(xx, 2*len, int) ;

 for (k=0; k<len; ++k) { 
  splitlong((unsigned long) a[k], &u, &v) ;  
  xx[k] = u ; 
  xx[k+len] = v ;
 }


 thash = hashintarr(xx, 2*len) ; ;  


 free(xx) ; 

 return thash ; 


}
int hashintarr(int *a, int len) 
{
 int thash, k ;  

 thash = len ;  

 for (k=0; k<len; ++k) {  
  thash += xhash1(a[k]) ; 
  thash = xcshift(thash, 3) ;
 }


 return thash ; 


}
void pspos(SNP *cupt) 
{

  printf(" %s %d %12.0f ", cupt -> ID, cupt -> chrom, cupt -> physpos) ;

}
void checkduppos(SNP **snpm, int nsnps) 
{
  SNP *cupt, *cupt2 ;
  int k, p1, p2 ; 

  for (k=1; k< nsnps; ++k) { 
   cupt = snpm[k] ; 
   cupt2 = snpm[k-1] ; 
   p1 = nnint(cupt -> physpos) ;
   p2 = nnint(cupt2 -> physpos) ;
   if (p1 != p2) continue ;
   if (cupt -> chrom != cupt2 -> chrom) continue ;
   printf("*** duplicate position: ") ;
   pspos(cupt2) ; pspos(cupt) ; printnl() ;
  }
   return ;
}
