#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>

#include "admutils.h"
#include "mcio.h"  
#include "mcmcpars.h"  
#include "egsubs.h"  

#define WVERSION   "100:3021" 
/** 
 reformats files.             
 pedfile junk (6, 7 cols, ACGT added)
 packped (.bed) added  input and output 
 pedfile can now be in random order
 tersemode supported 
 ped -> ancestrymap outputs alleles
 various fixups for ped files in computing variance, reference alleles 
 alleles can (say) be A X if second allele unknown
 allele flipped for bed format
 huge ascii file trapped now
 fastdup added
 r2thresh added
 killr2 added to convertf (default is NO)  
 bug in processing C0 SNPS in ped files fixed

 outputall added  (all indivs, snps output)
 genolistname added (list of packed ancestrymap files one per line)
 support for enhanced eigenstrat format and genolist can also be enhanced eigenstrat
 Chromosome # 1-89 now supported.  Improved treatment of funny chromosomes. 
 MT -> 90 
 XY -> 91

 Various bugfixes and checks for files size  
 inpack2 no longer needs big memory 
 hash table lookup for snpindex
 check size more carefully
 update to mcio.c  (snprawindex)
 maxmissfrac added

 poplistname added 
 remap (newsnpname) added
 flipsnpname, flipreference added  
 flipsnpname flips genotype values 
 if flipped we can flip reference allele too if we wish (default YES)
 badpedignore.  Crazy bases flagged as ignore
 remapind (newindivname) added
*/


#define MAXFL  50   
#define MAXSTR  512

char *trashdir = "/var/tmp" ;
int verbose = NO ;
int qtmode = NO ;
Indiv **indivmarkers, **indm2;
SNP **snpmarkers ;
SNP **snpm2 ;

int numsnps, numindivs, numind2 ; 
int nums2 ;

char  *genotypename = NULL ;
char  *genotypelist = NULL ;

char  *snpname = NULL ;
char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;
char  *indivname = NULL ;
char  *newindivname = NULL ;
char *badsnpname = NULL ;
char *flipsnpname = NULL ;
int flipreference = YES ;

char *poplistname = NULL ;
char *zapind1 = "Clint" ;
char *zapind2 = "Xwoman" ;
int zapmatch = YES ;

double r2thresh = -1.0 ;  
double r2genlim = 0.01 ; // Morgans 
double r2physlim = 5.0e6 ;     
double maxmissfrac = 1000.0 ; // no thresh
int killr2 = NO ;  
int mkdiploid = NO ;

int packout = -1 ;
int tersem  = YES ;
extern enum outputmodetype outputmode  ;
extern int checksizemode ;
char *omode = "packedancestrymap" ;
extern int packmode ;
int ogmode = NO ;
int fastdup = NO ;
int fastdupnum = 0  ;
int phasedmode = NO ;
int badpedignore = NO ;

int xchrom = -1 ;
int lopos = -999999999 ; 
int hipos = 999999999 ;
int minchrom = 1 ; 
int maxchrom = 97 ;

int deletedup = YES ; // only one marker at a position
char *newsnpname = NULL ;  // new map  


char  unknowngender = 'U' ;

void setomode(enum outputmodetype *outmode, char *omode)  ;
void readcommands(int argc, char **argv) ;
void outfiles(char *snpname, char *indname, char *gname, SNP **snpm, 
  Indiv **indiv, int numsnps, int numind, int packem, int ogmode) ;
void remap(SNP **s1, int nums1, SNP **s2, int nums2) ;
void remapind(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, Indiv **indm2, int numindivs, int numind2) ;
void pickx(SNP *c1, SNP *c2, SNP **px1, SNP **px2) ;
void dedupit(SNP **snpmarkers, int numsnps) ;
void flipsnps(char *fsname, SNP **snpm, int numsnps, int phasedmode) ;
void cswap(char *c1, char *c2) ;
int  mkindh2d(Indiv **indivmarkers, Indiv ***pindm2, int numindivs) ;
void remaph2d(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, Indiv **indm2, int numindivs, int numind2) ;

int main(int argc, char **argv)
{

  int **snppos ;
  int *snpindx ;
  char **snpnamelist, **indnamelist ;
  char **eglist ;
  int  lsnplist, lindlist, numeg ;
  int i,j; 
  SNP *cupt, *cupt1, *cupt2, *cupt3 ;
  Indiv *indx ;
  double gpos1,gpos2,cpos1,cpos2,gd, cd, gd100 ;
  double rthresh, zt ;
  int mpflag, ret, numvalidind, nvalid, numvalidsnps ;

  int ch1, ch2 ;
  int fmnum , lmnum ;
  int num, n1, n2 ;
  int nkill = 0 ;
  int t, k ;

  int nindiv = 0, e, f, lag=1  ;
  double xc[9], xd[4], xc2[9] ;
  double ychi, zscore, zthresh = 20.0 ;
  double y1, y2 ; 
  int nignore, numrisks = 1 ;

  char **genolist ;
  int numgenolist ;
  int maxmiss ; 
  int g1, g2, k1, k2 ;

  malexhet = YES ;    // convertf default is don't change the data
  tersem = YES ;     // no snp counts

  readcommands(argc, argv) ;

  setomode(&outputmode, omode) ;
  packmode = YES ;
  settersemode(tersem) ;

  if (r2thresh > 0.0) killr2 = YES ;
  if (badpedignore) setbadpedignore() ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks) ;


  for (i=0; i<numsnps; i++)  {  
   if (xchrom == -1) break ;  
   cupt = snpmarkers[i] ; 
   if (cupt -> chrom != xchrom) cupt -> ignore = YES ; 
   if (cupt -> ignore) continue ; 
   t = nnint(cupt -> physpos) ; 
   if ( (t< lopos) || (t >hipos)) cupt -> ignore = YES ;
  }

  nignore = 0 ;
  for (i=0; i<numsnps; i++)  {  
   cupt = snpmarkers[i] ; 
   if (cupt -> chrom > maxchrom) cupt -> ignore = YES ;  
   if (cupt -> chrom < minchrom) cupt -> ignore = YES ;  
   if (cupt -> ignore) ++nignore ;
  }

  if (numsnps == nignore) fatalx("no valid snps\n") ;


  numindivs = getindivs(indivname, &indivmarkers) ;

  if (genotypelist!= NULL) {  
    getgenos_list(genotypelist, snpmarkers, indivmarkers, 
     numsnps, numindivs, nignore) ;
  }

  else {
   setgenotypename(&genotypename, indivname) ;
   getgenos(genotypename, snpmarkers, indivmarkers, 
     numsnps, numindivs, nignore) ;
  }

  if (newsnpname != NULL) { 
    nums2 = 
     getsnps(newsnpname, &snpm2, 0.0, NULL, &nignore, numrisks) ;
     remap(snpmarkers, numsnps, snpm2, nums2) ;
     snpmarkers = snpm2 ; 
     numsnps = nums2 ;
  }

  if (newindivname != NULL) { 
    numind2 = getindivs(newindivname, &indm2) ;
    remapind(snpmarkers, numsnps, indivmarkers, indm2, numindivs, numind2) ;
    indivmarkers = indm2 ;
    numindivs = numind2 ;
  }

  if (mkdiploid) { 

    numind2 = mkindh2d(indivmarkers, &indm2, numindivs) ;
    remaph2d(snpmarkers, numsnps, indivmarkers, indm2, numindivs, numind2) ;

    indivmarkers = indm2 ;
    numindivs = numind2 ;

  }


  if (deletedup) dedupit(snpmarkers, numsnps) ; // only one marker per position
  flipsnps(flipsnpname, snpmarkers, numsnps, phasedmode) ;

  if (outputall) {
   outfiles(snpoutfilename, indoutfilename, genooutfilename, 
    snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

   printf("##end of convertfz run (outputall mode)\n") ;
   return 0 ;
  }

    k1 = indindex(indivmarkers, numindivs, zapind1) ;
    k2 = indindex(indivmarkers, numindivs, zapind2) ;
    if (k1<0) fatalx("bad zapind %s\n", zapind1) ;
    if (k2<0) fatalx("bad zapind %s\n", zapind2) ;
    for (k=0; k<numsnps; ++k) { 
     cupt = snpmarkers[k] ;
     g1 = getgtypes(cupt, k1) ;
     if (g1<0) cupt -> ignore = YES ;
     g2 = getgtypes(cupt, k2) ;
     if (g2<0) cupt -> ignore = YES ;
     if ((zapmatch) && (g1 == g2)) cupt -> ignore = YES ;
     if ((zapmatch == NO) && (g1 != g2)) cupt -> ignore = YES  ;
    }


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
  else 
  setstatus(indivmarkers, numindivs, "Case") ;

  numsnps = rmsnps(snpmarkers, numsnps) ;
  numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;

  if (killr2) {
   nkill = killhir2(snpmarkers, numsnps, numindivs, r2physlim, r2genlim, r2thresh) ;
   if (nkill>0) printf("killhir2.  number of snps killed: %d\n", nkill) ;
  }

  numvalidind = 0 ;
  for (i=0; i<numindivs; ++i)  { 
   indx = indivmarkers[i] ;
   if (indx -> ignore) continue ; 
   if (numvalidgtind(snpmarkers, numsnps, i) ==0) { 
    indx -> ignore = YES ; 
    printf("no data for individual: %s\n", indx -> ID) ;
   }
   if (indx -> ignore == NO) ++numvalidind ;
  }

  maxmiss  = (int) (maxmissfrac * (double) numvalidind) ;
  printf("numvalidind:  %5d  maxmiss: %5d\n", numvalidind, maxmiss)  ;
  if (numvalidind  == 0) fatalx("no valid samples!\n") ;

  for (k=0; k<numsnps; ++k) {  
   if (maxmiss>numvalidind) break ;
   cupt = snpmarkers[k] ;
   t = numvalidind - numvalidgtypes(cupt) ;
// printf("zz %20s %4d %4d\n", cupt -> ID, t, numvalidind-t) ;
   if (maxmiss <= t) { 
    cupt -> ignore = YES ;
   }
/**
   if (numvalidind ==  t) { 
    printf("no data for snp: %s\n", cupt -> ID) ;
    cupt -> ignore = YES ;
   }
*/

  }

  if (fastdup)  {  
   if (fastdupnum > 0) setfastdupnum(fastdupnum) ;
   fastdupcheck(snpmarkers, indivmarkers, numsnps, numindivs) ;  
  }

  if (decim>0) {  
   snpdecimate(snpmarkers, numsnps, decim, dmindis, dmaxdis) ;
  }

  outfiles(snpoutfilename, indoutfilename, genooutfilename, 
   snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode) ;

  printf("##end of convertfz run\n") ;
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

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
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
   getstring(ph, "genotypelist:", &genotypelist) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "flipsnpname:", &flipsnpname) ;
   getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "indivoutname:", &indoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "snpoutname:", &snpoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "genooutfilename:", &genooutfilename) ; 
   getstring(ph, "genotypeoutname:", &genooutfilename) ; /* changed 11/02/06 */
   getstring(ph, "outputformat:", &omode) ;
   getstring(ph, "outputmode:", &omode) ;
   getint(ph, "checksizemode:", &checksizemode) ;
   getint(ph, "badpedignore:", &badpedignore) ;

   getint(ph, "outputgroup:", &ogmode) ;
   getint(ph, "malexhet:", &malexhet) ;
   getint(ph, "nomalexhet:", &malexhet) ; /* changed 11/02/06 */
   getint(ph, "tersemode:", &tersem) ;
   getint(ph, "familynames:", &familynames) ;
   getint(ph, "packout:", &packout) ; /* now obsolete 11/02/06 */
   getint(ph, "decimate:", &decim) ; 
   getint(ph, "dmindis:", &dmindis) ; 
   getint(ph, "dmaxdis:", &dmaxdis) ; 
   getint(ph, "fastdup:", &fastdup) ;
   getint(ph, "flipreference:", &flipreference) ;
   getint(ph, "fastdupnum:", &fastdupnum) ;
   getint(ph, "killr2:",  &killr2) ;
   getint(ph, "hashcheck:", &hashcheck) ;
   getint(ph, "outputall:", &outputall) ;
   getint(ph, "sevencolumnped:", &sevencolumnped) ;
   getint(ph, "phasedmode:", &phasedmode) ;

   getdbl(ph, "r2thresh:", &r2thresh) ;
   getdbl(ph, "r2genlim:", &r2genlim) ;
   getdbl(ph, "r2physlim:", &r2physlim) ;

   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "lopos:", &lopos) ;
   getint(ph, "hipos:", &hipos) ;

   getint(ph, "minchrom:", &minchrom) ;
   getint(ph, "maxchrom:", &maxchrom) ;
   getdbl(ph, "maxmissfrac:", &maxmissfrac) ;
   
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "newsnpname:", &newsnpname) ;
   getstring(ph, "newindivname:", &newindivname) ;
   getint(ph, "deletedup:", &deletedup) ;
   getint(ph, "mkdiploid:", &mkdiploid) ;
   getint(ph, "zapmatch:", &zapmatch) ;

   writepars(ph) ;
   closepars(ph) ;

}
void remap(SNP **s1, int nums1, SNP **s2, int nums2) 
{
  SNP *cupt1, *cupt2 ; 
  SNP tcupt, *cupt ;
  int i, k ;
  
  for (i=0; i<nums2; ++i)  { 
   cupt2 = s2[i] ; 
   k = snpindex(s1, nums1, cupt2 -> ID) ; 
   if (k<0) { 
    printf("%20s not found\n") ; 
    cupt2 -> ignore = YES ;
    continue ;
   }
   cupt1 = s1[k] ; 
   tcupt = *cupt2 ;
   *cupt2 = *cupt1 ;
   cupt = &tcupt ;
   cupt2 -> chrom =   cupt -> chrom ; 
   cupt2 -> genpos =  cupt -> genpos ; 
   cupt2 -> physpos = cupt -> physpos ; 
  }
  freesnpindex() ;
}


void 
remapind(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, Indiv **indm2, int numindivs, int numind2) 

{

  int *g1, *g2 ;
  int *tind, t, i, j, k ; 
  Indiv *indx ;
  SNP *cupt ;

  ZALLOC(tind, numind2, int) ;
  ZALLOC(g2, numind2, int) ;
  ZALLOC(g1, numindivs, int) ;

  for (k=0; k<numind2; ++k) { 
   indx = indm2[k] ;   
   t = tind[k] = indindex(indivmarkers, numindivs, indx -> ID) ;
   if (t<0) fatalx("bad newindiv: %s\n", indx -> ID) ;
  }

  for (i=0; i<numsnps; i++)  { 
   cupt = snpmarkers[i]  ;

   for (j=0; j<numindivs; ++j) { 
    g1[j] = getgtypes(cupt, j) ;
   }

   for (k=0; k< numind2; ++k)  { 
     g2[k] = g1[tind[k]]  ;
   }

   ivclear(g1, -1, numindivs) ;
   copyiarr(g2, g1, numind2) ;

   for (k=0; k<numind2; ++k) { 
    putgtypes(cupt, k, g1[k]) ;
   }

  }

  free(g1) ; 
  free(g2) ;
  free(tind) ;

}



void  dedupit(SNP **snpmarkers, int numsnps) 
{
  SNP *cupt1, *cupt2 ;   
  SNP *x1, *x2 ;   
  int k ;

  cupt1 = NULL ;

  for (k=0; k<numsnps; ++k) { 
   cupt2 = snpmarkers[k] ;
   if (cupt2 -> ignore) continue ;
   if (cupt1 == NULL) { 
    cupt1 = cupt2 ;
    continue   ;
   }
   if (cupt1 -> chrom != cupt2 -> chrom) { 
    cupt1 = cupt2 ;
    continue   ;
   }
   if (cupt1 -> physpos != cupt2 -> physpos) { 
    cupt1 = cupt2 ;
    continue   ;
   }
   pickx(cupt1, cupt2, &x1, &x2) ;  // x2 bad
   x2 -> ignore = YES ;
   cupt1 = x1 ; 
  }
}
void pickx(SNP *c1, SNP *c2, SNP **px1, SNP **px2) 
{
// *px1 is retained *px2 dropped; Try and keep shorter rsnames (strip AFFX for instance)

  char *ch1, *ch2 ;
  int l1, l2, t ;
   
   ch1 = c1 -> ID ;
   ch2 = c2 -> ID ;
   l1 = strlen(ch1) ;
   l2 = strlen(ch2) ;

   if (l1<l2) { 
    *px2 = c2 ; 
    *px1 = c1 ; 
    return ;
   }
   if (l2<l1) { 
    *px2 = c1 ; 
    *px1 = c2 ; 
    return ;
   }
   t = strcmp(ch1, ch2) ;
   if (t>=0) { 
    *px2 = c2 ; 
    *px1 = c1 ; 
    return ;
   }
   *px2 = c1 ; 
   *px1 = c2 ; 
   return ;
}
void flipsnps(char *fsname, SNP **snpm, int numsnps, int phasedmode) 
{
  FILE *fff ;
  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *ss ;
  int nsplit, n, k ;
  SNP *cupt ; 

  if (fsname == NULL) return ;
  openit (fsname, &fff, "r") ;

  freesnpindex() ;

  while (fgets(line, MAXSTR, fff) != NULL)  {
    nsplit = splitup(line, spt, MAXFF) ; 
    if (nsplit==0) continue ;
    if (spt[0][0] == '#') { 
     freeup(spt, nsplit) ;
     continue ;
    }
    
     k = snpindex(snpm, numsnps, spt[0]) ;      
     if (k>=0) {
      cupt = snpm[k] ;
      if (phasedmode == NO)  flipalleles(cupt) ;
      if (phasedmode == YES)  flipalleles_phased(cupt) ;
// just flips genotypes  
      if (flipreference) cswap(&cupt -> alleles[0], &cupt -> alleles[1]) ;
     }
    freeup(spt, nsplit) ;
  }
  fclose (fff) ;
}
void cswap(char *c1, char *c2) 
{
  char cc ; 
  
  cc = *c1  ; 
  *c1 = *c2 ;  
  *c2 = cc  ; 

}
