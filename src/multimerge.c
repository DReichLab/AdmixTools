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
#include <xsearch.h>

#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "egsubs.h"

int maxsnpnum = 1000 ; 
int maxindnum = 500 ; 


#define WVERSION   "400"


char *target = "rs3094315" ;
int debug ; 

#define MAXFL  50
#define MAXSTR  512

extern int packmode;

char *trashdir = "/var/tmp";
extern int verbose;
int qtmode = NO;
extern enum outputmodetype outputmode;

Indiv **indivmarkers, **ind1;
SNP **snpmarkers, **snp1 = NULL, **snp1x = NULL;
int numsnps, numindivs;
int singlesnpname = NO ; 
char *dirname = NULL ;

char *globalsnpname = NULL ;
char *genotypename = NULL;
char *snpname = NULL;
char *genooutfilename = NULL;
char *indoutfilename = NULL;
char *snpoutfilename = NULL;
char *indivname = NULL;
char *badsnpname = NULL;
char *goodsnpname = NULL;
char *outstem = NULL ; 
char *badpairsname = NULL;
char *markername = NULL;
char *idname = NULL;
char *mergelistname = NULL ;

char **biglist ; 

char *outputname = NULL;
FILE *ofile;

char *omode = "packedancestrymap";
extern int packmode;
int ogmode = NO;

double fakespacing = 0.0;

char unknowngender = 'U';

void readcommands (int argc, char **argv);
void dophyscheck (SNP ** snpm, int numsnps);
int loadxlist(char *fname, char **names, int *plen, int xtype) ;  
void copysnpx(SNP **snpa, SNP **snpb, int numsnp) ;

Indiv **bigindivs = NULL ; 
SNP **bigsnps = NULL ; 

char **indnames,  **snpnames ; 

void setbig(char **inames, char **snames, int ilen, int slen) ;

 typedef struct { 
   int index ; 
   int xtype ; 
 } XITEM ; 

int
main (int argc, char **argv)
{

  int i, j, k, g, x, t;
  SNP *cupt, *cupt2, *cuptt = NULL;
  Indiv *indx, *indx2 ;
  char **inames, **snames, **stemnames, **gnames ; 
  int numstem ; 
  char ss[MAXSTR], ss2[MAXSTR], *sx ; 
  int rlen ;
  double y ;
  int nums1, numi1, nignore, nums1x = -1 ;
  int numrisks = 1 ; 
  int aind, asnp ;
  int a, b, xa, xb, g1, g2 ;

  unsigned char * packgg, *pp, *px ;
  int numout, packout = -1 ;
  double ymem ; 

  readcommands(argc, argv) ;

  printf ("## multimerge version: %s\n", WVERSION);

  cputime(0) ;
  calcmem(0) ;

  if (mergelistname == NULL) fatalx("no mergelistname!\n") ;
  numstem = numlines(mergelistname) ; 
  ZALLOC(stemnames, numstem, char *) ; 
  numstem = getss(stemnames, mergelistname) ;
  ZALLOC(inames, numstem, char *) ;
  ZALLOC(snames, numstem, char *) ;
  ZALLOC(gnames, numstem, char *) ;
//  printstrings(stemnames, numstem) ;

  maxsnpnum = MAX(maxsnpnum, 20*1000*1000) ; 
  maxindnum = MAX(maxindnum, 200*1000) ; 

  setdupcheck(NO) ;
  for (i=0; i<numstem; ++i) { 
   sx = stemnames[i] ; 

   makedfn(dirname, sx, ss2, MAXSTR) ;
   strcpy(ss, ss2) ; strcat(ss, ".ind") ; inames[i] = strdup(ss) ; 
   strcpy(ss, ss2) ; strcat(ss, ".snp") ; snames[i] = strdup(ss) ; 
   strcpy(ss, ss2) ; strcat(ss, ".geno") ; gnames[i] = strdup(ss) ; 


  }
  if (outstem == NULL) fatalx("*** no outputfiles!\n") ;
   sx = outstem ; 
   makedfn(dirname, sx, ss2, MAXSTR) ;
    strcpy(ss, ss2) ; strcat(ss, ".ind") ; indoutfilename = strdup(ss) ; 
    strcpy(ss, ss2) ; strcat(ss, ".snp") ; snpoutfilename   = strdup(ss) ; 
    strcpy(ss, ss2) ; strcat(ss, ".geno") ; genooutfilename = strdup(ss) ; 


/**
  printnl() ;
  printstrings(inames, numstem) ; 
  printnl() ;
  printstrings(snames, numstem) ; 
*/
 
  ZALLOC(biglist, maxsnpnum  + maxindnum, char *) ;  

  setbig(inames, snames, numstem, numstem) ;  

  printf("numsnps, numindivs: %d %d\n", numsnps, numindivs) ;

  ZALLOC (bigindivs, numindivs, Indiv *);

  for (i = 0; i < numindivs; i++) {
    ZALLOC (bigindivs[i], 1, Indiv);
  }
  clearind (bigindivs, numindivs);
  for (i = 0; i < numindivs; i++) {
   indx = bigindivs[i] ; 
   strcpy(indx -> ID, indnames[i]) ;
   indx -> affstatus = -99; // mark as unset
  }

  ZALLOC (bigsnps, numsnps, SNP *);
  y = (double) (numindivs * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));

  for (i = 0; i < numsnps; i++) {
    ZALLOC (bigsnps[i], 1, SNP);
    clearsnp (bigsnps[i]) ;
  }
  
  setomode (&outputmode, omode);
  packmode = YES;


//  if (singlesnpname) putsnpordered(YES) ;
  
/**
    if (globalsnpname != NULL) { 
     nums1 = getsnps (globalsnpname, &snp1, 0.0, NULL, &nignore, numrisks);
     for (a = 0; a< nums1; ++a) { 
      cupt2 = snp1[a] ;
      sx  = cupt2 -> ID ; 
      xa = xfindit(sx) ;
      cupt = bigsnps[xa] ;
      if (cupt -> ngtypes == 0) { // SNP unset
       *cupt = *cupt2 ; 
       cupt -> ngtypes = numindivs ;
       ZALLOC(cupt -> pbuff, rlen, unsigned char) ; 
       cclear(cupt -> pbuff, 0XFF, rlen) ; 
      } 
     }
    
     freesnps(&snp1, nums1) ; 
   }
*/

  for (j=0; j<numstem; ++j) { 
   snpname = snames[j] ; 
   indivname = inames[j] ;
   genotypename = gnames[j] ;

   if (snp1x == NULL)  {
    nums1 = getsnps (snpname, &snp1, 0.0, NULL, &nignore, numrisks);
   }

   if (snp1x != NULL) { 
    for (k=0; k<nums1x; ++k) { 
     snp1[k] = snp1x[k] ; // just pointer
    }
    nums1 = nums1x ;
   }


   if ((snp1x == NULL) && (singlesnpname == YES))  {   
    ZALLOC(snp1x, nums1,  SNP *) ;
    copysnpx(snp1, snp1x, nums1) ;
    nums1x = nums1 ;
   }
 


   numi1 = getindivs (indivname, &ind1);
   getgenos (genotypename, snp1, ind1, nums1, numi1, nignore);

   for (k = 0; k< nums1; ++k) { 
    sx = snp1[k] -> ID ;  asnp = xfindit(sx) ;
 // printf("snp: %s %d\n", sx, asnp) ;
   }
   for (k = 0; k< numi1; ++k) { 
    sx = ind1[k] -> ID ;  aind = xfindit(sx) - maxsnpnum;
   break ;
//  printf("ind: %s %d\n", sx, aind) ;
   }

   for (a = 0; a< nums1; ++a) { 
    debug = NO ;
    cupt2 = snp1[a] ;
    sx  = cupt2 -> ID ; 
    t = strcmp(sx, target) ; 
    t = -1 ;
    if (t==0) debug = YES ;
    xa = xfindit(sx) ;
    if (xa<0) fatalx("bad search %s\n", sx) ;
    cupt = bigsnps[xa] ;
    if (debug) cuptt = cupt ;
    if (cupt -> ngtypes == 0) { // SNP unset
     *cupt = *cupt2 ;
 //  printf ("allocating %s\n", cupt ->  ID) ;
     cupt -> ngtypes = numindivs ;
     ZALLOC(cupt -> pbuff, rlen, char) ; 
     cclear((unsigned char *) cupt -> pbuff, 0XFF, rlen) ; 
    } 
    
   for (b = 0; b< numi1; ++b) { 
    indx2 = ind1[b] ;  
    sx  = indx2 -> ID ; 
    xb = xfindit(sx) - maxsnpnum ;
    if (xb<0) fatalx("bad search %s\n", sx) ;
    indx = bigindivs[xb] ;
//  printf("zz %d %d %s %s %s\n", b, xb, indx -> ID, indx2 -> ID, indx2 -> egroup) ;
    if (indx -> affstatus == -99) { 
      *indx = *indx2 ; 
    }
     g1 = getgtypes(cupt, xb) ; 
     g2 = getgtypes(cupt2, b) ;
     putgtypes(cupt2, b, -1) ;
    if (g1>=0) continue ;
     putgtypes(cupt, xb, g2) ;
    }
   }

    freepackgenos() ;
     
    for (a = 0; a< nums1; ++a) { 
     cupt = snp1[a] ; 
     cupt -> ngtypes = 0 ; 
     if (cupt -> gtypes != NULL) { 
      free(cupt -> gtypes) ; 
      cupt -> gtypes = NULL ;
     }
    }
   
   
   if (singlesnpname == NO)  freesnps(&snp1, nums1) ; 
   freeinds(&ind1, numi1) ; 

   printf("stem %d processed -- topheap: %p\n", j, topheap()) ;

  if (cuptt != NULL) { 
   for (b=0; b<numindivs; ++b) {
    g1 = getgtypes(cuptt, b) ; 
//  printf("yy2 %20s %d %d\n", cuptt -> ID, b, g1) ;
   }
  }

 }

  settersemode(YES) ; 

  numout = outfiles (snpoutfilename, indoutfilename, genooutfilename,
     bigsnps, bigindivs, numsnps, numindivs, packout, ogmode);


  ymem = calcmem(1)/1.0e6 ; 
  printf("##end of multimerge: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;
}

void setbig(char **inames, char **snames, int ilen, int slen) 
{
  char *fname ; 
  FILE *fff ;
  int k, xilen = 0, xslen = 0 ; 

// xtype 1 for Indiv 2 for SNP

  int tlen, thlen ; 

  tlen = maxindnum + maxsnpnum ; 
  thlen = 2*tlen ;
  xhinit(tlen) ; 
  ZALLOC(indnames, maxindnum, char *) ; 
  ZALLOC(snpnames, maxsnpnum, char *) ; 
  
  ZALLOC(biglist, tlen, char *) ;

  xslen = xilen = 0 ;
  for (k=0; k<ilen; ++k) { 
   fname = inames[k] ; 
   loadxlist(fname, indnames, &xilen, 1) ; 
// printnl() ; 
  }
  numindivs = xilen ;
  if (globalsnpname != NULL) { 
     loadxlist(globalsnpname, snpnames, &xslen, 2) ; 
  }
 
  for (k=0; k<slen; ++k) { 
   fname = snames[k] ; 
    if ((k==0) || (singlesnpname == NO)) {
     loadxlist(fname, snpnames, &xslen, 2) ; 
    }

   }

//  printstrings(snpnames, xslen) ; 
  numsnps = xslen ;

 return  ;



}
int loadxlist(char *fname, char **names, int *plen, int xtype) 
{ 
   int blen, len, nraw, n, k, t, x  ; 
   XITEM xitem, *xpt ;
   char **lnames, *sx ; 
   ENTRY zentry, *zpt ; 

   nraw  = numlines(fname) ; 
   ZALLOC(lnames, nraw, char *) ;
   n = getss(lnames, fname) ;

   zpt = &zentry ;
   xpt = &xitem ; 
   xpt -> xtype = xtype ;
   blen = len = *plen ; 
   if (xtype == 1) blen += maxsnpnum ;  

   for (k=0; k<n; ++k) { 
    sx = lnames[k] ;  
    xpt -> index = len ;      
    zentry.key = sx ;  
    zentry.data = xpt ; 
    if (xlookup(sx, FIND) < 0) { ;  
     xstore(sx, blen) ;  
     names[len] = strdup(sx) ;
     biglist[blen] = strdup(sx) ;
     ++len ;
     ++blen ; 
    }
//  printf("zzz\n") ; 
//  dumpxs() ;
   }

   *plen = len ;
   freeup(lnames, nraw) ;
   return len ; 

}

void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  char *parname = NULL;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n;

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i) {

    case 'p':
      parname = strdup (optarg);
      break;

    case 'v':
      printf ("version: %s\n", WVERSION);
      break;

    case 'V':
      verbose = YES;
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }


  pcheck (parname, 'p');
  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);
  dostrsub (ph);

/**
DIR2:  /fg/nfiles/admixdata/ms2
SSSS:  DIR2/outfiles 
genotypename: DIR2/autos_ccshad_fakes
eglistname:    DIR2/eurlist  
output:        eurout
*/
  getint (ph, "packmode:", &packmode);  // controls internals 

  getstring (ph, "mergelist:", &mergelistname);
  getstring (ph, "mergelistname:", &mergelistname);
  getstring (ph, "outstem:", &outstem);
  getstring (ph, "globalsnpname:", &globalsnpname);
  getint (ph, "singlesnpname:", &singlesnpname);
  getint (ph, "maxsnpnum:", &maxsnpnum);
  getint (ph, "maxindnum:", &maxindnum);
  getstring (ph, "outputformat:", &omode);
  getstring (ph, "outputmode:", &omode);
  getstring (ph, "dirname:", &dirname);
  writepars (ph);
  closepars (ph);

}

void
dophyscheck (SNP ** snpm, int numsnps)
{
// catch places where physpos genpos are in opposite order
  SNP *cupt, *cuptold;
  int i;

  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (i == 0)
      cuptold = cupt;
    if (cupt->isfake)
      continue;
    if (cupt->ignore)
      continue;
    if (cupt->chrom == cuptold->chrom) {
      if (cupt->physpos < cuptold->physpos) {
        printf ("physcheck %20s %15s %12.3f %12.3f %13.0f %13.0f\n",
                cuptold->ID, cupt->ID,
                cuptold->genpos, cupt->genpos,
                cuptold->physpos, cupt->physpos);
      }
    }
    cuptold = cupt;
  }
}
void copysnpx(SNP **snpa, SNP **snpb, int numsnp) 
{
    int k ;  
    SNP *cupt ;

    for (k=0; k<numsnp; ++k) { 
     cupt = snpb[k] ; 
     if (cupt != NULL) freecupt(&cupt) ; 
     ZALLOC(snpb[k], 1, SNP) ;  
     cupt = snpb[k] ; 
     clearsnp(cupt) ; 
     *cupt = *snpa[k] ;  
     cupt -> pbuff = NULL ; 
     cupt -> ngtypes = 0 ; 
    }
}
     



