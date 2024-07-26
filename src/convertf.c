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
#include "egsubs.h"
#include "exclude.h"
#include "h2d.h"

#define WVERSION   "8600"
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

 Various bugfixes and checks for file size  
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
 flipstrandname moves alleles to opposite strand
 maxmissing added (counts alleles like smartpca)
 minor bug fixed for plink homozygous files
 polarize added (force homozygotes to 2 if possible)
   if not possible SNP set -> ignore

 bigread now set 
 map files output 23, 24 for X, Y
 fastdupthresh, fastdupkill added
 remap now allows allele changes when genotypes are killed

 zerodistance added 

 mcio now calls toupper on alleles
 bugfix for newsnpname + phasedmode

 downsample added (make pseudo homozygotes)
 chimpmode added support for chr2a etc

 support for snps out of order in packed format
 (pordercheck: NO) 
 remapcheck NO (useful for moving reference sequence)

 mcio.c imported from eig 5.0.1
 randommode, seed added (randomizes fastdup)
 newignore: NO => newsnpname with new snps will have genotypes filled in with unknown (default is don't output) 
 inddupcheck added (compulsory) 
 usesamples added => poplistname = NULL)
 copyalleles (if newsnpname alleles are copied from old file) 
 familypopnames (use egroup for family in .fam output
 dupcheck contains iter number

 minvalpop :: every pop must have at least this number of valids (default not set) 
 fillmissing: added 
 better handling of seed

 hiresgendis added

 O2 version
 support for .prob files

 wipeoutchrom added 
 probprintname added 

 imputeusingpops added  
 mkhaploid added 

 forcemiss added 
 slowdup added  

 .map output files no longer print alleles (use .pedsnp if wanted)  
 transpose packed format (TGEN)
 memorymapping implemented

 transform hets (1 -> 2) 
 shuffle 
 randomhetfix
 indivlistname added (overrules poplistname) 
 instem, outstem added
*/


#define MAXFL  50
#define MAXSTR  512

char *trashdir = "/var/tmp";
int qtmode = NO;
Indiv **indivmarkers, **probindivs, **indm2 ;
int numsnps, numindivs, numprobindivs, numind2 ;
SNP **snpmarkers;
SNP **snpm2;
int zerodistance = NO;		// YES => force gdis 0
int downsample = NO;		// make pseudo homozygotes
int pordercheck = YES;
int familypopnames = NO;
int hiresgendis = NO ;
int memorymap = NO ; 
int transformhets = NO ; 
int randomhetfix = NO ; 
int shuffle = NO ;

int nums2;

char *genotypename = NULL;
char *genotypelist = NULL;

char *snpname = NULL;
char *indoutfilename = NULL;
char *snpoutfilename = NULL;
char  *genooutfilename = NULL;
char *tgenooutfilename = NULL;
char *indivname = NULL;
char *newindivname = NULL;

char *instem = NULL, *outstem = NULL ;

char *probfilename = NULL;
char *proboutfilename = NULL;
char *probindivname = NULL;
char *probprintname = NULL ; 
char *probprintid = NULL ; 



char *badsnpname = NULL;
char *xregionname = NULL;
char *deletesnpoutname = NULL;
char *flipsnpname = NULL;
char *flipstrandname = NULL;
int flipreference = YES;
int remapcheck = YES;

char *poplistname = NULL;
char *indivlistname = NULL ; 

char *fillmissingpoplistname = NULL ; 
int haploidfill = NO ;
char **fpops ;
int nfpops = 0 ;
double makemiss = -1 ;

double r2thresh = -1.0;
double r2genlim = 0.01;		// Morgans 
double r2physlim = 5.0e6;
double maxmissfrac = 1000.0;	// no thresh
int maxmiss = -1;		// no thresh
int minvalpop = -1 ;
int killr2 = NO;
int mkdiploid = NO;
int mkhaploid = NO;

int packout = -1;
int tersem = YES;
int randommode = NO;
int seed = 0;

extern enum outputmodetype outputmode;
extern int checksizemode;
char *omode = "packedancestrymap";
extern int packmode;
int ogmode = NO;
int fastdup = NO;
int slowdup = NO;
int fastdupnum = 10;
double fastdupthresh = .75;
double fastdupkill = .75;
char *polarid = NULL;
int polarindex = -1;

int phasedmode = NO;
int badpedignore = NO;
int chimpmode = NO;

int xchrom = -1;
int lopos = -999999999;
int hipos = 999999999;
int minchrom = 1;
int maxchrom = 97;
int wipeoutchrom = -1 ;

int deletedup = YES;		// only one marker at a position
char *newsnpname = NULL;	// new map  
int newignore = YES;		// default ignore snps not in old list
int polarcheck = NO;
int copyalleles = NO;
int rmcompress = YES ; 

char *usesamples = NULL;

char unknowngender = 'U';
double nhwfilter = -1;


void setomode (enum outputmodetype *outmode, char *omode);
void readcommands (int argc, char **argv);
void remap (SNP ** s1, int nums1, SNP ** s2, int nums2);
void remapind (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	       Indiv ** indm2, int numindivs, int numind2);
void pickx (SNP * c1, SNP * c2, SNP ** px1, SNP ** px2);
void dedupit (SNP ** snpmarkers, int numsnps);
void flipsnps (char *fsname, SNP ** snpm, int numsnps, int phasedmode);
void flipstrand (char *fsname, SNP ** snpm, int numsnps);
int mkindh2d (Indiv ** indivmarkers, Indiv *** pindm2, int numindivs);
void remaph2d (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	       Indiv ** indm2, int numindivs, int numind2);
void flip1 (SNP * cupt, int phasedmode, int flipreference);

void fixaa (SNP * cupt1, SNP * cupt2);
void fvalg (SNP * cupt, int val);
char cxx (char *c1, char *c2);
void downsamp (SNP * cupt);
void forcemiss (double yprob);
int setsamp (Indiv ** indivmarkers, int numindivs, char *usesamples);
int testmisspop(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, int numindivs, int minvalpops)   ;
int fillmiss(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char **fpops, int nfpops)  ;
int fixsnpdistance(SNP **snpm, int numsnps)  ; 
long loadprobpack(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char *bigbuff) ;
void doshuffle(SNP **snpm, int numsnps, int numindivs)   ; 



int
main (int argc, char **argv)
{

  int **snppos;
  int *snpindx;
  char **snpnames, **indnames;
  char **eglist;
  int lsnplist, lindlist, numeg;
  int i, j;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;
  double gpos1, gpos2, cpos1, cpos2, gd, cd, gd100;
  double rthresh, zt;
  int mpflag, ret, numvalidind, nvalid, numvalidsnps;

  int ch1, ch2;
  int fmnum, lmnum;
  int num, n1, n2;
  int nkill = 0;
  int t, k, g;
  int k1, k2 ;

  int nindiv = 0, e, f, lag = 1;
  double xc[9], xd[4], xc2[9];
  double ychi, zscore, zthresh = 20.0;
  double y1, y2, ymem;
  int nignore, numrisks = 1;

  char **genolist;
  int numgenolist;
  char c1, c2;
  int t1, t2, x;

  unsigned char *packp, *packp2 ; 
  long plen, plen2, numx ;  
  int rl2 = 4 ; 

  int numout = -1 ;

  malexhet = YES;		// convertf default is don't change the data
  tersem = YES;			// no snp counts

  int jlast = -1 ; 
  FILE *probfile = NULL ;
  double *pp ; 

  readcommands (argc, argv);

  cputime(0) ;
  calcmem(0) ;
  
  printf("## %s version: %s\n", argv[0], WVERSION) ;

  if (seed == 0)  {
   seed = seednum() ; 
  }

  SRAND(seed) ; 

  if (fastdup) { 
    randommode = YES;
  }

  if (randommode)
    printf("seed: %d\n", seed) ;

  if (chimpmode) {
    setchimpmode (YES);
    setchr (YES);
  }
  if (familypopnames) {
    setfamilypopnames (YES);
  }
  if (outputall) memorymap = NO ;  // need to read in all data 
  if (memorymap) setmemorymap(YES) ;

  if (instem != NULL) setinfiles(&indivname, &snpname, &genotypename, instem) ;
  if (outstem != NULL) setoutfiles(&indoutfilename, &snpoutfilename, &genooutfilename, outstem) ;
  if (indivlistname != NULL) poplistname = NULL ;

  if (rmcompress == NO) { 
   fatalx("rmcompress obsolete -- use outputall!\n") ; 
  }

  setomode (&outputmode, omode);
  packmode = YES;
  settersemode (tersem);

  if (hiresgendis) sethiressnp() ;

  if (r2thresh > 0.0)
    killr2 = YES;
  if (badpedignore)
    setbadpedignore ();

  setpordercheck (pordercheck);

  numsnps =
    getsnps (snpname, &snpmarkers, 0.0, badsnpname, &nignore, numrisks);

  for (i = 0; i < numsnps; i++) {
    if (xchrom == -1)
      break;
    cupt = snpmarkers[i];
    if (cupt->chrom != xchrom)
      cupt->ignore = YES;
    if (cupt->ignore)
      continue;
    t = nnint (cupt->physpos);
    if ((t < lopos) || (t > hipos))
      cupt->ignore = YES;
  }

  nignore = 0;
  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    if (cupt->chrom > maxchrom)
      cupt->ignore = YES;
    if (cupt->chrom < minchrom)
      cupt->ignore = YES;
    if (cupt->ignore)
      ++nignore;
  }

  if (numsnps == nignore)
    fatalx ("no valid snps\n");


  numindivs = getindivs (indivname, &indivmarkers);
  if (polarid != NULL) {
    polarindex = indindex (indivmarkers, numindivs, polarid);
    if (polarindex < 0)
      fatalx ("polarid %s not found\n", polarid);
  }

  if (poplistname != NULL) {
    ZALLOC (eglist, numindivs, char *);
    numeg = loadlist (eglist, poplistname);
    seteglist (indivmarkers, numindivs, poplistname);
    for (i = 0; i < numindivs; ++i) {
      indx = indivmarkers[i];
      if (indx->affstatus == NO)
	indx->ignore = YES;
    }
  }
  else {
    setstatus (indivmarkers, numindivs, NULL) ; 
  }

  
  setgk(indivmarkers, numindivs, poplistname, NULL, NULL) ; 

  inddupcheck (indivmarkers, numindivs);

  if (indivlistname != NULL) { 
   lindlist = numlines(indivlistname) ; 
   ZALLOC(indnames, lindlist, char *) ; 
   lindlist = getlist(indivlistname, indnames) ;
   for (k=0; k<numindivs; ++k) { 
    indx = indivmarkers[k] ; 
    indx -> ignore = YES ; 
    indx -> gkode = -1 ; 
    t = indxstring(indnames, lindlist, indx -> ID) ; 
    if (t>=0) { 
     indx -> ignore = NO ; 
     indx -> gkode = 1 ; 
    }
   }
  }

  if (genotypelist != NULL) {
    getgenos_list (genotypelist, snpmarkers, indivmarkers,
		   numsnps, numindivs, nignore);
  }

  else {
    setgenotypename (&genotypename, indivname);
    getgenos (genotypename, snpmarkers, indivmarkers,
	      numsnps, numindivs, nignore);
  }

  for (i=0; i< numsnps; ++i) { 
   if (transformhets == NO) break ;
// het -> countallele
   cupt = snpmarkers[i] ; 
   for (j=0; j<numindivs; ++j) { 
    g = getgtypes(cupt, j) ; 
    if (g==1) putgtypes(cupt, j, 2) ;  
   }
  }

  for (i=0; i< numsnps; ++i) { 
   if (randomhetfix == NO) break ;
// het -> random homozygote
   cupt = snpmarkers[i] ; 
   for (j=0; j<numindivs; ++j) { 
    g = getgtypes(cupt, j) ; 
    if (g==1) {
      g = 2 * ranmod(2) ; 
      putgtypes(cupt, j, g) ;  
    }
   }
  }

// read in probs
  numprobindivs = 0 ; 
  plen = 0 ; 
  packp = NULL ;

  if ((probfilename != NULL) && (probindivname == NULL)) probindivname = indivname ; 

  if (probindivname != NULL)  { 
   numprobindivs = getindivs (probindivname, &probindivs);
    for (j=0; j<numprobindivs; ++j) { 
      indx = probindivs[j] ; 
      k = indindex(indivmarkers, numindivs, indx -> ID) ;
      indx -> idnum = k ;
      if (k<0) {
       indx -> ignore = YES ; 
       printf("*** warning *** ID %s missing in indivs\n", indx -> ID) ;
     }
    } 
    plen = numsnps*numprobindivs*rl2 ;  
    ZALLOC(packp, plen, unsigned char) ;
    printf("calling inprobx: hashcheck: %d\n",  hashcheck) ; 
    inprobx (probfilename,  snpmarkers, probindivs, numsnps, numprobindivs, (char *) packp) ;
  }


   
 if (numprobindivs > 0) {
  for (i=0; i<numsnps; i++) {
   cupt = snpmarkers[i] ;
   cupt -> probbuff = (char *) packp + i*numprobindivs*rl2 ;
   ZALLOC(cupt -> diplike, numindivs, double *) ; 
   for (j=0; j<numindivs; ++j) { 
     indx = indivmarkers[j] ; 
     if (indx -> ignore) continue ;  
     ZALLOC(cupt -> diplike[j], 3, double) ;
   } 
  }
  for (i=0; i<numsnps; i++) {
   cupt = snpmarkers[i] ;
   for (j=0; j<numprobindivs; ++j) { 
     indx = probindivs[j] ;  
     if (indx -> ignore) continue ;
     k = indx -> idnum ;
     if (k<0) continue ;  
     indx = indivmarkers[k] ; 
     if (indx -> ignore) continue ; 
     x = loaddiplike(cupt -> diplike[k], (unsigned char *) cupt -> probbuff + j*rl2) ; 
     jlast = MAX(j, jlast) ;
   } 
  }
 }
 if (makemiss > 0) printf("forcing missing with prob: %9.f\n", makemiss) ;

  if (newsnpname != NULL) {
    numindivs = rmindivs (snpmarkers, numsnps, indivmarkers, numindivs);
// clean up before funky stuff
    clearsnpord ();
    nums2 = getsnps (newsnpname, &snpm2, 0.0, NULL, &nignore, numrisks);
    remap (snpmarkers, numsnps, snpm2, nums2);
    snpmarkers = snpm2;
    numsnps = nums2;
  }

  if (newindivname != NULL) {
    numind2 = getindivs (newindivname, &indm2);
    remapind (snpmarkers, numsnps, indivmarkers, indm2, numindivs, numind2);
    indivmarkers = indm2;
    numindivs = numind2;
    if (polarid != NULL) {
      polarindex = indindex (indivmarkers, numindivs, polarid);
    }
    inddupcheck (indivmarkers, numindivs);
  }

  if (mkdiploid) {

    numindivs = rmindivs (snpmarkers, numsnps, indivmarkers, numindivs);
    numind2 = mkindh2d (indivmarkers, &indm2, numindivs);
    remaph2d (snpmarkers, numsnps, indivmarkers, indm2, numindivs, numind2);

    indivmarkers = indm2;
    numindivs = numind2;

  }

  if (mkhaploid) {

    numindivs = rmindivs (snpmarkers, numsnps, indivmarkers, numindivs);
    numind2 = mkindd2h (indivmarkers, &indm2, numindivs);
    remapd2h (snpmarkers, numsnps, indivmarkers, indm2, numindivs, numind2);

    indivmarkers = indm2;
    numindivs = numind2;


  }


  if (deletedup)
    dedupit (snpmarkers, numsnps);	// only one marker per position

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    if (zerodistance)
      cupt->genpos = 0.0;

    c1 = cupt->alleles[0];
    c2 = cupt->alleles[1];
    t1 = pedval (&c1) % 5;
    t2 = pedval (&c2) % 5;	// 0 and 5 are no good
    if ((t1 == 0) && (t2 > 0))
      flip1 (cupt, phasedmode, YES);
  }
  t = fixsnpdistance(snpmarkers, numsnps) ; 
  if (t>0) printf("%12d SNP positions adjusted\n", t) ;  

  if ((proboutfilename != NULL) && (numprobindivs > 0)) {
   plen2 = numsnps*numindivs*rl2 ; 
   ZALLOC(packp2, plen2, unsigned char) ;

   numx = loadprobpack(snpmarkers, indivmarkers, numsnps, numindivs,  (char *) packp2) ; 
  }
  

  if (deletedup)
    dedupit (snpmarkers, numsnps);	// only one marker per position

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    if (zerodistance)
      cupt->genpos = 0.0;

    c1 = cupt->alleles[0];
    c2 = cupt->alleles[1];
    t1 = pedval (&c1) % 5;
    t2 = pedval (&c2) % 5;	// 0 and 5 are no good
    if ((t1 == 0) && (t2 > 0))
      flip1 (cupt, phasedmode, YES);
  }

  flipstrand (flipstrandname, snpmarkers, numsnps);
  flipsnps (flipsnpname, snpmarkers, numsnps, phasedmode);

  if (polarindex >= 0) {
    for (i = 0; i < numsnps; i++) {
      cupt = snpmarkers[i];
      g = getgtypes (cupt, polarindex);
      if (g == 0) {
	printf ("polarizing %s", cupt->ID);
	printf (" %3d %12.0f", cupt->chrom, cupt->physpos);
	printnl ();
	fflush (stdout);
	flip1 (cupt, NO, YES);
	g = getgtypes (cupt, polarindex);
	if (g != 2)
	  fatalx ("badbug\n");
      }
      if (g != 2)
	cupt->ignore = YES;
    }
  }
  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    if (downsample)
      downsamp (cupt);
  }

  if (outputall) {
    numout = outfiles (snpoutfilename, indoutfilename, genooutfilename,
	      snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode);

    printf("numsnps output: %d\n", numout) ;

   if (proboutfilename != NULL) {
    printf("calling outprob\n") ; fflush(stdout) ; 
    outprobx(proboutfilename, snpmarkers, indivmarkers, numsnps, numindivs, (char *) packp2) ; 
    printf("PROB file %s written: %ld records\n", proboutfilename, numx) ;
   }

    printf ("##end of convertf run (outputall mode)\n");
    return 0;
  }


  if (usesamples != NULL) {
    poplistname = NULL;
    setsamp (indivmarkers, numindivs, usesamples);
  }

  if (fillmissingpoplistname != NULL) {  
   ZALLOC(fpops, numindivs, char *) ; 
   nfpops = loadlist(fpops, fillmissingpoplistname) ;
   t = fillmiss(snpmarkers, indivmarkers, numsnps, numindivs, fpops, nfpops) ;
   printf("%10d missing genotypes filled\n", t) ;
  }

  if (rmcompress) { 
   printf("before compress: snps: %d indivs: %d\n", numsnps, numindivs) ;
   numsnps = rmsnps (snpmarkers, numsnps, deletesnpoutname);
   numindivs = rmindivs (snpmarkers, numsnps, indivmarkers, numindivs);
   printf("after compress: snps: %d indivs: %d\n", numsnps, numindivs) ;
  }
  fflush(stdout) ;
  if (shuffle) doshuffle(snpmarkers, numsnps, numindivs) ;
  forcemiss(makemiss) ;
//  printf("got here! 2\n") ; fflush(stdout) ;

// force missing on wipeoutchrom.  retain snp.
  if (wipeoutchrom > 0)  printf("wiping out chrom: %d\n", wipeoutchrom) ; 
  for (i=0; i<numsnps; i++) { 
   if (wipeoutchrom<0) break ; 
   cupt = snpmarkers[i] ; 
   if (cupt -> chrom != wipeoutchrom) continue ; 
   fvalg(cupt, 999) ; // wipe out
  }

  if (killr2) {
    nkill =
      killhir2 (snpmarkers, numsnps, numindivs, r2physlim, r2genlim,
		r2thresh);
    if (nkill > 0)
      printf ("killhir2.  number of snps killed: %d\n", nkill);
  }


  if (nhwfilter > 0) {
    hwfilter (snpmarkers, numsnps, numindivs, nhwfilter, deletesnpoutname);
  }

  if (xregionname) {
    excluderegions (xregionname, snpmarkers, numsnps, deletesnpoutname);
  }


  numvalidind = 0;
  cputimes(0, 1) ;
  for (i = 0; i < numindivs; ++i) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    if (numvalidgtind (snpmarkers, numsnps, i) == 0) {
      if (probprintname == NULL) indx->ignore = YES;
      printf ("no data for individual: %s\n", indx->ID);
    }
    if (indx->ignore == NO)
      ++numvalidind;
  }

  if (maxmiss < 0)
    maxmiss = (int) (maxmissfrac * (double) numvalidind + 1);
  printf ("numvalidind:  %5d  maxmiss: %5d\n", numvalidind, maxmiss);

//   printf("clock 1: %9.3f\n", cputimes(1, 1)) ;

  if ((numvalidind == 0) && (probprintname ==  NULL))  
    fatalx ("no samples with valid genotypes!\n");

  cputimes(0, 2) ;
  t = testmisspop(snpmarkers, numsnps, indivmarkers, numindivs, minvalpop) ; 
  if (minvalpop>0) printf("minvalpop: deleted %d   retained %d\n", t, numsnps-t) ;

//   printf("clock 2: %9.3f\n", cputimes(1, 2)) ;

  cputimes(0, 3) ;
  for (k = 0; k < numsnps; ++k) {
    if (maxmiss > numvalidind)
      break;
    cupt = snpmarkers[k];
    t = numvalidind - numvalidgtypes (cupt);
    if (maxmiss < t) {
      cupt->ignore = YES;
    }

/**
   if (numvalidind ==  t) { 
    printf("no data for snp: %s\n", cupt -> ID) ;
    cupt -> ignore = YES ;
   }
*/

  }

//   printf("clock 3: %9.3f\n", cputimes(1, 3)) ;
  if (slowdup) { 
   fastdup = NO ; 
   fastdupthresh = 0.5 ;  fastdupkill = 2.0 ;
    printf ("calling slowdupcheck\n");
     setfastdupthresh (fastdupthresh, fastdupkill);
    for (k1=0; k1<numindivs; ++k1) { 
     if (indivmarkers[k1] -> ignore) continue ; 
    for (k2=k1+1; k2<numindivs; ++k2) { 
     if (indivmarkers[k2] -> ignore) continue ; 
     slowdupcheck (snpmarkers, indivmarkers, numsnps,  k1,  k2) ;
  }}} 


  if (fastdup) {

    printf ("fastdup set %d\n", fastdupnum);
    fflush (stdout);
    if (fastdupnum > 0) {
      setfastdupnum (fastdupnum);
      setfastdupthresh (fastdupthresh, fastdupkill);
      printf ("calling fastdupcheck\n");
      fflush (stdout);
      fastdupcheck (snpmarkers, indivmarkers, numsnps, numindivs);
    }
  }


  if (decim > 0) {
    snpdecimate (snpmarkers, numsnps, decim, dmindis, dmaxdis);
  }

  cputimes(0, 4) ;
  printf("callng outfiles\n") ;
  fflush(stdout) ;

  numout = outfiles (snpoutfilename, indoutfilename, genooutfilename,
	    snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode);

  if (tgenooutfilename != NULL) { 
    settrans(YES) ;
    outpack (tgenooutfilename, snpmarkers, indivmarkers, numsnps, numindivs);
  }


//   printf("clock 4: %9.3f\n", cputimes(1, 4)) ;
  printf("numsnps output: %d\n", numout) ;
  fflush(stdout) ;

/**
  x = snpindex(snpmarkers, numsnps, "rs12626123") ;
  cupt = snpmarkers[x] ; 
  printf("zzq2 %s ", cupt -> ID) ; 
  printmatl(cupt -> diplike[0], 1, 3) ;
*/

   if (proboutfilename != NULL) {
    if (numprobindivs <= 0) fatalx("proboutfilename set but not probindivname\n") ; 
    numx = loadprobpack(snpmarkers, indivmarkers, numsnps, numindivs, (char *) packp2) ; 
    outprobx(proboutfilename, snpmarkers, indivmarkers, numsnps, numindivs, (char *) packp2) ; 
    printf("PROB file %s written: %ld records\n", proboutfilename, numx) ;
    fflush(stdout) ;
   }
 if (probprintname != NULL)    {   
  if (numprobindivs <= 0) fatalx("probprintname set but not probindivname\n") ; 
  if (numprobindivs == 1) probprintid = probindivs[0] -> ID ;
  if (probprintid == NULL) fatalx("probprintid not set\n") ;  
  
  x = indindex(probindivs, numprobindivs, probprintid) ; 
  if (x<0) fatalx("probprintid: %s not in probindiv file\n", probprintid) ;  
  indx = probindivs[x] ;  
  openit(probprintname, &probfile, "w") ;  
  fprintf(probfile, "## prob for ID: %s file: %s\n", probprintid, indivname) ;

  for (i=0; i<numsnps; i++) {
   cupt = snpmarkers[i] ;
   if (cupt -> ignore) continue ;  
   fprintf(probfile, "%20s ", cupt -> ID) ; 
   fprintf(probfile, "%2d ",  cupt -> chrom) ; 
   fprintf(probfile, "%12.0f ", cupt -> physpos) ;
     k = indx -> idnum ;
     if (k<0) fatalx("badbug\n") ;  
     pp = cupt -> diplike[k] ; 
     printmatwxfile(pp, 1, 3, 3, probfile) ; 
     fprintf(probfile, "\n") ; 
  }
  fclose(probfile) ; 
 }

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of convertf: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0;
}

void
readcommands (int argc, char **argv)
{
  int i ;
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

  if (parname==NULL) { 
   printf("no parameter file (-p)\n") ; 
   exit(1) ;
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
  getstring (ph, "genotypename:", &genotypename);
  getstring (ph, "genotypelist:", &genotypelist);
  getstring (ph, "snpname:", &snpname);
  getstring (ph, "indivname:", &indivname);
  getstring (ph, "badsnpname:", &badsnpname);
  getstring (ph, "flipsnpname:", &flipsnpname);
  getstring (ph, "flipstrandname:", &flipstrandname);
  getstring (ph, "indoutfilename:", &indoutfilename);
  getstring (ph, "indivoutname:", &indoutfilename);	/* changed 11/02/06 */
  getstring (ph, "snpoutfilename:", &snpoutfilename);
  getstring (ph, "snpoutname:", &snpoutfilename);	/* changed 11/02/06 */
  getstring (ph, "genooutfilename:", &genooutfilename);
  getstring (ph, "tgenooutfilename:", &tgenooutfilename);
  getstring (ph, "genotypeoutname:", &genooutfilename);	/* changed 11/02/06 */
  getstring (ph, "tgenotypeoutname:", &tgenooutfilename) ;
  getstring (ph, "outputformat:", &omode);
  getstring (ph, "outputmode:", &omode);
  getstring (ph, "polarize:", &polarid);
  getstring (ph, "usesamples:", &usesamples);
  getint (ph, "zerodistance:", &zerodistance);
  getint (ph, "memorymap:", &memorymap);
  getint (ph, "mmap:", &memorymap);
  getint (ph, "checksizemode:", &checksizemode);
  getint (ph, "badpedignore:", &badpedignore);
  getint (ph, "downsample:", &downsample);
  getint (ph, "chimpmode:", &chimpmode);
  getint (ph, "pordercheck:", &pordercheck);
  getint (ph, "remapcheck:", &remapcheck);
  getint (ph, "seed:", &seed);
  getint (ph, "randommode:", &randommode);
  getint (ph, "familypopnames:", &familypopnames);
  getint (ph, "wipeoutchrom:", &wipeoutchrom);
  getint (ph, "rmcompress:", &rmcompress);

  getint (ph, "numchrom:", &numchrom);
  getstring (ph, "xregionname:", &xregionname);
  getdbl (ph, "hwfilter:", &nhwfilter);
  getstring (ph, "deletesnpoutname:", &deletesnpoutname);

  getint (ph, "outputgroup:", &ogmode);
  getint (ph, "malexhet:", &malexhet);
  getint (ph, "nomalexhet:", &malexhet);	/* changed 11/02/06 */
  getint (ph, "tersemode:", &tersem);
  getint (ph, "familynames:", &familynames);
  getint (ph, "packout:", &packout);	/* now obsolete 11/02/06 */
  getint (ph, "decimate:", &decim);
  getint (ph, "dmindis:", &dmindis);
  getint (ph, "dmaxdis:", &dmaxdis);
  getint (ph, "flipreference:", &flipreference);
  getint (ph, "fastdup:", &fastdup);
  getint (ph, "slowdup:", &slowdup);
  getint (ph, "fastdupnum:", &fastdupnum);
  getdbl (ph, "fastdupthresh:", &fastdupthresh);
  getdbl (ph, "fastdupkill:", &fastdupkill);
  getint (ph, "killr2:", &killr2);
  getint (ph, "hashcheck:", &hashcheck);
  getint (ph, "outputall:", &outputall);
  getint (ph, "sevencolumnped:", &sevencolumnped);
  getint (ph, "phasedmode:", &phasedmode);
  getint (ph, "polarcheck:", &polarcheck);
  getint (ph, "copyalleles:", &copyalleles);
// we assume with newsnpname we are (A,B) in S1 and (A, B) or (rev(A), rev(B)) in S2

  getdbl (ph, "r2thresh:", &r2thresh);
  getdbl (ph, "r2genlim:", &r2genlim);
  getdbl (ph, "r2physlim:", &r2physlim);

  getint (ph, "chrom:", &xchrom);
  getint (ph, "lopos:", &lopos);
  getint (ph, "hipos:", &hipos);

  getint (ph, "minchrom:", &minchrom);
  getint (ph, "maxchrom:", &maxchrom);
  getdbl (ph, "maxmissfrac:", &maxmissfrac);
  getint (ph, "maxmissing:", &maxmiss);
  getint (ph, "minvalpop:", &minvalpop);
  getdbl (ph, "makemiss:", &makemiss) ;

  getstring (ph, "poplistname:", &poplistname);
  getstring (ph, "indivlistname:", &indivlistname);
  getstring (ph, "instem:", &instem) ;
  getstring (ph, "outstem:", &outstem) ;
  getstring (ph, "fillmissingpoplistname:", &fillmissingpoplistname);
  getstring (ph, "newsnpname:", &newsnpname);
  getint (ph, "newignore:", &newignore);
  getstring (ph, "newindivname:", &newindivname);
  getint (ph, "deletedup:", &deletedup);
  getint (ph, "mkdiploid:", &mkdiploid);
  getint (ph, "mkhaploid:", &mkhaploid);
  getint (ph, "hiresgendis:", &hiresgendis);
  getstring (ph, "probfilename:", &probfilename);
  getstring (ph, "probname:", &probfilename);
  getstring (ph, "proboutfilename:", &proboutfilename);
  getstring (ph, "proboutname:", &proboutfilename);
  getstring (ph, "probtypeoutname:", &proboutfilename);
  getstring (ph, "probfileoutname:", &proboutfilename);
  getstring (ph, "probindivname:", &probindivname);
  getstring (ph, "probprint:", &probprintname) ;
  getstring (ph, "printprob:", &probprintname) ;
  getstring (ph, "probprintname:", &probprintname) ;
  getint (ph, "haploidfill:", &haploidfill) ;
  getint (ph, "transformhets:", &transformhets) ;
  getint (ph, "randomhetfix:", &randomhetfix) ;
  getint (ph, "shuffle:", &shuffle) ;

  writepars (ph);
  closepars (ph);

}

void
remap (SNP ** s1, int nums1, SNP ** s2, int nums2)
{
  SNP *cupt1, *cupt2;
  SNP tcupt, *cupt;
  int i, k;

  for (i = 0; i < nums2; ++i) {
    cupt2 = s2[i];
    k = snpindex (s1, nums1, cupt2->ID);
    if (k < 0) {
      printf ("%20s not found\n", cupt2->ID);
      if (newignore) {
	cupt2->ignore = YES;
      }
      continue;
    }
    cupt1 = s1[k];

    if (copyalleles) {
      cupt2->alleles[0] = cupt1->alleles[0];
      cupt2->alleles[1] = cupt1->alleles[1];
    }

    fixaa (cupt1, cupt2);
    if (cupt1->alleles[1] == 'X') {
      cupt1->alleles[1] = cxx (cupt1->alleles, cupt2->alleles);
    }
    tcupt = *cupt2;
    *cupt2 = *cupt1;
    cupt = &tcupt;
    cupt2->chrom = cupt->chrom;
    cupt2->genpos = cupt->genpos;
    cupt2->physpos = cupt->physpos;
    cupt2->alleles[0] = cupt->alleles[0];
    cupt2->alleles[1] = cupt->alleles[1];
  }
  freesnpindex ();
}

char
cxx (char *c1, char *c2)
{
  if (c1[0] == c2[0])
    return c2[1];
  if (c1[0] == c2[1])
    return c2[0];

  return 'X' ; 

}

void
fixaa (SNP * cupt1, SNP * cupt2)
{
  char *c1, *c2;
  int t, ok;
  char cc1, cc2;

  if (remapcheck == NO)
    return;
  c1 = cupt1->alleles;
  c2 = cupt2->alleles;

  if ((c1[1] == 'X') && (c1[0] == c2[0]))
    c1[1] = c2[1];
  if ((c1[1] == 'X') && (c1[0] == c2[1]))
    c1[1] = c2[0];
  if ((c2[1] == 'X') && (c1[0] == c2[0]))
    c2[1] = c1[1];
  if ((c2[1] == 'X') && (c1[1] == c2[0]))
    c2[1] = c1[0];
  if ((c1[0] == c2[0]) && (c1[1] == c2[1]))
    return;
  if (polarcheck) {
    ok = YES;
    cc1 = toupper (c1[0]);
    cc2 = toupper (c2[0]);
    if (cc2 != revchar (cc1))
      ok = NO;
    cc1 = toupper (c1[1]);
    cc2 = toupper (c2[1]);
    if (cc2 != revchar (cc1))
      ok = NO;
    if (ok == NO) {
      printf ("forcing all genos invalid for %s %c %c     %c %c\n", cupt1->ID,
	      c1[0], c1[1], c2[0], c2[1]);
      fvalg (cupt1, 999);
      cupt1->ignore = YES;
      return;
    }
    c1[0] = c2[0];
    c1[1] = c2[1];
  }

  if ((c1[0] == c2[1]) && (c1[1] == c2[0])) {
    flip1 (cupt1, phasedmode, YES);
    return;
  }
  t = 999;
  if ((c1[0] == c2[0]) && (c1[1] != c2[1])) {
    t = 2;
    if (phasedmode)
      t = 1;
  }
  if ((c1[1] == c2[1]) && (c1[0] != c2[0])) {
    t = 0;
    fvalg (cupt1, 0);		// only valid genotype ;
    return;
  }
  fvalg (cupt1, t);		// force all snps invalid
  if (t == 999) {
    printf ("forcing all genos invalid for %s %c %c     %c %c\n", cupt1->ID,
	    c1[0], c1[1], c2[0], c2[1]);
  }
}

void
forcemiss (double yprob)
{
  int i, k, g, t, g2;
  static int ncall = 0;
  SNP *cupt ; 

  if (yprob<=0) return ; 
  ++ncall;

 for (i=0; i<numsnps; ++i) { 
  cupt = snpmarkers[i] ; 
  if (cupt -> ignore) continue ;
  for (k = 0; k < numindivs; ++k) {
   t = prob1(yprob) ; 
   if (t==1) putgtypes(cupt, k, -1) ;
  }
 }
}

void
downsamp (SNP * cupt)
{
  int k, g, t, g2;
  static int ncall = 0;

  ++ncall;

  for (k = 0; k < numindivs; ++k) {
    g = getgtypes (cupt, k);
    if (g == 1) {
      t = ranmod (2);
      putgtypes (cupt, k, 2 * t);
    }
  }
}

void
fvalg (SNP * cupt, int val)
{
  int k, g;

  for (k = 0; k < numindivs; ++k) {
    g = getgtypes (cupt, k);
    if (g != val)
      putgtypes (cupt, k, -1);
  }
}


void
remapind (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	  Indiv ** indm2, int numindivs, int numind2)
{

  int *g1, *g2, *w1;
  int *tind, t, i, j, k;
  Indiv *indx;
  SNP *cupt;
  double **probbuff ; 

  if (numprobindivs>0) { 
   probbuff = initarray_2Ddouble(numind2, 3, 0) ; 
  }


  if ((numindivs != numind2) && (remapcheck == YES)) 
    fatalx ("different remapind sizes %d %d\n", numindivs, numind2);
  ZALLOC (tind, numind2, int);
  ZALLOC (g2, numindivs, int);
  ZALLOC (g1, numindivs, int);
  ZALLOC (w1, numindivs, int);

  for (k = 0; k < numind2; ++k) {
    indx = indm2[k];
    t = tind[k] = indindex (indivmarkers, numindivs, indx->ID);
    if (t < 0)
      fatalx ("bad newindiv: %s\n", indx->ID);
  }

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];

    for (j = 0; j < numindivs; ++j) {
      g1[j] = getgtypes (cupt, j);
    }
    copyiarr (g1, w1, numindivs);

    for (k = 0; k < numind2; ++k) {
      g2[k] = g1[tind[k]];
      if (numprobindivs>0) copyarr(cupt -> diplike[tind[k]], probbuff[k], 3) ; 
    }

    copyiarr (g2, g1, numind2);


    

    for (k = 0; k < numind2; ++k) {
      putgtypes (cupt, k, g1[k]);
      if (numprobindivs>0) copyarr(probbuff[k], cupt -> diplike[k], 3) ; 
    }

  }

  free (w1);
  free (g1);
  free (g2);
  free (tind);

  if (numprobindivs>0) { 
   free2D(&probbuff, numind2) ;
  }


}



void
dedupit (SNP ** snpmarkers, int numsnps)
{
  SNP *cupt1, *cupt2;
  SNP *x1, *x2;
  int k;

  cupt1 = NULL;

  for (k = 0; k < numsnps; ++k) {
    cupt2 = snpmarkers[k];
    if (cupt2->ignore)
      continue;
    if (cupt1 == NULL) {
      cupt1 = cupt2;
      continue;
    }
    if (cupt1->chrom != cupt2->chrom) {
      cupt1 = cupt2;
      continue;
    }
    if (cupt1->physpos != cupt2->physpos) {
      cupt1 = cupt2;
      continue;
    }
    pickx (cupt1, cupt2, &x1, &x2);	// x2 bad
    x2->ignore = YES;
    cupt1 = x1;
  }
}

void
pickx (SNP * c1, SNP * c2, SNP ** px1, SNP ** px2)
{
// *px1 is retained *px2 dropped; Try and keep shorter rsnames (strip AFFX for instance)

  char *ch1, *ch2;
  int l1, l2, t;

  ch1 = c1->ID;
  ch2 = c2->ID;
  l1 = strlen (ch1);
  l2 = strlen (ch2);

  if (l1 < l2) {
    *px2 = c2;
    *px1 = c1;
    return;
  }
  if (l2 < l1) {
    *px2 = c1;
    *px1 = c2;
    return;
  }
  t = strcmp (ch1, ch2);
  if (t >= 0) {
    *px2 = c2;
    *px1 = c1;
    return;
  }
  *px2 = c1;
  *px1 = c2;
  return;
}

void
flipstrand (char *fsname, SNP ** snpm, int numsnps)
// move alleles to opposite strand
{
  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXFF];
  char *ss;
  int nsplit, n, k;
  SNP *cupt;

  if (fsname == NULL)
    return;
  openit (fsname, &fff, "r");

  freesnpindex ();

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    if (spt[0][0] == '#') {
      freeup (spt, nsplit);
      continue;
    }

    k = snpindex (snpm, numsnps, spt[0]);
    if (k >= 0) {
      cupt = snpm[k];
      cupt->alleles[0] = compbase (cupt->alleles[0]);
      cupt->alleles[1] = compbase (cupt->alleles[1]);
    }
    freeup (spt, nsplit);
  }
  fclose (fff);
}

void
flipsnps (char *fsname, SNP ** snpm, int numsnps, int phasedmode)
{
  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXFF];
  char *ss;
  int nsplit, n, k;
  SNP *cupt;

  if (fsname == NULL)
    return;
  openit (fsname, &fff, "r");

  freesnpindex ();

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    if (spt[0][0] == '#') {
      freeup (spt, nsplit);
      continue;
    }

    k = snpindex (snpm, numsnps, spt[0]);
    if (k >= 0) {
      flip1 (snpm[k], phasedmode, flipreference);
    }
    freeup (spt, nsplit);
  }
  fclose (fff);
}

void
flip1 (SNP * cupt, int phasedmode, int flipreference)
{
  if (phasedmode == NO)
    flipalleles (cupt);
  if (phasedmode == YES)
    flipalleles_phased (cupt);
// just flips genotypes  
  if (flipreference)
    cswap (&cupt->alleles[0], &cupt->alleles[1]);
}

void
setsamplist (Indiv ** indivmarkers, int numindivs, char **samplist,
	     int nsplit)
{
  int i, j, k;
  Indiv *indx;
  int t = 0;

  setstatusv (indivmarkers, numindivs, NULL, NO);	// set affstatus to NO

  for (j = 0; j < nsplit; ++j) {
    k = indindex (indivmarkers, numindivs, samplist[j]);
    if (k < 0) {
      printf ("*** warning: sample ID: %s not found\n", samplist[j]);
      continue;
    }
    indx = indivmarkers[k];
    indx->affstatus = YES;
  }

  for (i = 0; i < numindivs; ++i) {
    indx = indivmarkers[i];
    if (indx->affstatus == NO)
      indx->ignore = YES;
    if (indx->ignore == NO)
      ++t;
  }

  if (t == 0)
    fatalx ("(setsamplist) no valids\n");
}

int
setsamp (Indiv ** indivmarkers, int numindivs, char *usesamples)
{
  char *spt[MAXFF];
  int nsplit;

  nsplit = splitupx (usesamples, spt, MAXFF, ':');

  setsamplist (indivmarkers, numindivs, spt, nsplit);

  freeup (spt, nsplit);

  return nsplit;

}
int mkpops(char **pops, Indiv **indm, int numind) 
// make a list of pops from indm
{
  int n = 0 ;
  Indiv *indx ; 
  int k, t ; 
  char *sx ; 
 
  for (k=0; k<numind; ++k) { 
   indx = indivmarkers[k] ;
   if (indx -> ignore) continue ;
   sx = indx -> egroup ; 
   t = indxstring(pops, n, sx) ; 
   if (t>=0) continue ; 
   pops[n] = strdup(sx) ;
   ++n ; 
  }
// make a list of pops from indm

  return n ; 


}
int
testmisspop(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, int numindivs, int minvalpops)  
{
    
  char **pops ; 
  int  npops, ndelete = 0 ;
  int **valcnt, i, k, x, g ; 
  Indiv *indx ; 
  SNP *cupt ;

  if (minvalpops<1) return 0 ; 
  ZALLOC(pops, numindivs, char *) ; 
  npops = mkpops(pops, indivmarkers, numindivs) ;
  if (npops == 0) fatalx("no valid data\n") ;
  if (verbose) {
   printf("zzpops:\n") ; printstrings(pops, npops) ;
  }

  valcnt = initarray_2Dint(numsnps, npops, 0) ; 
  for (i=0; i<numindivs; ++i) { 
   indx = indivmarkers[i] ; 
   if (indx -> ignore) continue ; 
   x = indxstring(pops, npops, indx -> egroup) ; 
   if (x<0) continue ; 
   for (k=0; k<numsnps; ++k) { 
    cupt = snpmarkers[k] ; 
    if (cupt -> ignore) continue ; 
    g = getgtypes(cupt, i) ; 
    if (g<0) continue ; 
    ++valcnt[k][x] ; 
   }
  }
  for (k=0; k<numsnps; ++k) { 
   if (verbose && (k<10))  printimat(valcnt[k], 1, npops) ;
   ivmaxmin(valcnt[k], npops, NULL, &x) ; 
   if (x<minvalpops) { 
     snpmarkers[k] -> ignore = YES ; 
     ++ndelete ; 
   }
  }

 freeup(pops, numindivs) ; 
 free2Dint(&valcnt, numsnps) ;
 return ndelete ; 

}

int count_column(SNP *cupt, int *xtypes, int numindivs, int **ncount, int *nmiss, int npops) 
// ncount npops x 3     nmiss ncount long   must be pre-allocated
{
    int i, k, g  ;

    if (cupt -> ignore)  return -1 ;

    iclear2D(&ncount, npops, 3, 0) ;
    ivzero(nmiss, npops) ;

    for (i=0; i<numindivs; ++i) { 
     k = xtypes[i]  ; 
     if (k<0) continue ; 
     if (k>= npops) fatalx("(count_column) overflow from xtypes\n") ;  
     g = getgtypes(cupt, i) ;
     if (g<0) {    
      ++nmiss[k] ;
      continue ;
     }
     ++ncount[k][g] ; 
    }

   return intsum(nmiss, npops) ; 

}

int fillmiss(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char **fpops, int nfpops) 

{

   int *xtypes, i, j, k, g, tmiss, g1, g2, l ; 
   int top, bot ;
   Indiv *indx ; 
   SNP *cupt ;
   int **ncount, *nmiss, *nc ;
   int nfill = 0 ;
   double y ; 

   if (nfpops == 0) return 0 ; 
   if (fpops == NULL) return 0 ; 

   ZALLOC(xtypes, numindivs, int) ; 
   ZALLOC(nmiss, nfpops, int) ; 
   ivclear(xtypes, -3, numindivs) ;  

   for (i=0; i<numindivs; ++i) {  
    indx = indivmarkers[i] ; 
    k = indxindex(fpops, nfpops, indx -> egroup) ; 
    xtypes[i] = k ; 
   }

   ncount = initarray_2Dint(nfpops, 3, 0) ;  // 0, 1 2  count

   for (j=0; j<numsnps; ++j) { 
    cupt = snpmarkers[j] ; 
    if (cupt -> ignore) continue ;

    tmiss  = count_column(cupt, xtypes, numindivs, ncount, nmiss, nfpops) ;

    for (k=0; k<nfpops; ++k) {  
     if (cupt -> ignore) break ; 
     if (tmiss==0) break ; 
     nc = ncount[k] ; 
     bot = 2*intsum(nc, 3) ;         
     if (bot==0) {  
      continue ;
     }
    if (nmiss[k] <= 0) continue ;  // no fill in needed
     top = 2*nc[2] + nc[1] ; 
     y = (double) top / (double) bot ; 
// now fill in 
     for (i=0; i<numindivs; ++i) { 
      if (xtypes[i] != k) continue ; 
      g = getgtypes(cupt, i) ;
      if (g>=0) continue ; 
      g1 = prob1(y) ; 
      if (haploidfill) g2 = g1 ; 
      else  g2 = prob1(y) ; 
      putgtypes(cupt, i, g1 + g2) ;
      ++nfill ;
     }
    } 
   }

   free2Dint(&ncount, nfpops) ; 
   free(xtypes) ;
   free(nmiss) ;

   return nfill ;
}
int fixsnpdistance(SNP **snpm, int numsnps) 
{
  int k, n=0  ; 
  SNP *cupt1, *cupt2 ;
  double dis ; 

  if (hiresgendis == NO) return 0 ; 
  
  for (k=1; k<numsnps; ++k) { 
   cupt1 = snpm[k-1] ; 
   cupt2 = snpm[k] ; 
   if (cupt2 -> chrom != cupt1 -> chrom) continue ;  
   dis = cupt2 -> genpos - cupt1 -> genpos ; 
   if (fabs(dis) < 1.0e-8) { 
    cupt2 -> genpos = cupt1 -> genpos + 1.0e-9 ; 
    ++n ; 
   }
  }
  return n ; 
}

long loadprobpack(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char *bigbuff) 

{
  SNP *cupt ;
  Indiv *indx ; 
  int rl2, i, j, x ;  
  int sval;
  double yy, ww[3] ;   
  unsigned char *buff ; 
  unsigned short bb[2] ; 
  long numx = 0 ;

  
  buff = (unsigned char *) bigbuff ;
  rl2 = 4 ; 
  sval = (1 << 16 ) -1 ; 

  for (i=0; i<numsnps; i++) { 
   cupt = snpmarkers[i] ; 
   if (cupt -> ignore) continue ; 
   cupt -> scount = 1 ; 
// printf("zz %3d %s %d %x\n", i, cupt -> ID, numprobindivs, cupt -> diplike) ;  fflush(stdout) ; 
   for (j=0; j<numindivs; ++j) {
     indx = indivmarkers[j] ;  
     if (indx -> ignore) continue ; 
//   printf("zzq %d %s %x %x %9.3f\n", j, indx -> ID, cupt -> diplike, cupt-> diplike[j], cupt->diplike[j][0]) ;  fflush(stdout) ;
     copyarr(cupt -> diplike[j], ww, 3) ;
     bb[0] = bb[1] = sval ; 
     if (ww[0] > -0.5) {     
      bal1(ww, 3) ; 
      yy = (double) sval * ww[0] ;  x = nnint(yy) ; bb[0] = (unsigned short) x ;
      yy = (double) sval * ww[2] ;  x = nnint(yy) ; bb[1] = (unsigned short) x ;
     }
     memcpy(buff, bb, rl2) ;
     buff += rl2 ;  
     ++numx ;
  }}

  return numx ; 

}

void doshuffle(SNP **snpm, int numsnps, int numindivs)  
// permute genetypes destroying LD + pop structure.  Retain inbreed 
{
  int *a, n = numindivs, i, k, g ; 
  SNP *cupt ;

  printf("doing shuffle!\n") ;

  ZALLOC(a, n, int) ;
  for (i=0; i< numsnps; ++i) { 
   cupt = snpm[i] ; 

   for (k=0; k<n; ++k) { 
    a[k] = getgtypes(cupt, k) ;
   }
   ranperm(a, n) ; 
   for (k=0; k<n; ++k) { 
     putgtypes(cupt, k, a[k]) ;
   }
 
  }

  free(a) ;
}
