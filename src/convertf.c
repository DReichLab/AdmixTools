#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h>

#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "egsubs.h"
#include "exclude.h"

#define WVERSION   "5000"
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
*/


#define MAXFL  50
#define MAXSTR  512

char *trashdir = "/var/tmp";
int qtmode = NO;
Indiv **indivmarkers, **indm2;
SNP **snpmarkers;
SNP **snpm2;
int zerodistance = NO;		// YES => force gdis 0
int downsample = NO;		// make pseudo homozygotes
int pordercheck = YES;
int familypopnames = NO;

int numsnps, numindivs, numind2;
int nums2;

char *genotypename = NULL;
char *genotypelist = NULL;

char *snpname = NULL;
char *indoutfilename = NULL;
char *snpoutfilename = NULL;
char *genooutfilename = NULL;
char *indivname = NULL;
char *newindivname = NULL;
char *badsnpname = NULL;
char *xregionname = NULL;
char *deletesnpoutname = NULL;
char *flipsnpname = NULL;
char *flipstrandname = NULL;
int flipreference = YES;
int remapcheck = YES;

char *poplistname = NULL;
char *fillmissingpoplistname = NULL ; 
char **fpops ;
int nfpops = 0 ;

double r2thresh = -1.0;
double r2genlim = 0.01;		// Morgans 
double r2physlim = 5.0e6;
double maxmissfrac = 1000.0;	// no thresh
int maxmiss = -1;		// no thresh
int minvalpop = -1 ;
int killr2 = NO;
int mkdiploid = NO;

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

int deletedup = YES;		// only one marker at a position
char *newsnpname = NULL;	// new map  
int newignore = YES;		// default ignore snps not in old list
int polarcheck = NO;
int copyalleles = NO;

char *usesamples = NULL;

char unknowngender = 'U';
double nhwfilter = -1;


void setomode (enum outputmodetype *outmode, char *omode);
void readcommands (int argc, char **argv);
void outfiles (char *snpname, char *indname, char *gname, SNP ** snpm,
	       Indiv ** indiv, int numsnps, int numind, int packem,
	       int ogmode);
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
int setsamp (Indiv ** indivmarkers, int numindivs, char *usesamples);
int testmisspop(SNP **snpmarkers, int numsnps, Indiv **indivmarkers, int numindivs, int minvalpops)   ;
int fillmiss(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, char **fpops, int nfpops)  ;

int usage (char *prog, int exval); 

int
usage (char *prog, int exval)
{
  
  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");

  exit(exval);
};

int
main (int argc, char **argv)
{

  int **snppos;
  int *snpindx;
  char **snpnamelist, **indnamelist;
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

  int nindiv = 0, e, f, lag = 1;
  double xc[9], xd[4], xc2[9];
  double ychi, zscore, zthresh = 20.0;
  double y1, y2;
  int nignore, numrisks = 1;

  char **genolist;
  int numgenolist;
  char c1, c2;
  int t1, t2;

  malexhet = YES;		// convertf default is don't change the data
  tersem = YES;			// no snp counts

  readcommands (argc, argv);

  
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

  setomode (&outputmode, omode);
  packmode = YES;
  settersemode (tersem);

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

  inddupcheck (indivmarkers, numindivs);

  if (genotypelist != NULL) {
    getgenos_list (genotypelist, snpmarkers, indivmarkers,
		   numsnps, numindivs, nignore);
  }

  else {
    setgenotypename (&genotypename, indivname);
    getgenos (genotypename, snpmarkers, indivmarkers,
	      numsnps, numindivs, nignore);
  }

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
    outfiles (snpoutfilename, indoutfilename, genooutfilename,
	      snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode);

    printf ("##end of convertf run (outputall mode)\n");
    return 0;
  }

  if (usesamples != NULL) {
    poplistname = NULL;
    setsamp (indivmarkers, numindivs, usesamples);
  }

  if (fillmissingpoplistname != NULL) {  
   printf("seed: %d\n", seed) ;
   ZALLOC(fpops, numindivs, char *) ; 
   nfpops = loadlist(fpops, fillmissingpoplistname) ;
   t = fillmiss(snpmarkers, indivmarkers, numsnps, numindivs, fpops, nfpops) ;
   printf("%10d missing genotypes filled\n", t) ;
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
  else
    setstatus (indivmarkers, numindivs, "Case");


  numsnps = rmsnps (snpmarkers, numsnps, deletesnpoutname);
  numindivs = rmindivs (snpmarkers, numsnps, indivmarkers, numindivs);
//  printf("got here! 2\n") ; fflush(stdout) ;

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
  for (i = 0; i < numindivs; ++i) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    if (numvalidgtind (snpmarkers, numsnps, i) == 0) {
      indx->ignore = YES;
      printf ("no data for individual: %s\n", indx->ID);
    }
    if (indx->ignore == NO)
      ++numvalidind;
  }

  if (maxmiss < 0)
    maxmiss = (int) (maxmissfrac * (double) numvalidind + 1);
  printf ("numvalidind:  %5d  maxmiss: %5d\n", numvalidind, maxmiss);
  if (numvalidind == 0)
    fatalx ("no valid samples!\n");

  t = testmisspop(snpmarkers, numsnps, indivmarkers, numindivs, minvalpop) ; 
  if (minvalpop>0) printf("minvalpop: deleted %d   retained %d\n", t, numsnps-t) ;


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

//  printf("got here! 3\n") ; fflush(stdout) ;

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

  outfiles (snpoutfilename, indoutfilename, genooutfilename,
	    snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode);

  printf ("##end of convertf run\n");
  return 0;
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
  if (argc == 1) { usage(basename(argv[0]), 1); }
  while ((i = getopt (argc, argv, "hp:vV")) != -1) {

    switch (i) {

    case 'h':
      usage(basename(argv[0]), 0);

    case 'p':
      parname = strdup (optarg);
      break;

    case 'v':
      printf ("version: %s\n", WVERSION);
      exit(0);

    case 'V':
      verbose = YES;
      break;

    default:
        usage(basename(argv[0]), 1);
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
  getstring (ph, "genotypeoutname:", &genooutfilename);	/* changed 11/02/06 */
  getstring (ph, "outputformat:", &omode);
  getstring (ph, "outputmode:", &omode);
  getstring (ph, "polarize:", &polarid);
  getstring (ph, "usesamples:", &usesamples);
  getint (ph, "zerodistance:", &zerodistance);
  getint (ph, "checksizemode:", &checksizemode);
  getint (ph, "badpedignore:", &badpedignore);
  getint (ph, "downsample:", &downsample);
  getint (ph, "chimpmode:", &chimpmode);
  getint (ph, "pordercheck:", &pordercheck);
  getint (ph, "remapcheck:", &remapcheck);
  getint (ph, "seed:", &seed);
  getint (ph, "randommode:", &randommode);
  getint (ph, "familypopnames:", &familypopnames);

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

  getstring (ph, "poplistname:", &poplistname);
  getstring (ph, "fillmissingpoplistname:", &fillmissingpoplistname);
  getstring (ph, "newsnpname:", &newsnpname);
  getint (ph, "newignore:", &newignore);
  getstring (ph, "newindivname:", &newindivname);
  getint (ph, "deletedup:", &deletedup);
  getint (ph, "mkdiploid:", &mkdiploid);

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
downsamp (SNP * cupt)
{
  int k, g, t, g2;
  static int ncall = 0;

  ++ncall;
  if (ncall == 1) {
    SRAND (77);
    printf ("downsample set\n");
  }

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
    }

    copyiarr (g2, g1, numind2);

    for (k = 0; k < numind2; ++k) {
      putgtypes (cupt, k, g1[k]);
    }

  }

  free (w1);
  free (g1);
  free (g2);
  free (tind);

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
      cupt -> ignore = YES;  
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
      g2 = prob1(y) ; 
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
