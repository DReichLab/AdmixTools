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

#define WVERSION   "2450"
#define MAXFL  50
#define MAXSTR  512

// phasedmode added
// bugfix for pedcols
// mcio (bigread)
// mergeit (overflow bug in clear) 
// and better code for mergeit routine (overflow bug in cclear)
// malexhet added
// polarmode added both input files assumed polarized

char *trashdir = "/var/tmp";
extern int verbose;
int qtmode = NO;
Indiv **indm1, **indm2;
SNP **snpm1, **snpm2;
int nums1, nums2;
int numi1, numi2;

char *snp1 = NULL;
char *ind1 = NULL;
char *geno1 = NULL;

char *snp2 = NULL;
char *ind2 = NULL;
char *geno2 = NULL;

char *indoutfilename = NULL;
char *snpoutfilename = NULL;
char *genooutfilename = NULL;
char *badsnpname = NULL;

int packout = -1;
int tersem = YES;
extern enum outputmodetype outputmode;
extern int checksizemode;
char *omode = "packedancestrymap";
extern int packmode;
extern long packlen;

int ogmode = NO;
int docheck = YES;
int strandcheck = YES;
int polarmode = NO;
int phasedmode = NO;
int allowdups = NO;

int xchrom = -1;
int lopos = -999999999;
int hipos = 999999999;
int minchrom = 1;
int maxchrom = 97;

/** 
 docheck:  YES  allele flipping check
 hashcheck:  YES  ... NO => snp names, Indiv names changed MUST retain order and number
 strandcheck: (default YES) if NO then alleles are assumed on same strand
 allowdups: YES (second set of data ... dup individuals are killed otherwise fatal error
 rewrote mergeit()  removed call to rmindivs (bug??)
*/

char unknowngender = 'U';

void setomode (enum outputmodetype *outmode, char *omode);
void readcommands (int argc, char **argv);
void outfiles (char *snpname, char *indname, char *gname, SNP ** snpm,
               Indiv ** indiv, int numsnps, int numind, int packem,
               int ogmode);
int checkmatch (SNP * cupt1, SNP * cupt2);

int
mergeit (SNP ** snpm1, SNP ** snpm2, Indiv *** pindm1, Indiv ** indm2,
         int nums1, int nums2, int numi1, int numi2);


int
main (int argc, char **argv)
{
  SNP **snpmarkers;
  Indiv **indivmarkers;
  int numsnps, numindivs;
  unsigned char *packg1, *packg2;

  int **snppos;
  int *snpindx;
  int lsnplist, lindlist, numeg;
  int i, j;
  SNP *cupt, *cupt1, *cupt2, *cupt3;
  Indiv *indx;

  int ch1, ch2;
  int fmnum, lmnum;
  int num, n1, n2;
  int nkill = 0;
  int t, k, x;

  int nignore, numrisks = 1;

  char **genolist;
  int numgenolist;
  int maxmiss;
  int sorder[2];

  tersem = YES;                 // no snp counts

  readcommands (argc, argv);

  setomode (&outputmode, omode);
  packmode = YES;
  settersemode (tersem);
  if (phasedmode)
    malexhet = YES;
  if (polarmode) {
    strandcheck = NO;
    printf ("polarmode set!\n");
  }

  nums1 = getsnps (snp1, &snpm1, 0.0, NULL, &nignore, numrisks);

  sorder[0] = getsnpordered ();
  putped (1);
  freeped ();

  sorder[1] = getsnpordered ();
  nums2 = getsnps (snp2, &snpm2, 0.0, NULL, &nignore, numrisks);

  putped (2);
  freeped ();

  for (x = 0; x < nums1; ++x) {
    cupt1 = snpm1[x];
    cupt1->tagnumber = -1;
  }
  for (x = 0; x < nums2; ++x) {
    cupt2 = snpm2[x];
    t = x % 1000;
// if (t==0) printf("zz %d %d\n", x, nums2) ;

    k = snpindex (snpm1, nums1, cupt2->ID);
    if (k < 0) {
      cupt2->ignore = YES;
      continue;
    }
    cupt1 = snpm1[k];
    cupt1->tagnumber = x;
    t = checkmatch (cupt1, cupt2);
    if (t == 1)
      continue;
    if (t == 2) {
      cupt2->isrfake = YES;
      continue;
    }
    if (t < 0) {
      cupt1->ignore = cupt2->ignore = YES;
      continue;
    }
    printf ("allele funny: %s", cupt1->ID);
    printalleles (cupt1, stdout);
    printalleles (cupt2, stdout);
    printnl ();
    cupt1->ignore = cupt2->ignore = YES;
    continue;
  }
  freesnpindex ();
  numi1 = getindivs (ind1, &indm1);
  numi2 = getindivs (ind2, &indm2);

  for (x = 0; x < numi2; ++x) {
    k = indindex (indm1, numi1, indm2[x]->ID);
    if ((k >= 0) && (allowdups == NO))
      fatalx ("dup ind: %s\n", indm2[x]->ID);   // fix later?  
    if ((k >= 0) && (allowdups) && (indm1[k]->ignore == NO))
      indm2[x]->ignore = YES;
  }


  setgenotypename (&geno1, ind1);
  getped (1);
  putsnpordered (sorder[0]);
  getgenos (geno1, snpm1, indm1, nums1, numi1, nignore);

  packg1 = (unsigned char *) getpackgenos ();
  clearpackgenos ();

  setgenotypename (&geno2, ind2);
  getped (2);
  putsnpordered (sorder[1]);
  getgenos (geno2, snpm2, indm2, nums2, numi2, nignore);

/**
  numi1 = rmindivs(snpm1, nums1, indm1, numi1)  ; 
  numi2 = rmindivs(snpm2, nums2, indm2, numi2)  ; 
*/

  packg2 = (unsigned char *) getpackgenos ();
  numindivs =
    mergeit (snpm1, snpm2, &indm1, indm2, nums1, nums2, numi1, numi2);

  snpmarkers = snpm1;
  numsnps = nums1;
  indivmarkers = indm1;

  free (packg1);
  free (packg2);

//  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;

  printf ("numsnps: %d  numindivs: %d\n", numsnps, numindivs);

  outfiles (snpoutfilename, indoutfilename, genooutfilename,
            snpmarkers, indivmarkers, numsnps, numindivs, packout, ogmode);

  printf ("##end of mergeit run\n");
  return 0;
}

int
checkmatch (SNP * cupt1, SNP * cupt2)
{

  char a1, a2, b1, b2;

  if (docheck == NO)
    return 1;
  if (cupt1->alleles == NULL)
    return -1;
  if (cupt2->alleles == NULL)
    return -1;

  a1 = cupt1->alleles[0];
  a2 = cupt1->alleles[1];

  b1 = cupt2->alleles[0];
  b2 = cupt2->alleles[1];

  a1 = toupper (a1);
  a2 = toupper (a2);

  b1 = toupper (b1);
  b2 = toupper (b2);

  if ((a1 == 'X') && (a2 == 'X') && (strandcheck))
    return -1;                  // flipcheck impossible

  if (strandcheck) {

    if ((a1 == 'A') && (a2 == 'T'))
      return -1;
    if ((a1 == 'T') && (a2 == 'A'))
      return -1;
    if ((a1 == 'C') && (a2 == 'G'))
      return -1;
    if ((a1 == 'G') && (a2 == 'C'))
      return -1;

    if ((b1 == 'A') && (b2 == 'T'))
      return -1;
    if ((b1 == 'T') && (b2 == 'A'))
      return -1;
    if ((b1 == 'C') && (b2 == 'G'))
      return -1;
    if ((b1 == 'G') && (b2 == 'C'))
      return -1;
  }

  if (polarmode) {
    if ((a1 == b1) && (a2 == b2))
      return 1;
    if ((a1 == b1) && (b2 == 'X'))
      return 1;
    b1 = compbase (b1);
    b2 = compbase (b2);
    if ((a1 == b1) && (a2 == b2))
      return 1;
    return -1;
  }

  if ((a1 == b1) && (a2 == 'X')) {
    cupt1->alleles[1] = b2;
    return 1;
  }

  if ((a1 == b2) && (a2 == 'X')) {
    cupt1->alleles[1] = b1;
    return 2;
  }


  if ((a1 == b1) && (a2 == b2))
    return 1;
  if ((a1 == b2) && (a2 == b1))
    return 2;

  if ((a1 == b1) && (b2 == 'X'))
    return 1;
  if ((a2 == b1) && (b2 == 'X'))
    return 2;

  if (strandcheck == NO)
    return 0;

  b1 = compbase (b1);
  b2 = compbase (b2);

  if ((a1 == b1) && (a2 == 'X')) {
    cupt1->alleles[1] = b2;
    return 1;
  }

  if ((a1 == b2) && (a2 == 'X')) {
    cupt1->alleles[1] = b1;
    return 2;
  }

  if ((a1 == b1) && (a2 == b2))
    return 1;
  if ((a1 == b2) && (a2 == b1))
    return 2;

  if ((a1 == b1) && (b2 == 'X'))
    return 1;
  if ((a2 == b1) && (b2 == 'X'))
    return 2;

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

  while ((i = getopt (argc, argv, "p:vVf")) != -1) {

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

    case 'f':
      phasedmode = YES;
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

  getstring (ph, "geno1:", &geno1);
  getstring (ph, "snp1:", &snp1);
  getstring (ph, "ind1:", &ind1);

  getstring (ph, "geno2:", &geno2);
  getstring (ph, "snp2:", &snp2);
  getstring (ph, "ind2:", &ind2);

  getstring (ph, "indoutfilename:", &indoutfilename);
  getstring (ph, "snpoutfilename:", &snpoutfilename);
  getstring (ph, "genooutfilename:", &genooutfilename);
  getstring (ph, "outputformat:", &omode);

  getint (ph, "malexhet:", &malexhet);
  getint (ph, "nomalexhet:", &malexhet);        /* changed 11/02/06 */

  getint (ph, "docheck:", &docheck);
  getint (ph, "hashcheck:", &hashcheck);
  getint (ph, "strandcheck:", &strandcheck);
  getint (ph, "phasedmode:", &phasedmode);
  getint (ph, "numchrom:", &numchrom);
  getint (ph, "allowdups:", &allowdups);
  getint (ph, "polarmode:", &polarmode);


  writepars (ph);
  closepars (ph);

}

int
mergeit (SNP ** snpm1, SNP ** snpm2, Indiv *** pindm1, Indiv ** indm2,
         int nums1, int nums2, int numi1, int numi2)
{
  SNP *cupt1, *cupt2;
  int k, x, g, t, tt;
  double y;
  long rlen;
  static unsigned char *packg;
  unsigned char *buff;
  Indiv **indm1;
  static Indiv **indivmarkers;
  int numindivs, numsnps;
  int newnumi1, newnumi2;

  indm1 = *pindm1;
  numindivs = numi1 + numi2;
  numsnps = nums1;
  ZALLOC (indivmarkers, numindivs, Indiv *);

  t = 0;
  for (x = 0; x < numi1; ++x) {
    if (indm1[x]->ignore)
      continue;
    indivmarkers[t] = indm1[x];
    ++t;
  }
  newnumi1 = t;
  for (x = 0; x < numi2; ++x) {
    if (indm2[x]->ignore)
      continue;
    indivmarkers[t] = indm2[x];
    ++t;
  }
  newnumi2 = t - newnumi1;
  numindivs = t;
// we don't bother with a destructor here.   Sloppy code

  y = (double) (numindivs * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  rlen = MAX (rlen, 48);
  packlen = numsnps * rlen;
  ZALLOC (packg, packlen, unsigned char);
  clearepath (packg);
// wipe to invalid

  buff = packg;
  for (k = 0; k < nums1; k++) {
    cupt1 = snpm1[k];
    x = cupt1->tagnumber;
    if (x < 0)
      cupt1->ignore = YES;
    if (cupt1->ignore)
      continue;
    cupt2 = snpm2[x];
    if (cupt2->isrfake) {
      if (phasedmode == NO)
        flipalleles (cupt2);
      if (phasedmode == YES)
        flipalleles_phased (cupt2);
    }
    tt = 0;
    for (t = 0; t < numi1; ++t) {
      if (indm1[t]->ignore)
        continue;
      g = getgtypes (cupt1, t);
      if (g < 0)
        g = 3;
      wbuff ((unsigned char *) buff, tt, g);
      ++tt;
    }
    for (t = 0; t < numi2; ++t) {
      if (indm2[t]->ignore)
        continue;
      g = getgtypes (cupt2, t);
      if (g < 0)
        g = 3;
      wbuff ((unsigned char *) buff, tt, g);
      ++tt;
    }
    cupt1->ngtypes = numindivs;
    cupt1->pbuff = (char *) buff;
    buff += rlen;
  }
  *pindm1 = indivmarkers;
  return numindivs;
}
