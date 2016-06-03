#include <fcntl.h>
#include <ctype.h>
#include  <mcio.h>
#include  <xsearch.h>
#include  <ranmath.h>

/*! \file mcio.c 
 *
 * \brief Input/Output Library
*/



/* global data */

extern int numchrom;
int usecm = NO;			//!< genetic distances are in cMorgans
int plinkinputmode = NO;
static int snprawtab = NO;
static int debug = NO;
extern char *trashdir;
extern int qtmode;		//!< user parameter (phenotype is quantitative)
extern int verbose;		//!< user parameter (print additional output to stdout)
extern int familynames;		//!< user parameter (prepend PLINK family names with colon to individual names)
extern double lp1, lp2;
extern double a1, b1;

extern int packmode;		//!< flag - input {is not,is} in packed mode
extern char *packgenos;		//!< packed genotype data (packit.h)
extern char *packepath;
extern long packlen;		//!< allocated size of packgenos data space
extern long rlen;		//!< number of bytes in packgenos space that each SNP's data occupies
extern int malexhet;		//!< user parameter (retain het genotype data on male X chromosome)
extern int hashcheck;		//!< user parameter (check input file hashes against input data)
extern int outputall;
extern int sevencolumnped;
static int dofreeped = YES;

int tempnum = 0;
int tempfake = 0;

static int *snpord = NULL;	//!< snpord[i] == j if and only if  snpm[j] is ith SNP in input file
static int numsnpord = 0;	//!< current size of array snpord
static int *snporda[3];		//!< Copies of snpord for various data sets (used by mergeit)
static int numsnporda[3];	//!< Number of elements of snporda in use

static int badpedignore = NO;	//!< flag - ignore bad allele symbols in PED file 

static int maxgenolinelength = -1;
static int tersemode = NO;
int checksizemode = YES;
int pedignore = YES;
enum outputmodetype outputmode = PACKEDANCESTRYMAP;
static double maxgpos[MAXCH];
static int chrmode = NO;
static int chimpmode = NO;
static int pordercheck = YES;
static int snpordered;
static int isgdis = YES;	// no means gsid 0 in input 
// fails if packed and out of order 
static int familypopnames = NO;
// in .fam output use popnames (egroup) 


SNPDATA *tsdpt;

/* local function prototypes */

int getbedgenos (char *gname, SNP ** snpmarkers, Indiv ** indivmarkers,
		 int numsnps, int numindivs, int nignore);

void freeped ();

static char x2base (int x);
static void gtox (int g, char *cvals, int *p1, int *p2);

int ancval (int x);
static int setskipit (char *sx);	// ignore lines in snp, map files
int calcishash (SNP ** snpm, Indiv ** indiv, int numsnps, int numind,
		int *pihash, int *pshash);

/* ---------------------------------------------------------------------------------------------------- */

void
setfamilypopnames (int fpop)
{

  familypopnames = fpop;

}


void
clearsnpord ()
{

  free (snpord);
  snpord = NULL;
  numsnpord = 0;

}

void
snpsortit (int **spos, int *indx, int n)
{
  long *lkode;
  int i, base[3];

  base[0] = 1;
  base[1] = 10 ^ 8;
  base[2] = 10 ^ 9;

  ZALLOC (lkode, n, long);
  for (i = 0; i < n; i++) {
    lkode[i] = lkodeitbb (spos[i], 3, base);
  }

  lsortit (lkode, indx, n);

  free (lkode);
  return;

}

int
getsnps (char *snpfname, SNP *** snpmarkpt, double spacing,
	 char *badsnpname, int *numignore, int numrisks)
{
  // returns number of SNPS
  // numrisks 
  /* read file of real SNPS store in temporary structure */

  SNPDATA **snpraw, *sdpt;
  static SNP **snpmarkers;
  SNP *cupt;
  int **snppos;
  int nreal, nfake, numsnps = 0, i, t, j;
  int *snpindx;
  double xspace;
  int failx = 0;

  if (snpfname == NULL)
    fatalx ("(getsnps) null snpname");
  xspace = spacing;
  nreal = getsizex (snpfname);
  if (nreal <= 0)
    fatalx ("no snps found: snpfname: %s\n", snpfname);
  ZALLOC (snpraw, nreal, SNPDATA *);

  if (snpord == NULL) {
    ZALLOC (snpord, nreal, int);
    ivclear (snpord, -1, nreal);
    numsnpord = nreal;
  }
  for (i = 0; i < nreal; i++) {
    ZALLOC (snpraw[i], 1, SNPDATA);
    cclear ((unsigned char *) snpraw[i]->cchrom, CNULL, 7);
    snpraw[i]->inputrow = -1;
    snpraw[i]->alleles[0] = '1';
    snpraw[i]->alleles[1] = '2';
  }
  nreal = readsnpdata (snpraw, snpfname);
  dobadsnps (snpraw, nreal, badsnpname);

  ZALLOC (snppos, nreal, int *);
  for (i = 0; i < nreal; i++) {
    ZALLOC (snppos[i], 3, int);
  }

  for (i = 0; i < nreal; i++) {
    sdpt = snpraw[i];
    snppos[i][0] = sdpt->chrom;
    if ((sdpt->ignore) && (plinkinputmode)) {
      snppos[i][0] = 99;
      if (pordercheck == YES) {
	pordercheck = NO;
	printf ("PLINK input. No check on SNP order\n");
      }
    }
    t = snppos[i][1] = nnint ((sdpt->gpos) * GDISMUL);
    if (isgdis)
      snppos[i][1] = 0;
    snppos[i][2] = nnint (sdpt->ppos);
    // sdpt -> gpos = ((double) t)/ GDISMUL ;
  }

/**
  for (i=nreal-10; i<nreal; i++) { 
   printf("zzyy: %d ", i) ; printimat(snppos[i], 1, 3) ;
  }
*/

  ZALLOC (snpindx, nreal, int);
//snpsortit(snppos, snpindx, nreal) ;
  ipsortit (snppos, snpindx, nreal, 3);

  snpordered = YES;

  for (i = 0; i < nreal; ++i) {
    j = snpindx[i];
    sdpt = snpraw[j];

    //printf("zzz %d %d %s ",   i, j, sdpt -> ID) ;
    //printimat(snppos[j], 1, 3) ;

    if (j != i) {
      snpordered = NO;
      ++failx;
      if (failx < 10) {
	printf
	  ("snp order check fail; snp list not ordered: %s (processing continues)",
	   snpfname);
	printimat (snppos[i], 1, 3);
	printf ("zzz %d %d\n", i, j);
      }
    }
  }

  if ((usecm) && (xspace > 0.5)) {
    printf ("*** warning fake spacing given in cM\n");
    xspace /= 100.0;
  }

  // get number of fakes

  nfake = numfakes (snpraw, snpindx, nreal, xspace);
  numsnps = nreal + nfake;

  tempnum = numsnps;
  tempfake = nfake;

  // allocate storage

  ZALLOC (snpmarkers, numsnps, SNP *);
  for (i = 0; i < numsnps; i++) {
    ZALLOC (snpmarkers[i], 1, SNP);
    cupt = snpmarkers[i];
    clearsnp (cupt);
    ZALLOC (cupt->modelscores, numrisks, double);
    ZALLOC (cupt->totmodelscores, numrisks, double);
  }
  tsdpt = snpraw[0];
  *snpmarkpt = snpmarkers;
  numsnps = loadsnps (snpmarkers, snpraw, snpindx, nreal, xspace, numignore);

/**
  for (i=numsnps-10; i<numsnps; i++) { 
   cupt = snpmarkers[i] ;
   printf("zzyy3: %d %d %12.0f\n", i, cupt -> chrom, cupt -> physpos) ;
  }
*/




  // and free up temporary storage
  for (i = 0; i < nreal; i++) {
    free (snpraw[i]);
    free (snppos[i]);
  }
  free (snpraw);
  free (snppos);
  free (snpindx);

  /* printf("numsnps: %d\n", numsnps) ; */

  /* 
     if (snpord != NULL) { 
     printimat(snpord, 1, MIN(100, numsnps)) ;
     }
   */
  cupt = snpmarkers[0];
  if (isnumword (cupt->ID))
    printf
      ("*** warning: first snp %s is number.  perhaps you are using .map format\n",
       cupt->ID);

  return numsnps;
}




/* ---------------------------------------------------------------------------------------------------- */
int
getsizex (char *fname)
{
  char line[MAXSTR + 1], c;
  char *spt[MAXFF], *sx;
  int nsplit, num = 0;
  int skipit;
  int len;

  FILE *fff;
  openit (fname, &fff, "r");
  line[MAXSTR] = '\0';
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = setskipit (sx);	// comment line
    if (skipit == NO) {
      ++num;
    }

    // now flush the rest of the line if necessary.
    len = strlen (line);
    c = line[len - 1];
    if (c != '\n') {
      while ((c = fgetc (fff)) != EOF) {
	if (c == '\n')
	  break;
      }
    }
    freeup (spt, nsplit);
    continue;
  }
  fclose (fff);
  fflush (stdout);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
ismapfile (char *fname)
{
  // PLINK map file ? 
  // just look at file name (perhaps should look at format) 

  char *sx;
  int len;
  len = strlen (fname);
  if (len < 4)
    return NO;
  sx = fname + len - 4;

  if (strcmp (sx, ".map") == 0)
    return YES;
  if (strcmp (sx, ".bim") == 0)
    return YES;

  if (len < 7)
    return NO;
  sx = fname + len - 7;
  if (strcmp (sx, ".pedsnp") == 0)
    return YES;

  return NO;

}



/* ---------------------------------------------------------------------------------------------------- */
int
ispedfile (char *fname)
{
  // PLINK ped file ? 
  // just look at file name (perhaps should look at format) 
  char *sx;
  int len;
  len = strlen (fname);
  if (len < 4)
    return NO;
  sx = fname + len - 4;

  if (strcmp (sx, ".ped") == 0)
    return YES;
  if (strcmp (sx, ".fam") == 0)
    return YES;

  if (len < 7)
    return NO;
  sx = fname + len - 7;
  if (strcmp (sx, ".pedind") == 0)
    return YES;

  return NO;
}


/* ---------------------------------------------------------------------------------------------------- */
int
isbedfile (char *fname)
{
  // PLINK ped file ? 
  // just look at file name (perhaps should look at format) 

  char *sx;
  int len;
  len = strlen (fname);
  if (len < 4)
    return NO;
  sx = fname + len - 4;

  if (strcmp (sx, ".bed") == 0)
    return YES;
  return NO;

}

/* ---------------------------------------------------------------------------------------------------- */
int
readsnpdata (SNPDATA ** snpraw, char *fname)
{
  char line[LONGSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k;
  int skipit;
  SNPDATA *sdpt;

  double maxg = -9999.0;

  FILE *fff;
  int chrom;
  int nbad = 0;

  plinkinputmode = NO;
  // if this is a PLINK file, call PLINK input routine
  if (ismapfile (fname)) {
    plinkinputmode = YES;
    return readsnpmapdata (snpraw, fname);
  }
  usecm = NO;

  vclear (maxgpos, -9999.0, MAXCH);
  openit (fname, &fff, "r");
  while (fgets (line, LONGSTR, fff) != NULL) {

    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = setskipit (sx);
    if (skipit == NO) {
      if (nsplit < 4)
	fatalx ("(readsnpdata) bad line: %s\n", line);
      sdpt = snpraw[num];
      sdpt->inputrow = num;

      if (strlen (spt[0]) >= IDSIZE)
	fatalx ("ID too long\n", spt[0]);
      strcpy (sdpt->ID, spt[0]);

      sdpt->chrom = chrom = str2chrom (spt[1]);
      strncpy (sdpt->cchrom, spt[1], 6);

      if ((chrom >= MAXCH) || (chrom <= 0)) {
	if (nbad < 10)
	  printf ("warning: bad chrom: %s", line);
	++nbad;

	sdpt->chrom = MIN (chrom, BADCHROM);
	sdpt->chrom = MAX (chrom, 0);
	sdpt->ignore = YES;
      }

      // the genetic positions will be converted to Morgans (assumed to be in cM) if and only if
      // any genetic position is greater than 100

      sdpt->gpos = atof (spt[2]);
      if (sdpt->gpos > 100) {
	if (sdpt->gpos > 1.0e6)
	  fatalx ("absurd genetic distance:\n%s\n", line);
	if (!usecm) {
	  printf ("*** warning.  genetic distances are in cM not Morgans\n");
	  printf ("%s\n", line);
	}
	usecm = YES;		// set flag to connvert to Morgans
      }

      maxgpos[chrom] = MAX (maxgpos[chrom], sdpt->gpos);
      maxg = MAX (maxg, maxgpos[chrom]);

      setsdpos (sdpt, atoi (spt[3]));
      if (nsplit < 8) {
	ivzero (sdpt->nn, 4);
	if (nsplit == 6) {
	  sx = spt[4];
	  sdpt->alleles[0] = toupper (sx[0]);
	  sx = spt[5];
	  sdpt->alleles[1] = toupper (sx[0]);
	}
      }
      else {			// QUESTION:  when does a SNP file have more than seven columns?
	for (k = 0; k < 4; k++) {
	  sdpt->nn[k] = atoi (spt[4 + k]);
	}
	if (nsplit == 10) {
	  sx = spt[8];
	  sdpt->alleles[0] = toupper (sx[0]);
	  sx = spt[9];
	  sdpt->alleles[1] = toupper (sx[0]);
	}
      }
      ++num;
    }
    freeup (spt, nsplit);
    continue;
  }				// elihw

  // if all genetic positions are set to zero, set from physical position 
  if (maxg <= 0.00001) {
    isgdis = NO;
    printf ("%s: genetic distance set from physical distance\n", fname);
    usecm = NO;
    for (k = 0; k < num; ++k) {
      snpraw[k]->gpos = 1.0e-8 * snpraw[k]->ppos;
    }
  }

  // convert to Morgans 
  if (usecm) {
    for (k = 0; k < num; ++k) {
      snpraw[k]->gpos /= 100.0;
    }
  }

  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
readsnpmapdata (SNPDATA ** snpraw, char *fname)
{
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k, t;
  int skipit, len;
  SNPDATA *sdpt;
  int nbad = 0;

  FILE *fff;
  int chrom;
  double maxg = -9999.0;

  vclear (maxgpos, -9999.0, MAXCH);
  openit (fname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = setskipit (sx);
    if (skipit == NO) {
      if (nsplit < 4)
	fatalx ("(readsnpmapdata) bad line: %s\n", line);
      sdpt = snpraw[num];
      if (strlen (spt[1]) >= IDSIZE)
	fatalx ("ID too long\n", spt[1]);
      strcpy (sdpt->ID, spt[1]);

      if (nsplit >= 6) {	// alleles in .map file are optional
	sx = spt[4];
	sdpt->alleles[0] = sx[0];
	sx = spt[5];
	sdpt->alleles[1] = sx[0];
	if (sdpt->alleles[0] == '0')
	  sdpt->alleles[0] = 'X';	// unknown
	if (sdpt->alleles[1] == '0')
	  sdpt->alleles[1] = 'X';
      }
      else {
	cclear ((unsigned char *) sdpt->alleles, CNULL, 2);
      }

      sx = spt[0];
      sdpt->chrom = chrom = str2chrom (sx);
      strncpy (sdpt->cchrom, sx, 6);

      if ((chrom >= MAXCH) || (chrom <= 0)) {
	if (nbad < 10)
	  printf ("warning (mapfile): bad chrom: %s", line);
	++nbad;

	sdpt->chrom = MIN (chrom, BADCHROM);
	sdpt->chrom = MAX (chrom, 0);
	sdpt->chrom = 99;
	strcpy (sdpt->cchrom, "99");
	sdpt->ignore = YES;
      }

      // the genetic positions will be converted to Morgans (assumed to be in cM) if and only if
      // any genetic position is greater than 100

      sdpt->gpos = atof (spt[2]);
      if (sdpt->gpos > 100) {
	if (sdpt->gpos > 1.0e6)
	  fatalx ("absurd genetic distance:\n%s\n", line);
	if (!usecm) {
	  printf ("*** warning.  genetic distances are in cM not Morgans\n");
	  printf ("%s\n", line);
	}
	usecm = YES;
      }
      maxgpos[chrom] = MAX (maxgpos[chrom], sdpt->gpos);
      maxg = MAX (maxg, maxgpos[chrom]);
      sdpt->ppos = atof (spt[3]);
      if (nsplit < 8) {
	ivzero (sdpt->nn, 4);
      }
      else {
	for (k = 0; k < 4; k++) {
	  sdpt->nn[k] = atoi (spt[4 + k]);
	}
      }
      sdpt->inputrow = num;
//    printf("zz %d %d %s %12.0f\n", num, sdpt -> chrom, sdpt -> ID, sdpt -> ppos) ;
      ++num;
    }
    freeup (spt, nsplit);
    continue;
  }

  if (maxg <= 0.00001) {
    printf ("genetic distance set from physical distance\n");
    usecm = NO;
    isgdis = NO;
    for (k = 0; k < num; ++k) {
      snpraw[k]->gpos = 1.0e-8 * snpraw[k]->ppos;
    }
  }

  if (usecm) {
    for (k = 0; k < num; ++k) {
      snpraw[k]->gpos /= 100.0;
    }
  }

  if (snpord == NULL) {
    ZALLOC (snpord, num, int);
    ivclear (snpord, -1, num);
    numsnpord = num;
  }

  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
numfakes (SNPDATA ** snpraw, int *snpindx, int nreal, double spacing)
{

  // it seems better for this internal routine 
  // to use the precomputed values

  int nignore, numsnps;
  int nfake = 0, i, k, indx;
  int num = 0;
  SNP *cupt;
  SNPDATA *sdpt;
  char *sname;
  int *sp;
  int xc = 0, chrom;
  double fakedis, realdis;	// gpos for fake marker 
  double yf, yr;
  double physpos;

  if (spacing <= 0.0)
    fakedis = 1.0e20;

  for (k = 0; k < nreal; k++) {

    indx = snpindx[k];
    sdpt = snpraw[indx];

    chrom = sdpt->chrom;
    realdis = sdpt->gpos;

    if (chrom != xc) {
      fakedis = nextmesh (realdis, spacing);
      xc = chrom;
    }
    while (fakedis < realdis) {
      fakedis += spacing;
      ++nfake;
    }
  }

  // nfake is number of multiples of fakedis in chromosome

  return nfake;
}



/* ---------------------------------------------------------------------------------------------------- */
double
nextmesh (double val, double spacing)
{
  double y;

  if (spacing == 0.0)
    return 1.0e8;
  y = ceil (val / spacing) * spacing;
  if (y < val)
    y += spacing;
  return y;
}


/* ---------------------------------------------------------------------------------------------------- */
/*! \fn int loadsnps(SNP **snpm, SNPDATA **snpraw, 
            int *snpindx, int nreal, double spacing, int *numignore)    
    \brief  Store raw SNP data in final array of type SNP *
    \param  snpm    Pointer to array of type SNP * in which to store data
    \param  snpraw  Pointer to array of type SNPDATA * in which preliminary data was stored
    \param  snpindx  On entry, kth element of snpindx is index of the kth SNP in snpraw (This is not
      the same as the value k itself if the SNPs were out of order in the file.)
    \param  nreal   Number of SNPs stored in snpraw
    \param  spacing    Maximum spacing between SNPs (not relevant to EIGENSOFT)
    \param  numignore  Return number of SNPs to ignore here
*/


int
loadsnps (SNP ** snpm, SNPDATA ** snpraw,
	  int *snpindx, int nreal, double spacing, int *numignore)
{

  // snppos, snpindx could be recalculated but 
  // it seems better for this internal routine 
  // to use the precomputed values
  // do NOT call externally

  int nignore, numsnps;
  int nfake = 0, i, k, indx;
  int num = 0, tnum;
  SNP *cupt = NULL, *lastcupt = NULL, *tcupt;
  SNPDATA *sdpt;
  char *sname;
  int *sp;
  int xc = 0, chrom;
  double fakedis, realdis, xrealdis;	// gpos for fake marker 
  double yf, yr;
  double physpos;
  double xl, xr, xmid, al, ar, fraw;
  double y;
  int nn[2], n0, n1;
  int cnum, t;
  int inputrow, chimpfudge, xchimpfudge;
  int ischimp = NO;
  char ss[6];

  if (spacing <= 0.0)
    fakedis = 1.0e20;
  strcpy (ss, "??");

  for (k = 0; k < nreal; k++) {

    indx = snpindx[k];
    sdpt = snpraw[indx];

    chrom = sdpt->chrom;
// defensive programming;  should not be needed:
    if (sdpt->cchrom[0] == CNULL) {
      sprintf (sdpt->cchrom, "%d", chrom);
    }
    sname = sdpt->ID;
    realdis = sdpt->gpos;
    physpos = sdpt->ppos;
    inputrow = sdpt->inputrow;
    if (sdpt->chimpfudge)
      ischimp = YES;

/**
    if (k>(nreal-10)) { 
     printf("zzyy2b %d %d %12.0f %d\n", k, chrom, physpos, inputrow) ;
    }
*/

    t = strcmp (ss, sdpt->cchrom);
    if (t != 0) {
      fakedis = nextmesh (realdis, spacing);
      xc = chrom;
      cnum = 0;
      strcpy (ss, sdpt->cchrom);
    }

    yf = fakedis;
    yr = realdis;

    // insert fake SNPs so the distance between SNPs is no greater than spacing
    while (fakedis < realdis) {

      if (cnum == 0)
	break;			// first SNP on chromosome
      if (sdpt->ignore)
	break;

      if (nfake >= tempfake)
	fatalx (" too many fake markers (bug) %d %d\n", num, nfake);
      if (num >= tempnum)
	fatalx (" too many markers (bug) %d %d\n", num, nfake);

      cupt = snpm[num];
      if (cupt == NULL)
	fatalx ("bad loadsnps\n");
      sprintf (cupt->ID, "fake-%d:%d", xc, nfake);
      cupt->estgenpos = cupt->genpos = fakedis;
      tcupt = lastcupt;
      for (;;) {
	xl = tcupt->genpos;
	if (xl < fakedis)
	  break;
	tnum = tcupt->markernum;
	--tnum;
	if (tnum < 0)
	  fatalx ("verybadbug\n");
	tcupt = snpm[tnum];
	if (tcupt->chrom != chrom)
	  fatalx ("badbug\n");
      }
      al = tcupt->physpos;
      xr = realdis;;
      ar = physpos;
      y = cupt->physpos = interp (xl, xr, fakedis, al, ar);
      if (chrom == -199) {
	printf
	  ("zzinterp %12.6f  %12.6f  %12.6f  %12.0f  %12.0f     %12.6f\n", xl,
	   xr, fakedis, al, ar, y);
      }
      cupt->markernum = num;
      cupt->isfake = YES;
      cupt->chrom = xc;
      strncpy (cupt->cchrom, ss, 6);
      fakedis += spacing;
      ++num;
      ++nfake;
    }

    cupt = snpm[num];
    if (cupt == NULL)
      fatalx ("bad loadsnps\n");
    strcpy (cupt->ID, sname);
    sdpt->cuptnum = num;
    cupt->estgenpos = cupt->genpos = realdis;
    cupt->physpos = physpos;
    cupt->markernum = num;
    cupt->isfake = NO;
    cupt->ignore = sdpt->ignore;
    // if ((cupt -> ignore == NO) && (cupt -> isfake == NO))   
    if (cupt->isfake = NO) {
      lastcupt = cupt;
      ++cnum;
    }
    cupt->isrfake = sdpt->isrfake;
    cupt->chrom = xc;
    strncpy (cupt->cchrom, ss, 6);
    cupt->tagnumber = inputrow;	// just used for pedfile 
    if (inputrow >= 0) {
      if (inputrow >= numsnpord)
	fatalx ("snpord overflow\n");
      snpord[inputrow] = num;
    }

    n0 = sdpt->nn[0];
    n1 = sdpt->nn[1];
    fraw = mknn (nn, n0, n1);
    copyiarr (nn, cupt->af_nn, 2);
    cupt->aftrue = cupt->af_freq = fraw;
    cupt->aa_aftrue = cupt->aa_af_freq = fraw;

    if (sdpt->alleles != NULL) {
      cupt->alleles[0] = sdpt->alleles[0];
      cupt->alleles[1] = sdpt->alleles[1];
    }
    else {
      cupt->alleles[0] = '1';
      cupt->alleles[1] = '2';
    }

    n0 = sdpt->nn[2];
    n1 = sdpt->nn[3];
    fraw = mknn (nn, n0, n1);
    copyiarr (nn, cupt->cauc_nn, 2);
    cupt->cftrue = cupt->cauc_freq = fraw;
    cupt->aa_cftrue = cupt->aa_cauc_freq = fraw;
    ++num;
  }

  // now make list of ignored snps used by loadgeno for check
  numsnps = num;
  for (k = 0; k < nreal; k++) {
    indx = snpindx[k];
    sdpt = snpraw[indx];
    if (sdpt->ignore == NO)
      continue;
    inputrow = sdpt->inputrow;
    chrom = sdpt->chrom;
    sname = sdpt->ID;
    realdis = sdpt->gpos;
    physpos = sdpt->ppos;
    cupt = snpm[sdpt->cuptnum];
    cupt->tagnumber = inputrow;	// just used for pedfile 
    /* 
       strncpy(cupt -> ID, sname, IDSIZE-1) ; 
       cupt -> genpos = realdis ;                          
       cupt -> physpos = physpos ;
       cupt -> markernum = num ;
       cupt -> isfake = NO ;
       cupt -> ignore = YES ;
       cupt -> chrom  = chrom ;
     */
    ++num;
  }
  nignore = 0;
  for (k = 0; k < numsnps; ++k) {
    cupt = snpm[k];
    if (ischimp && (cupt->chrom == 2))
      cupt->chimpfudge = YES;
    if (cupt->ignore)
      ++nignore;
  }
  *numignore = nignore;
  return numsnps;
}

/* ---------------------------------------------------------------------------------------------------- */
double
interp (double l, double r, double x, double al, double ar)
{
  // linearly interp ;
  double y, y1, y2;
  y = (r - l);
  if (y == 0.0)
    return 0.5 * (al + ar);
  y1 = (r - x) / y;
  y2 = (x - l) / y;
  return y1 * al + y2 * ar;
}

/* ---------------------------------------------------------------------------------------------------- */
int
getindivs (char *indivfname, Indiv *** indmarkpt)
{
  static Indiv **indivmarkers;
  int nindiv, i;

  if (indivfname == NULL)
    fatalx ("(getindivs) NULL indivfname\n");
  nindiv = getsizex (indivfname);
  if (nindiv <= 0)
    fatalx ("no indivs found: indivname: %s\n", indivfname);
  ZALLOC (indivmarkers, nindiv, Indiv *);

  for (i = 0; i < nindiv; i++) {
    ZALLOC (indivmarkers[i], 1, Indiv);
  }
  clearind (indivmarkers, nindiv);
  *indmarkpt = indivmarkers;
  readinddata (indivmarkers, indivfname);
  return nindiv;
}


/* ---------------------------------------------------------------------------------------------------- */
int
readinddata (Indiv ** indivmarkers, char *fname)
{
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k;
  int skipit;
  Indiv *indx;

  FILE *fff;

  // Call routine to read PLINK format file
  if (ispedfile (fname)) {
    plinkinputmode = YES;
    return readindpeddata (indivmarkers, fname);
  }

  // Read ANCESTRYMAP/EIGENSTRAT format individual file
  openit (fname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = setskipit (sx);
    if (skipit == NO) {
      if (nsplit < 3)
	fatalx ("%s bad line: %s", fname, line);
      indx = indivmarkers[num];
      if (strlen (sx) >= IDSIZE)
	fatalx ("ID too long: %s\n", sx);
      strcpy (indx->ID, sx);
      indx->idnum = num;
      sx = spt[1];
      indx->gender = sx[0];
      indx->affstatus = indx->ignore = NO;
      sx = spt[2];
      if (strcmp (sx, "Ignore") == 0)
	indx->ignore = YES;
      if ((qtmode) && (!indx->ignore)) {	// store quantitative phenotype in qval
	indx->egroup = strdup ("Case");
	indx->qval = indx->rawqval = atof (sx);
      }
      else {
	indx->egroup = strdup (sx);	// store discrete phenotype in egroup
      }
      // affstatus set by setstatus  
      ++num;
    }
    freeup (spt, nsplit);
    continue;
  }
  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
readindpeddata (Indiv ** indivmarkers, char *fname)
{
  char *line;
  char *spt[MAXFF], *sx, *sx0, gender;
  int nsplit, num = 0, k, i;
  int skipit;
  Indiv *indx;
  int nindiv;
  int maxnsplit = 0;
  char nnbuff[IDSIZE];
  int nok = 0;

  FILE *fff;

  maxgenolinelength = maxlinelength (fname);
  ZALLOC (line, maxgenolinelength + 1, char);
  openit (fname, &fff, "r");

  while (fgets (line, maxgenolinelength, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx0 = sx = spt[0];
    skipit = NO;
    if (sx[0] == '#')
      skipit = YES;
    if (skipit == NO) {
      if (nsplit < 6)
	fatalx ("%s bad line: %s", fname, line);
      indx = indivmarkers[num];
      if (strlen (sx) >= IDSIZE)
	fatalx ("ID too long: %s\n", sx);
      maxnsplit = MAX (maxnsplit, nsplit);

      sx = spt[1];
      pedname (nnbuff, sx0, sx);
      strcpy (indx->ID, nnbuff);
      indx->idnum = num;
      sx = spt[4];
      k = atoi (sx);
      gender = 'U';
      if (k == 1)
	gender = 'M';
      if (k == 2)
	gender = 'F';
      indx->gender = gender;
      indx->affstatus = indx->ignore = NO;

      sx = spt[5];
      if (qtmode) {
	indx->egroup = strdup ("Case");
	indx->qval = indx->rawqval = atof (sx);
      }
      else {
	k = 99;
	if (strcmp (sx, "-9") == 0)
	  k = -9;
	if (strcmp (sx, "9") == 0)
	  k = 9;
	if (strcmp (sx, "0") == 0)
	  k = 0;
	if ((pedignore == NO) && (k == 0))
	  k = 3;
	if (strcmp (sx, "1") == 0)
	  k = 1;
	if (strcmp (sx, "2") == 0)
	  k = 2;
	switch (k) {
	case 9:
	  indx->ignore = YES;
	  printf ("%s ignored\n", indx->ID);
	  break;
	case -9:
	  indx->ignore = YES;
	  printf ("%s ignored\n", indx->ID);
	  break;
	case 0:
	  indx->ignore = YES;
	  printf ("%s ignored\n", indx->ID);
	  break;
	case 1:
	  indx->egroup = strdup ("Control");
	  break;
	case 2:
	  indx->egroup = strdup ("Case");
	  break;
	case 3:
	  indx->egroup = strdup ("???");
	  break;
	default:
	  indx->egroup = strdup (sx);
	}
      }

      // affstatus set by setstatus  
      if (indx->ignore == NO)
	++nok;
      ++num;
    }
    freeup (spt, nsplit);
    continue;
  }

  if (nok == 0) {
    printf ("all individuals set ignore.  Likely input problem (col 6)\n");
    printf ("resetting all individual...\n");
    for (i = 0; i < num; i++) {
      indx = indivmarkers[i];
      indx->ignore = NO;
      indx->egroup = strdup ("???");
    }
  }

  if (maxnsplit < 8)
    maxgenolinelength = -1;
  free (line);

  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
void
pedname (char *cbuff, char *sx0, char *sx1)
{
  int l0, l1, ll;

  l0 = strlen (sx0);
  l1 = strlen (sx1);
  ll = l0 + l1 + 1;
  if (familynames == NO)
    ll = l1;
  if (ll >= IDSIZE) {
    fatalx ("idnames too long %s %s ll: %d limit: %d\n", sx0, sx1, ll,
	    IDSIZE - 1);
  }
  if (familynames == YES) {	// prepend family name to individual name
    strcpy (cbuff, sx0);
    cbuff[l0] = ':';
    strcpy (cbuff + l0 + 1, sx1);
    return;
  }
  strcpy (cbuff, sx1);

}

/* ---------------------------------------------------------------------------------------------------- */
int
readtldata (Indiv ** indivmarkers, int numindivs, char *inddataname)
{
  // warning printed if theta/lambda not in file
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k, ind, i;
  int skipit;
  Indiv *indx;
  double y;
  double gg[3];
  int *xcheck;

  FILE *fff;
  ZALLOC (xcheck, numindivs, int);
  openit (inddataname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = NO;
    if (strcmp (sx, "Indiv_Index") == 0) {
      // hack.  thetafile should be output with leading ## 
      freeup (spt, nsplit);
      continue;
    }
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < 8)
      fatalx ("%s bad line: %s", inddataname, line);
    sx = spt[1];
    ind = indindex (indivmarkers, numindivs, sx);
    if (ind < 0)
      fatalx ("(readtldata) indiv: %s not found \n", sx);
    indx = indivmarkers[ind];

    indx->theta_mode = atof (spt[3]);
    indx->lambda_mode = atof (spt[7]);
    indx->Xtheta_mode = atof (spt[5]);
    indx->Xlambda_mode = atof (spt[9]);
    xcheck[ind] = 1;

    freeup (spt, nsplit);
    continue;
  }
  for (i = 0; i < numindivs; ++i) {
    indx = indivmarkers[i];
    if (indx->ignore)
      continue;
    if (xcheck[i] == 1)
      continue;
    printf ("*** warning (readtldata) ");
    printf ("%s not found in tlname file", indx->ID);
    printnl ();
  }

  free (xcheck);
  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
readfreqdata (SNP ** snpm, int numsnps, char *inddataname)
{
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k, ind;
  int skipit;
  SNP *cupt;

  FILE *fff;
  openit (inddataname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = NO;
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < 6)
      fatalx ("%s bad line: %s", inddataname, line);
    sx = spt[2];
    ind = snpindex (snpm, numsnps, sx);
    if (ind < 0)
      fatalx ("(readfreqdata) snp %s not found \n", sx);
    cupt = snpm[ind];
    cupt->aa_af_freq = cupt->af_freq = atof (spt[3]);
    cupt->aa_cauc_freq = cupt->cauc_freq = atof (spt[5]);

    freeup (spt, nsplit);
    ++num;
    continue;
  }

  fclose (fff);
  return num;
}



/* ---------------------------------------------------------------------------------------------------- */
int
setstatus (Indiv ** indm, int numindivs, char *smatch)
{
  // return number set 
  // smatch = NULL => set everything 
  return setstatusv (indm, numindivs, smatch, YES);
}

/* ---------------------------------------------------------------------------------------------------- */
int
setstatusv (Indiv ** indm, int numindivs, char *smatch, int val)
{
  // return number set 
  // smatch = NULL => set everything 
  int i, n = 0;
  Indiv *indx;
  char *sx;
  for (i = 0; i < numindivs; i++) {
    indx = indm[i];
    if (indx->ignore)
      continue;
    sx = indx->egroup;
    if (smatch == NULL) {
      ++n;
      indx->affstatus = val;
      continue;
    }
    if (strcmp (sx, smatch) == 0) {
      ++n;
      indx->affstatus += val;
    }
    if (indx->affstatus > 1)
      fatalx ("aff2bug\n");
  }
  return n;
}


/* ---------------------------------------------------------------------------------------------------- */
long
getgenos (char *genoname, SNP ** snpmarkers, Indiv ** indivmarkers,
	  int numsnps, int numindivs, int nignore)
{
  // read genofile.  Use hashtable to improve search 
  // if genofile is gzipped decompress to trashdir
  char *gname, *genotmp = NULL;
  ENTRY *hashlist, *iteml;
  ENTRY item1;
  int k, num, indiv, lgt;
  int val;
  void *basept = 0;
  int bigoff;
  int tcheck;

  // we use a trick: want to store k 
  // store basept + k instead 
  // basept + k + bigoff for individual ID 

  SNP *cupt;
  Indiv *indx;

  char line[MAXSTR], cmd[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, nsnp;
  int skipit, kret, tpackmode, teigenstratmode;

  FILE *fff;
  int gnlen, ngenos = 0;

  double y;
  char *pbuff;

  item1.key = NULL;
  item1.data = NULL;

  if (genoname == NULL)
    fatalx ("(getgenos) NULL genoname\n");
  gname = genoname;
  gnlen = strlen (genoname);

  // Unzip file if necessary
  if (strcmp (genoname + gnlen - 3, ".gz") == 0) {
    makedir (trashdir);
    sprintf (line, "%s/genotmp:%d", trashdir, getpid ());
    genotmp = strdup (line);
    sprintf (cmd, "gunzip -c %s > %s", genoname, genotmp);
    printf ("unzip cmd: %s\n", cmd);
    system (cmd);
    kret = system (cmd);
    if (kret < 0) {
      perror ("gunzip failed\n");
      fatalx ("gunzip failed... probably out of disk space\n");
    }
    printf ("geno file unzipped\n");
    gname = genotmp;
  }

  // Enforce data size limits
  tcheck = checksize (numsnps, numindivs, outputmode);
  if (tcheck == -2)
    fatalx
      ("Data sets with more than 8 billion genotypes are not permitted\n");
  if (tcheck == -1)
    fatalx
      ("Output files of size >2GB are not permitted: use a more compact output data format. Also see documentation of chrom, badsnpname and checksizemode parameters.\n");

  // Call routine to read PLINK format unpacked genotype file
  if (ispedfile (gname)) {

    if (snpord == NULL)
      fatalx ("snpord not allocated (no map file ?)");
    getpedgenos (genoname, snpmarkers, indivmarkers, numsnps, numindivs,
		 nignore);
    freeped ();
    return numsnps * numindivs;
  }

  // Call routine to read PLINK format packed genotype file
  if (isbedfile (gname)) {
    return getbedgenos (genoname, snpmarkers, indivmarkers, numsnps,
			numindivs, nignore);
  }

  // Check whether file is packed ANCESTRYMAP format (packed EIGENSTRAT does not exist)
  tpackmode = ispack (gname);
  nsnp = numsnps;

  // Call routine to read packed ANCESTRYMAP format
  if (tpackmode) {
    inpack (gname, snpmarkers, indivmarkers, nsnp, numindivs);
    for (k = 0; k < nsnp; k++) {
      cupt = snpmarkers[k];
      if (cupt->ignore)
	continue;
      if ((cupt->isfake) && (!(cupt->isrfake)))
	continue;
      if (cupt->gtypes == NULL)
	ZALLOC (cupt->gtypes, 1, int);
      cupt->ngtypes = numindivs;
    }
    packmode = YES;
    return nsnp * numindivs;
  }

  teigenstratmode = iseigenstrat (gname);

  // Call routine to read EIGENSTRAT format
  if (teigenstratmode) {
    packmode = YES;
    ineigenstrat (gname, snpmarkers, indivmarkers, nsnp, numindivs);
    freeped ();
    return nsnp * numindivs;
  }


  // (If execution reaches here, the file is unpacked ANCESTRYMAP format)

  // rlen is number of bytes needed to store each SNP's genotype data
  y = (double) (numindivs * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  packlen = rlen * numsnps;
  if (packlen < 0)
    fatalx ("yuckk\n");
  if (packmode) {
    ZALLOC (packgenos, packlen, char);
    pbuff = packgenos;
    clearepath (packgenos);
  }

  // instantiate hash table
  num = nsnp + numindivs;
  xhcreate (5 * num);
  ZALLOC (hashlist, num, ENTRY);
  bigoff = nsnp + 100;

  // hash SNPs (key=name, value=index in snpmarkers)
  for (k = 0; k < nsnp; k++) {
    cupt = snpmarkers[k];
    if ((cupt->isfake) && (!(cupt->isrfake)))
      continue;
    iteml = hashlist + k;
    iteml->key = cupt->ID;
    iteml->data = basept + k;
    if (xhsearch (*iteml, FIND) != NULL)
      fatalx ("duplicate ID: %s\n", iteml->key);
    (void) xhsearch (*iteml, ENTER);
  }

  // hash individuals (key=name, value=index in indivmarkers)
  for (k = 0; k < numindivs; k++) {

    indx = indivmarkers[k];
    iteml = hashlist + numsnps + k;
    iteml->key = indx->ID;
    iteml->data = basept + k + bigoff;
    if (xhsearch (*iteml, FIND) != NULL)
      fatalx ("duplicate ID: %s\n", iteml->key);
    (void) xhsearch (*iteml, ENTER);
  }

  // read genotype file
  openit (gname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = NO;
    if (sx[0] == '#')
      skipit = YES;
    skipit = setskipit (sx);
    if (skipit == NO) {
      if (nsplit < 3)
	fatalx ("bad geno line.  missing field?\n", line);

      // Look up SNP and individual indices in hash table
      item1.key = spt[0];
      iteml = xhsearch (item1, FIND);
      if (iteml == NULL) {
	fatalx ("(genotypes) bad ID (SNP): %s\n", line);
      }
      k = (int) (iteml->data - basept);

      if (k >= numsnps) {
	fatalx ("bad genotype line: `snp' may be Indiv Id\n%s\n", line);
      }

      cupt = snpmarkers[k];
      if (cupt->ignore) {
	freeup (spt, nsplit);
	continue;
      }
      item1.key = spt[1];
      iteml = xhsearch (item1, FIND);
      if (iteml == NULL) {
	fatalx ("(genotypes) bad ID: (Indiv) %s\n", line);
      }
      indiv = (int) (iteml->data - basept);
      indiv -= bigoff;
      val = atoi (spt[2]);

      indx = indivmarkers[indiv];
      if (indx->ignore)
	val = -1;
      if (checkxval (cupt, indx, val) == NO)
	val = -1;
      if (val > 2) {
	printf ("*** warning invalid genotype: %s %s %d\n",
		cupt->ID, indx->ID, val);
	val = -1;
      }

      if (cupt->ngtypes == 0) {

	// If this is the first datum for this SNP, initialize
	// Set cupt->puff to point to the SNP's data in the genotype array.
	// Set cupt->gtypes to the number of individuals stored in the genotype.

	if (packmode == NO) {
	  ZALLOC (cupt->gtypes, numindivs, int);
	}
	else {
	  ZALLOC (cupt->gtypes, 1, int);
	  cupt->pbuff = pbuff;
	  pbuff += rlen;
	}
	cupt->ngtypes = numindivs;
	for (k = 0; k < numindivs; ++k) {
	  putgtypes (cupt, k, -1);	// initialize all individuals to "missing data"
	}
      }
      putgtypes (cupt, indiv, val);	// store this individual's genotype at this SNP
      ++ngenos;
    }
    freeup (spt, nsplit);
  }
  fclose (fff);

  // destroy hash table
  free (hashlist);
  xhdestroy ();


  // if this is a temporary file (gunzipped), delete it
  if (genotmp != NULL) {
    unlink (gname);
  }
  /* printf("genotype file processed\n") ; */
  freeped ();
  return ngenos;
}


/* ---------------------------------------------------------------------------------------------------- */
void
freeped ()
{
  // destructor for snpord 
  if (snpord == NULL)
    return;
  if (dofreeped == NO)
    return;
  free (snpord);
  snpord = NULL;
  numsnpord = 0;
  maxgenolinelength = -1;
}

/* ---------------------------------------------------------------------------------------------------- */
int
checkxval (SNP * cupt, Indiv * indx, int val)
{
  // check Male X marker not het
  if (cupt->chrom != numchrom + 1)
    return YES;
  if (indx->gender != 'M')
    return YES;
  if (val != 1)
    return YES;
  if (malexhet)
    return YES;
  return NO;
}

/* ---------------------------------------------------------------------------------------------------- */
void
clearsnp (SNP * cupt)
{

  cupt->af_freq = cupt->cauc_freq = -1;
  cupt->aa_af_freq = cupt->aa_cauc_freq = -1;
  cupt->estgenpos = 0;
  cupt->genpos = 0;
  cupt->physpos = 0;
  cupt->ngtypes = 0;
  cupt->pbuff = NULL;
  cupt->ebuff = NULL;
  cupt->gtypes = NULL;
  cupt->modelscores = NULL;
  cupt->totmodelscores = NULL;
  cupt->score = cupt->weight = 0.0;
  cupt->isfake = NO;
  cupt->ignore = NO;
  cupt->isrfake = NO;
  cupt->estdis = 0;
  cupt->dis = 0;
  cupt->esum = 0;
  cupt->lsum = 0;
  cupt->gpsum = 0;
  cupt->gpnum = 0;
  cupt->pcupt = NULL;
  cupt->tagnumber = -1;
  cclear (cupt->cchrom, CNULL, 7);
  strcpy (cupt->cchrom, "");
  cupt->chimpfudge = NO;
  cclear ((unsigned char *) cupt->alleles, CNULL, 2);
}

/* ---------------------------------------------------------------------------------------------------- */
int
rmindivs (SNP ** snpm, int numsnps, Indiv ** indivmarkers, int numindivs)
{
  // squeeze out ignore
  // dangerous bend.  Of course indivmarkers indexing will change
  int n = 0, g, i, k;
  int x;
  Indiv *indx;
  SNP *cupt;

  // n is index of next unused array element 

  for (k = 0; k < numindivs; ++k) {
    if (indivmarkers[k]->ignore == YES)
      continue;			// don't store
    if (n == k) {		// if no ignored found yet, 
      ++n;			//   next unused is next element
      continue;			//   and no need to copy
    }

    // copy k -> n 
    indx = indivmarkers[n];	// if kth element is not ignored, put it 
    indivmarkers[n] = indivmarkers[k];	//   into next unused element
    indx->idnum = n;
    for (i = 0; i < numsnps; i++) {
      cupt = snpm[i];
      if (cupt->gtypes == NULL)
	break;
      if (cupt->ignore)
	continue;		// copy only genotypes of non-ignored SNPs
      g = getgtypes (cupt, k);
      putgtypes (cupt, n, g);
    }
    ++n;
  }

  for (i = 0; i < numsnps; i++) {	// reset number of individuals 
    cupt = snpm[i];
    cupt->ngtypes = n;
  }
  return n;

}

/* ---------------------------------------------------------------------------------------------------- */
int
rmsnps (SNP ** snpm, int numsnps, char *deletesnpoutname)
{

  int i, x;
  SNP *cupt;
  int lastc, chrom;

  freesnpindex ();		// clear hash table

  // wipe out fakes not between real markers
  lastc = -1;
  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (cupt->ignore)
      continue;
    chrom = cupt->chrom;
    if ((cupt->isfake) && (chrom != lastc)) {
      cupt->ignore = YES;	// precedes first real SNP
      logdeletedsnp (cupt->ID, "isfake", deletesnpoutname);
    }
    if (!cupt->isfake)
      lastc = chrom;
  }

  lastc = -1;
  for (i = numsnps - 1; i >= 0; i--) {
    cupt = snpm[i];
    if (cupt->ignore)
      continue;
    chrom = cupt->chrom;
    if ((cupt->isfake) && (chrom != lastc)) {
      cupt->ignore = YES;	// follows last real SNP
      logdeletedsnp (cupt->ID, "isfake", deletesnpoutname);
    }
    if (!cupt->isfake)
      lastc = chrom;
  }

  x = 0;			// index of next retained SNP in the array 
  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (cupt->ignore) {
      freecupt (&cupt);
      continue;
    }
    snpm[x] = snpm[i];
    ++x;
  }

  // reset own-index field
  for (i = 0; i < x; i++) {
    cupt = snpm[i];
    cupt->markernum = i;
  }

  return x;
}


/* ---------------------------------------------------------------------------------------------------- */
void
freecupt (SNP ** cuppt)
{
  SNP *cupt;
  cupt = *cuppt;
  if (cupt->modelscores != NULL) {
    free (cupt->modelscores);
  }
  if (cupt->totmodelscores != NULL) {
    free (cupt->totmodelscores);
  }
  free (cupt);
  cupt = NULL;
}


/* ---------------------------------------------------------------------------------------------------- */
void
clearind (Indiv ** indm, int numind)
{
  Indiv *indx;
  double theta;
  int i;

  for (i = 0; i < numind; i++) {
    indx = indm[i];
    indx->egroup = NULL;
    indx->affstatus = indx->ignore = NO;
    indx->gender = 'U';
    indx = indm[i];
    indx->Xtheta_mode = indx->theta_mode = a1 / (a1 + b1);
    indx->Xlambda_mode = indx->lambda_mode = lp1 / lp2;
    indx->thetatrue = -1.0;	// silly value
    indx->qval = indx->rawqval = 0.0;
  }
  cleartg (indm, numind);
}

/* ---------------------------------------------------------------------------------------------------- */
void
cleartg (Indiv ** indm, int nind)
{
  int i;
  Indiv *indx;

  for (i = 0; i < nind; i++) {
    indx = indm[i];
    vzero (indx->totgamms, 3);
    indx->totscore = 0.0;
  }
}


/* ---------------------------------------------------------------------------------------------------- */
double
mknn (int *nn, int n0, int n1)
{
  double x;
  int t;

  nn[0] = n0 + 1;
  nn[1] = n1 + 1;

  // no clipping.  (Old code clipped here)
  t = intsum (nn, 2);
  x = ((double) nn[0]) / (double) t;

  return x;
}

/* ---------------------------------------------------------------------------------------------------- */
void
setug (Indiv ** indm, int numind, char gender)
{
  Indiv *indx;
  double theta;
  int i;

  for (i = 0; i < numind; i++) {
    indx = indm[i];
    if (indx->gender == 'U')
      indx->gender = gender;
  }
}


/* ---------------------------------------------------------------------------------------------------- */
void
dobadsnps (SNPDATA ** snpraw, int nreal, char *badsnpname)
{

  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXFF];
  char *ss;
  int indx, nsplit, n;

  if (badsnpname == NULL)
    return;
  openit (badsnpname, &fff, "r");

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    if (spt[0][0] == '#') {
      freeup (spt, nsplit);
      continue;
    }

    // look up index in snpraw
    indx = snprawindex (snpraw, nreal, spt[0]);
    if (indx >= 0) {
      snpraw[indx]->ignore = YES;
      if ((nsplit >= 2) && (checkfake (spt[1]))) {
	snpraw[indx]->ignore = NO;
	snpraw[indx]->isrfake = YES;
      }
    }
    freeup (spt, nsplit);
  }
  fclose (fff);
}


/* ---------------------------------------------------------------------------------------------------- */
int
checkfake (char *ss)
{
  // yes if string ss is "Fake" 
  // ss is overwritten

  ss[0] = tolower (ss[0]);
  if (strcmp (ss, "fake") == 0)
    return YES;
  return NO;

}

void
mkchrom (char *ss, int chrom, double *ppos, int fudge, int chrmode)
{
  char *sx;
  int big = 200 * 1000 * 1000;

  sx = ss;
  if (chrmode) {
    strcpy (ss, "chr");
    sx += 3;
  }
  if ((chrom != 2) || (fudge == NO)) {
    sprintf (sx, "%d", chrom);
    return;
  }

  if (*ppos <= big) {
    sprintf (sx, "2a");
  }

  if (*ppos > big) {
    sprintf (sx, "2b");
    *ppos -= big;
  }

}


/* ---------------------------------------------------------------------------------------------------- */
void
printsnps (char *snpoutfilename, SNP ** snpm, int num, Indiv ** indm,
	   int printfake, int printvalids)
{

  int i, chrom;
  double ppos;
  SNP *cupt;
  char ss[10];
  FILE *xfile;
  int numvcase, numvcontrol;
  char c;

  if ((snpoutfilename != NULL) && (strcmp (snpoutfilename, "NULL") == 0))
    return;
  if (snpoutfilename != NULL) {
    openit (snpoutfilename, &xfile, "w");
  }
  else
    xfile = stdout;

  if (tersemode == NO) {
    fprintf (xfile, "\n");
    fprintf (xfile, "###DETAILS ABOUT THE MARKERS\n");
    fprintf (xfile,
	     "##Gen_Pos: genetic position, Phys_pos: Physical position\n");
    fprintf (xfile,
	     "##Afr_vart: Parental African variant allele count, Afr_ref: Parental African reference allele count\n");
    fprintf (xfile,
	     "##Eur_vart: Parental European variant allele count, Eur:ref:Parental European reference allele count\n");

    fprintf (xfile, "\n");
    fprintf (xfile, "%20s %5s %10s %18s", "#SNP_Id", "Chr_Num", "Gen_Pos",
	     "Phys_Pos");
    fprintf (xfile, " %9s %9s %9s %9s", "Afr_vart", "Afr_ref", "Eur_vart",
	     "Eur_ref");
    fprintf (xfile, "\n");
  }
  for (i = 0; i < num; ++i) {
    cupt = snpm[i];
    if (outputall == NO) {
      if (!printfake && (ignoresnp (cupt)))
	continue;
      if (!printfake && (cupt->isrfake))
	continue;
    }

    ppos = cupt->physpos;

    mkchrom (ss, cupt->chrom, &ppos, cupt->chimpfudge, chrmode);
    fprintf (xfile, "%20s %5s ", cupt->ID, ss);

    if (cupt->genpos == 0.0) {
      fprintf (xfile, "%15.0f %15.0f", cupt->genpos, ppos);
    }
    else {
      fprintf (xfile, "%15.6f %15.0f", cupt->genpos, ppos);
    }

    if (tersemode) {
      printalleles (cupt, xfile);
      fprintf (xfile, "\n");
      continue;
    }

    fprintf (xfile, " %8d ", cupt->af_nn[0]);
    fprintf (xfile, "%8d ", cupt->af_nn[1]);
    fprintf (xfile, "%8d ", cupt->cauc_nn[0]);
    fprintf (xfile, "%8d", cupt->cauc_nn[1]);
    if (!printvalids) {
      printalleles (cupt, xfile);
      fprintf (xfile, "\n");
      continue;
    }
    numvcase = numvalidgtx (indm, cupt, 1);
    numvcontrol = numvalidgtx (indm, cupt, 0);
    fprintf (xfile, "   %6d %6d", numvcase, numvcontrol);
    fprintf (xfile, "    %d %d %d", cupt->ignore, cupt->isfake,
	     cupt->isrfake);
    printalleles (cupt, xfile);
    fprintf (xfile, "\n");
  }
  if (snpoutfilename != NULL)
    fclose (xfile);
}


/* ---------------------------------------------------------------------------------------------------- */
void
printalleles (SNP * cupt, FILE * fff)
{
  char c;
  if ((c = cupt->alleles[0]) != CNULL)
    fprintf (fff, " %c", c);
  if ((c = cupt->alleles[1]) != CNULL)
    fprintf (fff, " %c", c);

}

/* ---------------------------------------------------------------------------------------------------- */
void
printdata (char *genooutfilename, char *indoutfilename,
	   SNP ** snpm, Indiv ** indiv, int numsnps, int numind, int packem)
{

  FILE *gfile, *ifile;
  int i, j, t;
  SNP *cupt;
  Indiv *indx;
  char ss[MAXSTR];
  char *gfilename;
  int dogenos = YES;

  if (packem)
    printf ("packedancestrymap output\n");
  else
    printf ("ancestrymap output\n");

  if ((genooutfilename != NULL) && (strcmp (genooutfilename, "NULL") == 0))
    dogenos = NO;
  if (genooutfilename == NULL)
    dogenos = NO;

  if (dogenos) {
    gfilename = genooutfilename;
    if (packem) {
      outpack (genooutfilename, snpm, indiv, numsnps, numind);
      gfilename = NULL;
    }

    // print unpacked genotype output
    if (gfilename != NULL) {
      openit (gfilename, &gfile, "w");
      if (tersemode == NO)
	fprintf (gfile, "#SNP_ID,INDIV_ID,VART_ALLELE_CNT\n");
    }

    for (i = 0; i < numsnps; i++) {
      if (gfilename == NULL)
	break;
      cupt = snpm[i];

      if (outputall == NO) {
	if (ignoresnp (cupt))
	  continue;
	if (cupt->isrfake)
	  continue;
      }

      for (j = 0; j < cupt->ngtypes; j++) {
	indx = indiv[j];
	if (indx->ignore)
	  continue;
	fprintf (gfile, "%20s %20s %3d\n", cupt->ID, indx->ID,
		 getgtypes (cupt, j));
      }
    }

    if (gfilename != NULL)
      fclose (gfile);

  }

  if (indoutfilename == NULL)
    return;
  if ((indoutfilename != NULL) && (strcmp (indoutfilename, "NULL") == 0))
    return;
  if (indoutfilename != NULL)
    openit (indoutfilename, &ifile, "w");

  /* fprintf(ifile,"#INDIV,GENDER,POPULATION\n"); */
  for (i = 0; i < numind; i++) {
    indx = indiv[i];
    if (indx->ignore)
      continue;
    strcpy (ss, indx->egroup);
    if ((qtmode) && (!indx->ignore)) {
      sprintf (ss, "%9.3f", indx->rawqval);
    }
    if (tersemode) {
      fprintf (ifile, "%20s %c %10s", indx->ID, indx->gender, ss);
      fprintf (ifile, "\n");
      continue;
    }
    t = numvalids (indx, snpm, 0, numsnps - 1);
    fprintf (ifile, "%20s %c %10s %5d\n", indx->ID, indx->gender, ss, t);
  }

  if (indoutfilename != NULL)
    fclose (ifile);
}




/* ---------------------------------------------------------------------------------------------------- */
int
readindval (Indiv ** indivmarkers, int numindivs, char *inddataname)
{

  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k, ind;
  int skipit;
  Indiv *indx;
  double y;
  double gg[3];

  FILE *fff;
  openit (inddataname, &fff, "r");
  for (k = 0; k < numindivs; ++k) {
    indx = indivmarkers[k];
    indx->affstatus = NO;
    indx->qval = -999.0;
  }
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit < 2) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (strcmp (sx, "Indiv_Index") == 0) {
      freeup (spt, nsplit);
      continue;
    }
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    ind = indindex (indivmarkers, numindivs, sx);
    if (ind < 0)
      fatalx ("(readindval) indiv: %s not found \n", sx);
    indx = indivmarkers[ind];
    indx->qval = atof (spt[1]);
    indx->affstatus = YES;
    freeup (spt, nsplit);
    continue;
  }

  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
readgdata (Indiv ** indivmarkers, int numindivs, char *gname)
     // only needed for logreg 
     // not correct for X chromosome 
     // Needs correction  for males 
{
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k, ind;
  int skipit;
  Indiv *indx;
  double y;
  double gg[3];

  FILE *fff;

  cleartg (indivmarkers, numindivs);
  openit (gname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = NO;
    skipit = setskipit (sx);
    if (skipit == NO) {
      if (nsplit < 4)
	fatalx ("%s bad line: %s", gname, line);
      ind = indindex (indivmarkers, numindivs, sx);
      if (ind < 0)
	fatalx ("(readgdata) indiv: %s not found \n", sx);
      indx = indivmarkers[ind];
      for (k = 0; k < 3; k++) {
	gg[k] = atof (spt[k + 1]);
      }
      y = asum (gg, 3);
      vst (gg, gg, 1.0 / y, 3);
      y = 0.5 * (gg[1] + 2.0 * gg[2]);	/* est caucasian ancestry */
      indx->thetatrue = y;
      copyarr (gg, indx->totgamms, 3);
    }
    freeup (spt, nsplit);
    continue;
  }

  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
putweights (char *fname, SNP ** snpm, int numsnps)
{
  int num = 0, k;
  SNP *cupt;
  double weight;

  FILE *fff;
  openit (fname, &fff, "w");

  for (k = 0; k < numsnps; ++k) {
    cupt = snpm[k];
    if (cupt->ignore)
      continue;
    fprintf (fff, "%20s ", cupt->ID);
    fprintf (fff, "%15.9f ", cupt->weight);
    fprintf (fff, "\n");
    ++num;
  }
  fclose (fff);
  return num;
}


/* ---------------------------------------------------------------------------------------------------- */
int
getweights (char *fname, SNP ** snpm, int numsnps)
{
  // number of real lines 
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0;
  int skipit, k;
  double weight;

  FILE *fff;
  for (k = 0; k < numsnps; ++k) {
    snpm[k]->weight = 1.0;
  }
  openit (fname, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0) {
      continue;
    }
    sx = spt[0];
    skipit = NO;
    skipit = setskipit (sx);
    k = snpindex (snpm, numsnps, sx);
    if (k < 0)
      skipit = YES;
    if (skipit == NO) {
      if (nsplit > 1) {
	sx = spt[1];
	weight = atof (sx);
	snpm[k]->weight = weight;
	printf ("weight set: %20s %9.3f\n", snpm[k]->ID, weight);
	++num;
      }
    }
    freeup (spt, nsplit);
    continue;
  }
  fclose (fff);
  fflush (stdout);
  return num;
}

/* ---------------------------------------------------------------------------------------------------- */
void
outpack (char *genooutfilename, SNP ** snpm, Indiv ** indiv, int numsnps,
	 int numind)
{
  char **arrx;
  int n, num, ihash, shash, i, g, j, k;
  int nind, nsnp, irec;
  Indiv *indx;
  SNP *cupt;
  double y;
  unsigned char *buff;
  int fdes, ret;
  char *packit;

  n = numind;
  ZALLOC (arrx, n, char *);

  num = 0;
  for (i = 0; i < n; i++) {
    indx = indiv[i];
    if ((outputall == NO) && indx->ignore)
      continue;
    arrx[num] = strdup (indx->ID);
    ++num;
  }

  // compute hash on individuals
  ihash = hasharr (arrx, num);
  nind = num;
  freeup (arrx, num);
  free (arrx);

  n = numsnps;
  ZALLOC (arrx, n, char *);
  num = 0;
  for (i = 0; i < n; i++) {
    cupt = snpm[i];
    if (outputall == NO) {
      if (ignoresnp (cupt))
	continue;
      if (cupt->isrfake)
	continue;
    }
    arrx[num] = strdup (cupt->ID);
    ++num;
  }

  // compute hash on SNPs
  shash = hasharr (arrx, num);
  nsnp = num;
  freeup (arrx, num);
  free (arrx);
  // printf("ihash:  %x   shash: %x\n", ihash, shash) ; 

  // rlen is number of bytes each SNP will occupy in packed format
  y = (double) (nind * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  rlen = MAX (rlen, 48);
  // printf("nind: %d  rlen:  %d\n", nind, rlen) ;
  ZALLOC (buff, rlen, unsigned char);
  sprintf ((char *) buff, "GENO %7d %7d %x %x", nind, nsnp, ihash, shash);

  ridfile (genooutfilename);
  fdes = open (genooutfilename, O_CREAT | O_TRUNC | O_RDWR, 0666);

  if (fdes < 0) {
    perror ("bad genoout");
    fatalx ("open failed for %s\n", genooutfilename);
  }
  if (verbose)
    printf ("file %s opened\n", genooutfilename);

  ret = write (fdes, buff, rlen);
  if (ret < 0) {
    perror ("write failure");
    fatalx ("(outpack) bad write");
  }

  irec = 1;
  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (outputall == NO) {
      if (ignoresnp (cupt))
	continue;
      if (cupt->isrfake)
	continue;
    }
    cclear ((unsigned char *) buff, 0X00, rlen);
    num = 0;
    for (j = 0; j < numind; j++) {
      indx = indiv[j];
      if (indx->ignore)
	continue;
      g = getgtypes (cupt, j);
      if (g < 0)
	g = 3;
      wbuff (buff, num, g);	// store two-bit genotype in packed data buffer
      ++num;
    }
    ret = write (fdes, buff, rlen);	// print out all SNPs in packed data buffer
    if (ret < 0) {
      perror ("write failure");
      fatalx ("(outpack) bad write");
    }
    if (verbose) {
      printf ("record: %4d ", irec);
      for (k = 0; k < rlen; ++k) {
	printf (" %02x", (unsigned char) buff[k]);
      }
      printf ("\n");
    }
    ++irec;
  }
  close (fdes);
  free (buff);
  // printf("check: %s %d\n", genooutfilename, ispack(genooutfilename)) ;
}



/* ---------------------------------------------------------------------------------------------------- */
int
ispack (char *gname)
{
  // checks if file is packed gfile 
  int fdes, t, ret;
  char buff[8];

  fdes = open (gname, O_RDONLY);
  if (fdes < 0) {
    perror ("open failure");
    fatalx ("(ispack) bad open %s\n", gname);
  }
  t = read (fdes, buff, 8);
  if (t < 0) {
    perror ("read failure");
    fatalx ("(ispack) bad read");
  }
  close (fdes);
  buff[4] = '\0';
  ret = strcmp (buff, "GENO");
  if (ret == 0)
    return YES;
  return NO;

}



/* ---------------------------------------------------------------------------------------------------- */
int
iseigenstrat (char *gname)
{

  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k;


  openit (gname, &fff, "r");

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    freeup (spt, nsplit);
    fclose (fff);
    if (nsplit > 1)
      return NO;
    return YES;
  }
  fatalx ("(iseigenstrat) no genotyped data found\n");

}



/* ---------------------------------------------------------------------------------------------------- */
int
ineigenstrat (char *gname, SNP ** snpm, Indiv ** indiv, int numsnps,
	      int numind)
{
  // supports enhanced format fist character X => all missing data for SNP
  FILE *fff;
  char *line = NULL, c;
  char *spt[2], *sx;
  int nsplit, rownum = 0, k, num;
  int maxstr, maxff = 2;
  int nind, nsnp, len;
  double y;
  unsigned char *buff;
  char *packit, *pbuff;
  int g, g1, g2;
  SNP *cupt;
  Indiv *indx;
  int nbad = 0;


  packmode = YES;
  maxstr = numind + 10;
  ZALLOC (line, maxstr, char);

  nind = numind;
  nsnp = numsnps;

  // rlen is number of bytes used to store each SNP's genotype data
  y = (double) (nind * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  rlen = MAX (rlen, 48);
  ZALLOC (buff, rlen, unsigned char);

  packlen = rlen * nsnp;
  if (packgenos == NULL) {
    ZALLOC (packgenos, packlen, char);
    clearepath (packgenos);
  }

  openit (gname, &fff, "r");

  rownum = 0;
  pbuff = packgenos;
  while (fgets (line, maxstr, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0)
      continue;

    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }

    if (nsplit > 1)
      fatalx ("(ineigenstrat) more than 1 field\n");	// white space not expected

    if (rownum >= numsnps)
      fatalx ("(ineigenstrat) too many lines in file %d %d\n", rownum,
	      numsnps);
    num = snpord[rownum];
    cupt = snpm[num];
    ++rownum;
    if (cupt == NULL)
      continue;

    if (cupt->ngtypes == 0) {
      if (packmode == NO) {
	ZALLOC (cupt->gtypes, numind, int);
	ivclear (cupt->gtypes, -1, numind);
      }
      else {
	ZALLOC (cupt->gtypes, 1, int);
	cupt->pbuff = pbuff;
	pbuff += rlen;
      }
      cupt->ngtypes = numind;
    }

    if (sx[0] == 'X') {
      freeup (spt, nsplit);
      continue;
    }

    len = strlen (sx);
    if (len != nind) {
      printf ("(ineigenstrat) bad line %d ::%s\n", rownum, line);
      fatalx ("(ineigenstrat) mismatch line length %d %d\n", len, nind);
    }

    for (k = 0; k < len; k++) {
      sscanf (sx + k, "%c", &c);
      g = -2;
      if (c == '0')
	g = 0;
      if (c == '1')
	g = 1;
      if (c == '2')
	g = 2;
      if (c == '9')
	g = -1;

      if (g == -2)
	fatalx ("(ineigenstrat) bad character %c\n", c);
      if (indiv[k]->ignore)
	g = -1;
      if (checkxval (cupt, indiv[k], g) == NO)
	g = -1;

      indx = indiv[k];
      if (checkxval (cupt, indx, g) == NO)
	g = -1;
      g2 = g;
      if (g2 < 0)
	continue;
      g1 = getgtypes (cupt, k);
      if ((g1 >= 0) && (g1 != g2))
	++nbad;			// something is already stored there
      putgtypes (cupt, k, g2);
    }
    freeup (spt, nsplit);
  }
  if (rownum != numsnps)
    fatalx ("(ineigenstrat) mismatch in numsnps %d and numlines %d\n",
	    numsnps, rownum);
  fclose (fff);
  freestring (&line);

  return nbad;
}

/* ---------------------------------------------------------------------------------------------------- */
int
calcishash (SNP ** snpm, Indiv ** indiv, int numsnps, int numind, int *pihash,
	    int *pshash)
{
  char **arrx;
  int ihash, shash, n, num;
  int i;
  Indiv *indx;
  SNP *cupt;

  n = numind;
  ZALLOC (arrx, n, char *);

  num = 0;
  for (i = 0; i < n; i++) {
    indx = indiv[i];
    arrx[num] = strdup (indx->ID);
    ++num;
  }
  *pihash = hasharr (arrx, num);

  freeup (arrx, num);
  free (arrx);

  n = numsnps;
  ZALLOC (arrx, n, char *);
  num = 0;
  for (i = 0; i < n; i++) {
    cupt = snpm[i];
    if (cupt->isfake)
      continue;
    arrx[num] = strdup (cupt->ID);
    cupt->ngtypes = numind;
    ++num;
  }
  *pshash = hasharr (arrx, num);
  freeup (arrx, num);
  free (arrx);
  return num;

}


long
bigread (int fdes, char *packg, long numbytes)
{
  long x;
  int xx;
  char *pt;
  long nb, t, nr = 0;
  int pswitch = NO;


  pt = packg;

  x = nnint (pow (2, 30));

  nb = numbytes;
  if (nb > x)
    pswitch = YES;

  for (;;) {
    xx = MIN (x, nb);
    t = read (fdes, pt, xx);
    if (t < 0) {
      perror ("read failure");
      fatalx ("(bigread) bad data read");
    }
    if (t != xx) {
      perror ("read failure (length mismatch)");
      fatalx ("(bigread) bad data read (length mismatch) %ld %ld\n", t, xx);
    }
    nb -= xx;
    nr += xx;
    if (pswitch)
      printf ("read %ld bytes\n", nr);
    if (nb == 0)
      break;
    pt += xx;
  }
  return nr;
}

int
getsnpordered ()
{
  return snpordered;
}

void
putsnpordered (int mode)
{
  snpordered = mode;
}

void
setpordercheck (int mode)
{
  pordercheck = mode;
}

void
failorder ()
{
  fatalx
    ("snps out of order and packed format.  Run convertf with pordercheck: NO\n");
}

/* ---------------------------------------------------------------------------------------------------- */
void
inpack (char *gname, SNP ** snpm, Indiv ** indiv, int numsnps, int numind)
{

  char **arrx, junk[10];
  int n, num, ihash, shash, i, g, j, k;
  long t;
  int xihash, xshash, xnsnp, xnind;
  int nind, nsnp, irec;
  Indiv *indx;
  SNP *cupt;
  double y;
  unsigned char *buff;
  int fdes, ret;
  char *packit, *pbuff;

  nind = n = numind;
  nsnp = calcishash (snpm, indiv, numsnps, numind, &ihash, &shash);

  // rlen is the number of bytes needed to store one SNP's genotype data
  y = (double) (nind * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  rlen = MAX (rlen, 48);
  ZALLOC (buff, rlen, unsigned char);

  // open binary file and check readability
  fdes = open (gname, O_RDONLY);
  if (fdes < 0) {
    perror ("open failure");
    fatalx ("(ispack) bad open %s\n", gname);
  }
  t = read (fdes, buff, rlen);
  if (t < 0) {
    perror ("read failure");
    fatalx ("(inpack) bad read");
  }

  if (pordercheck && (snpordered == NO))
    failorder ();

  // check for file modification
  if (hashcheck) {
    sscanf ((char *) buff, "GENO %d %d %x %x", &xnind, &xnsnp, &xihash,
	    &xshash);
    if (xnind != nind)
      fatalx ("OOPS number of individuals %d != %d in input files\n", nind,
	      xnind);
    if (xnsnp != nsnp)
      fatalx ("OOPS number of SNPs %d != %d in input file: %s\n", nsnp, xnsnp,
	      gname);
    if (xihash != ihash)
      fatalx
	("OOPS indiv file has changed since genotype file was created\n");
    if (xshash != shash)
      fatalx ("OOPS snp file has changed since genotype file was created\n");
  }

  packlen = rlen * nsnp;
  ZALLOC (packgenos, packlen, char);
  clearepath (packgenos);

  // printf("packgenos: %x end: %x  len:  %d\n", packgenos, packgenos+packlen-1, packlen) ;

  t = bigread (fdes, packgenos, packlen);
  if (t < 0) {
    perror ("read failure");
    fatalx ("(inpack) bad data read");
  }
  if (t != packlen) {
    perror ("read failure (length mismatch)");
    printf ("numsnps: %d  nsnp (from geno file): %d\n", numsnps, nsnp);
    fatalx ("(inpack) bad data read (length mismatch) %ld %ld\n", t, packlen);
  }
  else
    printf ("packed geno read OK\n");

  // now set up pointers into packed data 
  pbuff = packgenos;
  for (i = 0; i < numsnps; i++) {
    j = snpord[i];
    if (snpordered == YES)
      j = i;
    if (j < 0)
      fatalx ("(inpack) bug\n");
    if (j > nsnp)
      fatalx ("(inpack) bug\n");
    cupt = snpm[j];
/**
    if ((i % 100000) == 0)  { 
     printf("zz %d %d %d\n", i, j, numsnps) ;  fflush(stdout) ;
    }
*/
    if (cupt->isfake)
      continue;
    cupt->pbuff = pbuff;
    pbuff += rlen;
    // now check xhets
    for (k = 0; k < numind; ++k) {
      indx = indiv[k];
      g = getgtypes (cupt, k);
      if (checkxval (cupt, indx, g) == NO) {
	putgtypes (cupt, k, -1);
      }
    }
  }

  free (buff);
  close (fdes);

  printf ("end of inpack\n");
  fflush (stdout);
}


/* ---------------------------------------------------------------------------------------------------- */
void
getsnpsc (char *snpscname, SNP ** snpm, int numsnps)
{

  FILE *fff;
  int score;
  SNP *cupt;
  char line[MAXSTR];
  char *spt[MAXFF], *sx;
  int nsplit, num = 0, k;
  double y;


  if (snpscname == NULL)
    fatalx ("no snpsc file\n");
  else
    openit (snpscname, &fff, "r");

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    k = snpindex (snpm, numsnps, sx);
    if (k < 0) {
      printf ("*** warning.  snp %s in snpscname but not in main snp file\n",
	      spt[0]);
      freeup (spt, nsplit);
      continue;
    }
    y = atof (spt[1]);
    y += .1 * gauss ();		// dither 
    cupt = snpm[k];
    cupt->score = y;
    freeup (spt, nsplit);
  }

  if (snpscname != NULL)
    fclose (fff);
}

/* ---------------------------------------------------------------------------------------------------- */
void
setepath (SNP ** snpm, int nsnps)
{

  int i;
  SNP *cupt;
  char *pbuff;

  if (packlen == 0)
    fatalx ("(setepath) packlen unset\n");
  ZALLOC (packepath, packlen, char);
  printf ("setepath. packlen: %ld  rlen: %ld\n", packlen, rlen);
  pbuff = packepath;
  for (i = 0; i < nsnps; i++) {
    cupt = snpm[i];
    if (cupt->isfake)
      continue;
    cupt->ebuff = pbuff;
    pbuff += rlen;
  }
  clearepath (packepath);
}


/* ---------------------------------------------------------------------------------------------------- */
void
clearepath (char *packp)
{
  cclear ((unsigned char *) packp, 0XFF, packlen);
}


/* ---------------------------------------------------------------------------------------------------- */
int
getpedgenos (char *gname, SNP ** snpmarkers, Indiv ** indivmarkers,
	     int numsnps, int numindivs, int nignore)
{
  int val;
  int ngenos = 0;

  SNP *cupt;
  Indiv *indx;

  char *line;
  char **spt, *sx;
  char c;
  int nsplit, num = 0;
  int skipit;
  int numf, snpnumber, nsnp;
  int k, n, t, i;
  FILE *fff;
  int **gcounts, *gvar, *gref;
  int xvar, xref;
  int parity, colbase, ncols;
  int snpnum;
  int markernum = -99;
  int n1, n2;

  /* 
     markernum = snpindex(snpmarkers, numsnps, "rs3002685") ;
     if (markernum <0) fatalx("qq1") ;
   */

  maxgenolinelength = MAX (maxgenolinelength, maxlinelength (gname));

  // printf("maxlinelen %d\n", maxlinelength(gname)) ;
  ZALLOC (line, maxgenolinelength + 1, char);

  cleargdata (snpmarkers, numsnps, numindivs);
  nsnp = numsnps;

  ZALLOC (gcounts, nsnp, int *);
  for (i = 0; i < nsnp; i++) {
    ZALLOC (gcounts[i], 5, int);
  }
  genopedcnt (gname, gcounts, nsnp);

  ZALLOC (gvar, nsnp, int);
  ZALLOC (gref, nsnp, int);

  // designate ref and var alleles from counts
  setgref (gcounts, nsnp, gvar, gref);

  // Override improvised ref and var designations if they were in the .map file
  for (i = 0; i < nsnp; ++i) {
    cupt = snpmarkers[i];
    if (cupt->alleles[0] != CNULL) {
      c = cupt->alleles[0];
      gvar[i] = xpedval (c);
      c = cupt->alleles[1];
      gref[i] = xpedval (c);
    }
    else {
      c = x2base (gvar[i]);
      cupt->alleles[0] = c;
      c = x2base (gref[i]);
      cupt->alleles[1] = c;
    }
  }

  numf = 2 * nsnp + 10;
  ZALLOC (spt, numf, char *);

  // Read genotype file, one line per individual
  openit (gname, &fff, "r");
  while (fgets (line, maxgenolinelength, fff) != NULL) {

    nsplit = splitup (line, spt, numf);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    skipit = NO;
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }

    // On first individual, set column base (6 or 7)
    if (num == 0) {
      parity = nsplit % 2;
      ncols = nsplit;
      colbase = 6 + parity;
    }
    if (nsplit != ncols)
      fatalx ("bad number of fields %d %d\n", ncols, nsplit);

    // Loop over SNPs 
    for (k = colbase; k < nsplit - 1; k += 2) {
      snpnumber = (k - colbase) / 2;

      if (snpnumber >= numsnpord)
	fatalx ("snpord overflow\n");
      snpnum = snpord[snpnumber];
      if (snpnum < 0)
	fatalx ("logic bug (bad snpord)\n");

      xvar = gvar[snpnum];
      xref = gref[snpnum];

      t = 0;

      n1 = n = pedval (spt[k]);
      n2 = pedval (spt[k + 1]);

      if ((n1 == 5) || (n2 == 5)) {	// Missing data or invalid
	val = -1;
	putgtypes (cupt, num, val);
	continue;
      }

      if ((n < 0) || (n > 4))
	fatalx ("(getpedgenos) %s bad geno %s\n", gname, spt[k]);
      if (n == xvar)
	++t;
      if ((n != xvar) && (n != xref))
	t = -10;

      n = n2;
      if ((n < 0) || (n > 4))
	fatalx ("(getpedgenos) %s bad geno %s\n", gname, spt[k + 1]);
      if (n == xvar)
	++t;
      if ((n != xvar) && (n != xref))
	t = -10;

      if (t < 0)
	t = -1;			// Any unexpected allele is stored as "missing"
      cupt = snpmarkers[snpnum];
      if (cupt->ignore)
	continue;
      val = t;
      if (checkxval (cupt, indivmarkers[num], val) == NO)
	val = -1;
      putgtypes (cupt, num, val);	// Store genotype 
      if (val >= 0)
	++ngenos;

    }				// rof (SNP)
    freeup (spt, nsplit);
    ++num;

  }				// elihw (individual)

  free (spt);
  fclose (fff);

  for (i = 0; i < nsnp; i++) {
    free (gcounts[i]);
  }

  free (gcounts);
  free (gref);
  free (gvar);
  free (line);

  printf ("genotype file processed\n");
  return ngenos;
}


/* ---------------------------------------------------------------------------------------------------- */
void
genopedcnt (char *gname, int **gcounts, int nsnp)
{
  char *line;
  char **spt, *sx;
  int nsplit, num = 0;
  int skipit;
  int numf, snpnumber, snpnum;
  int k, n;
  FILE *fff;
  int parity, ncols, colbase;

  // gcounts already zeroed 

  maxgenolinelength = MAX (maxgenolinelength, maxlinelength (gname));

  //  printf("maxlinelen %d\n", maxlinelength(gname)) ;
  ZALLOC (line, maxgenolinelength + 1, char);

  numf = 2 * nsnp + 10;
  ZALLOC (spt, numf, char *);

  openit (gname, &fff, "r");
  while (fgets (line, maxgenolinelength, fff) != NULL) {

    nsplit = splitup (line, spt, numf);
    if (nsplit == 0)
      continue;
    skipit = NO;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (num == 0) {
      parity = nsplit % 2;
      ncols = nsplit;
      colbase = 6 + parity;	// QUESTION:  what is the optional seventh column?
    }

    for (k = colbase; k < nsplit - 1; k += 2) {
      snpnumber = (k - colbase) / 2;
      if (snpnumber >= numsnpord)
	fatalx ("snpord overflow\n");
      snpnum = snpord[snpnumber];
      if (snpnum < 0)
	fatalx ("logic bug (bad snpord)\n");
      n = pedval (spt[k]);
      //  if ((n<0) || (n>4)) fatalx("(genopedcnt) %s bad geno %s\n", gname, spt[k]) ;
      if ((n < 0) || (n > 4))
	continue;
      if (n > 0) {
	++gcounts[snpnum][n];
	++num;
      }
      n = pedval (spt[k + 1]);
      //  if ((n<0) || (n>4)) fatalx("(genopedcnt) %s bad geno %s\n", gname, spt[k+1]) ;
      if ((n < 0) || (n > 4))
	continue;
      if (n > 0) {
	++gcounts[snpnum][n];
	++num;
      }
    }
    freeup (spt, nsplit);
    continue;
  }
  free (spt);
  free (line);
  fclose (fff);
  return;
}


/* ---------------------------------------------------------------------------------------------------- */
void
outfiles (char *snpname, char *indname, char *gname, SNP ** snpm,
	  Indiv ** indiv, int numsnps, int numindx, int packem, int ogmode)
{
  /*  call at end of main program usually  */

  int sizelimit = 10000000;
  int numind;

  // Squeeze out individuals with ignore flag set
  numind = rmindivs (snpm, numsnps, indiv, numindx);
  if (snpname == NULL) {
    printf ("*** warning output snpname NULL\n");
    printf ("snpname: %s %d\n", snpname, numsnps);
    printf ("indname:  %s %d\n", indname, numind);
    printf ("gname: %s\n", gname);
  }

  switch (outputmode) {

  case EIGENSTRAT:
    printf ("eigenstrat output\n");
    outeigenstrat (snpname, indname, gname, snpm, indiv, numsnps, numind);
    return;

  case PED:
    printf ("ped output\n");
    outped (snpname, indname, gname, snpm, indiv, numsnps, numind, ogmode);
    return;

  case PACKEDPED:
    printf ("packedped output\n");
    outpackped (snpname, indname, gname, snpm, indiv, numsnps, numind,
		ogmode);
    return;

  case PACKEDANCESTRYMAP:
    if (snpname != NULL)
      printsnps (snpname, snpm, numsnps, indiv, NO, NO);
    packem = YES;
    printdata (gname, indname, snpm, indiv, numsnps, numind, packem);
    return;

  case ANCESTRYMAP:
  default:
    if (snpname != NULL)
      printsnps (snpname, snpm, numsnps, indiv, NO, NO);
    packem = NO;
    if (numsnps > (sizelimit / numind))
      packem = YES;
    printdata (gname, indname, snpm, indiv, numsnps, numind, packem);
    return;
  }
}

/* ---------------------------------------------------------------------------------------------------- */
void
outeigenstrat (char *snpname, char *indname, char *gname, SNP ** snpm,
	       Indiv ** indiv, int numsnps, int numind)
{

  FILE *fff, *ifile;
  int g, i, k;
  SNP *cupt;
  Indiv *indx;
  char ss[MAXSTR];


  settersemode (YES);
  if (snpname != NULL)
    printsnps (snpname, snpm, numsnps, indiv, NO, NO);

  // Print individual data to .ind file
  if (indname != NULL) {
    openit (indname, &ifile, "w");
    for (i = 0; i < numind; i++) {
      indx = indiv[i];
      if (indx->ignore)
	continue;
      strcpy (ss, indx->egroup);
      if (qtmode) {
	sprintf (ss, "%9.3f", indx->rawqval);
      }
      fprintf (ifile, "%20s %c %10s", indx->ID, indx->gender, ss);
      fprintf (ifile, "\n");
      continue;
    }
    fclose (ifile);
  }

  if (gname == NULL)
    return;

  // Print genotypes to .geno file
  openit (gname, &fff, "w");
  for (k = 0; k < numsnps; k++) {
    cupt = snpm[k];
    if (outputall == NO) {
      if (ignoresnp (cupt))
	continue;
      if (cupt->isrfake)
	continue;
    }
    for (i = 0; i < numind; i++) {
      indx = indiv[i];
      if (indx->ignore)
	continue;
      g = getgtypes (cupt, i);
      if (g < 0)
	g = 9;
      fprintf (fff, "%1d", g);
    }
    fprintf (fff, "\n");
  }
  fclose (fff);
}



/* ---------------------------------------------------------------------------------------------------- */
void
setgref (int **gcounts, int nsnp, int *gvar, int *gref)
{
  int tt[5];
  int i, kmax;

  for (i = 0; i < nsnp; i++) {
    copyiarr (gcounts[i], tt, 5);
    tt[0] = -9999;		// Ensure "missing data" is not ref or var allele
    ivlmaxmin (tt, 5, &kmax, NULL);
    gvar[i] = kmax;		// designate major allele "variant"
    if (tt[kmax] == 0)
      gvar[i] = 5;
    tt[kmax] = -9999;
    ivlmaxmin (tt, 5, &kmax, NULL);
    gref[i] = kmax;		// designate minor allele "variant"
    if (tt[kmax] == 0)
      gref[i] = 5;
  }
}


/* ---------------------------------------------------------------------------------------------------- */
void
cleargdata (SNP ** snpmarkers, int numsnps, int numindivs)
{

  // wipe out all genotype data 
  int i, k;
  SNP *cupt;
  char *pbuff;
  double y;

  // rlen is number of bytes needed to store each SNP's genotype data in packed mode
  y = (double) (numindivs * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  rlen = MAX (rlen, 48);
  packlen = rlen * numsnps;

  if (packlen <= 0)
    fatalx ("bad packlen\n");

  if ((packmode) && (packgenos == NULL)) {
    ZALLOC (packgenos, packlen, char);
    clearepath (packgenos);
  }

  pbuff = packgenos;

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];
    //  if (cupt -> ignore) continue ;
    if (cupt->ngtypes == 0) {
      if (packmode == NO) {
	ZALLOC (cupt->gtypes, numindivs, int);
      }
      else {
	ZALLOC (cupt->gtypes, 1, int);
	cupt->pbuff = pbuff;
	pbuff += rlen;
      }
      cupt->ngtypes = numindivs;
      for (k = 0; k < numindivs; ++k) {
	putgtypes (cupt, k, -1);
      }
    }
  }
}


/* ---------------------------------------------------------------------------------------------------- */
void
setgenotypename (char **gname, char *iname)
{
  if (ispedfile (iname) == NO)
    return;
  if ((*gname != NULL) && strcmp (*gname, "NULL") == 0) {
    *gname = NULL;
    return;
  }
  if (*gname != NULL)
    return;
  *gname = strdup (iname);
}



/* ---------------------------------------------------------------------------------------------------- */
int
maxlinelength (char *fname)
{
  // linelength including \n 

  int len, maxlen;
  int nl, t;
  FILE *fff;

  maxlen = -1;

  len = 0;
  nl = (int) (unsigned char) '\n';

  openit (fname, &fff, "r");
  while ((t = fgetc (fff)) != EOF) {
    ++len;
    if (t == nl) {
      maxlen = MAX (maxlen, len);
      len = 0;
    }
  }
  return maxlen;
}


/* ---------------------------------------------------------------------------------------------------- */
void
settersemode (int mode)
{
  tersemode = mode;
}


/* ---------------------------------------------------------------------------------------------------- */
void
outindped (char *indname, Indiv ** indiv, int numind, int ogmode)
{

  FILE *fff, *ifile;
  int g, i, k;
  Indiv *indx;
  char c;
  int pgender, astatus;
  int dcode = 1;

  if (indname == NULL)
    return;

  openit (indname, &ifile, "w");
  for (i = 0; i < numind; i++) {
    indx = indiv[i];
    if (indx->ignore)
      continue;

    if (familypopnames != YES)
      fprintf (ifile, "%6d ", i + 1);
    if (familypopnames == YES)
      fprintf (ifile, "%20s ", indx->egroup);
    fprintf (ifile, "  %12s", indx->ID);
    fprintf (ifile, " %d %d", 0, 0);	// parents 
    c = indx->gender;
    pgender = 0;
    if (c == 'M')
      pgender = 1;
    if (c == 'F')
      pgender = 2;
    fprintf (ifile, " %d", pgender);
    if (ogmode == NO) {
      astatus = indx->affstatus + 1;
      if (qtmode) {
	fprintf (ifile, "%9.3f", indx->rawqval);
      }
      else {
	fprintf (ifile, " %d", astatus);
      }
    }
    if (ogmode == YES)
      fprintf (ifile, " %10s", indx->egroup);
    if (sevencolumnped)
      fprintf (ifile, " %d", dcode);
    fprintf (ifile, "\n");
  }
  fclose (ifile);
}


/* ---------------------------------------------------------------------------------------------------- */
void
outped (char *snpname, char *indname, char *gname, SNP ** snpm,
	Indiv ** indiv, int numsnps, int numind, int ogmode)
{

  FILE *fff, *ifile;
  int g, i, k;
  SNP *cupt;
  Indiv *indx;
  char c;
  int pgender, astatus;
  int g1, g2, dcode = 1;

  settersemode (YES);
  if (snpname != NULL)
    printmap (snpname, snpm, numsnps, indiv);	// print .map file

  if (indname != NULL)
    outindped (indname, indiv, numind, ogmode);	// print .pedind file

  // Here, printt the .ped file
  if (gname == NULL)
    return;
  openit (gname, &fff, "w");
  for (i = 0; i < numind; i++) {
    indx = indiv[i];
    if (indx->ignore)
      continue;
    fprintf (fff, "%6d %12s", i + 1, indx->ID);	// make up a family name (index) and print individual name
    fprintf (fff, " %d %d", 0, 0);	// set parents to "not in data set"  
    c = indx->gender;
    pgender = 0;
    if (c == 'M')
      pgender = 1;
    if (c == 'F')
      pgender = 2;
    fprintf (fff, " %d", pgender);
    if (ogmode == NO) {
      astatus = indx->affstatus + 1;
      if (qtmode) {
	fprintf (fff, "%9.3f", indx->rawqval);
      }
      else
	fprintf (fff, " %d", astatus);
    }
    if (ogmode == YES)
      fprintf (fff, " %10s", indx->egroup);
    if (sevencolumnped)
      fprintf (fff, " %d", dcode);
    for (k = 0; k < numsnps; k++) {
      cupt = snpm[k];
      if (outputall == NO) {
	if (ignoresnp (cupt))
	  continue;
	if (cupt->isrfake)
	  continue;
      }
      g = getgtypes (cupt, i);
      gtox (g, cupt->alleles, &g1, &g2);
      fprintf (fff, "  %d %d", g1, g2);
      if ((g1 > 4) || (g2 > 4)) {
	fprintf (fff, "\n");
	fflush (fff);
	fclose (fff);
	printf ("bad genotype for snp %s alleles: ", cupt->ID);
	printalleles (cupt, stdout);
	printf ("\n");
	fatalx ("trying to make invalid ped file %s\n", gname);
      }
    }
    fprintf (fff, "\n");
  }
  fclose (fff);
}


/* ---------------------------------------------------------------------------------------------------- */
void
gtox (int g, char *cvals, int *p1, int *p2)
{
  // output values for ped file using allele array
  int g1, g2;

  switch (g) {
  case -1:
    *p1 = *p2 = 0;
    return;
  case 0:
    g1 = 1;
    g2 = 1;
    break;
  case 1:
    g1 = 1;
    g2 = 2;
    break;
  case 2:
    g1 = 2;
    g2 = 2;
    break;
  default:
    fatalx ("(outped) bug %d\n", g);
  }

  if (cvals != NULL) {
    g1 = 3 - g1;
    g2 = 3 - g2;
    g1 = xpedval (cvals[g1 - 1]);
    g2 = xpedval (cvals[g2 - 1]);
  }

  *p1 = MIN (g1, g2);
  *p2 = MAX (g1, g2);

}


/* ---------------------------------------------------------------------------------------------------- */
void
outpackped (char *snpname, char *indname, char *gname, SNP ** snpm,
	    Indiv ** indiv, int numsnps, int numind, int ogmode)
{

  FILE *fff, *ifile;
  int g, i, k;
  SNP *cupt;
  Indiv *indx;
  char c;
  int pgender, astatus;
  int g1, g2, dcode = 1;
  unsigned char ibuff[3];
  unsigned char *buff;
  int fdes, ret, blen;
  int *gtypes;
  double y;

  settersemode (YES);
  if (snpname != NULL)
    printmap (snpname, snpm, numsnps, indiv);	// print .map (not .bim)

  if (indname != NULL)		// print .pedind file
    outindped (indname, indiv, numind, ogmode);

  if (gname == NULL)
    return;

  /*  magic constants for snp major bed file */
  ibuff[0] = 0x6C;
  ibuff[1] = 0x1B;
  ibuff[2] = 0x01;


  // blen is number of bytes each SNP's data requires
  y = (double) (numind * 2) / (8 * (double) sizeof (char));
  blen = nnint (ceil (y));
  ZALLOC (buff, blen, unsigned char);
  ZALLOC (gtypes, numind, int);

  // open output file and check readability
  fdes = open (gname, O_CREAT | O_TRUNC | O_RDWR, 0666);
  if (fdes < 0) {
    perror ("bad gname");
    fatalx ("open failed for %s\n", gname);
  }
  if (verbose)
    printf ("file %s opened\n", gname);

  if (fdes < 0) {
    perror ("bad genoout");
    fatalx ("open failed for %s\n", gname);
  }

  if (verbose)
    printf ("file %s opened\n", gname);

  ret = write (fdes, ibuff, 3);

  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (outputall == NO) {
      if (ignoresnp (cupt))
	continue;
      if (cupt->isrfake)
	continue;
    }
    for (k = 0; k < numind; ++k) {
      g = getgtypes (cupt, k);
      if (g >= 0)
	g = 2 - g;
      gtypes[k] = g;
    }
    setbedbuff ((char *) buff, gtypes, numind);
    ret = write (fdes, buff, blen);
    if (ret < 0) {
      perror ("write failure");
      fatalx ("(outpackped) bad write");
    }
  }

  free (buff);
  close (fdes);
}


/* ---------------------------------------------------------------------------------------------------- */
void
setbedbuff (char *buff, int *gtypes, int numind)
{
  int i, k;
  double y;
  int blen, wnum, wplace, bplace, t, g;
  unsigned char c;

  y = (double) (numind * 2) / (8 * (double) sizeof (char));
  blen = nnint (ceil (y));

  c = 0xAA;			// missing 
  cclear ((unsigned char *) buff, c, blen);	// initialize buffer to "missing"

  for (k = 0; k < numind; k++) {
    wnum = k / 4;
    t = k % 4;
    wplace = 3 - t;		// switch for bed 
    bplace = 4 * wnum + wplace;
    g = bedval (gtypes[k]);
    wbuff ((unsigned char *) buff, bplace, g);
  }
}


/* ---------------------------------------------------------------------------------------------------- */
int
bedval (int g)
{
  if (g < 0)
    return 1;
  if (g == 2)
    return 3;
  if (g == 1)
    return 2;
  if (g == 0)
    return 0;

  fatalx ("(bedval) bad g value %d\n", g);
}


/* ---------------------------------------------------------------------------------------------------- */
void
atopchrom (char *ss, int chrom)
{

  // ancestry chromosome -> map convention  

/**
  if ( chrom == numchrom+1 )  {
    strcpy(ss, "X") ;
    return ;
  }
  else if ( chrom == numchrom+2 )  {
    strcpy(ss, "Y") ;
    return ;
  }
*/
  sprintf (ss, "%d", chrom);
}

/* ---------------------------------------------------------------------------------------------------- */
int
ptoachrom (char *ss)
{
  // map -> ancestry  
  char c;
  c = ss[0];

  if (c == 'X')
    return (numchrom + 1);
  if (c == 'Y')
    return (numchrom + 2);
  return atoi (ss);
}


/* ---------------------------------------------------------------------------------------------------- */
void
printmap (char *snpname, SNP ** snpm, int numsnps, Indiv ** indiv)
{

  char ss[5];
  int i;
  FILE *fff;
  SNP *cupt;
  char c;

  if (snpname == NULL)
    return;
  openit (snpname, &fff, "w");
  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (outputall == NO) {
      if (ignoresnp (cupt))
	continue;
      if (cupt->isrfake)
	continue;
    }
    atopchrom (ss, cupt->chrom);
    fprintf (fff, "%-2s", ss);
    fprintf (fff, " %12s", cupt->ID);
    fprintf (fff, " %12.6f", cupt->genpos);
    fprintf (fff, " %12.0f", cupt->physpos);
    printalleles (cupt, fff);
    fprintf (fff, "\n");
  }
  fclose (fff);
}


/* ---------------------------------------------------------------------------------------------------- */
char
x2base (int x)
{
  // 12345 -> ACGTX
  char *blist = "?ACGT";
  if (x < 0)
    return '?';
  if (x > 4)
    return 'X';
  return blist[x];
}

/* ---------------------------------------------------------------------------------------------------- */
int
xpedval (char c)
{
  char bb[2];

  bb[1] = '\0';
  bb[0] = c;

  if (isdigit (c))
    return atoi (bb);
  return pedval (bb);
}


/* ---------------------------------------------------------------------------------------------------- */
int
pedval (char *sx)
{
  char c;

  c = sx[0];
  if (c == 'A')
    return 1;
  if (c == 'C')
    return 2;
  if (c == 'G')
    return 3;
  if (c == 'T')
    return 4;
  if (c == 'X')
    return 5;
  if (c == 'N')
    return 5;
  if (c == 'N')
    return 5;
  if (c == 'D')
    return 5;
  if (c == 'I')
    return 5;

  if (c == '1')
    return 1;
  if (c == '2')
    return 2;
  if (c == '3')
    return 3;
  if (c == '4')
    return 4;
  if (c == '0')
    return 0;

  if (badpedignore)
    return 5;

  return 9;

}


/* ---------------------------------------------------------------------------------------------------- */
int
getbedgenos (char *gname, SNP ** snpmarkers, Indiv ** indivmarkers,
	     int numsnps, int numindivs, int nignore)
{

  int val, i, k, x, j;
  int t, wnum, wplace;
  int nsnp;
  int ngenos = 0;

  SNP *cupt;
  Indiv *indx;

  unsigned char *buff, ibuff[3], jbuff[3];
  double y;
  int blen;
  int fdes;

  // magic numbers for BED identification
  ibuff[0] = 0x6C;
  ibuff[1] = 0x1B;
  ibuff[2] = 0x01;

  cleargdata (snpmarkers, numsnps, numindivs);
  nsnp = numsnps;

  if (pordercheck && (snpordered == NO))
    failorder ();

  // blen is number of bytes needed to store each SNP's genotype
  y = (double) (numindivs * 2) / (8 * (double) sizeof (char));
  blen = nnint (ceil (y));
  ZALLOC (buff, blen, unsigned char);

  // open binary file and check that it is readable
  fdes = open (gname, O_RDONLY);
  if (fdes < 0) {
    perror ("open failure");
    fatalx ("(getbedgenos) bad open %s\n", gname);
  }
  t = read (fdes, jbuff, 3);
  if (t < 0) {
    perror ("read failure");
    fatalx ("(getbedgenos) bad read");
  }

  // check magic
  for (k = 0; k < 3; k++) {
    if (ibuff[k] != jbuff[k]) {
      fprintf (stderr, "magic failure: ");
      fprintf (stderr, " %x %x %x", jbuff[0], jbuff[1], jbuff[2]);
      fprintf (stderr, " %x %x %x", ibuff[0], ibuff[1], ibuff[2]);
      fprintf (stderr, "\n");
      fatalx ("(getbedgenos) magic failure\n");
    }
  }

  // Read genotype data
  for (i = 0; i < nsnp; i++) {

    j = snpord[i];
    if (snpordered == YES)
      j = i;
    if (j < 0)
      fatalx ("(readbedgenos) bug\n");
    if (j > nsnp)
      fatalx ("(readbedgenos) bug\n");

    cupt = snpmarkers[j];
    t = read (fdes, buff, blen);

    if (t < 0) {
      perror ("read failure");
      fatalx ("(getbedgenos) bad read");
    }
    if (cupt->ignore)
      continue;

    for (k = 0; k < numindivs; k++) {
      indx = indivmarkers[k];
      wnum = k / 4;
      t = k % 4;
      wplace = 3 - t;		// switch for bed 
      wplace += 4 * wnum;
      x = rbuff (buff, wplace);
      val = ancval (x);
      if (checkxval (cupt, indx, val) == NO)
	val = -1;
      putgtypes (cupt, k, val);
      if (val >= 0)
	++ngenos;
    }
  }

  free (buff);
  printf ("genotype file processed\n");
  return ngenos;

}


/* ---------------------------------------------------------------------------------------------------- */
int
ancval (int x)
{
  // bed -> anc 
  // 1/22/07 allele flipped
  if (x == 1)
    return -1;
  if (x == 3)
    return 0;
  if (x == 2)
    return 1;
  if (x == 0)
    return 2;
  fatalx ("(ancval) bad value %d\n", x);
}



/* ---------------------------------------------------------------------------------------------------- */
void
setomode (enum outputmodetype *outmode, char *omode)
{

  char *ss;
  int len, i;

  if (outmode == NULL)
    return;
  *outmode = PACKEDANCESTRYMAP;
  if (omode == NULL)
    return;

  ss = strdup (omode);
  len = strlen (ss);
  for (i = 0; i < len; i++) {
    ss[i] = tolower (ss[i]);
  }

  if (strcmp (ss, "eigenstrat") == 0)
    *outmode = EIGENSTRAT;
  if (strcmp (ss, "ascii") == 0)
    *outmode = EIGENSTRAT;
  if (strcmp (ss, "alkes") == 0)
    *outmode = EIGENSTRAT;
  if (strcmp (ss, "ped") == 0)
    *outmode = PED;
  if (strcmp (ss, "packedped") == 0)
    *outmode = PACKEDPED;
  if (strcmp (ss, "packedancestrymap") == 0)
    *outmode = PACKEDANCESTRYMAP;
  if (strcmp (ss, "ancestrymap") == 0)
    *outmode = ANCESTRYMAP;

  free (ss);
}


/* ---------------------------------------------------------------------------------------------------- */
void
snpdecimate (SNP ** snpm, int nsnp, int decim, int mindis, int maxdis)
{
  int chrom = -1;
  SNP **cbuff, *cupt, *cupt2;
  int k, k2, n, t;

  printf ("snpdecimate called: decim: %d  mindis: %d  maxdis: %d\n", decim,
	  mindis, maxdis);
  ZALLOC (cbuff, nsnp, SNP *);
  for (k = 0; k < nsnp; ++k) {
    cupt = snpm[k];
    if (cupt->chrom != chrom) {
      chrom = cupt->chrom;
      n = 0;
      for (k2 = k; k2 < nsnp; ++k2) {
	cupt2 = snpm[k2];
	if (cupt2->chrom != chrom)
	  break;
	if (cupt2->ignore)
	  continue;
	if (cupt2->isfake)
	  continue;
	cbuff[n] = cupt2;
	++n;
      }
      if (n < decim)
	continue;
      decimate (cbuff, n, decim, mindis, maxdis);
    }
  }
}



/* ---------------------------------------------------------------------------------------------------- */
void
decimate (SNP ** cbuff, int n, int decim, int mindis, int maxdis)
{
  int k, t, u, dis, len;
  int *ttt;
  SNP *cupt;

  cupt = cbuff[0];
  if (n < 2)
    return;
  if (n < decim)
    return;
  ZALLOC (ttt, n, int);
  for (k = 1; k < n; ++k) {
    dis = (int) (cbuff[k]->physpos - cbuff[k - 1]->physpos);
    if (dis > maxdis) {
      decimate (cbuff, k - 1, decim, mindis, maxdis);
      decimate (cbuff + k, n - k, decim, mindis, maxdis);
      return;
    }
  }
  t = ranmod (decim);
  ttt[t] = 1;

  u = t + decim;

  for (;;) {
    if (u >= n)
      break;
    dis = (int) (cbuff[u]->physpos - cbuff[t]->physpos);
    if (dis < mindis) {
      ++u;
      continue;
    }
    len = u - t - 1;
    ivclear (ttt + t + 1, 1, len);
    t = u;
    u = t + decim;
  }
  for (k = 0; k < n; ++k) {
    if (ttt[k] == 1)
      cbuff[k]->ignore = YES;
  }
// debug 
  if (verbose) {
    for (k = 0; k < n; ++k) {
      printf ("zz %6d %20s %20d  %3d\n", k, cbuff[k]->ID,
	      (int) cbuff[k]->physpos, ttt[k]);
    }
  }

  free (ttt);

}


/* ---------------------------------------------------------------------------------------------------- */
int
killhir2 (SNP ** snpm, int numsnps, int numind, double physlim, double genlim,
	  double rhothresh)
{
  // physlim = genlim = 0 => kill monomorphs 
  double *badbuff;
  int *xbadbuff;
  SNP *cupt, *cupt1, *cupt2;

#define BADBUFFSIZE 100000 ;

  int badbuffsize = BADBUFFSIZE;
  int i, j, k, nbad, kmax, kmin, t, j1, j2, lo, hi;
  int *gtypes;
  double *x1, *x2, mean, dis, *p1;
  int nkill = 0, tj;
  double y1, y2, y, rho, smax;
  double **xx1, *yy1;
  SNP **snpxl;

  if (physlim < 0)
    return 0;
  if (genlim < 0)
    return 0;

  // step 1 give score to each SNP 
  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    if (cupt->ignore)
      continue;
    cupt->score = numvalidgtypes (cupt);
    cupt->score += DRAND ();	// jitter
  }

  ZALLOC (badbuff, badbuffsize, double);
  ZALLOC (xbadbuff, badbuffsize, int);
  ZALLOC (x1, numind, double);
  ZALLOC (x2, numind, double);
  ZALLOC (gtypes, numind, int);

  xx1 = initarray_2Ddouble (10000, numind, 0.0);
  ZALLOC (yy1, 10000, double);
  ZALLOC (snpxl, 10000, SNP *);

  for (i = 0; i < numsnps; i += 5000) {

    lo = i;
    hi = i + 10000 - 1;
    hi = MIN (hi, numsnps - 1);

    for (j = lo; j <= hi; ++j) {
      p1 = xx1[j - lo];
      cupt = snpm[j];
      snpxl[j - lo] = cupt;
      grabgtypes (gtypes, cupt, numind);
      floatit (p1, gtypes, numind);
      vvadjust (p1, numind, NULL);
      y1 = asum2 (p1, numind);
      yy1[j - lo] = y1;
      if (y1 < 0.01) {
	++nkill;
	cupt->ignore = YES;
      }
    }
    for (j1 = 0; j1 < 5000; ++j1) {
      if (j1 > (hi - lo))
	break;
      cupt1 = snpxl[j1];
      if (cupt1->ignore)
	continue;
      nbad = 0;
      tj = 0;
      for (j2 = j1 + 1; j2 <= hi - lo; ++j2) {
	cupt2 = snpxl[j2];
	if (cupt2->ignore)
	  continue;
	if (cupt2->chrom != cupt1->chrom)
	  break;

	dis = cupt2->genpos - cupt1->genpos;
	if (dis > genlim)
	  break;

	dis = cupt2->physpos - cupt1->physpos;
	if (dis > physlim)
	  break;
	++tj;

	y1 = yy1[j1];
	y2 = yy1[j2];

	y = vdot (xx1[j1], xx1[j2], numind) / sqrt (y1 * y2);	// compute correlation
	rho = y * y;
	if (rho < rhothresh)
	  continue;
	badbuff[nbad] = cupt2->score;
	xbadbuff[nbad] = j2 + lo;
	++nbad;

      }
      t = (j1 + lo) % 100;
      if (nbad == 0)
	continue;
      vlmaxmin (badbuff, nbad, &kmax, &kmin);
      smax = snpm[kmax]->score;
      if (smax > cupt1->score) {
	cupt1->ignore = YES;
	++nkill;
	continue;
      }
      for (k = 0; k < nbad; ++k) {
	j = xbadbuff[k];
	snpm[j]->ignore = YES;
	++nkill;
      }
    }
  }


  free2D (&xx1, 10000);
  free (yy1);
  free (snpxl);
  free (gtypes);
  free (badbuff);
  free (xbadbuff);
  free (x1);
  free (x2);

  // re-initialize scores
  for (i = 0; i < numsnps; i++) {
    cupt = snpm[i];
    cupt->score = 0.0;
  }

  printf ("killr2 complete\n");
  return nkill;
}


/* ---------------------------------------------------------------------------------------------------- */
int
vvadjust (double *cc, int n, double *pmean)
{
  // take off mean  force missing to zero 
  // simpler version of vadjust 

  double ynum, ysum, y, ymean;
  int i, nmiss = 0;

  ynum = ysum = 0.0;
  for (i = 0; i < n; i++) {
    y = cc[i];
    if (y < 0.0) {
      ++nmiss;
      continue;
    }
    ++ynum;
    ysum += y;
  }
  if (ynum <= 1.5) {
    // no data or monomorphic 
    vzero (cc, n);
    if (pmean != NULL)
      *pmean = ysum / (ynum + 1.0e-8);
    return nmiss;
  }
  ymean = ysum / ynum;
  for (i = 0; i < n; i++) {
    y = cc[i];
    if (y < 0.0)
      cc[i] = 0.0;
    else
      cc[i] -= ymean;
  }
  if (pmean != NULL)
    *pmean = ymean;
  return nmiss;

}




/* ---------------------------------------------------------------------------------------------------- */
static int
setskipit (char *sx)
{
  int skipit = NO;
  if (sx[0] == '#')
    skipit = YES;
  if (strcmp (sx, "SNP_ID") == 0)
    skipit = YES;
  if (strcmp (sx, "Indiv_ID") == 0)
    skipit = YES;
  if (strcmp (sx, "Chr") == 0)
    skipit = YES;
  return skipit;
}



/* ---------------------------------------------------------------------------------------------------- */
int
inpack2 (char *gname, SNP ** snpm, Indiv ** indiv, int numsnps, int numind)
{
  // load up packed genotype file for merge.

  char **arrx, junk[10];
  int n, num, ihash, shash, i, g, j, k, t, g1, g2;
  int xihash, xshash, xnsnp, xnind;
  int nind, nsnp, irec;
  Indiv *indx;
  SNP *cupt, *cupt2;
  SNP xsnp;
  double y;
  unsigned char *buff, *tbuff;
  int fdes, ret;
  char *packit, *pbuff;
  int nbad = 0;
  n = numind;

  ZALLOC (arrx, n, char *);

  // compute hashes to compare with file
  num = 0;
  for (i = 0; i < n; i++) {
    indx = indiv[i];
    arrx[num] = strdup (indx->ID);
    ++num;
  }
  ihash = hasharr (arrx, num);
  nind = num;

  freeup (arrx, num);
  free (arrx);

  n = numsnps;
  ZALLOC (arrx, n, char *);
  num = 0;
  for (i = 0; i < n; i++) {
    cupt = snpm[i];
    arrx[num] = strdup (cupt->ID);
    ++num;
  }
  shash = hasharr (arrx, num);
  nsnp = num;
  freeup (arrx, num);
  free (arrx);

  // rlen  is number of bytes each SNP's data requires in packed format
  y = (double) (nind * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));
  rlen = MAX (rlen, 48);
  ZALLOC (buff, rlen, unsigned char);
  ZALLOC (tbuff, rlen, unsigned char);

  // openfile and check readability
  fdes = open (gname, O_RDONLY);
  if (fdes < 0) {
    perror ("open failure");
    fatalx ("(inpack2) bad open %s\n", gname);
  }
  t = read (fdes, buff, rlen);
  if (t < 0) {
    perror ("read failure");
    fatalx ("(inpack2) bad read");
  }
  sscanf ((char *) buff, "GENO %d %d %x %x", &xnind, &xnsnp, &xihash,
	  &xshash);

  if (xnind != nind)
    fatalx ("(inpack2) nind mismatch %d %d \n", nind, xnind);
  if (xnsnp != nsnp)
    fatalx ("(inpack2) nsnp mismatch\n");
  if (xihash != ihash)
    fatalx ("(inpack2) ihash mismatch\n");
  if (xshash != shash)
    fatalx ("(inpack2) shash mismatch\n");


  // now copy genotypes 
  for (i = 0; i < n; i++) {
    t = read (fdes, tbuff, rlen);
    if (t != rlen) {
      perror ("read failure");
      fatalx ("(inpack2) bad data read");
    }
    cupt = snpm[i];
    if (cupt->isfake)
      continue;
    xsnp = *cupt;
    cupt2 = &xsnp;
    cupt2->pbuff = (char *) tbuff;
    for (k = 0; k < numind; ++k) {
      g2 = getgtypes (cupt2, k);	// store in temporary buffer
      if (g2 < 0)
	continue;
      g1 = getgtypes (cupt, k);
      if ((g1 >= 0) && (g1 != g2))
	++nbad;			// inconsistent data
      putgtypes (cupt, k, g2);
    }

    // now check xhets
    for (k = 0; k < numind; ++k) {
      if (cupt->chrom != (numchrom + 1))
	break;
      indx = indiv[k];
      g = getgtypes (cupt, k);
      if (checkxval (cupt, indx, g) == NO) {
	putgtypes (cupt, k, -1);
      }
    }
  }

  free (buff);
  free (tbuff);
  close (fdes);
  return nbad;
}

/* ---------------------------------------------------------------------------------------------------- */
void
getgenos_list (char *genotypelist, SNP ** snpmarkers, Indiv ** indivmarkers,
	       int numsnps, int numindivs, int nignore)
{

  char **fnames, *fn;
  int n;
  int k, nbad, isok;

  dofreeped = NO;
  n = numlines (genotypelist);
  ZALLOC (fnames, n, char *);

  // Read in list of genotype files
  n = getlist (genotypelist, fnames);

  // Load first one the ordinary way
  getgenos (fnames[0], snpmarkers, indivmarkers, numsnps, numindivs, nignore);

  // Load all others 
  for (k = 1; k < n; ++k) {
    fn = fnames[k];
    isok = NO;
    if (ispack (fn)) {
      nbad = inpack2 (fn, snpmarkers, indivmarkers, numsnps, numindivs);
      isok = YES;
    }
    if (iseigenstrat (fn)) {
      nbad = ineigenstrat (fn, snpmarkers, indivmarkers, numsnps, numindivs);
      isok = YES;
    }
    if (nbad > 0)
      printf ("%s genotypes mismatches: %d\n", fn, nbad);
    if (isok == NO)
      fatalx ("file %s must be packed or eigenstrat format\n");
  }

  dofreeped = YES;
  freeped ();

  free (fnames);
}



/* ---------------------------------------------------------------------------------------------------- */
int
setsdpos (SNPDATA * sdpt, int pos)
{
  int t;
  char ss[10], *sx;

  sdpt->ppos = pos;
  strcpy (ss, sdpt->cchrom);
  mkupper (ss);

  sdpt->chimpfudge = chimpmode;

  sx = strstr (ss, "CHR");
  if (sx != NULL)
    sx = ss + 3;
  else
    sx = ss;

  t = strcmp (sx, "2B");
  if (t == 0) {
    sdpt->ppos += 200000000;
    sdpt->chimpfudge = YES;
  }
  t = strcmp (sx, "2A");
  if (t == 0) {
    sdpt->chimpfudge = YES;
  }
  return sdpt->chimpfudge;
}

int
str2chrom (char *sss)
{
  char ss[6];
  if (strlen (sss) > 5)
    fatalx ("bad chrom: %s\n", sss);
  if (strstr (sss, "chr") != NULL) {
    strcpy (ss, sss + 3);
    setchr (YES);
  }
  else
    (strcpy (ss, sss));
  mkupper (ss);
  if (strcmp (ss, "X") == 0)
    return (numchrom + 1);
  if (strcmp (ss, "Y") == 0)
    return (numchrom + 2);
  if (strcmp (ss, "MT") == 0)
    return MTCHROM;
  if (strcmp (ss, "XY") == 0)
    return XYCHROM;
  if (strcmp (ss, "2A") == 0)
    return 2;
  if (strcmp (ss, "2B") == 0)
    return 2;
  if (!isnumword (ss))
    return -1;
  return atoi (ss);
}


/* ---------------------------------------------------------------------------------------------------- */
int
checksize (int numindivs, int numsnps, enum outputmodetype outputmode)
{
  // -1 try packed format

  double y;
  long z;

  if (sizeof (z) == 8)
    checksizemode = NO;
  if (checksizemode == NO)
    return 1;

  y = (double) numindivs;
  y *= (double) numsnps;

  if (y > 8.0e9)
    return -2;

  switch (outputmode) {

  case ANCESTRYMAP:
    if (y > 5.0e7)
      return -1;
    break;
  case EIGENSTRAT:
    if (y > 2.0e9)
      return -1;
    break;
  case PED:
    if (y > 4.0e8)
      return -1;
    break;
  case PACKEDPED:
    break;
  case PACKEDANCESTRYMAP:
    break;
  default:
    fatalx ("unknown outputmode\n");
  }
  return 1;
}

/* ---------------------------------------------------------------------------------------------------- */
int
snprawindex (SNPDATA ** snpraw, int nreal, char *sname)
{
  int k;
  char **ss;

  freesnpindex ();

  // if hash table is not set up, do it now
  if (snprawtab == NO) {
    snprawtab = YES;
    ZALLOC (ss, nreal, char *);
    for (k = 0; k < nreal; k++) {
      ss[k] = strdup (snpraw[k]->ID);
    }

    // hash SNP data (key=SNP name, data=index in snpraw)
    xloadsearch (ss, nreal);
    freeup (ss, nreal);
    free (ss);
  }

  // return index in snpraw    
  k = xfindit (sname);
  return k;
}



/* ---------------------------------------------------------------------------------------------------- */
void
freesnprawindex ()
{
  if (snprawtab == NO)
    return;
  snprawtab = NO;
  xdestroy ();
}


/* ---------------------------------------------------------------------------------------------------- */
void
cntpops (int *count, Indiv ** indm, int numindivs, char **eglist, int numeg)
{
  // count number of samples for each pop
  Indiv *indx;
  int t, k;

  ivzero (count, numeg);
  for (k = 0; k < numindivs; ++k) {
    indx = indm[k];
    if (indx->ignore)
      continue;
    t = indxindex (eglist, numeg, indx->egroup);
    if (t < 0)
      continue;
    ++count[t];
  }
}

/* ---------------------------------------------------------------------------------------------------- */
char *
getpackgenos ()
{
  return packgenos;
}


/* ---------------------------------------------------------------------------------------------------- */
void
clearpackgenos ()
{
  packgenos = NULL;
}


/* ---------------------------------------------------------------------------------------------------- */
void
genocloseit (genofile * gfile)
{

  genofile *gpt;
  SNP *cupt;
  int i;
  gpt = gfile;

  free (gpt->buff);
  for (i = 0; i < gpt->numsnps; i++) {
    cupt = gpt->snpm[i];
    freecupt (&cupt);
  }
  free (gpt->snpm);

  for (i = 0; i < gpt->numindivs; i++) {
    free (gpt->indivm[i]);
  }
  free (gpt->indivm);

  close (gpt->fdes);

}

/* ---------------------------------------------------------------------------------------------------- */
int
genoopenit (genofile ** gfile, char *geno2name, SNP ** snp2m,
	    Indiv ** indiv2m, int numsnp2, int numindiv2, int nignore)
{
  // only one gfile can be open 
  static genofile xfile;
  genofile *gpt;
  double y;
  int rlen, fdes, t;
  static unsigned char *buff;
  int xihash, xshash, xnsnp, xnind;
  int ihash, shash;
  char *gname;
  int nsnp, nind;


  if (geno2name == NULL)
    fatalx ("(genoopenit) null name\n");
  if (!ispack (geno2name))
    fatalx ("(genoopenit) not packed ancestrymap format\n");
  gpt = *gfile = &xfile;
  strcpy (gpt->gname, geno2name);
  gpt->snpm = snp2m;
  gpt->indivm = indiv2m;
  gpt->numsnps = numsnp2;
  gpt->numindivs = numindiv2;

  y = (double) (numindiv2 * 2) / (8 * (double) sizeof (char));
  rlen = nnint (ceil (y));

  gpt->rlen = rlen;
  rlen = MAX (rlen, 48);
  ZALLOC (buff, rlen, unsigned char);
  gpt->buff = buff;

  fdes = open (geno2name, O_RDONLY);
  if (fdes < 0)
    return fdes;
  gpt->fdes = fdes;
  gpt->snpindex = -1;

  t = read (fdes, buff, rlen);
  if (t < 0)
    fatalx ("(genoopenit) bad initial read\n");

  nsnp = numsnp2;
  nind = numindiv2;
  gname = geno2name;

  calcishash (snp2m, indiv2m, nsnp, nind, &ihash, &shash);
  if (hashcheck) {
    sscanf ((char *) buff, "GENO %d %d %x %x", &xnind, &xnsnp, &xihash,
	    &xshash);
    if (xnind != nind)
      fatalx ("OOPS number of individuals %d != %d in input files\n", nind,
	      xnind);
    if (xnsnp != nsnp)
      fatalx ("OOPS number of SNPs %d != %d in input file: %s\n", nsnp, xnsnp,
	      gname);
    if (xihash != ihash)
      fatalx
	("OOPS indiv file has changed since genotype file was created\n");
    if (xshash != shash)
      fatalx ("OOPS snp file has changed since genotype file was created\n");
  }

  return 0;


/* (Real definition is in admutils.h) 

typedef struct {  
 char gname[IDSIZE] ;  
 SNP **snpm ;
 Indiv **indivm ;  
 int numsnps; 
 int numindivs ; 
 int rlen ;
 int fdes ; 
 unsigned char *buff ;
 int snpindex ;
} genofile ;
*/

}

/* ---------------------------------------------------------------------------------------------------- */
int
genoreadit (genofile * gfile, SNP ** pcupt)
{
/*  
 return code 
 < 0 bad read  
 0 EOF 
 rlen good read 
*/
  genofile *gpt;
  SNP *cupt;
  int t, rlen, snum;
  int k;

  cupt = *pcupt = NULL;
  gpt = gfile;
  rlen = gpt->rlen;
  t = read (gpt->fdes, gpt->buff, rlen);
  if (t < 0)
    fatalx ("(genoreadit) bad read \n");
  if (t == 0)
    return 0;
  if (t < gpt->rlen)
    fatalx ("(genoopenit) premature EOF\n");
  ++gpt->snpindex;
  snum = gpt->snpindex;
  cupt = *pcupt = gpt->snpm[snum];
  cupt->tagnumber = snum;
  cupt->pbuff = (char *) gpt->buff;
  cupt->ngtypes = gpt->numindivs;
  if (cupt->gtypes == NULL)
    ZALLOC (cupt->gtypes, 1, int);
  return rlen;
}


/* ---------------------------------------------------------------------------------------------------- */
void
putped (int num)
{
  int *pp;
  int t;

  pp = snporda[num];
  if (pp != NULL)
    free (pp);
  pp = NULL;
  t = numsnporda[num] = numsnpord;
  if (t == 0)
    return;
  ZALLOC (snporda[num], t, int);
  pp = snporda[num];
  copyiarr (snpord, pp, t);
}

/* ---------------------------------------------------------------------------------------------------- */
void
getped (int num)
{
  int *pp;
  int t;

  pp = snpord;
  if (pp != NULL)
    free (pp);
  pp = NULL;
  t = numsnpord = numsnporda[num];
  if (t == 0)
    return;
  ZALLOC (snpord, t, int);
  pp = snpord;
  copyiarr (snporda[num], pp, t);
}

/* ---------------------------------------------------------------------------------------------------- */
void
setbadpedignore ()
{
  badpedignore = YES;
}

/* ---------------------------------------------------------------------------------------------------- */
void
logdeletedsnp (char *snpname, char *cmnt, char *deletesnpoutname)
{
  if (deletesnpoutname != NULL) {
    FILE *fid = fopen (deletesnpoutname, "a");
    fprintf (fid, "%-40s  %-40s\n", snpname, cmnt);
    fclose (fid);
  }
}


void
sortsnps (SNP ** snpa, SNP ** snpb, int n)
{
  SNP **tsnp, *cupt;
  int **snppos, *snpindx;
  int i, k;

  snppos = initarray_2Dint (n, 3, 0);
  ZALLOC (snpindx, n, int);
  ZALLOC (tsnp, n, SNP *);

  for (i = 0; i < n; i++) {
    cupt = snpa[i];
    snppos[i][0] = cupt->chrom;
    snppos[i][1] = nnint ((cupt->genpos) * GDISMUL);
    snppos[i][2] = nnint (cupt->physpos);
  }

  ZALLOC (snpindx, n, int);
  ipsortit (snppos, snpindx, n, 3);

  for (i = 0; i < n; ++i) {
    k = snpindx[i];
    tsnp[i] = snpa[k];
  }

  for (i = 0; i < n; ++i) {
    snpb[i] = tsnp[i];
  }

  free (snpindx);
  free2Dint (&snppos, n);
  free (tsnp);

}


int
setstatuslist (Indiv ** indm, int numindivs, char **smatchlist, int slen)
// return number set 
{
  int i, n, k;
  Indiv *indx;
  char *sx;

  n = 0;

  for (i = 0; i < numindivs; i++) {
    indx = indm[i];
    if (indx->ignore)
      continue;
    sx = indx->egroup;
    if (sx == NULL)
      continue;
    k = indxindex (smatchlist, slen, sx);
    if (k < 0)
      continue;
    indx->affstatus = k + 1;
    ++n;
  }
}










/* doxygen documentation */

/*! \fn int getsnps(char *snpfname, SNP ***snpmarkpt, double spacing,
                    char *badsnpname, int *numignore, int numrisks) 

    \brief Read SNP data from file
    \param snpfname   File name (.snp or .map)
    \param snpmarkpt  Pointer to array of type SNP * to store data in
    \param spacing 
    \param badsnpname Name of file with list of SNPs to ignore (or NULL for none)
    \param numignore  ???
    \param numrisks   ???

    Returns number of SNPs loaded
 */


/*!  \fn  int readsnpdata(SNPDATA **snpraw, char *fname)   
     \brief Read SNP file 
     \param snpraw  Array of (pointers to) type SNPDATA in which to temporarily store data
     \param fname   Name of SNP file

     For each SNP read in, stores data in one element (SNPDATA *) of snpraw.  
     Fills these elements of SNPDATA struct : inputrow (own index), chrom, gpos, ppos, alleles
     Also sets maxgpos[chrom] to highest genetic position in chromosome
 */


/*!  \fn  int readsnpmapdata(SNPDATA **snpraw, char *fname)   
     \brief Read PLINK format SNP file 
     \param snpraw  Array of (pointers to) type SNPDATA in which to temporarily store data
     \param fname   Name of SNP file

     For each SNP read in, stores data in one element (SNPDATA *) of snpraw.  
     Fills these elements of SNPDATA struct : inputrow (own index), chrom, gpos, ppos, alleles
     Also sets maxgpos[chrom] to highest genetic position in chromosome
 */


/*! \fn getsizex(char *fname)
    \brief Count number of non-comment lines in file
    This is the number of SNPs in a .map or .snp file
 */

/*! \fn int ismapfile(char *fname)
    \brief Look at file name to determine whether this is a PLINK .map file
    File is assumed to be PLINK if file extension is .map, .bim or .pedsnp
 */

/*! \fn int ispedfile(char *fname)
    \brief Look at file name to determine whether this is a PLINK .ped file.
    File is assumed to be PLINK if file extension is .ped or .fam
 */

/*! \fn int isbedfile(char *fname)
    \brief Look at file name to determine whether this is a PLINK .bed file.
    File is assumed to be PLINK if file extension is .ped or .fam
 */

/*! \fn static int setskipit(char *sx)   
    \brief Determine whether an input line from the SNP file should be skipped
    \param sx is the first token on the input line
    Skip if this is a comment or a line of column headers.
  /

/*! \fn int numfakes(SNPDATA **snpraw, int *snpindx, int nreal, double spacing)    
    \brief Return (approximate) number of fake SNPs that will be inserted

    Note:  for EIGENSOFT programs, spacing is always set to 0.0 (presumably not for 
    ANCESTRYMAP) which results in a return value of 0.  The fake SNPs are inserted
    so that the genetic distance between adjacent SNPs is not greater than spacing.
 */

/*! \fn  double nextmesh(double val, double spacing)  
    \brief Return least multiple of spacing greater than or equal to val
    (Used by numfakes and loadsnps to count number of fake SNPs.)
 */

/*! \fn double interp (double l, double r, double x, double al, double ar)   
    \brief Return linear interpolant a fractional x between the points (l,al) and (r, ar)
 */


/*! \fn int getindivs(char *indivfname, Indiv ***indmarkpt)    
    \brief Read individual data from file
    \param indivfname  File name
    \param indmarkpt   Pointer to array of type Indiv * in which data is to be stored.
 */

/*! \fn int readindpeddata(Indiv **indivmarkers, char *fname)   {
    \brief Read individual data from file
    \param indivfname  File name
    \param indmarkpt   Pointer to array of type Indiv * in which data is to be stored.
 */

/*! \fn void pedname(char *cbuff, char *sx0, char *sx1)    
    \brief Enforce name length requirements and prepend family names if desired.
 */


/*! \fn int readtldata(Indiv **indivmarkers, int numindivs, char *inddataname)    
    Not used in EIGENSOFT
 */

/*! \fn int readfreqdata(SNP **snpm, int numsnps, char *inddataname)     
    Not used in EIGENSOFT
 */

/*!  \fn int setstatus(Indiv **indm, int numindivs, char *smatch)    
     \brief Call setstatusv with value YES
 */

/*! \fn int setstatusv(Indiv **indm, int numindivs, char *smatch, int val)    
    \brief Set affstatus of all individuals with egroup equal to smatch to value val
    \param indm  Array in which individuals' data is stored
    \param numindivs  Number of individuals in the array
    \param  smatch  String in individual's field egroup to match
    \param  val     Value to set affstatus to
 */

/*! \fn int checksize(int numindivs, int numsnps, enum outputmodetype outputmode) 
    \brief Enforce size limits on genotype data file
  */

/*! \fn int ispack(char *gname)    
    \brief Open file and look for GENO at top.  If it's there, the file is packed (binary)
 */

/*! \fn int iseigenstrat(char *gname)    
    \brief If every line in the file is one "word" (i.e., no white space), the file is 
    assumed to be EIGENSTRAT format
 */

/*! \fn long getgenos(char *genoname, SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, int nignore)    
    \brief Read genotype data from file
    \param genoname      Name of genotype data file
    \param snpmarkers    Array in which SNP data is stored
    \param indivmarkers  Array in which individual data is stored
    \param numsnps       Number of SNPs in snpmarkers
    \param numindivs     Number of individuals in indivmarkers
    \param nignore       ???

    Returns number of genotypes read 

    Note:  Instantiates, uses and destroys the hash table.
 */

/*! \fn void genopedcnt(char *gname, int **gcounts, int nsnp)    
    \brief   Count number of alleles of each type in each SNP

    Return in gcounts[k][n] is number of "n" alleles in SNP k

    (This is used to discover and designate ref/var alleles) 

 */

/*! \fn void setgref(int **gcounts, int nsnp, int *gvar, int *gref)    
    \brief  Designate reference and variant alleles by looking at allele counts

    (This is for use with PED files)

 */

/*! \fn long getgenos(char *genoname, SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, int nignore)    
    \brief Read genotype data from PLINK .ped format file

    \param gname         Name of genotype data file
    \param snpmarkers    Array in which SNP data is stored
    \param indivmarkers  Array in which individual data is stored
    \param numsnps       Number of SNPs in snpmarkers
    \param numindivs     Number of individuals in indivmarkers
    \param nignore       ???

    Returns number of genotypes read 

    Note:  Instantiates, uses and destroys the hash table.
 */

/*! \fn int getbedgenos(char *gname, SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, int nignore)    
     \brief Read genotype data from PLINK .bed format file

    \param gname         Name of genotype data file
    \param snpmarkers    Array in which SNP data is stored
    \param indivmarkers  Array in which individual data is stored
    \param numsnps       Number of SNPs in snpmarkers
    \param numindivs     Number of individuals in indivmarkers
    \param nignore       ???

    Returns number of genotypes read 

    Note:  Instantiates, uses and destroys the hash table.
 */

 /*! \fn int inpack(char *gname, SNP **snpm, Indiv **indiv, int numsnps, int numind)    
    \brief Read genotype data from packed ANCESTRYMAP format file

    \param gname         Name of genotype data file
    \param snpmarkers    Array in which SNP data is stored
    \param indivmarkers  Array in which individual data is stored
    \param numsnps       Number of SNPs in snpmarkers
    \param numindivs     Number of individuals in indivmarkers

    Returns number of genotypes read 

    Note:  Instantiates, uses and destroys the hash table.
  */

/*! \fn void cleargdata(SNP **snpmarkers, int numsnps, int numindivs)
    \brief Wipe out all genotype data
 */

/*! \fn rmindivs(SNP ** snpm, int numsnps, Indiv **indivmarkers, int numindivs)
    \brief squeeze out ignored individuals

    Return number of individuals retained (not ignored)
 */


/*! \fn void freecupt(SNP **cuppt) 
    \brief Free memory associated with SNP *
 */

/*! \fn void clearind(Indiv **indm, int numind)
    \brief  Re-initialize all individuals
 */

/*! \fn void cleartg(indiv **indm, int nind 
    \brief  Zero out totgamms and totscore fields for all individuals
 */

/*! \fn void dobadsnps(SNPDATA **snpraw, int nreal, char *badsnpname)
    \brief  Read badsnps file and set ignore flag on all bad SNPs
    \param snpraw  Array of initial SNP data structs
    \param nreal   Number of elements in snpraw
    \param badsnpname  Name of badsnp file
 */

/*!  \fn void printsnps(char *snpoutfilename, SNP **snpm, int num, Indiv **indm, int printfake, int printvalids)
     \brief Print SNP output in EIGENSTRAT format
 */

/*! \fn void printalleles(SNP *cupt, FILE *fff) 
    \brief print SNP's alleles 
 */

/*!  \fn void outfiles(char *snpname, char *indname, char *gname, SNP **snpm, Indiv **indiv, 
                int numsnps, int numindx, int packem, int ogmode)
     \brief Determine which output function to call based on user parameter outputmode
     \param snpname   SNP output file name
     \param indname   Individual output file name
     \param gname     Genotype output file name
     \param snpm      Array of SNP data
     \param indiv     Array of individual data
     \param numsnps   Number of elements in snpm
     \param numindx   Number of elements in indiv
     \param packem    not used (used as local variable)
     \param ogmode    flag for PED, print quantitative or group phenotype
 */

/*! \fn void outpack(char *genooutfilename, SNP **snpm, Indiv **nindiv, int numsnps, int numind)
    \brief Print out genotype data in packed ANCESTRYMAP format
    \param genooutfilename  Genotype output file name
    \param snpm  Array of SNP data
    \param indiv Array of individual data
    \param numsnps  Number of elements in snpm
    \param numind   Number of elements in indiv
*/

/*! \fn void outeigenstrat(char *snpname, char *indname, char *gname, SNP **snpm, Indiv **indiv,
  int numsnps, int numind)
    \brief  Print output in EIGENSTRAT format
    \param snpname  SNP output file name
    \param indname  Individual output name
    \param gname    Genotype output name
    \param snpm     Array of SNP data
    \param indiv    Array of individual data
    \param numsnps  Number of elements in snpm
    \param numind   Number of elements in indiv
 */


/*!  \fn void outped(char *snpname, char *indname, char *gname, SNP **snpm, Indiv **indiv, int numsnps, int numind, int ogmode)
     \brief Print output in (unpacked) PED format

     \param snpname  SNP output file name
     \param indname  Individual output file name
     \param gname    Genotype output file name
     \param snpm     Array with SNP data
     \param indiv    Array with individual data
     \param numsnps  Number of elements in snpm
     \param numind   Number of individuals in indiv
     \param ogmode   phenotype output mode (quantitative or discrete)
 */



/*!  \fn void gtox(int g, char *cvals, int *p1, int *p2)
     \brief Get alleles in PED format
     \param g  Diploid genotype (0,1,2 or -1 for "missing")
     \param cvals  Array with char ref and var alleles
     \param p1  Output for first allele
     \param p2  Output for second allele

     If cvals is NULL, return alleles "1" and "2" (i.e., "A" and "C")
     Otherwise, look up actual alleles.  If the diploid is het, the alleles will be printed in the
     order (ref,var) not (var,ref)
 */

/*!  \fn void outindped(char *indname, Indiv **indiv, int numind, int ogmode)
     \brief  Print out individual names in PEDIND format (i.e., first six or seven columns of PED)
     \param indname  Individual output file name
     \param indiv    Array with individual data
     \param numind   Number of elements in indiv
     \param ogmode   Flag for phenotype type (quantitative or discrete)
 */

/*!  \fn void outpackped(char *snpname, char *indname, char *gname, SNP **snpm, Indiv **indiv, 
  int numsnps, int numind, int ogmode)
     \brief Print data in packed PED format (.bed)
     \param snpname  Output SNP file name
     \param indname  Output individual file name
     \param gname    Output genotype file name
     \param snpm     Array with SNP data
     \param indiv    Array with individual data
     \param numsnps  Number of elements in snpm
     \param numind   Number of individuals in indiv
     \param ogmode   Flag for phenotype type (quantitative or discrete)
 */


/*!  \fn void setbedbuff(char *buff, int *gtypes, int numind)  
     \brief Fill buffer with diploid genotypes in BED format
 */


/*!  \fn void bedval(int g)
     \brief Encode diploid genotype into its packed BED equivalent
 */


/*!  \fn void atopchrom(char *ss, int chrom)
     \brief Encode integer chromosome number to its MAP file equivalent
     \param  ss     output chromosome symbol (0-22 or "X" or "Y")  CHANGED 23 24
     \param  chrom  input chromosome symbol (0-24)
 */

/*!  \fn int ptoachrom(char *ss)
     \brief Encode MAP chromosome symbol to its numerical equivalent
     \param  ss     input chromosome symbol (0-22 or "X" or "Y")

     Return chromosome number (0-24)
 */


/*!  \fn void printmap(char *snpname, SNP **snpm, int numsnps, Indiv **indiv)
     \brief Print out SNP data in PLINK .map format
     \param snpname  Output SNP file name
     \param snpm     Array with SNP data
     \param numsnps  Number of elements in snpm
     \param indiv    not used
 */


/*!  \fn int calcishash(SNP **snpm, Indiv **indiv, int numsnps, int numind, int *pihash, int *pshash)
     \brief Calculate hashes on individuals and SNPs (to compare with file values.)
     \param snpm  Array of SNP data
     \param indiv Array if individual data
     \param numsnps  Number of elements in snpm
     \param numind  Number of elements in indiv 
     \param pihash  Output parameter for indiv hash
     \param pshash  Output parameter for SNP hash

     Return number of SNPs plus number if individuals
 */

/*!  \fn void freeped(void)
     \brief destructor for array snpord
 */

/*!  \fn int readinddata(Indiv **indivmarkers, char *fname)
     \brief  Read individual data from input file
     \param indivmarkers  Array to store data in
     \param fname         Individual input file
 */

/*! \fn int readtldata(Indiv **indivmarkers, int numindivs, char *inddataname)
    \brief Read theta/lambda data (for ANCESTRYMAP, not used in EIGENSOFT)
 */

/*! \fn int readfreqata(SNP **snpm, int numsnps, char *inddataname)
    \brief Read allele frequency data (for ANCESTRYMAP, not used in EIGENSOFT)
 */

/*! \fn int checkxval(SNP *cupt, Indiv *indx, int val)
    \brief Check that male X marker is not (invalidly) heterozygous
 */

/*! \fn void clearsnp(SNP *cupt)
    \brief Reinitialize all fields in SNP data structure
 */


/*! \fn int rmindivs(SNP **snpm, int numsnps, Indiv **indivmarkers, int numindivs)
    \brief Squeeze out individuals with ignore flag set.
    \param snpm      Array of SNP data
    \param numsnps   Number of elements in snpm
    \param indivmarkers  Array of individual data
    \param numindivs     Number of elements in indivmarkers
 */


/*! \fn int rmsnps(SNP **snpm, int numsnps)
    \brief Squeeze out SNPs with ignore flag set
    \param snpm      Array of SNP data
    \param numsnps   Number of elements in snpm
 */

/*! \fn void freecupt(SNP **cuppt)
    \brief Free memory associated with SNP data structure
 */

/*! \fn void clearind(Indiv **indm, int numind)
    \brief Free memory associated with all Indiv data structs in array
 */

/*! \fn void cleartg(Indiv **indm, int nind)
    \brief Free memory in two fields of all Indiv data structs in array
 */

/*! \fn void setug(Indiv **indm, int numind, char gender)
    \brief  Set all unknown gender fields to value passed in
 */

/*! \fn void dobadsnps(SNPDATA **snpraw, int nreal, char *badsnpname)
    \brief Remove bad SNPs from array
    \param snpraw     Array of (preliminary) SNP data
    \param nreal      Number of elements in snpraw
    \param badsnpname  Bad SNP file name
*/

/*! \fn int checkfake(char **ss)
    \brief Returns YES if and only if string ss is "fake"
 */

/*! \fn void printsnps(char *snpoutfilename, SNP **snpm, int num, Indiv **indm, int printfake, int printvalids)
    \brief Print ANCESTRYMAP format SNP file.
    \param snpoutfilename    Name of SNP output file
    \param snpm              Array with SNP data
    \param num               Number of SNPs in array
    \param indm              Array with individual data
    \param printfakes        Flag to print fake SNPs
    \param printvalids       Flag to print alleles
 */

/*! \fn void printalleles(SNP *cupt, FILE *fff)
    \brief  Print SNP alleles to file
 */

/*! \fn void printdata(char *genooutfilename, char *indoutfilename, SNP **snpm, 
   Indiv **indiv, int numsnps, int numind, int packem)
    \brief  Print ANCESTRYMAP format genotype file
    \param genooutfilename    Genotype output file name
    \param indoutfilename     Individual output file name
    \param snpm               Array of SNP data
    \param indiv              Array of individual data
    \param numsnps            Number of elements in snpm
    \param numind             Number of elements in indiv
    \param packem             Flag - print in packed mode if set
*/

/*! \fn void outpack(char *genooutfilename, SNP **snpm, Indiv **indiv, int numsnps, int numind)
    \brief  Print packed ANCESTRYMAP format genotype file
    \param genooutfilename    Genotype output file name
    \param snpm               Array of SNP data
    \param indiv              Array of individual data
    \param numsnps            Number of elements in snpm
    \param numind             Number of elements in indiv
*/

/*! \fn int ineigenstrat(char *gname, SNP **snpm, Indiv **indiv, int numsnps, int numind)
    \brief Read EIGENSTRAT genotype file
    \param gname    Genotype input file
    \param snpm     Array of SNP data
    \param indiv    Array of individual data
    \param numsnps            Number of elements in snpm
    \param numind             Number of elements in indiv

    Return number of errors encountered 

 */

/*! \fn void clearepath(char *packp)
    \brief Fill memory with 0xFF
 */


/*! \fn void getsnpsc(char *snpscname, SNP **snpm, int numsnps)
    \brief Read SNP score input file (not used in EIGENSOFT)
 */

/*! \fn void setepath(SNP **snpm, int nsnps)
    \brief Clear packed genotype memory (i.e., set to "missing") and point SNP buffers-pointers to the SNP's position in packed memory.
 */

/*! \fn int getpedgenos(char *gname, SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs, int nignore)
    \brief  Read PLINK format genotype file
    \param  gname   Name of genotype input file
    \param  snpmarkers  Array of SNP data
    \param  indivmarkers  Array of individual data
    \param  numsnps       Number of SNPs in snpmarkers
    \param  numindivs     Number of individuals in indivmarkers
    \param  nignore       (not used)
 */


/*! \fn void setgenotypename(char **gname, char *iname)
    \brief Copy PED genotypename from iname to *gname, checking that iname is not "NULL."
 */


/*! \fn int maxlinelength(char *fname)
    \brief Find and return length of longest line in file
 */


/*! \fn char x2base(int x)
    \brief Encode digit to char-type allele (PLINK convention)
 */

/*! \fn int xpedval(char c)
    \brief Encode char-type allele to digit (PLINK convention) 
 */

/*! \fn int pedval(char *sx)  
    \brief Encode char-type allele to digit (PLINK convention)
 */

/*! \fn int ancval(int x)
    \brief Encode BED allele digit to ANCESTRYMAP equivalent
 */

/*! \fn void setomode(enum outputmodetype *outmode, char *omode)
    \brief Set output mode from user parameter omode (default is packed ANCESTRYMAP)
 */

/*! \fn void decimate(SNP **cbuff, int n, int decim, int mindis, int maxdis)
    \brief (Undocumented feature)  Prune SNPs 
*/

/*! \fn void snpdecimate(SNP **snpm, int nsnp, int decim, int mindis, int maxdis)
    \brief (Undocumented feature)  Prune SNPs 
*/

/*! \fn int killhir2(SNP **snpm, int numsnps, int numind, double physlim, double genlim, double rhothresh)
    \brief Remove one of each pair of SNPs with r-squared greater than rhothresh 
    \param snpm     Array of SNP data
    \param numsnps  Number of SNPs in snpm
    \param numind   Number of individuals in each SNP's genotype data
    \param physlim  Only consider SNP pairs closer than this
    \param genlim   Only consider SNP pairs closer than this
    \param rhothresh  Maximum permissible r-squared value
*/


/*! \fn int vvadjust(double *cc, int n, double *pmean)
    \brief Mean-adjust data in array and force missing data to zero
    \param cc  Array of values to mean-adjust
    \param n   Number of values in array
    \param pmean Output parameter for returning the mean
 */


/*! \fn int inpack2(char *gname, SNP **snpm, Indiv **indiv, int numsnps, int numind)
    \brief Load packed genotype file for merge of genotype files (used by getgenos_list)
    \param gname    Name of input genotype file
    \param snpm     Array of SNP data
    \param indiv    Array of individual data
    \param numsnps  Number of SNPs in snpm
    \param numind   Number of individuals in indiv
 */


/*! \fn void getgenos_list(char *genotypelist, SNP **snpmarkers, Indiv **indivmarkers, int numsnps,
  int numindivs, int nignore)
    \brief (Undocumented feature) Read in data from all genotype files in a list 
    \param genotypelist  File with names of genotype files in it
    \param snpmarkers    Array of SNP data
    \param indivmarkers  Array of individual data
    \param numsnps       Number of SNPs in snpmarkers
    \param numindivs     Number of individuals in indivmarkers
    \param nignore       (not used)
 */


/*! \fn int str2chrom(char *ss)
    \brief Encode string representation of chromosome to digit equivalent
 */


/*! \fn int snprawindex(SNPDATA **snpraw, int nreal, char *sname)
    \brief Return index of SNP with name sname in array snpraw
 */


/*! \fn void freesnprawindex()
    \brief Free hash table used to look up indices in snpraw
 */

/*! \fn void cntpops(int *count, Indiv **indm, int numindivs, char **eglist, int numeg)
    \brief Count number of samples in each population
    \param count Array in which to store counts
    \param indm  Array of individual data
    \param numindivs  Number of individuals in indm
    \param eglist     Array of population names
    \param numeg      Number of individuals in eglist 
 */

/*! \fn int genoopenit(genofile **gfile, char *geno2name, SNP **snp2m, Indiv **indiv2m, int numsnp2,
   int numindiv2, int nignore)
    \brief Not used in EIGENSOFT (obsolete?)
 */

/*! \fn int genoreadit(genofile *gfile, SNP **pcupt)
    \brief Not used in EIGENSOFT (obsolete?)
 */

/*! \fn int putped(int num)
    \brief Store array snpord in snporda
    \param num   Index in snporda in which to store copy of array
 */

/*! \fn void getped(int num)
    \brief Copy array snpord from snporda
    \param num   Index in snporda from which to copy array
 */

/*! \fn int getweights(char *fname, SNP **snpm, int numsnps)
    \brief Read SNP weights from input file
    \param fname  Weight file name
    \param snpm   Array of SNP data
    \param numsnps  Number of SNPs in snpm
    \return Number of weights set
 */


void
setchr (int mode)
{
  chrmode = mode;
}

void
setchimpmode (int mode)
{
  chimpmode = mode;
}

void
ckdup (char **eglist, int n)
{
  if (checkdup (eglist, n)) {
    printdups (eglist, n);
    fatalx ("dup population found!\n");
  }
}
