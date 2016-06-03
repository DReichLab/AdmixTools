#include <exclude.h>
#include <math.h>
#include <mcio.h>

#define MAXRGN 1000

void
excluderegions (char *xregionname, SNP ** snps, int nsnps,
		char *deletesnpoutname)
{
  FILE *fp;

  int chr[MAXRGN];
  int lo[MAXRGN];
  int hi[MAXRGN];

  char line[MAXSTR];
  char *spt[MAXFF];
  int nsplit, nrgn, i, j;

  if ((fp = fopen (xregionname, "r")) == NULL) {
    printf ("excluderegions: can't open file %s\n", xregionname);
    return;
  }

  for (i = 0; i < MAXRGN; i++) {

    if (fgets (line, MAXSTR, fp) == NULL)
      break;

    nsplit = splitup (line, spt, MAXFF);
    if (nsplit != 3)
      continue;

    chr[i] = atoi (spt[0]);
    lo[i] = atoi (spt[1]);
    hi[i] = atoi (spt[2]);

  }
  fclose (fp);
  nrgn = i;

  for (i = 0; i < nsnps; i++) {
    SNP *cupt = snps[i];
    for (j = 0; j < nrgn; j++) {
      if (cupt->chrom == chr[j] && cupt->physpos >= lo[j]
	  && cupt->physpos <= hi[j]) {
	cupt->ignore = YES;
	if (deletesnpoutname != NULL) {
	  logdeletedsnp (cupt->ID, "xregion", deletesnpoutname);
	}
      }
    }
  }

  return;

}

void
hwfilter (SNP ** snps, int nsnps, int nindiv, double nhwfilter,
	  char *deletesnpoutname)
{

  int i, k;

  for (i = 0; i < nsnps; i++) {
    int num = 0, den = 0, het = 0, n0 = 0, n1 = 0, n2 = 0, nsamples;
    double p, Q, stdv;
    SNP *cupt = snps[i];

    for (k = 0; k < nindiv; k++) {
      int g = getgtypes (cupt, k);
      if (g >= 0) {
	num += g;
	den += 2;
      }
      if (g == 1) {
	het++;
	n1++;
      }
      else if (g == 0) {
	n0++;
      }
      else if (g == 2) {
	n2++;
      }
    }

    if ((nsamples = den / 2) == 0)
      continue;
    p = (double) num / den;
    Q = 2 * p * (1 - p);
    stdv = sqrt (Q * (1 - Q) / nsamples);
    if (fabs ((double) het / nsamples - Q) > nhwfilter * stdv) {
      printf ("SNP %s removed by Hardy-Weinberg filter\n", cupt->ID);
      cupt->ignore = YES;
      if (deletesnpoutname != NULL) {
	logdeletedsnp (cupt->ID, "hwfilt", deletesnpoutname);
      }
    }
  }


}
