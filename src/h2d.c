#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <nicklib.h>
#include <admutils.h>

extern int verbose;

int
mkindh2d (Indiv ** indivmarkers, Indiv *** pindm2, int numindivs)
{
  char ss[50];
  Indiv *indx, **indm2, *indp;
  int n, len, k;
  int numind2;

  numind2 = numindivs / 2;
  ZALLOC (*pindm2, numind2, Indiv *);
  indm2 = *pindm2;
  n = 0;
  for (k = 0; k < numindivs; k++) {
    indx = indivmarkers[k];
    strcpy (ss, indx->ID);
    len = strlen (ss);
    if (ss[len - 1] != 'A')
      continue;
    ss[len - 2] = CNULL;
    ZALLOC (indm2[n], 1, Indiv);
    indp = indm2[n];
    *indp = *indx;
    strcpy (indp->ID, ss);
    ++n;
  }
  if (n != numind2)
    fatalx ("(mkindh2d) bug\n");
  return n;
}

void
remaph2d (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	  Indiv ** indm2, int numindivs, int numind2)
{

  int *g1, *g2;
  int *x1, *x2;
  int *tind, tt, t, i, j, k, j1, j2;
  Indiv *indx;
  SNP *cupt;
  char s1[50], s2[50];

  ZALLOC (g2, numind2, int);
  ZALLOC (g1, numindivs, int);
  ZALLOC (x1, numindivs, int);
  ZALLOC (x2, numindivs, int);

  for (k = 0; k < numind2; ++k) {
    indx = indm2[k];
    sprintf (s1, "%s:A", indx->ID);
    sprintf (s2, "%s:B", indx->ID);
    t = x1[k] = indindex (indivmarkers, numindivs, s1);
    if (t < 0) {
      sprintf (s1, "%s_A", indx->ID);
      sprintf (s2, "%s_B", indx->ID);
      t = x1[k] = indindex (indivmarkers, numindivs, s1);
    }
    if (t < 0)
      fatalx ("bad newindiv: %s\n", indx->ID);
    t = x2[k] = indindex (indivmarkers, numindivs, s2);
    if (t < 0)
      fatalx ("bad newindiv: %s\n", indx->ID);
  }

  for (i = 0; i < numsnps; i++) {
    cupt = snpmarkers[i];

    for (j = 0; j < numind2; ++j) {
      t = x1[j];
      g1[j] = getgtypes (cupt, t);
      t = x2[j];
      g2[j] = getgtypes (cupt, t);
      tt = -1;
      if ((g1[j] >= 0) && (g2[j] >= 0))
	tt = g1[j] + g2[j];
      putgtypes (cupt, j, tt);
    }
  }

  free (g1);
  free (g2);
  free (x1);
  free (x2);

}
