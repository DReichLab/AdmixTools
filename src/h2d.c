#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <nicklib.h>
#include <admutils.h>

extern int verbose;
extern long rlen, packlen ; 
extern char *packgenos ; 

static char *pack2 ;   

static int pseudodip = YES ;  // for hap to dip.  

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
      if ((g1[j] >= 0) && (g2[j] >= 0)) {  
	tt = g1[j] + g2[j];
        if (pseudodip) tt /= 2 ; 
      }
      putgtypes (cupt, j, tt);
    }
  }

  free (g1);
  free (g2);
  free (x1);
  free (x2);

}

int mkindd2h (Indiv ** indivmarkers, Indiv *** pindm2, int numindivs)
{
  char ss[50], s1[50], s2[50] ;
  Indiv *indx, **indm2, *indp, *ind1, *ind2 ;
  int n, len, k;
  int numind2;

  numind2 = numindivs * 2;
  ZALLOC (*pindm2, numind2, Indiv *);
  indm2 = *pindm2;
  n = 0;
  for (k = 0; k < numindivs; k++) {
    indx = indivmarkers[k];
    strcpy (s1, indx->ID);
    strcat (s1, ":A") ; 
    strcpy (s2, indx->ID);
    strcat (s2, ":B") ; 
    if (strlen(s1) >=IDSIZE) fatalx("d2h overflow: %s\n", indx -> ID) ;
    ZALLOC (indm2[n], 1, Indiv);
    indp = indm2[n];
    *indp = *indx;
    strcpy (indp->ID, s1);
    ++n;
    ZALLOC (indm2[n], 1, Indiv);
    indp = indm2[n];
    *indp = *indx;
    strcpy (indp->ID, s2);
    ++n;
  }
  if (n != numind2) fatalx ("(mkindh2d) bug\n");
  return n;
}

void
remapd2h (SNP ** snpmarkers, int numsnps, Indiv ** indivmarkers,
	  Indiv ** indm2, int numindivs, int numind2)
{

  int g1, g2, gg;
  int *xx, *x2;
  int *tind, tt, t, i, j, k, j1, j2, x;
  Indiv *indx;
  SNP *cupt;
  char s1[50], s2[50];
  long rlen2, packlen2 ; 
  char *pbuff ;

  rlen2 = rlen  * 2 ; 
  packlen2 = packlen  * 2 ; 

  ZALLOC(pack2, packlen2, char) ;  
  pbuff = pack2 ; 

  ZALLOC (xx, numind2, int);
  ZALLOC (x2, numindivs, int);

  x = 0 ;
  for (i=0 ; i < numsnps; ++i) { 
   cupt = snpmarkers[i] ; 
   ivclear(xx, -1, numind2) ; 
   for (k=0 ; k < numindivs; ++k) { 
    gg = getgtypes(cupt, k) ; 
    if (gg<0) gg = 3 ; 
    g1 = g2 = gg ; 
    if (gg == 1) {  
     x = ranmod(2) ; 
     g1 = 2*x ; 
     g2 = 2*(1-x) ; 
    }
    xx[2*k] = g1 ; 
    xx[2*k+1] = g2 ; 
   }
   for (k=0; k<numind2; ++k) { 
    gg = xx[k] ; 
    wbuff( (unsigned char *) pbuff, k, gg) ; 
   }
   cupt -> pbuff = pbuff ; 
   pbuff += rlen2 ; 
  }

  free (xx);
  free (x2);
  free (packgenos) ; 
  packgenos = pack2 ; 
  rlen = rlen2 ; 
  packlen = packlen2 ;

}

