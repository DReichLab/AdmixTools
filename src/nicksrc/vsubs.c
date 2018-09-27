#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include "strsubs.h"
#include "vsubs.h"

/** 
 tiny routines BLAS? 
 a small library to do simple arithmetic
 on 1D vectors with no skips 
*/
void
vsp (double *a, double *b, double c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] + c;
}

void
vst (double *a, double *b, double c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] * c;
}

void
vvt (double *a, double *b, double *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] * c[i];
}

void
vvp (double *a, double *b, double *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] + c[i];
}

void
vvm (double *a, double *b, double *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] - c[i];
}

void
vvd (double *a, double *b, double *c, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    if (c[i] == 0.0)
      fatalx ("(vvd): zero value in denominator\n");
    a[i] = b[i] / c[i];
  }
}

void
vsqrt (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    if (b[i] < 0.0)
      fatalx ("(vsqrt): negative value %g\n", b[i]);
    if (b[i] == 0.0) {
      a[i] = 0.0;
      continue;
    }
    a[i] = sqrt (b[i]);
  }
}

void
vinvert (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    if (b[i] == 0.0)
      fatalx ("(vinvert): zero value\n");
    a[i] = 1.0 / b[i];
  }
}

void
vabs (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    a[i] = fabs (b[i]);
  }
}

void
vlog (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    if (b[i] <= 0.0)
      fatalx ("(vlog): negative or zero value %g\n", b[i]);
    a[i] = log (b[i]);
  }
}

void
vlog2 (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    if (b[i] <= 0.0)
      fatalx ("(vlog2): negative or zero value %g\n", b[i]);
    a[i] = NPlog2 (b[i]);
  }
}

void
vexp (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    a[i] = exp (b[i]);
  }
}

void
vclear (double *a, double c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = c;
}

void
vzero (double *a, int n)
{
  vclear (a, 0.0, n);
}

void
cpzero (char **a, int n)
{
  int i;
  for (i = 0; i < n; ++i) {
    a[i] = NULL;
  }
}

void
cclear (unsigned char *a, unsigned char c, long n)

/** 
 be careful nothing done about NULL at end
*/
{
  long i;
  for (i = 0; i < n; i++) {
    a[i] = c;
  }
}

void
charclear (char *a, unsigned char c, long n)
// fussy compiler warnigns about unsigned char conversions avoided
{

 cclear( (unsigned char *) a, c, n) ; 

}

void
ivvp (int *a, int *b, int *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] + c[i];
}

void
ivvm (int *a, int *b, int *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] - c[i];
}

void
ivsp (int *a, int *b, int c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] + c;
}

void
ivst (int *a, int *b, int c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] * c;
}

void
ivclear (int *a, int c, long n)
{
  long i;
  for (i = 0; i < n; i++)
    a[i] = c;
}

void
lvclear (long *a, long c, long n)
{
  long i;
  for (i = 0; i < n; i++)
    a[i] = c;
}

void
ivzero (int *a, int n)
{
  ivclear (a, 0, n);
}

void
lvzero (long *a, long n)
{
  lvclear (a, 0, n);
}

double
clip (double x, double lo, double hi)

/* clip off values to range [lo,hi] */
{
  if (x < lo)
    return lo;
  if (x > hi)
    return hi;
  return x;
}

void
ivclip (int *a, int *b, int loval, int hival, int n)
{

/* clip off values to range [loval,hival] */
  int i;
  int t;

  for (i = 0; i < n; i++) {
    t = MAX (b[i], loval);
    a[i] = MIN (t, hival);
  }
}

void
vclip (double *a, double *b, double loval, double hival, int n)
{

/* clip off values to range [loval,hival] */
  int i;
  double t;

  for (i = 0; i < n; i++) {
    t = MAX (b[i], loval);
    a[i] = MIN (t, hival);
  }
}

void
vmaxmin (double *a, int n, double *max, double *min)
{

  int i;
  double tmax, tmin;

  tmax = tmin = a[0];
  for (i = 1; i < n; i++) {
    tmax = MAX (tmax, a[i]);
    tmin = MIN (tmin, a[i]);
  }
  if (max != NULL)
    *max = tmax;
  if (min != NULL)
    *min = tmin;
}

void
vlmaxmin (double *a, int n, int *pmax, int *pmin)

/** 
 return location 
*/
{

  int i;
  double tmax, tmin;
  double lmax, lmin;

  tmax = tmin = a[0];
  lmax = lmin = 0;
  for (i = 1; i < n; i++) {
    if (a[i] > tmax) {
      tmax = a[i];
      lmax = i;
    }
    if (a[i] < tmin) {
      tmin = a[i];
      lmin = i;
    }
  }
  if (pmax != NULL)
    *pmax = lmax;
  if (pmin != NULL)
    *pmin = lmin;
}

void
ivmaxmin (int *a, int n, int *max, int *min)
{

  int i;
  int tmax, tmin;

  tmax = tmin = a[0];
  for (i = 1; i < n; i++) {
    tmax = MAX (tmax, a[i]);
    tmin = MIN (tmin, a[i]);
  }
  if (max != NULL)
    *max = tmax;
  if (min != NULL)
    *min = tmin;
}

int
minivec (int *a, int n)
{
  int t;

  ivmaxmin (a, n, NULL, &t);
  return t;

}

int
maxivec (int *a, int n)
{
  int t;

  ivmaxmin (a, n, &t, NULL);
  return t;

}

void
ivlmaxmin (int *a, int n, int *pmax, int *pmin)

/** 
 return location 
*/
{

  int i;
  int tmax, tmin;
  int lmax, lmin;

  tmax = tmin = a[0];
  lmax = lmin = 0;
  for (i = 1; i < n; i++) {
    if (a[i] > tmax) {
      tmax = a[i];
      lmax = i;
    }
    if (a[i] < tmin) {
      tmin = a[i];
      lmin = i;
    }
  }
  if (pmax != NULL)
    *pmax = lmax;
  if (pmin != NULL)
    *pmin = lmin;
}

double
vdot (double *a, double *b, int n)
{
  int i;
  double ans = 0.0;
  for (i = 0; i < n; i++)
    ans += a[i] * b[i];

  return ans;

}

double
corr (double *a, double *b, int n)
{
  double v12, v11, v22, y1, y2, y;
  double *aa, *bb;
  ZALLOC (aa, n, double);
  ZALLOC (bb, n, double);
  y1 = asum (a, n) / (double) n;
  y2 = asum (b, n) / (double) n;

  vsp (aa, a, -y1, n);
  vsp (bb, b, -y2, n);

  v12 = vdot (aa, bb, n);
  v11 = asum2 (aa, n);
  v22 = asum2 (bb, n);

  y = v11 * v22;
  if (y == 0.0)
    fatalx ("(corr) constant vector\n");


  free (aa);
  free (bb);
  return (v12 / sqrt (y));

}

double
corrx (double *a, double *b, int n)
// like corr but constant vec returns 0
{
  double v12, v11, v22, y1, y2, y;
  double *aa, *bb;

  ZALLOC (aa, n, double);
  ZALLOC (bb, n, double);
  y1 = asum (a, n) / (double) n;
  y2 = asum (b, n) / (double) n;

  vsp (aa, a, -y1, n);
  vsp (bb, b, -y2, n);

  v12 = vdot (aa, bb, n);
  v11 = asum2 (aa, n);
  v22 = asum2 (bb, n);

  free (aa);
  free (bb);

  y = v11 * v22;
  y += 1.0e-12;

  return (v12 / sqrt (y));

}


double
variance (double *a, int n)
{

  double *aa;
  double y1, y2;

  ZALLOC (aa, n, double);
  y1 = asum (a, n) / (double) n;
  vsp (aa, a, -y1, n);

  y2 = asum (aa, n) / (double) n;

  free (aa);
  return y2;

}

void
getdiag (double *a, double *b, int n)

/* extract diagonal */
{
  int i, k;

  for (i = 0; i < n; i++) {
    k = i * n + i;
    a[i] = b[k];
  }
}

void
setdiag (double *a, double *diag, int n)

/* set diagonal matrix */
{
  int i, k;

  vzero (a, n * n);
  for (i = 0; i < n; i++) {
    k = i * n + i;
    a[k] = diag[i];
  }
}

void
adddiag (double *a, double *diag, int n)

/* add diagonal matrix */
{
  int i, k;

  for (i = 0; i < n; i++) {
    k = i * n + i;
    a[k] += diag[i];
  }
}



void
copyarr (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void
revarr (double *b, double *a, int n)
{
  int i;
  double *x;
  ZALLOC (x, n, double);
  for (i = 0; i < n; i++) {
    x[n - i - 1] = a[i];
  }
  copyarr (x, b, n);
  free (x);
}

void
revlarr (long *b, long *a, int n)
{
  int i;
  long *x;
  ZALLOC (x, n, long);
  for (i = 0; i < n; i++) {
    x[n - i - 1] = a[i];
  }
  copylarr (x, b, n);
  free (x);
}


void
revuiarr (unsigned int *b, unsigned int *a, int n)
{
  int i;
  unsigned int *x;
  ZALLOC (x, n, unsigned int);
  for (i = 0; i < n; i++) {
    x[n - i - 1] = a[i];
  }
  for (i = 0; i < n; i++) {
    b[i] = x[i];
  }
  free (x);
}

void
reviarr (int *b, int *a, int n)
{
  int i;
  int *x;
  ZALLOC (x, n, int);
  for (i = 0; i < n; i++) {
    x[n - i - 1] = a[i];
  }
  copyiarr (x, b, n);
  free (x);

}

void
copyiarr (int *a, int *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void
copylarr (long *a, long *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void
copyiparr (int **a, int **b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void
dpermute (double *a, int *ind, int len)
{

  int i, k;
  double *rrr;

  ZALLOC (rrr, len, double);

  for (i = 0; i < len; i++) {
    rrr[i] = a[i];
  }

  for (i = 0; i < len; i++) {
    k = ind[i];
    a[i] = rrr[k];
  }

  free (rrr);
}

void
ipermute (int *a, int *ind, int len)
{

  int i, k;
  int *rrr;

  ZALLOC (rrr, len, int);

  copyiarr (a, rrr, len);

  for (i = 0; i < len; i++) {
    k = ind[i];
    a[i] = rrr[k];
  }

  free (rrr);
}

void
dppermute (double **a, int *ind, int len)
{

  int i, k;
  double **rrr;

  ZALLOC (rrr, len, double *);

  for (i = 0; i < len; i++) {
    rrr[i] = a[i];
  }

  for (i = 0; i < len; i++) {
    k = ind[i];
    a[i] = rrr[k];
  }

  free (rrr);
}

void
ippermute (int **a, int *ind, int len)
{

  int i, k;
  int **rrr;

  ZALLOC (rrr, len, int *);

  for (i = 0; i < len; i++) {
    rrr[i] = a[i];
  }

  for (i = 0; i < len; i++) {
    k = ind[i];
    a[i] = rrr[k];
  }

  free (rrr);
}

double
asum (double *a, int n)
{
  int i;
  double ans = 0.0;
  for (i = 0; i < n; i++)
    ans += a[i];

  return ans;
}

int
intsum (int *a, int n)
{
  int i;
  int ans = 0;
  for (i = 0; i < n; i++)
    ans += a[i];

  return ans;
}

long
longsum (long *a, int n)
{
  int i;
  long ans = 0;
  for (i = 0; i < n; i++)
    ans += a[i];

  return ans;
}

int
idot (int *a, int *b, int n)
{
  int i;
  int ans = 0.0;
  for (i = 0; i < n; i++)
    ans += a[i] * b[i];

  return ans;

}

int
iprod (int *a, int n)

/* overflow not checked */
{
  int i;
  int ans = 1;
  for (i = 0; i < n; i++)
    ans *= a[i];

  return ans;
}


double
aprod (double *a, int n)

/* overflow not checked */
{
  int i;
  double ans = 1.0;
  for (i = 0; i < n; i++)
    ans *= a[i];

  return ans;
}

double
asum2 (double *a, int n)
{
  int i;
  double ans = 0.0;
  for (i = 0; i < n; i++)
    ans += a[i] * a[i];

  return ans;
}

double
trace (double *a, int n)
{
  double *diags, t;
  ZALLOC (diags, n, double);
  getdiag (diags, a, n);        /* extract diagonal */
  t = asum (diags, n);
  free (diags);
  return t;
}

int
nnint (double x)
{
  long int lrint (double x);
// double round(double x) ;
  return (int) lrint (x);
}

void
countcat (int *tags, int n, int *ncat, int nclass)

/* simple frequency count of integer array */
{
  int i, k;
  ivzero (ncat, nclass);
  for (i = 0; i < n; i++) {
    k = tags[i];
    if ((k < 0) || (k >= nclass))
      fatalx ("(countcat) bounds error\n");
    ++ncat[k];
  }
}

void
rowsum (double *a, double *rr, int n)
// square matrix 
{
  int i, j;
  vclear (rr, 0.0, n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      rr[j] += a[i + j * n];
    }
  }
}

void
colsum (double *a, double *cc, int n)
// square matrix 
{
  int i, j;
  vclear (cc, 0.0, n);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      cc[i] += a[i + j * n];
    }
  }
}

void
rrsum (double *a, double *rr, int m, int n)
{
  int i, j;
  vclear (rr, 0.0, n);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      rr[j] += a[i + j * m];
    }
  }
}

void
ccsum (double *a, double *cc, int m, int n)
{
  int i, j;
  vclear (cc, 0.0, m);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      cc[i] += a[i + j * m];
    }
  }
}

void
printmatfile (double *a, int m, int n, FILE * fff)

/** 
 print a matrix n wide m rows  
*/
{
  printmatwfile (a, m, n, 5, fff);
}

void
printmatwfile (double *a, int m, int n, int w, FILE * fff)

/** 
 print a matrix n wide m rows  w to a row
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      fprintf (fff, "%9.3f ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        fprintf (fff, "  ...\n");
      }
    }
    fprintf (fff, "\n");
  }
}

void
printmatx (double *a, int m, int n)

/** 
 print a matrix n wide m rows   no final nl
*/
{
  printmatwx (a, m, n, 5);
}

void
printmat (double *a, int m, int n)

/** 
 print a matrix n wide m rows  
*/
{
  printmatw (a, m, n, 5);
}

void
printmatwx (double *a, int m, int n, int w)

/** 
 print a matrix n wide m rows  w to a row
 no final nl
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%9.3f ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    if (i < (m - 1))
      printf ("\n");
  }
}

void
printmatw (double *a, int m, int n, int w)

/** 
 print a matrix n wide m rows  w to a row
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%9.3f ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printmatl (double *a, int m, int n)

/** 
 print a matrix n wide m rows  
*/
{
  printmatwl (a, m, n, 5);
}

void
printmatlx (double *a, int m, int n)

/** 
 print a matrix n wide m rows  
 no final \n
*/
{
  printmatwlx (a, m, n, 5);
}

void
printmatwl (double *a, int m, int n, int w)

/** 
 print a matrix n wide m rows  w to a row
 15.9f format
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%15.9f ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printmatwlx (double *a, int m, int n, int w)

/** 
 print a matrix n wide m rows  w to a row
 15.9f format
 No final \n
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%15.9f ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
  }
}

void
printmatwf (double *a, int m, int n, int w, char *format)

/**
 print a matrix n wide m rows  w to a row with format
 no spacing introduced here.  User must supply
*/
{
  int i, j, jmod;
  if (format == NULL) {
    printmatw (a, m, n, w);
    return;
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf (format, a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printmat2D (double **a, int m, int n)
{
  int k;
  for (k = 0; k < m; ++k) {
    printf ("%3d: ", k);
    printmat (a[k], 1, n);
  }
}

void
int2c (char *cc, int *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    cc[i] = (char) b[i];
  }
  cc[n] = '\0';
}

void
floatit (double *a, int *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    a[i] = (double) b[i];
  }
}

void
printimatwfile (int *a, int m, int n, int w, FILE * fff)

/** 
 print a matrix n wide m rows  w to a row
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      fprintf (fff, "%5d ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        fprintf (fff, "  ...\n");
      }
    }
    fprintf (fff, "\n");
  }
}

void
printimatw (int *a, int m, int n, int w)

/** 
 print a matrix n wide m rows  w to a row
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%5d ", a[i * n + j]);
      jmod = (j + 1) % w;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printimatx (int *a, int m, int n)

/** 
 print a matrix n wide m rows  
 no final new line
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%5d ", a[i * n + j]);
      jmod = (j + 1) % 10;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
  }
}

void
printimatfile (int *a, int m, int n, FILE * fff)

/** 
 print a matrix n wide m rows  
*/
{
  printimatwfile (a, m, n, 10, fff);
}

void
printimat2D (int **a, int m, int n)
{
  int k;

  for (k = 0; k < m; ++k) {
    printimat (a[k], 1, n);
  }
}

void
printimat1 (int *a, int m, int n)

/** 
 print a matrix n wide m rows, %1d format  
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%1d", a[i * n + j]);
      jmod = (j + 1) % 50;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printimat (int *a, int m, int n)

/** 
 print a matrix n wide m rows  
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%5d ", a[i * n + j]);
      jmod = (j + 1) % 10;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printimatlfile (int *a, int m, int n, FILE * fff)

/** 
 print a matrix n wide m rows  %10d format
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      fprintf (fff, "%10d ", a[i * n + j]);
      jmod = (j + 1) % 10;
      if ((jmod == 0) && (j < (n - 1))) {
        fprintf (fff, "  ...\n");
      }
    }
    fprintf (fff, "\n");
  }
}

void
printimatl (int *a, int m, int n)

/** 
 print a matrix n wide m rows  %10d format
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%10d ", a[i * n + j]);
      jmod = (j + 1) % 10;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}

void
printimatlx (int *a, int m, int n)

/** 
 print a matrix n wide m rows  %10d format
 no final newline  
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%10d ", a[i * n + j]);
      jmod = (j + 1) % 10;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    if (i < (m - 1))
      printf ("\n");
  }
}

void
printstringf (char *ss, int w, FILE * fff)
{
  char *sss;
  char *sx;

  ZALLOC (sss, w + 1, char);
  cclear ((unsigned char *) sss, CNULL, w + 1);

  sx = ss;
  for (;;) {
    strncpy (sss, sx, w);
    if (strlen (sss) <= 0)
      break;
    sx += w;
    fprintf (fff, "%s\n", sss);
  }

  free (sss);
}


void
printstringbasepos (char *ss, int w, int basepos)
{
  char *sss;
  char *sx;
  int pos = basepos;

  ZALLOC (sss, w + 1, char);
  cclear ((unsigned char *) sss, CNULL, w + 1);

  sx = ss;
  for (;;) {
    strncpy (sss, sx, w);
    if (strlen (sss) <= 0)
      break;
    printf ("%12d ", pos);
    printf ("%s\n", sss);
    sx += w;
    pos += w;
  }

  free (sss);
}



void
printstring (char *ss, int w)
{
  printstringf (ss, w, stdout);
}


void
rndit (double *a, double *b, int n)
{
  int i;

  for (i = 0; i < n; ++i) {
    a[i] = nearbyint (b[i]);
  }
}


void
fixit (int *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    a[i] = nnint (b[i]);
  }
}
void
fixitl (long *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    a[i] = lrint (b[i]);
  }
}


int
findfirst (int *a, int n, int val)
{
  int i;
  for (i = 0; i < n; i++) {
    if (a[i] == val)
      return i;
  }
  return -1;
}

int
findfirstl (long *a, int n, long val)
{
  int i;
  for (i = 0; i < n; i++) {
    if (a[i] == val)
      return i;
  }
  return -1;
}

int
findfirstu (unsigned int *a, int n, unsigned int val)
{
  int i;
  for (i = 0; i < n; i++) {
    if (a[i] == val)
      return i;
  }
  return -1;
}

int
findlastu (unsigned int *a, int n, unsigned int val)
{
  int i;
  for (i = n - 1; i >= 0; i--) {
    if (a[i] == val)
      return i;
  }
  return -1;
}

int
findlast (int *a, int n, int val)
{
  int i;
  for (i = n - 1; i >= 0; i--) {
    if (a[i] == val)
      return i;
  }
  return -1;
}

int
binsearch (int *a, int n, int val)
// binary search.  a sorted in ascending order
{
#define TINYS 12
  int x, m, h, v;

  if (n <= TINYS)
    return findfirst (a, n, val);
  if (val < a[0])
    return -1;
  if (val > a[n - 1])
    return -1;
  h = n / 2;
  v = a[h];
  if (val < v)
    return binsearch (a, h, val);
  if (val == v)
    return h;
  m = (n - 1) - (h + 1) + 1;
  x = binsearch (a + h + 1, m, val);
  if (x < 0)
    return -1;
  return x + h + 1;
}

void
idperm (int *a, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = i;
}

double
NPlog2 (double y)
{
  if (y <= 0.0)
    fatalx ("(NPlog2) negative argument\n");
  return (log (y) / log (2.0));
}

double
logfac (int n)

/** 
 log (factorial n))
*/
{
  double y, x;
  x = (double) (n + 1);
  y = lgamma (x);
  return (y);
}

double
logbino (int n, int k)

/* log n choose k */
{
  double top, bot;

  top = logfac (n);
  bot = logfac (n - k) + logfac (k);

  return top - bot;
}

double
loghprob (int n, int a, int m, int k)
// http://www.math.uah.edu/stat/urn/Hypergeometric.xhtml
{

/** 
 n balls a black.  Pick m without replacement  
 return log prob (k black)
*/

  double ytop, ybot;

  if (k < 0)
    return -1.0e30;
  if (k > a)
    return -1.0e30;
  if (k > m)
    return -1.0e30;
  if ((m - k) > (n - a))
    return -1.0e30;

  ytop = logbino (a, k) + logbino (n - a, m - k);
  ybot = logbino (n, m);
  return ytop - ybot;

}

int hprobv(double *vprob, int n, int a, int m)
//  vprob[k] is hprob(n,a, m, k) ;  vprob preallocated m+1  
{
   int k, mp = m+1 ;
   double y ;  
 
   vclear(vprob, -1.0e40, mp) ;
   for (k=0; k<=m; ++k)    {
     vprob[k] = loghprob(n, a, m, k) ;
   }
  
   vexp(vprob, vprob, mp) ;
   y = bal1(vprob, mp) ;
   if (fabs(y-1.0) > 0.01) return -1 ;
/**
   printf("sum:: %15.9f\n", y) ;
   printmatl(vprob, 1, mp) ;
*/
   return 1 ;
}          

double
log2fac (int n)

/** 
 log base2 (factorial n))
*/
{
  double y, x;
  x = (double) (n + 1);
  y = lgamma (x);
  return (y / log (2.0));
}

double
addlog (double a, double b)
{
  /* given a = log(A)
     b = log(B)
     returns log(A+B) 
     with precautions for overflow etc
   */
  double x, y, z;

  x = MIN (a, b);
  y = MAX (a, b);

/** 
 answer is log(1 + A/B) + log (B)  
*/
  z = x - y;
  if (z < -50.0)
    return y;
  z = 1.0 + exp (z);
  z = log (z) + y;
  return (z);

}

double logsum(double *x, int n) 
{
// no test for 0
  
  double *w, y ; 
  ZALLOC(w, n, double) ; 

  vlog(w, x, n) ; 
  y = asum(w, n) ;

  free(w) ; 
  return y ; 


} 


double
vldot (double *x, double *y, int n)

/** 
 x. log(y) 
*/
{
  double *z, ans;
  double tiny = 1.0e-19;
  int i;

  ZALLOC (z, n, double);
  vsp (z, y, 1.0e-20, n);
  vlog (z, z, n);

  ans = 0.0;
  for (i = 0; i < n; i++) {
    if (x[i] > tiny)
      ans += x[i] * z[i];
  }
  free (z);
  return ans;
}

int
ipow2 (int l)
{
  return nnint (pow (2.0, l));
}

double
pow10 (double x)
{
  return exp (x * log (10.0));
}


void
vpow10 (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = exp (b[i] * log (10.0));
}

void
vlog10 (double *a, double *b, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = log10 (b[i]);
}

void
transpose (double *aout, double *ain, int m, int n)

/** 
 aout and ain must be identical or not overlap 
 does matrix transpose 

 input  m vectors of length n  (m x n) 
 output n vectors of length m 
*/
{
  double *ttt;
  int i, j, k1, k2;
  if (aout == ain) {
    ZALLOC (ttt, m * n, double);
  }
  else
    ttt = aout;

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) {
      k1 = i * n + j;
      k2 = j * m + i;
      ttt[k2] = ain[k1];
    }
  if (aout == ain) {
    copyarr (ttt, aout, m * n);
    free (ttt);
  }
}

int **
initarray_2Dint (int numrows, int numcolumns, int initval)
{
  int i, j;
  int **array;

  ZALLOC (array, numrows, int *);
  for (i = 0; i < numrows; i++) {
    ZALLOC (array[i], numcolumns, int);
    if (initval != 0)
      ivclear (array[i], initval, numcolumns);
  }
  return array;
}

long **
initarray_2Dlong (int numrows, int numcolumns, long initval)
{
  int i, j;
  long **array;


  ZALLOC (array, numrows, long *);
  for (i = 0; i < numrows; i++) {
    ZALLOC (array[i], numcolumns, long);
    if (initval != 0)
      lvclear (array[i], initval, numcolumns);
  }
  return array;
}


void
free2Dint (int ***xx, int numrows)
{
  int **array;
  int i;
  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    free (array[i]);
  }
  free (array);
  *xx = NULL;
}

void
free2Dlong (long ***xx, int numrows)
{
  long **array;
  int i;
  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    free (array[i]);
  }
  free (array);
  *xx = NULL;
}

void
free_darray (double **xx)
{
  free (*xx);
  *xx = NULL;
}

void
free_iarray (int **xx)
{
  free (*xx);
  *xx = NULL;
}


double **
initarray_2Ddouble (int numrows, int numcolumns, double initval)
{
  int i, j;
  double **array;


  ZALLOC (array, numrows, double *);
  for (i = 0; i < numrows; i++) {
    ZALLOC (array[i], numcolumns, double);
    if (initval != 0.0)
      vclear (array[i], initval, numcolumns);
  }
  return array;
}

long double **
initarray_2Dlongdouble (int numrows, int numcolumns, long double initval)
{
  int i, j;
  long double **array, *bb;


  ZALLOC (array, numrows, long double *);
  for (i = 0; i < numrows; i++) {
    ZALLOC (array[i], numcolumns, long double);
    if (initval != 0.0) {
      bb = array[i];
      for (j = 0; j < numcolumns; ++j) {
        bb[j] = initval;
      }
    }
  }
  return array;
}


void
clear2D (double ***xx, int numrows, int numcols, double val)
{
  double **array;
  int i;
  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    vclear (array[i], val, numcols);
  }

}

void
iclear2D (int ***xx, int numrows, int numcols, int val)
{
  int **array;
  int i;

  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    ivclear (array[i], val, numcols);
  }

}

void
lclear2D (long ***xx, int numrows, int numcols, long val)
{
  long **array;
  int i;

  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    lvclear (array[i], val, numcols);
  }

}

void
free2D (double ***xx, int numrows)
{
  double **array;
  int i;
  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    free (array[i]);
  }
  free (array);
  *xx = NULL;
}

void
free2Dlongdouble (long double ***xx, int numrows)
{
  long double **array;
  int i;

  array = *xx;

  for (i = numrows - 1; i >= 0; i--) {
    free (array[i]);
  }
  free (array);
  *xx = NULL;
}

void
addoutmul (double *mat, double *v, double mul, int n)
{
  int a, b;
  for (a = 0; a < n; ++a) {
    for (b = 0; b < n; ++b) {
      mat[a * n + b] += v[a] * v[b] * mul;
    }
  }
}




void
addouter (double *out, double *a, int n)

/* 
 add outerprod(a)  to out
 trival to recode to make ~ 2 * faster
*/
{

  addoutmul (out, a, 1.0, n);

}

void
subouter (double *out, double *a, int n)

/* 
 subtract outerprod(a)  to out
 trival to recode to make ~ 2 * faster
*/
{
  addoutmul (out, a, -1.0, n);

}

double
bal1 (double *a, int n)
// WARNING a is input and output
{
  double y;

  y = asum (a, n);
  if (y <= 0.0)
    fatalx ("bad bal1\n");
  vst (a, a, 1.0 / y, n);
  return y;

}

double
bal2 (double *a, int n)
// WARNING a is input and output
{
  double y;

  y = asum2 (a, n);

  if (y <= 0.0)
    fatalx ("bad bal2\n");

  y = sqrt(y) ; 
  vst (a, a, 1.0 / y, n);

  return y;
}

double
logmultinom (int *cc, int n)

/* log multinomial */
{
  int t, k, i;
  double y, ytot;

  if (n <= 1)
    return 0.0;
  t = intsum (cc, n);
  if (t == 0)
    return 0.0;
  ytot = 0;
  for (i = 0; i < n - 1; i++) {
    k = cc[i];
    y = logbino (t, k);
    ytot += y;
    t -= k;
  }
  return ytot;
}

void
flipiarr (int *a, int *b, int n)
// reverse array 
{
  int *x, k;
  ZALLOC (x, n, int);

  for (k = 0; k < n; ++k) {
    x[n - 1 - k] = b[k];
  }

  copyiarr (x, a, n);

  free (x);


}

void
fliparr (double *a, double *b, int n)
{
  double *x;
  int k;

  ZALLOC (x, n, double);

  for (k = 0; k < n; ++k) {
    x[n - 1 - k] = b[k];
  }

  copyarr (x, a, n);

  free (x);

}

void
vcompl (double *a, double *b, int n)
// a <- 1 - b 
{
  double *x;
  ZALLOC (x, n, double);

  vvm (x, x, b, n);
  vsp (x, x, 1.0, n);

  copyarr (x, a, n);

  free (x);
}

void
setidmat (double *a, int n)
// a <- identity matrix
{
  int i;
  vzero (a, n * n);
  for (i = 0; i < n; i++) {
    a[i * n + i] = 1.0;
  }
}


int
stripit (double *a, double *b, int *x, int len)
// copy b to a leave out elems where x < 0
{
  int k, n;

  n = 0;
  for (k = 0; k < len; ++k) {
    if (x[k] >= 0) {
      a[n] = b[k];
      ++n;
    }
  }
  return n;
}

int
istripit (int *a, int *b, int *x, int len)
// copy b to a leave out elems where x < 0
{
  int k, n;

  n = 0;
  for (k = 0; k < len; ++k) {
    if (x[k] >= 0) {
      a[n] = b[k];
      ++n;
    }
  }
  return n;

}

int
cstripit (char **a, char **b, int *x, int len)
// copy b to a leave out elems where x < 0
{
  int k, n;

  n = 0;
  for (k = 0; k < len; ++k) {
    if (x[k] >= 0) {
      a[n] = b[k];
      ++n;
    }
  }
  return n;
}

void
mapit (int *a, int *b, int n, int inval, int outval)
{
  int k;

  copyiarr (b, a, n);
  for (k = 0; k < n; ++k) {
    if (a[k] == inval)
      a[k] = outval;
  }
}

int
ifall (int n, int k)
// falling factorial
{

  int prod = 1, t = n, j;

  for (j = 0; j < k; ++j) {
    prod *= t;
    --t;
  }
  return prod;
}

double
hlife (double val)
{

  return -log (2.0) / log (val);

}

void *
topheap ()
// find top of heap (address).  Useful for finding memory leaks 
{

  return sbrk (0);
}

void
swap (double *pa, double *pb)
{
  double a, b, t;

  a = *pa;
  b = *pb;
  t = b;
  b = a;
  a = t;
  *pa = a;
  *pb = b;

}

void
iswap (int *pa, int *pb)
{
  int a, b, t;

  a = *pa;
  b = *pb;
  t = b;
  b = a;
  a = t;
  *pa = a;
  *pb = b;

}

void
cswap (char *c1, char *c2)
{
  char cc;

  cc = *c1;
  *c1 = *c2;
  *c2 = cc;


}

int
kodeitb (int *xx, int len, int base)
{
  int t = 0, i;

  for (i = 0; i < len; ++i) {
    t *= base;
    t += xx[i];
  }
  return t;

}

int
dekodeitb (int *xx, int kode, int len, int base)
{

  int i, t;

  t = kode;
  for (i = len - 1; i >= 0; --i) {
    xx[i] = t % base;
    t /= base;
  }
  return intsum (xx, len);      // weight

}

void
floatit2D (double **a, int **b, int nrows, int ncols)
{

  int x;

  for (x = 0; x < nrows; ++x) {
    floatit (a[x], b[x], ncols);
  }
}

void
copyarr2D (double **a, double **b, int nrows, int ncols)
{

  int x;

  for (x = 0; x < nrows; ++x) {
    copyarr (a[x], b[x], ncols);
  }
}

void
copyiarr2D (int **a, int **b, int nrows, int ncols)
{

  int x;

  for (x = 0; x < nrows; ++x) {
    copyiarr (a[x], b[x], ncols);
  }
}


// be very careful about rows and columns here. 
void
plus2Dint (int **a, int **b, int **c, int nrows, int ncols)
{
  int x;

  for (x = 0; x < nrows; ++x) {
    ivvp (a[x], b[x], c[x], ncols);
  }
}

void
minus2Dint (int **a, int **b, int **c, int nrows, int ncols)
{
  int x;

  for (x = 0; x < nrows; ++x) {
    ivvm (a[x], b[x], c[x], ncols);
  }
}

void
plus2D (double **a, double **b, double **c, int nrows, int ncols)
{
  int x;

  for (x = 0; x < nrows; ++x) {
    vvp (a[x], b[x], c[x], ncols);
  }
}

void
minus2D (double **a, double **b, double **c, int nrows, int ncols)
{
  int x;

  for (x = 0; x < nrows; ++x) {
    vvm (a[x], b[x], c[x], ncols);
  }
}

void
sum2D (double *a, double **b, int nrows, int ncols)
{
  int x;

  vzero (a, ncols);
  for (x = 0; x < nrows; ++x) {
    vvp (a, a, b[x], ncols);
  }
}

double
total2D (double **a, int nrows, int ncols)
{
  int x;
  double sum = 0;

  for (x = 0; x < nrows; ++x) {
    sum += asum (a[x], ncols);
  }

  return sum;

}

int
total2Dint (int **a, int nrows, int ncols)
{
  int x, sum = 0;

  for (x = 0; x < nrows; ++x) {
    sum += intsum (a[x], ncols);
  }

  return sum;

}


/** 
 mixed modulus coding (see .../popgen/kimfitdir     
*/
long
lkodeitbb (int *xx, int len, int *baselist)
{
  int i, base;
  long t = 0;

  for (i = 0; i < len; ++i) {
    base = baselist[i];
    t *= base;
    t += xx[i];
    if (t < 0)
      fatalx ("(lkodeitbb) overflow\n");
  }
  return t;
}

int
ldekodeitbb (int *xx, long kode, int len, int *baselist)
{
// return weight

  int i, base;
  long t;

  t = kode;
  for (i = len - 1; i >= 0; --i) {
    base = baselist[i];
    xx[i] = t % base;
    t /= base;
  }
  return intsum (xx, len);

}

int
kodeitbb (int *xx, int len, int *baselist)
{
  int t = 0, i, base;

  for (i = 0; i < len; ++i) {
    base = baselist[i];
    t *= base;
    t += xx[i];
    if (t < 0)
      fatalx ("(kodeitbb) overflow\n");
  }
  return t;
}

int
dekodeitbb (int *xx, int kode, int len, int *baselist)
{
// return weight

  int i, t, base;

  t = kode;
  for (i = len - 1; i >= 0; --i) {
    base = baselist[i];
    xx[i] = t % base;
    t /= base;
  }
  return intsum (xx, len);

}


long expmod(long a, long b, long n) 
{ 
 int t ; 
 long ax=1, bx, z, z2 ; 
 t = b % 2 ;  
 if (t==1) ax = a ; 
 bx = b/2; 
 if (bx == 0) return ax % n ; 
 z = expmod(a, bx, n) ; 
 z2 = (z*z) % n ; 
 z2 = (ax*z2) % n ; 
 
 return z2 ; 

}

long
nextprime (long num)
// return nextprime >= num
{
  long x, q;
  int t;

  for (x = num;; ++x) {
    q = expmod(2, x-1, x) ;  
    if (q != 1 ) continue ; 
    t = isprime (x);
    if (t == YES)
      return x;
  }
}

int
isprime (long num)
// naive algorithm.  Implement Pollard rho at some time
{
  int top, x, t;

  if (num < 2)
    return NO;
  if (num == 2)
    return YES;
  top = nnint (sqrt (num));

  t = num % 2 ; 
  if (t==0) return NO ;  

  for (x = 3; x <= top; x += 2) {
    t = num % x;
    if (t == 0)
      return NO;
  }

  return YES;

}


int
irevcomp (int xx, int stringlen)
// consists of stringlen "mininibbles" (2 bits) 
{
  int aa[32], xxx, k, t;

  if (stringlen > 16)
    fatalx ("stringlen > 16\n");
  xxx = xx;
  for (k = 0; k < stringlen; ++k) {
    aa[k] = (xxx & 3) ^ 3;
    xxx = xxx >> 2;
  }
  xxx = 0;
  for (k = 0; k < stringlen; ++k) {
    t = aa[k];
    xxx = (xxx << 2) | t;
  }
  return xxx;
}

long
lrevcomp (long xx, int stringlen)
// consists of stringlen "mininibbles" (2 bits) 
// could be rewritten to avoid array aa + 1 loop.
{
  int aa[32], k, t;
  long xxx;

  if (stringlen > 32)
    fatalx ("stringlen > 32\n");
  xxx = xx;
  for (k = 0; k < stringlen; ++k) {
    aa[k] = (xxx & 3) ^ 3;
    xxx = xxx >> 2;
  }
  xxx = 0;
  for (k = 0; k < stringlen; ++k) {
    t = aa[k];
    xxx = (xxx << 2) | t;
  }
  return xxx;
}

void
ismatch (int *a, int *b, int n, int val)
{

  int i;

  for (i = 0; i < n; i++) {
    if (b[i] == val)
      a[i] = YES;
    else
      a[i] = NO;

  }
}

int
pmult (double *a, double *b, double *c, int nb, int nc)
// polynomial multiplication 
{
  double *ww;
  int i, j;

  ZALLOC (ww, nb + nc + 1, double);

  for (i = 0; i <= nb; ++i) {
    for (j = 0; j <= nc; ++j) {
      ww[i + j] += b[i] * c[j];
    }
  }

  copyarr (ww, a, nb + nc + 1);
  free (ww);

  return nb + nc;

}

void
pdiff (double *a, double *b, int deg)
// differentiate univariate polynomial
{
  double *ww, y;
  int k;

  ZALLOC (ww, deg + 1, double);
  for (k = 1; k <= deg; ++k) {
    y = (double) k;
    ww[k - 1] = y * b[k];
  }

  copyarr (ww, a, deg + 1);
  free (ww);
}

int
mktriang (double *out, double *in, int n)
{
  int a, b, x = 0;
  for (a = 0; a < n; ++a) {
    for (b = a; b < n; ++b) {
      out[x] = in[a * n + b];
      ++x;
    }
  }
  return x;
}

int
mkfull (double *out, double *in, int n)
// inverse to mktriang
{
  int a, b, x = 0;
  for (a = 0; a < n; ++a) {
    for (b = a; b < n; ++b) {
      out[a * n + b] = out[b * n + a] = in[x];;
      ++x;
    }
  }
  return x;
}
void vswap(double *a, double *b, int n) 
{
  double *w ; 

  ZALLOC(w, n, double) ;

  copyarr(a, w, n) ; 
  copyarr(b, a, n) ; 
  copyarr(w, b, n) ;

  free(w) ;
}

void setlong(long *pplen, long a, long b)  
// *pplen is a*b with check for overflow
{
  long long int xx ; 

  xx = a*b ;  
  if (xx > LONG_MAX) fatalx("overflow:   Are you on a 32 bit machine?\n") ;
  *pplen = xx ;


}

long lmod (long x, long base) 
// like % but handles - numbers better
{
  long t ; 

  t = x % base ; 
  if (t < 0) return lmod(t+base, base) ; 
  return t ; 


}


void
printlmat (long *a, int m, int n)

/**
 print a matrix n wide m rows
*/
{
  int i, j, jmod;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      printf ("%10ld ", a[i * n + j]);
      jmod = (j + 1) % 5;
      if ((jmod == 0) && (j < (n - 1))) {
        printf ("  ...\n");
      }
    }
    printf ("\n");
  }
}



void
floatitl (double *a, long *b, int n)
{
  int i;
  for (i = 0; i < n; i++) {
    a[i] = (double) b[i];
  }
}

void
lvsp (long *a, long *b, long c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] + c;
}


void
lvvp (long *a, long *b, long *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] + c[i];
}

void
lvvm (long *a, long *b, long *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = b[i] - c[i];
}

long gcdx(long a, long b, long *x, long *y)
{
    long t, x1, y1; // To store results of recursive call
    if (a == 0)
    {
        *x = 0;
        *y = 1;
        return b;
    }
 
 
    // Update x and y using results of recursive
    // call
    t = gcdx(b%a, a, &x1, &y1) ;

    *x = y1 - (b/a) * x1;
    *y = x1;

    return t ;
    
    
 
}

// gcd = a*x + b*y 

long modinv(long a, long base) 
{
  long t, x, y, aa ; 

  aa = lmod(a, base) ; 
  if (aa==0) return 0 ;  // special case 
  t = gcdx(a, base, &x, &y) ; 
  if (t != 1) fatalx("(modinv) %ld %ld\n", a, base) ; 

  return x ; 


}

long lpow2(int n) 
{
 long x = 1 ; 

 return x << n ; 


}

double cputime (int mode) 
{
  static double ttt=0 ; 

 if (mode==0) {  
  ttt = clocktime() ;
  return 0 ;
 }

 return clocktime() - ttt ; 

}

double calcmem (int mode) 
{
  static void *ttt = 0 ; 

 if (mode==0) {  
  ttt = topheap() ;
  return 0 ;
 }

 return (double) (topheap()  - ttt) ; 

}

double exp1minus(double x) 
// 1 - exp(-x) to good precision
{
 double ans, top, bot, term ; 
 int n ; 

 if (fabs(x)>.001) return 1.0 - exp(-x) ; 
 
 ans = 0; top = x ; bot = 1 ; n = 1 ; 

 for (;;) { 
   term = top/bot ;
   ans += term ; 
   if (fabs(term) < 1.0e-20) break ; 
   top *= -x ; 
   ++n ; 
   bot *= (double) n ; 
 }
 return ans ;
}


