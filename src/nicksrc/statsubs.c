#include <stdio.h>
#include <math.h>
#include <string.h>

#include "statsubs.h"
#include "vsubs.h"
#include "linsubs.h"
#include "ranmath.h"
#include "sortit.h"

#include "twtable.h"

#define EPS1 .001
#define EPS2 1.0e-12
#define ZLIM 20
#define QSIZE 10


static double *bern;		/* bernouilli numbers */
static int bernmax = 0;
static double *ztable = NULL, *ptable = NULL;
static double ptiny;
static int numbox = QSIZE * ZLIM;



static double zzprob (double zval);
static double znewt (double z, double ptail);

static double ltlg1 (double a, double x);
static double ltlg2 (double a, double x);
static double rtlg1 (double a, double x);
static double rtlg2 (double a, double x);
static double pochisq (double x, int df);
static double pof (double F, int df1, int df2);
static double betacf (double a, double b, double x);
static void weightjackx (double *est, double *sig, double mean, double *jmean,
			 double *jwt, int g);

static int twtabsize = -1;
static double *twxval, *twxpdf, *twxtail;

static char *twxtable = NULL;

static double **bcotable = NULL;
static int bcosize = -1;	// max val of binomial coefficient  


double
nordis (double zval)
/* normal density */
{
  double pi, t;

  pi = 2.0 * acos (0.0);

  t = exp (-0.5 * zval * zval);
  t /= sqrt (2.0 * pi);

  return t;

}

double
ndens (double val, double mean, double sig)
{
  return nordis ((val - mean) / sig) / sig;
}

double
ntail (double zval)
/** normal distribution tail area 
 uses erfc 
*/
{
  double pi, t;
  double p, q, d;

  if (zval == 0.0)
    return 0.5;
  if (zval < 0.0)
    return (1.0 - ntail (-zval));
  if (zval < ZLIM)
    {
      t = zval / sqrt (2.0);
      q = erfc (t) / 2.0;
      return q;
    }

  pi = 2.0 * acos (0.0);

  t = exp (-0.5 * zval * zval);
  t /= (sqrt (2.0 * pi) * zval);

  return t;

}

double
zzprob (double pval)
{
  double x, dev, p, q, d, h, u;
  double pi;
  int iter;

/** approximate normal by 1/(sqrt 2 pi) * exp (-0.5*x*x) / x   */
/* Feller I  page 166 */

  if (pval == 0.0)
    return 50.0;

  pi = 2.0 * acos (0.0);
  u = -log (sqrt (2.0 * pi) * pval);
/* solve x*x/2 + log(x) = u */

  x = sqrt (2.0 * u);
  for (iter = 1; iter <= 10; ++iter)
    {
      q = (0.5 * x * x) + log (x);
      d = x + (1.0 / x);
      dev = u - q;
      h = dev / d;
      x += h;
      if (fabs (h) < 1.0e-7)
	return x;
    }
  return x;
}

double
medchi (int *cls, int len, int *n0, int *n1, double *xtail)
{
/* compute 2x2 chisq splitting at median */
  int i, m0, m1, n, m;
  double arr[4], y, ys, p, q, d;
  *n0 = *n1 = 0;
  for (i = 0; i < len; i++)
    {
      if (cls[i] > 1)
	continue;
      if (cls[i] < 0)
	continue;
      if (cls[i] == 0)
	++ * n0;
      if (cls[i] == 1)
	++ * n1;
    }
  if (MIN (*n0, *n1) == 0)
    {
      *xtail = 1.0;
      return 0;
    }
  m = (*n0 + *n1) / 2;
  m0 = m1 = 0;
  for (i = 0; i < len; i++)
    {
      if (cls[i] > 1)
	continue;
      if (cls[i] < 0)
	continue;
      if (cls[i] == 0)
	++m0;
      if (cls[i] == 1)
	++m1;
      if ((m0 + m1) == m)
	break;
    }

  arr[0] = (double) m0;
  arr[1] = (double) m1;
  arr[2] = (double) (*n0 - m0);
  arr[3] = (double) (*n1 - m1);

  y = conchi (arr, 2, 2);
  ys = sqrt (y + EPS2);
  q = ntail (ys);

  *xtail = q;

  return y;

}

double
ks2 (int *cls, int len, int *n0, int *n1, double *kstail)
{
/*
 compute KS statistic 
 cls should be 0 or 1.  if larger take as invalid
*/
  int i;
  double en0, en1, en;
  double y, ymax;
/* count class sizes */

  if (len <= 1)
    {
      *kstail = 1.0;
      return 0;
    }

  *n0 = *n1 = 0;
  for (i = 0; i < len; i++)
    {
      if (cls[i] > 1)
	continue;
      if (cls[i] < 0)
	continue;
      if (cls[i] == 0)
	++ * n0;
      if (cls[i] == 1)
	++ * n1;
    }
  if (MIN (*n0, *n1) == 0)
    {
/**
  printf("warning ks2 has only 1 class passed\n") ;
  for (i=0; i<len ; i++) {
   printf("zz1 %d %d\n",i,cls[i]) ;
  }
*/
      *kstail = 1.0;
      return 0;
    }

  en0 = (double) *n0;
  en1 = (double) *n1;


  ymax = y = 0.0; /* running stat */ ;
  for (i = 0; i < len; i++)
    {
      if (cls[i] > 1)
	continue;
      if (cls[i] < 0)
	continue;
      if (cls[i] == 0)
	y += 1.0 / en0;
      if (cls[i] == 1)
	y -= 1.0 / en1;
      ymax = MAX (ymax, fabs (y));
    }

/*  Numerical recipes p 626 */
  en = sqrt (en0 * en1 / (en0 + en1));
  y = en + .12 + (0.11 / en);
  y *= ymax;
  *kstail = probks (y);
  return y;
/** crude analysis:  
  variance of 1 step above is (1/n0 + 1/n1) / (n0+n1) 
  and so variance of y is brownian motion not bridge is (1/n0+1/n1) 
  We want to rescale y to correspond to Brownian bridge.  
  First order correction is en.  We actually use 
   a Bartlett correction of some sort 
  Normalized y seems like what to return.
*/

}

double
probks (double lam)
/* KS tail area: Numerical recipes p 626 */
{
  int j;
  double a2, fac = 2.0, sum = 0.0, term, termbf = 0.0;
  double t;

  a2 = -2.0 * lam * lam;
  for (j = 1; j <= 100; j++)
    {
      t = a2 * (double) (j * j);
      term = fac * exp (t);
      sum += term;
      t = fabs (term);
      if ((t <= EPS1 * termbf) || (t <= EPS2 * sum))
	return sum;
      fac = -fac;
      termbf = fabs (term);
    }
  return 1.0;
}

double
conchiv (double *a, int m, int n)
/* a is m rows n columns.  contingency chisq */
{
  double *rsum, *csum, ee, tot = 0, chsq = 0, y;
  int i, j, k;

  ZALLOC (rsum, m, double);
  ZALLOC (csum, n, double);

  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
	{
	  k = i * n + j;
	  rsum[i] += a[k];
	  csum[j] += a[k];
	  tot += a[k];
	}
    }
  if (tot < 0.001)
    fatalx ("(conchiv) no data\n");
  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
	{
	  k = i * n + j;
	  ee = rsum[i] * csum[j] / tot;
	  if (ee < EPS2)
	    {
	      printf ("bad conchi\n");
	      printmat (a, m, n);
	      fatalx ("(conchiv) zero row or column sum\n");
	    }
	  y = a[k] - ee;
	  printf ("conchiv: %4d %4d ", i, j);
	  printf ("%9.3f ", a[k]);
	  printf ("%9.3f   ", ee);
	  printf ("%9.3f ", y);
	  printf ("%9.3f ", y / sqrt (ee));
	  printnl ();
	  chsq += (y * y) / ee;
	}
    }
  free (rsum);
  free (csum);
  return chsq;
}

double
z2x2 (double *a)
/* a is 2 rows 2 columns.  contingency z-score */
{
  double *rsum, *csum, ee, tot = 0, chsq = 0, y, dev00, z;
  int i, j, k, m = 2, n = 2;

  ZALLOC (rsum, m, double);
  ZALLOC (csum, n, double);

  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
	{
	  k = i * n + j;
	  rsum[i] += a[k];
	  csum[j] += a[k];
	  tot += a[k];
	}
    }
  if (tot < 0.001)
    fatalx ("(z2x2) no data\n");
  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
	{
	  k = i * n + j;
	  ee = rsum[i] * csum[j] / tot;
	  if (ee < EPS2)
	    {
	      printf ("bad conchi\n");
	      printmat (a, m, n);
	      fatalx ("(conchi) zero row or column sum\n");
	    }
	  y = a[k] - ee;
	  if (k == 0)
	    dev00 = y;
	  chsq += (y * y) / ee;
	}
    }
  z = sqrt (chsq + 1.0e-12);
  if (dev00 < 0.0)
    z = -z;
  free (rsum);
  free (csum);
  return z;
}


double
conchi (double *a, int m, int n)
/* a is m rows n columns.  contingency chisq */
{
  double *rsum, *csum, ee, tot = 0, chsq = 0, y;
  int i, j, k;

  ZALLOC (rsum, m, double);
  ZALLOC (csum, n, double);

  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
	{
	  k = i * n + j;
	  rsum[i] += a[k];
	  csum[j] += a[k];
	  tot += a[k];
	}
    }
  if (tot < 0.001)
    fatalx ("(conchi) no data\n");
  for (i = 0; i < m; i++)
    {
      for (j = 0; j < n; j++)
	{
	  k = i * n + j;
	  ee = rsum[i] * csum[j] / tot;
	  if (ee < EPS2)
	    {
	      printf ("bad conchi\n");
	      printmat (a, m, n);
	      fatalx ("(conchi) zero row or column sum\n");
	    }
	  y = a[k] - ee;
	  chsq += (y * y) / ee;
	}
    }
  free (rsum);
  free (csum);
  return chsq;
}

double
chitest (double *a, double *p, int n)
/* a is n boxes.  Goodness of fit test to p */
{

  double *x, *b, *pp;
  double y1 = 0.0, y2 = 0.0;
  int i;

  ZALLOC (pp, n, double);
  if (p != NULL)
    copyarr (p, pp, n);
  else
    vclear (pp, 1.0, n);

  y1 = asum (pp, n);
  y2 = asum (a, n);

  if ((y1 == 0.0) || (y2 == 0.0))
    {
      free (pp);
      return 0.0;
    }

  ZALLOC (x, n, double);
  ZALLOC (b, n, double);


  vst (x, pp, y2 / y1, n);	/* expected */

  vsp (x, x, .0001, n);
  vvm (b, a, x, n);
  vvt (b, b, b, n);
  vvd (b, b, x, n);

  y1 = asum (b, n);

  free (x);
  free (b);

  return y1;

}

double
zprob (double ptail)
/** inverse normal */
{
  double z, p, t, plog;
  double ylo, yhi, ya, yb;
  int i, k;

  if (ztable == NULL)
    setzptable ();
  if (ptail == 0.5)
    return 0.0;
  if (ptail > 0.5)
    {
      z = zprob (1.0 - ptail);
      return -z;
    }
  if ((ptail > .40) && (ptail < .50))
    return znewt (0.0, ptail);
  if (ptail < ptiny)
    {
      z = zzprob (ptail);
      return znewt (z, ptail);
    }
/** replace by binary or interpolating search */
  plog = -log (ptail);
  k = firstgt (plog, ptable, numbox);
  if (k == 0)
    return ztable[0];
  if (k == numbox)
    return ztable[numbox - 1];
  ylo = ptable[k - 1];
  yhi = ptable[k];
  ya = (yhi - plog);
  yb = plog - ylo;
  z = (ya * ztable[k - 1] + yb * ztable[k]) / (ya + yb);
  if (isnan (z))
    fatalx ("zprob bug %15.9f %15.9f\n", z, ptail);
  t = znewt (z, ptail);
  if (isnan (t)) {
    fatalx ("zprob bug %15.9f %15.9f\n", z, ptail);
  }
  return t;
}

void
setzptable ()
{
  int i;
  double p, z;
  if (ptable != NULL)
    free (ptable);
  if (ztable != NULL)
    free (ztable);
  ZALLOC (ptable, numbox, double);
  ZALLOC (ztable, numbox, double);
  z = 0.0;
  for (i = 0; i < numbox; i++)
    {
      p = ntail (z);
      ztable[i] = z;
      ptable[i] = -log (p);
      z += 1.0 / (double) QSIZE;
    }
  ptiny = p;
}

double
znewt (double z, double ptail)
{
/** 
 newton step 
 z is very good approximation 
*/
  double p0, pder, h;
  double pi, zz;
  int iter;
  pi = 2.0 * acos (0.0);
  zz = z;
  for (iter = 1; iter <= 5; ++iter)
    {
      p0 = ntail (zz);
      pder = -exp (-0.5 * zz * zz) / sqrt (2 * pi);
      if (pder == 0.0)
	return zz;
      h = (ptail - p0) / pder;
      if (fabs (h) <= 1.0e-10)
	return zz;
      zz += h;
    }
  return zz;
}

int
ifirstgt (int val, int *tab, int n)
{
/* tab sorted in ascending order */
  int i;

  if (val >= tab[n - 1])
    return n;
  for (i = 0; i < n; i++)
    {
      if (val < tab[i])
	return i;
    }

    return n ; 
}

int
firstgt (double val, double *tab, int n)
{
/* tab sorted in ascending order */
  int i;

  if (val >= tab[n - 1])
    return n;
  for (i = 0; i < n; i++)
    {
      if (val < tab[i])
	return i;
    }

    return n ;
}

void
mleg (double a1, double a2, double *p, double *lam)
{
  int iter;
  double s, pp, ll;
  double top, bot, fval;
  int debug = NO;

/** 
 solve 
 p/lam = a1 ; psi(p) - log(lam) = a2 ;
 Thus psi(p) - log(p) = a2 - log(a1) 
*/
  s = a2 - log (a1);

  if (s >= 0.0)
    fatalx ("log E(x) < E(log (x)) \n");
  pp = -s;

  for (iter = 1; iter <= 30; ++iter)
    {
      fval = s - (psi (pp) - log (pp));
      if (debug)
	printf ("yy1 %3d %9.3f %9.3f\n", iter, pp, fval);
      if (fval < 0.0)
	break;
      pp *= 2.0;
    }

  for (iter = 1; iter <= 30; ++iter)
    {
      fval = s - (psi (pp) - log (pp));
      if (fval > 0.0)
	break;
      if (debug)
	printf ("yy2 %3d %9.3f %9.3f\n", iter, pp, fval);
      pp /= 2.0;
    }

  for (iter = 1; iter <= 10; ++iter)
    {
      fval = psi (pp) - log (pp);
      top = s - fval;
      bot = tau (pp) - (1.0 / pp);
      if (debug)
	printf ("%3d %12.6f %12.6f\n", iter, pp, top);
      pp += top / bot;
    }
  ll = pp / a1;
  *p = pp;
  *lam = ll;
}


double
psi (double x)
{
  double y, zz, term;
  int k;

  if (x <= 0.0)
    fatalx ("(psi) bad value:  %9.3f\n", x);
  bernload ();
  if (x < 10.0)
    return (psi (x + 1.0) - 1.0 / x);

  y = log (x) - 1.0 / (2.0 * x);
  zz = 1.0;
  for (k = 1; k <= bernmax / 2; k++)
    {
      zz /= (x * x);
      term = bernum (2 * k) / (double) (2 * k);
      term *= zz;
      y -= term;
    }
  return y;
}

double
tau (double x)
/*
 derivative of psi 
*/
{
  double y, zz, term;
  int k;

  if (x <= 0.0)
    fatalx ("(tau) bad value:  %9.3f\n", x);
  bernload ();
  if (x < 10.0)
    return (tau (x + 1.0) + 1.0 / (x * x));

  y = 1.0 / x + 1.0 / (2.0 * x * x);
  zz = 1.0 / x;
  for (k = 1; k <= bernmax / 2; k++)
    {
      zz /= (x * x);
      term = bernum (2 * k) / (double) (2 * k);
      term *= zz;
      term *= -(double) (2 * k);
      y -= term;
    }
  return y;
}

// Formulae from A + S p 374
double
logbessi1 (double x)
{
  double yans, logyans;
  double y;


  if (x < 0.0)
    fatalx ("bad bessi1\n");
  if (x < 3.75)
    {
      y = x / 3.75, y = y * y;
      yans = x * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
								 +
								 y *
								 (0.2658733e-1
								  +
								  y *
								  (0.301532e-2
								   +
								   y *
								   0.32411e-3))))));
      logyans = log (yans);
    }
  else
    {
      y = 3.75 / x;
      yans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
						      - y * 0.420059e-2));
      yans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
						    + y * (0.163801e-2 +
							   y *
							   (-0.1031555e-1 +
							    y * yans))));
      logyans = log (yans);
      logyans += x;
      logyans -= log (x) / 2;
    }
  return logyans;
}




double
bessi1 (double x)
// I1(x)
{
  double yans;
  double y;


  if (x < 0.0)
    fatalx ("bad bessi1\n");
  if (x < 3.75)
    {
      y = x / 3.75, y = y * y;
      yans = x * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
								 +
								 y *
								 (0.2658733e-1
								  +
								  y *
								  (0.301532e-2
								   +
								   y *
								   0.32411e-3))))));
    }
  else
    {
      y = 3.75 / x;
      yans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
						      - y * 0.420059e-2));
      yans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
						    + y * (0.163801e-2 +
							   y *
							   (-0.1031555e-1 +
							    y * yans))));
      yans *= exp (x) / sqrt (x);
    }
  return yans;
}

double
logbessi0 (double x)
{
// log I0(x)

  double yans, logyans;
  double y;

  if (x < 0.0)
    fatalx ("bad bessi0\n");

  if (x < 3.75)
    {
      y = x / 3.75, y = y * y;
      yans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
							  + y * (0.2659732 +
								 y *
								 (0.360768e-1
								  +
								  y *
								  0.45813e-2)))));
      logyans = log (yans);
    }
  else
    {
      y = 3.75 / x;
      yans = 0.39894228 + y * (0.1328592e-1
			       + y * (0.225319e-2 +
				      y * (-0.157565e-2 +
					   y * (0.916281e-2 +
						y * (-0.2057706e-1 +
						     y * (0.2635537e-1 +
							  y * (-0.1647633e-1 +
							       y *
							       0.392377e-2)))))));
      logyans = log (yans);
      logyans += x;
      logyans -= log (x) / 2;
    }
  return logyans;
}


double
bessi0 (double x)
{
// I0(x)

  double ans;
  double y;

  if (x < 0.0)
    fatalx ("bad bessi0\n");

  if (x < 3.75)
    {
      y = x / 3.75, y = y * y;
      ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
							 + y * (0.2659732 +
								y *
								(0.360768e-1 +
								 y *
								 0.45813e-2)))));
    }
  else
    {
      y = 3.75 / x;
      ans = 0.39894228 + y * (0.1328592e-1
			      + y * (0.225319e-2 +
				     y * (-0.157565e-2 +
					  y * (0.916281e-2 +
					       y * (-0.2057706e-1 +
						    y * (0.2635537e-1 +
							 y * (-0.1647633e-1 +
							      y *
							      0.392377e-2)))))));
      ans *= exp (x) / sqrt (x);
    }
  return ans;
}




void
bernload ()
{
  if (bernmax > 0)
    return;
  bernmax = 14;
  ZALLOC (bern, bernmax + 1, double);
  bern[0] = 1.0;
  bern[1] = -1.0 / 2.0;
  bern[2] = 1.0 / 6.0;
  bern[4] = -1.0 / 30.0;
  bern[6] = 1.0 / 42.0;
  bern[8] = -1.0 / 30.0;
  bern[10] = 5.0 / 66.0;
  bern[12] = -691.0 / 2730.0;
  bern[14] = 7.0 / 6.0;
}

double
bernum (int k)
{
  bernload ();
  if ((k < 0) || (k > bernmax))
    fatalx ("(bernum) bad arg: %s\n", k);
  return (bern[k]);
}

double
dilog (double x)
{
  return li2 (x);
}

double
li2 (double x)
{
  double pi, sum = 0.0, term, top, z;
  int k;

  pi = acos (0.0) * 2.0;
  if (x <= 0.0)
    return pi * pi / 6.0;
  if (x >= 1.0)
    return 0;
  if (x < 0.5)
    {
      return (-log (x) * log (1 - x) + (pi * pi / 6.0) - li2 (1.0 - x));
    }
  z = 1 - x;
  top = 1.0;
  for (k = 1; k <= 100; k++)
    {
      top *= z;
      term = top / (double) (k * k);
      sum += term;
      if (term <= 1.0e-20)
	break;
    }
  return sum;
}

double
hwstat (double *x)
/** Hardy-Weinberg equilibrium test 
    returns standard normal in null case. 
    +sign is excess heterozygosity 
    x[0] [1] [2] are counts for homozm hetero homo (alt allele)
*/
{

  double p, q, ysum, s1, y1, y2, ychi, sig;
  double a1[3], a2[3];

  ysum = asum (x, 3);
  if (ysum < 0.001)
    return 0.0;
  s1 = 2 * x[2] + x[1];
  p = 0.5 * s1 / ysum;
  q = 1.0 - p;

  a1[0] = q * q;
  a1[1] = 2 * p * q;
  a1[2] = p * p;

  vsp (a1, a1, 1.0e-8, 3);
  vst (a2, x, 1.0 / ysum, 3);
  vsp (a2, a2, 1.0e-8, 3);

  y2 = vldot (x, a2, 3);
  y1 = vldot (x, a1, 3);

  ychi = 2.0 * (y2 - y1);
  sig = sqrt (ychi + 1.0e-8);

  if (a2[1] < a1[1])
    sig = -sig;
/* negative => hets lo */

  return sig;


}

double
bprob (double p, double a, double b)
{
  double q, yl;
  q = 1.0 - p;
  yl = (a - 1) * log (p) + (b - 1) * log (q);
  if (!isfinite (yl))
    fatalx ("bad bprob %9.3f %9.3f %9.3f\n", p, a, b);
  yl -= lbeta (a, b);
  if (!isfinite (yl))
    fatalx ("bad bprob (lbeta) %9.3f %9.3f %9.3f\n", p, a, b);
  return yl;
}

double
gammprob (double x, double p, double lam)
/* gamma density */
{
  double xx, yl;
  xx = MAX (x, 1.0e-8);
  yl = (p - 1) * log (xx) - lam * xx;
  yl += p * log (lam);
  yl -= xlgamma (p);
  return yl;
}


double
lbeta (double a, double b)
{
  return (xlgamma (a) + xlgamma (b) - xlgamma (a + b));
}


double
dawson (double t)
/** 
 Dawson's Integral 
 [A + S 7.31]
 exp(-t*t) \int ( exp(x^2), x = 0..t) 
 loosely based on mcerror.for 
*/
{

  double z1, cs, cr, cl;
  double z1sq, cer;
  int k;

  z1 = fabs (t);
  if (z1 <= 1.0e-8)
    return t;
/* derivative is 1 at 0 */
  z1sq = -t * t;
  if (z1 < 4.5)
    {
      cs = cr = z1;
      for (k = 1; k <= 999; ++k)
	{
	  cr *= z1sq / ((double) k + 0.5);
	  cs += cr;
	  if (fabs (cr / cs) < 1.0e-15)
	    break;
	}
      cer = cs;
    }
  else
    {
      cl = 1 / z1;
      cr = cl;
      for (k = 1; k <= 13; ++k)
	{
	  cr *= -((double) k - 0.5) / z1sq;
	  cl += cr;
	  if (fabs (cr / cl) < 1.0e-15)
	    break;
	}
      cer = 0.5 * cl;
    }
  if (t < 0)
    cer = -cer;
  return cer;
}

double
binlogtail (int n, int t, double p, char c)
{

  double *bindis;
  double val, base;

  ZALLOC (bindis, n + 1, double);
  genlogbin (bindis, n, p);
  base = bindis[t];
  vsp (bindis, bindis, -base, n + 1);
  if (c == '+')
    {
      vexp (bindis + t, bindis + t, n - t + 1);
      val = asum (bindis + t, n - t + 1);
    }
  else
    {
      vexp (bindis, bindis, t);
      val = asum (bindis, t);
    }
  free (bindis);
  return (log (val) + base);

}

double
binomtail (int n, int t, double p, char c)
{
/** 
 c = '+':   P(S>=t) 
 c = '-':   P(S<t) 
 WARNING <= t use binomtail(n, t+1, ... 
*/
  double *bindis;
  double val;

  ZALLOC (bindis, n + 1, double);
  genlogbin (bindis, n, p);
  vexp (bindis, bindis, n + 1);
  if (c == '+')
    val = asum (bindis + t, n - t + 1);
  else
    val = asum (bindis, t);
  free (bindis);
  return val;
}

void
genbin (double *a, int n, double p)
{
  vzero (a, n + 1);
  if (p <= 0.0)
    {
      a[0] = 1.0;
      return;
    }
  if (p >= 1.0)
    {
      a[n] = 1.0;
      return;
    }

  genlogbin (a, n, p);
  vexp (a, a, n + 1);
  bal1 (a, n + 1);

}

void
genlogbin (double *a, int n, double p)
/* generate log prob for binomial distribution */
// a must be n+1 long
{
  double q, plog, qlog;
  double *lfac, x;
  int i, r, s;
  q = 1.0 - p;
  plog = log (p), qlog = log (q);

  ZALLOC (lfac, n + 1, double);
  for (i = 1; i <= n; i++)
    {
      x = (double) i;
      lfac[i] = lfac[i - 1] + log (x);
    }
  for (r = 0; r <= n; r++)
    {
      s = n - r;
      x = lfac[n] - lfac[r] - lfac[s];	/* log binom coeff */
      x += ((double) r) * plog;
      x += ((double) s) * qlog;
      a[r] = x;			/* log prob */
    }
  free (lfac);
}

double
xlgamma (double x)
// a version of lgamma since the linux version is 
// causing trouble ???
{
  static double l2pi = -1.0;
  double term, y, zz;
  int k, t;

  if (x <= 0.0)
    fatalx ("xlgamma: x <= 0  %9.3f\n", x);
// return lgamma(x) ;

  t = (int) sizeof (long);
  if (t >= 8)
    return lgamma (x);		// 64 bit machines 
  if (l2pi < 0.0)
    l2pi = log (4.0 * acos (0.0));
  bernload ();
  if (x < 10.0)
    return (xlgamma (x + 1.0) - log (x));
  y = (x - 0.5) * log (x) - x + 0.5 * l2pi;
  zz = 1.0 / x;
  for (k = 1; k <= bernmax / 2; k++)
    {
      t = 2 * k;
      term = bernum (2 * k) / (double) (t * (t - 1));
      term *= zz;
      y += term;
      zz /= (x * x);
    }
  if (!isfinite (y))
    fatalx ("bad xlgamma\n");
  return y;
}

double
rtlchsq (int df, double z)
{
  double a, x, y;
  if (df == 1)
    return 2.0 * ntail (sqrt (z));
  if (df == 2)
    return exp (-0.5 * z);

  if (df <= 0.0)
    return 1.0;
  y = pochisq (z, df);
  if (y < 1.0e-6)
    {
      a = 0.5 * (double) df;
      x = 0.5 * z;
      return rtlg (a, x);
    }
  return y;

}

/*ALGORITHM Compute probability of chi square value.
	Adapted from:
		Hill, I. D. and Pike, M. C.  Algorithm 299
		Collected Algorithms for the CACM 1967 p. 243
	Updated for rounding errors based on remark in
		ACM TOMS June 1985, page 185
         Perlman.  No copyright
*/
double
pochisq (double x, int df)
{
  double a, y, s;
  double e, c, z;
  double poz ();		/* computes probability of normal z score */
  int even;			/* true if df is an even number */

  if (x <= 0.0 || df < 1)
    return (1.0);

  a = 0.5 * x;
  even = (2 * (df / 2)) == df;
  if (df > 1)
    y = ex (-a);
  s = (even ? y : (2.0 * ntail (sqrt (x))));
  if (df > 2)
    {
      x = 0.5 * (df - 1.0);
      z = (even ? 1.0 : 0.5);
      if (a > BIGX)
	{
	  e = (even ? 0.0 : LOG_SQRT_PI);
	  c = log (a);
	  while (z <= x)
	    {
	      e = log (z) + e;
	      s += ex (c * z - a - e);
	      z += 1.0;
	    }
	  return (s);
	}
      else
	{
	  e = (even ? 1.0 : (I_SQRT_PI / sqrt (a)));
	  c = 0.0;
	  while (z <= x)
	    {
	      e = e * (a / z);
	      c = c + e;
	      z += 1.0;
	    }
	  return (c * y + s);
	}
    }
  else
    return (s);
}

double
critchi (int df, double p)
/* Perlman.  arguments interchanged */
{
  double minchisq = 0.0;
  double maxchisq = 100.0 * p;
  double chisqval, z, y;;

  if (p <= 0.0)
    return (maxchisq);
  else if (p >= 1.0)
    return (0.0);

  if (df == 1)
    {
      z = zprob (0.5 * p);
      return z * z;
    }

  if (df == 2)
    {
      y = -log (p);
      return 2 * y;
    }

  chisqval = df / sqrt (p);	/* fair first value */
  while (maxchisq - minchisq > CHI_EPSILON)
    {
      if (rtlchsq (df, chisqval) < p)
	maxchisq = chisqval;
      else
	minchisq = chisqval;
      chisqval = (maxchisq + minchisq) * 0.5;
    }
  return (chisqval);
}

/*
	Module:       f.c
	Purpose:      compute approximations to F distribution probabilities
	Contents:     pof(), critf()
	Programmer:   Gary Perlman
	Organization: Wang Institute, Tyngsboro, MA 01879
	Tester:       compile with -DFTEST to include main program
	Copyright:    none
	Tabstops:     4
*/


#ifndef	I_PI			/* 1 / pi */
#define	I_PI        0.3183098861837906715377675
#endif
#define	F_EPSILON     0.000001	/* accuracy of critf approximation */
#define	F_MAX      9999.0	/* maximum F ratio */


double
rtlf (int df1, int df2, double F)
{
  return pof (F, df1, df2);
}

static double
pof (double F, int df1, int df2)
{
  int i, j;
  int a, b;
  double w, y, z, d, p;

  if (F < F_EPSILON || df1 < 1 || df2 < 1)
    return (1.0);
  a = df1 % 2 ? 1 : 2;
  b = df2 % 2 ? 1 : 2;
  w = (F * df1) / df2;
  z = 1.0 / (1.0 + w);
  if (a == 1)
    if (b == 1)
      {
	p = sqrt (w);
	y = I_PI;		/* 1 / 3.14159 */
	d = y * z / p;
	p = 2.0 * y * atan (p);
      }
    else
      {
	p = sqrt (w * z);
	d = 0.5 * p * z / w;
      }
  else if (b == 1)
    {
      p = sqrt (z);
      d = 0.5 * z * p;
      p = 1.0 - p;
    }
  else
    {
      d = z * z;
      p = w * z;
    }
  y = 2.0 * w / z;
#ifdef	REMARK			/* speedup modification suggested by Tolman (wrong answer!) */
  if (a == 1)
    for (j = b + 2; j <= df2; j += 2)
      {
	d *= (1.0 + a / (j - 2.0)) * z;
	p += d * y / (j - 1.0);
      }
  else
    {
      double zk = 1.0;
      for (j = (df2 - 1) / 2; j; j--)
	zk *= z;
      d *= zk * df2 / b;
      p *= zk + w * z * (zk - 1.0) / (z - 1.0);
    }
#else /* original version */
  for (j = b + 2; j <= df2; j += 2)
    {
      d *= (1.0 + a / (j - 2.0)) * z;
      p = (a == 1 ? p + d * y / (j - 1.0) : (p + w) * z);
    }
#endif
  y = w * z;
  z = 2.0 / z;
  b = df2 - 2;
  for (i = a + 2; i <= df1; i += 2)
    {
      j = i + b;
      d *= y * j / (i - 2.0);
      p -= z * d / j;
    }
  /* correction for approximation errors suggested in certification */
  if (p < 0.0)
    p = 0.0;
  else if (p > 1.0)
    p = 1.0;
  return (1.0 - p);
}

// incomplete gamma function.  P(z<x) if  standard gamma  shape a
double
ltlg (double a, double x)
{
  if (x <= 0.0)
    return 0.0;
  if (x <= (1.0 + a))
    return ltlg1 (a, x);
  return ltlg2 (a, x);
}

// incomplete gamma function.  P(z>x) if  standard gamma  shape a
double
rtlg (double a, double x)
{
  if (x <= 0.0)
    return 1.0;
  if (x <= (1.0 + a))
    return rtlg1 (a, x);
  return rtlg2 (a, x);
}

double
ltlg1 (double a, double x)
{
  double r, s, tiny = 1.0e-14;
  double yk, y1, xam;
  int k;

  s = 1.0 / a;
  r = s;
  for (k = 1; k <= 60; ++k)
    {
      yk = (double) k;
      r *= (x / (a + yk));
      s += r;
      if (fabs (r / s) < tiny)
	break;
    }
  xam = (a * log (x)) - x;
  y1 = xam + log (s);
  y1 -= xlgamma (a);
  y1 = exp (y1);
  if (isnan (y1))
    fatalx ("bad ltlg1\n");
  return y1;

}

double
ltlg2 (double a, double x)
{
  return 1.0 - rtlg2 (a, x);
}


double
rtlg1 (double a, double x)
{
  return 1.0 - ltlg1 (a, x);
}

double
rtlg2 (double a, double x)
{
  double y1, y2;
  double yk, top, bot, t0, xam;
  int k;
  t0 = 0.0;

// ZJ p 64 ff
  for (k = 60; k >= 1; --k)
    {
      yk = (double) k;
      top = yk - a;
      bot = yk / (x + t0);
      ++bot;
      t0 = top / bot;
    }
  xam = (a * log (x)) - x;
  y1 = xam - log (x + t0);
  y1 -= xlgamma (a);
  y1 = exp (y1);
  if (isnan (y1))
    fatalx ("bad rtlg2\n");
  return y1;
}

void cinterp (double val, double x0, double x1,
	      double f0, double f0p, double f1, double f1p, double *fv,
	      double *fvp);
int firstgtx (double val, double *tab, int n);
static int gtx (double *tab, int lo, int hi, double val);
void gettw (double x, double *tailp, double *densp);


double
twdens (double twstat)
// Tracy-Widom prob density
{
  double dens, tail;

  gettw (twstat, &dens, &tail);
  return dens;

}

double
twtail (double twstat)
// Tracy-Widom right tail
{

  double dens, tail;
  static int ncall = 0;


  ++ncall;

  gettw (twstat, &tail, &dens);
/**
  printf("zz %9.3f %9.3f\n", twstat, tail) ;
  if (ncall==10) abort() ;
*/
  return tail;

}

double
twdensx (double tw)
//  Margetis-Edelman
{
  double bot, lbot, ltop, y1, y2;

  if (tw <= 0.0)
    return 0.0;
  bot = (SQRT_PI) * 4.0;
  lbot = log (bot);
  y1 = -0.25 * log (tw);
  y2 = -2.0 * pow (tw, 1.5) / 3.0;
  ltop = y1 + y2;
  return exp (ltop - lbot);
}

double
twtailx (double tw)
// right tail Margetis-Edelman
{
  double bot, lbot, ltop, y1, y2;

  if (tw <= 0.0)
    return 1.0;
  bot = (SQRT_PI) * 4.0;
  lbot = log (bot);
  y1 = -0.75 * log (tw);
  y2 = -2.0 * pow (tw, 1.5) / 3.0;
  ltop = y1 + y2;
  return exp (ltop - lbot);
}

void
twfree ()
// destructor.  Here for completeness
{
  if (twtabsize < 0)
    return;
  free (twxval);
  free (twxpdf);
  free (twxtail);
  twtabsize = -1;



}


double
twnorm (double lam, double p, double n)
// Ref Johnstone (2001) 
{
  double mu, phi, y1, y2;

  if (n < 0.0)
    return -10.0;
  if (p < 0.0)
    return -10.0;

  if (n < p)
    return twnorm (lam, n, p);
// not very important refinement as twnorm symmetric in p, n-1  NJP

  y1 = sqrt (n - 1) + sqrt (p);
  mu = y1 * y1;
  y2 = (1.0 / sqrt (n - 1)) + 1.0 / sqrt (p);
  phi = y1 * pow (y2, 1.0 / 3.0);
  return (lam - mu) / phi;
}

double
dotwcalc (double *lambda, int m, double *ptw, double *pzn, double *pzvar,
	  int minm)
{
  double nv, mv, tzn, tm;
  double *evals;
  double y, top, bot, zn, tw, ystat;
  double tail, lsum;

  if (m < minm)
    {
      *pzn = *pzvar = *ptw = -1;
      return -1.0;
    }
  lsum = asum (lambda, m);
  if (lsum <= 0.0)
    {
      *pzn = *pzvar = *ptw = -1;
      return -1.0;
    }

  tzn = *pzn;
  tm = (double) m;

  y = (double) m / lsum;
  ystat = lambda[0] * y * tzn;

  if (tzn > 0.0)
    {
      tw = twnorm (ystat, tm, tzn);
      *pzn = tzn;
      *ptw = tw;
      tail = twtail (tw);
      return tail;
    }
  ZALLOC (evals, m, double);
  vst (evals, lambda, y, m);
  top = (double) (m * (m + 2));
  bot = asum2 (evals, m) - (double) m;
  zn = top / bot;		// see appendix to eigenpaper  NJP
  y = evals[0] * zn;
  tw = twnorm (y, tm, zn);
  *pzn = zn;
  *ptw = tw;
  tail = twtail (tw);
  free (evals);
  return tail;
}

int
numgtz (double *a, int n)
{
#define THRESH  .000001

  int num, k;

  num = 0;
  for (k = 0; k < n; ++k)
    {
      if (a[k] > THRESH)
	++num;
    }
  return num;
}


static int
fgtx (double *tab, int lo, int hi, double val)
{

  int k;

  if (val >= tab[hi])
    return hi + 1;
  if (val < tab[lo])
    return lo;
  k = (lo + hi) / 2;
  if (val <= tab[k])
    return fgtx (tab, lo + 1, k, val);
  return fgtx (tab, k, hi - 1, val);
}

int
firstgtx (double val, double *tab, int n)
// tab sorted in ascending order 
{
  return fgtx (tab, 0, n - 1, val);
}

static int
jfgtx (int *tab, int lo, int hi, int val)
{

  int k;

  if (val >= tab[hi])
    return hi + 1;
  if (val < tab[lo])
    return lo;
  k = (lo + hi) / 2;
  if (val <= tab[k])
    return jfgtx (tab, lo + 1, k, val);
  return jfgtx (tab, k, hi - 1, val);
}

int
jfirstgtx (int val, int *tab, int n)
// tab sorted in ascending order 
// return first index with tab[x] >= val
{
  return jfgtx (tab, 0, n - 1, val);
}

int
settwxtable (char *table)
{
  FILE *fff;
  if (twxtable != NULL)
    return 1;
  if (table == NULL)
//    twxtable = strdup (TWXTABLE);
    return 1;
  else
    twxtable = strdup (table);
  fff = fopen (twxtable, "r");
  if (fff == NULL)
    return -1;
  fclose (fff);
  return 1;
}

void
gettw (double x, double *tailp, double *densp)
// main routine for accessing twtable
{
  int k, n;
  double x0, x1, f0, f1, f0p, f1p;
  double *xx[3];


  if (twtabsize == -1)
    {
      twtabsize = TWTABSIZE;
      ZALLOC(twxval,  twtabsize, double);
      ZALLOC(twxpdf,  twtabsize, double);
      ZALLOC(twxtail, twtabsize, double);
      memcpy(twxval,  TWXVAL,  twtabsize*sizeof(double));
      memcpy(twxpdf,  TWXPDF,  twtabsize*sizeof(double));
      memcpy(twxtail, TWXTAIL, twtabsize*sizeof(double));
    }
  n = twtabsize;

  k = firstgtx (x, twxval, n);

  if (k <= 0)
    {
      *tailp = 1.0;
      *densp = 0.0;
      return;
    }

  if (k >= n)
    {
      *tailp = twdensx (x);
      *densp = twtailx (x);
      return;
    }

  x0 = twxval[k - 1];
  x1 = twxval[k];
  f0 = twxtail[k - 1];
  f0p = twxpdf[k - 1];
  f1 = twxtail[k];
  f1p = twxpdf[k];

// now do cubic interpolation
  cinterp (x, x0, x1, f0, -f0p, f1, -f1p, tailp, densp);
  *densp = -*densp;

/**
     printf("zzz %9.3f %9.3f %9.3f\n", x0, x1, x) ;
     printf("zz1 %9.3f %9.3f %9.3f\n", f0, f1, *tailp) ;
     printf("zz2 %9.3f %9.3f %9.3f\n", f0p, f1p, *densp) ;
*/

}

void
cinterp (double val, double x0, double x1,
	 double f0, double f0p, double f1, double f1p, double *fv,
	 double *fvp)
// cubic interpolation val should be between x0 and x1
// fv is function at x 
// fvp is derivative
{

  double inc, yval, f, fp, a0, b0, a1, b1, a2, a3;
  double c0, c1, cc0, cc1;

  inc = x1 - x0;
  yval = (val - x0) / inc;
  b0 = f0;
  b1 = f0p * inc;
  c0 = f1;
  c1 = f1p * inc;

  a0 = b0;
  a1 = b1;
  cc0 = c0 - (a0 + a1);
  cc1 = c1 - a1;
  a2 = 3 * cc0 - cc1;
  a3 = cc1 - 2 * cc0;

  f = a3;
  f *= yval;
  f += a2;
  f *= yval;
  f += a1;
  f *= yval;
  f += a0;
  *fv = f;
  fp = 3 * a3;
  fp *= yval;
  fp += 2 * a2;
  fp *= yval;
  fp += a1;

  *fv = f;
  *fvp = fp / inc;

}

double
dirmult (double *pp, int *aa, int len)
{
  int t, i, m;
  double y1, y2, ysum;
  double top, bot;

  m = len;

  t = intsum (aa, m);
  if (t < 1)
    return 0.0;

  top = bot = 0.0;
  ysum = asum (pp, m);
  for (i = 0; i < m; i++)
    {
      top += lgamma (pp[i] + (double) aa[i]);
      bot += lgamma (pp[i]);
    }
  top += lgamma (ysum);
  bot += lgamma (ysum + (double) t);


  y1 = top - bot;

  return y1 - y2;
}


double
betaix (double a, double b, double lo, double hi)
{

  double y1, y2;

  y1 = betai (a, b, lo);
  y2 = betai (a, b, hi);

  if (!isfinite (y1))
    fatalx ("bad y1\n");
  if (!isfinite (y2))
    fatalx ("bad y2\n");

  return y2 - y1;

}

double
betai (double a, double b, double x)
{
  double betacf (double a, double b, double x);
  double bt;

  if (x < 0.0 || x > 1.0)
    fatalx ("Bad x in routine betai\n");
  if (x == 0.0)
    return 0.0;
  if (x == 1.0)
    return 1.0;
  /* Factors in front of the continued fraction. */
  bt =
    exp (lgamma (a + b) - lgamma (a) - lgamma (b) + a * log (x) +
	 b * log (1.0 - x));
  if (x < (a + 1.0) / (a + b + 2.0))	/* Use continued fraction directly. */
    return bt * betacf (a, b, x) / a;
  else				/* Use continued faction after making */
    return 1.0 - bt * betacf (b, a, 1.0 - x) / b;	/* the symmetry transformation. */
}

/*********************************************************************
   Continued fraction evaluation routine needed for the incomplete beta 
   function, I_x(a,b).	
   C.A. Bertulani        May/16/2000
*********************************************************************/

static double
betacf (double a, double b, double x)
/* Used by betai: Evaluates continued fraction for incomplete beta function 
   by modified Lentz's method.   */
{
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;			/* First step of Lentz's method.  */
  d = 1.0 - qab * x / qap;
  if (fabs (d) < FPMIN)
    d = FPMIN;
  d = 1.0 / d;
  h = d;
  for (m = 1; m <= MAXIT; m++)
    {
      m2 = 2 * m;
      aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1.0 + aa * d;		/* One step (the even one) of the recurrence. */
      if (fabs (d) < FPMIN)
	d = FPMIN;
      c = 1.0 + aa / c;
      if (fabs (c) < FPMIN)
	c = FPMIN;
      d = 1.0 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1.0 + aa * d;		/* Next step of the recurence the odd one) */
      if (fabs (d) < FPMIN)
	d = FPMIN;
      c = 1.0 + aa / c;
      if (fabs (c) < FPMIN)
	c = FPMIN;
      d = 1.0 / d;
      del = d * c;
      h *= del;
      if (fabs (del - 1.0) < EPS)
	break;			/* Are we done? */
    }
  if (m > MAXIT)
    fatalx ("a or b too big, or MAXIT too small in betacf\n");
  return h;
}

void
bpars (double *a, double *b, double mean, double var)
{

  double x2, g, xmean, x, m, v;
  double xa, xb, ym, yv;

  m = mean;
  v = var;
  x = (m * (1 - m) - v) / v;
  xa = x * m;
  xb = x * (1 - m);


/**
  ym = xa/(xa+xb) ;  
  yv = xa*xb/((xa+xb)*(xa+xb)*(xa+xb+1)) ;
  printf("%9.3f %9.3f\n", mean, ym) ;  
  printf("%9.3f %9.3f\n", var, yv) ;
*/

  *a = xa;
  *b = xb;

}

void
bmoments (double a, double b, double *mean, double *var)
{

  double x2, g, xmean, x, m, v;
  double xa, xb, ym, yv;

  x = a + b;
  *mean = a / x;
  *var = (a * b) / (x * x * (x + 1));

}

double
unbiasedest (int *ndx, int ndsize, int **counts)
{
/**
 ndx is ndsize array containing small integers coding pop index of each bracket (pop0 assumed)
 thus ndsize = 4   ndx = (1,1,1,1) codes (p_0-p_1)^4  
 thus ndsize = 4   ndx = (1,1,2,3) codes (p_0-p_1)^2  (p_0-p_2)  (p_0-p_3)

 counts [][] is integer array containing   counts[k][0] is count for variant allele for pop k
                                           counts[k][1] is count for reference allele for pop k

*/
  double xtop, xbot, yest, y;
  int popind[20];
  int popmax, j, k, n, nmax, a, t, s;
  int *tcounts;
  double **xmomest, yp;

  ivmaxmin (ndx, ndsize, &popmax, NULL);
  //printf("popmax: %d\n", popmax) ;
  ZALLOC (tcounts, popmax + 1, int);

  for (j = 0; j <= popmax; ++j)
    {
      tcounts[j] = counts[j][0] + counts[j][1];
    }

/** unbiased estimate of p_j^k */
  xmomest = initarray_2Ddouble (popmax + 1, ndsize, 0.0);
  for (j = 0; j <= popmax; ++j)
    {
      xmomest[j][0] = 1.0;
      for (k = 1; k <= ndsize; ++k)
	{
	  xtop = ifall (counts[j][0], k);
	  xbot = ifall (tcounts[j], k);
	  if (xbot <= 0.1)
	    xmomest[j][k] = -10000.0;
	  else
	    xmomest[j][k] = (double) xtop / (double) xbot;
	  //printf("zz %3d %3d %9.3f\n", j, k, xmomest[j][k] ) ;
	}
    }
  nmax = (1 << (ndsize)) - 1;
  yest = 0.0;
//printf("nmax: %d\n", nmax) ;

  for (n = 0; n <= nmax; ++n)
    {
      t = n;
      ivzero (popind, popmax + 1);
      for (k = 0; k < ndsize; ++k)
	{
	  a = 0;
	  s = t & 1;
	  t = t >> 1;
	  if (s == 1)
	    a = ndx[k];
	  ++popind[a];
	}
      yp = 1.0;
      for (j = 0; j <= popmax; ++j)
	{
	  t = popind[j];
	  s = 0;
	  if (j > 0)
	    s = t % 2;		// flags sign 
	  y = xmomest[j][t];

	  if (y < -1.0)
	    {
	      free (tcounts);
	      free2D (&xmomest, popmax + 1);
	      return (-10000.0);
	    }

	  if (s == 1)
	    y = -y;
	  yp *= y;
	}
      //printf(" %12.6f ", yp) ;
      //printimat(popind, 1, popmax+1) ; 
      yest += yp;
    }


  if (fabs (yest) >= 100)
    yest = -10000;


  free (tcounts);
  free2D (&xmomest, popmax + 1);
  return (yest);

}

void
weightjack (double *est, double *sig, double mean, double *jmean, double *jwt,
	    int g)
// test for jwt 0 
{
  double *jjmean, *jjwt;
  int i, n;

  ZALLOC (jjmean, g, double);
  ZALLOC (jjwt, g, double);
  n = 0;

  for (i = 0; i < g; ++i)
    {
      if (jwt[i] < 1.0e-6)
	continue;
      jjmean[n] = jmean[i];
      jjwt[n] = jwt[i];
      ++n;
    }

  weightjackx (est, sig, mean, jjmean, jjwt, n);
  free (jjmean);
  free (jjwt);
}

static void
weightjackx (double *est, double *sig, double mean, double *jmean,
	     double *jwt, int g)
// weighted jackknife see wjack.tex
// mean is natural estimate.  jmean[k] mean with block k removed.  jwt is weight for block (sample size)
{

  double *tdiff, *hh, *xtau, *w1, *w2;
  double jackest, yn, yvar;
  int k;

  if (g <= 1)
    fatalx ("(weightjack) number of blocks <= 1\n");
  ZALLOC (tdiff, g, double);
  ZALLOC (hh, g, double);
  ZALLOC (xtau, g, double);
  ZALLOC (w1, g, double);
  ZALLOC (w2, g, double);

  yn = asum (jwt, g);

  vsp (tdiff, jmean, -mean, g);
  vst (tdiff, tdiff, -1.0, g);
  jackest = asum (tdiff, g) + vdot (jwt, jmean, g) / yn;
// this is equation 2

  vclear (hh, yn, g);
  vvd (hh, hh, jwt, g);
/**
  for (k=0; k<g; ++k) {
   if (jwt[k] > 0.0) hh[k] /= jwt[k] ;  
   else hh[k] *= 1.0e20 ;
  }
*/
// jwt should be positive

  vst (xtau, hh, mean, g);
  vsp (w1, hh, -1.0, g);
  vvt (w2, w1, jmean, g);
  vvm (xtau, xtau, w2, g);

  vsp (xtau, xtau, -jackest, g);
  vvt (xtau, xtau, xtau, g);
  vvd (xtau, xtau, w1, g);
  yvar = asum (xtau, g) / (double) g;
  *est = jackest;
  *sig = 0.0;
  if (yvar > 0)
    *sig = sqrt (yvar);

  free (tdiff);
  free (hh);
  free (xtau);
  free (w1);
  free (w2);

}

int
modehprob (int n, int a, int m)
// mode of hypergeometric.  
{
  double v;
  int mode;

  v = (double) (a + 1) * (m + 1) / (double) (n + 2);
  mode = nnint (v - 0.5);

/** 
 v solves equation P(v) = P(v-1) 
 since P is log-concave it follows that mode is in 
 (v-1, v) 
*/

/**
 for (;;) { 
  if (loghprob(n, a, m, mode) < loghprob(n, a, m, mode+1)) { 
   mode = mode + 1 ;
   continue ;
  }
  if (loghprob(n, a, m, mode) < loghprob(n, a, m, mode-1)) { 
   mode = mode - 1 ;
   continue ;
  }
  break ; 
 }
*/

  return mode;

}



void
calcfc (double *c, int g, double rho)
// calculate spectral weights to correct correlation
{

  double *l, *d;
  int i, k;
  double pi;

  pi = 2.0 * acos (0.0);
  ZALLOC (l, g, double);
  ZALLOC (d, g, double);
  for (k = 0; k < g; k++)
    {
      l[k] = 1 + 2 * rho * cos ((2.0 * pi * k) / g);
      d[k] = 1.0 / sqrt (l[k]);
    }

  vzero (c, g);

  for (i = 0; i < g; i++)
    {
      c[i] = 0;
      for (k = 0; k < g; k++)
	{
	  c[i] += d[k] * cos ((2.0 * pi * k * i) / (double) g);
	}
    }
  vst (c, c, 1.0 / (double) g, g);

  free (l);
  free (d);

}

void
circconv (double *jp, double *c, double *jmean, int g)
{

  double *ww;
  int i, k;

  ZALLOC (ww, 2 * g, double);
  copyarr (c, ww, g);
  copyarr (c, ww + g, g);

  vzero (jp, g);
  for (i = 0; i < g; i++)
    {
      for (k = 0; k < g; k++)
	{
	  jp[i] += jmean[k] * ww[g + i - k];
	}
    }

  free (ww);
}

double
bino (int a, int b)
{
  if (b > a)
    return 0.0;
  if (a < 0)
    fatalx ("bad bino\n");
  if (b < 0)
    fatalx ("bad bino\n");
  if (a >= bcosize)
    fatalx ("bino overflow %d %d\n", a, bcosize - 1);

  return bcotable[a][b];


}

void
setbino (int maxbco)
{
  int i, j;
  if (maxbco < bcosize)
    return;
  if (bcosize >= 0)
    destroy_bino ();
  bcotable = initarray_2Ddouble (maxbco + 1, maxbco + 1, 0.0);
  bcosize = maxbco + 1;
  for (i = 0; i <= maxbco; ++i)
    {
      for (j = 0; j <= i; ++j)
	{
	  bcotable[i][j] = exp (logbino (i, j));
	}
    }
}

void
destroy_bino ()
{
  if (bcosize < 0)
    return;
  free2D (&bcotable, bcosize);
  bcosize = -1;
}

double scx(double *W, double *mean, double *x, int d) 
// v = x-mean compute v W v' 
// mean NULL ==> 0 
{
   double *ww, *w1  ;
   double y ;

   ZALLOC(w1, d, double) ; 
   ZALLOC(ww, d, double) ; 

   if (mean == NULL) { 
    copyarr(x, w1, d) ;
   }
   else {
    vvm(w1, x, mean, d) ;
   }
   mulmat(ww, W, w1, d, d, 1) ;
   y = vdot(ww, w1, d) ;

   free(w1) ; 
   free(ww) ;

   return y ; // -0.5 y is log likelihood

}



void
dither (double *xout, double *xin, int n)

/** add tiny random  xin  and xout can agree  */
{
  int i;
  double tail, y, top, bot;
  double *a;
  int *ind;

  ZALLOC (a, n, double);

  y = asum2 (xin, n) / (double) n;
  y = sqrt (y + 1.0e-10);
  y *= 1.0e-6;

  gaussa (a, n);
  vst (a, a, y, n);
  vvp (xout, a, xin, n);

  free (a);

}

void
probit (double *xout, double *xin, int n)

/** probit transform */
{
  int i;
  double tail, y, top, bot;
  double *a;
  int *ind;

  ZALLOC (a, n, double);
  ZALLOC (ind, n, int);

  vst (a, xin, -1.0, n);
  sortit (a, ind, n);
  bot = (double) (n + 1);
  for (i = 0; i < n; i++) {
    top = (double) (i + 1);
    tail = top / bot;
    y = zprob (tail);
    xout[ind[i]] = y;
  }

  free (a);
  free (ind);


}
double
exx (double x)
{
  double sum, term;
  int n;

  if (x < 0)
    return (1 - exp (-x)) / x;
  if (x == 0)
    return 1;                   // l'Hopital
  if (x > .01)
    return (1 - exp (-x)) / x;
  sum = term = 1;
  for (n = 2; n <= 100; ++n) {
    term *= (-x / (double) n);
    sum += term;
    if (fabs (term) < 1.0e-16)
      break;
  }
  return sum;

}

double
ubias (int a, int n, int k)
// a hits out of n
// return unbiased estimate of p^k where p is prob of hit
{

  int i;
  double u, v, y;

  if (k == 0)
    return 1;
  if (n < k)
    return -1;                  // impossible value should always check

  u = (double) a;
  v = (double) n;

  y = 1;
  for (i = 1; i <= k; ++i) {
    y = (y * u) / v;
    u -= 1.0;
    v -= 1.0;
  }
  return y;
}

