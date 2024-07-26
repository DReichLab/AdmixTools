#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>
#include "qpsubs.h"

#define WVERSION   "100"
#define MAXSTR  512

extern int gsldetails;
extern double gslprecision;

double scorit (double *www, int n, double *pfix, double *ans);
static int xnmix;
static double *xvmix;
static double *www;

gsl_multimin_fminimizer *s = NULL;
static gsl_vector *ss, *x;
gsl_multimin_function minex_func;

double fjunk (const gsl_vector *v, void *params);
static double fff (double *a, int n);

int
gslsetup (int nmix, double *vmix)
{
  int k;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;

  if (nmix == 0)
    return 0;
  xnmix = nmix;
  ZALLOC (xvmix, nmix, double);
  copyarr (vmix, xvmix, nmix);
  ZALLOC (www, 2 * nmix, double);


  minex_func.n = nmix;
  minex_func.f = fjunk;
  minex_func.params = NULL;

  x = gsl_vector_alloc (nmix);
  for (k = 0; k < nmix; ++k) {
    gsl_vector_set (x, k, vmix[k]);
  }

  /* Set initial step sizes to 0.1 */
  ss = gsl_vector_alloc (nmix);
  gsl_vector_set_all (ss, 0.1);

  /* Initialize method and iterate */

  s = gsl_multimin_fminimizer_alloc (T, nmix);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  gsl_set_error_handler_off ();

  printf ("gslsetup called\n");
  return 1;
}

int gslnewvals(int nmix, double *vmix) 
{
 int k ; 

  for (k = 0; k < nmix; ++k) {
    gsl_vector_set (x, k, vmix[k]);
  }

 return 1 ;

}

double
gslopt (double *wpars)
{
  size_t iter = 0, k;
  int status;
  double size;
  double q = 999999;
  double qbest = 1.0e40;

  /* Starting point */


  for (k = 0; k < xnmix; ++k) {
    gsl_vector_set (x, k, wpars[k]);
  }
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  gsl_set_error_handler_off ();

  do {
    iter++;
    status = gsl_multimin_fminimizer_iterate (s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, gslprecision);
    q = s->fval;
    if (gsldetails)
      printf ("gslopt: %3d %12.6f\n", (int) iter, q);
    if (q < qbest) {
      if (gsldetails)
	printf ("+++ new best\n");
      qbest = q;
      for (k = 0; k < xnmix; ++k) {
	wpars[k] = gsl_vector_get (s->x, k);
      }
    }
  }
  while (status == GSL_CONTINUE && iter < 10000);

  printf ("gslans: %4d %12.6f\n", (int) iter, q);
  fff (wpars, xnmix);
  printmat (wpars, 1, xnmix);

  fflush (stdout);
  return q;
}

double
fff (double *a, int n)
// make a cononical return penalty 
{
  double penalty = 0;
  int k;
  double xx, y, yy;
  for (k = 0; k < xnmix; ++k) {
    yy = xx = a[k];
    if ((xx >= 0) && (xx <= 1))
      continue;
    if (xx < 0) {
      yy = -yy;
    }
    if (yy > 1) {
      yy = 2.0 * modf (0.5 * yy, &y);	// periodicity is 2 
      if (yy > 1)
	yy = 2 - yy;
    }
    y = yy - xx;
    penalty += y * y;
    a[k] = xx = yy;
  }
  return penalty;
}

double
fjunk (const gsl_vector * v, void *params)
{
  double xx;
  double q;
  int k, t = 0;
  double *tt;
  double penalty = 0;

  ZALLOC (tt, xnmix, double);
  for (k = 0; k < xnmix; ++k) {
    tt[k] = gsl_vector_get (v, k);
  }
  penalty = fff (tt, xnmix);

  t = 0;
  for (k = 0; k < xnmix; ++k) {
    xx = tt[k];
    www[t] = xx;
    www[t + 1] = 1.0 - xx;
    t += 2;
  }
  q = scorit (www, t, NULL, NULL);
  q += 10 * penalty;
  if (gsldetails) {
    printf ("zzfjunk %9.3f\n", q);
    printmat (tt, 1, xnmix);
  }
  free (tt);
  return q;
}
