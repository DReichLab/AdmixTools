#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_errno.h>

#include <nicklib.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  
#include "regsubs.h"  
#include "egsubs.h"  
#include "qpsubs.h" 


#define WVERSION   "100"
#define MAXSTR  512

extern int gsldetails;
extern double gslprecision;
extern int gslmaxiter ; 

double scorit(double *chol, double *xtri, double *xvar, int n) ; 
static int xnpars, vn;
static double *xvpars;
static double *www;
static double *vxtri, *vxvar ; 

gsl_multimin_fminimizer *s = NULL;
static gsl_vector *ss, *x;
gsl_multimin_function minex_func;

double fjunk (const gsl_vector *v, void *params);
double mkcanon (double *b, double *a, int n) ;

int
gslsetup (int n, double *vpars, double initscale, double *xtri, double *xvar) 
{
  int k;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  int dv, npars ; 

  dv = n*(n+1)/2 ;
  ZALLOC(vxtri, dv, double) ;
  ZALLOC(vxvar, dv*dv, double) ;

  copyarr(xtri, vxtri, dv) ;
  copyarr(xvar, vxvar, dv*dv) ;

  vn = n ;

  npars = dv - 1 ; // bottom corner of cholesky matrix is zero

  if (npars <= 0)
    return 0;
  xnpars = npars;
  ZALLOC (xvpars, npars, double);
  copyarr (vpars, xvpars, npars);

  minex_func.n = npars;
  minex_func.f = fjunk;
  minex_func.params = NULL;

  x = gsl_vector_alloc (npars);
  for (k = 0; k < npars; ++k) {
    gsl_vector_set (x, k, vpars[k]);
  }

  if (gsldetails) {
   printf("zzsetup ") ;  
   printnl() ;
   printmatl(vpars, 1, npars) ;
   printnl() ;
   printmatw(xtri, 1, dv, dv) ;
   printnl() ;
   printmatw(xvar, dv, dv, dv) ;
   printnl() ;
   printnl() ;
  } 

  /* Set initial step sizes to initscale */
  ss = gsl_vector_alloc (npars);
  gsl_vector_set_all (ss, initscale);

  /* Initialize method and iterate */

  s = gsl_multimin_fminimizer_alloc (T, npars);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

  gsl_set_error_handler_off ();

//  printf ("gslsetup called\n");
  return 1;
}

void gslfree() 
{

  gsl_vector_free(ss) ; 
  gsl_vector_free(x) ; 
  gsl_multimin_fminimizer_free(s) ; 
  free(vxtri) ;
  free(vxvar) ;

}


double
gslopt (double *wpars)
{
  size_t iter = 0, k;
  int status;
  double size;
  double q = 999999;
  double qbest = 1.0e40;

  // Starting point add space for final zero


  for (k = 0; k < xnpars; ++k) {
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

    if (q < qbest) {
      qbest = q;
      for (k = 0; k < xnpars; ++k) {
	wpars[k] = gsl_vector_get (s->x, k);
      }
      wpars[xnpars] = 0 ;
      if (gsldetails) { 
	printf ("+++ new best %12.6f: ", qbest);
        printmatwl(wpars, 1, xnpars, xnpars) ;
//      printf("zzheap: %x\n", topheap()) ; 
      }
    }
  }
  while (status == GSL_CONTINUE && iter < gslmaxiter);

  mkcanon (wpars, wpars, vn);
/**
  printf ("gslans: %4d %12.6f\n", (int) iter, q);
  printmatwl (wpars, 1, xnpars, xnpars);
*/

  fflush (stdout);
  return qbest ; // 2 * log likelihood -- chi sq (1 d of freedom) 
}

double
fjunk (const gsl_vector * v, void *params)
{
  double xx;
  double q;
  int k, t = 0;
  double *tt;
  double penalty = 0;

  ZALLOC (tt, xnpars+1, double); // pass in last 0 
  for (k = 0; k < xnpars; ++k) {
    tt[k] = gsl_vector_get (v, k);
  }
  penalty = mkcanon (tt, tt, vn);

  q = scorit (tt, vxtri, vxvar, vn);
  q += 10 * penalty;
  if (gsldetails) {
    printf ("zzfjunk %12.6f: ", q);
    printmatl (tt, 1, xnpars);
  }
  free (tt);
  return q;
}

void setgsldetails(int val) 
{
 gsldetails = val ;
}

int getgsldetails() 
{
  return gsldetails ;

}
