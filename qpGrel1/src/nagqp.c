#include <nag.h>
#include <nag_stdlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <nage04.h>
#include <nicklib.h>  

extern int verbose ;


double scorit(double *www, int n, double *pfix, double *ans) ;
static int windex, pplen ;
static double *ppp, *ppbest, scbest ;
void normvec(double *www, int n) ;

double nagopt1(double lo, double hi, double *plam,  int wind);
void nagxfun(double xp, double *nagxval, Nag_Comm *comm) ;

void myfun(Integer n, double *xc, double *fc, Nag_Comm *comm) 
{  
   int nparams = n ;
   double *params ; 
   double y, yval, yfix ;
   int i  ;

   ZALLOC(params, nparams, double) ;
   copyarr(xc, params, n) ;

   yval = scorit(params, nparams, &yfix, NULL) ;                                         

   *fc = yval + 1000*yfix ;
 //printf("zz %15.9f %15.9f\n", *fc, yfix) ;

   free(params) ;

}
void simpest(double *params, int n) 
{
  double oldval, lam, val ; 
  int iter, maxiter = 50, t, k, ind, miniter=5 ;
  double xlimit = 1.0 ; // insist on this much improvement
  double xdiff, tmp ;

  ZALLOC(ppp, n, double) ; 
  ZALLOC(ppbest, n, double) ; 
  pplen = n ;

  copyarr(params, ppp, n) ;
  copyarr(params, ppbest, n) ;

  scbest = scorit(ppp, n, NULL, NULL) ;
  if (verbose) 
  printf("simpest: %4d %9.3f\n", 0, scbest) ;

  for (iter = 1; iter <=maxiter ; ++iter) { 
   oldval = scbest ;
   for (t=0; t<n; ++t)  {
/**
    getvmind(t, &k, &ind) ; 
    if (ind == 0) continue ;
*/ 
    lam = ppbest[t] ;
    val = nagopt1(0.0, 5.0, &lam, t) ;
    if (val<scbest) {
     ppbest[t] = lam ;
     normvec(ppbest, n) ;
     scbest = val ;
    }
   }
  if (verbose) 
   printf("simpest: %4d %9.3f\n", iter, scbest) ;
   xdiff = oldval - scbest ;
   if ((iter>=miniter) && (xdiff<xlimit)) break ;
  }

  copyarr(ppbest, params, n) ;
  normvec(params, n) ;
  val = scorit(params, n, NULL,  NULL) ;
  if (verbose) 
  printf("simpest final: %4d %9.3f\n", iter, val) ;


  free(ppp) ; 
  free(ppbest) ;

}

void nagopt(double *params, int n) 
{

  double val ;
  Nag_E04_Opt options;
  static NagError fail;
  double *ycoeffs ;
  static int ncall = 0 ; 
  int nparams = n, i ;
  double yval, oldyval, tmp ;

  if (n==0) return ;
  simpest(params, n) ;

  ++ncall ;
  e04xxc(&options);

  options.print_level = Nag_NoPrint;
  if (verbose) 
  options.print_level = Nag_Soln ;
  options.optim_tol = .00001 ;
  options.max_iter = n*200 ;

  options.list = FALSE ; 

  fail.print = TRUE;
  verbose = NO ;

  ZALLOC(ycoeffs, nparams, double) ;
  copyarr(params, ycoeffs, nparams) ;

  oldyval = yval = scorit(params, nparams, NULL, NULL) ;                                           

/**
  printf("nagopt ncall: %d  value:  %9.3f\n", ncall, yval) ;

  for (i=0; i<n; i++) { 
   printf("%3d %12.6fn", i, params[i]) ;  
  }
*/

  e04ccc(n, myfun, ycoeffs, &yval, &options, NAGCOMM_NULL, &fail);
  copyarr(ycoeffs, params, nparams) ;
  normvec(params, nparams) ;
  yval = scorit(params, nparams, &tmp, NULL) ;                                           
//  printf("nagopt ncall: %d  oldvalue:  %9.3f newvalue: %9.3f\n", ncall, oldyval, yval) ;


  free(ycoeffs) ;
}


double nagopt1(double lo, double hi, double *plam,  int wind)
{
    double e1, e2, minp, minval ;  
    double tlo, thi ;
    static NagError fail ;
    Nag_Comm comm ;
    
    Integer max_fun = 100 ;

    windex = wind ;
    copyarr(ppbest, ppp, pplen) ;

    tlo = lo ; 
    thi = hi ;
    minp = *plam ;
    e1 = 1.0e-5 ;
    e2 = 0.0 ;

    SET_FAIL(fail) ;
    fail.print = Nag_TRUE ;

    e04abc(nagxfun, e1, e2, &tlo, &thi, 
     max_fun, &minp, &minval, 
     &comm, &fail) ;

    if (fail.code != NE_NOERROR) { 
     fatalx("(e04abc) error\n") ;
    }

    *plam = minp ;
    return minval ;
}

void nagxfun(double xp, double *nagxval, Nag_Comm *comm) 
{
     double y ; 

     ppp[windex] = xp ;
     y = scorit(ppp, pplen, NULL, NULL) ;
//     printf("nagxfun: %9.3f %9.3f\n", xp, y) ;
     *nagxval =  y ;

}


