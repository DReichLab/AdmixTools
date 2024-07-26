#include <fcntl.h>
#include <ctype.h>
#include  <nicklib.h>
#include  "eigsubs.h"  


int gslsetup (int n, double *vpars, double initscale, double *xtri, double *xvar)  ; 
double gslopt(double *pars) ;
int initpars(double *pars, double *z, double *zv, int dim) ; 

void full2chol(double *mat, double *chol, int n) ;
void chol2full(double *mat, double *chol, int n) ; 
void chol2prod (double *prod, double *chol, int n) ; 
double rank1opt(double *ff, double *z, double *zv, int dim) ;
void solveco(double *coeff, double *mat, int dim)  ;
void initr1(double *outmat, double *mat, int n) ;
void xcheck(int dim) ; 

double pscale = 0.01 ;

extern int gsldetails ;
extern double gslprecision ;

static double *initvpars = NULL ; 


void destroypars() { 

 free(initvpars) ; 
 initvpars = NULL ;

}

void setpars(double *vv, int dim) {  
 int dv ; 

 dv = dim*(dim+1)/2 ; 
 if (initvpars != NULL) destroypars() ;
 ZALLOC(initvpars, dv, double) ;
 copyarr(vv, initvpars, dv) ;
}


double rank1opt(double *vv, double *z, double *zv, int dim) 
{

 int dv, t, s  ; 
 double **xx ; 
 int k, a, b ; 
 double *vpars,  *ww, *ww2 ;
 double *fmat ;
 double ychi ; 

 dv = dim*(dim+1)/2 ; 
 printmatw(z, 1, dv, dv) ; 
 printnl() ;
 printmatw(zv, dv, dv, dv) ; 

 ZALLOC(fmat, dim*dim, double) ;
 ZALLOC(vpars, dv, double) ;

 mkfull(fmat, z, dim) ; 
 printnl() ; 
 printmat(fmat, dim, dim) ;

 if (initvpars == NULL) initr1(vpars,  fmat, dim) ;
 else copyarr(initvpars, vpars, dv) ;
 printnl() ; 

 gslsetup(dim, vpars, pscale, z, zv) ;    

 ychi = gslopt(vpars) ; 
 printf("final score: %9.3f\n", ychi) ; 

 copyarr(vpars, vv, dv) ;

  printf("zzfmat:\n") ;
  chol2prod(fmat, vpars, dim) ; 
  printmatwl(fmat, dim, dim, dim) ; 


 free(fmat) ;

 return ychi ; 

}

int initpars(double *pars, double *z, double *zv, int dim) 
{
  double *ff, y ;  
  int x, a, b, dv ; 

  dv = dim*(dim+1)/2 ;
  vzero(pars, dv) ; 
  x = 0 ;
  for (a=0; a<dim; ++a) { 
   for (b=a; b<dim; ++b) { 
     pars[x] = 0 ; 
     if (a==b) { 
      y = fabs(z[x]) ; 
      pars[x] = sqrt(y) ; 
     }
     ++x ; 
   } 
  }
//  printf("zzpars ") ; printmat(pars, 1, dv) ;

  return 1 ;
}

void full2chol(double *mat, double *chol, int n) 
{
  int a, b, x ;

  x = 0 ; 
  for (a=0; a<n; ++a) { 
   for (b=a; b<n; ++b) { 
    chol[x] = mat[a*n+b] ; 
    ++x ;  
   }
  }
}
void chol2full(double *mat, double *chol, int n) 
{
  int a, b, x ;
  
  vzero(mat, n*n) ;
  x = 0 ; 
  for (a=0; a<n; ++a) { 
   for (b=a; b<n; ++b) { 
    mat[a*n+b] = chol[x] ; 
    ++x ;  
   }
  }
}

void chol2prod (double *prod, double *chol, int n) 

{

 double *mat ; 
 ZALLOC(mat, n*n, double) ;
 chol2full(mat, chol, n) ;
 txmulx(prod, mat, n, n) ;

 free(mat) ;

}

void prod2chol(double *prod, double *chol, int n) 
{
 
 double *ww  ; 
 ZALLOC(ww, n*n, double) ;

 cholesky(ww,  prod, n) ; 
 transpose(ww,  ww, n, n) ; 
 full2chol(ww, chol, n) ;

 free(ww) ;

}

double mkcanon(double *chol2, double *chol, int n) 
{

 double *full, *vv ; 
 double penalty = 0 ;
 int k ;
 double y ; 

 ZALLOC(full, n*n, double) ;
 chol2full(full, chol, n) ;
 
 for (k=0; k<n; ++k) { 
  y = full[k*n+k] ;  
  if (y<0) { 
   vv = full + k*n ; 
   vst(vv, vv, -1, n) ;
   penalty  -=  y ; 
  }
 }
 full2chol(full, chol2, n) ; 


 free(full) ; 
 return penalty ; 

}

double scorit(double *chol, double *xtri, double *xvar, int n) 
{

  int dv = n*(n+1)/2 ; 
  double *full, *tri ; 
  double y ; 

  ZALLOC(full, n*n, double) ; 
  ZALLOC(tri, dv, double) ;

  chol2prod(full, chol, n) ; 
  mktriang(tri, full, n) ; 
  y = scx(xvar, xtri, tri, dv) ; 

  free(full) ; 
  free(tri) ; 

  return y ; 
}

void solveco(double *coeff, double *mat, int dim) 
// corank1 a matrix ; we do some scaling here
{

 double *wco, *ww, *rhs, *r2 ; 

 ZALLOC(ww, dim*(dim+1), double) ; 
 ZALLOC(wco, dim*dim, double) ; 
 ZALLOC(r2, dim, double) ; 
 ZALLOC(rhs, dim+1, double) ; 

 copyarr(mat, ww, dim*dim) ; 
 vclear(ww+dim*dim, 1.0, dim) ; 
 rhs[dim] = 1.0 ; 
 txmulx(wco, ww, dim+1, dim) ; 
 mulmat(r2, rhs, ww, 1, dim+1, dim) ; 

 solvit(wco, r2, dim, coeff) ;
 printf("zzz\n") ;
 printmat(mat, dim, dim) ; 
 printnl() ;
 printmat(ww, dim+1, dim) ; 
 printnl() ;
 printmat(rhs, 1, dim+1) ; 
 printnl() ;
 printmat(wco, dim, dim) ; 
 printnl() ;
 printmat(r2, 1, dim) ; 

}


void xcheck(int dim) 
{

  int dv ;
  double *ff, *ch ;
  
  ZALLOC(ff, dim*dim, double) ; 
  dv = dim*(dim+1)/2 ;
  dv = dim*(dim+1)/2 ;
  ZALLOC(ch, dv, double) ; 

  setidmat(ff, dim) ; 
  full2chol(ff, ch, dim) ; 
  printf("check") ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  printnl() ;
  printmatw(ch, 1, dv, dv) ; 
  chol2full(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  chol2prod(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  vst(ch, ch, 2, dv) ; 
  chol2full(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  chol2prod(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 

 free(ff) ; 
 free(ch) ;
}

void initr1(double *outmat, double *mat, int n) 
{
   double *evecs, *evals, *full, *ww ; 
   int k, dim ;
  
   dim = n*(n+1)/2 ;
   ZALLOC(evecs, n*n, double) ;
   ZALLOC(evals, n, double) ;
   ZALLOC(full, n*n, double) ;
   ZALLOC(ww, n*n, double) ;

  eigvecs(mat, evals, evecs,  n) ; 
  if (gsldetails) {
   printf("mat: ") ; 
   printnl() ;
   printmat(mat, n, n) ;

   printf("evals: ") ; 
   printmat(evals, 1, n) ;
   printnl() ;
   printmat(evecs, n, n) ;
   printnl() ;
  }

   evals[n-1] = .0001 * asum(evals, n) ; 
   for (k=0; k<n; ++k) {  
    addoutmul(full, evecs+k*n, evals[k], n) ;
   }
  if (gsldetails) {
   printf("full: ") ; 
   printnl() ;
   printmat(full, n, n) ;
   printnl() ;
  }
   vzero(ww, n*n) ;
   prod2chol(full, ww, n) ; 
// printmatw(ww, n, n, n) ;
   copyarr(ww, outmat, dim) ;
   chol2prod(full, ww, n) ;

  if (gsldetails) {
   printf("check: ") ; 
   printnl() ;
   printmat(full, n, n) ;
   fflush(stdout) ;
  }
      
   free(evecs) ;
   free(evals) ;
   free(full) ;
   free(ww) ;
}
