#include "empr.h"  

extern int asiforce ; 

double emest(double *vv, double *mean, double **vdata, int ndim, int ndata, int niter, double **postest) 
/**  
 EM;  mm, vv should be set;  vdata is obs and upper triangle of noise covar
 return est mean and covariance of global dist. Aso mean and covar (upper triangle) of posterior
 postest can be NULL
*/
{

  double mm[10], var[100], mpost[10], vpost[100], moment[200], la[10], lb[10], mobs[10], vobs[100] ;  
  double *vd, s1, s10, ylike, ybase, yl, ybot, ytop, y ; 
  double qvar[100],  qobs[100], A, B, C, z=0, ww[200] ; 
  double aa, bb ;

  int iter, k, a, b, kk ;  
  int nmoment ;

   copyarr(mean, mm, ndim) ; 
   copyarr(vv, var, ndim*ndim) ;
   A = B = 0 ;
  
  printf("starting iterations\n") ;
  fflush(stdout) ;

  for (iter = 1; iter <=niter; ++iter) { 
   ylike = 0 ;
   vzero(moment, 200) ;
   for (k=0; k<ndata; ++k) { 
    ybot = v2q(var, mm, &z, qvar, la, &A, ndim) ;
    vd = vdata[k] ; 
    copyarr(vd, mobs, ndim) ; 
    kk = ndim ;
    for (a=0; a<ndim; ++a) { 
     for (b=a; b<ndim; ++b) { 
      vobs[a*ndim+b] = vd[kk] ; 
      vobs[b*ndim+a] = vd[kk] ; 
      ++kk ;
    }}
    ybot += v2q(vobs, mobs, &z, qobs, lb, &B, ndim) ;
    vvp(qobs, qobs, qvar, ndim*ndim) ;
    vvp(lb, lb, la, ndim) ;
    B += A ; 
    ytop = q2v(vobs, mobs, &C, qobs, lb, &B, ndim) ;
    ylike += (ytop-ybot) ;
    nmoment = v2m(ww, vobs, mobs, ndim) ;
    vvp(moment, moment, ww, nmoment) ;
//  if (iter==1) printf("zz: %d ", k) ; printmatl(moment, 1, 6) ;

    if (postest == NULL) continue ;
    vd = postest[k] ; 
    copyarr(mobs, vd, ndim) ; 

    vd = vd + ndim ; 
    mktriang(vd, vobs, ndim) ; 
   }

   if (iter == 1) ybase = ylike ; 
   ylike -= ybase ; 
// now reestimate 
   m2v(moment, var, mm, ndim) ; 
   printf("iter: %d %12.3f\n", iter, ylike) ; 
   printmat(mm, 1, ndim) ;  printnl() ;
   printmat(var, ndim, ndim) ;  printnl() ;
   fflush(stdout) ;
  }
  copyarr(mm, mean, ndim) ; 
  copyarr(var, vv, ndim*ndim) ; 
  return ybase + ylike ;
}
void m2v(double *mom, double *var, double *mean, int n) 
{
   int x, xx ;
   double *tt, *t2, *pt, y ; 
   double fpars[4], mm[6]  ; 

   ZALLOC(tt, n*n, double) ; 
   ZALLOC(t2, n*n, double) ; 

   x = n*(n+1)/2 ;
   xx = x + n + 1 ;
   y = mom[0] ; 
   if (y<.001) fatalx("no odata in m2v\n") ; 
   vst(mom, mom, 1.0/y, xx) ;  

   pt = mom+1 ; 
   if (asiforce) { 
    if (n!=2) fatalx("bad dimension") ;
    mtof(fpars, mom) ;  
    ftom(fpars, mm) ; 
    pt = mm+1 ;
   }

   copyarr(pt, mean, n) ; pt += n ;
   addouter(t2, mean, n) ;  
  
   mkfull(tt, pt, n) ;
   vvm(var, tt, t2, n*n) ;
   

   free(tt) ; 
   free(t2) ;


}

int v2m(double *mom, double *var, double *mean, int n) 
{
   int x ;
   double *tt, *tri, *pt ; 

   ZALLOC(tt, n*n, double) ; 
   copyarr(var, tt, n*n) ; 
   addouter(tt, mean, n) ;  
  
   x = n*(n+1)/2 ;
   ZALLOC(tri, x, double) ;
  
   pt = mom ;
   pt[0] = 1 ;  pt += 1 ;
   copyarr(mean, pt, n) ; 
   pt += n ; 
   mktriang(tri, tt, n) ; 
   copyarr(tri, pt, x) ; 

   free(tt) ; 
   free(tri) ;

   return x + n + 1 ;

}

double q2v(double *V, double *M, double *A, double *Q, double *L, double *B, int n)  
// Input different from output

{

  pdinv(V, Q, n) ; 
  solvit(Q, L, n, M) ; 
  *A = *B - vdot(L, M, n) ; 
  return lognormint(Q, L,  *B, n) ;

}

double v2q(double *V, double *M, double *A, double *Q, double *L, double *B, int n)  
// Input different from output

{

  pdinv(Q, V, n) ; 
  mulmat(L, Q, M, n, n, 1) ; 

  if (B==NULL) return 0 ; 

  *B = *A + vdot(L, M, n) ; 
  return lognormint(Q, L,  *B, n) ;

}

double rannormal(double m, double sig) 
{

  return sig*gauss() + m ; 

}
double lognormint(double *Q, double *L, double C, int n) 
// int exp -2 \sum_Q_{ij} x_i x_j -2 \sum L_i x_i + C 
{

  double y, yc, yl, yn, yldet ;
  double *mean ; 
  static double logtwopi = -1  ; 

  if (logtwopi < 0) { 
   y = pi() ; 
//  printf("pi: %15.10f", y) ;
   logtwopi = log(2.0*y) ;
  }

  ZALLOC(mean, n, double) ;
  solvit(Q, L, n, mean) ; 
  yc = C - vdot(mean, L, n) ; 
  yl = -0.5*yc ; 
  yn = (double) n ; 

  yl += (yn/2)*logtwopi ;  
  yldet = pdinv(NULL, Q, n) ; 
  yl -= -0.5 * yldet ; 
  
  free(mean) ;
  return yl ; 

}
void 
multgauss (double *rvec, int n, double *mean, double *covar)
// 1 sample
{
  genmultgauss(rvec, 1, n, covar) ; 
  vvp(rvec, rvec, mean, n) ;


}
void ftom(double *fpars, double *mpars)  
{
  double mm[6], s2[2], m0, m1, alpha, th1, y2, sigma[2] ; 

  th1 = fpars[0] ; 
  alpha = fpars[1] ; 
  copyarr(fpars+2, sigma, 2) ;
  mm[0] = 1 ; 
  
  vvt(s2, sigma, sigma, 2) ;
  m1 = mm[2] = th1 ;   
  m0 = mm[1] = -alpha + alpha*m1 ; 
  mm[5] = s2[1] ; 
  mm[3] = s2[0] ; 
  mm[4] = 0 ; 
  copyarr(mm, mpars, 6) ; 
  y2 = mpars[5] +=  m1*m1 ; 
  mpars[4] +=  -alpha*m1 + alpha*y2  ; 
  mpars[3] +=  alpha*alpha * (1 +  y2 - 2*m1 ) ;
}

void mtof(double *fpars, double *mpars)  
{
  double mm[6], s2[2], m0, m1, alpha, th1 ; 
  double scal, a1, a2, b1, b2  ; 

  copyarr(mpars, mm, 6) ;
  scal = mm[0] ; 
  vst(mm, mm, 1.0/scal, 6) ;

  m0 = mm[1] ; 
  m1 = mm[2] ; 

  a1 = mm[4] - m0 ; 
  a2 = 1 + mm[5] - 2*m1 ; 
  alpha = fpars[1]  = a1/a2 ; 

  b1 = mm[3] - 2 * alpha * a1 + alpha*alpha*a2 ;  
  
  s2[1]  = mm[5]  - m1*m1 ; 
  s2[0]  = b1 ; 
  printmat(fpars, 1, 2) ; 
  printmat(s2, 1, 2) ; 
  vsqrt(fpars+2, s2, 2) ;

}

