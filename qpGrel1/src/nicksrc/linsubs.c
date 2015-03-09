#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "vsubs.h"  
#include "strsubs.h" 
#include "linsubs.h" 
static int linsolvx(int nDim, double* pfMatr, double* pfVect, double* pfSolution) ;
static int calcdet = NO ; 

int setcalcdet (int val) 
{
  calcdet = val ;
}
void bal(double *a, double *b, int n) 
/** 
 normalize mean 0 s.d 1 
*/
{
    double t ;
    t = asum(b,n)/ (double) n ;
    vsp (a, b, -t, n) ;

    t = asum2(a,n)/ (double) n ;
    vst (a, a, 1.0/sqrt(t), n) ;
}

void mulmat(double *a, double *b, double *c, int a1, int a2, int a3) 
/* b is a1 x a2 , c a2 x a3 so a is a1 x a3  */  
{
    double *t ;
    int i,j,k ;
    ZALLOC(t, a1*a3, double) ;

    for (i=0; i<a1; i++)  
     for (j=0; j<a3; j++)  
      for (k=0; k<a2; k++)  
       t[i*a3+j] += b[i*a2+k]*c[k*a3+j] ;

    copyarr(t, a, a1*a3) ;

    free (t) ;

}
void imulmat(int *a, int *b, int *c, int a1, int a2, int a3) 
/* b is a1 x a2 , c a2 x a3 so a is a1 x a3  */  
{
    int *t ;
    int i,j,k ;
    ZALLOC(t, a1*a3, int) ;

    for (i=0; i<a1; i++)  
     for (j=0; j<a3; j++)  
      for (k=0; k<a2; k++)  
       t[i*a3+j] += b[i*a2+k]*c[k*a3+j] ;

    copyiarr(t, a, a1*a3) ;

    free (t) ;

}


double pdinv(double *cinv, double *coeff, int n) 
// cinv and coeff can be same
// cinv can be NULL 
// return log det (coeff) 
{
   double *tt;
   double *p ;
   double t, sum, y ;
   int i,j, k ;

/**
   pmat(coeff, n) ;
*/
   ZALLOC (tt, n*n, double);
   ZALLOC (p, n, double );
   

  copyarr(coeff,tt,n*n); 
  
  choldc (tt, n, p) ;

 
  for (i=0; i<n; i++) {
    tt[i*n+i] = 1.0/p[i] ;
    for (j=i+1; j<n; j++) {
      sum=0.0 ;
      for (k=i; k<j; k++) {
        sum -= tt[j*n+k]*tt[k*n+i] ;
      }
      tt[j*n+i] = sum/p[j] ;

    }
  }

   for (i=0; i<n; i++) 
    for (j=i; j<n; j++) {
     sum=0.0 ;
     if (cinv == NULL) break ;
     for (k=j; k<n; k++) {
      sum += tt[k*n+j]*tt[k*n+i] ;
     }
     cinv[i*n+j] = cinv[j*n+i] = sum ;
    }

    vlog(p, p, n) ; 
    y = 2.0*asum(p, n) ;


   free(tt) ;
   free(p) ;

   return y ;

}

int            
solvit (double *prod, double *rhs, int n, double *ans) 
{ 
  //The coefficient matrix should be positive definite

  /*AT : changed this code to take in matrix in a linear array form*/
  double *ttt;                   
  double *b; 
  double *p; 
  int i ; 
  int ret ;
 
 
  ZALLOC (ttt, n*n, double); 
  ZALLOC (p, n, double); 
  ZALLOC(b,n,double); 
 
  copyarr(prod,ttt,n*n); 
  copyarr(rhs,b,n); 
 
  ret = choldc (ttt, n, p); 
  if (ret<0) return -1 ;  // not pos def
  cholsl (ttt, n, p, b, ans); 
   
  free (ttt) ;  
  free(b);           
  free (p) ; 

//  printf("zzsol: %d %g %g %g\n", n, prod[0], rhs[0], ans[0]) ;
 
  return 1 ;
} 

void 
cholsl (double *a, int n, double p[], double b[], double x[])

/** 
 Numerical Recipes.  Must change 
*/

{
  /*AT: Changing the code*/

 int i, k; 
  double sum; 
     
 
  for (i = 0; i < n; i++) 
    { 
      sum = b[i]; 
      for (k = i - 1; k >= 0; k--) 
        sum -= a[i*n+k] * x[k]; 
      x[i] = sum / p[i]; 
    } 
 
  for (i = (n-1); i >= 0; i--) 
    { 
      sum = x[i]; 
      for (k = i + 1; k < n; k++) 
        sum -= a[k*n+i]* x[k]; 
      x[i] = sum / p[i]; 
    }        


}

int 
choldc (double *a, int n, double p[])
{
   int i, j,k;         
  double sum; 
     
  
  for (i = 0; i < n; i++) 
    {           
      for (j = i; j < n; j++) 
        {                 
          sum = a[i*n+j]; 
          for  (k = i - 1; k >= 0; k--) 
            sum -= a[i*n+k] * a[j*n+k]; 
          if (i == j) 
            {          
	      /**                                       
              printf("zzchol %d %20.10f %9.3f\n",i, sum, a[i][i]) ; 
	      */                                  
              if (sum <= 0.0) {                  
                return -1 ; // not pos def
              }  
              p[i] = sqrt (sum); 
            }   
          else  
            { 
              a[j*n+i] = sum / p[i]; 
               
            } 
        }                                    
    } 
  
    return 1 ;
                
}

void pmat(double *mat, int n)  

/** 
 print square matrix 
*/

{
 int i,j ;
 double *diag ;

 ZALLOC(diag,n, double) ;
 getdiag(diag, mat, n) ;

 printf("pmat:\n") ;
 
 for (i=0; i<n; i++) {
  printf("diag %5d %9.3f\n",i, diag[i]) ;
  for (j=0; j<n; j++) {
   if ((j%10) == 9)  printf("\n") ;
   if ((n%10) != 0) printf("%9.3f ",mat[i*n+j]) ;
  }
  printf("\n") ;
 }
 printf("\n") ;
 printf("\n") ;

 free(diag) ;
}

void cholesky(double *cf, double *a, int n) 
{  
  int i, j, k ;
  double *tt ;
  double *p ;
  
  ZALLOC(tt, n*n, double) ;
  ZALLOC(p, n, double) ; 
  copyarr(a,tt,n*n);


   choldc(tt, n, p ) ;

   vzero(cf, n*n) ;

   for (i = 0; i < n; i++) {
    tt[i*n+i] = p[i] ;
    
    for (j=0; j <= i ; j++) {  
     k = i*n+j ;
     cf[k] = tt[i*n+j] ;
    }
   }
  
   free(tt) ; 
   free(p) ;
}
//==============================================================================
//
// Linear System Solution by Gauss method
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================

int linsolv(int n, double* pfMatr, double* pfVect, double* sol)
// 1 on failure 
{

  int ret ; 
  double *a, *rhs ; 

  ZALLOC(a, n*n, double) ;
  ZALLOC(rhs, n, double) ;

  copyarr(pfMatr, a, n*n) ;
  copyarr(pfVect, rhs, n) ;

  ret = linsolvx(n, a, rhs, sol) ;

  free(a) ;
  free(rhs) ;


  return ret ;


}

int linsolvx(int nDim, double* pfMatr, double* pfVect, double* pfSolution)
{
  double fMaxElem;
  double fAcc;

  int i , j, k, m;


  for(k=0; k<(nDim-1); k++) // base row of matrix
  {
    // search of line with max element
    fMaxElem = fabs( pfMatr[k*nDim + k] );
    m = k;
    for(i=k+1; i<nDim; i++)
    {
      if(fMaxElem < fabs(pfMatr[i*nDim + k]) )
      {
        fMaxElem = pfMatr[i*nDim + k];
        m = i;
      }
    }
    
    // permutation of base line (index k) and max element line(index m)
    if(m != k)
    {
      for(i=k; i<nDim; i++)
      {
        fAcc               = pfMatr[k*nDim + i];
        pfMatr[k*nDim + i] = pfMatr[m*nDim + i];
        pfMatr[m*nDim + i] = fAcc;
      }
      fAcc = pfVect[k];
      pfVect[k] = pfVect[m];
      pfVect[m] = fAcc;
    }

    if( pfMatr[k*nDim + k] == 0.) return 1; // needs improvement !!!

    // triangulation of matrix with coefficients
    for(j=(k+1); j<nDim; j++) // current row of matrix
    {
      fAcc = - pfMatr[j*nDim + k] / pfMatr[k*nDim + k];
      for(i=k; i<nDim; i++)
      {
        pfMatr[j*nDim + i] = pfMatr[j*nDim + i] + fAcc*pfMatr[k*nDim + i];
      }
      pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
    }
  }

  for(k=(nDim-1); k>=0; k--)
  {
    pfSolution[k] = pfVect[k];
    for(i=(k+1); i<nDim; i++)
    {
      pfSolution[k] -= (pfMatr[k*nDim + i]*pfSolution[i]);
    }
    pfSolution[k] = pfSolution[k] / pfMatr[k*nDim + k];
  }

  return 0;
}
