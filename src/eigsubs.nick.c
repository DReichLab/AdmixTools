#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h> 
#include "eigsubs.h" 
/* ********************************************************************* */
void eigx_(double *pmat, double *ev, int *n) ;
void eigxv_(double *pmat, double *eval, double *evec, int *n) ;
void cdc_(double *pmat, int *n) ;
void inverse_(double *pmat, int *n);
void solve_(double *pmat, double *v, int *n) ;
void geneigsolve_(double *pmat, double *qmat, double *eval, int *n);

void packsym(double *pmat, double *mat, int n) ;


void eigvals(double *mat, double *evals, int n) 
{
	 double *pmat ;  
         int len ; 

         len = n*(n+1) ; len /= 2 ;
	 ZALLOC(pmat, len, double) ;

	 vst(mat, mat, -1.0, n*n) ;
	 packsym(pmat, mat, n) ;
         eigx_(pmat, evals, &n) ;
	 free(pmat) ;
	 vst(mat, mat, -1.0, n*n) ;
	 vst(evals, evals, -1.0, n) ;
}
void eigvecs(double *mat, double *evals, double *evecs, int n) 
{
	 double *pmat ;  
         int len ; 

         len = n*(n+1) ; len /= 2 ;
	 ZALLOC(pmat, len, double) ;

	 vst(mat, mat, -1.0, n*n) ;
	 packsym(pmat, mat, n) ;

         eigxv_(pmat, evals, evecs, &n) ;
	 free(pmat) ;
	 vst(mat, mat, -1.0, n*n) ;
	 vst(evals, evals, -1.0, n) ;
}

/*  note:  dpotrf requires the entire matrix, not packed lower-tri */
void chdecomp(double *mat, int n)  
{
         /* symetric matrix - don't need to 
	  * convert to column major order */

	 cdc_(mat, &n);
}

void inverse(double *mat, int n)  
{
         int i,j;

	 /* convert to column-major order */
	 for(i=0;i<n;i++)  {
	   for(j=0;j<i;j++)  {
	     double t = mat[n*i+j];
	     mat[n*i+j] = mat[n*j+i];
	     mat[n*j+i] = t;
           }
         }

         inverse_(mat, &n);

	 for(i=0;i<n;i++)  {
	   for(j=0;j<i;j++)  {
	     double t = mat[n*i+j];
	     mat[n*i+j] = mat[n*j+i];
	     mat[n*j+i] = t;
           }
         }
}

void solve(double *mat, double *b, double *v, int n)  
{
         int i,j;

	 double *mat2  ;
	 
         ZALLOC(mat2, n*n, double) ;

         /* lapack is going to put the lu-decomp into the matrix,
	  * so make a copy and convert to column-major order */
	 for(i=0;i<n;i++)  {
	   for(j=0;j<n;j++)  {
	     mat2[n*i+j] = mat[n*j+i];
           }
         }
	 printmat(mat, n, n) ;
	 printnl() ;

	 printmat(mat2, n, n) ; 
	 fflush(stdout) ;

	 /* copy b into v */
	 copyarr(b, v, n) ;

	 solvit(mat2, b, n, v) ;
//	 solve_(mat2, v, &n);

	 free(mat2);
	 return;
}
   
void
packsym(double *pmat, double *mat, int n) 
	//  lapack L mode (fortran)
{ 
	int i, j, k = 0 ;
	for (i=0; i<n; i++)  {  
          for (j=i; j<n; j++) { 
             pmat[k] = mat[i*n+j] ;
	     ++k ;
	  }
	}
}

void geneigsolve(double *pmat, double *qmat, double *evec, double *eval, int n)  {

  /* save copy of A and B, which LAPACK will overwrite */
  double *amat = (double *)malloc(n*n*sizeof(double));
  double *bmat = (double *)malloc(n*n*sizeof(double));

  int i, j;
  for(i=0;i<n*n;i++)  {
    amat[i] = pmat[i];
    bmat[i] = qmat[i];
  }

  /* matrices have to be symetric-definite, so don't need to convert to column-major order */
  geneigsolve_(pmat, qmat, eval, &n);

  /* copy eigenvectors to A and original A,B back */ 
  for(i=0;i<n;i++)  {
    for(j=0;j<n;j++)  {
      evec[i*n+j] = pmat[j*n+i];    /* back in row-major order */
    }
  }
  for(i=0;i<n*n;i++)  {
    pmat[i] = amat[i];
    qmat[i] = bmat[i];
  }

  free(amat);
  free(bmat);

}



