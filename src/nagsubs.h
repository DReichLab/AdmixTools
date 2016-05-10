#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  

double psi(double x) ;
double tau(double x) ;
void mleg(double a1, double a2, double *p, double *lam) ;
double chisqtail(double x, double df) ; 

void eigx(double *evec, double *eval, 
 double *mat, int *nval, int n)   ; 

void linsolv(double *ans, double *mat, 
 double *rhs, int n)   ; 


#define YES  1
#define NO   0

#ifdef __cplusplus
}
#endif
