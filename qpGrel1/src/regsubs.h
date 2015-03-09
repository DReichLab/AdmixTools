#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <nicklib.h>  

double regressit(double *ans, double *eq, double *rhs, int m, int n) ;
void regressitall(char **vname, double *eq, double *rhs, int m, int n) ;
void add1(int *a, int *b, int n)  ;

void ptoz(double *p, double *z, int n) ;
void ztop(double *p, double *z, int n) ;
double logregressit(double *ans, double *eq, double **rhs, int neq, int nv) ;
double logrscore(double *eq, double **rhs, int neq, int nv) ;

void calcgh(double *grad, double *hess, double *eq, double *z, 
 double *n0, double *n1, int neq, int nv) ;

double zlike(double *eq, double *n0, double *n1, 
  double *ans, int neq, int nv)  ;

void squish(double *xmat, double *mat, int nrow, int oldc, int newc)  ;

void
calcres(double *res, double *ans, double *eq, double *rhs, 
 int  neq, int nv) ;

#ifdef __cplusplus
}
#endif
