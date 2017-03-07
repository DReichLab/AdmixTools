#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

void bal(double *a, double *b, int n) ;

/* linear algebra */
void mulmat(double *a, double *b, double *c, int a1, int a2, int a3) ;
int solvit (double *prod, double *rhs,int n, double *ans);
int  solvitfix (double *prod, double *rhs, int n, double *ans, int *vfix, double *vvals, int nfix) ;
int  oldsolvitfix (double *prod, double *rhs, int n, double *ans, int *vfix, double *vvals, int nfix) ;
double pdinv(double *cinv, double *coeff, int n) ;

/* numer recipes p 97 */
double logdet(double *mat, int n) ;
int choldc (double *a, int n, double p[]);
void cholsl (double *a, int n, double p[], double b[], double x[]);
void cholesky(double *cf, double *a, int n) ;
void pmat(double *mat, int n)   ;
void imulmat(int *a, int *b, int *c, int a1, int a2, int a3) ;
int linsolv(int n, double* pfMatr, double* pfVect, double* sol) ; // Developer: Henry Guennadi Levkin

double qval(double *vv, double *q, double *l, int n) ; 
void qgrad(double *grad, double *vv, double *q, double *l, int n) ;
double mquad(double y0, double y1, double y2, double *pmx) ;
double qminpos(double *vv, double *q, double *l, int n)  ;
double qminposfix(double *vv, double *q, double *l, int n, int *fixlist, double *fixvals, int nfix)  ;
double qminposfixc(double *vv, double *q, double *l, int n, int *fixlist, double *fixvals, int nfix, int *constraint)  ;
double qmin(double *vv, double *q, double *l, int n) ;
double qminfix(double *vv, double *q, double *l, int n, int *fixlist, double *fixvals, int nfix)  ;;
double qmpc (double *vnew, double *vold, double *q, double *l, int *dead, int level, int *constraint, int n) ;
double qmp (double *vnew, double *vold, double *q, double *l, int *dead, int level, int n) ; 
