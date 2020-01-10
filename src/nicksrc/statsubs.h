#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  

#include "strsubs.h" 

#define	CHI_EPSILON     0.000001    /* accuracy of critchi approximation */
#define	CHI_MAX     99999.0         /* maximum chi square value */

#define	LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define	I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define	I_PI            0.3183098861837906715377675 /* 1 / pi  */
#define PI              3.1415926535897932384626434  
#define	BIGX           20.0         /* max value to represent exp (x) */
#define	ex(x)             (((x) < -BIGX) ? 0.0 : exp (x))
#define	SQRT_PI       (1.0/I_SQRT_PI)  /* sqrt (pi) */

double medchi(int *cls, int len, int *n0, int *n1, double *kstail)  ;
double ks2(int *cls, int len, int *n0, int *n1, double *kstail)  ;
double probks(double lam) ;

double nordis(double z) ;
double ndens(double val, double mean, double sig) ;
double ntail(double z) ;
void tailstats(double *x, double a, int isupper) ; 
double zprob(double p) ;
void setzptable() ;
double z2x2(double *a) ;
double conchi(double *a, int m, int n)  ;
double conchiv(double *a, int m, int n)  ;
double chitest(double *a, double *p, int n) ;
double pi() ;

double xlgamma(double x) ;
double psi(double x) ;
double tau(double x) ;
double logbessi0(double x) ;
double bessi0(double x) ;
double logbessi1(double x) ;
double bessi1(double x) ;
void bernload() ;
double bernum(int x) ;

void mleg(double a1, double a2, double *p, double *lam) ;

double dilog(double x) ;
double li2(double x) ;

double hwstat(double *x) ;

double gammprob(double x, double p, double lam)  ;
double bprob(double p, double a, double b)  ;
double lbeta(double a, double b)  ;
double poissloglike(int kk, double mean) ;
double poissexp(int kk, double mean) ;
double dirmult(double *pp, int *aa, int len) ;
double dawson(double t) ;

double binomtail(int n, int t, double p, char c) ;
double binlogtail(int n, int t, double p, char c) ;
void genbin(double *a, int n, double p) ;
void genlogbin(double *a, int n, double p) ;
int ifirstgt(int val, int *tab, int n)  ;
int firstgt(double val, double *tab, int n)  ;

void cinterp(double val, double x0, double x1, 
  double f0, double f0p, double f1, double f1p, double *fv, double *fvp)  ;
int firstgtx(double val, double *tab, int n) ;
int jfirstgtx(int val, int *tab, int n) ;

double rtlchsq(int df, double z) ;
double critchi(int df, double z) ;
double rtlf(int df1, int df2, double f) ;

double ltlg(double a, double x)  ;
double rtlg(double a, double x)  ;

double twdens(double twstat)  ;
double twtail(double twstat)  ;
double twtailx(double twstat)  ;
double twdensx(double twstat)  ;
double twnorm(double lam, double p, double n)  ; 
void twfree() ;
int settwxtable(char *table) ;
void gettw(double x,  double *tailp, double *densp)   ;
double dotwcalc(double *lambda, int m, double *ptw, double *pzn, double *pzvar, int minm) ;
int numgtz(double *a, int n) ;

double betaix(double a, double b, double lo, double hi)  ;
double betai(double a, double b, double x) ;
void bpars(double *a, double *b, double mean, double var)  ;
void bmoments(double a, double b, double *mean, double *var)  ;
double unbiasedest(int *ndx, int ndsize, int **counts) ;
void weightjack(double *est, double *sig, double mean, double *jmean, double *jwt, int g)  ;
int modehprob(int n, int a, int m) ;
void calcfc(double *c, int n, double rho) ;
void circconv(double *xout, double *xa, double *xb, int n) ;

double bino(int a, int b) ;
void setbino(int maxbco) ;
void destroy_bino() ; 
double exx(double x) ; 
double ubias(int a, int n, int k) ;  
double scx(double *W, double *mean, double *x, int d)  ; // maybe should be in vsubs  
void dither(double *xout, double *xin, int n)  ;
void probit(double *xout, double *xin, int n)  ;
// Berk-Jones 
double bjugauss(double *p, double *u, double *a, int n)  ; // return M_n  (Berk-Jones) 
void bjasympt(double *ptail, double *mtail, double *tail, double mplus, double mminus, int n)  ; 
void bj2(double *aa, double *bb, int a, int b,  double *plpv, double *prpv, double *ppv)  ;
double genhp(double **hp, int a, int b) ; 
double genhpt(int a, int b, int *lt, int *rt)  ; 
void gentail(double **ltail, double **rtail, double **hp, int a, int b)  ;
void setthresh(int *thresh, double **tail, int a, int b, double stat, int mode)  ; 
void bj2x(int *type, int a, int b, double *plpv, double *prpv, double *ppv)  ; 

void mlebeta(double *a, int n, double *p1, double *p2);
void estbpars(double *a, int n, double *p1, double *p2);
void mleb(double *p1, double *p2, double u, double v) ;

int  loadmptable(double ***mptable) ; 
double wynn(double *v, int n, double *acc, int *nacc) ;
double *vwynn(double **vv, int n, int dim, double **acc, int *nacc) ;

double rad2deg(double rad) ; 
double deg2rad(double deg) ; 

double quartile(double *x, int n, double q)  ;
int qinterp(double *a, double *b, int n, double val, double *ans)  ;
double truncexpmean(double m, double thresh, int isupper) ;
