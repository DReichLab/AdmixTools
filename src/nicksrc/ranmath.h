#include<math.h>
#include<stdlib.h>


#include <limits.h>

#define BIGINT INT_MAX   
#define SRAND  srandom
#define LRAND  random
#define DRAND() ( (double) (random() % BIGINT) / (double) (BIGINT)) 
#define DRAND2() ( drand2() ) 
/* random must return random integer in range 0 to BIGINT-1  */


#define NORMAL gauss

double  gauss()  ;  // standard normal
void    gaussa(double *a, int n)  ;  // array of standard normals
double  gds(double a)  ;  // obsolete
double  poidev(double mean) ;  // obsolete
double  ranpoiss(double mean) ;  // poisson Note double
double  ranpoissx(double mean) ;  // poisson | > 0
void  ranperm(int *a, int n) ;   // randomly permute a; if random permulation wanted : idperm(a,n) ; ranperm(a,n) 
double ranexp( void) ;  // exponential mean 1
double rangam(double a) ;  // standard gamma mean a
int randis(double *a, int n) ;  // element from discrete distribution a
void ransamp(int *samp, int nsamp, double *p, int plen) ; // sample nsamp samples from p
void pick2(int n, int *k1, int *k2) ;  // pick 2 elements from 0..n-1
int ranmod(int n)  ;  // random mod n 
int iranpick(int lo, int hi) ;  // random int in [lo, hi] 
double ranbeta(double a, double b) ;  //  beta 
int ranbinom(int n, double p) ;   // binomial  
void setrand(double *vv, int n) ;  // filll vv with U[0,1]
int ewens(int *a, int n, double theta) ;  // ewens sampling formula  
void genmultgauss(double *rvec, int num, int n, double *covar) ;  // multivariate 
double drand2() ;  
void ranmultinom(int *samp, int n, double *p, int len)  ;  // multinomial
double ranchi (int d)  ;   // chisq d dof.
void raninvwis(double *wis, int t, int d, double *s)  ;  // inverse wishart
double uniform(double lo, double hi) ;   // uniform (lo..hi)
void ransimplex(double *x, int n) ;   // uniform on n-simplex  
void randirichlet(double *x, double *pp, int n)  ;  // dirichlet parameter vector pp
void randirmult(double *pp, int *aa, int len, int m) ;  // dirichlet multinomial.  Output aa
int prob1(double p) ; // return YES with probability p 
double rant(double df) ;  // t distribution  
double samppow(double e, double a, double b) ;
double rejnorm(double lo, double hi) ;       // usually call ranboundnorm 
double ranboundnorm(double lo, double hi) ;  // sample standard normal in [lo, hi] 
double rtrunc2(double T) ;  // sample standard normal > T Rayleigh rejection         
double rantruncnorm(double T, int upper) ;  // sample standard normal > T (upper =1) < T (upper = 0) 
int ranhprob(int n, int a, int m) ;    // n balls a block sample m .  Fow many black
double rangeom (double theta) ;  // Geometric distribution, mean 1/theta
