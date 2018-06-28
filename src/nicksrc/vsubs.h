#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "strsubs.h" 

void vsp(double *a, double *b, double c, int n);
void vst(double *a, double *b, double c, int n);
void vvt(double *a, double *b, double *c, int n);
void vvp(double *a, double *b, double *c, int n);
void vvm(double *a, double *b, double *c, int n);
void vvd(double *a, double *b, double *c, int n);
void vsqrt(double *a, double *b,  int n) ;
void vinvert(double *a, double *b,  int n) ;
void vabs(double *a, double *b,  int n)  ;
void vlog(double *a, double *b,  int n)  ;
void vlog2(double *a, double *b,  int n)  ;
void vexp(double *a, double *b,  int n)  ;
void vclear(double *a,  double c, int n) ;
void vzero(double *a, int n) ;
void cpzero(char  **a, int n) ;
void ivvp(int *a, int *b, int *c, int n);
void ivvm(int *a, int *b, int *c, int n);
void ivsp(int *a, int *b, int c, int n);
void ivst(int *a, int *b, int c, int n);
void ivclear(int *a,  int c, long n) ;
void lvclear(long *a,  long c, long n) ;
void ivzero(int *a, int n) ;
void lvzero(long *a, long n) ;
void lvvp (long *a, long *b, long *c, int n) ;  
void lvvm (long *a, long *b, long *c, int n) ;  
void cclear(unsigned char *a,  unsigned char c, long n) ;
void charclear(char *a,  unsigned char c, long n) ;

double clip(double x, double lo, double hi)  ;
void ivclip(int *a, int *b,int loval, int hival,int n)  ; 
void vclip(double *a, double *b,double loval, double hival,int n)   ;

void vmaxmin(double *a, int n, double *max, double *min)  ;
void vlmaxmin(double *a, int n, int *max, int *min)  ;
void ivmaxmin(int *a, int n, int *max, int *min)  ;
int minivec(int *a, int n) ;
int maxivec(int *a, int n) ;
void ivlmaxmin(int *a, int n, int *max, int *min)  ;
void getdiag(double *a, double *b, int n)  ;
void setdiag(double *a, double *diag, int n)  ;
void adddiag(double *a, double *diag, int n)  ;
void flipiarr(int *a, int *b, int n) ;
void fliparr(double *a, double *b, int n)  ;
int ipow2 (int l) ;

void copyarr(double *a,double *b,int n) ;
void revarr(double *a, double *b,int n) ;
void reviarr(int *a,int *b,int n) ;
void revlarr(long *a,long *b,int n) ;
void revuiarr(unsigned int *a, unsigned int *b,int n) ;
void copyiarr(int *a,int *b,int n) ;
void copylarr(long *a, long *b, int n) ;
void copyiparr(int **a,int **b,int n) ;

void dpermute(double *a, int *ind, int len)  ;
void ipermute(int *a, int *ind, int len)  ;
void dppermute(double **a, int *ind, int len)  ;
void ippermute(int **a, int *ind, int len)  ;

double  asum(double *a, int n) ;
double  asum2(double *a, int n) ;
int     intsum(int *a, int n) ;
long  longsum(long *a, int n) ;
int     idot(int *a, int *b, int n)  ;
int     iprod(int *a, int n) ;
double  aprod(double *a, int n) ;
double  vdot(double *a, double *b, int n)  ;
double  corr(double *a, double *b, int n)  ;
double  corrx(double *a, double *b, int n)  ;
double  variance(double *a, int n)  ;
double trace(double *a, int n) ;
int  nnint(double a) ;
void countcat(int *tags, int n,int *ncat,int nclass)  ;
void rowsum(double *a, double *rr, int n) ;
void colsum(double *a, double *cc, int n) ;
void rrsum(double *a, double *cc, int m, int n)  ;
void ccsum(double *a, double *cc, int m, int n)  ;
void printmatfile(double *a, int m, int n, FILE *fff) ;
void printmatwfile(double *a, int m, int n, int w, FILE *fff) ;
void printmatx(double *a, int m, int n) ;
void printmat(double *a, int m, int n) ;
void printmatwx(double *a, int m, int n, int w) ;
void printmatw(double *a, int m, int n, int w) ;
void printmatl(double *a, int m, int n) ;
void printmatwl(double *a, int m, int n, int w) ;
void printmatwf(double *a, int m, int n, int w, char *format);
void int2c(char *cc, int *b, int n) ;
void floatit(double *a, int *b, int n) ;
void floatitl (double *a, long *b, int n) ;
void fixit(int  *a, double *b, int n) ;
void rndit(double  *a, double *b, int n) ;
void printimatw(int *a, int m, int n, int w) ;
void printimatx(int *a, int m, int n) ;
void printimat1(int *a, int m, int n) ;
void printimat1x(int *a, int m, int n) ;
void printimat(int *a, int m, int n) ;
void printimatl(int *a, int m, int n) ;
void printimatlfile(int *a, int m, int n, FILE *fff) ;
void printimatfile(int *a, int m, int n, FILE *fff) ;
void printimatwfile(int *a, int m, int n, int w, FILE *fff) ;
void printimat2D(int  **a, int m, int n)  ;
void printlmat (long *a, int m, int n) ;
void printmat2D(double **a, int m, int n)  ;
void printstring(char *ss, int width) ;
void printstringbasepos(char *ss, int w, int basepos) ;
void printstringf(char *ss, int width, FILE *fff) ;

int  findfirst(int *a, int n, int val) ;
int  findfirstl(long *a, int n, long val) ;
int  findfirstu(unsigned int *a, int n, unsigned int val) ;
int  findlastu(unsigned int *a, int n, unsigned int val) ;

int  findlast(int *a, int n, int val) ;
int  binsearch(int *a, int n, int val) ;
void idperm(int *a, int n)  ;
double NPlog2(double y) ;
double log2fac(int  n) ;
double logfac(int n)  ;
double logbino(int n, int k)  ;
double loghprob(int n, int a, int m, int k) ;  
int hprobv(double *vprob, int n, int a, int m) ;
/* hypergeometric probability */
double logmultinom(int *cc, int n) ;
double addlog(double a, double b) ;
double logsum(double *x, int n)  ;
double vldot(double *x, double *y, int n) ;
double pow10 (double x) ;
void vpow10 (double *a, double *b, int n) ;
void vlog10 (double *a, double *b, int n) ;
/* matrix transpose */
void transpose(double *aout, double *ain, int m, int n)  ;
void addoutmul(double *out, double *a, double mul, int n) ;
void addouter(double *out, double *a, int n) ;
void subouter(double *out, double *a, int n) ;
int mktriang(double *out, double *in, int n) ;
int mkfull(double *out, double *in, int n) ;

/* storage allocation */
int **initarray_2Dint(int numrows, int numcolumns, int initval);
long **initarray_2Dlong(int numrows, int numcolumns, long initval);
void free2Dint(int ***xx, int numrows) ;
void free2Dlong(long ***xx, int numrows) ;
double **initarray_2Ddouble(int numrows, int numcolumns, double initval);
long double **initarray_2Dlongdouble(int numrows, int numcolumns, long double initval);
void clear2D(double ***xx, int numrows, int numcols, double val)   ;
void iclear2D(int ***xx, int numrows, int numcols, int val)  ;
void lclear2D(long ***xx, int numrows, int numcols, long val)  ;
void free2D(double ***xx, int numrows) ;
void free2Dlongdouble(long double ***xx, int numrows) ;
void free_darray (double  **xx) ; 
void free_iarray (int  **xx)  ;          

double bal1 (double *a, int n)  ;
double bal2 (double *a, int n)  ;
void vcompl(double *a, double *b, int n) ;
void setidmat(double *a, int n) ; 

int stripit(double *a, double *b, int *x, int len) ;
int istripit(int *a, int *b, int *x, int len) ;
int cstripit(char **a, char **b, int *x, int len) ;

void mapit(int *a, int *b, int n, int inval, int outval) ;
int ifall(int n, int k)  ;  // falling factorial = n (n-1) (n-2) ... (n-k+1)
double hlife(double val) ;
void *topheap () ;

void swap (double *pa, double *pb) ;
void iswap (int *pa, int *pb) ;
void cswap(char *c1, char *c2) ;


void floatit2D(double **a, int **b, int nrows, int ncols)  ; 
void copyarr2D(double **a, double **b, int nrows, int ncols) ;  // a input b output
void copyiarr2D(int **a, int **b, int nrows, int ncols) ;  // a input b output
void plus2Dint(int **a, int **b, int **c, int nrows, int ncols) ;
void minus2Dint(int **a, int **b, int **c, int nrows, int ncols) ;

void plus2D(double **a, double **b, double **c, int nrows, int ncols) ;
void minus2D(double **a, double **b, double **c, int nrows, int ncols) ;
void sum2D(double *a, double **b, int nrows, int ncols) ;
double total2D(double **a, int nrows, int ncols) ;
int total2Dint(int **a, int nrows, int ncols) ;

int kodeitb(int *xx, int len, int base) ;
int dekodeitb(int *xx, int kode, int len, int base) ;
long lkodeitbb(int *xx, int len, int *baselist)  ;
int ldekodeitbb(int *xx, long kode, int len, int *baselist) ;
int kodeitbb(int *xx, int len, int *baselist)  ;
int dekodeitbb(int *xx, int kode, int len, int *baselist) ;

long expmod(long a, long b, long n) ;
int isprime(long num) ;
long nextprime(long num) ;

int irevcomp (int xx, int stringlen) ;
long lrevcomp (long xx, int stringlen) ;
void ismatch(int *a, int *b, int n, int val) ;
int pmult(double *a, double *b, double *c, int na, int nb) ;
void pdiff(double *a, double *b, int deg) ;
void vswap(double *a, double *b, int n)  ;
void setlong(long *pplen, long a, long b)   ;

long lmod (long x, long base)  ;
long gcdx(long b, long a, long *x, long *y) ; 
long modinv(long a, long base) ;  
long lpow2(int n) ; 
double exp1minus(double x) ;

double cputime (int mode) ;
double calcmem (int mode) ;

