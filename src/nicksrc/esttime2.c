#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <getpars.h> 
#include "regsubs.h"
#include "arms.h"

int verbose = NO ; 
#define  VERSION  "100"  

char *iname = "data2.txt" ; 
int *lomax, *himin, *lo, *hi ;
double *rr, *cc, *qq, *mu, *data ;  
double *rx, *rd, *cx, *cd ; 
double *ww ; 

double sp, slam, tp, tlam, qp, qlam ;  
int nspecies, ntimes, ncell ;
char **species ; 
char **pubage ; 
double *xleft, *xright, *xtot ; 
double *gpsi ; 

double b1, b2, yl1, yl2 ;

int bigiter = 100, numburn=100 ; 


struct xparm { 
 double lo, hi ; 
 double xa[3] ; 
 double ca[3] ; 
 double pnum; 
 double lval; 
 double blast[2] ; 
 double ysum; 
 double ylsum; 
 double yprior ; 
 int mode ; 
} ; 

double evgamma(double x, void *xdata)  ; 

void setlm(int *plo, int *phi, int lomax, int himin, double r, double *cc, double gp)   ; 
void mlegamm(double *a, double n, double *p, double *lam)  ;  
void readcommands(int argc, char **argv) ;
void rghyper(double *a, int n, double *psamp, double *plval, int mode)   ; 
void rbhyper(double *a, int n, double *p1, double *p2, int mode)   ; 
double loadrc(double *ww, double *rr, double *cc, int ntimes, int nspecies)  ;
void samprcq()  ;
void setqq( double *qq, int *lo, int *hi, int ntimes, int nspecies)  ;
void printit(char *str)   ;
void phist(char *sname, int *loh, int *hih, int num) ;  
double estgpsi(double *gp, int *hi, int *lo, int n)  ; 
void estbpars(double *a, int n, double *p1, double *p2)  ; 
double cputime (int mode)  ; 
void pcint(char *sname, double *ff, double *ll) ;
void calcxlr() ; 

int main(int argc, char **argv) 
{

  int nrows, ncols ; 
  int i, j, k, a ; 
  double **xx, y, ytop, ybot, ytotal, ya, yb ; 
  char **age,  ss[256] ; 
  int iter ; 
  int *lohist, *hihist ; 
  double *fhist, *lhist ; 
  int numiter = 0 ;
  double *specsum, *coltotal, *specave ; 
  double *ab1, *ab2 ;
  double  *ww0, *ww1, *ww2 ;
  double y1, y2 ;

  cputime(0) ; 
  readcommands(argc, argv) ;
  printf("## esttime2: iname: %s version: %s\n", iname, VERSION) ;

  nrows = numlines(iname) ;
  ncols = numcols(iname) ;  

//  printf("zznrows: %d  ncols: %d\n", nrows, ncols) ; 

  --ncols ;

  xx = initarray_2Ddouble(ncols, nrows, 0.0) ;
  ZALLOC(species, nrows, char *) ;
  nrows = getxxnames(&species, xx, nrows, ncols, iname) ;

  printf("nrows: %d  ncols: %d\n", nrows, ncols) ; 
  ZALLOC(specsum, nrows, double) ;
  ZALLOC(specave, nrows, double) ;
  ZALLOC(coltotal, ncols, double) ;
  ZALLOC(age, ncols, char *) ;
  ZALLOC(pubage, ncols, char *) ;

  nspecies = nrows ; 
  ntimes = ncols ; 

  sum2D(specsum, xx, ntimes, nspecies) ; 

  ZALLOC(gpsi, nspecies, double) ; b1 = 9; b2 = 1; 
  for (j=0; j<nspecies; ++j) { 
   gpsi[j] = ranbeta(b1, b2) ;
  } 

  for (k=0; k<ntimes; ++k) { 
   y = 142 - (double) k ; 
   ytop = y/10.0 ; 
   ybot = ytop - 0.1 ; 
   sprintf(ss, "%6.1f - %5.1f", ybot, ytop) ; 
   age[k] = strdup(ss) ; 
   substring(&age[k], " ", "") ; 
   sprintf(ss, "%6.1f", ybot) ; 
   if ((k==0) || (k==(ntimes-1))) sprintf(ss, "%6.1f*", ybot) ; 
   pubage[k] = strdup(ss) ; 
   substring(&pubage[k], " ", "") ; 
   coltotal[k] = asum(xx[k], nspecies) ;
  }
// 
  for (k=ntimes-1; k>=0; --k) { 
   if (coltotal[k] > 0) break ;  
  }

  ntimes = k+1 ; 
/**
  for (k=0; k<ntimes; ++k) { 
   printf("time: %12s %6.0f\n", age[k], coltotal[k]) ; 
  }
*/  


  vclip(specsum, specsum, 0.1, 1.0e6, nspecies) ; 
  vclip(coltotal, coltotal, 0.1, 1.0e6, ntimes) ; 
  

  ZALLOC(lomax, nspecies, int) ;
  ZALLOC(himin, nspecies, int) ;
  ivclear(lomax, 9999, nspecies) ;
  ivclear(himin, -9999, nspecies) ;
  ZALLOC(lo, nspecies, int) ;
  ZALLOC(hi, nspecies, int) ;

  ZALLOC(ab1, nspecies, double) ;
  ZALLOC(ab2, nspecies, double) ;
  ZALLOC(ww0, nspecies, double) ;
  ZALLOC(ww1, nspecies, double) ;
  ZALLOC(ww2, nspecies, double) ;
  ZALLOC(xleft, nspecies, double) ;
  ZALLOC(xright, nspecies, double) ;
  ZALLOC(xtot, nspecies, double) ;

  for (j=0; j<nspecies; ++j) { 
   for (k=0; k<ntimes; ++k) { 
     a = nnint(xx[k][j]) ; 
     if (a==0) continue ; 
     lomax[j] = MIN(lomax[j], k) ; 
     himin[j] = MAX(himin[j], k) ; 
   }
   y1 = himin[j] - lomax[j] + 1 ; 
   specave[k] = specsum[k] / y1 ; 
  }
  a = nnint(intsum(himin, ntimes)) ; 
  a -= nnint(intsum(lomax, ntimes)) ; 
  y = (double) a / (double) nspecies ; 
  printf("mean length of observed species occurrence: %9.3f\n", y+1) ; 
 

  for (k=0; k<nspecies; ++k) { 
   printf("species: %30s %6.0f", species[k], specsum[k]) ; 
   printf("%3d %3d", lomax[k], himin[k]) ;
   printnl() ; 
  }
  ncell = nspecies * ntimes ;   
  ZALLOC(rr, ntimes, double) ; 
  ZALLOC(rx, ntimes, double) ; 
  ZALLOC(rd, ntimes, double) ; 
  ZALLOC(cc, nspecies, double) ; 
  ZALLOC(cx, nspecies, double) ; 
  ZALLOC(cd, nspecies, double) ; 
  ZALLOC(qq, nspecies*ntimes, double) ; 
  ZALLOC(mu, nspecies*ntimes, double) ; 
  ZALLOC(data, nspecies*ntimes, double) ; 
  ZALLOC(ww, ncell, double) ; 
  ytotal = asum(specsum, nspecies) ; 
  vst(cc, specsum, 1/sqrt(ytotal), nspecies) ;
  vst(rr, coltotal, 1/sqrt(ytotal), ntimes) ;
  ZALLOC(lohist, ncell, int) ; 
  ZALLOC(hihist, ncell, int) ; 
  ZALLOC(fhist, ncell, double) ; 
  ZALLOC(lhist, ncell, double) ; 
// initialize 

  for (i=0; i<ntimes; ++i) {
   for (j=0; j<nspecies; ++j) {
    k = j*ntimes + i ; 
    qq[k] = rangam(2.0)/2.0 ; 
    mu[k] = rr[i]*cc[j]*qq[k] ;
    data[k] =  xx[i][j] ; 
  }} 
  qp = 2 ; qlam = 2 ; 
  for (j=0; j<nspecies; ++j) { 
   setlm(&lo[j], &hi[j], lomax[j], himin[j], cc[j], rr, gpsi[j]) ; 
   printf("%40s ", species[j]) ; 
   printf("%4d %4d   %4d %4d", lomax[j], himin[j], lo[j], hi[j]) ;
   printnl() ; 
  }
  
  mlegamm(rr, ntimes, &tp, &tlam) ; 
  mlegamm(cc, nspecies, &sp, &slam) ; 

   for (j=0; j<nspecies; ++j) { 
    setlm(&lo[j], &hi[j], lomax[j], himin[j], cc[j], rr, gpsi[j]) ; 
   }

   floatit(ab2, lomax, nspecies) ; // work area 
   floatit(ab1, himin, nspecies) ;
   vvm(ab1, ab1, ab2, nspecies) ;  
   vsp(ab1, ab1, 1, nspecies) ;
   vvd(ab2, specsum , ab1, nspecies) ; 
   copyarr(specsum, ab1, nspecies) ;
   copyarr(ab1, ww0, nspecies) ;

   setqq( qq, lo, hi, ntimes, nspecies) ;

   rghyper(qq,   ncell , &qp, &qlam, 12) ;
   estgpsi(gpsi, hi, lo, nspecies) ;
   rbhyper(gpsi,   nspecies , &b1, &b2, 12) ;
   printf("init hyper parameters(mle): %9.3f %9.3f   %9.3f %9.3f   %9.3f %9.3f", sp, slam, tp, tlam, qp, qlam) ;
   printf("   %9.3f %9.3f", b1, b2) ;
   printnl() ;


  rghyper(rr,   ntimes , &tp, &tlam, 10) ;
  rghyper(cc,   nspecies , &sp, &slam, 11) ;
  rghyper(qq,   ncell , &qp, &qlam, 12) ;

  for (iter= -numburn; iter <= bigiter; ++iter) { 
   rghyper(rr,   ntimes , &tp, &tlam, 0) ;
   rghyper(cc,   nspecies , &sp, &slam, 1) ;
   rbhyper(gpsi,   nspecies , &b1, &b2, 0) ;
   
   for (j=0; j<nspecies; ++j) { 
    setlm(&lo[j], &hi[j], lomax[j], himin[j], cc[j], rr, gpsi[j]) ; 

    if (iter >= 1) { 
      for (i=lo[j]; i<=lomax[j]; ++i)  { 
       k = j*ntimes + i ; 
       ++lohist[k] ; 
      }
      for (i=hi[j]; i>=himin[j]; --i)  { 
       k = j*ntimes + i ; 
       ++hihist[k] ; 
      }
      i = lo[j] ; k = j*ntimes+i; ++fhist[k] ;  
      i = hi[j] ; k = j*ntimes+i; ++lhist[k] ;  
    }
   }
   if (iter >= 1) ++numiter ; 
   setqq( qq, lo, hi, ntimes, nspecies) ;
   rghyper(qq,   ncell , &qp, &qlam, 2) ;
   estgpsi(gpsi, hi, lo, nspecies) ;
   rbhyper(gpsi,   nspecies , &b1, &b2, 0) ;
   loadrc(qq, rr, cc, ntimes, nspecies) ;
   printf("iter: %4d %9.3f %9.3f   %9.3f %9.3f   %9.3f %9.3f", iter, sp, slam, tp, tlam, qp, qlam) ;
   printf("   %9.3f %9.3f", b1, b2) ;
   printnl() ;


/**
    floatit(ww1, lo, nspecies) ;
    floatit(ww2, hi, nspecies) ;

    vvm(ww2, ww2, ww1, nspecies) ; 
    vvm(ww1, ww2, ww0, nspecies) ; 

    y1 = corr(ab1, ww1, nspecies) ;
    y2 = corr(ab2, ww1, nspecies) ;
*/

//  printf("corr: %6d %9.3f %9.3f", iter, y1, y2) ;
    printf("corr: %6d", iter) ;

    calcxlr()  ; 
    y1 = corr(xleft, specsum, nspecies) ;
    y2 = corr(xright, specsum, nspecies) ;
    printf("     %9.3f %9.3f",  y1, y2) ;

    y1 = corr(xleft, specave, nspecies) ;
    y2 = corr(xright, specave, nspecies) ;
    printf("     %9.3f %9.3f",  y1, y2) ;
    printnl() ;
    


  }


  for (j=0; j<nspecies; ++j) { 
   phist(species[j], lohist+j*ntimes, hihist+j*ntimes, numiter) ; 
  }

  for (j=0; j<nspecies; ++j) { 
   pcint(species[j], fhist+j*ntimes, lhist+j*ntimes) ; 
  }
 
  printf("end of esttime.  CPU :: %9.3f seconds\n", cputime(1)) ; 
  return 0 ; 

}
void samprcq() 
{
 int i, j, k ; 

 for (i=0; i<ntimes; i++) { 
  for (j=0; j<nspecies; j++) { 

     if (i<lo[j]) continue ; 
     if (i>hi[j]) continue ; 
     k = j*ntimes + i ; 
     if (ww[k] < -0.5) return ; 
     qq[k] = rangam(qp+data[k]) / (qlam+rr[i]*cc[j]) ; 


 }}

 for (i=0; i<ntimes; i++) { 
   rr[i] = rangam(tp+rd[i])/(tlam+rx[i])  ;
 }

 for (j=0; j<nspecies; j++) { 
   cc[j] = rangam(sp+cd[j])/(slam+cx[j]) ; 
 }
 
 


}
double loadrc(double *ww, double *rr, double *cc, int ntimes, int nspecies) 
// ww not altered
{
 int i, j, k ; 
 double w ; 

 vzero(rx, ntimes) ; 
 vzero(rd, ntimes) ; 
 vzero(cx, nspecies) ; 
 vzero(cd, nspecies) ; 


/**
  mlegamm(rr, ntimes, &tp, &tlam) ; 
  mlegamm(cc, nspecies, &sp, &slam) ; 
*/
 
 for (i=0; i<ntimes; i++) { 
  for (j=0; j<nspecies; j++) { 
     k = j*ntimes + i ; 
     if (i<lo[j]) continue ; 
     if (i>hi[j]) continue ; 
     w = ww[k] ; 
     if (w < -0.5) continue ; 
     rd[i] += data[k] ; 
     cd[j] += data[k] ; 
     rx[i] += w*cc[j] ; 
     cx[j] += w*rr[i] ; 
  }} 

 return asum(rx, ntimes) ; 

}
void estbpars(double *a, int n, double *p1, double *p2) 

{

// moments estimator for beta 

  double *w1, ymean, yvar ; 

  ZALLOC(w1, n, double) ;

  ymean = asum(a, n) / (double) n ; 
  vsp(w1, a, -ymean, n) ; 
  yvar = asum2(w1, n) / (double) n ; 

  bpars(p1, p2, ymean, yvar) ; 

  free(w1) ;


} 

void setlm(int *plo, int *phi, int lomax, int himin, double a, double *bb, double gp)  
{
  int k, tlo, thi,  t, x ; 
  double y,  p ; 
  double ww[1000], wend[1000] ;  
  static long ncall = 0 ;

  ++ncall ; 
  vclear(wend, -1, 1000) ; 
  tlo = lomax ; 
// tlo is probability we are alive at 
  t = 0 ; ww[t] = wend[t] = 1 ; ++t ; 
  for (k=lomax-1; k>=0; --k) { 
// 
    y = qp*(log(qlam) - log(qlam+a*bb[k])) ; 
    y = exp(y) ;   // probability of 0  | alive
    ww[t] = ww[t-1] * y * gp  ;    
    wend[t-1] = ww[t-1] - ww[t] ; 
    wend[t] = ww[t] ;  
    ++t ; 
  }
  if (ncall <= -1) {  
   printf("zzww ") ; printmat(wend, 1, t) ; 
  }

  bal1(wend, t) ;  x = randis(wend,t) ;  tlo = lomax - x ; 

  thi = himin ; 
  t = 0 ; ww[t] = wend[t] = 1 ; ++t ; 
  for (k=himin+1; k<ntimes; ++k) { 
    y = qp*(log(qlam) - log(qlam+a*bb[k])) ; 
    y = exp(y) ;   // probability of 0 
    ww[t] = ww[t-1] * y * gp ;   
    wend[t-1] = ww[t-1] - ww[t] ; 
    wend[t] = ww[t] ; 
    ++t ; 
  }

  bal1(wend, t) ;  x = randis(wend,t) ;  thi = himin + x  ; 

  *plo  = tlo ; 
  *phi  = thi ; 


}



void readcommands(int argc, char **argv) 

{
  int i;
  phandle *ph ;
  char str[512]  ;
  int n, kode ;
  int pops[2] ;

  while ((i = getopt (argc, argv, "i:n:b:")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'b':
	numburn = atoi(optarg) ;
	break;

      case 'n':
	bigiter = atoi(optarg) ;
	break;

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
}

void mlegamm(double *a, double n, double *p, double *lam) 
{

  double a1, a2 ; 
  int k, num=0 ; 

  a1 = a2 = 0 ;
  for (k=0; k<n; ++k) { 
   if (a[k] <= 0.0) continue ; 
   a1 += a[k] ; 
   a2 += log(a[k]) ; 
   ++num ; 
  }
  a1 /= (double) num ; 
  a2 /= (double) num ; 
  mleg(a1, a2, p, lam) ;

}

double evgamma(double x, void *xdata) 
{
  struct xparm *d;
  double y, z, xl, w ;
  int k ; 
  double ysum, ylsum, pnum, lval ; 

  d = xdata ; 
  if (x <= d -> lo) return -1.0e6 ;
  if (x >= d -> hi) return -1.0e6 ;

  pnum = d -> pnum ; 
  lval =  d -> lval ; 
//  if (d -> mode == 2)  lval = x ; 
  ysum =  d -> ysum ; 
  ylsum = d -> ylsum ;
  y = pnum*x*log(lval) ; 
  y -= pnum*lgamma(x)  ; 
  y += (x-1)*ylsum ; 
  y -= lval*ysum ;   
  if (!finite(y)) fatalx("evgamma bug\n") ;  
  return y ;
}


double evbeta(double x, void *xdata) 
// yl1, yl2 must be set from psi 
{
  struct xparm *d;
  double y, z, xl, w, pnum ;
  double a1, a2 ;

  d = xdata ; 
  if (x <= d -> lo) return -1.0e6 ;
  if (x >= d -> hi) return -1.0e6 ;
  pnum = d -> pnum ; 

  a1 = b1 ; a2 = b2 ; 

  if (d->mode == 1) a1 = x ; 
  if (d->mode == 2) a2 = x ; 

  y = -pnum*lbeta(a1, a2) ;
  
  y += (a1-1)*yl1 ; 
  y += (a2-1)*yl2 ; 

  return y ;

}
void rbhyper(double *a, int n, double *p1, double *p2, int mode)  
{

 int num = 0 ; 
 static double xxprev[3], xxlval[3], prev, lval, pval ; 
 static double xxinit[3][8] ; 
  double xinit [8] ;
  int lo, hi, i, k, t ;
  double y, yy, yinc, ylo, yhi, ylike ; 
  double yclip = .0001 ;

  double xsamp[10], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  double xprev, y0, y1, y2,  a1, a2, b1, b2, pnum ;
  int dometrop = 1;
  int err, ninit = 3, npoint = 101, nsamp = 1, ncent = 0 ;
  int neval, nhist;
  static struct xparm xdata, *dpt, *d ;
  static double xz[2][3] ;


 dpt = d = &xdata ; 

  if (mode >= 10) { 
   estbpars(a, n, p1, p2) ;


   b1 = d -> blast[0] = *p1 ;
   b2 = d -> blast[1] = *p2 ;

   xz[0][0] = b1*0.9 ; 
   xz[0][1] = b1 ; 
   xz[0][2] = b1*1.1 ; 

   xz[1][0] = b2*0.9 ; 
   xz[1][1] = b2 ; 
   xz[1][2] = b2*1.1 ; 
  }

 yl1 = yl2 = 0 ;
 d -> lo = 0 ; 
 d -> hi = 1.0e6 ;  
  for (k=0; k<n; ++k) { 
   y = clip(a[k], yclip, 1-yclip) ;
   yl1 += log(y) ;
   yl2 += log(1-y) ;
  }

  d -> mode = 1 ; 

  copyarr(xz[0], xinit, 3) ;

  d -> pnum = n ; 
  prev = d -> blast[0] ; 

  ylo = dpt -> lo = 0 ; 
  yhi = dpt -> hi = 1.0e6 ; 

    ninit = 3 ;
    err = arms(xinit,ninit,&ylo,&yhi, evbeta, &xdata, &convex,
           npoint,dometrop,&prev,xsamp,nsamp,
           qcent,xcent,ncent,&neval);

    *p1 = d -> blast[0] = xsamp[0] ;

    
  d -> mode = 2 ; 

  copyarr(xz[1], xinit, 3) ;

  d -> pnum = n ; 
  prev = d -> blast[1] ; 

  ylo = dpt -> lo = 0 ; 
  yhi = dpt -> hi = 1.0e6 ; 

    ninit = 3 ;
    err = arms(xinit,ninit,&ylo,&yhi, evbeta, &xdata, &convex,
           npoint,dometrop,&prev,xsamp,nsamp,
           qcent,xcent,ncent,&neval);

    *p2 = d -> blast[1] = xsamp[0] ;

}

void rghyper(double *a, int n, double *psamp, double *plval, int mode)  
{

 int num = 0 ; 
 static double xxprev[3], xxlval[3], prev, lval, pval ; 
 static double xxinit[3][8] ; 
  double xinit [8] ;
  int lo, hi, i, k, t ;
  double y, yy, yinc, ylo, yhi, ylike ; 

  double xsamp[10], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  double xprev, y0, y1, y2,  a1, a2, b1, b2, pnum ;
  int dometrop = 1;
  int err, ninit = 3, npoint = 101, nsamp = 1, ncent = 0 ;
  int neval, nhist;
  struct xparm xdata, *dpt, *d ;

 prev = *psamp ; 
 lval = *plval ;

 dpt = d = &xdata ; 
 d -> lo = 0 ; 
 d -> hi = 1.0e6 ;  
  a1 = a2 = 0 ;
  for (k=0; k<n; ++k) { 
   if (a[k] <= 0.0) continue ; 
   a1 += a[k] ; 
   a2 += log(a[k]) ; 
   ++num ; 
  }
  t = mode %10 ; 
  pnum = num ; 

  d -> mode = t ; 

  if (mode >= 10) { 
   b1 = a1/pnum ; 
   b2 = a2/pnum ; 
   mleg(b1, b2, &pval, &lval) ; 
   prev = pval ; 
   xxprev[t] = pval ;
   xxlval[t] = lval ; 
   y0 = MAX(pval-0.1, pval/2) ;
   y1 = pval ;     
   y2 = pval + 1.0 ;  
   xxinit[t][0] = y0; xxinit[t][1] = y1; xxinit[t][2] = y2 ;
  }

  d -> pnum = num ; 
  d -> lval = lval ; 
  
  d -> ysum = a1 ; 
  d -> ylsum = a2 ; 


  copyarr(xxinit[t], xinit, 3) ; 

  prev = xxprev[t] ; 

  ylo = dpt -> lo = 0 ; 
  yhi = dpt -> hi = 1.0e6 ; 

    ninit = 3 ;
    err = arms(xinit,ninit,&ylo,&yhi, evgamma, &xdata, &convex,
           npoint,dometrop,&prev,xsamp,nsamp,
           qcent,xcent,ncent,&neval);

    pval = *psamp  = xsamp[0] ;
    *plval = rangam(d->pnum*pval+1)/d->ysum ; 

    xxprev[t] = pval ; 
    xxlval[t] = *plval ; 
    

}

void setqq( double *qq, int *lo, int *hi, int ntimes, int nspecies) 

{
   int i, j, k ;

   vclear(qq, -1, ntimes*nspecies) ; 

   for (j=0; j< nspecies; ++j) { 
    for (i=lo[j]; i<=hi[j]; i++) { 
      k = j*ntimes + i ; 
      qq[k] = rangam(qp + data[k]) / (qlam+rr[i]*cc[j]) ; 
      qq[k] = clip(qq[k], .001, 1.e8) ;
    } 
   }
    
}

void printit(char *str)  
{
  int i,j,k ; 
  double y1 ;

  printf("%s  hyper parameters: %12.6f %12.6f   %12.6f %12.6f   %12.6f %12.6f\n", str, sp, slam, tp, tlam, qp, qlam) ;
  for (j=0; j<nspecies; ++j) { 
   printf("%s %40s ", str, species[j]) ; 
   printf("%4d %4d   %4d %4d", lomax[j], himin[j], lo[j], hi[j]) ;
   printnl() ; 
   for (i=lo[j]; i<=hi[j]; ++i) { 
    k = j*ntimes+i ; 
    y1 = rr[i]*cc[j] ;
    printf(" %9.3f %9.3f ", y1, qq[k]) ; 
    printf(" %9.3f %9.3f ", qq[k]*y1, data[k]) ; 
    printnl() ; 
   }
  }
}
void phist(char *sname, int *loh, int *hih, int num) 
{
 char ss1[4096], ss2[4096], *s1, *s2 ; 
 int k ; 
 double y ; 


 s1 = ss1 ; s2 = ss2 ; 
 s1 += sprintf(s1, "%30s ", sname) ; 
 s2 += sprintf(s2, "%30s ", "") ; 
 for (k=0; k<ntimes; ++k) { 
   if (loh[k] == 0) continue ; 
   y = (double) loh[k] / (double) num ;  
   s1 += sprintf(s1, "%6s ", pubage[k]) ; 
   s2 += sprintf(s2, "%6.3f ", y) ; 
 }
 s1 += sprintf(s1, " ... ") ;    
 s2 += sprintf(s2, " ... ") ;    
 for (k=0; k<ntimes; ++k) { 
   if (hih[k] == 0) continue ; 
   y = (double) hih[k] / (double) num ;  
   s1 += sprintf(s1, "%6s ", pubage[k]) ; 
   s2 += sprintf(s2, "%6.3f ", y) ; 
 }
 printf ("%s\n", ss1) ;
 printf ("%s\n", ss2) ;


}
void pcint(char *sname, double *ff, double *ll) 
{
  double *ww, y, thresh = .95 ; 
  int k, tlo, thi, x ; 
  ZALLOC(ww, ntimes, double) ; 

  x = indxstring(species, nspecies, sname) ; 
  copyarr(ff, ww, ntimes) ; bal1(ww, ntimes) ; 
  y = 0;  
  thi = -999 ; 
  for (k=ntimes; k>=0; --k) {
   if (ww[k] < .0001) continue ; 
   thi = MAX(k, thi) ; 
   y += ww[k] ; if (y>thresh) break ; 
  }
  
  tlo = k ; 
  printf("credible: %30s ", sname) ; 
  printf(" %6s %6s    ", pubage[tlo], pubage[thi]) ;
  
  copyarr(ll, ww, ntimes) ; bal1(ww, ntimes) ; 

  y = 0;  
  tlo = 99999 ; 

  for (k=0; k<ntimes; ++k) {
   if (ww[k] < .0001) continue ; 
   tlo = MIN(k, tlo) ; 
   y += ww[k] ; if (y>thresh) break ; 
  }
  
  thi = k ; 
  printf(" %6s %6s", pubage[tlo], pubage[thi]) ;
  printnl() ; 
  
  
 free(ww) ;

}
void calcxlr() 

{ 

 int k ; 
 double y1, qpsi ; 

 for (k=0; k<nspecies; ++k) {  
  qpsi = gpsi[k] ;
  y1 = lomax[k] - lo[k] ; 
  if (lo[k] ==0 )  { 
    y1 += (rangeom(1-qpsi) - 1) ;
  }
  xleft[k] = y1 ;

  y1 = hi[k] -  himin[k] ; 
  if (hi[k] == (nspecies-1))  { 
    y1 += (rangeom(1-qpsi) - 1) ;
  }
  xright[k] = y1 ;
 }

  vvp(xtot, xleft, xright, nspecies) ;


}
double estgpsi(double *gp, int *hi, int *lo, int n) 
{

  double ya, yb, qpsi ; 
  int k, nbonus=0, gg ; 
  double yclip = .0001, y ;

  yl1 = yl2 = 0 ; 
  for (k=0; k<n; ++k) { 
   qpsi = gp[k] ;
   ya = hi[k]-lo[k] ;  
   yb = 1 ;
   if ((lo[k]==0) || (hi[k]==(n-1))) {
    ya += (rangeom(1-qpsi) - 1) ;
   }
   y = gp[k] = clip(ranbeta(ya+b1, yb+b2), yclip, 1-yclip) ;
   yl1 += log(y) ;
   yl2 += log(1-y) ;
  }
}
double cputime (int mode) 
{
  static double ttt=0 ; 

 if (mode==0) {  
  ttt = clocktime() ;
  return 0 ;
 }

 return clocktime() - ttt ; 
   




}
