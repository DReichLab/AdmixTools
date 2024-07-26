#include <fcntl.h>
#include <ctype.h>
#include  <nicklib.h>
#include  <getpars.h>


int gslsetup (int n, double *vpars, double initscale, double *xtri, double *xvar)  ; 
double gslopt(double *pars) ;
int initpars(double *pars, double *z, double *zv, int dim) ; 

void full2chol(double *mat, double *chol, int n) ;
void chol2full(double *mat, double *chol, int n) ; 
void chol2prod (double *prod, double *chol, int n) ; 
double rank1opt(double *ff, double *z, double *zv, int dim) ;
void solveco(double *coeff, double *mat, int dim)  ;

char *iname = "zdump" ; 
double pscale = 0.01 ;

int gsldetails = YES;
double gslprecision = .001;

void xcheck(int dim) ; 
void readcommands (int argc, char **argv) ;


int main(int argc, char **argv)   
{
  
 int dim, dv, t, s  ; 
 double **xx ; 
 double *z, *zv ; 
 int k, a, b ; 
 double *vpars, *ff,  *ans ;
 double y, ytail ; 
 double *w1, *w2 ; 

 readcommands(argc, argv) ;
 
 printf("## qptestr1 -i %s\n", iname) ; 

 dv = s = numcols(iname) ;
 t = s + 2 ; 
 dim = (int) sqrt(2*s) ; 
 printf("zz %d %d\n", dim, dv) ;

 xx = initarray_2Ddouble(s, dv+1, 0) ; 
 t = getxx(xx, dv+1, dv, iname) ; 
 ZALLOC(z, dv, double) ;
 ZALLOC(zv, dv*dv, double) ;
 ZALLOC(ans, dim, double) ;
 ZALLOC(w1, dim*dim, double) ;
 ZALLOC(w2, dim*dim, double) ;


 for (k=0; k<dv; ++k) { 
   z[k] = xx[k][0] ; 
   for (a=0; a<dv; ++a) { 
     zv[a*dv+k] = xx[k][a+1] ;
   }
 }

 mkfull(w1, z, dim) ;  

 printf("zzinput\n") ;
 printmat(w1, dim, dim) ;

 ZALLOC(ff, dim*dim, double) ;

 y = rank1opt(ff, z, zv, dim) ;           
 ytail = rtlchsq(1, y) ; 
 printf("chisq: %12.6f  p-value: %9.3f\n", y, ytail) ;
 solveco(ans, ff,  dim)  ;

 

 printf("fullans: ") ; printmat(ans, 1, dim) ;
 printnl() ;
 printmat(ff, dim, dim) ;
 return 0 ;

}


double rank1opt(double *ff, double *z, double *zv, int dim) 
{

 int dv, t, s  ; 
 double **xx ; 
 int k, a, b ; 
 double *vpars,  *ww, *ww2 ;
 double y ; 

 dv = dim*(dim+1)/2 ; 
 printmat(z, 1, dv) ; 
 printnl() ;
 printmat(zv, dv, dv) ; 

 ZALLOC(vpars, dv, double) ;

 initpars(vpars, z, zv, dim) ;
 printnl() ; 

 vclear(vpars, 1.0, dv) ;
 gslsetup(dim, vpars, pscale, z, zv) ;    

 y = gslopt(vpars) ; 
 chol2prod(ff, vpars, dim) ; 
 printf("ans: %9.3f\n", y) ; 
 printmat(ff, dim, dim) ; 

// cholesky check
 ZALLOC(ww, dim*dim, double) ;
 ZALLOC(ww2, dim*dim, double) ;

  copyarr(ff, ww2, dim*dim) ; 

  choldc(ww2, dim, ww) ; 
  printmatl(ww, 1, dim) ;



 return y ; 

}

int initpars(double *pars, double *z, double *zv, int dim) 
{
  double *ff, y ;  
  int x, a, b, dv ; 

  dv = dim*(dim+1)/2 ;
  vzero(pars, dv) ; 
  x = 0 ;
  for (a=0; a<dim; ++a) { 
   for (b=a; b<dim; ++b) { 
     pars[x] = 0 ; 
     if (a==b) { 
      y = fabs(z[x]) ; 
      pars[x] = sqrt(y) ; 
     }
     ++x ; 
   } 
  }
  printf("zzpars ") ; printmat(pars, 1, dv) ;

}

void full2chol(double *mat, double *chol, int n) 
{
  int a, b, x ;

  x = 0 ; 
  for (a=0; a<n; ++a) { 
   for (b=a; b<n; ++b) { 
    chol[x] = mat[a*n+b] ; 
    ++x ;  
   }
  }
}
void chol2full(double *mat, double *chol, int n) 
{
  int a, b, x ;
  
  vzero(mat, n*n) ;
  x = 0 ; 
  for (a=0; a<n; ++a) { 
   for (b=a; b<n; ++b) { 
    mat[a*n+b] = chol[x] ; 
    ++x ;  
   }
  }
}

void chol2prod (double *prod, double *chol, int n) 

{

 double *mat ; 
 ZALLOC(mat, n*n, double) ;
 chol2full(mat, chol, n) ;
 txmulx(prod, mat, n, n) ;

 free(mat) ;

}
double mkcanon(double *chol2, double *chol, int n) 
{

 double *full, *vv ; 
 double penalty = 0 ;
 int k ;
 double y ; 

 ZALLOC(full, n*n, double) ;
 chol2full(full, chol, n) ;
 
 for (k=0; k<n; ++k) { 
  y = full[k*n+k] ;  
  if (y<0) { 
   vv = full + k*n ; 
   vst(vv, vv, -1, n) ;
   penalty  -=  y ; 
  }
 }
 full2chol(full, chol2, n) ; 


 free(full) ; 
 return penalty ; 

}

double scorit(double *chol, double *xtri, double *xvar, int n) 
{

  int dv = n*(n+1)/2 ; 
  double *full, *tri ; 
  double y ; 

  ZALLOC(full, n*n, double) ; 
  ZALLOC(tri, dv, double) ;

  chol2prod(full, chol, n) ; 
  mktriang(tri, full, n) ; 
  y = scx(xvar, xtri, tri, dv) ; 

  free(full) ; 
  free(tri) ; 

}

void solveco(double *coeff, double *mat, int dim) 
// corank1 a matrix ; we do some scaling here
{

 double *wco, *ww, *rhs, *r2 ; 

 ZALLOC(ww, dim*(dim+1), double) ; 
 ZALLOC(wco, dim*dim, double) ; 
 ZALLOC(r2, dim, double) ; 
 ZALLOC(rhs, dim+1, double) ; 

 copyarr(mat, ww, dim*dim) ; 
 vclear(ww+dim*dim, 1.0, dim) ; 
 rhs[dim] = 1.0 ; 
 txmulx(wco, ww, dim+1, dim) ; 
 mulmat(r2, rhs, ww, 1, dim+1, dim) ; 

 solvit(wco, r2, dim, coeff) ;
 printf("zzz\n") ;
 printmat(mat, dim, dim) ; 
 printnl() ;
 printmat(ww, dim+1, dim) ; 
 printnl() ;
 printmat(rhs, 1, dim+1) ; 
 printnl() ;
 printmat(wco, dim, dim) ; 
 printnl() ;
 printmat(r2, 1, dim) ; 

}


void xcheck(int dim) 
{

  int dv ;
  double *ff, *ch ;
  
  ZALLOC(ff, dim*dim, double) ; 
  dv = dim*(dim+1)/2 ;
  dv = dim*(dim+1)/2 ;
  ZALLOC(ch, dv, double) ; 

  setidmat(ff, dim) ; 
  full2chol(ff, ch, dim) ; 
  printf("check") ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  printnl() ;
  printmatw(ch, 1, dv, dv) ; 
  chol2full(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  chol2prod(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  vst(ch, ch, 2, dv) ; 
  chol2full(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 
  chol2prod(ff, ch, dim) ; 
  printnl() ;  printmat(ff, dim, dim) ; 

 free(ff) ; 
 free(ch) ;
}

void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n, t;

  while ((i = getopt (argc, argv, "i:")) != -1) {

    switch (i) {

    case 'i':
      iname = strdup (optarg);
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }

}
