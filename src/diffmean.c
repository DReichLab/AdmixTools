#include <stdio.h>
#include <stdlib.h>
#include <nicklib.h> 
#include <getpars.h> 

char *iname = NULL, *jname = NULL; 
int xflag = NO ;  
int terse = NO ;  

int getd (char *name, double *mm, double *vv, int n)  ;
void readcommands(int argc, char **argv) ;

int main(int argc, char **argv)
{

 int n1, n2, n, nx ; 
 double *m1, *var1, *m2, *var2, *dm, *dv, *dq, ychi, ytail ; 

 if (argc == 3) {
  iname = strdup(argv[1]) ; 
  jname = strdup(argv[2]) ; 
 }

 else { 
  readcommands(argc, argv) ; 
 }
   if (iname==NULL) printf ("Usage: diffmean -i in1 -j in2 [-x]\n")  ;
   if (jname==NULL) printf ("Usage: diffmean -i in1 -j in2 [-x]\n")  ;

   if ((iname == NULL) || (jname == NULL)) return -1 ; 
 n = n1 = numcols(iname) ; 
 n2 = numcols(jname) ; 
 if (n1 != n2) fatalx("mismatch source number\n") ;
 ZALLOC(m1, n,  double) ;
 ZALLOC(m2, n,  double) ;
 ZALLOC(var1, n*n,  double) ;
 ZALLOC(var2, n*n,  double) ;

 nx = getd(iname, m1, var1, n) ; 
 nx = getd(jname, m2, var2, n) ; 

if (terse == NO) {
 printmat(m1, 1, nx) ;
 printmat(m2, 1, nx) ;

 printnl() ;
 printmatl(var1, nx, nx) ;
 printmatl(var2, nx, nx) ;

}

 ZALLOC(dm, nx, double) ;
 ZALLOC(dv, nx*nx, double) ;
 ZALLOC(dq, nx*nx, double) ;
 
 vvm(dm, m1, m2, nx) ; 
 vvp(dv, var1, var2, nx*nx) ; 
 pdinv(dq, dv, nx) ; 

 ychi = scx(dq, NULL, dm, nx) ; 
 ytail = rtlchsq(nx, ychi) ;
 
 if (terse == NO) {
  printf("%15s ", iname) ; printmat(m1, 1, nx) ; 
  printf("%15s ", jname) ; printmat(m2, 1, nx) ; 
  printf("%15s ", "diff") ; printmat(dm, 1, nx) ; 
  printf("result: dof: %d chisq: %9.3f  tail: %12.6f\n", nx, ychi, ytail) ;
 }

 else { 
  printf(" %9.3f ", m1[0]) ;
  printf(" %9.3f ", m2[0]) ;
  printmatx(dm, 1, nx) ;  
  printf(" %12.6f", ytail) ;  
  printnl() ;
 }
  
 return 0; 

}

int getd (char *name, double *mm, double *vv, int n) 
{
 double **xx ; 
 int j, k, kk ;
 int nx ; 

 xx = initarray_2Ddouble(n+1, n+1, 0) ; 
 nx = n-1 ; 
 if (xflag) ++nx ;
 
 getxx(xx, n+1, n, name) ; 
 for (j=0; j<nx; ++j) { 
   mm[j] = xx[j][0] ; 
  for (k=1; k<=nx; ++k) { 
   kk = k-1 ; 
   vv[kk*nx+j] =  xx[j][k] ; 
  } 
 }
 vst(mm, mm, 1/1000.0, nx) ;  
 vst(vv, vv, 1/(1000.0*1000.0), nx*nx) ;  
 free2D(&xx, n+1) ;

 return nx ;

}
   
void readcommands(int argc, char **argv) 

{
  int i ;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "i:j:xt")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'j':
	jname = strdup(optarg) ;
	break;

      case 'x':
	xflag  = YES ;
	break; 

      case 't':
	terse  = YES ;
	break; 

       case '?': 
        printf ("Usage: diffmean -i in1 -j in2 [-x][-t]\n")  ;
	fatalx("bad params\n") ;
      }
  }


}

