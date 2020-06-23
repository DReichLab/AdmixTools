#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <getpars.h> 

static int verbose = NO ; 
#define  VERSION  "121"  

// -m flag now works  (first col 0-n, or 1-n + -m flag

#define N  2 
#define NX 3

char *iname = "outf3" ;
double xmean = -1.0e20 ;
int vflag = NO ; 

void readcommands(int argc, char **argv) ;

int main(int argc, char **argv)
{
  
#define MAXSTR 512  
#define MAXFF  10 
#define MAXNUM  60000 
#define MAXBLOCKS  10000 
#define BSIZE  5000000 

  char line[MAXSTR] ;
  int lastc, lastb ;
  char *sx ;
  char *hpos ;
  char *spt[MAXFF] ;
  int slo = 0, shi = 19, klo = 20, khi = 29, clo = 30, chi = 49 ;  
  int nsplit , n, t,  cc ;
  FILE *fff ;

  int chrom, pos, bnum  ;

  double wjack[MAXBLOCKS] ;
  double jmean[MAXBLOCKS] ;
  double jwt[MAXBLOCKS] ;

  int numblocks = 0, block ;
  int num, k, nline, x, xmax ; 
  double wt, val, twt ;
  double y1, y2, est, sig, mean, tnum, tden ;
  double **xx, *xvals, *xwt ;
  int *xind ; 

 readcommands(argc, argv) ;
/** 
 iname has 3 columns 
 block # weight jackknife value.   
 First row (block 0) is global mean 
*/

 openit(iname, &fff, "r") ;

 nline = numlines(iname)  ;  
 xx = initarray_2Ddouble(3, nline, 0) ;         
 bnum = nline = getxx(xx, nline, 3, iname) ;

 ZALLOC(xwt, nline+1, double) ; 
 ZALLOC(xvals, nline+1, double) ; 
 ZALLOC(xind,  nline+1, int) ; 

 xind[0] = -1 ;

 if (xmean > -1.0e15) { 
  xwt[0] = 0 ; 
  mean = xvals[0] = xmean ; 
  xind[0] = 0 ; 
 }  
 for (k=0; k<nline; ++k)  { 
  x = nnint(xx[0][k]) ; 
  if ((x<0) || (x>nline)) { 
   fatalx("bad index %d %d\n", x, nline) ;
  }
  xind[x] = x ; 
  xwt[x] = xx[1][k] ; 
  xvals[x] = xx[2][k] ; 
  xmax = x ; 
 }
 
  t = xind[0] ;         
  if (t!=0) fatalx("bad global row %d\n", t) ;

 mean = xvals[0] ;

 for (k=1; k<=xmax; ++k) { 
  jmean[k-1] = xvals[k] ;  
  jwt[k-1] =   xwt[k] ;
 }
 

  weightjack(&est, &sig, mean, jmean, jwt, xmax) ;

  if (vflag == NO) { 
   printf("## simpjack2: %s ", iname) ;
   printf("%12.6f %12.6f %12.6f  %9.3f", est, mean, sig, est/sig) ;
   printnl() ;
  }

  else { 
   printf("mean diff(a-b): %12.6f\n", mean) ;
   printf("Jackknife mean est: %12.6f\n", est) ;
   printf("std. err. %12.6f\n", sig) ;
   printf("Z: %9.3f\n", est/(sig+1.0e-20)) ; 
  }

  return 0 ; 
}


void readcommands(int argc, char **argv) 

{
  int i,haploid=0;

  while ((i = getopt (argc, argv, "i:m:v")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'm':
	xmean = atof(optarg) ;
	break;

      case 'v':
	vflag = YES ; 
	break;

      }

   }
}
