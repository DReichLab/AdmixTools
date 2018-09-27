#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <nicklib.h>  
#include <getpars.h>  

#include <globals.h>
#include <mcmcpars.h>
// #include "kimsubs.h"  
#include "egsubs.h" 

#include "badpairs.h"
#include "admutils.h"
#include "egsubs.h" 
#include "mcio.h"  
#include "qpsubs.h" 
char *trashdir = "/var/tmp" ;

int ztrans1mode = YES ;
int qtmode = NO ;

#define WVERSION 1216
#define MAXPOP 100
// #define MAXSTR 256 
//  haploid implemented  
//  actual counts written (last column) 

static int debug = 0 ;
int ranswitch = YES ;  // NO use full hypergeom, not random sample
int usedictnames = NO ;

char *specfilename = NULL ;       
char *specstring = NULL ;       
char *ascstring = NULL ;       

FILE *ofile ;
void readcommands(int argc, char **argv) ;

void loadssize(char *buff)  ;
void loadasc(char *buff) ;
int getstype(char *sx) ;
void countp(int ***counts, char **popnames, int numeg, SNP **xsnplist, Indiv **xindlist, int *xindex, int nrows, int ncols)  ;
void setspecp(double **xbsum, double *xbnum, int ***counts, int nblocks, int ncols, int npops, int *sampx1, SNP **xsnplist)  ;
void calcyp(double *ymean, double *ysig, double **ybsum, double *ybnum, int nblocks, int numx)  ;
void mkcv(int *cv, int *ascderiv, int *asctot, int kode, int *sampx1, int npops, int numeg)  ;
void setkodep(double *yprob, int **cc, int *sampx1, int npops)  ;
void dumpspec(char *specfilename,  double *specest, double *specsig)  ;  
void printabx(int abxmode)  ;

char **eglist ;  
char **popnames ;  
char *dictname = "dictfile" ;
char ***dictpops ; 
int numdict ;
int numeg = 0 ;
char **snames ; 
char **vnames ; 
int  numspec=0, xchrom = -1  ;
int abxmode = 0 ;  // 1 => AC CA AT
int abx(int a, int b)  ;
int abxok(int abx, int abxmode)  ;
int specsize = 0 ;
int flipmode = NO ; // ancestral counted  
double taest, tasig, *ttanum, *ttaden ;
double tbest, tbsig, *ttbnum, *ttbden ;

int *ascderiv, *asctot, *sampsize, *sampx1 ;
int *fascderiv, *fasctot, *tasctot, *maxacount, *maxacountp ;
int *ttspecpoly,  ttspecsize, *ww ;  
int *maxdeg ;

Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int npops ; 

char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *indivname = NULL ;
char *badsnpname = NULL ;
char *poplistname = NULL ;
char *haploidstring = NULL ; 
char **haploidpops ; 
int  numhaploidpops = 0 ; 

int **dekodetable ;      


double blgsize = 0.05 ;  // block size in Morgans */

int main  (int argc , char **argv) 
{

  SNP *cupt ;
  Indiv *indx ;

  int  i, j, k, t, u, a, b, n, x, tt, dd, numx, mx ; 
  double *specprobs, *zprobs, *yprobs, *ypvar, *ypsig ;
  int *specpoly, *ssize, numspec ;
  int *wexp ; 
  int numvind, nignore, numrisks = 1 ;

  SNP **xsnplist  ;
  Indiv **xindlist ;
  int *xindex ;
  int nrows, ncols, specnum, sval ;
  double *xsum, *xnum, *ysum, *ynum, *pprobs, *wwprobs, **xprobs, **y2probs ;
  int ***counts, **cc, ***counts2, **cc2 ;
  int *ktot, *kder, *vlist, kode, kode2, isbad ;
  int maxsize ;
  double y, y1, y2, ybase, yn, ysig ;

  int *blstart, *blsize, nblocks ;
  int  xnblocks ; /* for xsnplist */
  int *bcols, bnum ;
  double **xbsum, *xbnum ;
  double *wvar, *wa, *w1, *w2, yfixsig, yfixdiff ;  
  double dstat, dsig ;
  double *specest, *specsig ; 

  int xnumeg=0, t0 , t1 ; 
  char **xeglist ;
  double **xcsum, **xcnum ;
  char *sx ;
  double ya, yb, ya2, yb2, yy, tamean, tbmean ;
  double *ttanum, *ttaden, *ttbnum, *ttbden, *ttdiff, *ttnum ;
  int *specindex, *vsnps, *totsize ;
  int apoly[10], maxapoly, maxa[10] ; 
  double *apmeans, *aamean ;
  int numacount, *aapoly ;
  int *popsize ;

  int nbad = 0, ngood = 0 ;


  ZALLOC(fascderiv, MAXPOP, int) ;
  ZALLOC(fasctot, MAXPOP, int) ;
  ZALLOC(ascderiv, MAXPOP, int) ;
  ZALLOC(asctot, MAXPOP, int) ;
  ZALLOC(ww, MAXPOP, int) ;
  ZALLOC(sampsize, MAXPOP, int) ;
  ZALLOC(totsize, MAXPOP, int) ;
  ZALLOC(sampx1,  MAXPOP, int) ;
  ZALLOC(maxdeg,  MAXPOP, int) ;
  ZALLOC(eglist,  MAXPOP, char *) ;
  ZALLOC(dictpops,  2, char **) ;
  ZALLOC(dictpops[0],  MAXPOP, char *) ;
  ZALLOC(dictpops[1],  MAXPOP, char *) ;
  

  readcommands(argc, argv) ;

  cputime(0) ; 
  if (specstring == NULL) fatalx("no specstring!\n") ;
  if (ascstring != NULL) printf("ascertain: %s\n", ascstring) ;
  if (specstring != NULL) printf("sampsizes: %s\n", specstring) ;
  if (flipmode) printf("flipmode set!\n") ;
  if (ranswitch) printf("random subsampling!\n") ; 
  else printf("exact subsampling using hypergeometric\n") ;

  loadssize(specstring) ;  // place first
  npops = numeg ;  
  loadasc(ascstring) ;

/**
  for (x=0; x<numeg; ++x) { 
   printf("zspec: %3d %20s ", x, eglist[x]) ; 
   printf(" %3d %3d %3d", ascderiv[x], asctot[x], sampsize[x]) ; 
   printnl() ;
  }
*/

  if (haploidstring != NULL) { 
   ZALLOC(haploidpops, 10, char *) ; 
   subcolon(haploidstring) ; 
   numhaploidpops = splitup(haploidstring, haploidpops, 10) ; 
   printf("haploid pops:\n") ; 
   printstrings(haploidpops, numhaploidpops) ; 
  }

  numsnps = 
    getsnps(snpname, &snpmarkers, 0,  badsnpname, &nignore, numrisks) ;

// fakespacing 0.0 (default)

  numindivs = getindivs(indivname, &indivmarkers) ;
  setgenotypename(&genotypename, indivname) ;
  printf("genotypename:  %s\n", genotypename) ;

   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  if (abxmode>0) printabx(abxmode) ;


  for (i=0; i<numsnps; i++) { 
   cupt = snpmarkers[i] ;
   if (cupt -> chrom > 22) cupt -> ignore = YES ;
   if ((xchrom > 0) && (cupt -> chrom != xchrom)) cupt -> ignore = YES ;
   if (flipmode) flipalleles(cupt) ;
   if (abxmode == 99) continue ; 
   if (cupt -> alleles[0] == CNULL) continue ;
   a = base2num(cupt -> alleles[0]) ;
   b = base2num(cupt -> alleles[1]) ;
   t = abx(a, b) ; 
   if (t<0) cupt -> ignore = YES ;
   if (abxok(t, abxmode) == NO) cupt -> ignore = YES ;
// printf("zz2 %c %c %d %d\n", cupt -> alleles[0],  cupt -> alleles[1], t, cupt -> ignore) ;
  }

  numsnps = rmsnps(snpmarkers, numsnps, NULL) ;

  numdict =  numlines(dictname)  ;  // needed for next call
  numdict =  getnames(&dictpops, numdict, 2, dictname) ;

  ZALLOC(popnames, numeg, char *) ;
  for (k=0; k<numeg; ++k) { 
   t = indxstring(dictpops[0], numdict, eglist[k]) ;
   if (t<0) {  
    popnames[k] = strdup(eglist[k]) ;
   }
   else { 
    popnames[k] = strdup(dictpops[1][t]) ;
   }
  }
  for (k=0; k<numeg; ++k) { 
   printf("%3d: %15s %15s %3d %3d %3d\n", k, eglist[k], popnames[k], ascderiv[k], asctot[k], sampsize[k]) ;
  }

  ZALLOC(popsize, npops, int) ;

  for (k=0; k<npops; ++k) { 
    popsize[k] = t = setstatus(indivmarkers, numindivs, popnames[k]) ; 
    printf("popsize: %20s %4d\n", popnames[k], t) ; 
  }
  ivmaxmin(popsize, npops, NULL, &t) ; 
  if (t==0) fatalx("population with no data\\n") ; 

  for (i=0; i<numindivs; ++i) { 
    indx = indivmarkers[i] ; 
    if (indx -> affstatus == NO) {   
     indx -> ignore = YES ;
    } 
  }

  numindivs = rmindivs(snpmarkers, numsnps, indivmarkers, numindivs) ;

  numvind = numvalidind(indivmarkers, numindivs) ;
  printf("\n\n") ;
  printf("\n\n") ;
  printf("numindivs: %d valid: %d numsnps: %d nignore: %d\n" ,
    numindivs, numvind, numsnps, nignore) ; 

  printf("jackknife block size: %9.3f\n", blgsize) ;
  nblocks = numblocks(snpmarkers, numsnps, blgsize) ;
  printf("numsnps: %d\n", numsnps) ;
  nblocks = MIN(nblocks, 10000) ;
  ZALLOC(blstart, nblocks, int) ;
  ZALLOC(blsize, nblocks, int) ;
//  printf("number of blocks for block jackknife: %d\n", nblocks) ;

  ZALLOC(xindex, numindivs, int) ;
  ZALLOC(xindlist, numindivs, Indiv *) ;
  nrows = loadindx(xindlist, xindex, indivmarkers, numindivs) ;

  ZALLOC(xsnplist, numsnps, SNP *) ;
  ncols = loadsnpx(xsnplist, snpmarkers, numsnps, indivmarkers) ;

  printf("nrows: %d ncols: %d  blocks: %d\n", nrows, ncols, nblocks) ;

  setblocks(blstart, blsize, &xnblocks, xsnplist, ncols, blgsize)  ;

  ivvp(totsize, asctot, sampsize, numeg) ;

  numx = iprod(sampx1, npops) ;
  printf("numx: %d\n", numx) ;
  ZALLOC(vlist, npops, int) ;
  dekodetable = initarray_2Dint(numx, npops, 0) ;
  for (kode=0; kode<numx; ++kode) { 
    dekodeitbb(vlist, kode, npops, sampx1) ; 
    copyiarr(vlist, dekodetable[kode], npops) ;
  }
  
  ZALLOC(xsum, numx, double) ;
  ZALLOC(xnum, numx, double) ;
  ZALLOC(ysum, numx, double) ;
  ZALLOC(ynum, numx, double) ;
  ZALLOC(pprobs, numx, double) ;
  ZALLOC(wwprobs, numx, double) ;

  xbsum = initarray_2Ddouble(nblocks, numx, 0.0) ;
  ZALLOC(xbnum, nblocks, double) ;

  ZALLOC(counts, ncols, int **) ;          
  for (k=0; k<ncols; ++k) { 
   counts[k] = initarray_2Dint(npops, 2, 0) ;
  }

  countp(counts, popnames, npops, xsnplist, xindlist, xindex,  nrows, ncols) ;
  setspecp(xbsum, xbnum, counts,  nblocks, ncols, npops, sampx1, xsnplist)  ;
  sum2D(ysum, xbsum, nblocks, numx) ;

 ZALLOC(specest, numx, double) ;
 ZALLOC(specsig, numx, double) ;

 calcyp(specest, specsig, xbsum, xbnum, nblocks,  numx)  ;


 printf("sampsizes: ") ; 
 for (j=0; j<numeg; ++j)  { 
  printf(" %s :: %d ",  eglist[j], totsize[j]) ; 
  if (j<(numeg-1)) printf ("; ")  ;
 }
 printnl() ;

 fflush(stdout) ;

 for (k=0; k<numx; ++k) { 
  printf("specout: %3d: ", k) ; 
  mkcv(ww, ascderiv, asctot, k, sampx1, npops, numeg) ; 
  for (j=0; j<numeg; ++j) { 
   printf("%3d ", ww[j]) ;
  }
  printf(" %12.6f %12.6f",  specest[k], specsig[k]) ; 
  printf( "  %9.0f", ysum[k]) ;
  printnl() ;
  fflush(stdout) ;
 }

  dumpspec(specfilename,  specest, specsig) ; 

  printf("## end of gcount: cputime %9.3f seconds\n", cputime(1)) ;
//  printslurmenv() ;


  return 0 ;
  
}
void mkss(char *ss, int a, int sa, int b, int sb) 
{
  int t = a + b ; 

  sprintf(ss, "%d", t) ; 
  if (sa==0) return ;
  if (sb==0) return ;
  sprintf(ss, "%d:%d", a, b) ;


}
void mkcvs(char **ss, int *ascderiv, int *asctot, int kode, int *sampx1, int npops, int numeg) 
{

 int *vlist, j ; 

 ZALLOC(vlist, numeg, int) ; 
 dekodeitbb(vlist, kode, npops, sampx1) ; 

 for (j=0; j<numeg; ++j) { 
  mkss(ss[j], ascderiv[j], asctot[j], vlist[j], sampx1[j]-1) ;
 }

 free(vlist) ;

}

void mkcv(int *cv, int *ascderiv, int *asctot, int kode, int *sampx1, int npops, int numeg) 
{

 int *vlist ; 

 ZALLOC(vlist, numeg, int) ; 
 dekodeitbb(vlist, kode, npops, sampx1) ; 
 ivvp(cv, vlist, ascderiv, numeg) ; 

 free(vlist) ;

}
void dumpspec(char *specfilename,  double *specest, double *specsig)  
{

 FILE *fff ;
 int j, k, numx, *ww, *totsize ;
 char **sss, ss[20], *sx ;
 double  minsig, y  ;


 fflush(stdout) ;
 if (specfilename == NULL) return ; 
 openit(specfilename, &fff, "w") ;
 numx = iprod(sampx1, npops) ;
 ZALLOC(ww, numeg, int) ;
 ZALLOC(totsize, numeg, int) ;
 ZALLOC(sss, numeg, char *) ;

 for (j=0; j<numeg; ++j) { 
  ZALLOC(sss[j], 20, char) ; 
 }

 fprintf(fff, "newformat:\n") ;
 ivvp(totsize, asctot, sampsize, numeg) ;

 fflush(fff) ; 
 fprintf(fff, "sampsizes: ") ; 
 for (j=0; j<numeg; ++j)  { 
   mkss(ss, asctot[j], asctot[j], sampsize[j], sampsize[j]) ;
  sx = popnames[j] ; 
  if (usedictnames) sx = eglist[j] ; 
  fprintf(fff, " %s :: %s ",  sx, ss) ; 
  if (j<(numeg-1)) fprintf (fff, "; ")  ;
 }
 fprintf(fff, "\n") ;
 fflush(fff) ; 

 minsig = 1.0e10 ;
 for (k=0; k<numx; ++k) { 
  y = specsig[k] ; 
  if (y==0.0) continue ; 
  minsig = MIN(minsig, y) ;  
 }
// smalllest non-zero sigma
 
 for (k=0; k<numx; ++k) { 
  fprintf(fff, "spec: %5d: ", k) ; 
  mkcvs(sss, ascderiv, asctot, k, sampx1, npops, numeg) ; 
  for (j=0; j<numeg; ++j) { 
   fprintf(fff, "%5s ", sss[j]) ;
  }
  y = MAX(minsig, specsig[k]) ; 
  fprintf(fff, " %15.9f %15.9f\n",  specest[k], y) ; 
  fflush(fff) ;
 }

 fclose(fff) ;
 freeup(sss, numeg) ;


}
void setspecp(double **xbsum, double *xbnum, int ***counts, int nblocks, int ncols, int npops, int *sampx1, SNP **xsnplist) 
{  

 int j, k, t, bnum, *vlist ; 
 int kode, numkode, cder, canc, ctot, nsamp ; 
 SNP *cupt ; 
 int **cc ; 
 double y, *yprob ;  
 int nbl = 0,  nogood ; 

 numkode = iprod(sampx1, npops) ;
 ZALLOC(yprob, numkode, double) ;
 vzero(xbnum, nblocks) ; 
 clear2D(&xbsum, nblocks, numkode, 0) ; 
 for (j=0; j<ncols; ++j) {  
   cupt = xsnplist[j] ; 
   if (cupt -> ignore) continue ;  
   bnum = xsnplist[j] -> tagnumber ;
   nbl = MAX(nbl, bnum+1) ; 
   if (bnum < 0) continue ; 
   if (bnum >= nblocks) fatalx("logic bug: %d %d\n", bnum, nblocks) ;
   vzero(yprob, numkode) ; 
   cc = counts[j] ; 
   nogood = NO ;       
   for (k=0; k<npops; ++k) { 
    cder = cc[k][0] ;  
    canc = cc[k][1] ;  
    ctot = cder + canc ;
    nsamp = sampx1[k] - 1 ;  
    if (nsamp > ctot) nogood = YES  ; // infeasible
   } 
   if (nogood) continue ;       
  setkodep(yprob, cc, sampx1,  npops) ; // main work
  vvp(xbsum[bnum], yprob, xbsum[bnum], numkode) ;  
  xbnum[bnum] += asum(yprob, numkode) ; 
 }
 free(yprob) ; 
}
void oldsetkodeq(double *yprob, int **cc, int *sampx1, int npops) 
// obsolescent.  setkodeq calls hprobv
{

   int kode, k,  cder, canc, ctot, nsamp, t, *vlist ; 
   double **lprobs, *ww, y ;
   int n, a, m, j, numx ;

   numx = iprod(sampx1, npops) ; 
   t = intsum(sampx1, npops) ;   
   ZALLOC(ww, npops, double) ;
   lprobs = initarray_2Ddouble(npops, t, 0) ; 
   for (k=0; k<npops; ++k) { 
    cder = cc[k][0] ;  
    canc = cc[k][1] ;  
    n = ctot = cder + canc ;
    a = cder ; 
    m = sampx1[k]-1 ;    
    for (j=0; j<=m; ++j) { 
     lprobs[k][j] = loghprob(n, a, m, j) ;
    }
   }
   for (kode=0; kode<numx; ++kode) { 
    vlist = dekodetable[kode] ;
    for (j=0; j<npops; ++j) { 
     t = vlist[j] ; 
     ww[j] = lprobs[j][t] ;
    }
    y = asum(ww, npops) ; 
    yprob[kode] = exp(y) ;
   }

   free2D(&lprobs, npops) ;
   free(ww) ;

}

void setkodeq(double *yprob, int **cc, int *sampx1, int npops) 
{

   int kode, k,  cder, canc, ctot, nsamp, tnum, t, *vlist, isbad ; 
   double **vprobs, *ww, y ;
   int n, a, m, j, numx ;

   numx = iprod(sampx1, npops) ; 
   ivmaxmin(sampx1, npops, &tnum, NULL) ;
   ZALLOC(ww, npops, double) ;
   vprobs = initarray_2Ddouble(npops, tnum, 0) ; 
   isbad = NO ; 
   for (k=0; k<npops; ++k) { 
    cder = cc[k][0] ;  
    canc = cc[k][1] ;  
    n = ctot = cder + canc ;
    a = cder ; 
    m = sampx1[k]-1 ;    
    t = hprobv(vprobs[k], n, a, m) ;
    if (t<0) isbad = YES ;  
   }
   for (kode=0; kode<numx; ++kode) { 
    vlist = dekodetable[kode] ;
    for (j=0; j<npops; ++j) { 
     t = vlist[j] ; 
     ww[j] = vprobs[j][t] ;
    }
    y = aprod(ww, npops) ; 
    yprob[kode] = y ;  
   }

   free(ww) ;
   if (isbad) {
    vzero(yprob, numx) ;
   }
   else { 
    bal1(yprob, numx) ;
   }
   free2D(&vprobs, npops) ;

}

void setkodep(double *yprob, int **cc, int *sampx1, int npops) 
{

   int kode, k,  cder, canc, ctot, nsamp, t, *vlist ; 
   int numx ;
   int numcall=0, numbad = 0 ;
   double *w1, *w2, y1, y2 ; 

   ++numcall ;
   numx = iprod(sampx1, npops) ; 
   vzero(yprob, numx) ;
   if (ranswitch == NO) {
     setkodeq(yprob, cc, sampx1,  npops) ;
     return ;
   }
   ZALLOC(vlist, npops, int) ; 

   for (k=0; k<npops; ++k) { 
    cder = cc[k][0] ;  
    canc = cc[k][1] ;  
    ctot = cder + canc ;
    nsamp = sampx1[k] - 1 ;  
    if (nsamp > ctot) continue ; // infeasible
    
    t = ranhprob(ctot, cder, nsamp) ; 
    if (t<0) continue ; 
    vlist[k] = t  ;  
  }
  kode = kodeitbb(vlist, npops, sampx1) ; 
  yprob[kode] = 1 ;  // if probabilistic must normalize here

  free(vlist) ; 

}

void countp(int ***counts, char **popnames, int numeg, SNP **xsnplist, Indiv **xindlist, int *xindex, int nrows, int ncols)           
// status must be set
// counts initialized to zero
{
  int *xtypes ;
  int k, a, i, j  ; 
  Indiv *indx ;
  int *ishaploid ; 
  char *sx ; 

  ZALLOC(xtypes, nrows, int) ;
  ZALLOC(ishaploid, nrows, int) ;
  
  for (a=0; a<numeg; ++a) { 
    sx = popnames[a] ;   
    k = indxindex(haploidpops, numhaploidpops, sx) ;               
    if (k>=0) { 
     ishaploid[a] = 1 ; 
     printf("haploid %s set!\n", sx) ; 
    }
  }



  for (i=0; i<nrows; i++) {
   indx = xindlist[i] ;
   k = indxindex(popnames, numeg, indx -> egroup) ;
   xtypes[i] = k ;
  }

  countpops(counts, xsnplist, xindex, xtypes, nrows, ncols) ;

  for (a=0; a<numeg; ++a) { 
   if (ishaploid[a] == 0) continue ; 
   printf("setting counts for %s haploid\n", popnames[a]) ;  
   for (j=0; j<ncols; ++j) { 
    counts[j][a][0] /=2 ;
    counts[j][a][1] /=2 ;
   }
  }

  free(xtypes) ;
  free(ishaploid) ;
}


int addpop(char *pop)  
{
 int t ;
 if (numeg==0) t=-1 ;
 else t = indxstring(eglist, numeg, pop) ; 
 if (t>=0) return t ;

 eglist[numeg] = strdup(pop) ; 
 t = numeg ;
 ++numeg ;
 return t ;
}


void readcommands(int argc, char **argv) 

{
  int i ;
  char *parname = NULL ;
  phandle *ph ;
  int n ;

  while ((i = getopt (argc, argv, "p:a:s:fvV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'a':
	ascstring = strdup(optarg) ;
	break;

      case 's':
	specstring = strdup(optarg) ;
	break;

      case 'o':
	specfilename = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %d\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getlongstring(ph, "ascertain:", &ascstring) ;
   getlongstring(ph, "sampsizes:", &specstring) ;
   getstring(ph, "oname:", &specfilename) ;
   getstring(ph, "specfile:", &specfilename) ;
   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "badsnpname:", &badsnpname) ;
   getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "dictname:", &dictname) ;
   getint(ph, "chrom:", &xchrom) ;
   getint(ph, "abxmode:", &abxmode) ;  // may take more than 2 values 
   getdbl(ph, "blgsize:", &blgsize) ;
   getlongstring(ph, "haploidpops:", &haploidstring) ;
   getint(ph, "ranswitch:", &ranswitch) ;
   getint(ph, "flipmode:", &flipmode) ;
   getint(ph, "usedictnames:", &usedictnames) ; // on output Dictionary shorthand used 
   
   writepars(ph) ;
   closepars(ph) ;

}

void loadfixedasc(char *buff)  
{
  char *spt[MAXFF], *sx, *s2[10] ;
  int nsplit, nasc, k, z, xd, xt, n ; 

  nasc = nsplit = splitupx(buff, spt, MAXFF, ';') ;
  for (k=0; k<nasc; ++k) {  
   sx = strdup(spt[k]) ;
   substring(&sx, ":", " ") ;
   n = splitup(sx, s2, 10) ;
   if (n != 3) fatalx("bad loadfixedasc :: %s\n", buff) ;
   z = addpop(s2[0]) ;
   xd = atoi(s2[1]) ;
   xt = atoi(s2[2]) ;
   fascderiv[z] = xd ;
   fasctot[z] = xt ;
   freeup(s2, n) ;
   freestring(&sx) ;
  }

  freeup(spt, nsplit) ;
}
void loadasc(char *buff)  
{
  char *spt[MAXFF], *sx, *s2[10] ;
  int nsplit, nasc, k, z, xd, xt, n ; 

  if (buff==NULL) return ; 
  nasc = nsplit = splitupx(buff, spt, MAXFF, ';') ;
  for (k=0; k<nasc; ++k) {  
   sx = strdup(spt[k]) ;
   substring(&sx, ":", " ") ;
   n = splitup(sx, s2, 10) ;
   if (n != 3) fatalx("bad loadasc :: %s\n", buff) ;
   z = addpop(s2[0]) ;
   xd = atoi(s2[1]) ;
   xt = atoi(s2[2]) ;
   ascderiv[z] = xd ;
   asctot[z] = xt ;
   freeup(s2, n) ;
   freestring(&sx) ;
  }

  freeup(spt, nsplit) ;
}
void loadssize(char *buff)  
{
  char *spt[MAXFF], *sx, *s2[10] ;
  int nsplit, ns, k, z, siz,  n ; 

  ns = nsplit = splitupx(buff, spt, MAXFF, ';') ;
  for (k=0; k<ns; ++k) {  
   sx = strdup(spt[k]) ;
   substring(&sx, ":", " ") ;
   n = splitup(sx, s2, 10) ;
   if (n != 2) fatalx("bad loadasc :: %s\n", buff) ;
   z = addpop(s2[0]) ;
   siz = atoi(s2[1]) ;
   sampsize[z] = siz ;
   freeup(s2, n) ;
   freestring(&sx) ;
  }
  ivsp(sampx1, sampsize, 1, numeg) ;
  specsize = iprod(sampx1, numeg) ;

  freeup(spt, nsplit) ;
  

}


int getstype(char *sx) 
{
 int t ;
 t = strcmp(sx, "ascertain:") ; if (t==0) return 1 ;
 t = strcmp(sx, "sampsizes:") ; if (t==0) return 2 ;
 t = strcmp(sx, "fixedasc:") ; if (t==0) return 3 ;
 return -1 ;

}

void print_polyn(char **vnames, int *poly, double *cpoly, int n, int numterm)
{
 int i, j, t ;

 printf("%6s ", "") ;
 for (j=0; j<n; j++) {
  printf("%3s ", vnames[j]) ;
 }
 printnl() ;

 for (i=0; i<numterm; i++) {
  printf("%5d: ", i) ;
  for (j=0; j<n; j++) {
   t = poly[i*n+j] ;
   printf("%3d ", t) ;
  }
  printf ("  %12.6f", cpoly[i]) ;
  printnl() ;
 }

}

void printabx(int abxmode) 
{
  int t, a, b, abxkode, ca, cb ; 

  printf("printabx: %3d\n", abxmode) ; 
  for (a=0; a<=3; ++a) { 
   for (b=0; b<=3; ++b) { 
    if (a==b) continue ; 
    abxkode = abx(a,b) ; 
    printf("abxvalid: %c %c ", num2base(a), num2base(b)) ; 
    t = abxok(abxkode, abxmode) ; 
    if (t==YES) printf("YES") ;
    if (t==NO)  printf(" NO") ;
    printnl() ;  


  }} 


}

int abx(int a, int b) 
{

 if (a<0) return -1 ;
 if (b<0) return -1 ;
 if (a==b) return -1 ;
 if (a>1) return abx(3-a, 3-b) ;

 if (a==0) return b-1 ; 
 if (b==0) return 3 ;  // CA
 if (b==2) return 4 ;  // CG
 if (b==3) return 5 ;  // CT

 fatalx("badbug\n") ;
 
}

int abxok(int abx, int abxmode) { 

 int t ;

 if (abxmode == 99) return YES ;
 if (abxmode >= 10) { 
  t = abxmode - 10 ;
  if (abx == t) return YES ;
  return NO ;
 }
   
 switch (abxmode)  { 
 case 0:  
  return YES ;
 case 1:
  if (abx==0) return YES ; // AC
  if (abx==2) return YES ; // AT
  if (abx==3) return YES ; // CA
  return NO ;
 case 2: 
//No AG, CT
  if (abx==1) return NO ; 
  if (abx==5) return NO ; 

  return YES ;
 default:  
  fatalx("abxmode %d not implemented\n", abxmode) ;
 }
}


void calcyp(double *ymean, double *ysig, double **ybsum, double *ybnum, int nblocks, int numx) 
{
  double *jmean, **jjmean, *jwt, *ysum, ynum, *tsum, tnum ;
  double *ypvar, y; 
  int k ; 


  ZALLOC(jmean, numx, double) ; 
  ZALLOC(jwt, nblocks, double) ; 

  ZALLOC(ysum, numx, double) ; 
  ZALLOC(tsum, numx, double) ; 
  ZALLOC(ypvar, numx*numx, double) ;

  jjmean = initarray_2Ddouble(nblocks, numx, 0.0) ;
  fflush(stdout) ;

  sum2D(ysum, ybsum, nblocks, numx) ;
  ynum = asum(ybnum, nblocks) + 1.0e-20; 

  vst(jmean, ysum, 1.0/ynum, numx) ;
  bal1(jmean, numx) ;

  fflush(stdout) ;
  for (k=0; k<nblocks; ++k) { 
   vvm(tsum, ysum, ybsum[k], numx) ;
   tnum = ynum - ybnum[k] ; 
   jwt[k] = ybnum[k] ;
   vst(jjmean[k], tsum, 1.0/tnum, numx) ;
   bal1(jjmean[k], numx) ;
   fflush(stdout) ;
  }
  fflush(stdout) ;
 
  wjackvest(ymean, ypvar, numx,  jmean, jjmean, jwt, nblocks) ; 
  for (k=0; k<numx; ++k) { 
    y = ypvar[k*numx+k] ; 
    ysig[k] = sqrt(y) ;
  }

  free(jmean) ;
  free(jwt) ;
  free(ysum) ;
  free(tsum) ;
  free2D(&jjmean, nblocks) ;
  free(ypvar) ;

}

