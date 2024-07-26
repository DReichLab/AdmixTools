//nclude <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <libgen.h>

#include <nicklib.h>
#include <getpars.h>

#include "admutils.h"
#include "mcmcpars.h"
#include "regsubs.h"
#include "egsubs.h"
#include "qpsubs.h"
#include "eigsubs.h"
#include "globals.h"

#define WVERSION   "180"

// printsd added 
// clinetest added (as in qpdslow) Makes most sense for 2 f4 stats
// globaltest added  

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *trashdir = "/var/tmp";
int details = NO;
int hires = NO;			/* %10.4f */
int qtmode = NO;
int globaltest = YES ;  
int globalforce = NO ;  

int printsd = NO ; 
int clinetest = NO ; 

char *f3name = NULL;
char *fstatsname = NULL;

int seed = 0;
int popsizelimit = -1;

double diag = 0.0;

int isinit = NO;

char *poplistname = NULL;
char *popfilename = NULL;
char *fmvoutname = NULL;


char *outpop = NULL;
// list of outliers
char *basepop = NULL;
int basenum = -1;
double baseval = 0.0;

// outnum used for weights 
// basenum for f3 status  


double lambdascale = -1.0;
int *f2ind, *ind2f;
double *vest, *vvar, *vvinv, *xvvar, *xvvinv;
int nh2, numeg;
int *ezero = NULL;
double wtmin = .0001;
double minvar = 0.0;		// minvalue for variance term
int xnumeg;


char *outputname = NULL;
char *weightname = NULL;
FILE *ofile;
char **eglist;
char **egshort;
char **enames;

void readcommands (int argc, char **argv);
char *getshort (char *ss, int n);
void map4x (double *aa, double *bb, int n2, int *indx);
void map4y (double *aa, double *bb, int n2, int *indx);
void getmv (int a, int b, int c, int d, double *mean, double *var,
	    double *xest, double *xvar);
void bumpm (double *y, int a, int c, double val, double *xest);
void bumpw (double *w, int a, int c, double val);
double calcxx (double *xxans, double *qmat, double *ppwts, double *rhs,
	       int nrow, int ncol);
void setwww (double **tmix, double *www, int n);
void getwww (double **tmix, double *www, int n);
int iscanon (int a, int b, int c, int d);

void dump1 (FILE * dumpfile, double *ww, int n);
void loadpars (char *loadname, double *www, int nwts, double *xxans,
	       int nedge);
void read1 (FILE * loadfile, double *ww, int n);
void print4 (double *ff3, int a, int b, int c, int d, int numeg);
void printf3 (char *sss, FILE * fff, double *ww, char **eglist, int n);
int listsubset (int **x, int n, int k);
void balw (double **ww, int **vv, int n, int *nw);
void estff3 (double *fv, double *v, int numv, int *elist, int n);
void rcsquish (double *xmat, double *mat, int *cols, int oldn, int newn);
double ff4val (double *ff3, int a, int b, int c, int d, int numeg);
void dumpit (char *dumpname, double *ff3, double *ff3var, char **eglist,
	     int numeg);
void checkpd (double *a, int n);
void dumpf3 (char *dumpf3name, double **btop, double **bbot, int nblock);
int usage (char *prog, int exval);
void load4(int *x, int a, int b, int c, int d)  ;
void loadco(double *co, int *fs, int *fsb, int np) ;
void loadfsindex(int **fsindex, char ***qlist, char **eglist, int numfs, int numeg)  ;

int
main (int argc, char **argv)
{

  char sss[MAXSTR];
  int i, j, k, k1, k2, k3, k4, kk ;
  double y1, y2, y, sig, tail, yy1, yy2, ychi, ychi0;
  double *lambda, *evecs ; 

  int t, num, n, n1, n2, n3, n4;

  int nindiv = 0, e, f, lag = 1;
  int a, b, c, d, col, u, v, x;
  int *popsizes;
  double *qpscore;

  int popx[4];
  double *f2, *f2sig, *fst;
  double *ff3, *gg3, *ff3var, *ff3fit, *ff3sig;

  double *f3, *f3sig, *www, *ww2;
  int ng2, ng3, ng4;
  FILE *fff = NULL ;
  double scale, eps = 1.0e-20 ;

  int **findex, np, *fsnum, nfstats ; 
  int *tt, *jwork ; 
  int **basis, *basisnum,  *basisfn,  nbasis ; 
  double **fcoeffs ; 
  double *fsmean, *fssig ; 
  double *fbmean, *fbcovar ; 
  int numfs   ;       // number of entries in poplistname  
  int nqlist = 0  ; 
  char ***qlist ; 
  double *fsm, *fsv, *fsvinv, *fszero ; 
  int **fsindex ; 
  double *work, *yco ;     
  double *w1, *w2 ; 
  char *sx ; 
  int dof ; 
  char ***plists;
  int nplist ; 
  double *ww ; 
  double zmax = 0.0 ;
  double ymem ; 


  readcommands (argc, argv);
  printf ("## qpfmv version: %s\n", WVERSION);
  if (parname == NULL)
    return 0;


  cputime(0) ;
  calcmem(0) ;
  

  ZALLOC(eglist, MAXPOPS, char *) ; 
  numeg = np = fstats2popl(fstatsname, eglist) ; 
  if (numeg < 0) fatalx("no fstats file!\n") ;
  t = np*np ; 
 
  if (globalforce) printf("globalforce set\n") ;

  ZALLOC(ff3, np*np, double) ; 
  ZALLOC(ff3var, np*np*np*np, double) ; 

  loadfstats(fstatsname, ff3, ff3var, eglist, numeg) ; 
  printf("numeg: %d\n", numeg) ; 
  if (popfilename != NULL) nqlist = numlines (popfilename);
  if (nqlist==0) { 
   printf("no fstats to compute!\n") ; 
   return 0 ; 
  }
    ZALLOC (qlist, 4, char **);
    for (k = 0; k < 4; ++k) {
      ZALLOC (qlist[k], nqlist, char *);
    }
    numfs =  getnames(&qlist, nqlist, 4, popfilename) ; 
    printf ("number of fstats  %d\n", numfs);


  nh2 = numeg*(numeg-1) ; nh2 /= 2 ; 
  ZALLOC(vest, nh2, double) ; 
  ZALLOC(vvar, nh2*nh2, double) ;

  setvv(vest, vvar, ff3, ff3var, ind2f,  numeg) ; 
 
  fsindex = initarray_2Dint(numfs, 4, -1) ; 
  loadfsindex(fsindex, qlist, eglist, numfs, numeg) ; 
  ZALLOC(yco, np*np*numfs, double) ; 
  
 // printimat2D(fsindex, numfs, 4) ; 

  ZALLOC(fsm, numfs, double) ; 
  ZALLOC(fszero, numfs, double) ; 
  ZALLOC(fsv, numfs*numfs, double) ; 
  ZALLOC(fsvinv, numfs*numfs, double) ; 
  ZALLOC(work, numfs*numfs, double) ; 
  ZALLOC(jwork, numfs*numfs, int) ; 

  nbasis = mkcoeffs(yco, fsindex, numeg, numfs) ; 

  ZALLOC(w1, numfs, double) ;
  ZALLOC(w2, numfs, double) ;
  ZALLOC(lambda, numfs, double) ;
  ZALLOC(evecs, numfs * numfs, double) ;
 
  xmultx(work, yco, numfs, nbasis) ; 
  eigvals(work, lambda, numfs) ; 
  
  bal1(lambda, numfs) ;
  vst(lambda, lambda, (double) numfs, numfs) ;

/*
  printf("evals:\n") ; 
  printmat(lambda, 1, numfs) ; 
*/ 

  dof = 0 ; 
  for (k=0; k<numfs; ++k) { 
   if (lambda[k] < .0001) break ; 
   ++dof ; 
  }
  printf("dof:: %d\n", dof) ;

  fflush(stdout) ; 

  vv2ww(fsm, fsv, vest, vvar, numeg, fsindex, numfs) ; 

/* 

  vst(work, fsm, 1000, numfs) ; 
  fixit(jwork, work, numfs) ; 
  printimatw(jwork, 1, numfs, numfs) ; 
  printnl() ; 
  vst(work, fsv, 1000*1000, numfs*numfs) ; 
  fixit(jwork, work, numfs*numfs) ; 
  printimatw(jwork, numfs, numfs, numfs) ; 

  y = trace(fsv, numfs) ; 
  printf("trace: %15.9f\n", y) ; 
  vclear(work, y*.0001, numfs) ; 
  adddiag(fsv, work, numfs) ; 
  
  printmatwl(fsm, 1, numfs, numfs) ; 
  printnl() ;
  printmatwl(fsv, numfs, numfs, numfs) ; 


*/
  
  fflush(stdout) ; 

  printnl() ; 
  fff = stdout ; 
  if (fmvoutname != NULL) openit(fmvoutname, &fff, "w") ; 
  fprintf(fff,"%8s ", "" ) ;
  fprintf(fff,"%20s ", "P1" ) ;
  fprintf(fff,"%20s ", "P2" ) ;
  fprintf(fff,"%20s ", "P3" ) ;
  fprintf(fff,"%20s ", "P4" ) ;
  fprintf(fff,"%12s ", "fstat") ; 
  if (printsd) fprintf(fff,"%12s ", "s.err") ; 
  fprintf(fff,"%9s", "Z") ; 
  fprintf(fff,"\n") ; 
  for (k=0; k<numfs; ++k) { 
   y1 = fsm[k] ; 
   y2 = fsv[k*numfs+k] ; 
   y = y1/sqrt(y2+1.0e-20) ; 
   fprintf(fff, "%8s ", "result: ") ; 
   for (j=0; j<4; ++j) { 
     t = fsindex[k][j] ; 
     sx = eglist[t] ; 
     fprintf(fff, "%20s ", sx) ;
   }
// printimatxfile(fsindex[k], 1, 4, fff)  ;     
   fprintf(fff,"%12.6f ", y1) ;
   if (printsd) fprintf(fff, "%12.6f ", sqrt(y2)) ;
   fprintf(fff, "%9.3f ", y) ;
   fprintf(fff, "\n") ; 
   zmax = MAX(zmax, fabs(y)) ; 
   t = k % 100 ; if (t==0) fflush(fff) ;  
  }
   fprintf(fff, "\n") ; 
   fflush(fff) ; 

 if ((zmax > 10.0) && (globalforce == NO))  {  
  printf("## There are large Z-scores:: no globaltest\n") ;
  globaltest = NO ;
 }

 if (globaltest) { 
  eigvecs(fsv, lambda, evecs, numfs) ; 
  copyarr(lambda, w1, numfs) ; 
  bal1(lambda, numfs) ;
  vst(lambda, lambda, (double) numfs, numfs) ;

  if (verbose) {
    printf("evals:\n") ; 
    printmat(lambda, 1, numfs) ; 
  }

  ychi = 0.0 ; 
  for (a=0; a<dof; ++a) { 
   y  = vdot(evecs+a*numfs, fsm, numfs) ; 
   ychi += y*y / w1[a] ; 
  }
  y1 = rtlchsq(dof, ychi) ; 
  fprintf(fff, "##Hotelling T2: %9.3f dof: %d  tail: %12.6f\n", ychi , dof, y1) ; 
 }

 nplist = numfs ; 
 ZALLOC(ww, nplist, double) ; 
 ZALLOC(plists, nplist, char**) ;  
  for (a=0; a< nplist ; ++a) {
    ZALLOC(plists[a], 4, char *) ; 
    for (j=0; j<4; ++j) { 
      t = fsindex[a][j] ; 
      sx = eglist[t] ; 
      plists[a][j] = strdup(sx) ;
    }
  }
 
  
 for (a=0; a< nplist ; ++a) {
    if (clinetest == NO) break ;
     for (b=a+1; b< nplist ; ++b) {
      printf("clinetest:: ") ;
      printstringsx(plists[a], 4) ;
      printf(" :: ") ;
      printstringsx(plists[b], 4) ;
      printnl() ;

    for (k=0; k<=100; ++k) {
     vzero(ww, nplist) ; 
     y = ww[a] = (double) k / 100.0 ;
     ww[b] = 1.0-y ;
     y1 = vdot(fsm, ww, nplist) ;
     y2 = scx(fsv, NULL, ww, nplist) ;  y2 = sqrt(y2) ;

     printf("mixtable: %9.3f ", y) ;
     printf(" %12.6f %12.6f %9.3f\n", y1, y2, y1/y2) ;

    }
  }} 



  ymem = calcmem(1)/1.0e6 ;
  printf("##end of qpfmv: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;

  return 0 ; 
  np = numeg ; 

  return 0;

}


void
estff3 (double *fv, double *v, int numv, int *elist, int n)
// extract covarance for pops
{
  int a, b, xa, xb;

  for (a = 0; a < n; a++) {
    for (b = 0; b < n; b++) {
      xa = elist[a];
      xb = elist[b];
      fv[a * n + b] = v[xa * numv + xb];
    }
  }
}

void
print4 (double *ff3, int a, int b, int c, int d, int numeg)
{
  double y;

  y = ff4val (ff3, a, b, c, d, numeg);
  printf ("%4s ", get3 (egshort[a]));
  printf ("%4s ", get3 (egshort[b]));
  printf ("   ");
  printf ("%4s ", get3 (egshort[c]));
  printf ("%4s ", get3 (egshort[d]));
  printf ("%6d", nnint (1000.0 * y));
  printnl ();
}

double
ff4val (double *ff3, int a, int b, int c, int d, int numeg)
{
  double y1, y2, y3, y4;

  y1 = dump2 (ff3, a, c, numeg);
  y2 = dump2 (ff3, a, d, numeg);
  y3 = dump2 (ff3, b, c, numeg);
  y4 = dump2 (ff3, b, d, numeg);

  return y1 + y4 - (y2 + y3);

}


void
read1 (FILE * loadfile, double *ww, int n)
{

  int k;
  vclear (ww, -1.0, n);
  for (k = 0; k < n; ++k) {
    fscanf (loadfile, "%lf\n", &ww[k]);
    printf ("zzfs %d %9.3f\n", k, ww[k]);
  }

}


void
printf3 (char *sss, FILE * fff, double *ww, char **eglist, int n)
{
  int i, j, x;

  if (fff == NULL)
    return;
  fprintf (fff, "%s\n", sss);
  fprintf (fff, " %4s", "   ");
  for (i = 0; i < n; i++) {
    fprintf (fff, " %4s", get3 (eglist[i]));
  }
  fprintf (fff, "\n");
  for (i = 0; i < n; i++) {
    fprintf (fff, "%4s", get3 (eglist[i]));
    for (j = 0; j < n; j++) {
      x = nnint (1000 * ww[i * n + j]);
      fprintf (fff, " %4d", x);
    }
    fprintf (fff, "\n");
  }
  fprintf (fff, "\n");
}

void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n, t;

  while ((i = getopt (argc, argv, "p:o:vVh")) != -1) {

    switch (i) {

    case 'h':
      usage(basename(argv[0]), 0);
      break ; 

    case 'p':
      parname = strdup (optarg);
      break;

    case 's':
      seed = atoi (optarg);
      break;

    case 'v':
      printf ("version: %s\n", WVERSION);
      break;


    case 'V':
      verbose = YES;
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }


  if (parname == NULL) {
    fprintf (stderr, "no parameters\n");
    return;
  }

  pcheck (parname, 'p');
  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);
  dostrsub (ph);

//  getstring (ph, "poplistname:", &poplistname);
  getstring (ph, "fstatsname:", &fstatsname);
  getstring (ph, "popfilename:", &popfilename);
  getstring (ph, "fmvoutname:", &fmvoutname);

  getint (ph, "seed:", &seed);
  getint (ph, "printsd:", &printsd);
  getint (ph, "clinetest:", &clinetest) ; 
  getint (ph, "globaltest:", &globaltest) ; 
  getint (ph, "globalforce:", &globalforce) ; 

  printf ("### THE INPUT PARAMETERS\n");
  printf ("##PARAMETER NAME: VALUE\n");
  writepars (ph);

}
void
map4y (double *aa, double *bb, int n2, int *indx)
// map 4d array (n1 x n1 x n1 x n1  -> b  n2 x n2 x n2 x n2 
{
  int u, v, a, b, c, d;
  int x;
  double y1;
  int nh2;

  nh2 = n2 * (n2 - 1);
  nh2 /= 2;
  vzero (aa, nh2 * nh2);
  for (u = 0; u < nh2; ++u) {
    for (v = u; v < nh2; ++v) {

      x = indx[u];
      a = x / n2;
      b = x % n2;
      x = indx[v];
      c = x / n2;
      d = x % n2;

      y1 = dump4 (bb, a, b, c, d, n2);
      aa[u * nh2 + v] = y1;
      aa[v * nh2 + u] = y1;
    }
  }
}

void
bumpm (double *y, int a, int c, double val, double *xest)
{

  int k;

  if (a == basenum)
    return;
  if (c == basenum)
    return;
  k = f2ind[a * numeg + c];
  *y += val * xest[k];

}

void
bumpw (double *w, int a, int c, double val)
{

  int k;

  if (a == basenum)
    return;
  if (c == basenum)
    return;
  k = f2ind[a * numeg + c];
  w[k] += val;

}

void
dumpf3 (char *dumpf3name, double **btop, double **bbot, int nblock)
{
#define SQ  40
  int fdes, ret;
  int k, ng2, ng4;
  char sss[SQ];

  if (dumpf3name == NULL)
    return;

  ridfile (dumpf3name);
  fdes = open (dumpf3name, O_CREAT | O_TRUNC | O_RDWR, 006);
  cclear ((unsigned char *) sss, CNULL, SQ);
  sprintf (sss, "numeg: %d nblock: %d nh2: %d", numeg, nblock, nh2);
  ret = write (fdes, sss, SQ * sizeof (char));
  if (ret < 0) {
    perror ("write failure");
    fatalx ("bad write:  %s", sss);
  }
  for (k = 0; k < numeg; ++k) {
    cclear ((unsigned char *) sss, CNULL, SQ);
    strncpy (sss, eglist[k], SQ);
    ret = write (fdes, sss, SQ * sizeof (char));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  %s", sss);
    }
  }
  for (k = 0; k < nblock; ++k) {
    ret = write (fdes, btop[k], nh2 * sizeof (double));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  block %d", k);
    }
    ret = write (fdes, bbot[k], nh2 * sizeof (double));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  block %d", k);
    }
  }
  close (fdes);
  chmod (dumpf3name, 0644);
}

void
dumpit (char *dumpname, double *ff3, double *ff3var, char **eglist, int numeg)
{
#define SQ  40
  int fdes, ret;
  int k, ng2, ng4;
  char sss[SQ];
  if (dumpname == NULL)
    return;

  ridfile (dumpname);
  fdes = open (dumpname, O_CREAT | O_TRUNC | O_RDWR, 006);
  cclear ((unsigned char *) sss, CNULL, SQ);
  sprintf (sss, "numeg: %d", numeg);
  ret = write (fdes, sss, SQ * sizeof (char));
  if (ret < 0) {
    perror ("write failure");
    fatalx ("bad write:  %s", sss);
  }
  for (k = 0; k < numeg; ++k) {
    cclear ((unsigned char *) sss, CNULL, SQ);
    strncpy (sss, eglist[k], SQ);
    ret = write (fdes, sss, SQ * sizeof (char));
    if (ret < 0) {
      perror ("write failure");
      fatalx ("bad write:  %s", sss);
    }
  }
  ng2 = numeg * numeg;
  ng4 = ng2 * ng2;
  ret = write (fdes, ff3, ng2 * sizeof (double));
  if (ret < 0) {
    perror ("write failure");
  }
  ret = write (fdes, ff3var, ng4 * sizeof (double));
  close (fdes);
  chmod (dumpname, 0644);

}

void
checkpd (double *a, int n)
{
  double *b, *d;
  ZALLOC (b, n * n, double);
  ZALLOC (d, n, double);
  vsp (b, a, 1.0e-20, n * n);
  vst (b, b, 1.0e10, n * n);
  getdiag (d, b, n);
  printmat (d, 1, n);
  choldc (b, n, b);
  printnl ();
  printnl ();
  getdiag (d, b, n);
  printmat (d, 1, n);
  free (b);
  free (d);
}

int usage (char *prog, int exval)
{

  (void)fprintf(stderr, "Usage: %s [options] <file>\n", prog);
  (void)fprintf(stderr, "   -h          ... Print this message and exit.\n");
  (void)fprintf(stderr, "   -f <nam>    ... use <nam> sa fixname.\n");
  (void)fprintf(stderr, "   -b <val>    ... use <va> as base value.\n");
  (void)fprintf(stderr, "   -p <file>   ... use parameters from <file> .\n");
  (void)fprintf(stderr, "   -g <>   ... .\n");
  (void)fprintf(stderr, "   -s <val>   ... use <val> as seed.\n");
  (void)fprintf(stderr, "   -o <>   ... .\n");
  (void)fprintf(stderr, "   -l <val>    ... use <val> as lambda scale.\n");
  (void)fprintf(stderr, "   -v          ... print version and exit.\n");
  (void)fprintf(stderr, "   -V          ... toggle verbose mode ON.\n");
  (void)fprintf(stderr, "   -x          ... toggle doAnalysis ON.\n");

  exit(exval);
};

void load4(int *x, int a, int b, int c, int d)  
{
 x[0] = a ; 
 x[1] = b ; 
 x[2] = c ; 
 x[3] = d ; 
}

void loadco(double *co, int *fs, int *fsb, int np) 

{
  int a, b, c, d, t  ;  

  a = fs[0] ; 
  b = fs[1] ; 
  c = fs[2] ; 
  d = fs[3] ; 
 
  t = fsb[a*np+c] ;  
  if (t>=0) co[t] += 1 ; 
  t = fsb[b*np+d] ;  
  if (t>=0) co[t] += 1 ; 
  t = fsb[a*np+d] ;  
  if (t>=0) co[t] -= 1 ; 
  t = fsb[b*np+c] ;  
  if (t>=0) co[t] -= 1 ; 


}
void
loadfsindex(int **fsindex, char ***qlist, char **eglist, int numfs, int numeg)  
{
  int k, j, t ; 
  char *sx ; 

  for (k=0; k<numfs; ++k) { 
   for (j=0; j<4; ++j) { 
    sx = qlist[j][k] ; 
    t = indxstring(eglist, numeg, sx) ; 
    if (t<0) fatalx("(loadfsindex) bad pop: %s :: %d %d\n", sx, k, j) ; 
    fsindex[k][j] = t ; 
   }
// printf("zzfs %3d ", k) ; 
// printimat(fsindex[k], 1, 4) ; 
  }
}
