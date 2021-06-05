#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <getpars.h>
#include <globals.h>
#include <xsearch.h>  

#include "admutils.h"
#include "mcio.h"
#include "mcmcpars.h"
#include "regsubs.h"
#include "egsubs.h"
#include "qpsubs.h"

int gsldetails = NO;
double gslprecision = .001;

int  gslsetup (int npars, double *vpars, double initscale) ; 
void gslfree() ; 
double gslopt(double *pars) ; 

double  maxtau = 10.0 ;



/** 
 includes calendar time calculation
*/

int hashsize = 1000 * 1000 ;
ENTRY *hashlist ;   

int npars = 0 ;
int paramindex[1000] ;  
double *vpars ;
int newformat = NO ;

enum actiontype {TASC, TSAMP, TSPEC, TLIFT, TROOT, TADMIX, TNEWFORMAT} ;  

typedef struct { 
 int *cl ; 
 int lev ; 
 double bprob ; 
 int largeliter ;
}  BSTORE ;  

int **clstore ; 
char **clcstore ;
int numclstore ;  
int maxclstore ; 
int xcload = 0  ;
int xchit = 0 ;
double lastkimscore = -100 ;
double minsig = 0 ; // after normalization

BSTORE  **bstlist ; 
int maxbstnum ; 

typedef struct { 
 int count[200] ; 
 int *acount ; 
 int *bcount ;  
 int *ac ; 
 int kode ;  
 int level ;
 int number ;  
 double aprob ; 
 double bprob ;  
 double dprob ;  // external probability (data) 
 double dsig ;   // std. error
 double bmult ;
 int isroot ;  
} CONFIG ;  

typedef struct {  
  enum actiontype type ;    //  1-> 1 parent, 2-> 2 parents, 10-> 0 parents (root) 
  int pop1 ;    // Population index of this vertex (specfied in variables)
  int pop2 ;    // admixture vertex : index of parent;  for lift index of kid
  int pop3 ;    // Index of parent / kid
  int leaf ;  
  int n1 ;
  int n2 ;
  double tau1 ;
  double tau2 ;
  double tmix ;
  double thetaval1 ;
  double thetaval2 ;
  int padmix ;
  int ptau ;
  int ptheta ;
  int pnumber ;
  int mutate ;
} LOGENTRY ; 



#define WVERSION   "811"
// timescore, timegmul added  but not tested
// minsig addes
// iterate.  pass scale to gslsetup 
// seed, ranstart
// tau capped off (setptrans)
//  big bug fixed in mutlike
// default theta  = 1 
// script created if needed 

#define MAXFL  50
#define MAXSTR  512
#define MAXPOPS 100

char *parname = NULL;
char *trashdir = "/var/tmp";
int details = NO;
int qtmode = NO;
int fstdmode = NO;
int hires = NO;
int inbreed = NO;

int isinit = NO;

int seed = 0;
int ranstart = 100 ; 
double ransd = 0.01 ;  

char *graphname = NULL;
char *graphoutname = NULL;
char *graphdotname = NULL;
char *poplistname = NULL;

char *dumpname = NULL;
char *loadname = NULL;
char *rootname;

char *outpop = NULL;
// list of outliers
char *basepop = NULL;
int basenum = -1;
int calchash = NO;
int gsimp = NO;

int bigiter = 1 ;

// outnum used for weights 
// basenum for f3 status  

char **poplist ; 
int *popvertex ; 

double **vmix;
int *lmix, nmix;
int nh2, numeg;
int *ezero = NULL;
double wtmin = .0001;
int xnumeg;
char *newnodename = NULL;
char *edgename = NULL;

double *amut = NULL, bmut, **cmut, *dmut ; // see note.  These are coefficients for mutation along edge


char *outputname = NULL;
char *weightname = NULL;
char *countname = NULL ; 
char *scriptname = NULL ; 

FILE *ofile;
char **eglist;
char **egshort;
char **enames;

void readcommands (int argc, char **argv);
char *mytemp(char *string) ; 

void getvpars(double *vp)   ;  // vp -> log 
void setvpars(double *vp)   ;  // log -> vp 
void putgpars(double *vpars)  ; // log -> graph 
void printpars()    ;

void initlog(LOGENTRY *logpt) ;   
void printlog(LOGENTRY *logpt) ;   
void genrcoeffs(double ****rr, int n)  ;
void print4(char *title, double ****rr, int n)  ;

void genwcoeffs(double **ww, int n)  ;
void printwcoeffs(double **ww, int n)  ;
void setatrans(double **atrans, int n, double tmix)  ;
void setptrans(double **ptrans, int n, double tau)  ;
void setltrans(double **ltrans, double **ptrans, int u, int v) ;
void free4(double *****pvv, int a, int b, int c, int d)  ;
void alloc4(double *****pvv, int a, int b, int c, int d)  ; 
void printv(int *v)  ;
double bprob1(int *cl, int lev)  ; 
void setmut(int n)  ;
double mutlike(double theta, double tau, int u, int v, double **xtrans) ;
double mutlike0(double tau, int u) ;

void setcaltime()  ;
int calctime(int s) ;
double calctscore() ;

int numvertex, numalleles = -1 ;  
char **vnames ;  
// int *acount, *bcount, *tcount ;
int *sampsize, *sampx1,  *tt, specsize ;  
char **sampstr ; 
int *popindex, npop=0 ;  
double *caltime, **calss ;  
int timecheck = -1 ;  
double timesig = 30, timescore , timegscore, timegmul= 0.01 ; 
double ****ttrans, **wtrans ; 



void free_cl()  ;
void alloc_cl()  ;
void set_cl()  ;
void printbcl()  ;

double configprob(int *cfig)  ;
double kimscore() ;
double grscore(double *vp)  ; 
void printconfig(CONFIG *cpt, int mode, char *title) ;
int printtime()   ; 

  CONFIG **configlist ; 
  CONFIG *cpt ; 
  int numconfig =0 , maxconfig = 1000 ;

  LOGENTRY **loglist ;  
  LOGENTRY *logpt ;
  LOGENTRY **loglev ;  
  int numlogs= 0 , maxlog = 100 ;   
  double ****rcoeffs = NULL ;       
  double *****ttable1 ;  // level 1 u v a b 
  double *****ttable2 ;  // level 2 u v a b 
  double *** ttabmix ;  // level  u  a  or level v b 

  int **activelist, **maxlevel, *vact, *vmax, *maxtarget  ; 
  int *isslice, maxlevels ;  
  int *amax, *bmax, *tmax ;

// configurations by level
  int ***clist, *nclist = NULL ; 
  double **wlist ; 

  int maxbst, maxcl  ;
  int largeliter = 0 ;

int
main (int argc, char **argv)
{

  int nedge, t, s, i, j, k, l, z  ;
  int a, b, c, d, u, v, x ; 

  char *psname, *pdfname;
  char sss[MAXSTR];
  char *stmp ; 

  int **lconf; 
  double *lprob ; 
  double *lsig  ; 
  double tau ; 

  int nlconf = 0 ;  
  int level, numlevels ;
  int *tfig, *xpt ;

  double ylbase, y1, y2, yscore, y, yscale ;  
  double *w1, *dprob ;
  int xamax, xbmax ; 
  double *sciter ; 
  int iter ;
  double *tpars, *wpars ;
  double initscore, tscore ;

  readcommands (argc, argv);
  printf ("## kimf version: %s\n", WVERSION);
  printf("graph: %s  counts: %s\n", graphname, countname) ;
  xhcreate(hashsize*2) ; 
  ZALLOC(hashlist, hashsize, ENTRY) ;
  cputime(0) ; 
  if (seed == 0) seed = seednum() ; 
  printf("seed: %d\n", seed) ; 
  SRAND(seed) ;
  printf("maxtau: %9.3f\n", maxtau) ;
    

  if (graphdotname == NULL) graphdotname = strdup("ttt.dot") ;

  numeg = loadgraph (graphname, &eglist);

  numvertex = getnumvertex() ;  
  printf("number of vertices: %d\n", numvertex) ; 
  ZALLOC(configlist, maxconfig, CONFIG *) ; 
  for (i=0; i<maxconfig; ++i) {  
   ZALLOC(configlist[i], 1, CONFIG) ;
   cpt = configlist[i] ; 
   cpt -> acount = cpt -> count ; 
   cpt -> bcount = cpt -> count + numvertex ;
   cpt -> ac = cpt -> count + 2*numvertex ; 
   cpt -> number = i ;  
  }
  ZALLOC(vnames, numvertex, char *) ; 
  t = numlines(countname);
  lconf = initarray_2Dint(t, numvertex, 0) ;

  ZALLOC(loglist, maxlog, LOGENTRY *) ;
  ZALLOC(loglev, maxlog, LOGENTRY *) ;
  for (i=0; i<maxlog; ++i) {  
   ZALLOC(loglist[i], 1, LOGENTRY) ;
   logpt  = loglist[i] ; 
   initlog(logpt) ; 
   loglev[i] = NULL ;
  } 
  
  maxbst = hashsize ;

  ZALLOC(lprob, t, double) ;
  ZALLOC(lsig, t, double) ;

  getvnames(vnames) ;  


  ZALLOC(sampsize, numvertex, int) ;
  ZALLOC(poplist, numvertex, char *) ;
  ZALLOC(popvertex, numvertex, int) ;
  ZALLOC(sampstr, numvertex, char *) ; 
  ZALLOC(sampx1, numvertex, int) ;
  ZALLOC(popindex, numvertex, int) ;  
  ZALLOC(caltime, numvertex, double) ;  
  calss = initarray_2Ddouble(numvertex,3, 0) ;
  npop = 0 ;

  ZALLOC(amax, numvertex, int) ; 
  ZALLOC(bmax, numvertex, int) ; 
  ZALLOC(tmax, numvertex, int) ; 

  for (k=0; k<numvertex; ++k) { 
   sampstr[k] = NULL ; 
  }

  ZALLOC(tt, numvertex, int) ;

  nlconf = loadspec(countname, lconf, vnames, numvertex) ;  // load menu 
  if (scriptname == NULL)  { 
   stmp = mytemp("qqspectra") ; 
   scriptname = strdup(stmp) ; 
   sprintf(sss, "qpreroot -g %s -s %s > /dev/null", graphname, scriptname) ; 
   system(sss) ; 
  }

  printnl() ; 
  printf("scriptname: %s\n", scriptname) ; 
  sprintf(sss, "cat %s", scriptname) ;
  system(sss) ;               
  fflush(stdout) ; 
  
  nlconf = loadspec(scriptname, lconf, vnames, numvertex) ;  // load menu 
  printnl() ; 
  printnl() ; 
  printf("numalleles: %d", numalleles) ;
  printf( " npars: %d ", npars) ;
  printf( " nlconf: %d ", nlconf) ;
  printf( " numlogs: %d ", numlogs) ;
  printnl() ;
  fflush(stdout) ; 

  for (k=0; k<npop; ++k) { 
   z = popindex[k] ;
   printf("pop: %s %s %d\n", poplist[k], sampstr[z], sampsize[z]) ;
  }

 
  ZALLOC(vpars, npars, double) ; 
  setvpars(vpars) ; 
// sets numalleles
  maxlevels = MAX(numlogs + 1, nlconf)  ; 

  activelist = initarray_2Dint(maxlevels, numvertex,  0) ;  
  maxlevel = initarray_2Dint(maxlevels, numvertex,  -1) ;  

  setmut(numalleles) ;
  printpars()  ;

  ZALLOC(isslice, maxlevels, int) ; 
  ZALLOC(maxtarget, numvertex, int) ; 
  ivclear(isslice, NO, maxlevels) ;
   
  level = 0 ; 
  vact = activelist[level] ; 
  vmax = maxlevel[level] ; 
  for (k=0; k<numvertex; ++k) { 
   if (sampsize[k]>0) vact[k] = 1 ;
   if (vact[k] == 1) vmax[k] = level ; 
  }

  xamax = xbmax = 0 ;
  y = 0.0 ; 
  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ;
   y += cpt -> dprob ;
  }
  printf("config scale:  %9.3f\n", y ) ;
  yscale = 1.0/y ; 
  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ;
   cpt -> dprob *= yscale ;
   cpt -> dsig  *= yscale ;
   cpt -> dsig = MAX(cpt -> dsig, minsig) ;
  }
  
  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ;
  
    t = intsum(cpt-> acount, numvertex) ; 
    xamax = MAX(xamax, t) ;
    t = intsum(cpt-> bcount, numvertex) ; 
    xbmax = MAX(xbmax, t) ;

   for (j=0; j<numvertex; ++j) { 

    a = cpt -> acount[j] ; 
    b = cpt -> bcount[j] ; 
    
    a = MIN(a, numalleles) ;
    b = MIN(b, numalleles) ;
    t = a + b ; 
    t = MIN(t, numalleles) ;

    amax[j] = MAX(amax[j], a) ;
    bmax[j] = MAX(bmax[j], b) ;
    tmax[j] = MAX(tmax[j], t) ;

   }
  }
  for (j=0; j<numvertex; ++j) { 
    amax[j] = MIN(amax[j], xamax) ;
    bmax[j] = MIN(bmax[j], xbmax) ;
    tmax[j] = MIN(tmax[j], numalleles) ;
  }

  int a1, b1, t1, a2, b2, t2 ;
  int tact = NO ;

  for (l=0; l<numlogs; ++l) { 
   logpt = loglist[l] ;
   printf("log: %d\n", l) ; printlog(logpt) ; 

   if  (logpt -> type == TROOT) { 
    ++level ;  
    loglev[level] = logpt ; 

    copyiarr(activelist[level-1], activelist[level], numvertex) ; 
    vact = activelist[level] ; 
    s = logpt -> pop1 ; 
    if (vact[s] <= 0) fatalx("rootbug\n") ;

    for (k=0; k<numvertex; ++k) { 
     if (k==s) continue ;
     if (vact[k] != 0) fatalx("rootbug2\n") ;
    }

   }  

   if  (logpt -> type == TADMIX) { 
   ++level ;  
   loglev[level] = logpt ; 
   tact = YES ;
   t = logpt -> pop1 ;  if (t<0) fatalx(" TADMIX badbug\n") ; 
   s = logpt -> pop2 ; 
   maxtarget[s] = level ;

   copyiarr(activelist[level-1], activelist[level], numvertex) ; 
   vact = activelist[level] ; 
   vmax = maxlevel[level] ; 
   a1 = b1 = t1 = 0 ;
   a2 = b2 = t2 = 0 ;

    vact[t] = 0 ; 
    vact[s] = 1 ; 
    a1= amax[t] ; b1=bmax[t] ; t1=tmax[t] ;  

   amax[s] = MAX(amax[s], a1) ;
   bmax[s] = MAX(bmax[s], b1) ;
   tmax[s] = MAX(tmax[s], t1) ;

    s = logpt -> pop3 ; 
    maxtarget[s] = level ;

    vact[t] = 0 ; 
    vact[s] = 1 ; 
    a2= amax[t] ; b2=bmax[t] ; t2=tmax[t] ;  

   amax[s] = MAX(amax[s], a2) ;
   bmax[s] = MAX(bmax[s], b2) ;
   tmax[s] = MAX(tmax[s], t2) ;

  }

   if (logpt -> type == TLIFT) { 
   ++level ;  
   loglev[level] = logpt ; 
   tact = YES ;  
   t = logpt -> pop1 ;  if (t<0) fatalx("badbug\n") ; 
   s = logpt -> pop2 ; 
   maxtarget[t] = level ;
   copyiarr(activelist[level-1], activelist[level], numvertex) ; 
   vact = activelist[level] ; 
   vmax = maxlevel[level] ; 
   a1 = b1 = t1 = 0 ;
   a2 = b2 = t2 = 0 ;
   if (s>=0) { 
    vact[s] = 0 ; 
    vact[t] = 1 ; 
    a1= amax[s] ; b1=bmax[s] ; t1=tmax[s] ;  
   }
   s = logpt -> pop3 ; 
   if (s>=0) { 
    vact[s] = 0 ; 
    vact[t] = 1 ; 
    a2= amax[s] ; b2=bmax[s] ; t2=tmax[s] ;  
   }

   amax[t] = MAX(amax[t], a1+a2) ;
   bmax[t] = MAX(bmax[t], b1+b2) ;
   tmax[t] = MAX(tmax[t], t1+t2) ;

 }

  for (j=0; j<numvertex; ++j) { 
    amax[j] = MIN(amax[j], xamax) ;
    bmax[j] = MIN(bmax[j], xbmax) ;
    tmax[j] = MIN(tmax[j], numalleles) ;
  }

   if (tact == NO) continue ;
   for (k=0; k<numvertex; ++k) { 
    if (vact[k] == 1) vmax[k] = level ; 
   }
  }

  numlevels = level+1 ;
  for (level=0; level < numlevels;  ++level) { 
   t = YES ; 
   vact = activelist[level] ; 
   for (k=0; k<numvertex; ++k) { 
    if (vact[k] == 0) continue ; 
    if (maxtarget[k] > level) t = NO ;
   }
   isslice[level] = t ;
   if (t==YES) { 
    printf("slice on level: %3d\n", level) ;
    printv(activelist[level]) ;  
   }
  }

  printf("amax:\n") ; printv(amax) ; 
  printf("bmax:\n") ; printv(bmax) ; 
  printf("tmax:\n") ; printv(tmax) ; 
  t = numalleles + 1 ; 

  alloc4(&rcoeffs, t, t, t, t) ;
  genrcoeffs(rcoeffs, numalleles)  ;
  wtrans = initarray_2Ddouble(t, t, 0) ;  
  
  ZALLOC(ttable1, maxlevels, double ****) ;
  ZALLOC(ttable2, maxlevels, double ****) ;
  ZALLOC(ttabmix, maxlevels, double **) ;

  for (k=0; k<= numlogs; ++k) { 
   alloc4(&ttable1[k], t, t, t, t) ;
   alloc4(&ttable2[k], t, t, t, t) ;
   ttabmix[k] = initarray_2Ddouble(t, t, 0) ;
  }


  ZALLOC(bstlist, maxbst, BSTORE *) ;  
  for (i=0; i<maxbst; ++i) { 
   ZALLOC(bstlist[i], 1, BSTORE)  ;  
   bstlist[i] -> largeliter = 1 ;  
   bstlist[i] -> bprob = -1 ;
  }

  clstore = initarray_2Dint(hashsize, numvertex*2, 0) ;
  ZALLOC(clcstore, hashsize, char *) ;
  maxclstore = hashsize ; 
  numclstore = 0 ;
  
  alloc_cl() ;

  ++largeliter ;
  

  setvpars(vpars) ; 
  yscore = - grscore(vpars) ;

  printf("init: %9.3f\n", yscore) ;
  printbcl() ;

  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ; 
   printconfig(cpt, 2, "config:") ;
  }
// random start 
  ZALLOC(tpars, npars, double) ;
  ZALLOC(wpars, npars, double) ;
  initscore = yscore ;
  if (bigiter <= 0) ranstart = 0 ;
  for (iter = 1; iter <= ranstart; ++iter) { 
   gaussa(tpars, npars) ;     
   vst(tpars, tpars, ransd, npars) ;
   vvp(tpars, tpars, vpars, npars) ; 
   vabs(wpars, tpars, npars) ; 
   tscore = -grscore(wpars) ; 
   printf("raninit: %3d %9.3f", iter, tscore) ; 
   if (tscore > yscore) {  
    printf(" *** ") ;
    yscore = tscore ; 
    copyarr(wpars, vpars, npars) ;
   }
   printnl() ;
  }

  iter = 0 ; 
  printf("iter: %3d %12.6f\n", iter, yscore) ;
//  printf("bigiter:: %d\n", bigiter) ;  fflush(stdout) ; 

 if (bigiter>0) {
  ZALLOC(sciter, bigiter+1, double) ;
  sciter[0] = yscore ; 
  for (iter = 1; iter <= bigiter; ++iter) { 
   gslsetup(npars, vpars, 0.1) ;
   gslopt(vpars) ;
   getvpars(vpars) ; 
   yscore = -grscore(vpars) ;
   putgpars(vpars) ;
   for (k=0; k<numconfig; ++k) { 
    cpt = configlist[k] ;
    printconfig(cpt, 2, "config:") ;
   }
   printf("iter: %3d %12.6f\n", iter, yscore) ;
   sciter[iter] = yscore ; 
   if (yscore < (sciter[iter-1] + .001)) break ; 
   gslfree() ;
  }
  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ;
   printconfig(cpt, 2, "finalconfig:") ;
  }
  printf("finalscore: %12.6f\n", yscore) ;
 }

  setgtime(caltime) ;  
  dumpgraph (graphoutname);
  dumpdotgraph (graphdotname);

  printpars(vpars) ;
  calctscore() ; 
  printtime() ;
  printf ("## end of kimf time: %9.3f seconds\n", cputime(1));

  return 0 ;

}

void setvpars(double *vp)  
{
   int k, t, tt ; 
   LOGENTRY *logpt ;

  for (k=0; k<npars; ++k) { 
   t = paramindex[k] ; 
   tt = t/1000 ;
   switch (tt) { 

    case 0:  
      logpt = loglist[t] ;  
      vp[k] = logpt -> tau1 ;
      break ; 

    case 1:  
      t -= 1000 ;
      logpt = loglist[t] ;  
      vp[k] = logpt -> tau2 ;
      break ; 

    case 2:  
      t -= 2000 ;
      logpt = loglist[t] ;  
      vp[k] = logpt -> tmix ;
      break ; 

    case 3:  
      t -= 3000 ;
      logpt = loglist[t] ;  
      vp[k] = logpt -> thetaval1;
      break ; 

    case 4:  
      t -= 4000 ;
      logpt = loglist[t] ;  
      vp[k] = logpt -> thetaval2;
      break ; 

     default: fatalx("(setvpars)\n") ;
   }
  }

}

void printbcl() 
// dummy
{

}
void getvpars(double *vp)  
{
   int k, t, tt ; 
   LOGENTRY *logpt ;

  for (k=0; k<npars; ++k) { 
   t = paramindex[k] ; 
   tt = t/1000 ;

   if (tt==0)  {
      logpt = loglist[t] ;  
      logpt -> tau1 = vp[k]  ;
   }

   if (tt==1)  {
      t -= 1000 ;
      logpt = loglist[t] ;  
      logpt -> tau2 = vp[k]  ;
   }

   if (tt==2)  {
      t -= 2000 ;
      logpt = loglist[t] ;  
      logpt -> tmix = vp[k]  ;
   }

   if (tt==3)  {
      t -= 3000 ;
      logpt = loglist[t] ;  
      logpt -> thetaval1 = vp[k]  ;
   }
   if (tt==4)  {
      t -= 4000 ;
      logpt = loglist[t] ;  
      logpt -> thetaval2 = vp[k]  ;
   }
 }
}


int getstype(char *sx) { 

 int t ;

 t = strcmp(sx, "ascertain:") ; 
 if (t==0) return TASC ; 

 t = strcmp(sx, "sampsizes:") ; 
 if (t==0) return TSAMP ; 

 t = strcmp(sx, "spec:") ; 
 if (t==0) return TSPEC ; 

 t = strcmp(sx, "lift:") ; 
 if (t==0) return TLIFT ; 

 t = strcmp(sx, "root:") ; 
 if (t==0) return TROOT ; 

 t = strcmp(sx, "admix:") ; 
 if (t==0) return TADMIX ; 

 t = strcmp(sx, "newformat:") ; 
 if (t==0) { 
   return TNEWFORMAT ; 
 }

 return -1 ;
}

int vnum(char *sx) 
{
  int t ;  

  t = indxstring(vnames, numvertex, sx) ; 
  if (t<0) fatalx("vertex %s not found\n", sx) ;

  return t ; 
}


void printpars()   
{
   int k, t, tt, a, b, c, q ;
   LOGENTRY *logpt ;
   double y, theta ; 

  for (k=0; k<npars; ++k) {
   q = t = paramindex[k] ;
   tt = t/1000 ;
   switch (tt) {

    case 0:
      logpt = loglist[t] ;
      y = logpt -> tau1 ;
      printf("printp: %3d %6d %12.6f\n", k, q, y) ; 
      break ;

    case 1:
      t -= 1000 ;
      logpt = loglist[t] ;
      y = logpt -> tau2 ;
      printf("printp: %3d %6d %12.6f\n", k, q, y) ; 
      break ;

    case 2:
      t -= 2000 ;
      logpt = loglist[t] ;
      y = logpt -> tmix ;
      printf("printp: %3d %6d %12.6f\n", k, q, y) ; 
      break ;

    case 3:
      t -= 3000 ;
      logpt = loglist[t] ;
      theta = logpt -> thetaval1 ;
      printf("printp: %3d %6d %12.6f\n", k, q, theta) ; 
      break ;

    case 4:
      t -= 4000 ;
      logpt = loglist[t] ;
      theta = logpt -> thetaval2 ;
      printf("printp: %3d %6d %12.6f\n", k, q, theta) ; 
      break ;

     default: 
       if (t>4000) fatalx("(printpars)\n") ;
   }
  }
}

void putgpars(double *vpars)   
{
   int k, t, tt, a, b, c, q ;
   LOGENTRY *logpt ;
   double y, theta ; 

  for (k=0; k<npars; ++k) {
   q = t = paramindex[k] ;
   tt = t/1000 ;
   switch (tt) {

    case 0:
      logpt = loglist[t] ;
      y = logpt -> tau1 ;
      theta = logpt -> thetaval1 ;
      a = logpt -> pop1 ; 
      b = logpt -> pop2 ; 
      setepar(a, b, y, theta) ;
      break ;

    case 1:
      t -= 1000 ;
      logpt = loglist[t] ;
      y = logpt -> tau2 ;
      theta = logpt -> thetaval2 ;
      a = logpt -> pop1 ; 
      b = logpt -> pop3 ; 
      setepar(a, b, y, theta) ;
      break ;

    case 2:
      t -= 2000 ;
      logpt = loglist[t] ;
      y = logpt -> tmix ;
      a = logpt -> pop1 ;
      b = logpt -> pop2 ;
      c = logpt -> pop3 ;
      setapar(a, b, c, y) ; 
      break ;

     default: 
       if (tt>4) fatalx("(putgpars)\n, k, q") ;
   }
  }
}

void loadroot(char *buff) 
{
  char *spt[MAXFF], *sx, *s2[20] ;
  int nsplit, nasc, k, z, xd, xt, n ;
  LOGENTRY *logpt ;
  int p1, p2 ; 
  double tau1, tau2 ; 

   nsplit = splitup(buff, spt, MAXFF) ;

   sx = spt[0] ; 
   mkupper(sx) ; 

   z  = vnum(sx) ;

   logpt = loglist[numlogs] ;
   logpt -> type = TROOT ;  
   logpt -> pop1 = z ; 
   logpt -> pop2 = logpt -> pop3 = -1 ;

   ++numlogs ;
}

void loadadmix(char *buff) 
{
  char *spt[MAXFF], *sx, *s2[20] ;
  int nsplit, nasc, k, z, xd, xt, n ;
  LOGENTRY *logpt ;
  int p1, p2 ; 
  double tmix  ; 

   nsplit = splitup(buff, spt, MAXFF) ;

   sx = spt[0] ; 
   mkupper(sx) ; 

   z  = vnum(sx) ;
   getmixinfo(z, &p1, &p2, &tmix) ;

   logpt = loglist[numlogs] ;
   logpt -> type = TADMIX ;  
   logpt -> pop1 = z ; 
   logpt -> pop2 = p1 ; 
   logpt -> pop3 = p2 ; 
   logpt -> tmix = tmix ; 

    paramindex[npars] = numlogs + 2000 ; 
    ++npars ;

   ++numlogs ;
}
void loadlift(char *buff) 
{
  char *spt[MAXFF], *sx, *s2[20] ;
  int nsplit, nasc, k, z, xd, xt, n ;
  LOGENTRY *logpt ;
  int p1, p2 ; 
  double tau1, tau2 ; 
  double th1, th2 ;

   nsplit = splitup(buff, spt, MAXFF) ;

   sx = spt[0] ; 
   mkupper(sx) ; 

   z  = vnum(sx) ;
   getkidinfo(z, &p1, &p2, &tau1, &tau2, &th1, &th2) ;

   logpt = loglist[numlogs] ;
   logpt -> type = TLIFT ;  
   logpt -> pop1 = z ; 
   logpt -> pop2 = p1 ; 
   logpt -> pop3 = p2 ; 
   logpt -> tau1 = tau1 ; 
   logpt -> tau2 = tau2 ; 
   if (th1==0.0) th1 = -1 ;
   if (th2==0.0) th2 = -1 ;
   logpt -> thetaval1 = th1 ; 
   logpt -> thetaval2 = th2 ; 

   if (p1>=0)  { 
    paramindex[npars] = numlogs ; 
    ++npars ;
   }
   if (p2>=0)  { 
    paramindex[npars] = numlogs + 1000 ; 
    ++npars ;
   }

   if ((th1>0) && (p1>=0)){ 
    paramindex[npars] = numlogs + 3000 ; 
    ++npars ;
   }

   if ((th2>0) && (p2>=0)) { 
    paramindex[npars] = numlogs + 4000 ; 
    ++npars ;
   }

   ++numlogs ;
}


/**
void loadasc(char *buff)
{
  char *spt[MAXFF], *sx, *s2[20] ;
  int nsplit, nasc, k, z, xd, xt, n ;

  nasc = nsplit = splitupx(buff, spt, MAXFF, ';') ;
  ivzero(acount, numvertex) ;
  ivzero(tcount, numvertex) ;

  for (k=0; k<nasc; ++k) {
   sx = strdup(spt[k]) ;
   substring(&sx, ":", " ") ;
   n = splitup(sx, s2, 20) ;
   if (n==0) continue ;
   if (n != 3) fatalx("bad loadasc :: %s\n%s %d\n", buff, sx, n ) ;

   mkupper(s2[0]) ; 
   z  = vnum(s2[0]) ; 

   xd = atoi(s2[1]) ;
   xt = atoi(s2[2]) ;
   acount[z] = xd ;
   tcount[z] = xt ;
   freeup(s2, n) ;
   freestring(&sx) ;
  }
  freeup(spt, nsplit) ;
}
*/

int nsize (char *ss) 
// add integer separated by : 
{
  char *spt[MAXFF] ;
  int  k, t, n ;

  n = splitupx(ss, spt, MAXFF, ':') ; 
  t = 0 ; 
  for (k=0; k<n; ++k) { 
    t += atoi(spt[k]) ;
  }
  freeup(spt, n) ; 

  return t ; 
 
}

void loadssize(char *buff)
// sampsizes:  d   :: 1  ;  y   :: 2:20
{
  char *spt[MAXFF], *sx, *sy,  *syy, *s2[10] ;
  int nsplit, ns, k, z, siz,  n, t ;

  if (numvertex == 0) fatalx("bad parse\n") ;
  ns = nsplit = splitupx(buff, spt, MAXFF, ';') ;

  for (k=0; k<ns; ++k) {
   sx = strdup(spt[k]) ;
   substring(&sx, "::", " ") ;
   n = splitup(sx, s2, 10) ;
   if (n==0) continue ;
   if (n != 2) fatalx("bad loadssize :: %s\n", buff) ;

   sy = s2[0] ; 
   syy = strdup(sy) ; 
   
   mkupper(sy)  ; 

   z  = indxstring(vnames, numvertex, sy) ; 
   if (z<0) { 
    z = findlabel(syy) ;
   }
   if (z<0) {
     fatalx("vertex: %s not found\n", syy) ; 
   }

   sx = s2[1] ; 
   siz = nsize(sx) ;  
   if (siz == 0) { 
    freeup(s2, n) ;
    continue ;
   } 
   sampsize[z] = siz ;
   sampstr[z] = strdup(sx) ; 
   poplist[npop] = strdup(syy) ; 
   freestring(&syy) ;

   if (siz==0) continue ;  
   popindex[npop] = z ;
   ++npop ;
  }
  ivsp(sampx1, sampsize, 1, numvertex) ;
  specsize = iprod(sampx1, numvertex) ;

  numalleles = intsum(sampsize, numvertex) ;
  setbino(numalleles+1) ;

  
  freeup(spt, nsplit) ;
 

}

int jcrackit(char *ss, int *a) 
{

  char *spt[MAXFF]  ;
  int nsplit, k ;

  nsplit = splitupx(ss, spt, MAXFF, ':') ; 

  for (k=0; k<nsplit; ++k) { 
   a[k] = atoi(spt[k]) ;
  }

  freeup(spt, nsplit) ;

  return nsplit ;

}

double bcoeff(char *sa, char *sb) 
// product of binomial coeffs;  sa sb are strings : separated which -> int arrays
{
  int *aa, *bb, k, na, nb ; 
  double y, *yy ; 

  ZALLOC(aa, numvertex, int) ; 
  ZALLOC(bb, numvertex, int) ; 
  ZALLOC(yy, numvertex, double) ; 

  na = jcrackit(sa, aa) ; 
  nb = jcrackit(sb, bb) ; 

  if (na != nb) fatalx("specmissmatch: %s %s\n", sa, sb) ;  

  
  for (k=0; k<na; ++k) { 
   yy[k] = bino(aa[k], bb[k]) ;
  }


  y = aprod(yy, na) ;

  if (y<0.1) { 
   printf("bcoeff error\n") ; 
   printmat(yy, 1, na) ;
   fatalx("speczero: %s %s\n", sa, sb) ;  
  }

  free(aa) ; 
  free(bb) ; 
 
  return y ; 

}

void loaddata(char *buff)
{
  char *spt[MAXFF], *sx, *sy  ;
  int nsplit, ns, k, z, t, n, a, b, kk  ;
  double y, ysig, ybprob ; 
  CONFIG *cpt ; 
  int *tt ; 

  if (numvertex == 0) fatalx("bad parse\n") ;
  ZALLOC(tt, numvertex, int) ;  
  nsplit = splitup(buff, spt, MAXFF) ; 

  cpt = configlist[numconfig] ;
  y = atof(spt[nsplit-2]) ; 
  ysig = atof(spt[nsplit-1]) ; 

  ybprob = 1 ;
  for (k=0; k<npop; ++k)  { 
    kk = k ; 
    if (newformat) ++kk ; 
    z = popindex[k] ; 
    sx = spt[kk] ; 
    cpt -> ac[k] = tt[z] = nsize(sx) ;  
    sy = sampstr[z] ;  
    ybprob *= bcoeff(sy, sx) ;
  }

  copyiarr(tt, cpt -> acount, numvertex) ; 
  ivvm(tt, sampsize, tt, numvertex) ; 
  copyiarr(tt, cpt -> bcount, numvertex) ; 
// store configuration  

  cpt -> aprob = 1 ; 
  cpt -> bprob = -1 ; 
  cpt -> dprob = y ; 
  cpt -> dsig = ysig ; 
  cpt -> level = 0 ; 
//set bmult ;
  cpt -> bmult = ybprob ;
  ++numconfig ;  

  freeup(spt, nsplit) ;
  free(tt) ;
 

}

void free4(double *****pvv, int a, int b, int c, int d) 
{
  double ****vv ; 
  int i, j ;

  if (*pvv==NULL) return ;  
  vv = *pvv ; 

  for (i=0; i<a; i++)  {
   for (j=0; j<b; j++)  {
    free2D(&vv[i][j], c) ;
   }
   free(vv[i]) ;
  }
  free(vv) ; 
  *pvv = NULL ;
}
  




void alloc4(double *****pvv, int a, int b, int c, int d) 
{
  double ****vv ; 
  int i, j ;

  free4(pvv, a, b, c, d) ;  
  
  ZALLOC(vv, a, double ***) ; 

  for (i=0; i<a; i++)  {
   ZALLOC(vv[i], b, double **) ; 
   for (j=0; j<b; j++)  {
    vv[i][j] = initarray_2Ddouble(c, d, 0) ;
   }
  }
  *pvv = vv ;
}


void setconf(CONFIG *pt, int *acount, int *bcount, int level) 
{
   int *cc ;; 
   cc = pt -> count ; 
   copyiarr(acount, cc, numvertex) ; 
   cc += numvertex ; 
   copyiarr(bcount, cc, numvertex) ; 
   pt -> level = level ;  

   return ; 
}



int loadspec(char *countname, int **lconf, char **vnames, int numvertex)  
{
#define MAXFF  50

 FILE *fff ;
 char line[MAXSTR+1], c ;
 char buff[MAXSTR] ;
 char *spt[MAXFF], *sx ;
 int nsplit, num=0 ;
 int skipit ;
 int t, k, z, k1, k2, k3, oldmaxterms  ;
 int *w ;
 static int numl = 0 ;  

  openit(countname, &fff, "r") ;
  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ;
   if (nsplit == 0) continue ; // blank line
   sx = spt[0] ;
   if (sx[0] == '#') {
     freeup(spt, nsplit) ;
     continue ;
   }
   catxx(buff, spt+1, nsplit-1) ; // params
   t = strlen(line) ; 

   line[t-1] = CNULL ;
   t = getstype(spt[0]) ;
   printf("%s type: %d\n", line, t) ; 
   if (t<0) fatalx("bad parse:\n%s\n", line) ;

   switch (t) { 

    case TASC: 
     fatalx("ascertain: no longer supported\n") ;

    case TSAMP: 
     loadssize(buff) ; 
     printf("specsize: %d\n", specsize) ;
     break ;

    case TSPEC: 
     loaddata(buff) ; 
     ++numl ; 
     break ; 

    case TLIFT: 
     loadlift(buff) ; 
     break ; 

    case TROOT: 
     loadroot(buff) ; 
     break ; 

    case TADMIX: 
     loadadmix(buff) ; 
     break  ; 

    case TNEWFORMAT: 
     newformat = YES ;
     printf("newformat set!\n")  ;
     break  ; 

    default:  
     printf("???\n") ; 
   }
 }
 fclose(fff) ; 
 return numl ;
}

void initlog(LOGENTRY *logpt)
{
  logpt -> type = -1 ;
  logpt -> pop1 = -1 ;
  logpt -> pop2 = -1 ;
  logpt -> pop3 = -1 ;
  logpt -> n1 = 0 ;
  logpt -> leaf = NO ;
  logpt -> tau1 = logpt -> tau2 = 0 ;
  logpt -> tmix =  -1 ;
  logpt -> ptheta  = -1 ;
  logpt -> thetaval1 = -1 ;
  logpt -> thetaval2 = -1 ;
  logpt -> pnumber = -1 ;
}

char *vn(int k) 
{
 if (k<0) return "?" ; 
 return vnames[k] ;
}

void printlog(LOGENTRY *logpt)
{

  printf("scripttype: %d   %s %s %s",  logpt-> type, vn(logpt -> pop1), vn(logpt -> pop2), vn(logpt -> pop3)) ;
  if (logpt -> type == 1) printf("%9.3f %9.3f", logpt -> tau1, logpt -> tau2) ;
  printnl() ;
  fflush(stdout) ; 

}

void printwcoeffs(double **ww, int n)  
{
   int u, v, a, b, t ; 

   printf("wcoeffs:\n") ;
   for (t=2; t<=n; ++t)  {  
    for (u=1; u<t; ++u)  {  
     v = t-u ; 
       printf(" %d %d ", u, v) ;
       printf("%12.6f", ww[u][v]) ;
       printnl() ;  
    }
   printnl() ;
  }
}
void print4(char *title, double ****rr, int n)  
{
   int u, v, a, b, t ; 
   double y ; 

   printf("coeffs4: %s\n", title) ;
   for (t=0; t<=n; ++t)  {  
    for (u=0; u<=t; ++u)  {  
     v = t-u ; 
     for (a=0; a<=u; ++a)  {  
      for (b=0; b<=v; ++b)  {  
       printf(" %d %d %d %d ", u, v, a, b) ;
       y =  rr[u][v][a][b] ;
       printf("%12.6f", y) ; 
       if ((a==1) && (b==1))  printf(" %12.6f ***", bino(u+v, v) * y) ; 
       printnl() ;  
      }
     }
    printnl() ;
   }
   printnl() ;
  }
}

void setltrans(double **ltrans, double **ptrans, int u, int v) 
// genrcoeff must have been called    
{
  int t = u+v, tt, a, b    ; 
  clear2D(&ltrans, u+1, v+1, 0) ;  

  for (a=0; a<=u; ++a)  { 
   for (b=0; b<=v; ++b)  { 
    tt = a + b ; 
    ltrans[a][b] = rcoeffs[u][v][a][b]*ptrans[t][tt] ;
   }
  }
}

void setatrans(double **atrans, int n, double tmix)     
// set up coeffs of [(1-tmix) + tmix * x] ^ k 
{
  double ya, yb, *ww  ; 
  int k ; 

  ya = tmix; 
  yb = 1.0 - ya ;

  ZALLOC(ww, n+1, double) ;
  clear2D(&atrans, n+1, n+1, 0) ;
  atrans[0][0] = 1 ; 

  for (k=1; k<=n; ++k) {  
   vst(atrans[k] + 1, atrans[k-1], ya, k) ; // multiply by a x 
   vst(ww, atrans[k-1], yb, n) ;   // multiply by b 
   vvp(atrans[k], atrans[k], ww, k) ;  
   if (verbose) printf("pacheck: %d %12.6f\n", k, asum(atrans[k], n+1)) ; 
  }  

  free(ww) ;
}

void setptrans(double **ptrans, int n, double tau)     
{
  int  i, j, k  ; 
  double *dd,  yt, yb, y, ttau ; 
  
 //  ttau = tau ; 
  ttau = MIN(tau, maxtau) ;

  ZALLOC(dd, n+1, double) ;
  clear2D(&ptrans, n+1, n+1, 0) ;
  ptrans[0][0] = 1 ; 

  for (k=1; k<= n; ++k) {
   y = dd[k] = (double) (k*(k-1)) / 2.0 ;
   ptrans[k][k]  = exp(-y*ttau) ;
  }
   
   for (i=1; i<=n; ++i) { 
    for (j=i-1; j>=1; --j)  { 
     if (tau<=0.0) break ;
     yt = dd[i]*ptrans[i-1][j] - dd[j+1]*ptrans[i][j+1] ; 
     yb = dd[i] - dd[j] ;  
     ptrans[i][j] = yt/yb ;  
    }
    if (verbose) printf("ptcheck: %3d  %9.3f %3d %12.6f\n", n, tau, i, asum(ptrans[i], n+1) ) ;
   }
  free(dd) ; 
}

void genwcoeffs(double **ww, int n)  
// recursion from "considering the first coalescent" .  At root  and polymorphic
{
   int u, v,  t ; 
   double c1, c2, y1, yt, yb ; 

   ww[1][1] = 1 ;  
  
   for (t=3; t<=n; ++t) {  
    u = 1 ;  // special case.  Mutation can occur  
     v = t - u ; 

     y1 = 2.0 / (double) (v*(v+1))  ; 
     c1 = (double) (v-1) / (double) (v+1) ;  
     ww[u][v] = y1 + c1*ww[u][v-1] ; 

      for (u=2; u<t; ++u) {  
       v = t - u ; 
       yb =  (double) (t*(t-1))  ; 
       c1  = (double) (u*(u-1)) / yb ; 
       c2  = (double) (v*(v-1)) / yb ; 
       ww[u][v] = c1*ww[u-1][v] + c2*ww[u][v-1] ;  
      }
   } 
}


void genrcoeffs(double ****rr, int n)  
// recursion from "considering the first coalescent" 
{
   int u, v, a, b, t ; 
   double c1, c2, yb ; 
   rr[0][0][0][0] = 1 ; 
   for (u=1; u<=n; ++u) {  
    for (a=1; a<=u; ++a) {  
     rr[u][0][a][0] = 1 ;  
     rr[0][u][0][a] = 1 ;  
    }
   }
  
   for (t=1; t<=n; ++t) {  
    for (u=1; u<t; ++u) {  
     v = t - u ; 
     rr[u][v][u][v] = 1 ; 

     yb =  (double) (t*(t-1))  ; 
     c1  = (double) (u*(u-1)) / yb ; 
     c2  = (double) (v*(v-1)) / yb ; 

      for (a=1; a<=u; ++a) {  
       for (b=1; b<=v; ++b) {  
        if (t==(a+b)) continue ; 
         rr[u][v][a][b] = c1*rr[u-1][v][a][b] + c2*rr[u][v-1][a][b] ; 
       }
      }
    }
   } 
}

void free_cl() 
// call when parameters are reset 
{
 int level ; 

 if (nclist == NULL) return ;
  
 free(nclist) ; 
 free2D(&wlist, maxlevels+1) ;
  for (level=0; level <= maxlevels; ++level) { 
   free2Dint(&clist[level], maxcl) ;
  }

  free(clist) ;
  nclist = NULL ;


}
void alloc_cl() 
{
  int level, i ; 
  int a, b, u, v ;
  int s, t, maxa, maxb, maxu, maxv, maxt ;
  double **ttr ; 
  double **** ttab ; 
  LOGENTRY *logpt ; 
  
  ++largeliter ; 
  if (verbose) printf("largeliter: %d  xcload xchit: %d %d %x score: %12.6f\n", largeliter, xcload, xchit, topheap(), lastkimscore) ;
  fflush(stdout) ;

  xcload = xchit = 0 ;

  free_cl() ;

  ZALLOC(clist, maxlevels+1, int **) ; 
  for (level=0; level <= maxlevels; ++level) { 
   clist[level]  = initarray_2Dint(maxcl, numvertex*2, 0) ; 
  }
  wlist = initarray_2Ddouble(maxlevels+1, maxcl, 0) ;
  ZALLOC(nclist, maxlevels+1, int) ;      

  ttr = initarray_2Ddouble(numalleles+1, numalleles+1, 0) ;  
  for (level=0; level<=maxlevels+1; ++level) {  
   if (loglev[level] == NULL) continue ; 
   logpt = loglev[level] ; 
   if (logpt -> type == TADMIX) { 
    s = logpt -> pop1 ;
    maxt = tmax[s]  ; 
    setatrans(ttabmix[level], maxt, logpt -> tmix) ; 
    continue ; 
   }
   t = logpt -> pop1 ; 
   maxa = amax[t] ; 
   maxb = bmax[t] ; 
   s = logpt -> pop2 ;
   if (s>=0) {  
     maxu = amax[s] ; 
     maxv = bmax[s] ; 
     maxt = tmax[s] ;  
     setptrans(ttr, maxt, logpt -> tau1) ;
     ttab = ttable1[level] ; 
     for (u=0; u<=maxu; u++)  { 
      for (v=0; v<=maxv; v++)  { 
       if ((u+v) > maxt) continue ; 
        setltrans(ttab[u][v], ttr, u, v) ; 
       }
     }
   }

   s = logpt -> pop3 ;
   if (s>=0) {  
     maxu = amax[s] ; 
     maxv = bmax[s] ; 
     maxt = tmax[s] ;  
     setptrans(ttr, maxt, logpt -> tau2) ;
     ttab = ttable2[level] ; 
     for (u=0; u<=maxu; u++)  { 
      for (v=0; v<=maxv; v++)  { 
       if ((u+v) > maxt) continue ; 
        setltrans(ttab[u][v], ttr, u, v) ; 
       }
     }
   }
  }
}

void set_cl() 
{
  int level ; 
  int **cl ;

  for (level=0; level <= maxlevels; ++level) { 
   cl = clist[level] ;
   iclear2D(&cl, maxcl, numvertex, 0) ;
   vzero(wlist[level], maxcl) ;
  }
  ivzero(nclist, maxlevels+1) ;

}
void liftit(int lev, int u, int v, double **ans, int mode)  
{
 double **source ; 
 int t ; 
  
 t = MAX(u, v) ;  

 if (mode == 3) { 
  source = ttabmix[lev] ; 
  copyarr2D(source, ans, t+1, t+1) ;  
  return ; 
 }

 if (mode == 1) source = ttable1[lev][u][v] ; 
 if (mode == 2) source = ttable2[lev][u][v] ; 

 copyarr2D(source, ans, u+1, v+1) ; 
}

void putab (int *cl, int z, int pa, int pb) 
// store pa, pb in config     
{
 if (z<0) return ; 
 cl[z]   = pa ; 
 cl[z+numvertex] = pb ;
 return ;
}

void getab (int *cl, int z, int *pa, int *pb) 
//extract a, b from config  
{
 if (z<0) { 
  *pa = *pb = 0 ; 
 }
 else {
  *pa = cl[z] ; 
  *pb = cl[z+numvertex] ;
 }
 return ; 
}
double top_prob(int a, int b) 
{
  double yfac = 1.0, y, y1, ans ; 

    yfac = 1.0 ;
    y = bino(a+b, a) ;  
    y1 = 2.0/ (double) a ;
    ans = yfac*y1/y ;  
    return ans ;
}

double bprob1_root(int *cl, int lev) 
{
  LOGENTRY *logpt ;
  int p1, a, b ;  


//  yfac = pow(2, numalleles) ;
  logpt = loglev[lev] ;
 
   p1 = logpt -> pop1 ;
   getab(cl, p1, &a, &b) ; 

    return top_prob(a, b) ;

}


/**
typedef struct { 
 int *cl ; 
 double bprob ; 
}  BSTORE ;  

int **clstore ; 
int numclstore ;  
int maxclstore ; 

BSTORE  **bstlist ; 

*/


int hashiarr(char *str, int *a, int n) 
{
  char *sx ; 
  static char **spt  ; 
  char ss[10] ; 
  int k, t ; 
  static long ncall = 0 ;

  ++ncall ;  
  if (ncall==1) { 
   ZALLOC(spt, 100, char *) ;  
  }
  sx = str ;
  for (k=0; k<n; ++k) { 
   sprintf(ss, "%d", a[k]) ; 
   spt[k] = strdup(ss) ;
  }
  catxc(str, spt, n, ':') ;
  freeup(spt, n) ;
  return(strlen(str)) ;
}

int seekbst(int lev, int *cl) 
{
 ENTRY item1, *ept, *fpt ; 
 static void *base ;
 int x, i, *act ; 

 BSTORE *bstpt ;
 int *tcl, k, t, ll ; 
 int w1[100], *jpt ; 
 char s1[200] ; 
 static int ncall = 0 ;

 ++ncall ;
 base = 0  ; 
 x = 0 ;
 act = activelist[lev] ;  
 for (i=0; i<2*numvertex; ++i) { 
 // if (act[i] == 0) continue ;  
  w1[x] = cl[i] ; 
  ++x ;
 }
 w1[x] = lev ; 
 ++x ; 

 hashiarr(s1, w1, x) ;  
  
 ept = &item1 ; 
 ept -> key = s1 ; 
 ept -> data   = NULL ;
 
 fpt = xhsearch(*ept, FIND) ;

 if (fpt == NULL) { 
  copyiarr(cl, clstore[numclstore], 2*numvertex) ;
  ept -> key = clcstore[numclstore] = strdup(s1) ;  
  ept -> data = base +  numclstore ;
  ++numclstore ;
  ++xcload ;  
  fpt = xhsearch(item1, ENTER) ;

  t =  fpt -> data - base ; 
  bstpt = bstlist[t] ; 
  bstpt -> largeliter = -1; 
  bstpt -> bprob = -2 ;
 }

 t =  fpt -> data - base ; 
 return t ; 

}

double bprob1_admix(int *cl, int lev) 
{

 int lp, p1, p2, p3, t, x, k, u, v, tmax ; 
 int olda1, oldb1, olda2, oldb2 ; 
 int maxa, maxb ; 
 int a1, a2, b1, b2 ; 
 int xa1, xa2, xb1, xb2 ; 
 double y1, y2, cprob, ytot, ans ; 
 double z1, z2 ; 
 LOGENTRY *logpt ; 
 double **x1trans, **ttrans ; 
 double *xprobs ; 
 int **xlist, nx ; 
 int *wcl ; 
 BSTORE *bstpt ; 

 verbose = NO ;
 lp = lev + 1 ; 
 logpt = loglev[lp] ;

// admix 
 p1 = logpt -> pop1 ; // source  
 p2 = logpt -> pop2 ; 
 p3 = logpt -> pop3 ; 



 ZALLOC(wcl, numvertex*2, int) ;
 copyiarr(cl, wcl, numvertex*2) ;  

 if (verbose) { 
  printf("zzadmix:\n") ; 
  printimat(cl, 2, numvertex) ;  
 }

 
 getab(cl, p1, &xa1, &xb1) ;
 getab(cl, p2, &olda1, &oldb1) ;
 getab(cl, p3, &olda2, &oldb2) ;

 u = maxa = xa1  ; 
 v = maxb = xb1  ; 

  tmax = t = MAX(maxa, maxb) + 1 ;

  x1trans = initarray_2Ddouble(t, t, 0) ; 
  ttrans = initarray_2Ddouble(maxa+1, maxb+1, 0) ; 

  putab(wcl, p1, 0, 0) ; 

  liftit(lp, cl[p1],  cl[p1+numvertex], x1trans, 3) ; 

 for (a1=0; a1<=u; ++a1)  { 
  y1 = x1trans[u][a1] ;  
  for (b1=0; b1<=v; ++b1)  { 
   y2 = x1trans[v][b1] ; 
   ttrans[a1][b1] = y1*y2 ; 
   if (verbose == NO) continue ;
   printf("zzqtrans: %d %d %d %d ", u, v, a1, b1) ;  
   printf("%12.6f ", y1) ;
   printf("%12.6f ", y2) ;
   printf("%12.6f ", y1*y2) ;
   printnl() ;
  }}

  ytot = 0 ;
  for (a1=0; a1<=u; ++a1)  { 
   for (b1=0; b1<=v; ++b1)  { 
    z1 = y1 = ttrans[a1][b1] ; 
    if (y1<=0) continue ;  
    if (verbose) printf("zza %d  %3d %3d %15.9f %d %d\n", lev, a1, b1, y1, xa1, xb1)  ; 
    putab(wcl, p2, a1+olda1, b1+oldb1) ;
    a2 = u-a1  ; 
    b2 = v-b1  ;
    putab(wcl, p3, a2+olda2, b2+oldb2) ;
    y2 = bprob1(wcl, lev+1) ;
    if (verbose) {
     printf("zz3 %15.9f %12.6f %12.6f %12.6f %12.6f\n", y2, y1,  y1*y2) ;  
     printimat(wcl, 2, numvertex) ;
    }
    ytot += y1*y2 ;
   }}


 free2D(&x1trans, tmax) ; 
 free2D(&ttrans, maxa+1) ; 

 verbose = NO ;

 free(wcl) ;
 ans = ytot ;  
 return ans ;

// bprob1 calls addbst

}

int numnonz(int *a, int n) 
{
 int k, t ;

 t = 0 ;
 for (k=0; k<n; ++k) { 
  if (a[k] != 0) ++t ;  
 } 

 return t ; 

}

void addbst(int bstindex,  double ans) 
{
 BSTORE *bstpt ;

 bstpt = bstlist[bstindex] ; 
 bstpt -> bprob = ans ; 
 bstpt -> largeliter = largeliter ;

}

double bprob1(int *cl, int lev) 
{

 int lp, p1, p2, p3, t, x, k, olda, oldb; 
 int maxa, maxb ; 
 int a1, a2, b1, b2, u, v ; 
 int xa1, xa2, xb1, xb2 ; 
 double y1, y2, cprob, ytot, ans ; 
 LOGENTRY *logpt ; 
 double **x1trans, **x2trans, **ttrans, **xtrans ; 
 double *xprobs ; 
 int **xlist, nx ; 
 int *wcl ; 
 int bstindex ; 
 BSTORE *bstpt ; 
 double checkprob = -1  ;  
 double theta, tau ;

 lp = lev + 1 ; 
 logpt = loglev[lp] ;

 bstindex = k = seekbst(lev, cl) ; 

 if (k >= 0) { 
  bstpt = bstlist[k] ;
  if (bstpt -> largeliter == largeliter) {
   checkprob = bstpt -> bprob ;
   ++xchit ; 
   return bstpt -> bprob ;
  }
 }

 if (k<0) fatalx("seekbst failed\n") ;



 if (logpt -> type == TROOT) {
   ans =  bprob1_root(cl, lev) ;
   addbst(bstindex,  ans) ;
   return ans ; 
 }

 if (logpt -> type == TADMIX) {
   ans =  bprob1_admix(cl, lev) ;
   addbst(bstindex,  ans) ;
   return ans ; 
 }
// lift 
 p1 = logpt -> pop1 ; // target for lift 
 p2 = logpt -> pop2 ; 
 p3 = logpt -> pop3 ; 



 ZALLOC(wcl, numvertex*2, int) ;
 copyiarr(cl, wcl, numvertex*2) ;  

 
 getab(cl, p1, &olda, &oldb) ;
 getab(cl, p2, &xa1, &xb1) ;
 getab(cl, p3, &xa2, &xb2) ;

 maxa = xa1 + xa2 + olda ; 
 maxb = xb1 + xb2 + oldb ; 


 x1trans = initarray_2Ddouble(maxa+1, maxb+1, 0) ; 
 x2trans = initarray_2Ddouble(maxa+1, maxb+1, 0) ; 
 ttrans = initarray_2Ddouble(maxa+1, maxb+1, 0) ; 


 putab(wcl, p2, 0, 0) ; 
 putab(wcl, p3, 0, 0) ; 

 x1trans[0][0] = x2trans[0][0] = 1 ;  
 if (p2>=0)  liftit(lp, cl[p2], cl[p2+numvertex], x1trans, 1) ; 
 if (p3>=0)  liftit(lp, cl[p3], cl[p3+numvertex], x2trans, 2) ; 


 for (a1=0; a1<=xa1; ++a1)  { 
  for (b1=0; b1<=xb1; ++b1)  { 
   y1 = x1trans[a1][b1] ; 
   if (verbose) printf("zzqq1 %d  %3d %3d %15.9f %d %d\n", lev, a1, b1, y1, xa1, xb1) ; 
   if (y1<=0) continue ;
 for (a2=0; a2<=xa2; ++a2)  { 
  for (b2=0; b2<=xb2; ++b2)  { 
   y2 = x2trans[a2][b2] ; 
   if (verbose)  printf("zzqq2 %d %d %12.6f %d %d\n", a2, b2, y2, xa2, xb2)  ; 
   if (y2<=0) continue ;
   ttrans[a1+a2+olda][b1+b2+oldb] += y1*y2 ; 
 }}}}

  ytot = 0 ;
  for (a1=0; a1<=maxa; ++a1)  { 
   for (b1=0; b1<=maxb; ++b1)  { 
    y1 = ttrans[a1][b1] ; 
    if (y1<=0) continue ;  
    putab(wcl, p1, a1, b1) ;
    y2 = bprob1(wcl, lev+1) ;
    ytot += y1*y2 ;
    if (verbose) printf("zztran: %3d %3d %3d  %12.6f %12.6f %15.9f %d %d %d\n", lev, a1, b1, y1, y2, y1*y2, p1, p2, p3) ; 
   }}

 ans = ytot ;  

 t = numnonz(cl, numvertex) ;
// number of pops with derived count non zero
 xtrans = NULL ;
 if ((t==1) && (xa1>0))  {  
   u = xa1 ; v = xb1 ; xtrans = x1trans, tau = logpt -> tau1 ;
   theta = logpt -> thetaval1 ;
   if (theta<0) theta = 1.0 ;
 }
  if ((t==1) && (xa2>0)) { 
   u = xa2 ; v = xb2 ; xtrans = x2trans ; tau = logpt -> tau2 ;
   theta = logpt -> thetaval2 ;
   if (theta<0) theta = 1.0 ;
  }
  if (xtrans != NULL) ans += mutlike(theta, tau, u, v, xtrans) ;

/** 
 We just add probabilities here as the model has an eps in mutation probability 
 So prob of mutation on edge is O(eps) and remaining prob is multiplied by 1 - eps. 
 Thus to order eps in total prob we just add 
*/


 free2D(&x1trans, maxa+1) ; 
 free2D(&x2trans, maxa+1) ; 
 free2D(&ttrans, maxa+1) ; 


 free(wcl) ;
 if (checkprob >= 0) { 
  printf("zzcheck %12.6f %12.6f", ans, checkprob) ;
  if (ans != checkprob) {
    printf (" ***") ;
    printnl() ; 
    printf ("lev: index: %d %d\n", lev, bstindex) ; 
    printimat(cl, 2, numvertex) ; 
    printnl() ;
    printimat(clstore[bstindex], 2, numvertex) ; 
    printnl() ; 
  }
  printnl() ;
 }

 addbst(bstindex, ans) ;
 return ans ; 

}



double configprob(int *cfig) 
{
  int level ; 

  set_cl() ; // check 

  return bprob1(cfig, 0) ;
  
}
void printv(int *v) 
{
 int k ; 
 for (k=0; k<numvertex; ++k)  {  
  printf(" %3s", vnames[k]) ;
 }
 printnl() ; 
 for (k=0; k<numvertex; ++k)  {  
  printf(" %3d", v[k]) ;
 }
 printnl() ; 
}



void
readcommands (int argc, char **argv)
{
  int i, haploid = 0;
  phandle *ph;
  char str[5000];
  char *tempname;
  int n;

  while ((i = getopt (argc, argv, "s:p:r:c:g:o:d:n:m::vVx")) != -1) {

    switch (i) {

    case 'p':
      parname = strdup (optarg);
      break;

    case 's':
      scriptname = strdup (optarg);
      break;

    case 'r':
      rootname = strdup (optarg);
      break;

    case 'g':
      graphname = strdup (optarg);
      break;

    case 'c':
      countname = strdup (optarg);
      break;

    case 'o':
      graphoutname = strdup (optarg);
      break;

    case 'd':
      graphdotname = strdup (optarg);
      break;

    case 'n':
      bigiter = atoi (optarg);
      break;

    case 'm':
      hashsize = atoi (optarg);
      break;

    case 'v':
      printf ("version: %s\n", WVERSION);
      break;

    case 'V':
      verbose = YES;
      break;

    case 'x':
      bigiter = 0 ;
      break;

    case '?':
      printf ("Usage: bad params.... \n");
      fatalx ("bad params\n");
    }
  }


  if (parname == NULL) {
    return;
  }

  printf ("parameter file: %s\n", parname);
  ph = openpars (parname);
  dostrsub (ph);

  getstring (ph, "graphname:", &graphname);
  getstring (ph, "graphoutname:", &graphoutname);
  getstring (ph, "graphdotname:", &graphdotname);
  getstring (ph, "output:", &outputname);

  getstring (ph, "dumpname:", &dumpname);
  getstring (ph, "loadname:", &loadname);
  getint(ph, "gsldetails:", &gsldetails) ;
  getint(ph, "hires:", &hires) ;
  getdbl(ph, "timegmul:", &timegmul) ;
  getdbl(ph, "timesig:", &timesig) ;
  getdbl(ph, "minsig:", &minsig) ;

  getint(ph, "seed:", &seed) ;
  getint(ph, "ranstart:", &ranstart) ;
  getdbl(ph, "ransd:", &ransd) ;
  getdbl(ph, "maxtau:", &maxtau) ;

  writepars (ph);

}
double grscore(double *vp) 
{
  double y ; 
  double bigtau = 0  ;  
  int k ;  

  getvpars(vp) ; 
  
  for (k=0; k<npars; ++k) { 
   if (vp[k] > maxtau) bigtau += (vp[k] - bigtau) ; 
  }


  alloc_cl() ; // free_cl called here

  y = kimscore() ;
  y -= 0.01 * bigtau ; 
  return -y ; 

}
double kimscore () 
{
  double *w1, *dprob ; 
  int *tfig, *xpt ;
  double y1, y2, yscore, ylbase ; 
  CONFIG *cpt ;
  int k ; 
  double ytim ; 

  ZALLOC(w1, numconfig, double) ;
  ZALLOC(dprob, numconfig, double) ;
  ZALLOC(tfig, 2*numvertex, int) ;
  if (timecheck != 0) setcaltime() ; 

  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ;
   xpt = tfig ; 
   copyiarr(cpt -> acount, xpt, numvertex) ;
   xpt += numvertex ; 
   copyiarr(cpt -> bcount, xpt, numvertex) ;

   y1 = configprob(tfig)  ; 
   if (y1 == 0.0) {
    printconfig(cpt, 2, "bad configuration") ; 
    fatalx("bad score\n") ;
   }
   y2 = y1 * cpt -> bmult ;
   y2 += 1.0e-10 ; 
   w1[k] = y2 ; 
   dprob[k] = cpt -> dprob ; 
  }

  y1 = bal1(w1, numconfig) ;  
  bal1(dprob, numconfig) ;  

  for (k=0; k<numconfig; ++k) { 
   cpt = configlist[k] ;
   cpt -> aprob = w1[k] ;
  }

  vst(w1, w1, numconfig, numconfig) ; // score against random . ie all configs equiprobable 
  yscore = 100 * vldot(dprob, w1, numconfig) ;

  free(dprob) ; 
  free(w1) ; 
  free(tfig) ;

  ytim= calctscore() ; 
  y1 = yscore ; // old score 
  yscore -= ytim ; 
  lastkimscore = yscore ;     
  if (verbose) printf("zzsc %12.6f %12.6f %12.6f\n" , y1, ytim, yscore) ;
  return yscore ;

}
void setmut(int n) 
{
 int k, j ; 
 double y ; 

 ZALLOC(amut, n+1, double) ;
 bmut = 1 ; 
 cmut = initarray_2Ddouble(n+1, n+1, 0) ; 
 ZALLOC(dmut, n+1, double) ; 
 for (k=2; k<=n; ++k) { 
  dmut[k] = (double) (k*(k-1)) / 2.0 ; 
 }
 for (k=2; k<=n; ++k) { 
  amut[k] = amut[k-1]  - (1.0 / dmut[k]) ;  
 }
 for (k=2; k<=n; ++k) { 
  for (j=2; j<k; ++j) { 
   y = dmut[k]-dmut[j] ; 
   cmut[k][j] = (dmut[k]/y) *cmut[k-1][j] ; 
  }
  y = amut[k] + asum(cmut[k], k) ; 
  cmut[k][k] = -y ; 
 }
/**
 for (k=2; k<=n; ++k) { 
  printf("mutdiag: %15.9f %15.6f\n", cmut[k][k], 1/cmut[k][k]) ;
 }
*/

}

double mutlike(double theta, double tau, int u, int v, double **xtrans) 
{
  double y0, y1, y2, trprob, y ;
  int a, b ;

  if (theta<=0) return 0 ;
  if (amut == NULL) setmut(numalleles) ;  

  if (v==0) return  theta*mutlike0(tau, u) ;

  y1 = top_prob(u, v) ;
  y2 = 0 ;

/**
  for (a=1; a<=u; ++a)  { 
   for (b=1; b<=v; ++b)  { 
    trprob = xtrans[a][b] ; 
    y = trprob*top_prob(a, b) ;
    y2 += y ;
    printf("zzq %d %d %d %d %12.6f %12.6f\n", a, b, u, v, trprob, y) ; 
  }} 
  y0 = y1 - y2 ;
*/

  y1 = (1.0 - xtrans[u][v]) * top_prob(u,v) ;  
  for (a=1; a<=u; ++a)  { 
   for (b=1; b<=v; ++b)  { 
    if ((a==u) && (b==v)) continue ; 
    trprob = xtrans[a][b] ; 
    y = trprob*top_prob(a, b) ;
    y2 += y ;
//  printf("zzq %d %d %d %d %12.6f %12.6f\n", a, b, u, v, trprob, y) ; 
  }} 
  y0 = y1 - y2 ;

//printf("zzmutlike %d %d %12.6f %12.6f %12.6f\n", u, v, y1, y2, y0) ;  

 return y0 * theta ;

}


double mutlike0(double tau, int u)

{ 

  double *ww, y  ;  

  if (u==0) return 0 ;
  y = amut[u] + bmut*tau ;
  ZALLOC(ww, u+1, double) ; 

  vst(ww, dmut, -tau, u+1) ; 
  vexp(ww, ww, u+1) ; 
  ww[0] = 0 ;

  y += vdot(ww, cmut[u], u+1) ;

 free(ww) ;
 return y ;

} 

double mkcanon (double *vp)
// make up canonical return penalty
{
  double penalty = 0;
  int k, s, t, tt;
  double xx, y, yy;

  for (k=0; k<npars; ++k) { 
   t = paramindex[k] ; 
   tt = t/1000 ;

   s = 0 ;

   if (tt==2)  {
     s = 1 ;  
   }


   yy = xx = vp[k] ; 
   if ((xx>=0) && (xx<=1)) continue ;
    if (xx<0) yy = -yy ; 
    if ((s==1) & (yy>1)) {  // bad mixing coeff
     yy = 2.0*modf(0.5*yy, &y) ;  
     if (yy>1) yy = 2-yy ; 
   }
     
  y  = yy - xx;
  penalty += y * y;
  vp[k] = xx = yy;
 }
 
  return penalty;
}

void printconfig(CONFIG *cpt, int mode, char *title) 
{
   double y1, y2 ; 

   printf("%15s %3d npop: %d", title, cpt -> number, npop) ; 
   if (mode == 2) { 
    printimatx(cpt -> ac, 1, npop) ;
   }

   if (mode >= 1) { 
    printf("  %12.6f", cpt -> dprob) ;
    printf(" %12.6f", cpt -> aprob) ;
    y1 = cpt -> aprob - cpt -> dprob ;  
    y2 = cpt -> dsig ;
    printf(" %9.3f", y1/y2) ;   // Z score 
    printnl() ; 
    return ; 
   }

   printnl() ; 
   printimat(cpt -> acount, 1, numvertex) ; 
   printimat(cpt -> bcount, 1, numvertex) ; 
   printf("bmult: %12.0f",  cpt -> bmult) ;
   printnl() ;

}

// caltime code

void inittime() 
{
  double y ; 
  int k, z ;


  vclear(caltime, -999., numvertex) ; 
  clear2D(&calss, numvertex, 3, 0) ; 

  for (k=0; k < npop; ++k)  { 
   z = popindex[k] ;
   y = caltime[z] = 0 ;  // will generalise
   calss[z][0] = 1 ; 
   calss[z][1] = y ; 
   calss[z][2] = y * y  ; 
  }

}

double  callength(double tau, double theta)  
// time in generations,  10K nominal eff. pop. size
{
  double y  ; 
  if (tau<0) return -999 ;  
  if (theta < 0)  return -888 ; 
 
  y = tau*theta*2*10000 ;  // effect pop 10000 time in generations

  return y  ; 

}

int printtime()  
{
 int k ; 

 if (timecheck==0) return 0 ; 
  for (k=0; k<numvertex; ++k) { 
   printf("vtime: %20s %9.0f", vnames[k], caltime[k]) ; 
   printnl() ;
  }
  printf("timescore: %9.3f  timegscore: %9.3f\n", timescore, timegscore) ;
  return timecheck ;
}
int settime1(int target, int source, double len) 
// set time array (just calss)  
{
  double ys, yt ;

  if (target < 0) return -99 ; 
  if (source < 0) return -9 ; 

  calctime(source) ; 
  ys = caltime[source] ; 
  if (ys<0) return -8 ; 
  if (len <0) return -7 ; 

  yt = ys + len ; 
  calss[target][0] += 1 ;
  calss[target][1] += yt ;
  calss[target][2] += yt*yt ;

  return 1  ; 

}
int calctime(int s) 
// set caltime from calss (just mean)  
{
  double y0, y1 ;

  if (s<0) return -2 ;
  y0 = calss[s][0] ;  
  if (y0<=0.1) return -1 ; 

  y1 = calss[s][1] ;  
  caltime[s] = y1/y0 ; 

  return 1 ;  

}
// **caltime 

void setcaltime() 
{
  int t, l, a, b, s, tt, t2, x ; 
  double y ; 
  LOGENTRY *logpt ;

  inittime() ;
  for (l=0; l<numlogs; ++l) { 
   logpt = loglist[l] ;

   if  (logpt -> type == TROOT) { 
    s = logpt -> pop1 ;  
    calctime(s) ;
    continue ;
   }

   if  (logpt -> type == TADMIX) { 

   a = logpt -> pop1 ;  if (a<0) fatalx(" TADMIX badbug\n") ; 
   b = logpt -> pop2 ;  if (b<0) fatalx(" TADMIX badbug\n") ; 
   s = logpt -> pop3 ; 
   tt = calctime(s) ;   // evaluate time 
   if (tt<0) continue ;  
   
   settime1(b, s, 0) ; 
   settime1(a, s, 0) ; 
   calctime(a) ; 
   calctime(b) ; 
   continue ;   
  }  

   if (logpt -> type == TLIFT) { 
   t = logpt -> pop1 ;  if (t<0) fatalx("badbug\n") ; 
   s = logpt -> pop2 ; 
   tt = calctime(s) ; 
   y = -77 ; 
   if (tt>0) {
    y = callength(logpt -> tau1, logpt -> thetaval1) ; 
    t2  = settime1(t, s, y) ; 
    x = strcmp(vnames[t], "R") ; 
   }

   s = logpt -> pop3 ;   
   tt = calctime(s) ; 
   if (tt>0) { 
    y = callength(logpt -> tau2, logpt -> thetaval2) ; 
    t2 = settime1(t, s, y) ; 
    x = strcmp(vnames[t], "R") ; 
   }
   calctime(t) ; 
  }
  continue ; 
 }
 timecheck = 0 ; 

 for (s=0; s<numvertex; ++s) {  
  t = nnint(calss[s][0]) ; 
  if (t>1) timecheck += (t-1) ;
 }

}

double calctscore() 
{   

// double timesig = 30, timescore , timegscore, timegmul=1.0 ;
   

 int k ; 
 double tsc, y0, y1, y2, ym ;

 timescore = timegscore = 0 ; 
//  printf("zztch %d\n", timecheck) ; 
 if (timecheck <= 0) return 0 ;  
 
  tsc = 0 ; 
        

  for (k=0; k<numvertex; ++k) {  
   y0 = calss[k][1] ; 
   if (y0<1.5) continue ; 
    y1 = calss[k][1] ; 
    y2 = calss[k][2]  ;
     ym = y1/y0 ; 
     y2 -= y0*ym*ym ; 
     y2 /= y0 ;  // variance
     tsc += y2 ; 
  }

  timescore = timegscore = tsc ; 
  timescore /=  timesig ;   
  tsc = timescore * timegmul ; 

  return tsc  ;   

}
char *mytemp (char *qqq) 
{
  char ss[MAXSTR] ; 
  int t ; 

  t = (int) getpid() ; 
  sprintf(ss, "/tmp/%s.%d", qqq, t) ; 
  return strdup(ss) ;
}




