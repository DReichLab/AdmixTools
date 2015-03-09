#ifdef __cplusplus
extern "C" {
#endif
#include <ctype.h>
#include <nicklib.h>
#define IDSIZE 40

#ifndef ADMUTILS

typedef struct {
  char ID[IDSIZE];
  char *egroup ;
  char gender;     /* 'M' or 'F' */
  double theta_mode;     /* most likely theta on mesh */
  double lambda_mode; /* mean of log10(lambda) from probability distribution */

  double Xtheta_mode;     /* most likely theta on mesh */
  double Xlambda_mode;  /* most likely lambda on mesh */

  int idnum ;
  int affstatus;      /* affected status */
  int ignore ;        /* YES => do not use */
  double thetatrue ;
  double Xthetatrue ;
  double lambdatrue ;
  double Xlambdatrue ;
  double totgamms[3] ;
  double totscore ;
  double rawqval ;
  double qval ;
} Indiv; 

typedef struct {
  char ID[IDSIZE];
  int chrom;
  double genpos;
  double physpos;
  double aa_cauc_freq;    /* frequencies in ancestral pop to AA */
  double aa_af_freq;
  double cauc_freq;    /* frequencies of variant allele */
  double af_freq;
  double  cftrue ;
  double  aftrue ;
  double  aa_cftrue ;
  double  aa_aftrue ;
  int  markernum ; /* marker number */
  char *pbuff ; 
  char *ebuff ;   /* for random ethnic path */
  int  isfake ;   /* 1 if fake marker  else 0 */
  int  isrfake ;  
/* real marker treated as fake.  Used for 2D scoring */
  int ignore ;
  int  ngtypes ;  /* 	number of gtypes */
  int  *gtypes ;
  int  af_nn[2] ;  
  int  cauc_nn[2] ;  
  double *modelscores ;
  double *totmodelscores ;
  double score ;
  double weight ;
  double estgenpos ; 
  double estdis ; 
  double dis ; 
  double esum ;  
  double lsum ;
  double gpsum ; 
  double gpnum ;
  void *pcupt ;
  int tagnumber ;
  char alleles[2] ;
} SNP; 

typedef struct {
 char ID[IDSIZE] ;
 Indiv *father ; 
 Indiv *mother ; 
 Indiv *child ;
 int findex ;
 int mindex ;
 int cindex ; 
 int base ;
}   TRIO ;
// ?index into original Indiv array.  base is work variable used in phasetrio to store init index in new array

typedef struct {  
 char gname[IDSIZE] ;  
 SNP **snpm ;
 Indiv **indivm ;  
 int numsnps; 
 int numindivs ; 
 int rlen ;
 int fdes ;
 int snpindex ;
 unsigned char *buff ;
} genofile ;

typedef struct {  
 double xd[4] ;
 double xc[9] ;
 double ap1 ;
 double ap2 ;
 double cp1 ;
 double cp2 ;
 double rpowersum ;
 double crpowersum ;
 double gammasum[2] ;
 double gammanum[2] ;
 int  pubfmodern ;
}  SNPMC ;
// gammasum for cases/controls
#endif
#define ADMUTILS

void loadstats(FILE *statsfile, Indiv *indiv_array, int *numindivs);
void loadXstats(FILE *Xstatsfile, Indiv *indiv_array, int numindivs, int *numloaded);

void sett1(double *tt, double theta, int numstates) ;
void sett1r(double *tt, double theta, int numstate, double risks) ;
void gettln(SNP *cupt, Indiv *indx, 
  double *ptheta, double *plambda, int *pnumstates, int *pignore)  ;

void puttln(SNP *cupt, Indiv *indx, 
  double theta, double lambda) ;


/* UTILITY FUNCTIONS */

int countcol (char *fname);
int countcolumns (FILE *fp);

void fataly(const char *name);
int compare_doubles (const void *a, const void *b);

void pcheck (char *name, char x) ;
void printm(double **M, int numstates) ;
int numvalids(Indiv *indx, SNP **snpmarkers, int fc, int lc)  ; 
void gethpos(int *fc, int *lc, SNP **snpm, int numsnps,  
 int xchrom, int lo, int hi)  ; 
int numvalidgtypes(SNP *cupt) ;
double malefreq(Indiv **indivmarkers, int numindivs)  ;
int isimatch(int a, int b) ;   
void makedir(char *dirname) ;   
int indxindex(char **namelist, int len, char *strid)  ;
int indindex(Indiv **indivmarkers, int numindivs, char *indid)  ;
int snpindex(SNP **snpmarkers, int numsnps, char *snpid) ;
void freesnpindex() ;
int ignoresnp(SNP *cupt) ;
double entrop(double *a, int n)  ;
double xxlog2(double t)  ;
void testnan(double *a, int n) ;
void hap2dip(SNP *cupt) ; 
void flipalleles(SNP *cupt) ;
void flipalleles_phased(SNP *cupt) ;
int   getgtypes(SNP *cupt, int k) ;
void  putgtypes(SNP *cupt, int k, int val) ;
int   getep(SNP *cupt, int k) ;
void  putep(SNP *cupt, int k, int val) ;
int hasharr(char **xarr, int nxarr)   ;
void wbuff(unsigned char *buff, int num, int g) ;  
int rbuff(unsigned char *buff, int num)  ;
int ridfile(char *fname) ;
double hwcheck(SNP *cupt, double *cc) ;
double hwcheckx(SNP *cupt, Indiv **indm, double *cc) ;
void cntit(double *xc, SNP *cupt1, SNP *cupt2) ;
// dup routines 
void setfastdupnum(int num) ;
void setfastdupthresh(double thresh, double kill) ;
void killxhets(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs) ;
void fastdupcheck(SNP **snpmarkers, Indiv **indivmarkers, int numsnps, int numindivs) ;
int grabgtypes(int *gtypes, SNP *cupt, int numindivs) ;
int kcode(int *w, int len, int base) ;
void cdup(SNP **snpm, Indiv **indm, int nsnp, int *buff, int lbuff) ;
void printdup(SNP **snpm, int nsnp, Indiv *inda, Indiv *indb, int nmatch, int nnomatch); 
void killdup(Indiv *inda, Indiv *indb, SNP **snpm, int nsnp) ;
double kurtosis(double *a, int n) ;
int getlist(char *name, char **list) ; 
void printvers(char *progname, char *vers) ;
int numvalidind(Indiv **indivmarkers, int  numind)   ;
void numvalidgtallind(int *x, SNP **snpm, int numsnps, int numind) ; 
int numvalidgtind(SNP **snpm, int numsnps, int ind)  ;
int numvalidgt(Indiv **indivmarkers, SNP *cupt)   ;
int numvalidgtx(Indiv **indivmarkers, SNP *cupt, int affst)  ;
int isxmale(SNP *cupt, Indiv *indx) ;

void printmatz(double *ww, char **eglist, int n) ;
void printmatz5(double *ww, char **eglist, int n) ;
void printmatz10(double *ww, char **eglist, int n) ;
char *get3(char *ss) ;
char *getshort(char *ss, int n) ;


#undef max 
#define max(A,B)  ((A) > (B) ? (A) : (B))

#define MAXNUMR  200 
// max number models

#define CNULL '\0' 
#ifdef __cplusplus
}
#endif
