#  include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <globals.h>
#include <mcmcpars.h>
#include <getpars.h>

#include "badpairs.h"
#include "admutils.h"
#include "mcio.h"  

#define WVERSION   "420" 

// badpairsname added 

#define MAXFL  50   
#define MAXSTR  512

extern int packmode ;

char *trashdir = "/var/tmp" ;
extern int verbose  ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 
int randommode = NO ;

char *instem = NULL ; 
char  *genotypename = NULL ;
char  *snpname = NULL ;
char  *indivname = NULL ;

FILE *ofile ;

char *weightname = "weights.txt" ;
int seed = 89 ;

double fakespacing = 0.0 ;

char  unknowngender = 'U' ;

void readcommands(int argc, char **argv) ;
int getmsd(double *pmean, double *psd, SNP **snpmarkers, int numsnps, char *psname)  ;
void msd(double *mean, double *sd, double *a, int n)  ;
int sampleg (double p)  ;


double jcorr(double *pjest, double *pjsig, 
 double *pest, double *a, double *b, int n) ; 
double jvar(double *pjest, double *pjsig, 
 double *pest, double *a,  int n) ; 
double jbcorr(double *pjest, double *pjsig, 
 double *pest, double *a, double *b, int n) ; 
double javar(double *pjest, double *pjsig, 
 double *pest, double *a,  int n) ; 


int main(int argc, char **argv)
{

  int i, j, k, g, a, b, c, t, x, gf, gm ; 
  SNP *cupt ;
  Indiv *indx ;

  int numvind, nignore, numrisks = 1 ;
  char ***trios ; 
  int ntrios ; 
  int *findex, *mindex, *cindex ; 
  char *sx ; 
  double y, y1, y2, ysd, yvar ;
  double **gg, **g1, **g2, *wts, *ww, *w1, *w2, *wts2 ;
  double *pmean, *psd, *pp ;
  double *scor, *fsc, *msc, *csc, *sdu ;
  double *fz, *mz, *fmtt, *truett, xa, xb, z1, z2 ;
  double *ymul ;
  double ymem ; 
  double  *weven, *wodd ;
  double  *ceven, *codd ;
  double  *seven, *sodd ;
  int dof ; 
  double ytail ; 
  double rtrue, aest, vest, best, beta ;
  double ajest, ajsig, az, bjest, bjsig, bz ;

  ofile = stdout; 
  packmode = YES ;
  readcommands(argc, argv) ;

  cputime(0) ;
  calcmem(0) ;

  SRAND(seed) ;

// fakespacing 0.0 (default)

  if (instem != NULL) { 
   setinfiles(&indivname, &snpname, &genotypename, instem) ; 
  } 

  numindivs = getindivs(indivname, &indivmarkers) ;

  setgenotypename(&genotypename, indivname) ;

  printf("genotypename:  %s\n", genotypename) ;

   getgenos(genotypename, snpmarkers, indivmarkers, 
    numsnps, numindivs, nignore) ;

  numvind = numvalidind(indivmarkers, numindivs) ;
  if (randommode) printf("randommode set!\n") ;
  printf("\n\n") ;
  printf("numindivs: %d valid: %d numsnps: %d nignore: %d\n" ,
    numindivs, numvind, numsnps, nignore) ; 

  getweights(weightname, snpmarkers, numsnps) ; 

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of ascore: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
  return 0 ;
}

int sampleg (double p) 
{
  int g1, g2 ;  
  
  g1 = prob1(p) ; 
  g2 = prob1(p) ; 

  return g1 + g2 ;



}
void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:s:rvV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

      case 'r':
	randommode = YES ;
	break; 

      case 's':  
       seed = atoi(optarg) ;
       break ;


      case '?':
	printf ("Usage: bad params.... \n") ;
	fatalx("bad params\n") ;
      }
  }

         
   pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getstring (ph, "instem:", &instem);
   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;

   writepars(ph) ;
   closepars(ph) ;

}

