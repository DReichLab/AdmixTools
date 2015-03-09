#ifdef __cplusplus
extern "C" {
#endif
#include <nicklib.h>
#include <math.h>
#include "admutils.h" 

extern int verbose ;

double xpest(double **gg, int *gobs, int *na, int *nb,  
 int neq, double *ppa, double *ppb) ;

void mk2from3ml(double *xd, double *xc, double p, double pp)  ;
void mk2from2(double *xd, double *xc ) ;

void loadpprob(double *pprob, double pa, double pb) ;
void gen3(double *ww, double a, double b) ;

double xpest2like(double **gg, int *gobs, int *na, int *nb,  
 int *iscasearr, 
 int neq, double ppa, double ppb, double risk)  ;


// clean up SANS when finalized 

typedef struct {
  SNP *cupt ;
  int numsamps ;
  double admbayessc ; 
  double *baymodelsc ;
  double admfsc ; 
  double admzscore ; 
  double admbsc ;
  double admfccsc ;
  double simpsc[2] ;
  double admyl[3] ;
  double lrmax ;
  double lrsig ; 
  double maxlod ;
/* now start of fine-mapping scores */
  double gscore ;
  double gcheck ; 
  double gbayes ;
} SANS ; 

#ifdef __cplusplus
}
#endif
