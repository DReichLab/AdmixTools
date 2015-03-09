#ifdef __cplusplus
extern "C" {
#endif
#include <nicklib.h>
#include "admutils.h" 

void dobadpairs(char *badpairsname, SNP **snpm, int numsnps) ;
void dogood(char *goodsnpname, SNP **snpm, int numsnps)  ;
void getsnpsc(char *snpscname, SNP **snpm, int numsnps) ;
int killsnps(Indiv **indivmarkers, SNP **snpmarkers, int numsnps, int mincasenum)  ;
void loadbadpsc(SNP **snpm, int numsnps, int rmode, char *gname)  ;

double entrop(double *a, int n)  ;
double xxlog2(double t) ;
double mutx(double *dd) ;
double mutxx(double *dd, int m , int n) ;


#ifdef __cplusplus
}
#endif
