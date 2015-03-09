
typedef struct {
  double S0 ; 
  double S1 ; 
  double S2 ; 
  double S11 ; 
  double S12 ; 
  double S22 ; 
  double m1 ;
  double m2 ;
  double v11 ;
  double v12 ;
  double v22 ;
  double corr ;
  double Z ;
} CORR; 

int calccorr(CORR *corrpt, int mode, int ztrans)  ;
void printcorr(CORR *corrpt)  ;
void clearcorr(CORR *corrpt)  ;
void addcorr(CORR *corrpt, double x1, double x2) ; 
void addcorrn(CORR *corrpt, double x1, double x2, double yn) ; 
void minuscorr(CORR *out, CORR *c1, CORR *c2) ;

double lddip(double *xc) ; 
double zdip(double *xc) ; 
void setzdipmode(int mode) ;
void setzdphasedmode(int mode) ;
