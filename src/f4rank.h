typedef struct {
  int nl ; 
  int nr ; 
  char *lbase ;
  char *rbase ;
  char **lpops ;
  char **rpops ;
  int rank ; 
  double dofjack ; 
  double dof ; 
  double dofdiff ; 
  double chisq ;
  double chisqdiff ;
  double *A ; 
  double *B ;
  double *mean ; 
  double *resid ; 
  struct F4INFO *bestparent ;
  struct F4INFO *bestchild ;
} F4INFO  ;
// dofjack is dof from jackknife blocks 
// dof is dof from approx chisq     


void doranktest(double *mean, double *var, int m, int n, int rank, F4INFO *f4pt) ;
void doranktestfix(double *mean, double *var, int m, int n, int rank, F4INFO *f4pt, int *vfix) ;
double ranktest(double *mean, double *var, int m, int n, int rank, double *pA, double *pB) ;
double ranktestfix(double *mean, double *var, int m, int n, int rank, double *pA, double *pB, int *vfix) ;
void normab(double *A, double *B, int m, int n, int rank)  ;
int dofrank(int m, int n, int rank) ;
void f4info_init(F4INFO *f4pt, int nl, int nr, char **popllist, char **poprlist, int rank) ;
void printf4info(F4INFO *f4pt)  ;


