#ifdef __cplusplus
extern "C" {
#endif
#include "admutils.h"  


int makeeglist(char **eglist, int maxnumeg, Indiv **indivmarkers, int numindivs)  ; 
int mkeglist(Indiv **indm, int numindivs, char **eglist) ;
void seteglist(Indiv **indm, int nindiv, char *eglistname)  ;
void seteglistv(Indiv **indm, int nindiv, char *eglistname, int val) ;
int loadlist(char **list,  char *listname)  ;
int  loadlist_type(char **list, char *listname, int *ztypes, int off)   ;
#ifdef __cplusplus
}
#endif
