#ifdef __cplusplus
extern "C" {
#endif
#include <search.h>  

int xfindit(char *ss)   ; 
void xloadsearch(char **ss, int n)   ; 
void xdestroy()  ;

void xhcreate (int n)  ;
void xhdestroy()  ;
ENTRY *xhsearch(ENTRY item, ACTION act)  ;

int xlookup(char *key, ACTION act)  ;
int xhash (char *key)  ;
int xhash1(int ww)  ;
int xhash2 (int x)  ;
int xcshift(int x, int shft)  ;
int stringhash(char *key) ;

#ifdef __cplusplus
}
#endif
