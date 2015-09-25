#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>

int splitup (char *strin,  char *strpt[],int maxpt) ;
int splitupx(char *strin, char **spt, int maxpt, char splitc)  ;
int splitupwxbuff(char *strin, char **spt, int maxpt, char *bigbuff, int bigbufflen)  ;
int splitupxbuff(char *strin, char **spt, int maxpt, char splitc, char *bigbuff, int bigbufflen)  ;
int oldsplitup (char *strin,  char *strpt[],int maxpt) ;
void freeup (char *strpt[],int numpt) ;
int split1 (char *strin,  char *strpt[], char splitc);
int first_word(char *string, char *word, char *rest) ;
char *fnwhite (char *ss) ;
char *fwhite (char *ss) ;
char *ftab (char *ss) ;
int NPisnumber (char c) ;
int isnumword (char *str)  ;
void fatalx( char *fmt, ...) ;
long seednum() ;
void printbl(int n) ;
void printnl() ;
void striptrail(char *sss, char c) ;
void catx(char *sout, char **spt, int n) ;
void catxx(char *sout, char **spt, int n) ;
void catxc(char *sout, char **spt, int n, char c) ;
void makedfn(char *dirname, char *fname, char *outname, int maxstr) ;
int substring (char **ap, char *inx, char *outx) ;
int numcols (char *name) ;
int numlines(char *name) ;
void openit(char *name, FILE **fff, char *type)  ;
int getxx(double **xx, int maxrow, int numcol, char *fname) ;
int getss(char  **ss, char *fname) ;
double clocktime() ;  // cpu time in seconds
void crevcomp(char *sout, char *sin)  ;  
int indxstring(char **namelist, int len, char *strid)  ;
char *strstrx(char *s1, char *s2)  ;  // case insensitive strstr
int getxxnames(char ***pnames, double **xx, int maxrow, int numcol, char *fname);
int getjjnames(char ***pnames, int **xx, int maxrow, int numcol, char *fname);
int getxxnamesf(char ***pnames, double **xx, int maxrow, int numcol, FILE *fff) ;
int getnames(char ****pnames, int maxrow, int numcol, char *fname) ;
char num2base (int num) ; 
int base2num(char c) ;
char *int_string(int a, int len, int base) ;
char *binary_string(int a, int len) ;
int string_binary(char *sx) ;
void freestring (char **ss) ;
void copystrings(char **sa, char **sb, int n) ;
void printstrings(char **ss, int n)  ;
int ridfile(char *fname) ; 
char compbase(char x) ;
void mkupper(char *sx) ;
void mklower(char *sx) ;
int iubdekode(char *a, char iub) ;  
int char2int(char cc) ;
char int2char(int x) ;




#define ZALLOC(item,n,type)      if ((item = (type *)calloc((n),sizeof(type))) == NULL) \
                                        fatalx("Unable to allocate %d unit(s) for item \n",n)

#undef MAX
#undef MIN

#define MAX(a,b)   ( (a) < (b) ?  (b) : (a) ) 
#define MIN(a,b)   ( (a) < (b) ?  (a) : (b) ) 
#define YES  1
#define NO   0
#define TRUE   1
#define FALSE  0
#define CNULL  '\0' 
#define CNL  '\n' 
#define CTAB  '\t' 

#ifdef __cplusplus
}
#endif
