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
void ffprint (FILE *fff, char *fmt, ...) ; 
void enuf( char *fmt, ...) ;
void fatalx( char *fmt, ...) ;
int docommand( char *fmt, ...) ;
long seednum() ;
void printbl(int n) ;
void printnl() ;
void striplead(char *sss, char c) ;
void striptrail(char *sss, char c) ;
void catx(char *sout, char **spt, int n) ;
void catxx(char *sout, char **spt, int n) ;
void catxc(char *sout, char **spt, int n, char c) ;
void makedfn(char *dirname, char *fname, char *outname, int maxstr) ;
int substring (char **ap, char *inx, char *outx) ;
int mapstrings(char **pstr, char **insub, char **outsub, int n)  ;
int upstring (char *ss)  ; 
int numcols (char *name) ;
int numlines(char *name) ;
int openit_trap (char *name, FILE ** fff, char *type); 
void openitntry (char *name, FILE ** fff, char *type, int ntry) ;
void openit(char *name, FILE **fff, char *type)  ;
int  ftest(char *aname) ;
void fcheckr(char *name) ;
void fcheckw(char *name) ; 
int getxx(double **xx, int maxrow, int numcol, char *fname) ;
int getss(char  **ss, char *fname) ;
int  loadlist(char **list, char *listname)    ;  // with dup check
void printdups(char **list, int n)  ;
int checkdup(char **list, int n)  ;
double clocktime() ;  // cpu time in seconds
void crevcomp(char *sout, char *sin)  ;  
int indxstring(char **namelist, int len, char *strid)  ;
int indxstringr(char **namelist, int len, char *strid)  ;
char *strstrx(char *s1, char *s2)  ;  // case insensitive strstr
int getxxnames(char ***pnames, double **xx, int maxrow, int numcol, char *fname);
int getjjnames(char ***pnames, int **xx, int maxrow, int numcol, char *fname);
int getxxnamesf(char ***pnames, double **xx, int maxrow, int numcol, FILE *fff) ;
int getnameslohi(char ****pnames, int maxrow, int numcol, char *fname, int lo, int hi) ;
int getnamesstripcolon(char ****pnames, int maxrow, int numcol, char *fname, int lo, int hi) ;
int getnames(char ****pnames, int maxrow, int numcol, char *fname) ;
char num2iub (int num) ; 
char revchar(char c) ;
int  iub2num(char c) ;
char num2base (int num) ; 
int base2num(char c) ;
char *int_string(int a, int len, int base) ;
char *binary_string(int a, int len) ;
int string_binary(char *sx) ;
void freestring (char **ss) ;
void copystrings(char **sa, char **sb, int n) ;
void printstringsw(char **ss, int n, int slen, int width)  ;
void printstrings(char **ss, int n)  ;
void printstringsx(char **ss, int n)  ;
int ridfile(char *fname) ; 
char compbase(char x) ;
void mkupper(char *sx) ;
void mklower(char *sx) ;
int iubdekode(char *a, char iub) ;  
int isiub(char iub) ;  
int isiub2(char iub) ;  
int iubcbases(char *cbases, char iub) ;
int ishet(char c) ;
int cttype(char cc) ;
int char2int(char cc) ;
char int2char(int x) ;
void chomp(char *str) ;

int numcmatch(char *cc, int len, char c) ;
int numcnomatch(char *cc, int len, char c) ;
char *strnotchar(char *s, char c) ;
char *findupper(char *s) ;
char *fgetstrap(char *buff, int maxlen, FILE *fff, int *ret)  ;
char readtonl(FILE *fff) ; 
int  filehash(char *name) ;
char *mytemp (char *qqq) ; 
void printslurmenv ()  ; 
int getfline(char *ss, char *fname, int maxstr) ;
int copyfs(char *infile, FILE *fff)  ;
int getxxq(double **xx, int maxrow, int numcol, char *fname) ; 
int numcolsq (char *name) ;
int getdata(char *buff, int nbytes, char *fname)  ;
int putdata(char *buff, int nbytes, char *fname)  ;
void writestrings(char *fname, char **ss, int n)  ;


#define ZALLOC(item,n,type)      if ((item = (type *)calloc((n),sizeof(type))) == NULL) \
                                        fatalx("Unable to allocate %ld unit(s) for item \n", (long) n)

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

