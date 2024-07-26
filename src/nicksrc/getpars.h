typedef struct {
 int numpars ;
 FILE *fx ;
 char **ppars ;
 char **pdata ;
}  phandle ;

void writepars(phandle  *pp) ;
void closepars(phandle  *pp) ;
phandle *openpars(char *fname) ;

int getlongstring(phandle *pp, char *parname, char **kret)  ;
// whole of line
int getstring(phandle *pp, char *parname, char **kret)  ;
int getint(phandle *pp, char *parname, int *kret)  ;
int getlong(phandle *pp, char *parname, long *kret)  ;
int getints(phandle *pp, char *parname, int *aint, int nint) ;
int getintss(phandle *pp, char *parname, int *aint, int *xint) ;

int getdbl(phandle *pp, char *parname, double *dbl)  ;
int getdbls(phandle *pp, char *parname, double *dbl, int ndbl) ;
int getdblss(phandle *pp, char *parname, double *dbl, int *ndbl) ;
int subst(char *outstr, char *instr, char *ins, char *outs)  ;
void dostrsub(phandle *pp)  ; 
int upstring (char *ss)  ;
void subcolon(char *ss)  ;

#define MAXFIELD 1000
