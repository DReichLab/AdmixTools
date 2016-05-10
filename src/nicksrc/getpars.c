#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>

#define MAXSTR  5000
#define MAXPARS 200

#include "nicklib.h"
#include "getpars.h"

/** 
 a very simple keyword parameter cracker 
 comments with line starting #  
 syntax:  keyword param(s)   
 not designed to read large databases (needs hash code for that)
*/
void subcolon (char *ss);
void stripcomment (char *str);
int findpname (phandle * pp, char *parname);
int indxstring (char **ppars, int npars, char *ww);

static int parchange = 0;
static int debug = NO;

phandle *
openpars (char *fname)

/* constructor */
{
  phandle *pp;
  FILE *ff;

  char line[MAXSTR + 1];
  char str[MAXSTR];
  char ww[MAXSTR];
  char rest[MAXSTR];

  char *ppars[MAXPARS];
  char *pdata[MAXPARS];

  int npars = 0, i;
  int len, plen, t;

  pp = (phandle *) malloc (sizeof (phandle));
  ff = pp->fx = fopen (fname, "r");
  if (ff == NULL) {
    perror ("Can't open file\n");
    fatalx ("can't open %s\n", fname);
  }

  line[MAXSTR] = '\0';          // defensive programmming
  while (fgets (line, MAXSTR, ff) != NULL) {

    len = strlen (line);

    if (isspace (line[len - 1]))
      line[len - 1] = '\0';

    if (first_word (line, ww, rest) > 0) {

      if (ww[0] == '#')
        continue;
      /*AT: 12/2/04: Adding check to make sure that the parameter name ends in : */
      plen = strlen (ww);
      if (ww[plen - 1] != ':')
        printf
          ("**warning: dubious parameter, please check the parameter name %s\n",
           ww);
      t = indxstring (ppars, npars, ww);
      if (t >= 0)
        fatalx ("duplicate parameter: %s\n", ww);
      ppars[npars] = strdup (ww);

      striptrail (rest, ' ');   /* no trailing blanks */
      stripcomment (rest);
      pdata[npars] = strdup (rest);
      ++npars;

      if (debug)
        printf ("param %d %s %s\n", npars, ppars[npars - 1],
                pdata[npars - 1]);

    }
  }
  pp->numpars = npars;

  if (npars > 0) {
    ZALLOC (pp->ppars, npars, char *);
    ZALLOC (pp->pdata, npars, char *);
  }
  else {
    fprintf (stderr, "***warning: no parameters in %s\n", fname);
  }


  for (i = 0; i < npars; i++) {
    pp->ppars[i] = strdup (ppars[i]);
    pp->pdata[i] = strdup (pdata[i]);


    /*  printf("zz: %d %s %s\n",i,ppars[i],pp->ppars[i]) ; */

  }

  for (i = 0; i < npars; i++) {
    free (ppars[i]);
    free (pdata[i]);
  }

  return pp;

}

void
closepars (phandle * pp)

/* destructor */
{
  int n, i;

  fclose (pp->fx);
  n = pp->numpars;
  for (i = 0; i < n; i++) {
    free (pp->ppars[i]);
    free (pp->pdata[i]);
  }

  free (pp->ppars);
  free (pp->pdata);

  free (pp);
  pp = NULL;


}

void
stripcomment (char *str)
{
  int i, len;

  len = strlen (str);

  for (i = 0; i < len; i++) {
    if (str[i] == '#') {
      str[i] = '\0';
      return;
    }
  }
}

#define MAXFIELD 1000

int
getstring (phandle * pp, char *parname, char **strng)
{

  char *field[MAXFIELD];
  int n, kode;

  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;
  n = splitup (pp->pdata[kode], field, MAXFIELD);
  *strng = strdup (field[0]);
  freeup (field, n);
  return 1;
}


int
getint (phandle * pp, char *parname, int *kret)
{

  char *field[MAXFIELD];
  char str[MAXSTR];
  int n, kode;

  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;
  n = splitup (pp->pdata[kode], field, MAXFIELD);
  strcpy (str, field[0]);
  freeup (field, n);
  if (strcmp (str, "YES") == 0) {
    *kret = YES;
    return 1;
  }
  if (strcmp (str, "NO") == 0) {
    *kret = NO;
    return 1;
  }

  *kret = atoi (str);
  return 1;
}

int
getints (phandle * pp, char *parname, int *aint, int nint)
{
  char *field[MAXFIELD];
  int n, kode, i;
  char str[MAXSTR];

  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;

  strcpy (str, pp->pdata[kode]);
  subcolon (str);
  n = splitup (str, field, MAXFIELD);

  if (nint < n)
    fatalx ("(getints) wrong number of ints in line\n");

  for (i = 0; i < n; i++) {
    aint[i] = atoi (field[i]);
  }
  freeup (field, n);
  return n;
}

int
getintss (phandle * pp, char *parname, int *aint, int *xint)

/* white space separated */
{
  char *field[MAXFIELD];
  int n, kode, i;
  char str[MAXSTR];

  *xint = 0;
  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;
  strcpy (str, pp->pdata[kode]);
  subcolon (str);
  n = splitup (str, field, MAXFIELD);
  for (i = 0; i < n; i++) {
    aint[i] = atoi (field[i]);
  }
  freeup (field, n);
  *xint = n;
  return 1;
}

int
getdbls (phandle * pp, char *parname, double *dbl, int ndbl)
{
  char *field[MAXFIELD];
  int n, kode, i;
  char str[MAXSTR];

  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;
  strcpy (str, pp->pdata[kode]);
  subcolon (str);
  n = splitup (str, field, MAXFIELD);
  if (ndbl != n)
    fatalx ("(getdbls) wrong number of dbls in line\n");
  for (i = 0; i < n; i++) {
    dbl[i] = atof (field[i]);
  }
  freeup (field, n);
  return 1;
}

int
getdblss (phandle * pp, char *parname, double *dbl, int *ndbl)

/* separated no check */
{
  char *field[MAXFIELD];
  int n, kode, i;
  char str[MAXSTR];

  *ndbl = 0;
  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;
  strcpy (str, pp->pdata[kode]);
  subcolon (str);
  n = splitup (str, field, MAXFIELD);
  for (i = 0; i < n; i++) {
    dbl[i] = atof (field[i]);
  }
  freeup (field, n);
  *ndbl = n;
  return 1;
}

void
subcolon (char *ss)

/* substitute ' ' for ':' or ',' */
{
  int i, l;
  l = strlen (ss);
  for (i = 0; i < l; i++) {
    if (ss[i] == ':')
      ss[i] = ' ';
    if (ss[i] == ',')
      ss[i] = ' ';
  }
}


int
getdbl (phandle * pp, char *parname, double *dbl)
{
  char *field[MAXFIELD];
  int n, kode;

  kode = findpname (pp, parname);
  if (kode < 0)
    return kode;
  n = splitup (pp->pdata[kode], field, MAXFIELD);
  *dbl = atof (field[0]);
  freeup (field, n);
  return 1;
}

int
findpname (phandle * pp, char *parname)
{
  int k;

  for (k = 0; k < pp->numpars; k++) {
    if (strcmp (parname, pp->ppars[k]) == 0)
      return k;
  }
  return -1;
}

void
writepars (phandle * pp)
{
  int k;
  if (pp == NULL)
    fatalx ("(writepars) phandle not open\n");
  for (k = 0; k < pp->numpars; k++) {
    printf ("%s %s\n", pp->ppars[k], pp->pdata[k]);
  }
}

void
fwritepars (FILE * fff, phandle * pp)
{
  int k;
  if (pp == NULL)
    fatalx ("(writepars) phandle not open\n");
  for (k = 0; k < pp->numpars; k++) {
    fprintf (fff, "%s %s\n", pp->ppars[k], pp->pdata[k]);
  }
}


void
dostrsub (phandle * pp)
{

/** 
 massage phandle data structure  
 to do string substitutions of 
 UPPER case paramewters
*/

  int n, nu, nl;
  int *isupp, *islow, *slenm, *indx;
  int i, j, k, l, ind;
  int nchange = 0;
  char *inx, *outx;
  n = pp->numpars;

  ZALLOC (isupp, n, int);
  ZALLOC (islow, n, int);
  ZALLOC (slenm, n, int);
  ZALLOC (indx, n, int);

  nu = nl = 0;
  for (k = 0; k < pp->numpars; k++) {
    if (upstring (pp->ppars[k])) {
      isupp[nu] = k;
      l = strlen (pp->ppars[k]);
      slenm[nu] = -(10000 * l - k);

/** 
 after sorting longest first then order in file
 ensures sort is stable
*/

      ++nu;
    }
    else {
      islow[nl] = k;
      ++nl;
    }
  }
  if (nu == 0)
    return;
  isortit (slenm, indx, nu);
  for (i = 0; i < nu; i++) {
    ind = indx[i];
    k = isupp[ind];
    inx = strdup (pp->ppars[k]);
    l = strlen (inx);
    if (inx[l - 1] == ':')
      inx[l - 1] = '\0';
    outx = strdup (pp->pdata[k]);
    for (j = 0; j < nl; j++) {
      k = islow[j];
      nchange += substring (&(pp->pdata[k]), inx, outx);
      if (nchange > 0)
        break;
    }
    free (inx);
    free (outx);
  }
  free (isupp);
  free (islow);
  free (slenm);
  free (indx);

  parchange += nchange;
  if (parchange > 10000)
    fatalx ("(getpars) dostrsub looping\n");
  if (nchange > 0)
    dostrsub (pp);

}

int
upstring (char *ss)

/* 
 YES if at least one upper case character 
 and no lower case  
*/
{
  int nupper = 0;
  int i;
  for (i = 0; i < strlen (ss); i++)
  {
    if (islower (ss[i]))
      return NO;
    if (isupper (ss[i]))
      ++nupper;
  }
  if (nupper > 0)
    return YES;
  return NO;

}
