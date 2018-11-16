#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdarg.h>
#include <sys/times.h>
#include <time.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <xsearch.h>   


#define MAXSTR 10000
#define MAXFF  500 

#include "strsubs.h"
#include "vsubs.h"
#include "getpars.h"

extern int errno;


int
oldsplitup (char *strin, char **spt, int maxpt)

/**
 retained in case there are compatibility problems 
*/
{
  char *s1, *s2, *sx;
  char *str;
  int i, len, num;

  len = strlen (strin);
  if (len == 0)
    return 0;
  ZALLOC (str, 2 * len, char);
  num = 0;
  sx = strin;
  for (i = 0; i < maxpt; i++) {
    s1 = fnwhite (sx);
    if (s1 == NULL) {
      break;
    }
    s2 = fwhite (s1);
    if (s2 == NULL) {
      s2 = s1 + strlen (s1);
    }
    s2--;                       /* now points at last character of next word */
    len = s2 - s1 + 1;
    strncpy (str, s1, len);
    str[len] = '\0';
    spt[num] = strdup (str);
    ++num;
    sx = s2 + 1;
  }
  freestring (&str);
  return num;
}

void
freeup (char *strpt[], int numpt)

/** free up array of strings */
{
  int i;
  for (i = numpt - 1; i >= 0; i--) {
    if (strpt[i] != NULL)
      freestring (&strpt[i]);
  }
}

int
first_word (char *string, char *xword, char *xrest)

/*  first_word(string, *word, *rest)

        Break the string into the first word and the rest.  Both word and
rest begin with non-white space, unless rest is null.
        Return:
                0 means string is all white
                1 means word is non-white, but rest is white
                2 means word and rest are non-white
 
        If string and rest coincide, string will be overwritten
 
*/
{
  char *spt, x;
  char *ss = NULL, *sx;
  int l1, l2;

  ss = strdup (string);
  if (ss == NULL) {
    printf ("strdup fails\n");
    printf ("%s\n", string);
    fatalx ("first_word... strdup fails\n");
  }
  fflush (stdout);
  spt = ss;
  xword[0] = xrest[0] = '\0';
  if ((spt = fnwhite (ss)) == NULL) {
    free (ss);
    return 0;
  }
  sx = fwhite (spt);
  if (sx == NULL) {
    strcpy (xword, spt);
    free (ss);
    return 1;
  }
  l1 = sx - spt;
  l2 = strlen (sx) - 1;
  *sx = '\0';
  strcpy (xword, spt);
  if (l2 <= 0) {
    free (ss);
    return 1;
  }

  sx = fnwhite (sx + 1);
  if (sx == NULL) {
    free (ss);
    return 1;
  }
  strcpy (xrest, sx);
  free (ss);
  return 2;
}

char *
fnwhite (char *ss)

/* return first non white space */
{
  char *x;
  if (ss == NULL)
    fatalx ("fnwhite: logic bug\n");
  for (x = ss; *x != '\0'; ++x) {
    if (!isspace (*x))
      return x;
  }
  return NULL;
}

char *
ftab (char *ss)

/* return first tab  */
{
  char *x;
  int n;
  for (x = ss; *x != '\0'; ++x) {
    if (*x == CTAB)
      return x;
  }
  return NULL;
}

char *
fwhite (char *ss)

/* return first white space */
{
  char *x;
  int n;
  for (x = ss; *x != '\0'; ++x) {
    if (isspace (*x))
      return x;
  }
  return NULL;
}

static char Estr[MAXSTR];

void
fatalx (char *fmt, ...)
{
  va_list args;

  va_start (args, fmt);
  vsprintf (Estr, fmt, args);
  va_end (args);
  fflush (stdout);

  fprintf (stderr, "fatalx:\n%s", Estr);
  fflush (stderr);
#ifdef USER_EXIT
  exit(1);
#else
  abort ();
#endif
}

int
NPisnumber (char c)

/**
 returns 1 if - + or digit 
*/
{
  if (isdigit (c))
    return 1;
  if (c == '+')
    return 1;
  if (c == '-')
    return 1;

  return 0;
}

int
isnumword (char *str)
{

  int i, len, numpt;
  char c;
  len = strlen (str);

  numpt = 0;
  for (i = 0; i < len; i++) {
    c = str[i];

    if ((c == '.') && (numpt == 0)) {
      ++numpt;
      continue;
    }

    if (!NPisnumber (c))
      return NO;
    if (!isdigit (c) && (i > 0))
      return NO;
  }
  return YES;

}

long
seednum ()
{
  long a, b, c, d;
  struct tms tbuff;

  a = (long) getpid ();
  b = (long) getuid ();
  d = times (&tbuff);

  c = d ^ ((a + b) << 15);


  return c;

}

int
splitupwxbuff (char *strin, char **spt, int maxpt, char *bigbuff,
               int bigbufflen)
// splits by white space; No zero length strings
{
  char *sx, *sy;
  int num, len, k, klo;
  int empty = YES;

  len = strlen (strin);
  if (len >= bigbufflen)
    fatalx ("(splitupwxbuff) overflow\n%s", strin);
  strcpy (bigbuff, strin);
  num = 0;
  for (k = 0; k < len; ++k) {
    if (!isspace (bigbuff[k])) {
      empty = NO;
      klo = k;
      sx = bigbuff + k;
      break;
    }
  }
  if (empty)
    return 0;
  for (k = klo; k < len; ++k) {
    if (isspace (strin[k])) {
      bigbuff[k] = CNULL;
      if (num >= maxpt)
        break;
      spt[num] = sx;
      if (strlen (sx) > 0)
        ++num;
      sx = bigbuff + k + 1;
    }
  }
  if (num >= maxpt)
    return num;
  spt[num] = sx;
  if (strlen (sx) > 0)
    ++num;
  return num;
}

int
splitupxbuff (char *strin, char **spt, int maxpt, char splitc, char *bigbuff,
              int bigbufflen)
{
  char *sx, *sy;
  int num, len, k, klo;
  int empty = YES;

  len = strlen (strin);
  if (len >= bigbufflen)
    fatalx ("(splitupxbuff) overflow \n%s\n", strin);
  strcpy (bigbuff, strin);
  num = 0;
  for (k = 0; k < len; ++k) {
    if (strin[k] != splitc) {
      empty = NO;
      klo = k;
      sx = bigbuff + k;
      break;
    }
  }
  if (empty)
    return 0;
  for (k = klo; k < len; ++k) {
    if (strin[k] == splitc) {
      bigbuff[k] = CNULL;
      if (num >= maxpt)
        fatalx ("overflow\n");
      spt[num] = sx;
      sx = bigbuff + k + 1;
      ++num;
    }
  }
  if (num >= maxpt)
    fatalx ("overflow\n");
  spt[num] = sx;
  ++num;
  return num;
}

int
splitup (char *strin, char **spt, int maxpt)
{
  char *bigb, **qpt;
  int num, len, k;

  if (strin == NULL)
    return 0;
  len = strlen (strin);
  ZALLOC (bigb, len + 1, char);
  ZALLOC (qpt, maxpt + 10, char *);
  num = splitupwxbuff (strin, qpt, maxpt, bigb, len + 1);
  for (k = 0; k < num; ++k) {
    spt[k] = strdup (qpt[k]);
  }
  free (bigb);
  free (qpt);
  return num;
}

int
splitupx (char *strin, char **spt, int maxpt, char splitc)
{
  char *bigb, **qpt;
  int num, len, k;

  if (strin == NULL)
    return 0;
  len = strlen (strin);
  ZALLOC (bigb, len + 1, char);
  ZALLOC (qpt, maxpt + 10, char *);
  num = splitupxbuff (strin, qpt, maxpt, splitc, bigb, len + 1);
  for (k = 0; k < num; ++k) {
    spt[k] = strdup (qpt[k]);
  }
  free (qpt);
  free (bigb);
  return num;
}

int
split1 (char *strin, char *strpt[], char splitc)

/*
take a string and break it into 2 substrings separated by splitc ;
numpt is number of words returned  (1 or 2) 
*/
{
  char rest[MAXSTR], str[MAXSTR], ww[MAXSTR];
  int len, i, l;

  strncpy (str, strin, MAXSTR);
  len = strlen (strin);
  for (i = 0; i < len; i++) {
    if (str[i] == splitc) {
      l = i;
      strncpy (ww, str, l);
      ww[l] = '\0';
      strpt[0] = strdup (ww);
      l = len - (i + 1);
      if (l <= 0)
        return 1;
      strncpy (rest, str + i + 1, l);
      rest[l] = '\0';
      strpt[1] = strdup (rest);
      return 2;
    }
  }
  strpt[0] = strdup (strin);
  strpt[1] = NULL;
  return 1;
}

void
printbl (int n)
{
  int i;
  for (i = 0; i < n; i++) {
    printf (" ");
  }
}

void
printnl ()
{
  printf ("\n");
}

void
striptrail (char *sss, char c)

/** 
 strip out trailing characters 
 c will usually be ' '
*/
{
  int len, i;
  len = strlen (sss);
  for (i = len - 1; i >= 0; --i) {
    if (sss[i] != c)
      return;
    sss[i] = '\0';
  }
}

void
catx (char *sxout, char **spt, int n)
{
  int i;
  sxout[0] = CNULL;

  for (i = 0; i < n; i++) {
    strcat (sxout, spt[i]);
  }


}

void
catxx (char *sxout, char **spt, int n)

/** 
 like catx but with space between items 
*/
{
  int i;
  sxout[0] = CNULL;

  for (i = 0; i < n; i++) {
    strcat (sxout, spt[i]);
    if (i < (n - 1))
      strcat (sxout, " ");
  }
}

void
catxc (char *sxout, char **spt, int n, char c)

/** 
 like catx but with char c between items 
*/
{
  int i;
  char cc[2];

  sxout[0] = CNULL;

  cc[0] = c;
  cc[1] = CNULL;

  for (i = 0; i < n; i++) {
    strcat (sxout, spt[i]);
    if (i < (n - 1))
      strcat (sxout, cc);
  }
}

void
makedfn (char *dirname, char *fname, char *outname, int maxstr)

/** makes full path name.    
  If fname starts with '/' or dirname = NULL we 
  so nothing. 
  outname MUST be allocated of length at least maxstr 
*/
{
  char *ss;
  int len;

  if ((dirname == NULL) || (fname[0] == '/')) {

/* if fname starts with / we assume absolute pathname */
    len = strlen (fname);
    if (len >= maxstr)
      fatalx ("(makedfn) maxstr too short\n");
    strcpy (outname, fname);
    return;
  }
  len = strlen (dirname) + strlen (fname) + 1;
  if (len >= maxstr)
    fatalx ("(makedfn) maxstr too short\n");

  ss = outname;
  strcpy (ss, dirname);
  ss = ss + strlen (dirname);
  ss[0] = '/';
  ++ss;
  strcpy (ss, fname);
}

int
substringx (char **ap, char *inx, char *outx, int niter)

/** 
 *ap is original string 
 all occurrences of inx are substituted with outx 
 can loop so be careful !!  

 NB.  ap must be on heap.  Fixed allocation not supported 
*/
{
  char *a, *pt;
  char *str;
  int len, off, x;

  if (niter > 50)
    fatalx ("bad string replacement\n %s\n", *ap);

  a = *ap;
  len = strlen (a) + strlen (inx) + strlen (outx) + 1;
  pt = strstr (a, inx);
  if (pt == NULL) {
    return 0;
  }
  ZALLOC (str, len, char);
  off = pt - a;
  strncpy (str, a, off);
  strcpy (str + off, outx);
  x = strlen (outx);
  pt += strlen (inx);
  strcpy (str + off + x, pt);

  freestring (&a);
  *ap = strdup (str);
  free (str);
  return (1 + substringx (ap, inx, outx, niter + 1));
}

int
substring (char **ap, char *inx, char *outx)

/** 
 *ap is original string 
 all occurrences of inx are substituted with outx 
 can loop so be careful !!  
 

 NB.  ap must be on heap.  Fixed allocation not supported 
*/
{
  return (substringx (ap, inx, outx, 0));
}

int mapstrings(char **pstr, char **insub, char **outsub, int n)  
{
  char *ss ;  
  int k, t = 0 ; 
  ss = strdup(*pstr) ;
  
  for (k=0; k<n; ++k) {  
    t += substring(&ss, insub[k], outsub[k]) ; 
  }

  *pstr = strdup(ss) ; 
  freestring(&ss) ;
  return t ; 

}
  


int upstring (char *ss)

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

int
numcols (char *name)
// number of cols 
#define MAXCOLS 1000 
{
  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXSTR];
  char *sx;
  int nsplit, num = 0;

  if (name == NULL)
    fatalx ("(numlines)  no name");
  openit (name, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXCOLS);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    freeup (spt, nsplit);
    fclose (fff);
    return nsplit;
  }
  return -1 ;
}

int
numlines (char *name)
// number of lines   no comments or blanks
{
  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXSTR];
  char *sx;
  int nsplit, num = 0;

  num = 0;
  if (name == NULL)
    fatalx ("(numlines)  no name");
  openit (name, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, MAXFF);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    ++num;
    freeup (spt, nsplit);
  }
  fclose (fff);
  return num;
}

int
ftest (char *sss)
// can we open file for reading
{
  FILE *fdummy;

  fdummy = fopen (sss, "r");

  if (fdummy == NULL)
    return NO;
  fclose (fdummy);
  return YES;
}

void
openit (char *name, FILE ** fff, char *type)
{
  char *ss;
  if (name == NULL)
    fatalx ("\n(openit) null name\n");
  *fff = fopen (name, type);
  if (*fff == NULL) {
    ss = strerror (errno);
    printf ("bad open %s\n", name);
// system("lsof | fgrep np29") ;
    fatalx ("(openit) can't open file %s of type %s\n error info: %s\n", name, type,
            ss);
  }
}

void fcheckr(char *name)
// like ftest but just calls fatalx on errir
{
 FILE *fff ;

 openit(name, &fff,  "r") ;
 fclose(fff) ;


}

void fcheckw(char *name)
{
 FILE *fff ;

 openit(name, &fff,  "w") ;
 fclose(fff) ;


}


int
getxx (double **xx, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff;
  FILE *fff;
  int nbad = 0;

  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  maxff = MAX (MAXFF, numcol);

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < numcol) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      xx[i][num] = atof (spt[i]);
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

double
clocktime ()
{
  double xtime;
  double y;

  xtime = (double) clock ();
  y = xtime / (double) CLOCKS_PER_SEC;
  return y;
}

int
indxstring (char **namelist, int len, char *strid)
// look for string in list.  Was called indxindex
{
  int k;
  for (k = 0; k < len; k++) {
    if (namelist[k] == NULL)
      continue;
    if (strcmp (namelist[k], strid) == 0)
      return k;
  }
  return -1;
}

int
indxstringr (char **namelist, int len, char *strid)
// look for string in list.  Searches array in reverse ;
{
  int k;
  for (k = len - 1; k >= 0; k--) {
    if (namelist[k] == NULL)
      continue;
    if (strcmp (namelist[k], strid) == 0)
      return k;
  }
  return -1;
}

int
getnamesstripcolon (char ****pnames, int maxrow, int numcol, char *fname,
                    int lo, int hi)
{

// count is base 1
  char line[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff, numcolp, lcount = 0;
  FILE *fff;
  int nbad = 0;
  char ***names;

  names = *pnames;
  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  numcolp = numcol + 1;
  maxff = MAX (MAXFF, numcolp);

  while (fgets (line, MAXSTR, fff) != NULL) {
    subcolon (line);
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < numcol) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    ++lcount;
    if ((lcount < lo) || (lcount > hi)) {
      freeup (spt, nsplit);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      names[i][num] = strdup (spt[i]);
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

int
getnameslohi (char ****pnames, int maxrow, int numcol, char *fname, int lo,
              int hi)
{

// count is base 1
  char line[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff, numcolp, lcount = 0;
  FILE *fff;
  int nbad = 0;
  char ***names;

  names = *pnames;
  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  numcolp = numcol + 1;
  maxff = MAX (MAXFF, numcolp);

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < numcol) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    ++lcount;
    if ((lcount < lo) || (lcount > hi)) {
      freeup (spt, nsplit);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      names[i][num] = strdup (spt[i]);
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

int
getnames (char ****pnames, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff, numcolp;
  FILE *fff;
  int nbad = 0;
  char ***names;

  names = *pnames;
  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  numcolp = numcol + 1;
  maxff = MAX (MAXFF, numcolp);

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < numcol) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      names[i][num] = strdup (spt[i]);
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

int
getxxnames (char ***pnames, double **xx, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR];
  char sstt[MAXSTR] ;  
  char **spt ;
  char *sx;
  int nsplit, i, j, num = 0, maxff, numcolp;
  FILE *fff;
  int nbad = 0;
  char **names = NULL;

  if (pnames != NULL)
    names = *pnames;
  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  numcolp = numcol + 1;
  maxff = MAX (MAXFF, numcolp);
  ZALLOC (spt, maxff, char *) ;

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (names != NULL)
      names[num] = strdup (sx);
    if (nsplit < numcolp) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
// NA -> 0 
      strcpy(sstt, spt[i+1]) ; 
      if (strcmp(sstt, "NA") == 0) strcpy(sstt, "0") ;
      xx[i][num] = atof (sstt) ;      
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

int
getxxnamesf (char ***pnames, double **xx, int maxrow, int numcol, FILE * fff)

/** 
like getxxnames but file already open 
*/
{


  char line[MAXSTR];
  char sstt[MAXSTR] ; 
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff, numcolp;
  int nbad = 0;
  char **names;

  if (pnames != NULL)
    names = *pnames;

  numcolp = numcol + 1;
  maxff = MAX (MAXFF, numcolp);

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (names != NULL)
      names[num] = strdup (sx);
    if (nsplit < numcolp) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      strcpy(sstt, spt[i+1]) ; 
      if (strcmp(sstt, "NA") == 0) strcpy(sstt, "0") ;
      xx[i][num] = atof (sstt) ;      
    }
    freeup (spt, nsplit);
    ++num;
  }
  return num;
}


int
getss (char **ss, char *fname)

/** 
 get list of names 
*/
{

  char line[MAXSTR];
  char qqq[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff;
  FILE *fff;


  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  maxff = MAXFF;

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < 1) {
      continue;
    }
    ss[num] = strdup (spt[0]);
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

int
checkdup (char **list, int n)
{
  int a, b, t;
  for (a = 0; a < n; ++a) {
    for (b = a + 1; b < n; ++b) {
      t = strcmp (list[a], list[b]);
      if (t == 0)
        return YES;             // dup
    }
  }
  return NO;
}

void
printdups (char **list, int n)
{
  int a, b, t;
  for (a = 0; a < n; ++a) {
    for (b = a + 1; b < n; ++b) {
      t = strcmp (list[a], list[b]);
      if (t == 0)
        printf ("dup: %s\n", list[a]);
    }
  }
}


int
loadlist (char **list, char *listname)
// listname is just a list of names ... 
// dup check made
{
  int n, a, b, t;
  n = getss (list, listname);
  for (a = 0; a < n; ++a) {
    for (b = a + 1; b < n; ++b) {
      t = strcmp (list[a], list[b]);
      if (t == 0) {
        printf ("(loadlist) duplicate in list: %s\n", list[a]);
        fflush (stdout);
        fatalx ("(loadlist) duplicate in list: %s\n", list[a]);
      }
    }
  }
  return n;
}

char
revchar (char c)
{
  char cc;

  cc = toupper (c);
  if (cc == 'A')
    return 'T';
  if (cc == 'C')
    return 'G';
  if (cc == 'G')
    return 'C';
  if (cc == 'T')
    return 'A';

  return c;
}

void
crevcomp (char *sout, char *sin)
{
  char *sss, c, cout;
  int len;
  int i, j, t;

  len = strlen (sin);
  ZALLOC (sss, len + 1, char);
  sss[len] = CNULL;

  for (i = 0; i < len; ++i) {
    j = len - i - 1;
    c = sin[i];
    t = base2num (c);
    if (t < 0) {
      sss[j] = c;
      continue;
    }
    cout = num2base (3 - t);
    if (islower (c))
      cout = tolower (cout);
    sss[j] = cout;
  }
  strcpy (sout, sss);
  free (sss);
}

char *
int_string (int a, int len, int base)
{
  static char ss[100];
  int t = a, k, i;
  char *binary = "01";

  ss[len] = CNULL;
  for (i = 0; i < len; i++) {
    k = t % base;
    ss[len - i - 1] = '0' + k;
    t = t / base;
  }
  return ss;
// fragile
}

char *
binary_string (int a, int len)
{
  static char ss[100];
  int t = a, k, i;
  char *binary = "01";

  ss[len] = CNULL;
  for (i = 0; i < len; i++) {
    k = t % 2;
    ss[len - i - 1] = binary[k];
    t = t / 2;
  }
  return ss;
// fragile
}

char
num2iub (int num)
{

  char *iubstring = "ACGTMRWSYKVHDBX";
  char c;

  c = '?';
  if (num < 0)
    return c;
  if (num > 14)
    return c;

  return iubstring[num];

}

int
iub2num (char c)
{

  char *iubstring = "ACGTMRWSYKVHDBX";
  int t;
  char *sx;

  sx = strchr (iubstring, c);
  if (sx == NULL)
    return -1;
  return sx - iubstring;

}

char
num2base (int num)
{

  char *bases = "ACGT", c;
  c = '?';
  if (num < 0)
    return c;
  if (num > 3)
    return c;
  return bases[num];

}

int
base2num (char c)
{
  char cc;

  cc = toupper (c);

  switch (cc) {
  case 'A':
    return 0;
    break;
  case 'C':
    return 1;
    break;
  case 'G':
    return 2;
    break;
  case 'T':
    return 3;
    break;
  default:
    return -1;
  }
}

int
string_binary (char *sx)
{
  int *aa, len, i, t;
  char c;

  len = strlen (sx);
  ZALLOC (aa, len, int);

  for (i = 0; i < len; i++) {

    c = sx[i];
    if (c == '0')
      continue;
    if (c != '1')
      fatalx ("bad string: %s\n", sx);
    aa[i] = 1;
  }
  t = kodeitb (aa, len, 2);
  free (aa);
  return t;

}

void
freestring (char **ss)

/* note extra indirection */
{
  if (*ss == NULL)
    return;
  free (*ss);
  *ss = NULL;
}

void
copystrings (char **sa, char **sb, int n)
{
  int i;
  for (i = 0; i < n; ++i) {
    sb[i] = strdup (sa[i]);
  }
}

void
printstringsw (char **ss, int n, int slen, int width)
{
  int k, kmod;
  char fmt[10], s1[5];

  sprintf (s1, "%ds ", slen);
  strcpy (fmt, "%");
  strcat (fmt, s1);

  for (k = 0; k < n; ++k) {
    if (ss[k] != NULL)
      printf (fmt, ss[k]);
    else
      printf (fmt, "NULL");
    kmod = (k + 1) % width;
    if ((kmod == 0) && (k < (n - 1))) {
      printnl ();
    }
  }
  printnl ();
}

void
printstrings (char **ss, int n)
{
  int k;

  for (k = 0; k < n; ++k) {
    if (ss[k] != NULL)
      printf ("%s", ss[k]);
    else
      printf ("%s", "NULL");
    printnl ();
  }
}

int
ridfile (char *fname)
{
  int t;

  chmod (fname, 0777);
  t = unlink (fname);
  return t;
}

char
compbase (char x)
// upper case !!
// return complement
{
  if (x == 'A')
    return 'T';
  if (x == 'C')
    return 'G';
  if (x == 'G')
    return 'C';
  if (x == 'T')
    return 'A';

  return x;

}

void
mkupper (char *sx)
{
  int len, k;

  len = strlen (sx);
  for (k = 0; k < len; ++k) {
    sx[k] = toupper (sx[k]);
  }
}

void
mklower (char *sx)
{
  int len, k;

  len = strlen (sx);
  for (k = 0; k < len; ++k) {
    sx[k] = tolower (sx[k]);
  }
}


char *
strstrx (char *s1, char *s2)
// like strstr but case insensitive
// see also strcasestr
{
  char *ss1, *ss2, *spt;

  ss1 = strdup (s1);
  ss2 = strdup (s2);


  mkupper (ss1);
  mkupper (ss2);

  spt = strstr (ss1, ss2);
  if (spt != NULL) {
    spt = s1 + (spt - ss1);
  }

  freestring (&ss1);
  freestring (&ss2);

  return spt;

}


int
getjjnames (char ***pnames, int **jj, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff, numcolp;
  FILE *fff;
  int nbad = 0;
  char **names;

  names = *pnames;
  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  numcolp = numcol + 1;
  maxff = MAX (MAXFF, numcolp);

  while (fgets (line, MAXSTR, fff) != NULL) {
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    names[num] = strdup (sx);
    if (nsplit < numcolp) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      jj[i][num] = atoi (spt[i + 1]);
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}

int
isiub (char iub)
{
  char ss[5];
  int t;

  t = iubdekode (ss, iub);
  if (t == 0)
    return NO;
  return YES;

}

int
isiub2 (char iub)
{
// base or het iub should be upper case
  char ss[5];
  int t;

  t = iubdekode (ss, iub);
  if (t == 1)
    return YES;
  if (t == 2)
    return YES;
  return NO;

}

int
iubdekode (char *aa, char iub)
// a should be 5 long 
{

  char a[5];

  switch (iub) {

  case 'A':
    strcpy (a, "A");
    break;
  case 'C':
    strcpy (a, "C");
    break;
  case 'G':
    strcpy (a, "G");
    break;
  case 'T':
    strcpy (a, "T");
    break;
  case 'M':
    strcpy (a, "AC");
    break;
  case 'R':
    strcpy (a, "AG");
    break;
  case 'W':
    strcpy (a, "AT");
    break;
  case 'S':
    strcpy (a, "CG");
    break;
  case 'Y':
    strcpy (a, "CT");
    break;
  case 'K':
    strcpy (a, "GT");
    break;
  case 'V':
    strcpy (a, "ACG");
    break;
  case 'H':
    strcpy (a, "ACT");
    break;
  case 'D':
    strcpy (a, "AGT");
    break;
  case 'B':
    strcpy (a, "CGT");
    break;
  case 'X':
    strcpy (a, "ACGT");
    break;
  case 'N':
    strcpy (a, "ACGT");
    break;

  default:
    a[0] = CNULL;
  }
  if (aa != NULL)
    strcpy (aa, a);
  return strlen (a);
}

int
iubcbases (char *cbases, char iub)
// crack iub into 2 bases (which may agree) 
// return number of bases (1 or 2) or -1. 
{

  char uu[5];
  int nuu;

  nuu = iubdekode (uu, iub);

  if (nuu < 1)
    return -1;
  if (nuu > 2)
    return -1;
  if (nuu == 1)
    uu[1] = uu[0];

  cbases[0] = uu[0];
  cbases[1] = uu[1];

  return nuu;

}

int
ishet (char c)
{
  char aa[5];
  int n;

  n = iubdekode (aa, c);
  if (n == 2)
    return YES;
  return NO;
}

int
cttype (char cx)
// useful for CPG detection
{
  char cc;
  cc = toupper (cx);

  if (cc == 'A')
    return 2;
  if (cc == 'C')
    return 1;
  if (cc == 'G')
    return 2;
  if (cc == 'T')
    return 1;
  return -1;
}

char *
lastff (char *sss)
{
  char *sx;
  sx = strrchr (sss, '/');
  if (sx == NULL)
    return sss;
  return sx + 1;
}

int
char2int (char cc)
{

  int x;
  x = (int) (cc - '0');
  return x;

}

char
int2char (int x)
{

  char c;
  c = (char) ('0' + x);

  return c ; 
}

void
chomp (char *cc)
{
  int len;
  len = strlen (cc);

  if (len == 0)
    return;
  if (cc[len - 1] == CNL)
    cc[len - 1] = CNULL;

}

int
numcmatch (char *cc, int len, char c)
{
  int k, t = 0;

  for (k = 0; k < len; ++k) {
    if (cc[k] == c)
      ++t;
  }

  return t;

}

int
numcnomatch (char *cc, int len, char c)
{
  int k, t = 0;

  for (k = 0; k < len; ++k) {
    if (cc[k] != c)
      ++t;
  }

  return t;

}

char *
findupper (char *s)
// CNULL not tested
{
  char x;
  int len, k;

  len = strlen (s);

  for (k = 0; k < len; ++k) {
    if (isupper (s[k]))
      return s + k;
  }

  return NULL;

}

char *
strnotchar (char *s, char c)
// CNULL not tested
// return pointer to first char NOT c 
{
  char x;
  int len, k;

  len = strlen (s);

  for (k = 0; k < len; ++k) {
    if (s[k] != c)
      return s + k;
  }

  return NULL;

}

char
readtonl (FILE * fff)
{
  char c;

  for (;;) {
    c = fgetc (fff);
    if (c == EOF)
      break;
    if (c == CNL)
      break;
  }
  return c;
}

char *
fgetstrap (char *buff, int maxlen, FILE * fff, int *ret)
// fgets with long lines trapped
{
  int len;
  char c;

  if (fgets (buff, maxlen, fff) == NULL)
    return NULL;

  *ret = 1;
  len = strlen (buff);
  if (buff[len - 1] == CNL)
    return buff;

  *ret = 0;

  c = readtonl (fff);
  buff[0] = c;
  buff[1] = CNULL;
  return buff;

}
int filehash(char *name) 
{
#define MAXKL 256   
  FILE *fff;
  char line[MAXKL];
  int num = 0;
  int hash = 0, thash ; 

  num = 0;
  if (name == NULL)
    fatalx ("(filehash)  no name");
  openit (name, &fff, "r");
  while (fgets (line, MAXKL, fff) != NULL) {
    thash = stringhash(line) ; 
    hash += xhash1(thash ^ num) ;
    ++num ; 
  }
  fclose (fff);
  return abs(hash) ;
}

char *mytemp (char *qqq) 
// make temporary file name.   qqq is header string, eg "junk1" 
{
  char ss[MAXSTR] ; 
  int t ; 

  t = (int) getpid() ; 
  sprintf(ss, "/tmp/%s.%d", qqq, t) ; 
  return strdup(ss) ;
}

void printslurmenv () 
{
 char *ss ; 
 char sss[256] ;  

 ss = getenv("SLURM_JOBID") ; 
 if (ss==NULL) return ; 
 sprintf(sss, "sstat -j %s --format=jobid,avecpu,averss", ss) ;  
 fflush(stdout) ;
 printf("\n") ; 
 printf("====================================================================================\n") ;
 system(sss) ; 
 printf("\n") ; 
 fflush(stdout) ;

}

int getfline(char *ss, char *fname, int maxstr)
{
  FILE *fff ;
  int n ;
  char *pt ;
  openit(fname, &fff, "r") ;

  pt = fgets(ss, maxstr, fff) ;

  fclose(fff) ;
  if (pt == NULL) ss[0] = CNULL ;
  else {
   n = strlen(ss) ;
   ss[n-1] = CNULL ;
  }
  return strlen(ss) ;

}

int
numcolsq (char *name)
// number of cols : and , stripped 
#define MAXCOLS 1000 
{
  FILE *fff;
  char line[MAXSTR];
  char *spt[MAXSTR];
  char *sx;
  int nsplit, num = 0;

  if (name == NULL)
    fatalx ("(numlines)  no name");
  openit (name, &fff, "r");
  while (fgets (line, MAXSTR, fff) != NULL) {
    subcolon(line) ;
    nsplit = splitup (line, spt, MAXCOLS);
    if (nsplit == 0)
      continue;
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    freeup (spt, nsplit);
    fclose (fff);
    return nsplit;
  }
  return -1 ;
}

int
getxxq (double **xx, int maxrow, int numcol, char *fname)
// : and , stripped
{

  char line[MAXSTR];
  char *spt[MAXFF];
  char *sx;
  int nsplit, i, j, num = 0, maxff;
  FILE *fff;
  int nbad = 0;

  if (fname == NULL)
    fff = stdin;
  else {
    openit (fname, &fff, "r");
  }
  maxff = MAX (MAXFF, numcol);

  while (fgets (line, MAXSTR, fff) != NULL) {
    subcolon(line) ;
    nsplit = splitup (line, spt, maxff);
    if (nsplit == 0) {
      freeup (spt, nsplit);
      continue;
    }
    sx = spt[0];
    if (sx[0] == '#') {
      freeup (spt, nsplit);
      continue;
    }
    if (nsplit < numcol) {
      ++nbad;
      if (nbad < 10)
        printf ("+++ bad line: nsplit: %d numcol: %d\n%s\n", nsplit, numcol,
                line);
      continue;
    }
    if (num >= maxrow)
      fatalx ("too much data\n");
    for (i = 0; i < numcol; i++) {
      xx[i][num] = atof (spt[i]);
    }
    freeup (spt, nsplit);
    ++num;
  }
  if (fname != NULL)
    fclose (fff);
  return num;
}


int copyfs(char *infile, FILE *fff) 
// copy file  to stream
{
  char line[MAXSTR];
  int num = 0;
  FILE *ggg ; 

  openit(infile, &ggg, "r") ;
  while (fgets (line, MAXSTR, ggg) != NULL) {
    fprintf(fff, "%s", line) ; 
    ++num ; 
  }

  fclose(ggg) ;
  return num ;
}

