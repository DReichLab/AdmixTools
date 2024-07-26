#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <nicklib.h>
#include "xsearch.h"
#include "sortit.h"

static ENTRY *xentry;
static ENTRY **xxee;
static int xxeenum = -1;
static int xnum, xloaded;

static int debug = NO;

static int fancyhash = NO;



/* ********************************************************************* */


void
xhcreate (int n)
{
  int t, i;
  if (xentry != NULL) {
    free (xentry);
  }
  if (n == 0)
    fatalx ("(xhcreate) zero length\n");
  xnum = n;
  t = xnum % 17;
  if (t == 0)
    ++xnum;                     // for crude hash below
  ZALLOC (xentry, xnum, ENTRY);
  for (i = 0; i < xnum; i++) {
    xentry[i].key = NULL;
  }
  xloaded = 0;
}

void xhinit(int n) 
{
 int i ; 
 xhcreate(2*n) ; 
 ZALLOC(xxee, n, ENTRY *) ;
 for (i=0; i<n; ++i) { 
  ZALLOC (xxee[i], 1, ENTRY);
 }
 xxeenum = n ; 
}  

void
xhdestroy ()
{
  free (xentry);
  xentry = NULL;
  xnum = xloaded = 0;
}

ENTRY *
xhsearch (ENTRY item, ACTION act)
{

  ENTRY *itempt, *xtempt;
  int x;
  char *ccc;
  double yload;

  itempt = &item;
  ccc = itempt->key;
  x = xlookup (itempt->key, act);
  if (debug)
    printf ("lookup: %s %d\n", itempt->key, x);
  if ((x < 0) && (act == FIND))
    return NULL;
  if ((x >= 0) && (act == FIND))
    return xentry + x;
  if ((x < 0) && (act == ENTER))
    fatalx ("duplicate key %s\n", itempt->key);

  xtempt = xentry + x;
  xtempt->key = itempt->key;
  xtempt->data = itempt->data;
  ++xloaded;
  yload = (double) xloaded / (double) xnum;
  if (yload > 0.9)
    fatalx ("excessive xsearch load\n");
  return xtempt;
}

void xstore(char *key, int val) 
{
  ENTRY *pitem ;  
  char xx[10] ; 
  static int xxnum = 0 ; 
  
  pitem = xxee[xxnum] ;
  ++xxnum ;
  pitem->key = strdup(key) ;
  sprintf (xx, "%d", val);
  pitem -> data = strdup(xx) ;
  xhsearch (*pitem, ENTER);
  ++xloaded;
}

int
xlookup (char *key, ACTION act)
{
  ENTRY *xpt;
  int xbase, x, k;

  xbase = x = xhash (key);

  for (;;) {
    xpt = xentry + x;
    if (xpt->key == NULL) {
      if (act == FIND)
        return -1;
      return x;
    }
    k = strcmp (key, xpt->key);
    if (k == 0) {
      if (act == FIND)
        return x;
      return -1;
    }
    ++x;
    if (x >= xnum)
      x = 0;
  }
}

void dumpxs()  
{
 int x ; 
 ENTRY *xpt ; 

 for (x=0; x<xnum; ++x) { 
  xpt = xentry + x ; 
  if (xpt -> key == NULL) continue ;
  printf ("zzdumpxs: %8d %s\n", x, xpt -> key) ;
 }

}

int
xhash (char *key)
{
  int t;
  t = stringhash (key);
  return abs (t) % xnum;
}


int
stringhash (char *key)
{

#define MAXKEYLEN 512
  int xpack[MAXKEYLEN];
  int len, wlen, w;
  unsigned char t;
  int thash = 7;
  int jmax, jmin;
  int i, j;

  if (key == NULL)
    return 13;
  len = strlen (key);
  if (len == 0)
    return 17;
  if (len >= MAXKEYLEN)
    fatalx ("key too long\n");

  if (fancyhash==2) return fnv_hash(key) ; 

  wlen = (len - 1) / 4;
  ++wlen;

  thash += wlen ; 
  thash += key[0] ;  

  for (i = 0; i < wlen; ++i) {
    jmin = 4 * i;
    jmax = MIN (len - 1, jmin + 3);
    w = 0;
    for (j = jmin; j <= jmax; ++j) {
      t = (unsigned char) key[j];
      w = (w << 8) ^ t;
    }
    xpack[i] = xcshift (w, i);
  }
  if (debug)
    printf ("zz %s %x %x\n", key, w, xpack[0]);

  for (i = 0; i < wlen; i++) {
    thash += xhash1 (xpack[i]);
    thash = xcshift (thash, 3);
  }

  return thash ; 

}

void setfancyhash(int val)  

{
 fancyhash = val ;
}

int getfancyhash()  
{
 return fancyhash ; 
}

int
xhash1 (int ww)
{

  int k, w, w1, w2;
  w = xcshift (ww, 17);
  if (fancyhash == NO)
    return 17 * w;
  for (k = 0; k < 3; ++k) {
    w1 = w >> 16;
    w2 = w << 16;
    w = w2 ^ xhash2 (w1) ^ (w2 >> 16);
  }
  return w;
}

int
xhash2 (int x)
{

  int xmax = 65535;
  int t;

  if (x == 0)
    return xmax;
  if (x == xmax)
    return 0;

  t = x * 11;
  return t % xmax;

}

int
xcshift (int x, int shft)
{
  int a, b;

  if (shft == 0)
    return x;
  a = x << shft;
  b = x >> (32 - shft);

  return a ^ b;

}

void
xdestroy ()
{
  int i, num;
  ENTRY *pitem;

  if (xxee == NULL)
    return;
  num = xxeenum;
  for (i = 0; i < num; i++) {
    pitem = xxee[i];
    if (pitem == NULL)
      continue;
    free (pitem->key);
    free (pitem->data);
    free (pitem);
  }
  free (xxee);
  xhdestroy ();
}

int
xloadsearchx (char **ss, int n)
{

  ENTRY item, *pitem;
  char xx[10];
  int i, t;

  xhcreate (2 * n);
  ZALLOC (xxee, n, ENTRY *);
  xxeenum = n;
  for (i = 0; i < n; i++) {
    t = xlookup (ss[i], FIND);
    if (t >= 0)
      return i;
    ZALLOC (xxee[i], 1, ENTRY);
    pitem = xxee[i];
    pitem->key = strdup (ss[i]);
    sprintf (xx, "%d", i);
    pitem->data = strdup (xx);
    xhsearch (*pitem, ENTER);
  }
  return -1;
}

void
xloadsearch (char **ss, int n)
{

  ENTRY item, *pitem;
  char xx[10];
  int i;

  xhcreate (2 * n);
  ZALLOC (xxee, n, ENTRY *);
  xxeenum = n;
  for (i = 0; i < n; i++) {
    ZALLOC (xxee[i], 1, ENTRY);
    pitem = xxee[i];
    pitem->key = strdup (ss[i]);
    sprintf (xx, "%d", i);
    pitem->data = strdup (xx);
    xhsearch (*pitem, ENTER);
  }
}

int
xfindit (char *ss)
{

  ENTRY item, *pitem;
  int k;

  item.key = ss;
  pitem = xhsearch (item, FIND);
  if (pitem == NULL)
    return -1;
  sscanf (pitem->data, "%d", &k);
  return k;

}

int
finddup (char **ss, int n)
{

  int t;

  t = xloadsearchx (ss, n);
  xdestroy ();
  return t;

}
int
finddupalt (char **sss, int n)
{
// no hash table 

  int t, k, t1, t2, ttt ; 
  char **ss, *s1, *s2  ; 
  int *indx ; 
  
  ZALLOC(ss, n, char *) ;
  ZALLOC(indx, n, int) ;

  copystrings(sss, ss, n) ;
  sortstrings(ss, indx, n) ;  
  t = -1 ;
 // printf("qqq\n") ; printstrings(sss, n) ; printimat(indx, 1, n) ; printnl() ;

  for (k=1; k<n; ++k) { 
   t1 = indx[k-1] ; 
   t2 = indx[k] ; 
   s1 = sss[t1] ;
   s2 = sss[t2] ;
   ttt = strcmp(s1, s2)  ;
//  printf("zzz %d %d %d %s %s\n", t1, t2, ttt, s1, s2) ;
    if (ttt == 0) { 
    t = k ; 
    break ;
   }
  }
  free(indx) ; 
  freeup(ss, n) ; 
  free(ss) ;

  return t;

}


int xshash(int x) 
// slow but good 32 bit hash
{

 long w ; 
 int a, b ; 

 w = xlhash(x) ; 

 a = w >> 32 ; 
 return a ; 
 


}

long xlhash (long x) 
// slow but good hash.  
{
 static long bigp = 0 ;
 int k, nround = 7 ; 
 long a1, b1, a2, b2, w ; 
 static long *addpat ; 

 if (bigp == 0) {
   bigp = lpow2(61) - 1 ; 
   ZALLOC(addpat, nround+1, long) ; 
   for (k=1; k<=nround; ++k) { 
    w = lpow2(41) + 89 + k ; 
    addpat[k] = modinv(w, bigp) ; 
   }
 }

 w = x ^ 101010101 ; ; 
 for (k=1; k<=nround; ++k) { 
  a1 = w & bigp ; 
  b1 = w >> 61 ; 
  a2 = modinv(a1, bigp) ; 
  b2 = (a2 >> 7) & 7 ; 
  w  = a2 << 3 ; 
  b2 = b2 ^ b1 ; 
  w  = w ^ b2 ; 
  w ^= addpat[k] ; 
 }

 return w ; 
}

int fnv_hash(char *strng) 
{
// from fnv.c fnv.h on net.  Modified
    int h ; 
    char *p  ; 

    h = FNV_OFFSET_BASIS;

    p = strng ; 

    for (;;) { 
        if (*p == CNULL) break ; 
        h *= FNV_PRIME;
        h ^= *p;
        ++p ; 
    }
    return h ;
}

