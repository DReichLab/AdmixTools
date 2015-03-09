#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include "xsearch.h" 

static ENTRY *xentry ;
static ENTRY **xxee  ;
static int xxeenum = -1 ;
static int xnum, xloaded  ;

static int debug = NO ;

static int fancyhash = NO ;


/* ********************************************************************* */

void xhcreate (int n) 
{  
 int t, i ;
 if (xentry != NULL) { 
   free(xentry) ;
 }
 if (n==0) fatalx("(xhcreate) zero length\n") ;
 xnum = n ;
 t = xnum % 17 ; 
 if (t==0) ++xnum ; // for crude hash below
 ZALLOC (xentry, xnum , ENTRY) ;
 for (i=0; i<xnum; i++) {  
  xentry[i].key = NULL ;
 }
 xloaded = 0 ;
}

void xhdestroy() 
{ 
 free(xentry) ;
 xentry = NULL ;
 xnum = xloaded = 0 ;
}

ENTRY *xhsearch(ENTRY item, ACTION act) 
{
   
  ENTRY *itempt, *xtempt ;
  int x ;
  char *ccc ;
  double yload ;

  itempt = &item ;
  ccc = itempt -> key ;
  x = xlookup(itempt -> key, act) ;
  if (debug)
   printf("lookup: %s %d\n", itempt -> key, x) ;
  if ((x < 0)  && (act == FIND)) return NULL ;
  if ((x>=0) && (act == FIND)) return xentry + x ;
  if ((x<0) && (act == ENTER)) fatalx("duplicate key %s\n", itempt ->key) ;

  xtempt = xentry + x ; 
  xtempt -> key = itempt -> key ;  
  xtempt -> data = itempt -> data ;
  ++xloaded ;  
  yload = (double) xloaded / (double) xnum ; 
  if (yload>0.9) fatalx("excessive xsearch load\n") ;
  return xtempt ;
}

int xlookup(char *key, ACTION act) 
{
  ENTRY *xpt ;
  int xbase, x, k  ; 

  xbase = x = xhash(key) ;  
  for (;;)  {  
   xpt = xentry + x ;
   if (xpt -> key == NULL) {  
    if (act == FIND) return -1 ;
    return x ;
   }
   k = strcmp(key, xpt -> key) ;
   if (k==0) {
    if (act == FIND) return x  ;
    return -1 ;
   }
   ++x ; 
   if (x>=xnum) x=0 ;
  }
}
int xhash (char *key) 

{
 int t ; 
 t = stringhash(key) ;
 return abs(t) % xnum ;
}

int stringhash(char *key) 
{ 

#define MAXKEYLEN 512
 int xpack[MAXKEYLEN] ;
 int len, wlen, w ;
 unsigned char t  ;
 int thash = 7 ;
 int jmax, jmin ;
 int i, j ; 

 if (key == NULL) return 13 ;
 len = strlen(key) ; 
 if (len ==0) return 17  ;
 if (len >= MAXKEYLEN) fatalx("key too long\n") ;

 wlen = (len-1)/4 ; 
 ++wlen ; 
 
 for (i=0; i<wlen; ++i)  {  
  jmin = 4*i ; 
  jmax = MIN(len-1, jmin+3) ;
  w = 0 ;
  for (j= jmin; j <= jmax ; ++j) {  
   t = (unsigned char) key[j] ;
   w = (w << 8) ^ t ;      
  }
  xpack[i] = xcshift(w, i) ;
 }
 if (debug)
   printf("zz %s %x %x\n", key, w, xpack[0]) ;
 for (i=0; i<wlen; i++) {  
  thash += xhash1(xpack[i]) ;
  if (debug)
   printf("zz2  %x\n", thash) ;
  thash = xcshift(thash, 3) ;
 }
 if (debug)
  printf("key: %s  hash: %x\n", key, thash) ;

 return thash ;


}

int xhash1(int ww) 

{ 
  
  int k, w, w1, w2 ;
  w = xcshift(ww, 17) ;  
  if (fancyhash == NO) return 17*w ; 
  for (k=0; k<3; ++k) {
   w1 = w >> 16 ; 
   w2 = w << 16 ;  
   w = w2 ^ xhash2(w1) ^ (w2 >> 16);
  } 
  return w ;
}

int xhash2 (int x)  
{

 int xmax = 65535 ; 
 int t ;

 if (x==0) return xmax ;  
 if (x==xmax) return 0 ;

 t = x * 11 ; 
 return t % xmax ;

}

int xcshift(int x, int shft) 
{
 int a, b  ;

 if (shft==0) return x ;
 a = x << shft ;
 b = x >> (32 - shft) ;

 return a ^ b ;

}

void xdestroy()
{
 int i, num ;
 ENTRY *pitem ;

 if (xxee == NULL) return ;
 num = xxeenum ;
 for (i=0; i<num; i++) {  
  pitem = xxee[i] ;
  if (pitem == NULL) continue ; 
  free(pitem -> key) ;
  free(pitem -> data) ;
  free(pitem) ;
 }
 free(xxee) ;
 xhdestroy() ;
}

int xloadsearchx(char **ss, int n)  

{  

 ENTRY item, *pitem ; 
 char  xx[8] ;
 int i, t ;

 xhcreate(2*n) ;
 ZALLOC(xxee, n, ENTRY *) ;
 xxeenum = n ;
 for (i=0; i<n; i++) {  
  t = xlookup(ss[i], FIND) ;
  if (t>=0) return i ;
  ZALLOC(xxee[i], 1, ENTRY) ;
  pitem = xxee[i] ;
  pitem -> key = strdup(ss[i]) ;
  sprintf(xx, "%d", i) ;
  pitem -> data = strdup(xx) ; 
  xhsearch(*pitem, ENTER) ;
 }
 return -1 ;
}

void xloadsearch(char **ss, int n)  

{  

 ENTRY item, *pitem ; 
 char  xx[8] ;
 int i ;

 xhcreate(2*n) ;
 ZALLOC(xxee, n, ENTRY *) ;
 xxeenum = n ;
 for (i=0; i<n; i++) {  
  ZALLOC(xxee[i], 1, ENTRY) ;
  pitem = xxee[i] ;
  pitem -> key = strdup(ss[i]) ;
  sprintf(xx, "%d", i) ;
  pitem -> data = strdup(xx) ; 
  xhsearch(*pitem, ENTER) ;
 }
}

int xfindit(char *ss)    
{

  ENTRY item, *pitem ; 
  int k ;
 
  item.key = ss ;
  pitem = xhsearch(item, FIND) ; 
  if (pitem == NULL) return -1 ;
  sscanf(pitem -> data, "%d", &k) ;
  return k ;

}
 
int finddup(char **ss, int n) 
{

  int t ;

  t = xloadsearchx(ss, n) ;
  xdestroy() ;
  return t ;

}
