#include  <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdlib.h>
#include "strsubs.h"  
#include "sortit.h"  
#include "vsubs.h"  

/** 
 a simple sort routine
*/  

static double *ttt ;
static int *ittt ;
static int **pttt ;
static int plen=0 ;
static int *porder = NULL ;

void setorder(int *pp, int rlen) 
{
  int *tt ;

  if (plen > 0) { 
   if (porder != NULL) free(porder) ;
  }

   if (pp==NULL) { 
    porder = NULL ;
    plen = rlen ; 
    return ;
   }
   ZALLOC(porder, rlen, int) ;
   ZALLOC(tt, rlen, int) ;
   copyiarr(pp, tt, rlen) ;
   isortit(tt, porder, rlen) ;
   free(tt) ;
   plen = rlen ;
}
void sortit(double *a, int *ind, int len) 
{
  int i,k  ;
  int *inda ;

  if (len==0) fatalx("(sortit) len = 0\n") ;
  ZALLOC(ttt, len, double) ;
  ZALLOC(inda, len, int) ;

  for (i=0; i<len; i++) {
   inda[i] = i ;
  }

  copyarr(a,ttt,len) ;
  qsort((int *) inda, len, sizeof(int), (int (*) (const void *, const void *)) compit);

  for (i=0; i<len; i++) {
   k = inda[i] ;
   a[i] = ttt[k] ;
  }
  free (ttt) ;
  if (ind != NULL) copyiarr(inda, ind, len) ;
  free(inda) ;
}

int compit (int *a1, int *a2) 
{
 if (ttt[*a1] < ttt[*a2]) return -1 ;
 if (ttt[*a1] > ttt[*a2]) return 1 ;
 return 0 ;
}

void isortit(int *a, int *ind, int len) 
{
  int i,k  ;
  int *inda ;

  if (len==0) fatalx("(isortit) len = 0\n") ;
  ZALLOC(ittt, len, int) ;
  ZALLOC(inda, len, int) ;

  for (i=0; i<len; i++) {
   inda[i] = i ;
  }

  copyiarr(a,ittt,len) ;
  qsort((int *) inda, len, sizeof(int), (int (*) (const void *, const void *)) icompit);

  for (i=0; i<len; i++) {
   k = inda[i] ;
   a[i] = ittt[k] ;
  }
  free (ittt) ;
  if (ind != NULL) copyiarr(inda, ind, len) ;
  free(inda) ;
}

int icompit (int *a1, int *a2) 
{
 if (ittt[*a1] < ittt[*a2]) return -1 ;
 if (ittt[*a1] > ittt[*a2]) return 1 ;
 return 0 ;
}

void invperm(int *a, int *b, int n) {
/** 
 a, b can be same 
*/
     int i, j ;
     int *x ;

     if (n==0) return ;
     ZALLOC(x, n, int) ;

     ivclear(x,-1,n) ;
     for (i=0; i<n; i++)  {
      j=b[i] ;
      x[j]=i ;
     }
     copyiarr(x, a, n) ;
     free(x) ;
}

void ipsortit(int **a, int *ind, int len, int rlen) 
{
 ipsortitp(a, ind, len, rlen, NULL) ;

}
void ipsortitp(int **a, int *ind, int len, int rlen, int *order) 
/** 
 sort integer array pointers 
 rows of array are sorted in lexicographical order 

 compiarr can be called outside the sort
*/
{
  int i,k  ;
  int *inda ;

  if (len==0) fatalx("(ipsortit) len = 0\n") ;
  ZALLOC(pttt, len, int *) ;
  ZALLOC(inda, len, int) ;

  setorder(order, rlen) ; // order defines order as sorted in ascending order.  

  for (i=0; i<len; i++) {
   if (a[i] == NULL) fatalx("(ipsortit) array pointer %d NULL\n",i) ;
   inda[i] = i ;
  }

  copyiparr(a,pttt,len) ;
  qsort((int *) inda, len, 
   sizeof(int), (int (*) (const void *, const void *)) ipcompit);

  for (i=0; i<len; i++) {
   k = inda[i] ;
   a[i] = pttt[k] ;
  }
  if (ind != NULL) copyiarr(inda, ind, len) ;
  free(inda) ;
  free (pttt) ;
}

int ipcompit (int *a1, int *a2) 
{
 int l ;  
 l = compiarr(pttt[*a1], pttt[*a2], plen) ;
 return l ;
}

int compiarr(int *a, int *b, int len) 
{
      int i, k ;
      for (i=0; i<len; i++) {
       k = i ; if (porder != NULL) k = porder[i] ;
       if (a[k] < b[k]) return -1 ;
       if (a[k] > b[k]) return 1 ;
      }
      return 0 ;
}

