#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "vsubs.h"
#include "strsubs.h"
#include "statsubs.h"
#include "linsubs.h"
static int linsolvx (int nDim, double *pfMatr, double *pfVect,
                     double *pfSolution);
static int calcdet = NO;

void
setcalcdet (int val)
{
  calcdet = val;
}

void
bal (double *a, double *b, int n)

/** 
 normalize mean 0 s.d 1 
 no error checking 
*/
{
  double t;
  t = asum (b, n) / (double) n;
  vsp (a, b, -t, n);

  t = asum2 (a, n) / (double) n;
  vst (a, a, 1.0 / sqrt (t), n);
}


void xmultx(double *a, double *b, int m, int n) 
{
 // M x M' : a is m x m 
 double *tb ; 

 ZALLOC(tb, m*n, double) ; 
 
 transpose(tb, b, m, n) ;  
 mulmat(a, b, tb, m, n, m) ; 
 

 free(tb) ;


}
void txmulx(double *a, double *b, int m, int n) 
{
 // M' x M : a is n x n 
 double *tb ; 

 ZALLOC(tb, m*n, double) ; 
 
 transpose(tb, b, m, n) ;  
 mulmat(a, tb, b, n, m, n) ; 
 

 free(tb) ;



}
void
mulmat (double *a, double *b, double *c, int a1, int a2, int a3)

/* b is a1 x a2 , c a2 x a3 so a is a1 x a3  */
{
  double *t;
  int i, j, k;
  ZALLOC (t, a1 * a3, double);

  for (i = 0; i < a1; i++)
    for (j = 0; j < a3; j++)
      for (k = 0; k < a2; k++)
        t[i * a3 + j] += b[i * a2 + k] * c[k * a3 + j];

  copyarr (t, a, a1 * a3);

  free (t);

}

void
imulmat (int *a, int *b, int *c, int a1, int a2, int a3)

/* b is a1 x a2 , c a2 x a3 so a is a1 x a3  */
{
  int *t;
  int i, j, k;
  ZALLOC (t, a1 * a3, int);

  for (i = 0; i < a1; i++)
    for (j = 0; j < a3; j++)
      for (k = 0; k < a2; k++)
        t[i * a3 + j] += b[i * a2 + k] * c[k * a3 + j];

  copyiarr (t, a, a1 * a3);

  free (t);

}

double logdet(double *mat, int n) 
{
  return pdinv(NULL, mat, n) ;
}


double
pdinv (double *cinv, double *coeff, int n)
// cinv and coeff can be same
// cinv can be NULL 
// return log det (coeff) 
{
  double *tt;
  double *p;
  double t, sum, y;
  int i, j, k;

/**
   pmat(coeff, n) ;
*/
  ZALLOC (tt, n * n, double);
  ZALLOC (p, n, double);


  copyarr (coeff, tt, n * n);

  choldc (tt, n, p);


  for (i = 0; i < n; i++) {
    tt[i * n + i] = 1.0 / p[i];
    for (j = i + 1; j < n; j++) {
      sum = 0.0;
      for (k = i; k < j; k++) {
        sum -= tt[j * n + k] * tt[k * n + i];
      }
      tt[j * n + i] = sum / p[j];

    }
  }

  for (i = 0; i < n; i++)
    for (j = i; j < n; j++) {
      sum = 0.0;
      if (cinv == NULL)
        break;
      for (k = j; k < n; k++) {
        sum += tt[k * n + j] * tt[k * n + i];
      }
      cinv[i * n + j] = cinv[j * n + i] = sum;
    }

  vlog (p, p, n);
  y = 2.0 * asum (p, n);


  free (tt);
  free (p);

  return y;

}

void
vzclear (double *mat, double *rhs, int dim, int varfix, double vval)
{

  int k, n = dim, v = varfix;

  for (k = 0; k < dim; ++k) {
    rhs[k] -= mat[k * n + v] * vval;    // was serious bug. this line was omitted
    mat[k * n + v] = mat[v * n + k] = 0;
  }

  mat[v * n + v] = 1;
  rhs[v] = vval;
}

int  solvitfix(double *prod, double *rhs, int n, double *ans, int *vfix, double *vvals, int nfix) 
{
  double *ww, *wa, *ww2, *lam, *cfix, *rfix, *qfix  ; 
  double *pmat ; 
  double *qco, *rco, *pp ;
  double y ;
  int k, j, t ; 

    t = solvit(prod, rhs, n, ans) ;
    if (t<0) return -2 ;  
    if (nfix==0) return t ;

  ZALLOC(pmat, n*n, double) ; 
  ZALLOC(ww, n, double) ;
  ZALLOC(wa, n, double) ;
  ZALLOC(ww2, n, double) ;
  ZALLOC(cfix, nfix*n, double) ;
  ZALLOC(qfix, nfix*n, double) ;
  ZALLOC(rfix, nfix, double) ;
  ZALLOC(lam, nfix, double) ;

  for (k=0; k<nfix; ++k) { 
    j = vfix[k] ;  
    qfix[k*n+j] = 1 ;   // qfix  is nfix * n 
    rfix[k] = vvals[k] ;  
  }
  pdinv(pmat, prod, n) ;
  mulmat(wa, rhs, pmat, 1, n, n) ;  // u 
  mulmat(cfix, qfix, pmat, nfix, n, n) ;  // d  

  ZALLOC(qco, nfix*nfix, double) ;
  ZALLOC(rco, nfix, double) ;

  for (k=0; k<nfix; ++k) {  
   pp = qfix+k*n ;    // c_k  
   rco[k] = -vdot(wa, pp, n) + rfix[k]; 
   for (j=0; j<nfix; ++j) { 
    qco[k*nfix+j] = vdot(pp, cfix+j*n, n) ;
   }
  }

   t = linsolv(nfix, qco, rco, lam) ;
   if (t<0) fatalx("(solvitfix)\n") ;

   mulmat(ww2, lam, qfix, 1, nfix, n) ; 
   vvp(ww, rhs, ww2, n) ; 

   mulmat(ans, pmat, ww, n, n, 1) ;

  for (k=0; k<nfix; ++k) { 
    j = vfix[k] ;  
    ans[j] = vvals[k] ; 
  }

  free(pmat) ; 
  free(wa) ; 
  free(ww) ; 
  free(ww2) ;
  free(cfix) ; 
  free(qfix) ;
  free(lam) ;

  return t ;

}

int
oldsolvitfix (double *prod, double *rhs, int n, double *ans, int *vfix,
           double *vvals, int nfix)
// force variables in vfix list to vvals) 
{
  //The coefficient matrix should be positive definite

  /*AT : changed this code to take in matrix in a linear array form */
  double *ttt;
  double *b;
  double *p;
  int i, k, t;
  int ret;
  double y, ymul;


  if (nfix == 0)
    return solvit (prod, rhs, n, ans);

  ZALLOC (ttt, n * n, double);
  ZALLOC (b, n, double);

  copyarr (prod, ttt, n * n);
  copyarr (rhs, b, n);

  y = trace (ttt, n) / (double) n;
  ymul = 1.0 / y;
  vst (ttt, ttt, ymul, n * n);
  vst (b, b, ymul, n);

  for (k = 0; k < nfix; ++k) {
    vzclear (ttt, b, n, vfix[k], vvals[k]);
  }

  ret = solvit (ttt, b, n, ans);


  free (ttt);
  free (b);

  return ret;
}

int
solvit (double *prod, double *rhs, int n, double *ans)
{
  //The coefficient matrix should be positive definite

  /*AT : changed this code to take in matrix in a linear array form */
  double *ttt;
  double *b;
  double *p;
  int i;
  int ret;

  ret = visnan(prod, n*n) ; if (ret==YES) return -4 ; 
  ret = visnan(rhs, n) ; if (ret==YES) return -3 ; 

  ZALLOC (ttt, n * n, double);
  ZALLOC (p, n, double);
  ZALLOC (b, n, double);

  copyarr (prod, ttt, n * n);
  copyarr (rhs, b, n);

  ret = choldc (ttt, n, p);
  if (ret < 0)
    return -1;                  // not pos def
  cholsl (ttt, n, p, b, ans);

  free (ttt);
  free (b);
  free (p);


  return 1;
}

void
cholsl (double *a, int n, double p[], double b[], double x[])

/** 
 Numerical Recipes.  Must change 
*/
{
  /*AT: Changing the code */

  int i, k;
  double sum;


  for (i = 0; i < n; i++) {
    sum = b[i];
    for (k = i - 1; k >= 0; k--)
      sum -= a[i * n + k] * x[k];
    x[i] = sum / p[i];
  }

  for (i = (n - 1); i >= 0; i--) {
    sum = x[i];
    for (k = i + 1; k < n; k++)
      sum -= a[k * n + i] * x[k];
    x[i] = sum / p[i];
  }


}

int
choldc (double *a, int n, double *p)
{
  int i, j, k;
  double sum;


  vzero (p, n);

  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      sum = a[i * n + j];
      for (k = i - 1; k >= 0; k--)
        sum -= a[i * n + k] * a[j * n + k];
      if (i == j) {

              /**                                       
              printf("zzchol %d %20.10f %9.3f\n",i, sum, a[i][i]) ; 
	      */
        if (sum <= 0.0) {
          return -1;            // not pos def
        }
        p[i] = sqrt (sum);
      }
      else {
        a[j * n + i] = sum / p[i];

      }
    }
  }

  return 1;

}
int isposdef (double *a, int n) 
{
  double *aa, *p ; 
  int ret, k ; 
 
/** 
 quick check on diagonal
*/  
  for (k=0; k<n; ++k) { 
   if (a[k*n+k] <= 0.0) return NO ;
  }

  ZALLOC(aa, n*n, double) ; 
  ZALLOC(p, n, double) ; 

  copyarr(a, aa, n*n) ;
  ret = choldc(aa, n, p) ;


 free(aa) ; 
 free(p) ;

 if (ret>0) return YES ; 
 return NO ; 

}

void
pmat (double *mat, int n)

/** 
 print square matrix 
*/
{
  int i, j;
  double *diag;

  ZALLOC (diag, n, double);
  getdiag (diag, mat, n);

  printf ("pmat:\n");

  for (i = 0; i < n; i++) {
    printf ("diag %5d %9.3f\n", i, diag[i]);
    for (j = 0; j < n; j++) {
      if ((j % 10) == 9)
        printf ("\n");
      if ((n % 10) != 0)
        printf ("%9.3f ", mat[i * n + j]);
    }
    printf ("\n");
  }
  printf ("\n");
  printf ("\n");

  free (diag);
}

void
cholesky (double *cf, double *a, int n)
{
  int i, j, k;
  double *tt;
  double *p;

  ZALLOC (tt, n * n, double);
  ZALLOC (p, n, double);
  copyarr (a, tt, n * n);


  choldc (tt, n, p);

  vzero (cf, n * n);

  for (i = 0; i < n; i++) {
    tt[i * n + i] = p[i];

    for (j = 0; j <= i; j++) {
      k = i * n + j;
      cf[k] = tt[i * n + j];
    }
  }

  free (tt);
  free (p);
}

//==============================================================================
//
// Linear System Solution by Gauss method
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================

int
linsolv (int n, double *pfMatr, double *pfVect, double *sol)
// 1 on failure 
{

  int ret;
  double *a, *rhs;

  ZALLOC (a, n * n, double);
  ZALLOC (rhs, n, double);

  copyarr (pfMatr, a, n * n);
  copyarr (pfVect, rhs, n);

  ret = linsolvx (n, a, rhs, sol);

  free (a);
  free (rhs);


  return ret;


}

int
linsolvx (int nDim, double *pfMatr, double *pfVect, double *pfSolution)
{
  double fMaxElem;
  double fAcc;

  int i, j, k, m;


  for (k = 0; k < (nDim - 1); k++)      // base row of matrix
  {
    // search of line with max element
    fMaxElem = fabs (pfMatr[k * nDim + k]);
    m = k;
    for (i = k + 1; i < nDim; i++) {
      if (fMaxElem < fabs (pfMatr[i * nDim + k])) {
        fMaxElem = pfMatr[i * nDim + k];
        m = i;
      }
    }

    // permutation of base line (index k) and max element line(index m)
    if (m != k) {
      for (i = k; i < nDim; i++) {
        fAcc = pfMatr[k * nDim + i];
        pfMatr[k * nDim + i] = pfMatr[m * nDim + i];
        pfMatr[m * nDim + i] = fAcc;
      }
      fAcc = pfVect[k];
      pfVect[k] = pfVect[m];
      pfVect[m] = fAcc;
    }

    if (pfMatr[k * nDim + k] == 0.)
      return 1;                 // needs improvement !!!

    // triangulation of matrix with coefficients
    for (j = (k + 1); j < nDim; j++)    // current row of matrix
    {
      fAcc = -pfMatr[j * nDim + k] / pfMatr[k * nDim + k];
      for (i = k; i < nDim; i++) {
        pfMatr[j * nDim + i] =
          pfMatr[j * nDim + i] + fAcc * pfMatr[k * nDim + i];
      }
      pfVect[j] = pfVect[j] + fAcc * pfVect[k]; // free member recalculation
    }
  }

  for (k = (nDim - 1); k >= 0; k--) {
    pfSolution[k] = pfVect[k];
    for (i = (k + 1); i < nDim; i++) {
      pfSolution[k] -= (pfMatr[k * nDim + i] * pfSolution[i]);
    }
    pfSolution[k] = pfSolution[k] / pfMatr[k * nDim + k];
  }

  return 0;
}

double
mquad (double y0, double y1, double y2, double *pmx)
{

// max quadratic with max in interval (-1, 1) 3. pt interpolation
// max value return *pmx = arg max 

  double a, b, c;
  double mx;

/**
  a - b + c = y0 ;
  c  =  y1 ;
  a + b + c =  y2;
*/
  c = y1;
  b = (y2 - y0) / 2.0;
  a = y2 - (b + c);

  if (a == 0.0) {
    *pmx = mx = 0;
  }
  else
    *pmx = mx = -b / (2 * a);
  return a * mx * mx + b * mx + c;

}

void
qgrad (double *grad, double *vv, double *q, double *l, int n)
{
  double y;
  double *w1;

  mulmat (grad, vv, q, 1, n, n);
  vst (grad, grad, 2, n);
  vvp (grad, grad, l, n);

}

double
qval (double *vv, double *q, double *l, int n)
{
  double y;
  y = scx (q, NULL, vv, n);
  y += vdot (vv, l, n);
  return y;
}

double
seekzz (double *vv, double *vtarget, double *vsource, int *dead, int *constraint, int n)
// vsource feasible vtarget maybe ;  constraint = 1 => positive
{

  double *vmove, fmove, y;
  int k, kfix = -1;


  ZALLOC (vmove, n, double);
  vvm (vmove, vtarget, vsource, n);

  fmove = 1;

  for (k = 0; k < n; ++k) {
    if ((vtarget[k] < 0) && (dead[k] == 1) && (constraint[k] == 1)) fatalx ("(seekz) bug\n");
    if (dead[k] == 1) continue;
    if (vtarget[k] >= 0) continue;
    if (constraint[k] == 0) continue ;
    y = -vsource[k] / vmove[k];
    if (y < fmove) {
      fmove = y;
      kfix = k;
    }
  }
  if (kfix < 0) {
    copyarr (vtarget, vv, n);
    fmove = 1;
  }
  else {
    vst (vmove, vmove, fmove, n);
    vvp (vv, vsource, vmove, n);
    vv[kfix] = 0;
  }

  free (vmove);
  return fmove;
}

double
seekz (double *vv, double *vtarget, double *vsource, int *dead, int n)
// vsource feasible vtarget maybe 
{

  double *vmove, fmove, y;
  int k, kfix = -1;


  ZALLOC (vmove, n, double);
  vvm (vmove, vtarget, vsource, n);

  fmove = 1;

  for (k = 0; k < n; ++k) {
    if ((vtarget[k] < 0) && (dead[k] == 1)) fatalx ("(seekz) bug\n");
    if (dead[k] == 1)
      continue;
    if (vtarget[k] >= 0)
      continue;
    y = -vsource[k] / vmove[k];
    if (y < fmove) {
      fmove = y;
      kfix = k;
    }
  }
  if (kfix < 0) {
    copyarr (vtarget, vv, n);
    fmove = 1;
  }
  else {
    vst (vmove, vmove, fmove, n);
    vvp (vv, vsource, vmove, n);
    vv[kfix] = 0;
  }

  free (vmove);
  return fmove;
}

double qmp (double *vnew, double *vold, double *q, double *l, int *dead, int level, int n)
{
  int *constraint ; 
  double y ;

   ZALLOC(constraint, n, int) ;
   ivclear(constraint, 1, n) ;

   y = qmpc (vnew, vold, q, l, dead, level, constraint, n) ;

   free(constraint) ; 
   return y ; 

}


double
qmin (double *vv, double *q, double *l, int n)
{
  double *rhs;
  int t;

  ZALLOC (rhs, n, double);
  vst (rhs, l, -0.5, n);
  t = solvit (q, rhs, n, vv);

  if (t < 0) {
    printmatl (q, n, n);
    printmatl (l, 1, n);
    printmatl (vv, 1, n);
    fatalx ("bad qmin\n");
  }

  free (rhs);
  return qval (vv, q, l, n);
}

double
qminfix (double *vv, double *q, double *l, int n, int *fixlist,
         double *fixvals, int nfix)
{

  double *rhs;
  int t;

  if (nfix == 0)
    return qmin (vv, q, l, n);

  ZALLOC (rhs, n, double);
  vst (rhs, l, -0.5, n);
  t = oldsolvitfix (q, rhs, n, vv, fixlist, fixvals, nfix);

  if (t < 0) {
    printmatl (q, n, n);
    printmatl (l, 1, n);
    printmatl (vv, 1, n);
    fatalx ("bad qminfix\n");
  }

  free (rhs);
  return qval (vv, q, l, n);



}

double
qminposc (double *vv, double *q, double *l, int *constraint, int n)
// q is pos def matrix.  Find min v q v' + l.v ;  vv is minimum.  returns minimum value
// q is scaled sensibly to improve numerics.   Method constrained hill climb
// constraint[k] = 1 => pos constraint
{
  double *qq, *ll;
  double *best, *w1;
  int *dead;

  double *scale, *s2, y;
  int k;

  ZALLOC (qq, n * n, double);
  ZALLOC (ll, n, double);
  ZALLOC (w1, n, double);

  ZALLOC (scale, n, double);
  ZALLOC (s2, n * n, double);

  ZALLOC (best, n, double);
  ZALLOC (dead, n, int);

  vclear (best, 0, n);          // start at boundary qmpc will move away on positive gradient
  getdiag (w1, q, n);
// crude check for pos def. 
  for (k = 0; k < n; ++k) {
    y = w1[k];
    if (y > 0)
      continue;
    if (y < 0)
      fatalx ("(qminpos) negative diag element\n");
    if (l[k] != 0)
      fatalx ("(minpos) zero diag\n");
// k is silly variable
    w1[k] = 1;
    dead[k] = 1;
    best[k] = 0;
  }

  vsqrt (scale, w1, n);
  addouter (s2, scale, n);

  vvd (qq, q, s2, n * n);
  vvd (ll, l, scale, n);

  qmpc (best, best, qq, ll, dead, 0, constraint, n);
  vvd (vv, best, scale, n);


  free (qq);
  free (ll);
  free (best);
  free (dead);
  free (scale);
  free (s2);
  free (w1);



  return qval (vv, q, l, n);

}

double
qminpos (double *vv, double *q, double *l, int n)
// q is pos def matrix.  Find min v q v' + l.v ;  vv is minimum.  returns minimum value
// q is scaled sensibly to improve numerics.   Method constrained hill climb
{

   int *constraint ;
   double y ;

   ZALLOC(constraint, n, int) ;
   ivclear(constraint, 1, n) ;

   y = qminposc (vv,  q, l, constraint, n) ;

   free(constraint) ; 

   return y ; 

}

double
qminposfixc (double *vv, double *q, double *l, int n, int *fixlist,
            double *fixvals, int nfix, int *constraint)
{

  double *wfix, *w1, *w2, *qq, *ll;
  int i, j, k;
  double y;

  if (nfix == n) {
    for (j = 0; j < nfix; ++j) {
      k = fixlist[j];
      vv[k] = fixvals[j];
    }
    return qval (vv, q, l, n);
  }

  if (nfix == 0) {
    return qminposc (vv, q, l, constraint, n);
  }

  ZALLOC (wfix, n, double);
  ZALLOC (w1, n, double);
  ZALLOC (w2, n, double);
  ZALLOC (qq, n * n, double);
  ZALLOC (ll, n, double);


  copyarr (q, qq, n * n);
  copyarr (l, ll, n);

  vclear (w2, 1, n);

  for (j = 0; j < nfix; ++j) {
    k = fixlist[j];
    wfix[k] = fixvals[j];
    for (i = 0; i < n; ++i) {
      qq[i * n + k] = qq[k * n + i] = 0;
      ll[i] += 2 * q[i * n + k] * wfix[k];
    }
  }
  for (j = 0; j < nfix; ++j) {
    k = fixlist[j];
    qq[k * n + k] = 1;
    ll[k] = -2 * wfix[k];
  }

  qminposc (vv, qq, ll, constraint, n);

  free (wfix);
  free (w1);
  free (w2);
  free (qq);
  free (ll);

  return qval (vv, q, l, n);

}
double
qminposfix (double *vv, double *q, double *l, int n, int *fixlist,
            double *fixvals, int nfix)
{

  double *wfix, *w1, *w2, *qq, *ll;
  int i, j, k;
  double y;
  int *constraint ; 

  if (nfix == n) {
    for (j = 0; j < nfix; ++j) {
      k = fixlist[j];
      vv[k] = fixvals[j];
    }
    return qval (vv, q, l, n);
  }

  ZALLOC(constraint, n, int) ; 
  ivclear(constraint, 1, n) ;

  y = qminposfixc(vv, q, l, n, fixlist, fixvals, nfix, constraint) ;

  free(constraint) ; 
  return y ;  

}

double qmpc (double *vnew, double *vold, double *q, double *l, int *dead, int level,  int *constraint, int n)
// constraint[] == 1 => positive
{
  int *vfix, nfix, k;
  int *dd;
  double *vvals, *v2, *xpt, *bpt, *grad, *vv, *rhs, *oldgrad;
  double y, y1, y2;
  int t, olddead, newdead, stuck = 0;
  double pp[3];

  ZALLOC (vvals, n, double);
  ZALLOC (v2, n, double);
  ZALLOC (xpt, n, double);
  ZALLOC (bpt, n, double);
  ZALLOC (vfix, n, int);
  ZALLOC (dd, n, int);
  ZALLOC (vv, n, double);
  ZALLOC (grad, n, double);
  ZALLOC (oldgrad, n, double);

  ZALLOC (rhs, n, double);

  copyarr (vold, vv, n);
  vst (rhs, l, -0.5, n);

  if (level > 200) {
    fatalx (" (qmpc) looping\n");
  }

  for (;;) {

    for (k = 0; k < n; ++k) {
      if (dead[k])
        vv[k] = 0;
    }
    qgrad (grad, vv, q, l, n);
    for (k = 0; k < n; ++k) {
      if (dead[k])
        grad[k] = 0;
    }

    copyiarr (dead, dd, n);
    nfix = 0;
    for (k = 0; k < n; ++k) {
      if ((vv[k] <= 0) && (constraint[k]==1)) vv[k] = 0;
      if ((vv[k] == 0) && (grad[k] >= 0)&& (constraint[k]==1)) {
        vfix[nfix] = k;
        ++nfix;
        dd[k] = 1;
      }
    }

    y1 = qval (vv, q, l, n);

    if (level<100) { 
     t = oldsolvitfix (q, rhs, n, xpt, vfix, vvals, nfix);
    }
    else { 
     t = oldsolvitfix (q, rhs, n, xpt, vfix, vvals, nfix);
    }
    if (t<0) fatalx("(qmpc) solvitfix failure kode: %d\n", t) ;

    seekzz (bpt, xpt, vv, dd, constraint, n);
    y2 = qval (bpt, q, l, n);

 
 if (level>100) {
  printf("zzloop  %d %d %12.6f %12.6f\n", intsum(dd, n), level, y1, y2) ;
  printf("zz1: ") ; printmatw(vv, 1, n, 10) ;
  printf("zz2: ") ; printmatw(xpt, 1, n, 10) ;
  printimatw(dd, 1, n, 10) ;
  printf("zz3: ") ; printmatw(bpt, 1, n, 10) ;
  printf("zz4: ") ; printmatw(grad, 1, n, 10) ;
  fflush(stdout) ;
 }

    if (y2 < y1) {
      copyarr (bpt, vv, n);
      continue;
    }

    copyiarr (dead, dd, n);
    qgrad (grad, bpt, q, l, n);
    copyarr (grad, oldgrad, n);
    for (k = 0; k < n; k++) {
      if (dead[k] == 1) {
        bpt[k] = 0;
        grad[k] = oldgrad[k] = 0;
      }

      if ((bpt[k] <= 0) && (grad[k] >= 0) && (constraint[k]==1)) {
        grad[k] = oldgrad[k] = 0;
      }
      if ((xpt[k] <= 0.0) && (constraint[k]==1)) {
        grad[k] = 0;
        dd[k] = 1;
      }
    }

    y = asum2 (grad, n) / (double) n;
    if (y > 1.0e-8)
      y2 = qmpc (vv, bpt, q, l, dd, level + 1, constraint, n);
    if (y2 < y1) {
      continue;
    }

// now are there coords at boundary 

    y = asum2 (oldgrad, n) / (double) n;
    if (y < 1.0e-8)
      break;

    vvm (v2, bpt, oldgrad, n);
    pp[0] = -qval (v2, q, l, n);
    pp[1] = -qval (bpt, q, l, n);
    vvp (v2, bpt, oldgrad, n);
    pp[2] = -qval (v2, q, l, n);
    mquad (pp[0], pp[1], pp[2], &y);
//  printf ("zzmquad: %9.3f\n", y);
    vst (v2, oldgrad, y, n);
    vvp (xpt, v2, bpt, n);
    seekzz (bpt, v2, xpt, dead, constraint, n) ;
    y2 = qval (bpt, q, l, n);
    if (y2 < y1) {
      copyarr (bpt, vv, n);
      continue;
    }
    break;

  }

  copyarr (vv, vnew, n);

  free (vvals);
  free (v2);
  free (xpt);
  free (bpt);
  free (vfix);
  free (grad);
  free (dd);
  free (rhs);
  free (vv);

  return y1;

}

