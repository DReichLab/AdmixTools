#include <nag.h>
#include <nag_stdlib.h>
#include <nagg01.h>
#include <nagf02.h>
#include <nagf04.h>

#include "nagsubs.h"
#include <nicklib.h>

static int debug = NO;

double
chisqtail (double x, double df)
{
  static NagError fail;
  double q;
  q = nag_prob_chi_sq (Nag_UpperTail, x, df, &fail);
  return q;
}

#define NMAX 1100

void
mleg (double a1, double a2, double *p, double *lam)
{
  int iter;
  double s, pp, ll;
  double top, bot, fval;

/** 
 solve 
 p/lam = a1 ; psi(p) - log(lam) = a2 ;
 Thus psi(p) - log(p) = a2 - log(a1) 
*/
  s = a2 - log (a1);

  if (s >= 0.0)
    fatalx ("log E(x) < E(log (x)) \n");
  pp = -s;

  for (iter = 1; iter <= 30; ++iter) {
    fval = s - (psi (pp) - log (pp));
    if (debug)
      printf ("yy1 %3d %9.3f %9.3f\n", iter, pp, fval);
    if (fval < 0.0)
      break;
    pp *= 2.0;
  }

  for (iter = 1; iter <= 30; ++iter) {
    fval = s - (psi (pp) - log (pp));
    if (fval > 0.0)
      break;
    if (debug)
      printf ("yy2 %3d %9.3f %9.3f\n", iter, pp, fval);
    pp /= 2.0;
  }

  for (iter = 1; iter <= 10; ++iter) {
    fval = psi (pp) - log (pp);
    top = s - fval;
    bot = tau (pp) - (1.0 / pp);
    if (debug)
      printf ("%3d %12.6f %12.6f\n", iter, pp, top);
    pp += top / bot;
  }
  ll = pp / a1;
  *p = pp;
  *lam = ll;
}

double
psi (double x)
{
  static NagError fail;
  Integer k;
  double q;
  double nag_real_polygamma (double x, Integer k, NagError * fail);

  k = 0;
  q = nag_real_polygamma (x, k, &fail);
  return q;
}

double
tau (double x)
{
  static NagError fail;
  Integer k;
  double q;
  double nag_real_polygamma (double x, Integer k, NagError * fail);

  k = 1;
  q = nag_real_polygamma (x, k, &fail);
  return q;
}

void
eigxsym (double *evec, double *eval, double *mat, int *nval, int n)
{

/** 
 returns eigenvalues, eigenvectors of 
 symmetric matrix mat.  

  real eigenvectors returned.
*/

#define TDA NMAX
#define TDV NMAX
  double amat[NMAX][TDV];
  double r[NMAX];
  double v[NMAX][TDV];
  Integer iter[NMAX];

  static NagError fail;
  int j, k, l, nv;
  int *ind;


  if (n > NMAX)
    fatalx ("(eigxsym) n too large\n");

  for (k = 0; k < n; k++) {
    copyarr (mat + k * n, amat[k], n);
  }

  f02abc (n, (double *) amat, (Integer) TDA, r,
          (double *) v, (Integer) TDV, &fail);

  copyarr (r, eval, n);

  *nval = nv = n;

  ZALLOC (ind, nv, int);
  sortit (eval, ind, nv);

  if (evec != NULL) {
    for (k = 0; k < nv; ++k) {
      j = ind[k];
      for (l = 0; l < n; ++l) {
        if (evec != NULL)
          evec[k * n + l] = v[l][j];
      }
    }
  }
  free (ind);
}


void
eigx (double *evec, double *eval, double *mat, int *nval, int n)
{

/** 
 returns eigenvalues, eigenvectors of 
 unsymmetric matrix mat.  

 Only real eigenvectors returned.
*/

#define TDA NMAX
#define TDV NMAX
  double amat[NMAX][TDA];
  Complex r[NMAX];
  Complex v[NMAX][TDV];
  Integer iter[NMAX];

  static NagError fail;
  int j, k, l, nv;
  int *ind;

#define COMPLEX(A)  A.re, A.im
#define REAL(A)  A.re
#define IMAG(A)  A.im

  if (n > NMAX)
    fatalx ("(eigx) n too large\n");

  printmat (mat, n, n);

  for (k = 0; k < n; k++) {
    copyarr (mat + k * n, amat[k], n);
  }

  f02agc (n, (double *) amat, (Integer) TDA, r,
          (Complex *) v, (Integer) TDV, iter, &fail);
  j = 0;

  for (k = 0; k < n; k++) {
    printf ("iter: %d\n", iter[k]);
    printf ("%9.3f %9.3f\n", COMPLEX (r[k]));
    if (fabs (IMAG (r[k])) < 1.0e-8) {
      eval[j] = REAL (r[k]);
      ++j;
    }
  }
  *nval = nv = j;
  ZALLOC (ind, nv, int);
  sortit (eval, ind, nv);


  for (k = 0; k < nv; ++k) {
    j = ind[k];
    for (l = 0; l < n; ++l) {
      if (evec != NULL)
        evec[k * n + l] = REAL (v[l][j]);
    }
  }
  free (ind);
}

void
svdx (double *X, int m, int n, double *U, double *S, double *VT, int *nval)
     /* X is input m x n matrix.  SVD is X = U S VT.
        nval is number of positive singular values: at most min(m,n)
        S is diagonal, hence represented as array of length min(m,n) */
{

  /* returns SVD */

#define MMAX 11555
#define TDA NMAX
#define TDB 0
#define TDQ 0
#define TDPT NMAX

  double aX[MMAX][NMAX];
  double aVT[NMAX][NMAX];
  double dummyaU[1][1];
  double dummyaB[1][1];
  double dummyE[1][1];
  Integer iter[1];
  static NagError fail;
  Integer failinfo;
  int k, a, b, nvalval;

  if (n > NMAX)
    fatalx ("(svdx) n too large\n");
  if (m > MMAX)
    fatalx ("(svdx) m too large\n");

  for (k = 0; k < m; k++)
    copyarr (X + k * n, aX[k], n);

  f02wec (m, n, (double *) aX, (Integer) TDA, 0,
          (double *) dummyaB, (Integer) TDB, (Boolean) TRUE,
          (double *) dummyaU, (Integer) TDQ, (double *) S, (Boolean) TRUE,
          (double *) aVT, (Integer) TDPT, iter, (double *) dummyE, &failinfo,
          &fail);

  printf ("f02wec ran and produced fail parameter %d\n", fail);

  for (a = 0; a < m; a++) {
    for (b = 0; b < n; b++)
      U[a * n + b] = aX[a][b];
  }

  for (a = 0; a < n; a++) {
    for (b = 0; b < n; b++)
      VT[a * n + b] = aVT[a][b];
  }

  nvalval = 0;
  for (k = 0; (k < m) && (k < n); k++) {
    if (S[k] > 1.0e-8)
      nvalval++;
  }
  *nval = nvalval;

  /* may still want to sort singular values */
}

void
linsolv (double *ans, double *mat, double *rhs, int n)
{

/** 
 returns solution of mat.ans = rhs 
 unsymmetric matrix mat.  

*/

#define TDA NMAX
#define TDV NMAX
  double amat[NMAX][TDA];

  static NagError fail;
  int j, k, l, nv;
  int *ind;

  double *trhs;

#define COMPLEX(A)  A.re, A.im
#define REAL(A)  A.re
#define IMAG(A)  A.im

  if (n > NMAX)
    fatalx ("(linsolv) n too large\n");

  for (k = 0; k < n; k++) {
    copyarr (mat + k * n, amat[k], n);
  }
  ZALLOC (trhs, n, double);
  copyarr (rhs, trhs, n);

  f04arc (n, (double *) amat, (Integer) TDA, trhs, ans, &fail);

  free (trhs);

}
