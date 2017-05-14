#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <nicklib.h>
#include "f4rank.h"
#include "eigsubs.h"

extern int verbose;
extern double yscale ;

void addscaldiag(double *mat, double scal, int n) ;

int
dofrank (int m, int n, int rank)
{
  int x, t, i;

  if (rank >= m)
    return 0;
  if (rank >= n)
    return 0;
  if (rank == 0)
    return m * n;
  t = (m + n) * rank;
  t -= rank * rank;
  x = MAX (m * n - t, 0);
  return x;

}

int
dofrankfix (int m, int n, int rank, int nfix)
{
  int x, t, i;


  if (nfix > rank)
    return 0;
  if (rank >= m)
    return 0;
  if (rank >= n)
    return 0;
  if (rank == 0)
    return m * n;
  t = 0;

  t += m * (rank - nfix);
  t += n * rank;
  t -= rank * (rank - nfix);

  x = MAX (m * n - t, 0);	// m*n is saturated model
  return x;

}


void
doranktestfix (double *mean, double *var, int m, int n, int xrank,
	       F4INFO * f4pt, int *vfix)
{

  double y, tail = 1.0;
  int dof, t, nfix;
  double *ww;
  int rank = xrank;


  nfix = intsum (vfix, m);
  if (nfix == m) {
    y = ranktest (mean, var, m, n, 0, NULL, NULL);
    dof = m * n;
    rank = 0;
  }

  else {
    y = ranktestfix (mean, var, m, n, rank, f4pt->A, f4pt->B, vfix);
    dof = dofrankfix (m, n, rank, nfix);
  }
  

  f4pt -> chisq = 0 ; 

  if (dof > 0) { 
    tail = rtlchsq (dof, y);
    f4pt->chisq = y;
  }

   f4pt->dof = (double) dof;

  copyarr (mean, f4pt->mean, m * n);
  if (rank > 0) {
    mulmat (f4pt->resid, f4pt->A, f4pt->B, m, rank, n);
  }
  else {
    vzero (f4pt->resid, m * n);
    f4pt->A = f4pt->B = NULL;
  }
  vvm (f4pt->resid, f4pt->mean, f4pt->resid, m * n);
}

void
doranktest (double *mean, double *var, int m, int n, int rank, F4INFO * f4pt)
{

  double y, tail = 1.0;
  int dof;
  double *ww;

  y = ranktest (mean, var, m, n, rank, f4pt->A, f4pt->B);
  dof = dofrank (m, n, rank);

  f4pt -> chisq = 0 ; 
  if (dof > 0) {
    tail = rtlchsq (dof, y);
    f4pt->chisq = y;
  }

  f4pt->dof = (double) dof;

  copyarr (mean, f4pt->mean, m * n);
  if (rank > 0) {
    mulmat (f4pt->resid, f4pt->A, f4pt->B, m, rank, n);
  }
  else {
    vzero (f4pt->resid, m * n);
  }
  vvm (f4pt->resid, f4pt->mean, f4pt->resid, m * n);

//  printf4info(f4pt) ;
//  printf("dorank: %3d dof: %3d chisq: %9.3f  tail: %12.6f\n", rank, dof, y, tail) ;


}

void
printstrmat (char **ss, double *mat, int m, int n)
{
  int a;
  double *x, *scale, y;

  if (m == 0)
    return;
  if (n == 0)
    return;
  ZALLOC (x, m * n, double);
  ZALLOC (scale, n, double);

  transpose (x, mat, m, n);
  for (a = 0; a < n; ++a) {
    y = asum2 (x + a * m, m);
    if (y >= 0.0) {
      y = sqrt (y);
      y /= sqrt ((double) m);
      scale[a] = 1.0 / y;
    }
  }
  printf ("%15s ", "scale");
  printmat (scale, 1, n);


  for (a = 0; a < m; ++a) {
    printf ("%15s ", ss[a]);
    vvt (x, mat + a * n, scale, n);
    printmat (x, 1, n);
  }

  free (x);
  free (scale);
}

void
printf4info (F4INFO * f4pt)
{
  int m, n, rank;
  double *BT;

  m = f4pt->nl;
  n = f4pt->nr;
  rank = f4pt->rank;
  printf ("f4info: \n");
  printf ("f4rank: %d", rank);
  printf (" dof: %6.0f", f4pt->dof);
  printf (" chisq: %9.3f", f4pt->chisq);
  printf (" tail: %20.9g", rtlchsq (f4pt->dof, f4pt->chisq));
  printf (" dofdiff: %6.0f", f4pt->dofdiff);
  printf (" chisqdiff: %9.3f", f4pt->chisqdiff);
  printf (" taildiff: %20.9g", rtlchsq (f4pt->dofdiff, f4pt->chisqdiff));
//printf(" basepops: %s %s", f4pt -> lbase, f4pt -> rbase) ; 
  printnl ();

  if (rank > 0) {
    ZALLOC (BT, n * rank, double);
    transpose (BT, f4pt->B, f4pt->rank, n);
    printf ("B:\n");
    printstrmat (f4pt->rpops, BT, n, rank);
    free (BT);
    printf ("A:\n");
    printstrmat (f4pt->lpops, f4pt->A, m, rank);
  }
  printnl ();
}

int
solvitforcez (double *coeffs, double *rhs, int dim, double *ans, int *vl,
	      int nf)
{
  int ret ; 
  double *vals;
  double *tco ;

  ZALLOC (vals, dim, double);
  ZALLOC (tco, dim*dim, double);

  copyarr(coeffs, tco, dim*dim) ; 
  addscaldiag(tco, yscale, dim) ;
  
  ret = solvitfix(tco, rhs, dim, ans, vl, vals, nf) ; 

  if (ret<0) fatalx("bad solvitforcez\n") ; 
  free(vals) ; 
  free(tco) ;

  return ret ; 
  

}

void addscaldiag(double *mat, double scal, int n) 
{
 double y ;  
 int i, k ; 

  y = scal * trace(mat, n) ;
  for (i=0; i<n; ++i) { 
   k = i*n + i ; 
   mat[k] += y ; 
  }
}
 
double
ranktestfix (double *mean, double *var, int m, int n, int rank, double *pA,
	     double *pB, int *vfix)
// vfix is m long
{
  int d = m * n, dd;
  int f, i, j, k, l, a, b, s, t, r1, r2, k1, k2, u1, u2;
  int iter, numiter = 20, ret;
  double *ww, *varinv;
  double T2, tail;
  double *A, *B, *wleft, *wright, *mt, *evecs, *xmean;
  double y, y0, y1;
  int adim, bdim, tdim;
  double *coeffs, *rhs, *ans;

  int *vl;
  double *vfl;

  int nf;

  t = intsum (vfix, m);
  if (t > rank)
    fatalx ("(ranktestfix) too many fixed variables\n");
  ZALLOC (vfl, m * rank, double);
  nf = 0;
  ZALLOC (vl, m * rank, int);
  k = 0;
  for (i = 0; i < m; i++) {
    if (vfix[i] == 0)
      continue;
    for (j = 0; j < m; ++j) {
      if (i == j)
	continue;
      l = vl[nf] = j * rank + k;	// column k
// variables to force to zero
      if (verbose)
	printf ("zzvl %d %9.3f\n", nf, vl[nf], vfl[nf]);
      ++nf;
    }
    ++k;
  }



  dd = d * d;
  ZALLOC (varinv, dd, double);
  ZALLOC (ww, d, double);

  copyarr(var, varinv, d*d) ; 
  addscaldiag(varinv, yscale, d) ;
  pdinv (varinv, varinv, d);

  adim = rank * m;
  bdim = rank * n;
  tdim = MAX (adim, bdim);
  ZALLOC (A, adim, double);
  ZALLOC (B, bdim, double);
  ZALLOC (wright, n * n, double);
  ZALLOC (evecs, n * n, double);
  ZALLOC (mt, m * n, double);
  ZALLOC (xmean, m * n, double);

  ZALLOC (coeffs, tdim * tdim, double);
  ZALLOC (rhs, tdim, double);
  ZALLOC (ans, tdim, double);

/** 
 initialize 
*/

/**
  transpose(mt, mean, m, n) ;
  mulmat(wright, mt, mean, n, m, n) ;
  eigvecs(wright, ww, evecs, n) ;

  for (i=0; i<rank; ++i) {  
   copyarr(evecs+i*n, B+i*n, n) ;
  }
  for (a=0; a<m; ++a) {  
   for (i=0; i<rank; ++i) {  
    A[a*rank+i] = vdot(mean+a*n, B+i*n, n) ; 
   }
  }
*/

  gaussa (A, adim);
  gaussa (B, bdim);

  normab (A, B, m, n, rank);

/**    
  printf ("B:\n") ;  printmat(B, rank, n) ;
  printf ("A:\n") ;  printmat(A, m, rank) ; 
*/

  y0 = scx (varinv, mean, xmean, d);
  mulmat (xmean, A, B, m, rank, n);
  y1 = scx (varinv, mean, xmean, d);
  if (verbose)
    printf ("init scores:  %9.3f %9.3f\n", y0, y1);

  for (iter = 1; iter <= numiter; ++iter) {
    vzero (coeffs, tdim * tdim);
    vzero (rhs, tdim);

    for (i = 0; i < m; ++i) {
      for (s = 0; s < n; ++s) {
	a = i * n + s;
	for (j = 0; j < m; ++j) {
	  for (t = 0; t < n; ++t) {
	    b = j * n + t;
	    y = varinv[a * d + b];
	    for (r1 = 0; r1 < rank; ++r1) {
	      k1 = i * rank + r1;
	      u1 = r1 * n + s;
	      rhs[u1] += y * A[k1] * mean[b];
	      for (r2 = 0; r2 < rank; ++r2) {
		k2 = j * rank + r2;
		u2 = r2 * n + t;
		coeffs[u1 * bdim + u2] += y * A[k1] * A[k2];
	      }
	    }
	  }
	}
      }
    }

    addscaldiag(coeffs, yscale, bdim) ;
    solvit (coeffs, rhs, bdim, ans);

    copyarr (ans, B, bdim);
    normab (A, B, m, n, rank);
    mulmat (xmean, A, B, m, rank, n);
    y1 = scx (varinv, mean, xmean, d);
    if (verbose)
      printf ("iter B:  %9.3f\n", y1);

    vzero (coeffs, tdim * tdim);
    vzero (rhs, tdim);
    for (i = 0; i < m; ++i) {
      for (s = 0; s < n; ++s) {
	a = i * n + s;
	for (j = 0; j < m; ++j) {
	  for (t = 0; t < n; ++t) {
	    b = j * n + t;
	    y = varinv[a * d + b];
	    for (r1 = 0; r1 < rank; ++r1) {
	      k1 = i * rank + r1;
	      u1 = r1 * n + s;
	      rhs[k1] += y * B[u1] * mean[b];
	      for (r2 = 0; r2 < rank; ++r2) {
		k2 = j * rank + r2;
		u2 = r2 * n + t;
		coeffs[k1 * adim + k2] += y * B[u1] * B[u2];
	      }
	    }
	  }
	}
      }
    }
    if (verbose) {
      printnl ();
      printmat (coeffs, adim, adim);
      printmat (rhs, m, rank);
    }

    vzero (ans, adim);
    solvitforcez (coeffs, rhs, adim, ans, vl, nf);
// ret = solvitfix(coeffs, rhs, adim, ans, vl, vfl, nf) ; 

    copyarr (ans, A, adim);
    normab (A, B, m, n, rank);

    if (verbose) {
      printnl ();
      printmatl (A, m, rank);
    }

    mulmat (xmean, A, B, m, rank, n);
    y1 = scx (varinv, mean, xmean, d);
    if (verbose)
      printf ("iter A: %d %9.3f\n", iter, y1);
  }

  if (pA != NULL)
    copyarr (A, pA, adim);
  if (pB != NULL)
    copyarr (B, pB, bdim);

  free (A);
  free (B);
  free (wright);
  free (evecs);
  free (mt);
  free (ww);
  free (varinv);
  free (xmean);
  free (coeffs);
  free (rhs);
  free (ans);
  free (vfl);
  free (vl);


  return y1;
}

double
ranktest (double *mean, double *var, int m, int n, int rank, double *pA,
	  double *pB)
{
  int d = m * n, dd;
  int i, j, a, b, s, t, r1, r2, k1, k2, u1, u2;
  int iter, numiter = 20;
  double *ww, *varinv;
  double T2, tail;
  double *A, *B, *wleft, *wright, *mt, *evecs, *xmean;
  double y, y0, y1;
  int adim, bdim, tdim;
  double *coeffs, *rhs, *ans;



  dd = d * d;
  ZALLOC (varinv, dd, double);
  ZALLOC (ww, d, double);

  copyarr(var, varinv, d*d) ; 
  addscaldiag(varinv, yscale, d) ;
  pdinv (varinv, varinv, d);

  if (rank == 0) {
    vzero (ww, d);
    T2 = scx (varinv, mean, ww, d);
    tail = rtlchsq (d, T2);
    if (verbose)
      printf ("T2: %9.3f  tail: %15.9g\n", T2, tail);
    free (ww);
    free (varinv);
    return T2;
  }

  adim = rank * m;
  bdim = rank * n;
  tdim = MAX (adim, bdim);
  ZALLOC (A, adim, double);
  ZALLOC (B, bdim, double);
  ZALLOC (wright, n * n, double);
  ZALLOC (evecs, n * n, double);
  ZALLOC (mt, m * n, double);
  ZALLOC (xmean, m * n, double);

  ZALLOC (coeffs, tdim * tdim, double);
  ZALLOC (rhs, tdim, double);
  ZALLOC (ans, tdim, double);

/** 
 initialize 
*/
  transpose (mt, mean, m, n);
  mulmat (wright, mt, mean, n, m, n);
  eigvecs (wright, ww, evecs, n);

  for (i = 0; i < rank; ++i) {
    copyarr (evecs + i * n, B + i * n, n);
  }
  for (a = 0; a < m; ++a) {
    for (i = 0; i < rank; ++i) {
      A[a * rank + i] = vdot (mean + a * n, B + i * n, n);
    }
  }

  normab (A, B, m, n, rank);

/**    
  printf ("B:\n") ;  printmat(B, rank, n) ;
  printf ("A:\n") ;  printmat(A, m, rank) ; 
*/

/**
  debug 
  gaussa(rhs, tdim) ; vst(rhs, rhs, 1.0, tdim) ; vvp(B, B, rhs, bdim) ;
  gaussa(rhs, tdim) ; vst(rhs, rhs, 1.0, tdim) ; vvp(A, A, rhs, adim) ;
  normab(A, B, m, n, rank) ;
*/

  y0 = scx (varinv, mean, xmean, d);
  mulmat (xmean, A, B, m, rank, n);
  y1 = scx (varinv, mean, xmean, d);
  if (verbose)
    printf ("init scores:  %9.3f %9.3f\n", y0, y1);

  for (iter = 1; iter <= numiter; ++iter) {
    vzero (coeffs, tdim * tdim);
    vzero (rhs, tdim);

    for (i = 0; i < m; ++i) {
      for (s = 0; s < n; ++s) {
	a = i * n + s;
	for (j = 0; j < m; ++j) {
	  for (t = 0; t < n; ++t) {
	    b = j * n + t;
	    y = varinv[a * d + b];
	    for (r1 = 0; r1 < rank; ++r1) {
	      k1 = i * rank + r1;
	      u1 = r1 * n + s;
	      rhs[u1] += y * A[k1] * mean[b];
	      for (r2 = 0; r2 < rank; ++r2) {
		k2 = j * rank + r2;
		u2 = r2 * n + t;
		coeffs[u1 * bdim + u2] += y * A[k1] * A[k2];
	      }
	    }
	  }
	}
      }
    }
    
    addscaldiag(coeffs, yscale, bdim) ;
    solvit (coeffs, rhs, bdim, ans);

    copyarr (ans, B, bdim);
    normab (A, B, m, n, rank);
    mulmat (xmean, A, B, m, rank, n);
    y1 = scx (varinv, mean, xmean, d);
    if (verbose)
      printf ("iter B:  %9.3f\n", y1);

    vzero (coeffs, tdim * tdim);
    vzero (rhs, tdim);
    for (i = 0; i < m; ++i) {
      for (s = 0; s < n; ++s) {
	a = i * n + s;
	for (j = 0; j < m; ++j) {
	  for (t = 0; t < n; ++t) {
	    b = j * n + t;
	    y = varinv[a * d + b];
	    for (r1 = 0; r1 < rank; ++r1) {
	      k1 = i * rank + r1;
	      u1 = r1 * n + s;
	      rhs[k1] += y * B[u1] * mean[b];
	      for (r2 = 0; r2 < rank; ++r2) {
		k2 = j * rank + r2;
		u2 = r2 * n + t;
		coeffs[k1 * adim + k2] += y * B[u1] * B[u2];
	      }
	    }
	  }
	}
      }
    }
    if (verbose) {
      printnl ();
      printmat (coeffs, adim, adim);
      printmat (rhs, m, rank);
    }

    vzero (ans, adim);
    addscaldiag(coeffs, yscale, adim) ;
    solvit(coeffs, rhs, adim, ans);

    copyarr (ans, A, adim);
    normab (A, B, m, n, rank);

    if (verbose) {
      printnl ();
      printmatl (A, m, rank);
    }

    mulmat (xmean, A, B, m, rank, n);
    y1 = scx (varinv, mean, xmean, d);
    if (verbose)
      printf ("iter A: %d %9.3f\n", iter, y1);
  }

  if (pA != NULL)
    copyarr (A, pA, adim);
  if (pB != NULL)
    copyarr (B, pB, bdim);

  free (A);
  free (B);
  free (wright);
  free (evecs);
  free (mt);
  free (ww);
  free (varinv);
  free (xmean);
  free (coeffs);
  free (rhs);
  free (ans);


  return y1;
}


void
normab (double *A, double *B, int m, int n, int rank)
// each B vector forced to norm 1 and positive sum if possible
// not needed but seems like good practice.  Makes answer canonical
{
  int i, r;
  double y, *bpt;

  for (r = 0; r < rank; ++r) {
    bpt = B + r * n;
    y = asum2 (bpt, n);
    y = sqrt (y);
    y /= sqrt ((double) n);
    if (asum (bpt, n) < 0.0)
      y = -y;
    vst (bpt, bpt, 1.0 / y, n);
    for (i = 0; i < m; ++i) {
      A[i * rank + r] *= y;
    }
  }
}

void
f4info_init (F4INFO * f4pt, int nl, int nr, char **popllist, char **poprlist,
	     int rank)
{

  if (rank == 0) {
    f4pt->A = f4pt->B = NULL;
  }
  else {
    ZALLOC (f4pt->A, nl * rank, double);
    ZALLOC (f4pt->B, nr * rank, double);
  }
  ZALLOC (f4pt->mean, nl * nr, double);
  ZALLOC (f4pt->resid, nl * nr, double);
  f4pt->nl = nl;
  f4pt->nr = nr;
  f4pt->rank = rank;
  f4pt->lpops = popllist + 1;
  f4pt->rpops = poprlist + 1;
  f4pt->lbase = popllist[0];
  f4pt->rbase = poprlist[0];
  f4pt->dof = dofrank (nl, nr, rank);
  f4pt->dofjack = f4pt->chisq = 0.0;
  f4pt->dofdiff = f4pt->chisqdiff = 0.0;
  f4pt->bestparent = f4pt->bestchild = NULL;
}

