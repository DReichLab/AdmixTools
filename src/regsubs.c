#include "regsubs.h"

extern int verbose;
void squishx (double *xmat, double *mat, int nrow, int oldc, int *cols,
	      int newc);

double
regressit (double *ans, double *eq, double *rhs, int m, int n)
// return 2 log likelihood score
{
  double *co, *rr, *ww, *w2;
  double vres, vbase, ynum, y, trace;
  double *traneq;
  int i, j, k, ret;

  ZALLOC (co, n * n, double);
  ZALLOC (rr, n, double);
  ZALLOC (ww, m, double);
  ZALLOC (w2, m, double);

  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      rr[j] += eq[i * n + j] * rhs[i];
      for (k = j; k < n; k++) {
	co[j * n + k] = co[k * n + j] += eq[i * n + j] * eq[i * n + k];
      }
    }
  }

/**
 y = 1.00001 ;
 for (j=0; j<n ; j++) {  
  co[j*n+j] *= y ;
 }
*/

  if (verbose) {
    printf ("coeffs:\n");
    printmat (co, n, n);
    printf ("\n\n");
    printmat (rr, n, 1);

    for (i = 0; i < n; i++) {
      printf ("diag: %3d %9.3f\n", i, co[i * n + i]);
    }

    fflush (stdout);
  }


  ret = solvit (co, rr, n, ans);
  if (ret < 0) {
    printf("*** warning bad regress\n") ;
    printmatwl(eq, m, n, n) ;
    printnl() ; 
    printmatwl(co, n, n, n) ;
    printnl() ;
    printmatwl(rr, n, n, n) ;
    return -1.0e20 ;
  }
  for (i = 0; i < m; i++) {
    ww[i] = rhs[i] - vdot (ans, eq + i * n, n);
  }

  ynum = (double) m;
  vres = asum2 (ww, m) / ynum;
  vbase = asum2 (rhs, m) / ynum;

/**
 printf("zzreg  %15.9f  %15.9f\n", log(vbase), log(vres)) ;
 printmat(rr, 1, n) ;
 printmat(co, n, n) ;
 printf("\n") ;
*/

  free (co);
  free (rr);
  free (ww);
  free (w2);
  return ynum * log (vbase / vres);
}

void
regressitall (char **vname, double *eq, double *rhs, int m, int n)
{
  double *ans;
  int i, j, k, wt;
  int npow;
  int **tab, *tweight, *cols;
  double *teq;
  double yscore;

  npow = (int) pow (2.0, (double) n);
  ZALLOC (tab, npow, int *);
  ZALLOC (tweight, npow, int);
  ZALLOC (cols, n, int);

  for (k = 0; k < npow; k++) {
    ZALLOC (tab[k], n, int);
  }
  for (k = 1; k < npow; k++) {
    add1 (tab[k], tab[k - 1], n);
    tweight[k] = intsum (tab[k], n);
  }
  ZALLOC (ans, n, double);
  ZALLOC (teq, m * n, double);

  for (wt = 1; wt <= n; ++wt) {
    printf ("weight: %d\n", wt);
    for (k = 0; k < npow; ++k) {
      if (tweight[k] != wt)
	continue;
      for (i = 0, j = 0; i < n; i++) {
	if (tab[k][i] == 0)
	  continue;
	cols[j] = i;
	++j;
      }
      squishx (teq, eq, m, n, cols, wt);
      yscore = regressit (ans, teq, rhs, m, wt);
      printf ("chisq: %9.3f\n", yscore);
      for (i = 0, j = 0; i < n; i++) {
	if (tab[k][i] == 0)
	  continue;
	printf ("%15s %9.3f\n", vname[i], ans[j]);
	++j;
      }
      printf ("\n");
    }
  }


  free (ans);
  free (teq);

  for (k = 0; k < npow; k++) {
    free (tab[k]);
  }
  free (tab);
  free (tweight);
  free (cols);
}

void
add1 (int *a, int *b, int n)
// b is 0, 1 vector as base 2 integer.  a = b + 1
{
  if (n == 0)
    return;
  copyiarr (b, a, n);
  a[n - 1] = b[n - 1] + 1;
  if (a[n - 1] == 2) {
    a[n - 1] = 0;
    add1 (a, b, n - 1);
  }
}

// now logistic regression stuff

double
logregressit (double *ans, double *eq, double **rhs, int neq, int nv)
// return log likelihood NOT chi-sq
{
  double *p, *z, *q;
  double *n0, *n1, *tans;
  double *grad, *hess, rr[2];
  double y0, y1, y, ylike, ybase, yold;
  int i, j;
  int iter, numiter = 10;
  int ret;

  ZALLOC (p, neq, double);
  ZALLOC (q, neq, double);
  ZALLOC (z, neq, double);
  ZALLOC (n0, neq, double);
  ZALLOC (n1, neq, double);
  ZALLOC (tans, neq, double);
  ZALLOC (grad, nv, double);
  ZALLOC (hess, nv * nv, double);


  for (i = 0; i < neq; i++) {
    y0 = n0[i] = rhs[i][0];
    y1 = n1[i] = rhs[i][1];
    y0 += 1.0;
    y1 += 1.0;
    y = y1 / (y0 + y1);
    y = MIN (y, 0.75);
    y = MAX (y, 0.25);
// may need changing for some problems
    p[i] = y;
  }
  y0 = asum (n0, neq);
  y1 = asum (n1, neq);
  y = y1 / (y0 + y1);
  for (i = 0; i < neq; i++) {
    p[i] = (p[i] + y) / 2.0;
    if (p[i] < 0.0)
      fatalx ("bugbug\n");
    if (p[i] > 1.0)
      fatalx ("bugbug\n");
  }

  if (verbose) {
    vzero (rr, 2);
    for (j = 0; j < neq; j++) {
      vvp (grad, grad, eq + j * nv, nv);
      vvp (rr, rr, rhs[j], 2);
      addouter (hess, eq + j * nv, nv);
    }
    y = 1.0 / (double) neq;
    vst (grad, grad, y, nv);
    vst (rr, rr, y, 2);
    vst (hess, hess, y, nv * nv);
    printf ("## averages\n");
    printmat (grad, 1, nv);
    printmat (rr, 1, 2);
    printmat (hess, nv, nv);
  }

  ptoz (p, z, neq);
  regressit (ans, eq, z, neq, nv);
  for (j = 0; j < neq; j++) {
    z[j] = vdot (eq + j * nv, ans, nv);
  }
  ybase = zlike (eq, n0, n1, ans, neq, nv);
  ztop (p, z, neq);

  calcgh (grad, hess, eq, z, n0, n1, neq, nv);
  y = .001;
  for (i = 0; i < nv; i++) {
    if (!verbose)
      break;
    copyarr (ans, tans, nv);
    tans[i] += y;
    ylike = zlike (eq, n0, n1, tans, neq, nv);
    printf ("zzgrad %3d %12.6f %12.6f\n", i, ylike - ybase, grad[i] * y);
  }

  for (iter = 1; iter <= numiter; ++iter) {
    calcgh (grad, hess, eq, z, n0, n1, neq, nv);
    ret = solvit (hess, grad, nv, tans);
    if (ret < 0)
      return -1000.0;
    if (verbose) {
      printf ("zzzz\n");
      printmat (ans, 1, nv);
      printmat (grad, 1, nv);
      printmat (hess, nv, nv);
      printmat (tans, 1, nv);
      printf ("\n\n");
    }
    vvp (ans, ans, tans, nv);
    for (j = 0; j < neq; j++) {
      z[j] = vdot (eq + j * nv, ans, nv);
    }
    ylike = zlike (eq, n0, n1, ans, neq, nv);
/**
  if (verbose)  {
   printf("iter: %3d  llike: %15.9f incr: %15.9f\n", iter, ylike, ylike-ybase) ;
   printmat(ans, 1, nv) ;
  }
*/
    if ((iter > 1) && (ylike < (yold + .0001)))
      break;
    yold = ylike;
  }

  free (p);
  free (q);
  free (z);
  free (n0);
  free (n1);
  free (tans);
  free (grad);
  free (hess);

  return ylike;
}

void
ptoz (double *p, double *z, int n)
{
  double *w1, *w2;

  ZALLOC (w1, n, double);
  ZALLOC (w2, n, double);

  vst (w2, p, -1.0, n);
  vsp (w2, w2, 1.0, n);		// q 
  vvd (w1, p, w2, n);
  vlog (z, w1, n);
  free (w1);
  free (w2);
}

void
ztop (double *p, double *z, int n)
{
  double *ww, *w1;

  ZALLOC (ww, n, double);
  ZALLOC (w1, n, double);

  vexp (ww, z, n);
  vsp (w1, ww, 1.0, n);		// 1 + e^z 
  vvd (p, ww, w1, n);		// p 

  free (ww);
  free (w1);
}

void
calcgh (double *grad, double *hess, double *eq, double *z,
	double *n0, double *n1, int neq, int nv)
{


  double *ww, *w1, *w2, *x0, *x1;
  int j;

  ZALLOC (ww, neq, double);
  ZALLOC (w1, neq, double);
  ZALLOC (w2, neq, double);
  ZALLOC (x0, neq, double);
  ZALLOC (x1, neq, double);

  vexp (ww, z, neq);
  vsp (w1, ww, 1.0, neq);
  vvt (w2, w1, w1, neq);
  vvt (x0, n0, ww, neq);
  vvm (x0, x0, n1, neq);
  vvd (x0, x0, w1, neq);

  vvp (x1, n0, n1, neq);
  vvt (x1, x1, ww, neq);
  vvd (x1, x1, w2, neq);

  vzero (grad, nv);
  vzero (hess, nv * nv);

  for (j = 0; j < neq; j++) {
    vst (ww, eq + j * nv, x0[j], nv);
    vvm (grad, grad, ww, nv);
    vst (ww, eq + j * nv, sqrt (x1[j]), nv);
    addouter (hess, ww, nv);	// actually -hess 
  }
  free (ww);
  free (w1);
  free (w2);
  free (x0);
  free (x1);
}

double
zlike (double *eq, double *n0, double *n1, double *ans, int neq, int nv)
{
  double *z, *p, *q;
  double ylike, pprob, qprob, y0, y1, ybase;
  int j;

  ZALLOC (z, neq, double);
  ZALLOC (p, neq, double);
  ZALLOC (q, neq, double);

  y0 = asum (n0, neq);
  y1 = asum (n1, neq);

  y0 += 1.0e-10;
  y1 += 1.0e-10;

  pprob = y1 / (y0 + y1);
  qprob = y0 / (y0 + y1);

  ybase = y1 * log (pprob) + y0 * log (qprob);

  for (j = 0; j < neq; j++) {
    z[j] = vdot (eq + j * nv, ans, nv);
  }
  ztop (p, z, neq);
  vst (q, p, -1.0, neq);
  vsp (q, q, 1.0, neq);
  ylike = vldot (n1, p, neq) + vldot (n0, q, neq);
  ylike -= ybase;
  free (z);
  free (p);
  free (q);
  return ylike;
}

double
logrscore (double *eq, double **rhs, int neq, int nv)
// test significance of last regressor  
{
  double *teq, *ans;
  double y1, y2, ychi;
  int i;

  ZALLOC (teq, neq * nv, double);
  ZALLOC (ans, nv, double);

  squish (teq, eq, neq, nv, nv - 1);

  y1 = logregressit (ans, teq, rhs, neq, nv - 1);
  y2 = logregressit (ans, eq, rhs, neq, nv);

  ychi = 2.0 * (y2 - y1);

  free (teq);
  free (ans);

  return ychi;

}

void
rcsquish (double *xmat, double *mat, int *cols, int oldn, int newn)
// copy matrix picking indices
{
  int i, j, xi, xj;

  for (i = 0; i < newn; i++) {
    for (j = 0; j < newn; j++) {

      xi = cols[i];
      xj = cols[j];
      xmat[i * newn + j] = mat[xi * oldn + xj];

    }
  }

}

/**
// Now in vsubs
void
squish (double *xmat, double *mat, int nrow, int oldc, int newc)
// in place legal !!
{
  int i;
  double *ww;

  ZALLOC (ww, nrow * newc, double);

  for (i = 0; i < nrow; i++) {
    copyarr (mat + i * oldc, ww + i * newc, newc);
  }

  copyarr (ww, xmat, nrow * newc);
  free (ww);

}
*/

void
squishx (double *xmat, double *mat, int nrow, int oldc, int *cols, int newc)
// copy cols of mat to xmat 
{
  int i, j, k;
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < newc; ++j) {
      k = cols[j];
      xmat[i * newc + j] = mat[i * oldc + k];
    }
  }
}

void
calcres (double *res, double *ans, double *eq, double *rhs, int neq, int nv)
/**   
 calculate residual
*/
{

  int i;
  for (i = 0; i < neq; i++) {
    res[i] = rhs[i] - vdot (eq + i * nv, ans, nv);
  }
}

static double
qwmax1 (double *co, double *ans, int n, int col)
{
  double *rr, *ww, *aa, *tco;
  int *cc;
  int a, t;
  double y, ymin, amax, amin;

  if (n <= 0)
    fatalx ("(qwmax1) bad call\n");
  if (col < 0) {
    if (n == 1) {
      ans[0] = 1.0;
      return co[0];
    }
    ZALLOC (rr, n, double);
    ZALLOC (ww, n, double);
    ZALLOC (aa, n, double);
    vclear (rr, 1.0, n);
    solvit (co, rr, n, aa);
    vmaxmin (aa, n, &amax, &amin);
    if (amax * amin >= 0.0) {
      y = asum (aa, n);
      vst (ans, aa, 1.0 / y, n);
      mulmat (ww, co, ans, n, n, 1);
      y = vdot (ww, ans, n);
      free (rr);
      free (ww);
      free (aa);
      return y;
    }
  }
  if (col >= 0) {
    ZALLOC (tco, n * n, double);
    ZALLOC (cc, n, int);
    ZALLOC (aa, n, double);
    a = 0;
    for (t = 0; t < n; ++t) {
      if (t == col)
	continue;
      cc[a] = t;
      ++a;
    }
    rcsquish (tco, co, cc, n, n - 1);
    y = qwmax1 (tco, aa, n - 1, -1);

    free (cc);
    free (tco);

    a = 0;
    ans[col] = 0.0;
    for (t = 0; t < n; ++t) {
      if (t == col)
	continue;
      ans[t] = aa[a];
      ++a;
    }
    free (aa);
    return y;
  }

  ZALLOC (aa, n, double);
  ymin = 1.0e40;
  for (t = col + 1; t < n; ++t) {
    y = qwmax1 (co, aa, n, t);
    if (y < ymin) {
      ymin = y;
      copyarr (aa, ans, n);
    }
  }
  free (aa);
  return ymin;
}

double
qwmax (double *co, double *ans, int n)
// return min of pos def q form co, coeffs,  constrained to sum to 1 and be non-negative
{
  return qwmax1 (co, ans, n, -1);
}

static void setx(double *xout, double *xin, int n, int *svec, double step) 
{
 double *ww ; 
 ZALLOC(ww, n, double) ;
 floatit(ww, svec, n) ; 
 vst(ww, ww, step, n) ; 
 vvp(xout, xin, ww, n) ; 
 free(ww) ; 
}
static void setco(double *f0, double *f1, double *f2, double *xb, double *xa, int n) 

{
 double *ww ;
 ZALLOC(ww, n, double) ;
 vvm(ww, xb, xa, n) ;
 *f0 = 1 ; 
 copyarr(ww, f1, n) ; 
 vzero(f2, n*n) ; 
 addouter(f2, ww, n) ; 
 free(ww) ;
}
static int setcc(double *cc, double f0, double *f1, double *f2, int n) 
{
  int i, j, a ; 
  double y ; 

  a = 0 ; 
  cc[a] = f0; ++a ; 
  copyarr(f1, cc+a, n) ;  a += n ; 
  for (i=0; i<n; ++i) {  
   for (j=i; j<n; ++j) {  
    y = f2[i*n+j] +f2[j*n+i] ;
    cc[a] = 0.5*y ; 
    ++a ; 
   }
  }
  return a ; 
}
static int unsetcc(double *cc, double *f0, double *f1, double *f2, int n) 
{
  int i, j, a ; 
  double y ; 

  a = 0 ; 
  *f0 = cc[a] ;  ++a ; 
  copyarr(cc+a, f1,  n) ;  a += n ; 
  for (i=0; i<n; ++i) {  
   for (j=i; j<n; ++j) {  
    y = cc[a] ; 
    f2[i*n+j] = f2[j*n+i]  = y ; 
    ++a ; 
   }
  }
  return a ; 
}
void hgrad(double ff(double *xx, int nn), double *x, int n, double step, double *fval, double *grad, double *hess)  
{
// evaluates ff at 3^n points and solves least squares
   int **tab ; 
   int neval, dim, ddim ; 
   double *co, *rhs, *cc, *xx, *ww ;  
   double y ; 
   double f0, *f1, *f2 ; 
   int k ; 


   neval = pow(3, n) ;  
   tab = initarray_2Dint(neval, n, 0) ;
   ZALLOC(f1, n, double) ; 
   ZALLOC(xx, n, double) ; 
   ZALLOC(f2, n*n, double) ; 
   dim = n*(n-1) ; dim = dim/2 ; dim += (2*n+1) ;
   ZALLOC(co, dim*neval, double) ; 
   ZALLOC(rhs, neval, double) ; 
   ZALLOC(cc, dim, double) ; 
   for (k=0; k<neval; ++k) { 
    dekodeitb(tab[k], k, n, 3) ; 
    ivsp(tab[k], tab[k], -1, n) ;
   }
   for (k=0; k<neval; ++k) { 
     setx(xx, x, n, tab[k], step) ; 
     y = ff(xx,n) ;  
     setco(&f0, f1, f2, xx, x, n)  ; 
     ddim = setcc(cc, f0, f1, f2, n) ;
     if (ddim != dim) fatalx("(hgrad) bad bug\n") ;
     rhs[k] = y ; 
     copyarr(cc, co+k*dim, dim) ;
   }

   regressit(cc, co, rhs, neval, dim) ;   


   ddim = unsetcc(cc, fval, grad, hess, n) ;

   free2Dint(&tab, neval) ;
   free(f1) ;
   free(f2) ;
   free(co) ; 
   free(rhs) ; 
   free(cc) ; 
   free(xx) ; 

}

void xline(double *line, double *coeff, int m, int n) 
{
// coeff is ideally rank n-1;  line is linear relation among columns
 double *bigco, *rhs, y, z, *aa ;
 double ww[100] ; 
 
 if (m<(n-1)) fatalx("(xline) m is too small: %d %d\n", m, n) ;

 ZALLOC(bigco, (m+1)*n, double) ; 
 ZALLOC(rhs, m+1, double) ; 
 vmaxmin(coeff, m*n, &y, NULL) ;
 copyarr(coeff, bigco, m*n) ; 
 aa = bigco + n*n ; 
 vclear(bigco+m*n, 1, n) ; 
 rhs[m] = 1 ;  
 z = regressit(line, bigco, rhs, m+1, n) ;
 if (z < -1000*1000) vclear(line, 1, n) ;
 y = asum(line, n) ; 
 if (y==0.0) y = 1.0e20 ;
 vst(line, line, 1.0/y, n) ;

 printmat(bigco, m+1, n) ;
 printmat(line, 1, n) ;
 mulmat(ww, bigco, line, m+1, n, 1) ; 
 printf("wwdebug: ") ; printmat(ww, 1, m+1) ;  // ideally 0 0 0 1
 printf("wwsum: %12.5f\n", asum(line, n)) ;

 free(bigco) ; 
 free(rhs) ;

}
