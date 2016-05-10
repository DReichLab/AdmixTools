#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "ldsubs.h"

#include <nicklib.h>

extern int verbose;;

static double
mk2from4 (double *xe, double *xf, double *xc, double *pe, double *pf);

static double mk4from9 (double *xd, double *xc, double *p);

static int zdipmode = 0;
static int zdphasedmode = 0;

void
setzdipmode (int mode)
{
  zdipmode = mode;
}

void
setzdphasedmode (int mode)
{
  zdphasedmode = mode;
  if (mode == 1)
    zdipmode = mode;
}


double
lddip (double *xc)
{

  double p1, p2;
  double q1, q2;
  double ylike, ybase, yl, ychi, y;
  double ylast;

  double xe[2], xf[2];
  double pe[2], pf[2];
  double p4[4], p4t[4];
  double xd[4];

  int e1, f1, k1;
  int e2, f2, k2;
  int e3, f3, k3;
  int n, iter, i;

  vclear (p4, 0.25, 4);
  vclear (pe, 0.5, 2);
  vclear (pf, 0.5, 2);

  y = asum (xc, 9);

  if (y < 2.5) {
    printf ("*** zzwarning... small counts\n");
    return -1.0;
  }

  ylast = -1.0e40;
  for (iter = 1; iter <= 10; ++iter) {
// p4 is probability of chromosome if unlinked
    yl = ylike = mk4from9 (xd, xc, p4);
    mk2from4 (xe, xf, xd, pe, pf);
    if (iter == 1)
      ybase = ylike;
    ylike -= ybase;
    y = asum (xe, 2);
    vst (pe, xe, 1.0 / y, 2);
    y = asum (xf, 2);
    vst (pf, xf, 1.0 / y, 2);
    for (e1 = 0; e1 <= 1; ++e1) {
      for (f1 = 0; f1 <= 1; ++f1) {
        k1 = 2 * e1 + f1;
        p4[k1] = pe[e1] * pf[f1];
      }
    }
    if (verbose)
      printf ("itera: %3d %9.3f\n", iter, ylike);
    if (ylike < (ylast + .0001))
      break;
    ylast = ylike;
  }
  copyarr (p4, p4t, 4);

  ylast = -1.0e40;
  for (iter = 1; iter <= 10; ++iter) {
    ylike = mk4from9 (xd, xc, p4);
    ylike -= yl;
    y = asum (xd, 4);
    vst (p4, xd, 1.0 / y, 4);
    if (verbose)
      printf ("iterb: %3d %9.3f\n", iter, ylike);

    if (ylike < (ylast + .001))
      break;
    ylast = ylike;
  }
  ychi = 2.0 * ylike;
  if (verbose) {
    printf ("%12.6f\n", ychi);
    for (i = 0; i <= 3; i++) {
      printf ("%3d %9.3f %9.3f\n", i, p4[i], p4t[i]);
    }
  }
  return ychi;
}

double
lddipx (double *xc, double *xd)
// xc 3 x 3   xd  2 x 2 
{

  double p1, p2;
  double q1, q2;
  double ylike, ybase, yl, ychi, y;

  double xe[2], xf[2];
  double pe[2], pf[2];
  double p4[4], p4t[4];
  double yd[4];

  int e1, f1, k1;
  int e2, f2, k2;
  int e3, f3, k3;
  int n, iter, i;

  vclear (p4, 0.25, 4);
  vclear (pe, 0.5, 2);
  vclear (pf, 0.5, 2);

  y = asum (xc, 9);

  if (y < 2.5) {
    printf ("*** zzwarning... small counts\n");
    return 0.0;
  }

  for (iter = 1; iter <= 10; ++iter) {
    ylike = mk4from9 (yd, xc, p4);
    ylike += vldot (xd, p4, 4);
    yl = ylike;
    vvp (yd, yd, xd, 4);
    mk2from4 (xe, xf, yd, pe, pf);
    if (iter == 1)
      ybase = ylike;
    ylike -= ybase;
    y = asum (xe, 2);
    vst (pe, xe, 1.0 / y, 2);
    y = asum (xf, 2);
    vst (pf, xf, 1.0 / y, 2);
    for (e1 = 0; e1 <= 1; ++e1) {
      for (f1 = 0; f1 <= 1; ++f1) {
        k1 = 2 * e1 + f1;
        p4[k1] = pe[e1] * pf[f1];
      }
    }
    if (verbose)
      printf ("iterxa: %3d %9.3f\n", iter, ylike);
  }
  copyarr (p4, p4t, 4);

  for (iter = 1; iter <= 10; ++iter) {
    ylike = mk4from9 (yd, xc, p4);
    ylike += vldot (xd, p4, 4);
    vvp (yd, yd, xd, 4);
    ylike -= yl;
    y = asum (yd, 4);
    vst (p4, yd, 1.0 / y, 4);
    if (verbose)
      printf ("iterxb: %3d %9.3f\n", iter, ylike);
  }
  ychi = 2.0 * ylike;
  if (verbose) {
    for (i = 0; i <= 3; i++) {
      printf ("%3d %9.3f %9.3f\n", i, p4[i], p4t[i]);
    }
  }
  return ychi;
}

double
mk4from9 (double *xd, double *xc, double *p)
{
  double y, ylike;

  int c;
  int e1, e2, e3;
  int f1, f2, f3, k1, k2, k3;
  int k, t;
  double pp[9];

  vclear (pp, 1.0e-20, 9);

  for (e1 = 0; e1 <= 1; ++e1) {
    for (e2 = 0; e2 <= 1; ++e2) {
      for (f1 = 0; f1 <= 1; ++f1) {
        for (f2 = 0; f2 <= 1; ++f2) {
          k1 = 2 * e1 + f1;
          k2 = 2 * e2 + f2;
          e3 = e1 + e2;
          f3 = f1 + f2;
          k3 = 3 * e3 + f3;
          pp[k3] += p[k1] * p[k2];
        }
      }
    }
  }
  ylike = vldot (xc, pp, 9);

  vzero (xd, 4);
  for (e1 = 0; e1 <= 1; ++e1) {
    for (e2 = 0; e2 <= 1; ++e2) {
      for (f1 = 0; f1 <= 1; ++f1) {
        for (f2 = 0; f2 <= 1; ++f2) {
          k1 = 2 * e1 + f1;
          k2 = 2 * e2 + f2;
          e3 = e1 + e2;
          f3 = f1 + f2;
          k = k3 = 3 * e3 + f3;
          y = p[k1] * p[k2] / pp[k3];
          xd[k1] += y * xc[k];
          xd[k2] += y * xc[k];
        }
      }
    }
  }
  return ylike;
}

double
mk2from4 (double *xe, double *xf, double *xc, double *pe, double *pf)
{
  double y;

  int e1, e2, e3;
  int f1, f2, f3, k1, k2, k3;
  int k, t;
  double pp[4];

  vzero (xe, 2);
  vzero (xf, 2);

  for (e1 = 0; e1 <= 1; ++e1) {
    for (f1 = 0; f1 <= 1; ++f1) {
      k1 = 2 * e1 + f1;
      xe[e1] += xc[k1];
      xf[f1] += xc[k1];
      pp[k1] = pe[e1] * pf[f1];
    }
  }
  y = vldot (xc, pp, 4);
  return y;
}

double
lewont (double p, double q, double x)

/* d(ij) statistic */
{
  double d1, dmax;
  d1 = x - p * q;
  if (d1 < 0.0) {
    dmax = MIN (p * q, (1.0 - p) * (1.0 - q));
  }
  if (d1 > 0.0) {
    dmax = MIN (p * (1.0 - q), q * (1.0 - p));
  }
  if (d1 == 0.0)
    return 0.0;
  if (dmax == 0.0)
    fatalx ("(lewont) bad freq %9.3f %9.3f\n", p, q);
  return d1 / dmax;
}


double
lewontindprime (double *p4)
{
  double rs[2], cs[2], x, p, q, y;

  rowsum (p4, rs, 2);
  colsum (p4, cs, 2);

  y = asum (rs, 2);
  x = p4[0] / y;
  p = rs[0] / y;
  q = cs[0] / y;

  return fabs (lewont (p, q, x));

}



void
lewontinv (double *p4, double lew, double p, double pp, int ispos)

/* returns a 4 long distribution matching lewontins dprime lew 
ispos YES x = p4[0] >= p*pp 
      NO  x = p4[0] <= p*pp 
p is first row sum
pp is first col sum 
p, pp must NOT be 0 or 1 
*/
{
  double q, qq, dmax, d1, x;

  q = 1.0 - p;
  qq = 1.0 - pp;

  if ((p == 0) || (pp == 0.0))
    fatalx ("bad lewontinv %f %f\n", p, pp);
  if ((q == 0) || (qq == 0.0))
    fatalx ("bad lewontinv %f %f\n", p, pp);

  if (ispos == YES) {
    dmax = MIN (p * qq, pp * q);
  }
  else {
    dmax = -MIN (p * pp, q * qq);
  }
  d1 = dmax * lew;
  x = d1 + p * pp;
  p4[0] = x;
  p4[1] = p - x;
  p4[2] = pp - x;
  p4[3] = qq + x - p;
}

double
dprime (int *a1, int *a2, int n)
{

/** 
 a1, a2 are categories and contain small integers 
 computes Lewontin's dprime statistic for multialle locus
 Ref: Hedrick Genetics 117 (1987) pp 331-341 
*/
  int max1, max2, len, a, b, i;
  double *x, t, dd;
  double *f1, *f2;

  ivmaxmin (a1, n, &max1, NULL);
  ivmaxmin (a2, n, &max2, NULL);
  len = (max1 + 1) * (max2 + 1);
  ZALLOC (x, len, double);
  ZALLOC (f1, max1 + 1, double);
  ZALLOC (f2, max2 + 1, double);
  for (i = 0; i < n; i++) {
    a = a1[i];
    b = a2[i];
    ++x[a * (max2 + 1) + b];
    ++f1[a];
    ++f2[b];
  }
  t = asum (f1, max1 + 1);
  vst (f1, f1, 1.0 / t, max1 + 1);
  t = asum (f2, max2 + 1);
  vst (f2, f2, 1.0 / t, max2 + 1);
  t = asum (x, len);
  vst (x, x, 1.0 / t, len);

  t = 0.0;
  for (a = 0; a <= max1; ++a) {
    if (f1[a] == 0.0)
      continue;
    for (b = 0; b <= max2; ++b) {
      if (f2[b] == 0.0)
        continue;
      dd = lewont (f1[a], f2[b], x[a * (max2 + 1) + b]);

/**
   printf("zzlewont %9.3f %9.3f %9.3f %9.3f\n", 
   f1[a],f2[b],x[a*(max2+1)+b],dd ) ;
*/

      t += f1[a] * f2[b] * fabs (dd);
    }
  }
  if (verbose == YES) {
    for (i = 0; i < n; i++) {
      a = a1[i];
      b = a2[i];
      printf ("yy1 %4d %4d %4d\n", i, a, b);
    }
    for (a = 0; a <= max1; ++a) {
      if (f1[a] == 0.0)
        continue;
      for (b = 0; b <= max2; ++b) {
        if (f2[b] == 0.0)
          continue;
        dd = lewont (f1[a], f2[b], x[a * (max2 + 1) + b]);
        printf ("yy2 %d %d %9.3f %9.3f %9.3f %9.3f\n", a, b, f1[a], f2[b],
                x[a * (max2 + 1) + b], dd);
      }
    }

  }
  free (x);
  free (f1);
  free (f2);
  return t;
}


double
zdip0 (double *xc)
{

  double p1, p2;
  double q1, q2;
  double ylike, ybase, yl, ychi, y;
  double ylast, zscore;

  double xe[2], xf[2];
  double pe[2], pf[2];
  double p4[4], p4t[4];
  double xd[4];

  int e1, f1, k1;
  int e2, f2, k2;
  int e3, f3, k3;
  int n, iter, i;

  vclear (p4, 0.25, 4);
  vclear (pe, 0.5, 2);
  vclear (pf, 0.5, 2);

  y = asum (xc, 9);

  if (y < 2.5) {
    printf ("*** zzwarning... small counts\n");
    return -1.0;
  }

  ylast = -1.0e40;
  for (iter = 1; iter <= 10; ++iter) {
// p4 is probability of chromosome if unlinked
    yl = ylike = mk4from9 (xd, xc, p4);
    mk2from4 (xe, xf, xd, pe, pf);
    if (iter == 1)
      ybase = ylike;
    ylike -= ybase;
    y = asum (xe, 2);
    vst (pe, xe, 1.0 / y, 2);
    y = asum (xf, 2);
    vst (pf, xf, 1.0 / y, 2);
    for (e1 = 0; e1 <= 1; ++e1) {
      for (f1 = 0; f1 <= 1; ++f1) {
        k1 = 2 * e1 + f1;
        p4[k1] = pe[e1] * pf[f1];
      }
    }
    if (verbose)
      printf ("itera: %3d %9.3f\n", iter, ylike);
    if (ylike < (ylast + .001))
      break;
    ylast = ylike;
  }
  copyarr (p4, p4t, 4);

  ylast = -1.0e40;
  for (iter = 1; iter <= 10; ++iter) {
    ylike = mk4from9 (xd, xc, p4);
    ylike -= yl;
    y = asum (xd, 4);
    vst (p4, xd, 1.0 / y, 4);
    if (verbose)
      printf ("iterb: %3d %9.3f\n", iter, ylike);

    if (ylike < (ylast + .001))
      break;
    ylast = ylike;
  }
  ychi = 2.0 * ylike;

  if (verbose) {
    printf ("%12.6f\n", ychi);
    for (i = 0; i <= 3; i++) {
      printf ("%3d %9.3f %9.3f\n", i, p4[i], p4t[i]);
    }
  }

  if (ychi <= 0.0)
    return 0.0;
  zscore = sqrt (ychi);
  if (p4[0] < p4t[0])
    zscore *= -1.0;
  return zscore;
}

double
zdip (double *xc)
{

  CORR xcorr;
  CORR *corrpt;
  int a, b, t;
  static int ncount = 0, xran;
  double y;
  double ww[4];

  ++ncount;
  if (zdipmode == 0)
    return zdip0 (xc);

  if (zdphasedmode) {
    ww[0] = xc[0];
    ww[1] = xc[1];
    ww[2] = xc[3];
    ww[3] = xc[4];
    y = z2x2 (ww);

/**
  xran = ranmod(1000) ;
  if (xran == 0) {     
   printmat(ww, 2, 2) ;
   printf("ncount: %8d Z: %9.3f\n", ncount, y) ; 
  }
*/
    return y;
  }


  corrpt = &xcorr;
  clearcorr (corrpt);

  for (a = 0; a < 3; ++a) {
    for (b = 0; b < 3; ++b) {
      y = xc[3 * a + b];
      addcorrn (corrpt, a, b, y);
    }
  }
  t = calccorr (corrpt, 0, YES);
  if (t < 0)
    return 0;

/**
 xran = ranmod(1000) ;
 if (xran == 0) {     
  printmat(xc, 3, 3) ;
  printf("ncount: %8d Z: %9.3f\n", ncount, corrpt -> Z) ;
  printcorr(corrpt) ;
 }
*/
  return corrpt->Z;
}

int
calccorr (CORR * corrpt, int mode, int ztrans)
// mode = 1 => do NOT take off mean
{

  double y, yn, m1, m2, v11, v12, v22, r;

  corrpt->corr = corrpt->Z = 0.0;
  if (corrpt->S0 < 0.5)
    return -1;

  yn = corrpt->S0;
  corrpt->m1 = m1 = corrpt->S1 / corrpt->S0;
  corrpt->m2 = m2 = corrpt->S2 / corrpt->S0;

  if (mode == 1) {
    m1 = m2 = 0.0;
  }

  corrpt->v11 = v11 = (corrpt->S11 - yn * m1 * m1) / yn;
  corrpt->v12 = v12 = (corrpt->S12 - yn * m1 * m2) / yn;
  corrpt->v22 = v22 = (corrpt->S22 - yn * m2 * m2) / yn;


  y = corrpt->corr = v12 / sqrt (v11 * v22 + 1.0e-20);
  corrpt->Z = sqrt (yn) * y;

  if (ztrans) {

    if (yn < 4)
      return -1;

    y = MIN (y, 0.9);
    y = MAX (y, -0.9);

    r = 0.5 * log ((1 + y) / (1 - y));
    corrpt->Z = sqrt (yn - 3) * r;
  }
  return 1;
}

void
printcorr (CORR * corrpt)
{
  printf ("S0:   %12.3f\n", corrpt->S0);
  printf ("S1:   %12.3f\n", corrpt->S1);
  printf ("S2:   %12.3f\n", corrpt->S2);
  printf ("S11:  %12.3f\n", corrpt->S11);
  printf ("S12:  %12.3f\n", corrpt->S12);
  printf ("S22:  %12.3f\n", corrpt->S22);
  printf ("m1:   %12.3f\n", corrpt->m1);
  printf ("m2:   %12.3f\n", corrpt->m2);
  printf ("v11:  %12.3f\n", corrpt->v11);
  printf ("v12:  %12.3f\n", corrpt->v12);
  printf ("corr: %12.3f\n", corrpt->corr);
  printf ("Z:    %12.3f\n", corrpt->Z);
}

void
clearcorr (CORR * corrpt)
{
  corrpt->S0 = 0;
  corrpt->S1 = 0;
  corrpt->S2 = 0;               // was buggy
  corrpt->S11 = 0;
  corrpt->S12 = 0;
  corrpt->S22 = 0;
  corrpt->m1 = 0;
  corrpt->m2 = 0;
  corrpt->v11 = 0;
  corrpt->v12 = 0;
  corrpt->v22 = 0;
  corrpt->corr = 0;
  corrpt->Z = 0;
}

void
addcorr (CORR * corrpt, double x1, double x2)
{
  corrpt->S0 += 1;
  corrpt->S1 += x1;
  corrpt->S2 += x2;
  corrpt->S11 += x1 * x1;
  corrpt->S12 += x1 * x2;
  corrpt->S22 += x2 * x2;
}

void
addcorrn (CORR * corrpt, double x1, double x2, double yn)
{


  corrpt->S0 += yn;
  corrpt->S1 += x1 * yn;
  corrpt->S2 += x2 * yn;
  corrpt->S11 += x1 * x1 * yn;
  corrpt->S12 += x1 * x2 * yn;
  corrpt->S22 += x2 * x2 * yn;
}

void
minuscorr (CORR * out, CORR * c1, CORR * c2)
// subtract corr suff stats.  Used in jackknife
{
  out->S0 = c1->S0 - c2->S0;
  if (out->S0 < -1.0e-6)
    fatalx ("(minuscorr) S0: %15.9f\n", out->S0);
  out->S1 = c1->S1 - c2->S1;
  out->S2 = c1->S2 - c2->S2;
  out->S11 = c1->S11 - c2->S11;
  out->S12 = c1->S12 - c2->S12;
  out->S22 = c1->S22 - c2->S22;
}
