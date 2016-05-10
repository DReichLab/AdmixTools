#include  "ranmath.h"

double
gauss ()
{

/** 
 Numer alg. in C pp 289 ff
*/

  static int iset = 0;
  static double gset;
  double v1, v2, rsq, fac;

  if (iset == 1) {
    iset = 0;
    return gset;
  }

  do {
    v1 = 2.0 * DRAND2 () - 1.0;
    v2 = 2.0 * DRAND2 () - 1.0;
    rsq = v1 * v1 + v2 * v2;
  } while (rsq >= 1.0 || rsq == 0.0);

  fac = sqrt (-2.0 * log (rsq) / rsq);
  gset = v1 * fac;
  iset = 1;
  return v2 * fac;

}

void
gaussa (double *a, int n)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = gauss ();

}
