#include <stdio.h>
#include <limits.h>
#include <math.h>

#include <vsubs.h>
#include <eigsubs.h>


/* ********************************************************************* */
void eigx_ (double *pmat, double *ev, int *n);
void eigxv_ (double *pmat, double *eval, double *evec, int *n);
void cdc_ (double *pmat, int *n);
void inverse_ (double *pmat, int *n);
void solve_ (double *pmat, double *v, int *n);
void geneigsolve_ (double *pmat, double *qmat, double *eval, int *n);

void packsym (double *pmat, double *mat, int n);


void
eigvals (double *mat, double *evals, int n)
{
  double *pmat;
  int len;

  len = n * (n + 1);
  len /= 2;
  ZALLOC (pmat, len, double);

  vst (mat, mat, -1.0, n * n);
  packsym (pmat, mat, n);
  eigx_ (pmat, evals, &n);
  free (pmat);
  vst (mat, mat, -1.0, n * n);
  vst (evals, evals, -1.0, n);
}

void
eigvecs (double *mat, double *evals, double *evecs, int n)
{
  double *pmat;
  int len;

  len = n * (n + 1);
  len /= 2;
  ZALLOC (pmat, len, double);

  vst (mat, mat, -1.0, n * n);
  packsym (pmat, mat, n);

  eigxv_ (pmat, evals, evecs, &n);
  free (pmat);
  vst (mat, mat, -1.0, n * n);
  vst (evals, evals, -1.0, n);
}

/*  note:  dpotrf requires the entire matrix, not packed lower-tri */
void
chdecomp (double *mat, int n)
{
  /* symetric matrix - don't need to 
   * convert to column major order */

  cdc_ (mat, &n);
}

void
inverse (double *mat, int n)
{
  int i, j;

  /* convert to column-major order */
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      double t = mat[n * i + j];
      mat[n * i + j] = mat[n * j + i];
      mat[n * j + i] = t;
    }
  }

  /*** DEBUGGING: ***/
  {
    FILE *fid = fopen ("eigsubs.dbg", "a");
    fprintf (fid, "matrix U\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        fprintf (fid, "%8.4f ", mat[i * n + j]);
      }
      fprintf (fid, "\n");
    }
  }

  /*******************/

  inverse_ (mat, &n);

  /*** DEBUGGING: ***/
  {
    FILE *fid = fopen ("eigsubs.dbg", "a");
    fprintf (fid, "inverse of matrix U\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        fprintf (fid, "%8.4f ", mat[i * n + j]);
      }
      fprintf (fid, "\n");
    }
  }

  /*******************/

  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      double t = mat[n * i + j];
      mat[n * i + j] = mat[n * j + i];
      mat[n * j + i] = t;
    }
  }
}

void
solve (double *mat, double *b, double *v, int n)
{
  int i, j;

  double *mat2 = (double *) malloc (n * n * sizeof (double));

  /* lapack is going to put the lu-decomp into the matrix,
   * so make a copy and convert to column-major order */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      mat2[n * i + j] = mat[n * j + i];
    }
  }

  /* copy b into v */
  for (i = 0; i < n; i++) {
    v[i] = b[i];
  }

  solve_ (mat2, v, &n);

  free (mat2);
  return;
}

void
packsym (double *pmat, double *mat, int n)
        //  lapack L mode (fortran)
{
  int i, j, k = 0;
  for (i = 0; i < n; i++) {
    for (j = i; j < n; j++) {
      pmat[k] = mat[i * n + j];
      ++k;
    }
  }
}

void
geneigsolve (double *pmat, double *qmat, double *evec, double *eval, int n)
{

  /* save copy of A and B, which LAPACK will overwrite */
  double *amat = (double *) malloc (n * n * sizeof (double));
  double *bmat = (double *) malloc (n * n * sizeof (double));

  int i, j;
  for (i = 0; i < n * n; i++) {
    amat[i] = pmat[i];
    bmat[i] = qmat[i];
  }




  {
    FILE *fid = fopen ("eigsubs.dbg", "a");
    fprintf (fid, "matrix A\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        fprintf (fid, "%8.4f ", amat[i * n + j]);
      }
      fprintf (fid, "\n");
    }

    fprintf (fid, "matrix B\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        fprintf (fid, "%8.4f ", bmat[i * n + j]);
      }
      fprintf (fid, "\n");
    }
    fclose (fid);
  }


  /* matrices have to be symetric-definite, so don't 
   * need to convert to column-major order */
  geneigsolve_ (pmat, qmat, eval, &n);


  /* copy eigenvectors to A and original A,B back */
  /* ith eigenvector should be in row i */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      evec[i * n + j] = pmat[i * n + j];        /* don't put back in row-major order (?) */
    }
  }
  for (i = 0; i < n * n; i++) {
    pmat[i] = amat[i];
    qmat[i] = bmat[i];
  }

  /* LAPACK puts evals and evecs in ascending order           */
  /* reorder evals and evecs so evals are in descending order */

  for (i = 0; i < n / 2; i++) {
    double t = eval[i];
    eval[i] = eval[n - 1 - i];
    eval[n - 1 - i] = t;

    for (j = 0; j < n; j++) {   /* exchange row i and row(n-1-i) */
      t = evec[i * n + j];
      evec[i * n + j] = evec[(n - 1 - i) * n + j];
      evec[(n - 1 - i) * n + j] = t;
    }
  }

  free (amat);
  free (bmat);

}

void
mkorth (double *orth, double *ww, int n)
// special purpose.  Construct basis of vectors orthogonal to ww
{

  double *vv, *evec, *qq;
  double y;

  ZALLOC (vv, n * n, double);
  ZALLOC (evec, n * n, double);
  ZALLOC (qq, n, double);

  y = asum2 (ww, n);
  vst (qq, ww, 1.0 / sqrt (y), n);
  setidmat (vv, n);
  addouter (vv, qq, n);

  eigvecs (vv, orth, evec, n);
  copyarr (evec + n, orth, n * (n - 1));

  free (vv);
  free (qq);
  free (evec);

}
