// LAPACK-using version of eigensrc/eigx.f, easier to build on OS X
// Christopher Chang (chrchang@alumni.caltech.edu), BGI Cognitive Genomics Lab

#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else

#ifdef _WIN32

#define HAVE_LAPACK_CONFIG_H
#define LAPACK_COMPLEX_STRUCTURE
#include "lapack/lapacke/include/lapacke.h"
typedef int __CLPK_integer;

typedef double __CLPK_doublereal;

#else // begin !_WIN32
#if __LP64__
typedef int __CLPK_integer;
#else
typedef long int __CLPK_integer;
#endif
typedef double __CLPK_doublereal;

int dspev_ (char *jobz, char *uplo, __CLPK_integer * n,
	    __CLPK_doublereal * ap, __CLPK_doublereal * w,
	    __CLPK_doublereal * z__, __CLPK_integer * ldz,
	    __CLPK_doublereal * work, __CLPK_integer * info);

int dpotrf_ (char *uplo, __CLPK_integer * n, __CLPK_doublereal * a,
	     __CLPK_integer * lda, __CLPK_integer * info);

int dgetrf_ (__CLPK_integer * m, __CLPK_integer * n, __CLPK_doublereal * a,
	     __CLPK_integer * lda, __CLPK_integer * ipiv,
	     __CLPK_integer * info);

int dgetri_ (__CLPK_integer * n, __CLPK_doublereal * a, __CLPK_integer * lda,
	     __CLPK_integer * ipiv, __CLPK_doublereal * work,
	     __CLPK_integer * lwork, __CLPK_integer * info);

int dgetrs_ (char *trans, __CLPK_integer * n, __CLPK_integer * nrhs,
	     __CLPK_doublereal * a, __CLPK_integer * lda,
	     __CLPK_integer * ipiv, __CLPK_doublereal * b,
	     __CLPK_integer * ldb, __CLPK_integer * info);

int dsygv_ (__CLPK_integer * itype, char *jobz, char *uplo,
	    __CLPK_integer * n, __CLPK_doublereal * a, __CLPK_integer * lda,
	    __CLPK_doublereal * b, __CLPK_integer * ldb,
	    __CLPK_doublereal * w, __CLPK_doublereal * work,
	    __CLPK_integer * lwork, __CLPK_integer * info);
#endif // end !_WIN32
#endif // end !__APPLE__

void
mem_error ()
{
  fprintf (stderr, "CM\n");
  exit (1);
}

void
inverse_error (char *procname, int info)
{
  if (info < 0) {
    fprintf (stderr, "error (%s): illegal argument %d\n", procname, -info);
  }
  else {
    fprintf (stderr, "error (%s): singular matrix %d\n", procname, info);
  }
  exit (1);
}

void
eigx_ (double *pmat, double *ev, __CLPK_integer * n)
{
  char jobz = 'N';
  char uplo = 'L';
  __CLPK_integer ldz = *n;
  __CLPK_integer info;
  double *z;
  double *work;
  z = (double *) malloc (ldz * ldz * sizeof (double));
  if (!z) {
    mem_error ();
  }
  work = (double *) malloc (3 * ldz * sizeof (double));
  if (!work) {
    free (z);
    mem_error ();
  }
  dspev_ (&jobz, &uplo, n, pmat, ev, z, &ldz, work, &info);
  free (z);
  free (work);
  if (info) {
#if __LP64__ || _WIN32
    fprintf (stderr, "INFO: %d\n", info);
#else
    fprintf (stderr, "INFO: %ld\n", info);
#endif
    exit (1);
  }
}

void
eigxv_ (double *pmat, double *eval, double *evec, __CLPK_integer * n)
{
  char jobz = 'V';
  char uplo = 'L';
  __CLPK_integer ldz = *n;
  __CLPK_integer info;
  double *work = (double *) malloc (3 * ldz * sizeof (double));
  if (!work) {
    mem_error ();
  }
  dspev_ (&jobz, &uplo, n, pmat, eval, evec, &ldz, work, &info);
  free (work);
  if (info) {
#if __LP64__ || _WIN32
    fprintf (stderr, "INFO: %d\n", info);
#else
    fprintf (stderr, "INFO: %ld\n", info);
#endif
    exit (1);
  }
}

void
cdc_ (double *pmat, __CLPK_integer * n)
{
  char uplo = 'L';
  __CLPK_integer lda = *n;
  __CLPK_integer info;
  dpotrf_ (&uplo, n, pmat, &lda, &info);
  if (info) {
    if (info < 0) {
#if __LP64__ || _WIN32
      fprintf (stderr, "error (CDC): illegal argument %d\n", -info);
#else
      fprintf (stderr, "error (CDC): illegal argument %ld\n", -info);
#endif
    }
    else {
#if __LP64__ || _WIN32
      fprintf (stderr, "error (CDC): minor not positive definite %d\n", info);
#else
      fprintf (stderr, "error (CDC): minor not positive definite %ld\n",
	       info);
#endif
    }
    exit (1);
  }
}

void
inverse_ (double *pmat, __CLPK_integer * n)
{
  __CLPK_integer lwork = (*n) * (*n);
  __CLPK_integer info;
  __CLPK_integer *ipiv;
  double *work;
  ipiv = (__CLPK_integer *) malloc ((*n) * sizeof (__CLPK_integer));
  if (!ipiv) {
    mem_error ();
  }
  work = (double *) malloc (lwork * sizeof (double));
  if (!work) {
    free (ipiv);
    mem_error ();
  }
  dgetrf_ (n, n, pmat, n, ipiv, &info);
  if (info) {
    free (ipiv);
    free (work);
    inverse_error ("INVERSE", info);
    exit (1);
  }
  dgetri_ (n, pmat, n, ipiv, work, &lwork, &info);
  free (ipiv);
  free (work);
  if (info) {
    inverse_error ("INVERSE", info);
  }
}

void
solve_ (double *pmat, double *v, __CLPK_integer * n)
{
  __CLPK_integer ldb = *n;
  char trans = 'N';
  __CLPK_integer nrhs = 1;
  double *work;
  __CLPK_integer *ipiv;
  __CLPK_integer info;
  ipiv = (__CLPK_integer *) malloc (ldb * sizeof (__CLPK_integer));
  if (!ipiv) {
    mem_error ();
  }
  work = (double *) malloc (ldb * ldb * sizeof (double));
  if (!work) {
    free (ipiv);
    mem_error ();
  }
  dgetrf_ (n, n, pmat, n, ipiv, &info);
  if (info) {
    free (ipiv);
    free (work);
    inverse_error ("SOLVE", info);
  }
  dgetrs_ (&trans, n, &nrhs, pmat, n, ipiv, v, &ldb, &info);
  free (ipiv);
  free (work);
  if (info < 0) {
    inverse_error ("SOLVE", info);
  }
}

void
geneigsolve_ (double *pmat, double *qmat, double *eval, __CLPK_integer * n)
{
  __CLPK_integer lwork = (*n) * (*n);
  double *work = (double *) malloc (lwork * sizeof (double));
  __CLPK_integer wood_elf = 1;	// Sameer Merchant memorial temporary variable
  __CLPK_integer info;
  if (!work) {
    mem_error ();
  }
  dsygv_ (&wood_elf, "V", "U", n, pmat, n, qmat, n, eval, work, &lwork,
	  &info);
  free (work);
  if (info && (info <= 2 * (*n))) {
    if (info < 0) {
#if __LP64__ || _WIN32
      fprintf (stderr, "error (GENEIGSOLVE): illegal argument %d\n", -info);
#else
      fprintf (stderr, "error (GENEIGSOLVE): illegal argument %ld\n", -info);
#endif
    }
    else if (info <= (*n)) {
#if __LP64__ || _WIN32
      fprintf (stderr, "error (GENEIGSOLVE): failure to converge %d\n", info);
#else
      fprintf (stderr, "error (GENEIGSOLVE): failure to converge %ld\n",
	       info);
#endif
    }
    else {
#if __LP64__ || _WIN32
      fprintf (stderr, "error (GENEIGSOLVE): not positive definite %d\n",
	       info);
#else
      fprintf (stderr, "error (GENEIGSOLVE): not positive definite %ld\n",
	       info);
#endif
    }
    exit (1);
  }
}
