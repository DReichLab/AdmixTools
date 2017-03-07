//==============================================================================
//
// Linear System Solution by Gauss method
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================

#include <stdio.h>
#include <math.h>
#include <nicklib.h>

//==============================================================================
// return 1 if system not solving
// nDim - system dimension
// pfMatr - matrix with coefficients
// pfVect - vector with free members
// pfSolution - vector with system solution
// pfMatr becames trianglular after function call
// pfVect changes after function call
//
// Developer: Henry Guennadi Levkin
//
//==============================================================================
int LinearEquationsSolving(int nDim, double* pfMatr, double* pfVect, double* pfSolution)
{
  double fMaxElem;
  double fAcc;

  int i , j, k, m;


  for(k=0; k<(nDim-1); k++) // base row of matrix
  {
    // search of line with max element
    fMaxElem = fabs( pfMatr[k*nDim + k] );
    m = k;
    for(i=k+1; i<nDim; i++)
    {
      if(fMaxElem < fabs(pfMatr[i*nDim + k]) )
      {
        fMaxElem = pfMatr[i*nDim + k];
        m = i;
      }
    }
    
    // permutation of base line (index k) and max element line(index m)
    if(m != k)
    {
      for(i=k; i<nDim; i++)
      {
        fAcc               = pfMatr[k*nDim + i];
        pfMatr[k*nDim + i] = pfMatr[m*nDim + i];
        pfMatr[m*nDim + i] = fAcc;
      }
      fAcc = pfVect[k];
      pfVect[k] = pfVect[m];
      pfVect[m] = fAcc;
    }

    if( pfMatr[k*nDim + k] == 0.) return 1; // needs improvement !!!

    // triangulation of matrix with coefficients
    for(j=(k+1); j<nDim; j++) // current row of matrix
    {
      fAcc = - pfMatr[j*nDim + k] / pfMatr[k*nDim + k];
      for(i=k; i<nDim; i++)
      {
        pfMatr[j*nDim + i] = pfMatr[j*nDim + i] + fAcc*pfMatr[k*nDim + i];
      }
      pfVect[j] = pfVect[j] + fAcc*pfVect[k]; // free member recalculation
    }
  }

  for(k=(nDim-1); k>=0; k--)
  {
    pfSolution[k] = pfVect[k];
    for(i=(k+1); i<nDim; i++)
    {
      pfSolution[k] -= (pfMatr[k*nDim + i]*pfSolution[i]);
    }
    pfSolution[k] = pfSolution[k] / pfMatr[k*nDim + k];
  }

  return 0;
}
//==============================================================================
// testing of function
//==============================================================================

#define MATRIX_DIMENSION 4

int main(int nArgs, char** pArgs)
{
  int nDim = MATRIX_DIMENSION;
  double fMatr[MATRIX_DIMENSION*MATRIX_DIMENSION] =
  {
  1.0,  2.0, -1.0, -2.0,
  1.0,  3.0, -1.0, -2.0,
  2.0,  1.0,  1.0,  1.0,
  3.0,  1.0,  2.0,  1.0,
  };
  double fVec[MATRIX_DIMENSION] = {-6.0, -4.0, 11.0, 15.0};
  
  double fSolution[MATRIX_DIMENSION];
  int res;
  int i;

  printf("Using nicklib\n") ;
  res = linsolv(nDim, fMatr, fVec, fSolution); // !!!
  printf("Solution:\n");
  printmat(fSolution, 1, nDim);

  printf("Using Linear...\n") ;

  res = LinearEquationsSolving(nDim, fMatr, fVec, fSolution); // !!!

  if(res)
  {
    printf("No solution!\n");
    return 1;
  }
  else
  {
    printf("Solution:\n");
    printmat(fSolution, 1, nDim);

  }

  
  return 0;
}
