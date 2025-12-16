/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator

  for(size_t i = 0; i < *la; ++i)
  {
    for(size_t j = 0; j < *kv; ++j)
      AB[i * (*lab) + j] = 0.0;
  }

  for(size_t i = 1; i < *la; ++i)
  {
    AB[i * (*lab) + *kv] = -1.0;
  }

  for(size_t i = 0; i < *la; ++i)
  {
    AB[i * (*lab) + (*kv + 1)] = 2.0;
  }

  for(size_t i = 0; i < (*la - 1); ++i)
  {
    AB[i * (*lab) + (*lab - 1)] = -1.0;
  } 

}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  
  for(size_t i = 0; i < *la; ++i)
  {
    for(size_t j = 0; j < (*kv + *lab); ++j)
      AB[i * (*lab) + j] = j == (*lab - 2) ? 1.0 : 0.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector

  RHS[0] = (*BC0);
  RHS[(*la - 1)] = (*BC1);

  for(size_t i = 1; i < (*la - 1); ++i)
    RHS[i] = 0.0;

}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D

  for(size_t i = 1; i < (*la - 1); ++i)
    EX_SOL[i] = *BC0 + i *(*BC1 - *BC0);
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]

  const double dx = 1.0/((*la) + 1);
  x[0] = dx;

  for(size_t i = 1; i < (*la); ++i)
    x[i] = x[i - 1] + dx;

}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)

  for(size_t i = 0; i < (*la); ++i)
    y[i] = x[i] - y[i];

  return cblas_dnrm2(*la, y, 1)/cblas_dnrm2(*la, x, 1);
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return i * (*lab) + j;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices

  for(size_t k = 0; k < (*n - 1); ++k)
  {
    for(size_t i = k + 1; i < *n; ++i)
    {
      AB[indexABCol(*ku + i - k, k, lab)] = AB[indexABCol(*ku + i - k, k, lab)]/AB[indexABCol(*ku , k, lab)];
    }

    for(size_t i = k + 1; i < (*n); ++i)
    {
      for(size_t j = k + 1; j < (*n); ++j)
      {
        AB[indexABCol(*ku + i - j, j, lab)] = AB[indexABCol(*ku + i - j, j, lab)] - AB[indexABCol(*ku + i - k, k, lab)]/AB[indexABCol(*ku + k - j, j, lab)];
      }
    }
  }

  return *info;
}
