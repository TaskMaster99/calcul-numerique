/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#define PI 3.14

void eig_poisson1D(double* eigval, int *la){
  // TODO: Compute all eigenvalues for the 1D Poisson operator
  const double h = PI/(*la +1);
  for(size_t k = 0; k < *la ; ++k)
  {
    eigval[k] = 2.0 - 2.0 * cos( k * PI*h);
  }
}

double eigmax_poisson1D(int *la){
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
  const double h = PI/(*la +1);
  return 2.0 - 2.0 * cos(h * (*la + 1)) ;
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
 const double h = PI/(*la +1);
  return 2.0 - 2.0 * cos(0.0 * h) ;
}

double richardson_alpha_opt(int *la){
  // TODO: Compute alpha_opt
  const double eig = eigmax_poisson1D(la) + eigmin_poisson1D(la);
  return 2.0 / eig;
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat


  double* r = (double *) calloc(*la, sizeof(double));
  
  size_t it = 0;

  const double norm_b = 1/cblas_dnrm2(*la, RHS, 1);
  cblas_dcopy(*la, RHS, 1, r, 1);
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *lab, *kl, *ku, -1, AB, *la, X, 1, 1,r, 1);

  do
  {
    cblas_daxpy(*la, *alpha_rich, r, 1, X, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *lab, *kl, *ku, -1, AB, *la, X, 1, 1,r, 1);

    cblas_dcopy(*la, RHS, 1, r, 1);

    resvec[it] = cblas_dnrm2(*la, r, 1) * norm_b;

    if(resvec[it] < *tol)break;
    
    it++;
  }while(it < *maxit);

  free(r);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
}

