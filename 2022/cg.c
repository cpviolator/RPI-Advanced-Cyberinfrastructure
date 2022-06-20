#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

/*
1. #pragma acc data present array[...]
   This pragma tells the compiler that the data in `array` is already present on 
   the device. As such, PGI++ does not have to generate implicit copies
   of `array` each time a GPU accelerator loops is compiled
2. #pragma acc parallel loop collapse(N)
   Here we tell the compiler that there are two loops in a nested structure,
   and that each iteration of the loop is independent. We say that the loop
   may be 'flattened'.
3. #pragma acc parallel loop tile(N,M)
   Here we tell the compiler that there is spatial locality in the data that
   it can exploit. for example, one thread at (x,y) must fetch data from 
   site x+1 to perform a differential operator calculation in the +ve x-dim. 
   A different thread at (x+1,y-1) needs the exact same data to perform the
   differential operator in the +ve y-dim. The tile(N,M) clause instructs
   the compiler to store NxM lots of data in fast memory cache. You must 
   experiment with N and M.
4. #pragma acc update self (array[])
   #pragma acc update device (array[])
   These pragmas assume that `array` is defined on both the device and the host
   'update self' will copy the data on the device array to the host array.
   'update device' will copy the data on the host array to the device array.
5. Data causes. The 'data' clause in a pragma indicates that we are declaring memory 
   on the GPU. Three types of clause are listed:
    
   a) copy(A[0:N])      The 'copy' clause tells the compiler that it should copy 
                        the HOST data in the A array to the GPU at the start of 
                        the caluclation, and it is also to copy the data in the 
                        GPU back to the host at the end of the calculation.
    
   b) copyin(rhs[0:N])  The 'copyin' clause tells the compiler to copy the data 
                        from the HOST to the GPU at the start, but it should not 
                        copy the data in the GPU back to the host at the end. 
                        Can you guess what the opposite clause would be? To not 
                        copy the HOST data to the GPU at the start, but to copy 
                        from the GPU to the HOST at the end?
    
   c) create(Anew[0:N]) The 'create' clause simply creates memory on the device. 
                        No data transfer is performed for this array

*/

#define FLOAT double

int n = 1024;

// y = ax + y
void axpy(FLOAT a, FLOAT *x, FLOAT *y) {  
#pragma acc data present(x[0:n]) present(y[0:n])
#pragma acc parallel loop collapse(1)
  for(int i=0; i<n; ++i) y[i] += a*x[i];
}

//z = ax + y
void axpyz(FLOAT a, FLOAT *x, FLOAT *y, FLOAT *z) {  
#pragma acc data present(x[0:n]) present(y[0:n]) present(z[0:n])
#pragma acc parallel loop collapse(1)
  for(int i=0; i<n; ++i) z[i] = a*x[i] + y[i];
}

// result = < x | y >
FLOAT dotProd(FLOAT *x, FLOAT *y) {  
#pragma acc data present(x[0:n]) present(y[0:n])
  FLOAT result = 0.0;
#pragma acc parallel loop reduction(+:result)
  for(int i=0; i<n; ++i) result += x[i]*y[i];
  return result;
}

// Perform Mx = y
void matVec(FLOAT **mat, FLOAT *y, FLOAT *x){
#pragma acc data present(mat[0:n][0:n]) present(x[0:n]) present(y[0:n])
  for(int i=0; i<n; ++i) {
    FLOAT result = 0;
#pragma acc parallel loop reduction(+:result)
    for(int j=0; j<n; ++j) {
      result += mat[i][j] * x[j];
    }
    y[i] = result;
  }
#pragma acc update device (y[0:n])
}

// L2 norm
FLOAT norm2(FLOAT *x){
#pragma acc data present(x[0:n])
  FLOAT nrm2 = 0.0;
#pragma acc parallel loop reduction(+:nrm2)
  for(int i=0; i<n; ++i) nrm2 += x[i]*x[i];
  return nrm2;
}

// zero vector
void zero(FLOAT *x){
#pragma acc data present(x[0:n])
#pragma acc parallel loop collapse(1)
  for( int i = 0; i < n; ++i ) x[i] = 0;
}

// copy vector
void copy(FLOAT *x, FLOAT *y) {
#pragma acc data present(x[0:n]) present(y[0:n])
  {
#pragma acc parallel loop collapse(1)
    for(int i=0; i<n; i++)
      x[i] = y[i];
  }
}

//https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
int cg( FLOAT **mat, FLOAT *x, FLOAT *b, double tol, int maxiter){  
  // Temp objects
  FLOAT *temp = (FLOAT*)malloc( n*sizeof(FLOAT) );
  FLOAT *res = (FLOAT*)malloc( n*sizeof(FLOAT) );
  FLOAT *p = (FLOAT*)malloc( n*sizeof(FLOAT) );
  FLOAT *Ap = (FLOAT*)malloc( n*sizeof(FLOAT) );

  bool use_init_guess = false;
  int success = 0;
  int iter = 0;
  double eps = tol*tol;
  double rsq_new, rsq;
  double alpha, beta;
  double denom;
  
#pragma acc data create(temp[0:n], res[0:n], p[0:n], Ap[0:n])
  {    
    // Find norm of rhs.
    double bnorm = norm2(b);
    double bsqrt = sqrt(bnorm);
    
    // Sanity  
    if(bsqrt == 0 || bsqrt != bsqrt) {
      printf("Error in inverterCG: inverting on zero source or Nan\n");
      exit(0);
    }   
    
    // compute initial residual
    //---------------------------------------  
    if (norm2(x) > 0.0) {
      use_init_guess = true;
      // Initial guess supplied: res = b - A*x0    
      matVec(mat, temp, x);
      zero(res);
      axpyz(-1.0, temp, b, res);
      
      // temp contains the original guess
      copy(temp, x);
      
      printf("using initial guess, |x0| = %g, |b| = %g, |res| = %g\n",
	     sqrt(norm2(temp)), bsqrt, sqrt(norm2(res)));
      
    } else {
      // No initial guess supplied. Initial residual is the source.    
      copy(res, b);
      zero(temp);
    }    
    zero(x);
    copy(p, res);
    rsq = norm2(res);
    //---------------------------------------
    
    // Iterate until convergence
    //---------------------------------------
    for (iter = 0; iter < maxiter; iter++) {
      
      // Compute Ap.
      // This part is usually the bottle neck in any program
      // and we work very hard to optimise it!
      matVec(mat, Ap, p);
      denom = dotProd(p, Ap); // This is a reduction
      alpha = rsq/denom;
      
      axpy( alpha, p,    x);
      axpy(-alpha, Ap, res);
      
      // Exit if new residual is small enough
      rsq_new = norm2(res); // This is a reduction
      printf("CG iter %d, res = %g\n", iter+1, sqrt(rsq_new));
      if (rsq_new < eps*bnorm) {
	rsq = rsq_new;
	break;
      }
      
      // Update vec using new residual
      beta = rsq_new/rsq;
      rsq = rsq_new;
      
      axpyz(beta, p, res, p);
    } // End loop over iter
    //---------------------------------------
    
    if(iter == maxiter) {
      // Failed convergence 
      printf("CG: Failed to converge iter = %d, rsq = %.16e\n", iter+1, rsq); 
      success = 0; 
    } else {
      // Convergence 
      success = iter+1; 
    }
    
    // We must add back the exact part if using an initial guess
    if(use_init_guess) axpy(1.0, temp, x);  
    // x now contains the solution to the RHS b.
    
    //Sanity
    printf("source norm = %g\n", sqrt(norm2(b)));
    printf("sol norm = %g\n", sqrt(norm2(x)));
    
    matVec(mat, temp, x);
    axpyz(-1.0, temp, b, res);
    double truersq = norm2(res);
    printf("CG: Converged iter = %d, rsq = %.16e, truersq = %.16e\n", iter+1, rsq, truersq/(bsqrt*bsqrt));
  }
  return success;  
}

int main( int argc, char* argv[] ){

  int maxiter = 1000; // The maximal allowed CG iterations  
  double tol = 1e-20; // The tolerance on the quality of the solution vector  
  double diag = 150;  // We add this value to the diagonal elements
  // of the matrix to ensure it is positive semi-definite, i.e.
  // all eigenvalues are real and greater than zero. If you change
  // N and the CG fails to converge, you need to increase this
  // value.

  if( argc == 5 ) {
    n = atoi( argv[1] );
    maxiter = atoi( argv[2] );
    tol = atof( argv[3] );
    diag = atoi( argv[4] );
  } else {
    printf("Error in command line params\n");
    exit(0);
  }

  // Starting guess
  FLOAT *x = (FLOAT*)malloc( n*sizeof(FLOAT) );
  
  // For the matrix M, find x = M^{-1} * b
  FLOAT *b = (FLOAT*)malloc( n*sizeof(FLOAT) );  
  b[0] = 1.0;

  // rands
  FLOAT **rands = (FLOAT**)malloc( n*sizeof(FLOAT*) );
  for(int i=0; i<n; i++) {
    rands[i] = (FLOAT*)malloc( n*sizeof(FLOAT) );
    for(int k=0; k<n; k++) {
      rands[i][k] = drand48();
    }
  }
  
  // The row major matrix, enforced orthogonality
  FLOAT **mat = (FLOAT**)malloc( n*sizeof(FLOAT*) );
  for(int i=0; i<n; i++) {
    mat[i] = (FLOAT*)malloc( n*sizeof(FLOAT) );
    for(int k=0; k<n; k++) {
      mat[i][k] = rands[i][k] + rands[k][i];
      if(i==k) mat[i][k] += diag;
    }
  }
  
  // OpenACC Init 
#pragma acc init
  
  // Invert on the GPU
  FLOAT time = -(FLOAT)clock();
#pragma acc data copyin(mat[0:n][0:n], b[0:n]) copy(x[0:n])
  {
    cg(mat, x, b, tol, maxiter);
  }  
  time += (FLOAT)clock();
  printf( "GPU CG invert time %g sec\n", time*1e-6 );

  
  // Perturb the matrix, use previous solution as initial guess
  for(int k=0; k<n; k++) {
    mat[k][k] += 1e-4;
  }
  // Pertub the guess, use the same matrix.
  //x[0] =- 0.000001;
#pragma acc data copyin(mat[0:n][0:n], b[0:n]) copy(x[0:n])
  {
    cg(mat, x, b, tol, maxiter);  
  }  
  time += (FLOAT)clock();
  printf( "GPU CG invert with initial guess time %g sec\n", time*1e-6 );  
} 
