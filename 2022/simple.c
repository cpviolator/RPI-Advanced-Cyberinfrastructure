#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define FLOAT float

void vecfill( FLOAT *a, FLOAT *b, int n ){
#pragma present(a[0:n]) present(b[0:n]) 
  for( int i = 0; i < n; ++i ) {
    a[i] = (FLOAT)(i+1);
    b[i] = (FLOAT)(10*i);
  }
}

void vecaddgpu( FLOAT *r, FLOAT *a, FLOAT *b, int n ){
#pragma present(a[0:n]) present(b[0:n]) present(r[0:n])
  for( int i = 0; i < n; ++i ) r[i] = a[i] + b[i];
}

void vecmultgpu( FLOAT *r, FLOAT *a, FLOAT *b, int n ){
#pragma present(a[0:n]) present(b[0:n]) present(r[0:n])
  for( int i = 0; i < n; ++i ) r[i] = a[i] * b[i];
}

int main( int argc, char* argv[] ){

  int n; /* vector length */
  int k; /* operation repetitions */
  FLOAT * a; /* input vector 1 */
  FLOAT * b; /* input vector 2 */
  FLOAT * r; /* output vector */
  int i;

  if( argc > 2 ) {
    n = atoi( argv[1] );
    k = atoi( argv[2] );
  } else {
    n = 100000; /* default vector length */
    k = 1; /* one repetition */
  }

  a = (FLOAT*)malloc( n*sizeof(FLOAT) );
  b = (FLOAT*)malloc( n*sizeof(FLOAT) );
  r = (FLOAT*)malloc( n*sizeof(FLOAT) );
  
  // compute and data transfer on the GPU
  FLOAT time = -(FLOAT)clock();
#pragma acc data create(a[0:n], b[0:n], r[0:n])
  {
    vecfill( a, b, n );
    for (i=0; i<k; i++) {
      vecmultgpu( r, a, b, n );
      vecaddgpu( r, a, b, n );
    }
  }
  time += (FLOAT)clock();
  
  printf( "Time %g sec\n", time*1e-6 );
  return 0;
} 
