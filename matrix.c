#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "matrix.h"

#define PRECISION 1000000000000

// Allocate -------------------------------------------------------------

int *AllocIntVector (int n) {
  /* create n Vector */
  int *b;
  int i;
  b = malloc(n * sizeof(int));
  for (i=0; i<n; i++) {
    b[i] = 0.;
  }
  return b;
}

_Complex double *AllocComplexVector (int n) {
  /* create n Vector */
  _Complex double *b;
  int i;
  b = malloc(n * sizeof(_Complex double));
  for (i=0; i<n; i++) {
    b[i] = 0. + 0. * I;
  }
  return b;
}

int **AllocIntMatrix (int m, int n) {
  int **A;
  int i,j;
  A = malloc(m * sizeof(*A) +  m * n * sizeof(**A));
  for (i=0; i<m; i++) 
    A[i] = ((int *)(A + m)) + n*i;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      A[i][j] = 0;
    }
  }
  return A;
}

_Complex double ***Alloc3D(int a_max, int b_max, int c_max) {
  _Complex double ***A;
  int a, b, c, d;
  /* alloc space for n_max 3-pointers + ... */
  A = malloc(a_max * sizeof(*A) + 
	     a_max * b_max * sizeof(**A) + 
	     a_max * b_max * c_max * sizeof(***A));
  for (a=0; a<a_max; a++){
    A[a] = ((_Complex double **)(A + a_max)) + b_max*a;
    for (b=0; b<b_max; b++){
      A[a][b] = ((_Complex double *)(A + a_max + b_max*a_max)) + c_max*b_max*a + c_max*b;
      for (c=0; c<c_max; c++){
	A[a][b][c] = 0. + 0. * I;
      }
    }
  }
  return A;
}

_Complex double ****Alloc4D(int a_max, int b_max, int c_max, int d_max) {
  _Complex double ****A;
  int a, b, c, d;
  /* alloc space for n_max 3-pointers + ... */
  A = malloc(a_max * sizeof(*A) +
	     a_max * b_max * sizeof(**A) +
	     a_max * b_max * c_max * sizeof(***A) +
	     a_max * b_max * c_max * d_max * sizeof(****A));
  for (a=0; a<a_max; a++){
    A[a] = ((_Complex double ***)(A + a_max)) + b_max*a;
    for (b=0; b<b_max; b++){
      A[a][b] = ((_Complex double **)(A + a_max + b_max*a_max)) + c_max*b_max*a + c_max*b;
      for (c=0; c<c_max; c++){
	A[a][b][c] = ((_Complex double *)(A + a_max + b_max*a_max + c_max*b_max*a_max)) + d_max*c_max*b_max*a + d_max*c_max*b + d_max*c;
	for (d = 0; d<d_max; d++) {
	  A[a][b][c][d] = 0. + 0. * I;
	}
      }
    }
  }
  return A;
}


// Fill -----------------------------------------------------------------

_Complex double ***setZero3D(_Complex double ***A, int N) {
  int n, al, a;
  for (n=0; n<N; n++) {
    for (al=0; al<4; al++) {
      for (a=0; a<3; a++) {
	A[n][al][a] = 0. + 0. * I;
      }
    }
  }
  return A;
}

_Complex double ***setRandom3D(_Complex double ***A, int N) {
  int n, al, a;
  double x, y;
  for (n=0; n<N; n++) {
    for (al=0; al<4; al++) {
      for (a=0; a<3; a++) {
	x = (rand() * rand()) % PRECISION;
	x /= PRECISION;
	y = (rand() * rand()) % PRECISION;
	y /= PRECISION;
	A[n][al][a] = x + y * I;
      }
    }
  }
  return A;
}


// Calculate ------------------------------------------------------------

double normSquared(_Complex double ***vec, int m){
  int i,f,a,al;
  double kahan_rest = 0.;
  double result = 0.;
  double kahan_new, kahan_add;
    for (i=0; i<m; i++){
      for (al=0; al<4; al++){
	for (a=0; a<3; a++){
	  kahan_new = (conj(vec[i][al][a]) * vec[i][al][a]) 
	    - kahan_rest;
	  kahan_add = result + kahan_new;
	  kahan_rest = (kahan_add - result) - kahan_new;
	  result = kahan_add;
	}
      }
    }
  return result;
}

_Complex double scalarProductComplex(_Complex double ***vec1, _Complex double ***vec2, int m){
  int i,f,a,al;
  _Complex double kahan_rest = 0. + 0. * I;
  _Complex double result = 0. + 0. * I;
  _Complex double kahan_new, kahan_add;
    for (i=0; i<m; i++){
      for (al=0; al<4; al++){
	for (a=0; a<3; a++){
	  kahan_new = (conj(vec1[i][al][a]) * vec2[i][al][a]) 
	              - kahan_rest;
	  kahan_add = result + kahan_new;
	  kahan_rest = (kahan_add - result) - kahan_new;
	  result = kahan_add;
	}
      }
    }
  return result;
}


// Print ----------------------------------------------------------------

void fPrint3D(FILE *f, _Complex double ***A, int b_max, int c_max, int d_max) {
  int b, c, d, pr=0;

  for (b=0; b<b_max; b++) {
      for (c=0; c<c_max; c++) {
	for (d=0; d<d_max; d++) {


	  if ((creal(A[b][c][d]) != 0.) || (cimag(A[b][c][d]) != 0.)) {
	    fprintf(f, "A[n=%i][alpha=%i][a=%i] = %.20f +i*%.20f\n",
		    b, c, d, creal(A[b][c][d]), cimag(A[b][c][d]));
	  }

	}
      }
    }

}
