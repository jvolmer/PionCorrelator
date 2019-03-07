#ifndef MATRIX_H
#define MATRIX_H

// Allocate -------------------------------------------------------------
int *AllocIntVector (int n);
double *AllocDoubleVector (int n);
_Complex double *AllocComplexVector (int n);
int **AllocIntMatrix (int m, int n);
_Complex double ***Alloc3D(int a_max, int b_max, int c_max);
_Complex double ****Alloc4D(int n_max, int mu_max, int a_max, int b_max);
_Complex double *****Alloc5D(int a_max, int b_max, int c_max, int d_max, int e_max);

// Fill -----------------------------------------------------------------
_Complex double ***setZero3D(_Complex double ***A, int N);
_Complex double ***setRandom3D(_Complex double ***A, int N);

// Calculate ------------------------------------------------------------
double normSquared(_Complex double ***a, int m);
_Complex double scalarProductComplex(_Complex double ***a, _Complex double ***b, int m);

// Print ----------------------------------------------------------------

void fPrint3D(FILE *f, _Complex double ***A, int b_max, int c_max, int d_max);

#endif
