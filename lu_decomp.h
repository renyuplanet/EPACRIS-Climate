/**
 * lu_decomp.h
 * 
 * Header file for LU decomposition functions
 * Used for solving linear systems via LU decomposition
 */

#ifndef LU_DECOMP_H
#define LU_DECOMP_H

// LU decomposition: decomposes matrix a into L*U form
// a: input matrix (modified in-place), n: matrix size
// indx: output permutation vector, d: output determinant sign
void ludcmp(double **a, int n, int *indx, double *d);

// LU back substitution: solves linear system after LU decomposition
// a: LU-decomposed matrix, n: matrix size
// indx: permutation vector from ludcmp, b: right-hand side (modified to solution)
void lubksb(double **a, int n, int *indx, double b[]);

#endif // LU_DECOMP_H

