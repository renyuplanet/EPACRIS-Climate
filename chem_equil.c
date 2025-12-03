/*----------------------- chem_equil.c ---------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified:

------------------------------------------------------------------ */

/* Procedure to calculate chemical equilibrium.  Based on 
White, Johnson & Dantzig 1958

Inputs: P - Pressure
        g0 - Gibb's free energy of formation / RT
	a_ij - Matrix defining the atomic makeup of each molecule
	b_j - Relative abundances of each atom
        y - Initial solution for relative abundances.  This is
	    replaced by the actual solution upon output back to
	    main_chem.c
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "include.h"
#include "constant.h"
#include "nrutil.h"

#include "gaussj.c"

/* --- Global variables ------------------------------------------ */

/* --- Function prototypes --------------------------------------- */

 void gaussj(double **a, int n, double **b, int m);              

/* ------- begin ---------------- main --------------------------- */


void chem_equil(double P, double *g0, int **a_ij, double *b_j, double *y)
{
  double *x, *f_y, *delta, *pi, *b_test, y_sum, x_sum, 
    f_sum, delta_sum, delta_max, u, lambda, df_dlambda, tiny;
  double **matrix, **solution;
  int i, j, k, dum, count, max_count;
  FILE *file;

  count = 0;
  max_count = 100000;
  tiny = 1.0e-90;
  P *= 1.0/ATM_TO_PA;

/* Allocate memory */

  solution = dmatrix(0, NATOMS, 0, NATOMS);
  matrix = dmatrix(0, NATOMS, 0, NATOMS);
  x = dvector(0, NMOLECULES-1);
  f_y = dvector(0, NMOLECULES-1);
  delta = dvector(0, NMOLECULES-1);
  pi = dvector(0, NATOMS-1);
  b_test = dvector(0, NATOMS-1);

/* Start iteration scheme */
  
  do{

/* Calculate y_sum, f(y), and F(Y) */

    y_sum = 0.0;
    f_sum = 0.0;
    for(i=0; i<NMOLECULES; i++)
      y_sum += y[i];
    
    for(i=0; i<NMOLECULES; i++){
      f_y[i] = y[i]*(g0[i]+log(P)+log(y[i]/y_sum));
      f_sum += f_y[i];
    }
    
/* Set up matrix equations
                      --  later can save computing time by taking
                          advantage of the fact that matrix is symmetric*/
    
    for (j=0; j<NATOMS; j++){
      for (k=0; k<NATOMS; k++){
	matrix[j][k] = 0.0;
	for (i=0; i<NMOLECULES; i++)
	  matrix[j][k] += a_ij[i][j]*a_ij[i][k]*y[i];
      }
    }
    
    for (j=0; j<NATOMS; j++){
      b_test[j] = 0.0;
      for (i=0; i<NMOLECULES; i++)
	b_test[j] += y[i]*a_ij[i][j];
    }
	

    for (j=0; j<NATOMS; j++){
      matrix[NATOMS][j] = b_test[j];
      matrix[j][NATOMS] = b_test[j];
    }
    
    matrix[NATOMS][NATOMS] = 0.0;
    
    for (j=0; j<NATOMS+1; j++){
      solution[j][0] = 0.0;
      for (i=0; i<NMOLECULES; i++){
	if(j < NATOMS){
	  solution[j][0] += a_ij[i][j]*f_y[i];
	} else {
	  solution[j][0] += f_y[i];
	}
      }
    }

/* Solve matrix equation  -- later switch to a less computationally
                             intense solution */

    gaussj(matrix, NATOMS+1, solution, 1);
    
    for (j=0; j<NATOMS; j++)
      pi[j] = solution[j][0];
    u = solution[NATOMS][0];
    
/* Determine x_i */

    x_sum = (u+1)*y_sum;
    
    for (i=0; i<NMOLECULES; i++){
      x[i] = -f_y[i] + (y[i]/y_sum)*x_sum;
      for (j=0; j<NATOMS; j++)
	x[i] += pi[j]*a_ij[i][j]*y[i];
    }

/* Determine lambda */

    delta_sum = 0.0;
    for (i=0; i<NMOLECULES; i++){
      delta[i] = x[i] - y[i];
      delta_sum += delta[i];
    }

    if (y[0] == tiny)
      delta_max = 0.0;
    else
      delta_max = fabs(delta[0]/y[0]);
    for (i=1; i<NMOLECULES; i++)
      if (y[i] > tiny)
	delta_max = MAX(delta_max, fabs(delta[i]/y[i]));
    
    lambda = -1.0;
    for (i=0; i<NMOLECULES; i++){
      if(y[i]/delta[i] > lambda && delta[i] < 0.0 && y[i] > tiny)
	lambda = y[i]/delta[i];
    }
    if(lambda != -1.0 && lambda < -0.0001)
      lambda += 0.0001;
    if (count == 0 && lambda < -0.0001)
      lambda = -0.0001;
    
    df_dlambda = 0.0;
    for (i=0; i<NMOLECULES; i++)
      df_dlambda += delta[i] * (g0[i]+log(P)+
		    log((y[i]+lambda*delta[i])/(y_sum+lambda*delta_sum)));

/* Determine new values of y_i for next iteration */

    for (i=0; i<NMOLECULES; i++){
      y[i] = y[i] - lambda*delta[i];
      if (y[i] < tiny)
	y[i] = tiny;
    }

    count++;
    
  }while(count < max_count && (delta_max >= 0.001 || lambda+1.0 > 0.000001));

/* Check that convergence is achieved and mass-balance is preserved. */
  if(count == max_count){
    printf("\nWARNING: Convergence not achieved\n"); 
    printf("lambda %f\n", lambda);
  }

  for (j=0; j<NATOMS; j++){
    b_test[j] = 0.0;
    for (i=0; i<NMOLECULES; i++)
      b_test[j] += y[i]*a_ij[i][j];
    if (fabs((b_test[j]-b_j[j])/b_j[j]) > 0.1 && b_j[j]>1E-20 ){
      printf("ERROR: Mass balance not conserved\n"); 
      printf("%e %e %e %d\n", b_test[j], (b_test[j]-b_j[j])/b_j[j], b_j[j], j);
      exit(1);
    }
  }

/* Clean up */

  free_dmatrix(solution, 0, NATOMS, 0, NATOMS);
  free_dmatrix(matrix, 0, NATOMS, 0, NATOMS);
  free_dvector(x, 0, NMOLECULES-1);
  free_dvector(f_y, 0, NMOLECULES-1);
  free_dvector(delta, 0, NMOLECULES-1);
  free_dvector(pi, 0, NATOMS-1);
  free_dvector(b_test, 0, NATOMS-1);
  
  return;
}
