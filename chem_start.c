/*----------------------- chem_equil.c ---------------------------

Author: Eliza Miller-Ricci (emillerricci@cfa.harvard.edu)
Last modified:

------------------------------------------------------------------ */

/* Starting solution for calculating chemical equilibrium

Inputs: a_ij - Matrix defining the atomic makeup of each molecule
        b_j - Relative abundances of each atom
        y - Initial solution for relative abundances.  This is
	    filled in by chem_start and then returned to 
	    main_chem.c.
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "include.h"
#include "nrutil.h"

/* --- Global variables ------------------------------------------ */

/* --- Function prototypes --------------------------------------- */

/* ------- begin ---------------- main --------------------------- */


void chem_start(int **a_ij, double *b_j, double *y)
{
  double *b_test, small;
  int *sum, a_sum, i, j, k;

/* Allocate memory */
  
  sum = ivector(0, NATOMS-1); 
  b_test = dvector(0, NATOMS-1);

/* Partition atoms into the various molecules */
    
    for (j=0; j<NATOMS; j++) {
        /* printf("%s %d %e\n", "b_j", j, b_j[j]); */
    }

  for (j=0; j<NATOMS; j++){
    sum[j] = 0;
	  for (i=0; i<NMOLECULES; i++) {
          sum[j] += a_ij[i][j];
          /* printf("%s %d %d %d\n", "a_ij_new", i, j, a_ij[i][j]);*/
      }
    
      //printf("%s %d %d\n", "sum", j, sum[j]);
      
	  if (j == 0) { 
		  small = b_j[0]/sum[0];
		  /* small = MIN(b_j[0]/sum[0], b_j[1]/sum[1]); */
	  } else {
		  small = MIN(small, b_j[j]/sum[j]);
	  }

  }

    for (i=0; i<NMOLECULES; i++) {
        y[i] = small;
    }
    
    //printf("%s %e\n", "small", small);
  
  for (j=0; j<NATOMS; j++){
    b_test[j] = 0.0;
    for (i=0; i<NMOLECULES; i++){
      b_test[j] += a_ij[i][j]*y[i];
    }
  }

  for (j=0; j<NATOMS; j++){
    for (i=0; i<NMOLECULES; i++){
      a_sum = 0;

      for (k=0; k<NATOMS; k++){
      a_sum += a_ij[i][k];
      }
      
      if(a_sum == 1 && a_ij[i][j] == 1){
	y[i] += b_j[j] - b_test[j];
      }
							
    }
  }
  
/* Check that mass balance is preserved */

  for (j=0; j<NATOMS; j++){
    b_test[j] = 0.0;
    for(i=0; i<NMOLECULES; i++)
      b_test[j] += a_ij[i][j]*y[i];
    if ((b_test[j]-b_j[j])/b_j[j] > 0.0001){
      printf("ERROR: Mass balance not conserved\n"); 
      printf("%f %f\n", b_test[j], b_j[j]);
      exit(1);
    }
  }

/* Make sure that all of the initial values are positive */

  for (i=0; i<NMOLECULES; i++)
    if (y[i] < -1.0E-10){
      printf("%s\t%e\t%d\n","ERROR: Negative initial abundance",y[i],i);
      exit(1);
    } else if (y[i]<0.0) {
		y[i] = fabs(y[i]);
	}


/* Clean up */
  
  free_ivector(sum, 0, NATOMS-1);
  free_dvector(b_test, 0, NATOMS-1);

  return;
}

