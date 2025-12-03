/*----------------------- chemequil.c ----------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: February 10, 2010

------------------------------------------------------------------ */

/* A function to calculate chemical equilibrium for specied T and P.
Based on methods of White, Johnson & Dantzig 1958; and code of Eliza Miller-Ricci.
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "include.h"
#include "constant.h"
#include "nrutil.h"
#include "chem_start.c"
#include "chem_equil.c"
#include "config.h"

/* --- Global variables ------------------------------------------ */


/* --- Function prototypes --------------------------------------- */

void chem_start(int **a_ij, double *b_j, double *y);
void chem_equil(double P, double *g0, int **a_ij, double *b_j, double *y);

/* ------- begin ---------------- main --------------------------- */

void chemquil(double P[], double T[], int NL, int species[], int NS, double **mixequil, char atomfile[])
{
  double *g0, *y, *y_start, *a, *b, *c, *d, *e, *b_j, norm, **ysave;
  char molecule[NMOLECULES][8], atom[NATOMS][3];
  int **a_ij, i, j, k, l, dum, *id;
  FILE *file;

/* Allocate memory */

  a_ij = imatrix(0, NMOLECULES-1, 0, NATOMS-1);
  g0 = dvector(0, NMOLECULES-1);
  y = dvector(0, NMOLECULES-1);
  ysave = dmatrix(0, NL-1, 0, NMOLECULES-1);
  y_start = dvector(0, NMOLECULES-1);
  a = dvector(0, NMOLECULES-1);
  b = dvector(0, NMOLECULES-1);
  c = dvector(0, NMOLECULES-1);
  d = dvector(0, NMOLECULES-1);
  e = dvector(0, NMOLECULES-1);
  b_j = dvector(0, NATOMS-1);
  id = ivector(0, NMOLECULES-1);

  /* Read in relative abundances of each element */

  file = fopen(atomfile, "r");  
  for(j=0; j<NATOMS; j++){
    fscanf(file, "%d %s", &dum, atom[j]);
    fscanf(file, "%le\n", &b_j[j]);
  }
  fclose(file);
  
  // DEBUG: Print Fe abundance
  //printf("DEBUG: Fe abundance from file: %e (atom index %d)\n", b_j[20], 20);
  //printf("DEBUG: Fe abundance from file: %e (atom index %d)\n", b_j[21], 21);  
  for (j=0; j<NATOMS; j++){
    if (b_j[j] < 1.0e-20)
      b_j[j] = 1.0e-20;
  }

/* Read in molecular data --  
       g0 is initially in units of calories*mol^-1  */

  file = fopen(MOL_DATA_FILE, "r");

    for(i=0; i<NATOMS+7; i++){
        fscanf(file, "%s", molecule[0]);
    }/* The header of the file */

  for(i=0; i<NMOLECULES; i++){

    fscanf(file, "%s %d %le %le %le %le %le", molecule[i], &id[i], &a[i], &b[i], &c[i], 
	   &d[i], &e[i]);
      /*printf("%s %d %e %e %e %e %e\n", molecule[i], id[i], a[i], b[i], c[i],
             d[i], e[i]);*/
      for (j=0; j<NATOMS; j++){
          fscanf(file, "%d", &a_ij[i][j]);
          /*printf("%s %d %d %d %d\n", "a_ij", i, j, a_ij[i][j]);*/
      }
      

  }
    
  fclose(file);

  /* Chemical equilibrium calculation and return results*/

  /* Initial solution */

  chem_start(a_ij, b_j, y_start);
  

	
	//printf("%s\n","start chemequil");

  /* Loop over P-T */

  for (l=1; l<NL; l++){
      
      for (i=0; i<NMOLECULES; i++) {
          y[i] = y_start[i];
      }

      for (i=0; i<NMOLECULES; i++){
          g0[i] = a[i]/T[l] + b[i] + c[i]*T[l] + d[i]*SQ(T[l]) + e[i]*CUBE(T[l]);
          g0[i] *= CAL_TO_JOULE / (R_GAS * T[l]);
      }

      chem_equil(P[l], g0, a_ij, b_j, y);
	  
//ms2022	  printf("%s %d\n","compute chemequil for layer", l);
      
      /* Normalize the y's to 1 */
      
      norm = 0.0;
      for (i=0; i<NMOLECULES; i++) {norm += y[i];}
      
      for (i=0; i<NMOLECULES; i++) {
		  ysave[l][i] = y[i]/norm; 
		  /* printf("E %lf %lf %d %e\n", P[l], T[l], i, ysave[l][i]); */
	  }
      


  }
	
	//printf("%s\n","*** Chemical equilibrium achieved ***");
  
  	  /* Return the requested mixing ratios */
 
      for (i=0; i<NMOLECULES; i++){
          k=0;
          for (j=1; j<=NS; j++) {
              if (species[j]==id[i]) {k=j;}
		  }
          if (k!=0) {
			  for (l=1; l<NL; l++) {mixequil[l][k]=ysave[l][i];}
              

          } else {
			  for (l=1; l<NL; l++) {mixequil[l][k] = 0.0;}
		  }

      }
	
	//printf("%s\n","*** Chemical equilibrium achieved ***");

/* Clean up */
  
  free_imatrix(a_ij, 0, NMOLECULES-1, 0, NATOMS-1);
  free_dvector(g0, 0, NMOLECULES-1);
  free_dvector(y, 0, NMOLECULES-1);
  free_dvector(a, 0, NMOLECULES-1);
  free_dvector(b, 0, NMOLECULES-1);
  free_dvector(c, 0, NMOLECULES-1);
  free_dvector(d, 0, NMOLECULES-1);
  free_dvector(e, 0, NMOLECULES-1);
  free_dvector(y_start, 0, NMOLECULES-1);
  free_dvector(b_j, 0, NATOMS-1);
  free_ivector(id, 0, NMOLECULES-1);
  free_dmatrix(ysave, 0, NL-1, 0, NMOLECULES-1);
  
  return;
}

