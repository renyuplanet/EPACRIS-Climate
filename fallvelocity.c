/*Function to calculate fall velocity of aerosol */

/* Formulation from Seinfield & Pandis (2006) Book */
/* SI */
/* Except the output velocity has been converted into CGS */


#include <math.h>
#include "constant.h"

void fallvelocity(double g);

void fallvelocity(double g)
{
	int i;
	double lambda, drag;
	
	for (i=1; i<=zbin; i++) {
		lambda=2.0*AIRVIS/pl[i]/pow(8*meanmolecular[i]*1.0E-3/PI/R_GAS/tl[i],0.5); /* Mean free path */
		/* printf("%f %f %e\n",pl[i],tl[i],lambda); */
		drag=1.0+2.0*lambda/AERSIZE*(1.257+0.4*exp(-1.1*AERSIZE/2/lambda)); /* Drag coefficient */
		VFall[i]=1.0E+2/18.0*pow(AERSIZE,2)*AERDEN*g*drag/AIRVIS;
	}

}
