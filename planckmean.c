/*Function to calculate the infrared Planck mean opacity at each level */
/* input cross section in m^2 */
/* output mean cross section in cm^2 */
/* also output the cross section in cm^2 at each level in each wavelength */

#include <math.h>
#include "constant.h"

void planckmean(double Mean[], double MeanS[], double **xsc);

void planckmean(double Mean[], double MeanS[], double **xsc)
{		
    int i, j, k, s;
    double Fplanck[zbin+1], Cplanck[zbin+1], wave[NLAMBDA], h1, h2;
    double Fstar[zbin+1], Cstar[zbin+1];
    
    /* get wave */
    for (i=0; i<NLAMBDA; i++) {
        wave[i] = wavelength[i]*1.0E-9; /* convert to m */
    }
    
	/* Compute the total planck flux in this range */
	for (i=1; i<=zbin; i++) {
		Fplanck[i] = 0.0;
        Fstar[i] = 0.0;
		for (j=0; j<NLAMBDA-1; j++) {
			Fplanck[i] += 2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1)*(wave[j+1]-wave[j]);
            Fstar[i] += solar[j]*(wave[j+1]-wave[j]);
		}
	}
	
	/* Calculate the cross section */
	for (i=1; i<=zbin; i++) {
		Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
		for (j=0; j<NLAMBDA-1; j++) {
			h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
			h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
			if (fabs(xsc[i][j+1]-xsc[i][j])<1E-20) {
				Cplanck[i] += (wave[j+1]-wave[j])/2.0*(xsc[i][j]*h1+xsc[i][j+1]*h2);
			} else {
				Cplanck[i] += (wave[j+1]-wave[j])*(xsc[i][j+1]*h2-xsc[i][j]*h1)/(log(xsc[i][j+1]*h2*1.0E+30)-log(xsc[i][j]*h1*1.0E+30));
			}
		}
        for (j=0; j<NLAMBDA-1; j++) {
			h1=solar[j];
			h2=solar[j+1];
			if (fabs(h2*xsc[i][j+1]-h1*xsc[i][j])<1E-20) {
				Cstar[i] += (wave[j+1]-wave[j])/2.0*(xsc[i][j]*h1+xsc[i][j+1]*h2);
			} else {
				Cstar[i] += (wave[j+1]-wave[j])*(xsc[i][j+1]*h2-xsc[i][j]*h1)/(log(xsc[i][j+1]*h2 + 2.0E-80)-log(xsc[i][j]*h1 + 1.0E-80));
			}
		}
        Mean[i] = Cplanck[i]/Fplanck[i];
        MeanS[i] = Cstar[i]/Fstar[i];
        /*printf("%e %e %e\n", pl[i], Mean[i], MeanS[i]);*/
	}

}
