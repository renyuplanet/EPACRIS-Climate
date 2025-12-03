/*Function to calculate the infrared Planck mean opacity at each level */
/* input cross section in m^2 */
/* output mean cross section in cm^2 */
/* also output the cross section in cm^2 at each level in each wavelength */

#include <math.h>
#include "constant.h"

void planckmeanCIA();

void planckmeanCIA()
{		
	int i, j, k, s;
	double Fplanck[zbin+1], Cplanck[zbin+1], wave[NLAMBDA], h1, h2;
    double Fstar[zbin+1], Cstar[zbin+1];
	
	printf("Starting planckmeanCIA calculation...\n");
	
	/* get wave1 from wave */
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
		// Safety check to avoid division by zero
		if (Fplanck[i] < 1e-60) {
			printf("Warning: Very small Planck flux at layer %d, setting to 1e-60\n", i);
			Fplanck[i] = 1e-60;
		}
		if (Fstar[i] < 1e-60) {
			printf("Warning: Very small stellar flux at layer %d, setting to 1e-60\n", i);
			Fstar[i] = 1e-60;
		}
	}
	
	/* Calculate the cross section */
    
    /* H2-H2 */
	printf("Calculating H2-H2 CIA means...\n");
	for (i=1; i<=zbin; i++) {
		Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
		for (j=0; j<NLAMBDA-1; j++) {
			// Add checks to avoid NaN and Inf
			if (!isfinite(H2H2CIA[i][j]) || !isfinite(H2H2CIA[i][j+1])) {
				continue;
			}
			
			h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
			h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
			
			// Check for non-finite values in calculations
			if (!isfinite(h1) || !isfinite(h2)) {
				continue;
			}
			
			if (fabs(H2H2CIA[i][j+1]-H2H2CIA[i][j])<1.0E-60) {
				Cplanck[i] += (wave[j+1]-wave[j])/2.0*(H2H2CIA[i][j]*h1+H2H2CIA[i][j+1]*h2);
			} else {
				// Avoid log(0) and division by zero
				double denom = log(H2H2CIA[i][j+1]*h2 + 1E-80) - log(H2H2CIA[i][j]*h1 + 1E-80);
				if (fabs(denom) < 1e-80) {
					Cplanck[i] += (wave[j+1]-wave[j])/2.0*(H2H2CIA[i][j]*h1+H2H2CIA[i][j+1]*h2);
				} else {
					Cplanck[i] += (wave[j+1]-wave[j])*(H2H2CIA[i][j+1]*h2-H2H2CIA[i][j]*h1)/(denom);
				}
			}
		}
		
		// Safety check for mean calculation
		if (isfinite(Cplanck[i]) && isfinite(Fplanck[i]) && fabs(Fplanck[i]) > 1e-80) {
			MeanH2H2CIA[i] = Cplanck[i]/Fplanck[i];
		} else {
			printf("Warning: Invalid calculation for MeanH2H2CIA[%d], setting to 0.0\n", i);
			MeanH2H2CIA[i] = 0.0;
		}

		if (isfinite(Cstar[i]) && isfinite(Fstar[i]) && fabs(Fstar[i]) > 1e-80) {
			SMeanH2H2CIA[i] = Cstar[i]/Fstar[i];
		} else {
			printf("Warning: Invalid calculation for SMeanH2H2CIA[%d], setting to 0.0\n", i);
			SMeanH2H2CIA[i] = 0.0;
		}
	}
	
    /* H2-He */
    for (i=1; i<=zbin; i++) {
        Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
        for (j=0; j<NLAMBDA-1; j++) {
            h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
            h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
            if (fabs(H2HeCIA[i][j+1]-H2HeCIA[i][j])<1.0E-50) {
                Cplanck[i] += (wave[j+1]-wave[j])/2.0*(H2HeCIA[i][j]*h1+H2HeCIA[i][j+1]*h2);
            } else {
                Cplanck[i] += (wave[j+1]-wave[j])*(H2HeCIA[i][j+1]*h2-H2HeCIA[i][j]*h1)/(log(H2HeCIA[i][j+1]*h2 + 1E-80)-log(H2HeCIA[i][j]*h1 + 1E-80));
            }
        }
        for (j=0; j<NLAMBDA-1; j++) {
            h1=solar[j];
            h2=solar[j+1];
            if (fabs(H2HeCIA[i][j+1]-H2HeCIA[i][j])<1.0E-60) {
                Cstar[i] += (wave[j+1]-wave[j])/2.0*(H2HeCIA[i][j]*h1+H2HeCIA[i][j+1]*h2);
            } else {
                Cstar[i] += (wave[j+1]-wave[j])*(H2HeCIA[i][j+1]*h2-H2HeCIA[i][j]*h1)/(log(H2HeCIA[i][j+1]*h2 + 1E-80)-log(H2HeCIA[i][j]*h1 + 1E-80));
            }
        }
        MeanH2HeCIA[i] = Cplanck[i]/Fplanck[i];
        SMeanH2HeCIA[i] = Cstar[i]/Fstar[i];
        /* printf("%s\t%e\t%e\n", "H2HeCIA", MeanH2HeCIA[i], SMeanH2HeCIA[i]); */
    }
    
    /* H2-H */
    for (i=1; i<=zbin; i++) {
        Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
        for (j=0; j<NLAMBDA-1; j++) {
            h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
            h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
            if (fabs(H2HCIA[i][j+1]-H2HCIA[i][j])<1.0E-60) {
                Cplanck[i] += (wave[j+1]-wave[j])/2.0*(H2HCIA[i][j]*h1+H2HCIA[i][j+1]*h2);
            } else {
                Cplanck[i] += (wave[j+1]-wave[j])*(H2HCIA[i][j+1]*h2-H2HCIA[i][j]*h1)/(log(H2HCIA[i][j+1]*h2 + 1E-80)-log(H2HCIA[i][j]*h1 + 1E-80));
            }
        }
        for (j=0; j<NLAMBDA-1; j++) {
            h1=solar[j];
            h2=solar[j+1];
            if (fabs(H2HCIA[i][j+1]-H2HCIA[i][j])<1.0E-60) {
                Cstar[i] += (wave[j+1]-wave[j])/2.0*(H2HCIA[i][j]*h1+H2HCIA[i][j+1]*h2);
            } else {
                Cstar[i] += (wave[j+1]-wave[j])*(H2HCIA[i][j+1]*h2-H2HCIA[i][j]*h1)/(log(H2HCIA[i][j+1]*h2 + 1E-80)-log(H2HCIA[i][j]*h1 + 1E-80));
            }
        }
        MeanH2HCIA[i] = Cplanck[i]/Fplanck[i];
        SMeanH2HCIA[i] = Cstar[i]/Fstar[i];
        /* printf("%s\t%e\t%e\n", "H2HCIA", MeanH2HCIA[i], SMeanH2HCIA[i]); */
    }
    
    /* N2-H2 */
    for (i=1; i<=zbin; i++) {
        Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
        for (j=0; j<NLAMBDA-1; j++) {
            h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
            h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
            if (fabs(N2H2CIA[i][j+1]-N2H2CIA[i][j])<1.0E-60) {
                Cplanck[i] += (wave[j+1]-wave[j])/2.0*(N2H2CIA[i][j]*h1+N2H2CIA[i][j+1]*h2);
            } else {
                Cplanck[i] += (wave[j+1]-wave[j])*(N2H2CIA[i][j+1]*h2-N2H2CIA[i][j]*h1)/(log(N2H2CIA[i][j+1]*h2 + 1E-80)-log(N2H2CIA[i][j]*h1 + 1E-80));
            }
        }
        for (j=0; j<NLAMBDA-1; j++) {
            h1=solar[j];
            h2=solar[j+1];
            if (fabs(N2H2CIA[i][j+1]-N2H2CIA[i][j])<1.0E-60) {
                Cstar[i] += (wave[j+1]-wave[j])/2.0*(N2H2CIA[i][j]*h1+N2H2CIA[i][j+1]*h2);
            } else {
                Cstar[i] += (wave[j+1]-wave[j])*(N2H2CIA[i][j+1]*h2-N2H2CIA[i][j]*h1)/(log(N2H2CIA[i][j+1]*h2 + 1E-80)-log(N2H2CIA[i][j]*h1 + 1E-80));
            }
        }
        MeanN2H2CIA[i] = Cplanck[i]/Fplanck[i];
        SMeanN2H2CIA[i] = Cstar[i]/Fstar[i];
        /* printf("%s\t%e\t%e\n", "N2H2CIA", MeanN2H2CIA[i], SMeanN2H2CIA[i]); */
    }
    
    /* N2-N2 */
    for (i=1; i<=zbin; i++) {
        Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
        for (j=0; j<NLAMBDA-1; j++) {
            h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
            h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
            if (fabs(N2N2CIA[i][j+1]-N2N2CIA[i][j])<1.0E-60) {
                Cplanck[i] += (wave[j+1]-wave[j])/2.0*(N2N2CIA[i][j]*h1+N2N2CIA[i][j+1]*h2);
            } else {
                Cplanck[i] += (wave[j+1]-wave[j])*(N2N2CIA[i][j+1]*h2-N2N2CIA[i][j]*h1)/(log(N2N2CIA[i][j+1]*h2 + 1E-80)-log(N2N2CIA[i][j]*h1 + 1E-80));
            }
        }
        for (j=0; j<NLAMBDA-1; j++) {
            h1=solar[j];
            h2=solar[j+1];
            if (fabs(N2N2CIA[i][j+1]-N2N2CIA[i][j])<1.0E-60) {
                Cstar[i] += (wave[j+1]-wave[j])/2.0*(N2N2CIA[i][j]*h1+N2N2CIA[i][j+1]*h2);
            } else {
                Cstar[i] += (wave[j+1]-wave[j])*(N2N2CIA[i][j+1]*h2-N2N2CIA[i][j]*h1)/(log(N2N2CIA[i][j+1]*h2 + 1E-80)-log(N2N2CIA[i][j]*h1 + 1E-80));
            }
        }
        MeanN2N2CIA[i] = Cplanck[i]/Fplanck[i];
        SMeanN2N2CIA[i] = Cstar[i]/Fstar[i];
        /* printf("%s\t%e\t%e\n", "N2N2CIA", MeanN2N2CIA[i], SMeanN2N2CIA[i]); */
    }
    
    /* CO2-CO2 */
    for (i=1; i<=zbin; i++) {
        Cplanck[i] = 0.0;
        Cstar[i] = 0.0;
        for (j=0; j<NLAMBDA-1; j++) {
            h1=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j],5)/(exp(HPLANCK*CLIGHT/wave[j]/KBOLTZMANN/tl[i])-1);
            h2=2*HPLANCK*pow(CLIGHT,2)/pow(wave[j+1],5)/(exp(HPLANCK*CLIGHT/wave[j+1]/KBOLTZMANN/tl[i])-1);
            if (fabs(CO2CO2CIA[i][j+1]-CO2CO2CIA[i][j])<1.0E-60) {
                Cplanck[i] += (wave[j+1]-wave[j])/2.0*(CO2CO2CIA[i][j]*h1+CO2CO2CIA[i][j+1]*h2);
            } else {
                Cplanck[i] += (wave[j+1]-wave[j])*(CO2CO2CIA[i][j+1]*h2-CO2CO2CIA[i][j]*h1)/(log(CO2CO2CIA[i][j+1]*h2 + 1E-80)-log(CO2CO2CIA[i][j]*h1 + 1E-80));
            }
        }
        for (j=0; j<NLAMBDA-1; j++) {
            h1=solar[j];
            h2=solar[j+1];
            if (fabs(CO2CO2CIA[i][j+1]-CO2CO2CIA[i][j])<1.0E-60) {
                Cstar[i] += (wave[j+1]-wave[j])/2.0*(CO2CO2CIA[i][j]*h1+CO2CO2CIA[i][j+1]*h2);
            } else {
                Cstar[i] += (wave[j+1]-wave[j])*(CO2CO2CIA[i][j+1]*h2-CO2CO2CIA[i][j]*h1)/(log(CO2CO2CIA[i][j+1]*h2 + 1E-80)-log(CO2CO2CIA[i][j]*h1 + 1E-80));
            }
        }
        MeanCO2CO2CIA[i] = Cplanck[i]/Fplanck[i];
        SMeanCO2CO2CIA[i] = Cstar[i]/Fstar[i];
        /* printf("%s\t%e\t%e\n", "CO2CO2CIA", MeanCO2CO2CIA[i], SMeanCO2CO2CIA[i]); */
    }
	
}
