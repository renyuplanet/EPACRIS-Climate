/*Function to provide solar radiation at a given height*/
#include <math.h>
#include "constant.h"

void RadTransfer(double **radiation, double **OptD, int stdn[], double qysum[], double **cross, double **crosst, int NWave);

void RadTransfer(double **radiation, double **OptD, int stdn[], double qysum[], double **cross, double **crosst, int NWave)
{
/*Line by Line calculation of radiation field*/
	int i=0, j, k, jj, l, j1;	

	double mole2dust;
	mole2dust = PI*pow(AERSIZE,3)*AERDEN/6.0/AMU*2.0558;
	
	double w[zbin+1], wa[zbin+1], ws[zbin+1], g[zbin+1], tau[zbin+1], tc[zbin+1], gc[zbin+1], wc[zbin+1];
	double *gamma1, *gamma2, *gamma3, *gamma4, *lambda, *gamma, *e1, *e2, *e3, *e4, *CP0, *CP1, *CM0, *CM1;
	double B0, B1, miu1=0.5, xlll, surfaceref;
	int    ntau, nvector;
	double *A, *AS, *B, *D, *DS, *E, *ES, *Y, *Y1, *Y2, *aflux, dflux[zbin+1];
	double *tc1, *gc1, *wc1, *tau1, *w1, *g1;
	
	/* process temperature for cross sections */
	double temperature[zbin+1], crossl;
	for (j=1; j<=zbin; j++) {
		temperature[j] = tl[j];
		if (temperature[j] > TDEPMAX) {
			temperature[j] = TDEPMAX; 
		}
		if (temperature[j] < TDEPMIN) {
			temperature[j] = TDEPMIN; 
		}
	}
    
    double heliumnumber[zbin+1];
    for (j=1; j<=zbin; j++) {
        heliumnumber[j] = MM[j];
        for (i=1; i<=NSP; i++) {
            heliumnumber[j] -= xx[j][i];
        }
        if (heliumnumber[j]<0.0) {
            heliumnumber[j] = 0.0;
        }
    }
    
	
	for (i=0; i<NWave; i++) {
		
		/* Cumulative Optical Depth */
		for (j=0; j<=zbin; j++) {
			OptD[i][j]=0.0;
		}
		for (j=zbin-1; j>=0; j--) {
			OptD[i][j] = OptD[i][j+1];
			/* molecular absorption associated with photolysis */
			for (k=1; k<=nump; k++) { 
				crossl = cross[k][i] + crosst[k][i]* ( temperature[j+1] - 295.0 ) ;
				OptD[i][j] += xx[j+1][stdn[k]]*crossl*thickl/qysum[k];
			}
			/* rayleigh scattering */
			OptD[i][j] += crossr[i]*MM[j+1]*thickl;
			/* Mie scattering from H2SO4AER */
			OptD[i][j] += crossa[1][i]*xx[j+1][78]*98.0/mole2dust*thickl;
			/* Mie scattering from S8AER */
			OptD[i][j] += crossa[2][i]*xx[j+1][111]*256.0/mole2dust*thickl;
            
			/* molecular absorption not associated with photolysis */
			OptD[i][j] += opacCO2[j+1][i]*xx[j+1][52]*thickl;
			OptD[i][j] += opacO2[j+1][i]*xx[j+1][54]*thickl;
			OptD[i][j] += opacH2O[j+1][i]*xx[j+1][7]*thickl;
			OptD[i][j] += opacOH[j+1][i]*xx[j+1][4]*thickl;
			OptD[i][j] += opacH2CO[j+1][i]*xx[j+1][22]*thickl;
			OptD[i][j] += opacH2O2[j+1][i]*xx[j+1][6]*thickl;
			OptD[i][j] += opacHO2[j+1][i]*xx[j+1][5]*thickl;
            OptD[i][j] += opacCO[j+1][i]*xx[j+1][20]*thickl;
			OptD[i][j] += opacO3[j+1][i]*xx[j+1][2]*thickl;
			OptD[i][j] += opacCH4[j+1][i]*xx[j+1][21]*thickl;
            OptD[i][j] += opacC2H2[j+1][i]*xx[j+1][27]*thickl;
            OptD[i][j] += opacC2H4[j+1][i]*xx[j+1][29]*thickl;
            OptD[i][j] += opacC2H6[j+1][i]*xx[j+1][31]*thickl;
            OptD[i][j] += opacCH2O2[j+1][i]*xx[j+1][23]*thickl;
            OptD[i][j] += opacHCN[j+1][i]*xx[j+1][37]*thickl;
			OptD[i][j] += opacNH3[j+1][i]*xx[j+1][9]*thickl;
            OptD[i][j] += opacN2O[j+1][i]*xx[j+1][11]*thickl;
            OptD[i][j] += opacNO[j+1][i]*xx[j+1][12]*thickl;
            OptD[i][j] += opacNO2[j+1][i]*xx[j+1][13]*thickl;
            OptD[i][j] += opacHNO3[j+1][i]*xx[j+1][18]*thickl;
            OptD[i][j] += opacH2S[j+1][i]*xx[j+1][45]*thickl;
            OptD[i][j] += opacSO2[j+1][i]*xx[j+1][43]*thickl;
            OptD[i][j] += opacOCS[j+1][i]*xx[j+1][49]*thickl;
            
            OptD[i][j] += H2H2CIA[j+1][i]*xx[j+1][53]*xx[j+1][53]*thickl; /* H2H2CIA */
            OptD[i][j] += H2HeCIA[j+1][i]*xx[j+1][53]*heliumnumber[j+1]*thickl; /* H2HeCIA */
            OptD[i][j] += H2HCIA[j+1][i]*xx[j+1][53]*xx[j+1][3]*thickl; /* H2HCIA */
            OptD[i][j] += N2N2CIA[j+1][i]*xx[j+1][55]*xx[j+1][55]*thickl; /* N2N2CIA */
            OptD[i][j] += N2H2CIA[j+1][i]*xx[j+1][55]*xx[j+1][53]*thickl; /* N2H2CIA */
            OptD[i][j] += CO2CO2CIA[j+1][i]*xx[j+1][52]*xx[j+1][52]*thickl; /* CO2CO2CIA */
			
			/* if (wavelength[i]==121.6) {
				for (k=1; k<=nump; k++) { 
					crossl = cross[k][i] + crosst[k][i]* ( temperature[j+1] - 295.0 ) ;
					printf("%d %d %e\n", j,k, xx[j+1][stdn[k]]*crossl*thickl/qysum[k]);
				}
				printf("%d %s %e\n", j, "accu", OptD[i][j]);
			} */
			
		}

		
		/* Actinic Flux For Photolysis */
		for (j=0; j<=zbin; j++) {
			radiation[i][j] = solar[i]*exp(-OptD[i][j]/cos(THETAREF)); /* direct sunlight */
		}
		
		/* Diffuse Radiation */
		if (IFDIFFUSE == 1) {
			
            if (wavelength[i]>3000.0) {
                surfaceref = 1.0-PSURFEM;
            } else {
                surfaceref = PSURFAB;
            }
			
			/* Single Scattering Albedo */	
			for (j=1; j <= zbin; j++) {
				wa[j] = 0;
				ws[j] = 0;
				w[j]  = 0;
				j1    = zbin+1-j;
				for (k=1; k<=nump; k++) {
                    crossl = cross[k][i] + crosst[k][i]* ( temperature[j1] - 295.0 ) ;
                    wa[j] += xx[j1][stdn[k]]*crossl/qysum[k];
                }
				wa[j] += opacCO2[j1][i]*xx[j1][52];
                wa[j] += opacO2[j1][i]*xx[j1][54];
				wa[j] += opacH2O[j1][i]*xx[j1][7];
				wa[j] += opacOH[j1][i]*xx[j1][4];
				wa[j] += opacH2CO[j1][i]*xx[j1][22];
				wa[j] += opacH2O2[j1][i]*xx[j1][6];
				wa[j] += opacHO2[j1][i]*xx[j1][5];
				wa[j] += opacCO[j1][i]*xx[j1][20];
				wa[j] += opacO3[j1][i]*xx[j1][2];
				wa[j] += opacCH4[j1][i]*xx[j1][21];
                wa[j] += opacC2H2[j1][i]*xx[j1][27];
                wa[j] += opacC2H4[j1][i]*xx[j1][29];
                wa[j] += opacC2H6[j1][i]*xx[j1][31];
                wa[j] += opacCH2O2[j1][i]*xx[j1][23];
				wa[j] += opacHCN[j1][i]*xx[j1][37];
                wa[j] += opacNH3[j1][i]*xx[j1][9];
                wa[j] += opacN2O[j1][i]*xx[j1][11];
                wa[j] += opacNO[j1][i]*xx[j1][12];
                wa[j] += opacNO2[j1][i]*xx[j1][13];
                wa[j] += opacHNO3[j1][i]*xx[j1][18];
                wa[j] += opacH2S[j1][i]*xx[j1][45];
                wa[j] += opacSO2[j1][i]*xx[j1][43];
                wa[j] += opacOCS[j1][i]*xx[j1][49];
                
                wa[j] += H2H2CIA[j1][i]*xx[j1][53]*xx[j1][53]; /* H2H2CIA */
                wa[j] += H2HeCIA[j1][i]*xx[j1][53]*heliumnumber[j1]; /* H2HeCIA */
                wa[j] += H2HCIA[j1][i]*xx[j1][53]*xx[j1][3]; /* H2HCIA */
                wa[j] += N2N2CIA[j1][i]*xx[j1][55]*xx[j1][55]; /* N2N2CIA */
                wa[j] += N2H2CIA[j1][i]*xx[j1][55]*xx[j1][53]; /* N2H2CIA */
                wa[j] += CO2CO2CIA[j1][i]*xx[j1][52]*xx[j1][52]; /* CO2CO2CIA */
                
				wa[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*(1.0-sinab[1][i]);
				wa[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*(1.0-sinab[2][i]);
				ws[j] += crossr[i]*MM[j1];
				ws[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*sinab[1][i];
				ws[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*sinab[2][i];
				if (ws[j] > 0.0) {
					w[j]  = ws[j]/(wa[j]+ws[j]);
				} else {
					w[j]  = 0.0;
				}
                if (w[j] > 0.9999999999999) {
                    w[j] = 0.9999999999999;
                }
                if (w[j] < 0.0000000000001) {
                    w[j] = 0.0000000000001;
                }
			}
			
			/* Asymmetry Factor */
			for (j=1; j <= zbin; j++) {
				g[j] = 0.0;
				j1   = zbin+1-j;
				g[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*sinab[1][i]*asym[1][i];
				g[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*sinab[2][i]*asym[2][i];
				if (ws[j] > 0.0) {
					g[j] = g[j]/ws[j];
				} else {
					g[j] = 0;
				}
			}
			
			/* Optical Depth of Each Layer */
			for (j=1; j <= zbin; j++) {
				tau[j] = (wa[j] + ws[j])*thickl;
			}
			
			/* Cumulative Optical Depth */
			for (j=0; j <= zbin; j++) {
				tc[j] = 0.0;
				gc[j] = 0.0;
				wc[j] = 0.0;
			}
			for (j=1; j <= zbin; j++) {
				for (jj=1; jj <= j; jj++) {
					tc[j] += tau[jj];
				}
			}
			gc[0] = g[1];
			wc[0] = w[1];
			gc[zbin] = g[zbin];
			wc[zbin] = w[zbin];
			for (j=1; j<zbin; j++) {
				gc[j] = 0.5*(g[j]+g[j+1]);
				wc[j] = 0.5*(w[j]+w[j+1]);
			}
			
			if (tc[zbin] > TAUMAX) {
				surfaceref = 0.0;
			}
            
            if (tc[zbin] > 0.1) {
            
			/* Determine the optical depth scale for computation */;
			ntau = ceil((fmin(tc[zbin],TAUMAX) - tc[0])/TAUTHRESHOLD);
			nvector = 2*ntau;
			tc1  = dvector(0,ntau);
			gc1  = dvector(0,ntau);
			wc1  = dvector(0,ntau);
			tau1 = dvector(1,ntau);
			w1   = dvector(1,ntau);
			g1   = dvector(1,ntau);
			
			for (j=1; j<=ntau; j++) {
				tau1[j] = (fmin(tc[zbin],TAUMAX) - tc[0])/ntau;
			}
			for (j=0; j<=ntau; j++) {
				tc1[j] = tc[0] + j*tau1[1];
			}
			Interpolation(tc1, ntau+1, gc1, tc, gc, zbin+1, 0);
			Interpolation(tc1, ntau+1, wc1, tc, wc, zbin+1, 0);
			for (j=1; j<=ntau; j++) {
				w1[j] = 0.5*(wc1[j-1]+wc1[j]);
				g1[j] = 0.5*(gc1[j-1]+gc1[j]);
			}
			
			/* Delta Function Adjustment */
			if (DELADJUST == 1) {
				for (j=1; j <= ntau; j++) {
					tau1[j] *= (1-w1[j]*g1[j]*g1[j]);
					w1[j]   =  (1-g1[j]*g1[j])*w1[j]/(1-w1[j]*g1[j]*g1[j]);
					g1[j]   =  g1[j]/(1+g1[j]);
				}
			}
			
			gamma1 = dvector(1,ntau);
			gamma2 = dvector(1,ntau);
			gamma3 = dvector(1,ntau);
			gamma4 = dvector(1,ntau);
			gamma  = dvector(1,ntau);
			lambda = dvector(1,ntau);
			e1     = dvector(1,ntau);
			e2     = dvector(1,ntau);
			e3     = dvector(1,ntau);
			e4     = dvector(1,ntau);
			CP0    = dvector(1,ntau);
			CP1    = dvector(1,ntau);
			CM0    = dvector(1,ntau);
			CM1    = dvector(1,ntau);
			A      = dvector(1,nvector);
			B      = dvector(1,nvector);
			D      = dvector(1,nvector);
			E      = dvector(1,nvector);
			AS     = dvector(1,nvector);
			DS     = dvector(1,nvector);
			ES     = dvector(1,nvector);
			Y      = dvector(1,nvector);
			Y1     = dvector(1,ntau);
			Y2     = dvector(1,ntau);
			aflux  = dvector(0,ntau);
			
			/* Eddington Approximation */
			for (j=1; j <= ntau; j++) {
				gamma1[j] = (7.0-(4.0+3.0*g1[j])*w1[j])/4.0;
				gamma2[j] = -(1.0-(4.0-3.0*g1[j])*w1[j])/4.0;
				gamma3[j] = (2.0-3.0*g1[j]*cos(THETAREF))/4.0;
				gamma4[j] = 1.0-gamma3[j];
			}
			/* Two-stream computation scheme */
			/* some expressions */
			for (j=1; j <= ntau; j++) {
				lambda[j] = sqrt(gamma1[j]*gamma1[j]-gamma2[j]*gamma2[j]);
				gamma[j]  = (gamma1[j]-lambda[j])/gamma2[j];
				e1[j] = 1.0+gamma[j]*exp(-lambda[j]*tau1[j]);
				e2[j] = 1.0-gamma[j]*exp(-lambda[j]*tau1[j]);
				e3[j] = gamma[j]+exp(-lambda[j]*tau1[j]);
				e4[j] = gamma[j]-exp(-lambda[j]*tau1[j]);
			}
			/* source functions */
			for (j=1; j <= ntau; j++) {
				CP0[j] = 0.0;
				CP1[j] = 0.0;
				CM0[j] = 0.0;
				CM1[j] = 0.0;
				CP0[j] += solar[i]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
				CP1[j] += solar[i]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
				CM0[j] += solar[i]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
				CM1[j] += solar[i]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
			}
			/* equations */
			for (l=1; l<=nvector; l++) {
				A[l] = 0.0;
				B[l] = 0.0;
				D[l] = 0.0;
				E[l] = 0.0;
				Y[l] = 0.0;
				AS[l] = 0.0;
				DS[l] = 0.0;
				ES[l] = 0.0;
			}
			A[1] = 0;
			B[1] = e1[1];
			D[1] = -e2[1];
			E[1] = -CM0[1];
			for (j=1; j < ntau; j++) {
				l = 2*j+1;
				A[l] = e2[j]*e3[j] - e4[j]*e1[j];
				B[l] = e1[j]*e1[j+1] - e3[j]*e3[j+1];
				D[l] = e3[j]*e4[j+1] - e1[j]*e2[j+1];
				E[l] = e3[j]*(CP0[j+1]-CP1[j]) - e1[j]*(CM0[j+1]-CM1[j]);
			}
			for (j=1; j < ntau; j++) {
				l = 2*j;
				A[l] = e2[j+1]*e1[j] - e3[j]*e4[j+1];
				B[l] = e2[j]*e2[j+1] - e4[j]*e4[j+1];
				D[l] = e1[j+1]*e4[j+1] - e2[j+1]*e3[j+1];
				E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) - e4[j+1]*(CM0[j+1]-CM1[j]);
			}
			l = nvector;
			A[l] = e1[ntau] - surfaceref*e3[ntau];
			B[l] = e2[ntau] - surfaceref*e4[ntau];
			D[l] = 0;
			E[l] = -CP1[ntau] + surfaceref*CM1[ntau];
			E[l] += surfaceref*cos(THETAREF)*solar[i]*exp(-tc1[ntau]/cos(THETAREF));
			/* solve equations */
			l=nvector;
			AS[l] = A[l]/B[l];
			DS[l] = E[l]/B[l];
			for (l=nvector-1; l>0; l--) {
				xlll = 1.0/(B[l] - D[l]*AS[l+1]);
				AS[l] = A[l]*xlll;
				DS[l] = (E[l]-D[l]*DS[l+1])*xlll;
			}
			Y[1] = DS[1];
			for (l=2; l <= nvector; l++) {
				Y[l] = DS[l] - AS[l]*Y[l-1];
			}
			for (j=1; j <= ntau; j++) {
				Y1[j] = Y[2*j-1];
				Y2[j] = Y[2*j];
			}
			/* obtain diffuse actinic flux */
			aflux[0] = 1/miu1*(Y1[1]*e3[1] - Y2[1]*e4[1] + CP0[1]);
			for (j=1; j<=ntau; j++) {
				aflux[j] = 1/miu1*(Y1[j]*(e1[j]+e3[j])+Y2[j]*(e2[j]+e4[j])+CP1[j]+CM1[j]);
			}
			for (j=0; j<=zbin; j++) {
				dflux[j] = 0.0;
			}
			Interpolation(tc,zbin+1,dflux,tc1,aflux,ntau+1,0);
			
			/* final */
			for (j=0; j<=zbin; j++) {
				if (dflux[zbin-j]>0.0) {
					radiation[i][j] += dflux[zbin-j]; /* add diffuse radiation */
				}
			}
			
			free_dvector(tc1,0,ntau);
			free_dvector(wc1,0,ntau);
			free_dvector(gc1,0,ntau);
			free_dvector(w1,1,ntau);
			free_dvector(g1,1,ntau);
			free_dvector(tau1,1,ntau);
			free_dvector(gamma1,1,ntau);
			free_dvector(gamma2,1,ntau);
			free_dvector(gamma3,1,ntau);
			free_dvector(gamma4,1,ntau);
			free_dvector(gamma,1,ntau);
			free_dvector(lambda,1,ntau);
			free_dvector(e1,1,ntau);
			free_dvector(e2,1,ntau);
			free_dvector(e3,1,ntau);
			free_dvector(e4,1,ntau);
			free_dvector(CP0,1,ntau);
			free_dvector(CP1,1,ntau);
			free_dvector(CM0,1,ntau);
			free_dvector(CM1,1,ntau);
			free_dvector(Y1,1,ntau);
			free_dvector(Y2,1,ntau);
			free_dvector(Y,1,nvector);
			free_dvector(A,1,nvector);
			free_dvector(B,1,nvector);
			free_dvector(D,1,nvector);
			free_dvector(E,1,nvector);
			free_dvector(AS,1,nvector);
			free_dvector(DS,1,nvector);
			free_dvector(ES,1,nvector);
			free_dvector(aflux,0,ntau);
                
            }
			
			
		}
		
	}
	
}
