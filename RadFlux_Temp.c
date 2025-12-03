/*Function to calculate the reflected solar light */

#include <math.h>
#include "constant.h"

void RadFlux(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint);

void RadFlux(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint)
{

	int i=0, j, k, jj, l, j1, kk;
    
    double GA;
    GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET;

	double mole2dust, planck, xlll;
	mole2dust = PI*pow(AERSIZE,3)*AERDEN/6.0/AMU*2.0558;
	
	double lam1[NLAMBDA];
	for (i=0; i<NLAMBDA; i++) {
		lam1[i] = wavelength[i]*1.0E-9; /* convert to m */
	}
    
    /* Temperature variation */
    int nrl;
    nrl = zbin-ncl;
    double **Tvar;
    double tempvar=RJACOB_TEMPVAR;
    double tempvarc;
    Tvar = dmatrix(0,zbin,0,nrl+1);
    for (j=0; j<=zbin; j++) {
        Tvar[j][0] = T[zbin-j];
    }
    for (j=0; j<=zbin; j++) {
        Tvar[j][1] = T[zbin-j];
    }
    Tvar[0][1] = Tvar[0][1]+tempvar;
    int isconv1[zbin+1],radcount;
    double lapse1[zbin+1];
    for (k=1; k<=zbin; k++) {
        isconv1[k] = isconv[zbin+1-k]; /* reverse, from top to bottom */
        lapse1[k] = lapse[zbin+1-k];
    }
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            for (j=0; j<=zbin; j++) {
                Tvar[j][radcount+1] = T[zbin-j];
            }
            Tvar[k][radcount+1] = Tvar[k][radcount+1]+tempvar;
            kk=1;
            tempvarc=tempvar;
            /* printf("%s %e %f\n", "perturb p=", P[zbin-k], tempvarc); */
            while (k+kk<=zbin && isconv1[k+kk] == 1) { /* also perturb all subsequent convective layer */
                tempvarc  = tempvarc*pow(P[zbin-k-kk]/P[zbin-k-kk+1],lapse1[k+kk]);
                /* printf("%s %e %f\n", "convective perturb p=", P[zbin-k-kk], tempvarc); */
                Tvar[k+kk][radcount+1] = Tvar[k+kk][radcount+1]+tempvarc;
                kk = kk+1;
            }
        }
    }
    /*for (j=0; j<=zbin; j++) {
        for (k=0; k<=nrl+1; k++) {
            printf("%f\t",Tvar[j][k]);
        }
        printf("\n");
    }*/
	
	double w[zbin+1], wa[zbin+1], ws[zbin+1], g[zbin+1], tau[zbin+1], tc[zbin+1], gc[zbin+1], wc[zbin+1];
	double *gamma1, *gamma2, *gamma3, *gamma4, *lambda, *gamma, *e1, *e2, *e3, *e4, *CP0, *CP1, *CM0, *CM1, *CPh, *CMh;
	double B0, B1, miu1=0.5;
	int    ntau, nvector;
	double *A, *AS, *B, *D, *DS, *E, *ES, *Y, *Y1, *Y2, dflux[zbin+1];
	double *tc1, *gc1, *wc1, *tau1, *w1, *g1;
	double temp[zbin+1], *temp1;
	double surfaceref;
    
    double solarfraction;
    solarfraction = FADV/cos(THETAREF);
	
	/* process temperature for UV cross sections */
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
	
    /* Net flux */
    double **NetFlux;
    NetFlux = dmatrix(0,nrl,0,nrl+1);
    
	for (j=0; j<=nrl; j++) {
        for (k=0; k<=nrl+1; k++) {
            NetFlux[j][k] = 0.0;
        }
	}
    
    double radiationI0, radiationI1, radiationO;
    radiationI0=0;
    radiationI1=0;
    radiationO=0;
    
    FILE *fp, *fopa;
    fp=fopen("AuxillaryOut/checkemission.dat","w");
    fopa=fopen("AuxillaryOut/checkopacity.dat","w");


//========================================================    
//====== Wavelength loop =================================	
	for (i=0; i<(NLAMBDA-1); i++) {
//========================================================    
//========================================================    


		/* printf("%s %f %f %f %f\n", "crossa", crossa[1][i], cross[2][i], sinab[1][i], sinab[2][i]); */
		/* printf("%s\t%f\n","wavelength",wavelength[i]); */
		
		/* Single Scattering Albedo */	
		for (j=1; j <= zbin; j++) {
            wa[j] = 0;
            ws[j] = 0;
            w[j]  = 0;
            j1    = zbin+1-j;
            /*for (k=1; k<=nump; k++) {
                crossl = cross[k][i] + crosst[k][i]* ( temperature[j1] - 295.0 ) ;
                wa[j] += xx[j1][stdn[k]]*crossl/qysum[k];
            }*/
            /*wa[j] += opacCO2[j1][i]*xx[j1][52];
            wa[j] += opacO2[j1][i]*xx[j1][54];*/
            wa[j] += opacH2O[j1][i]*xx[j1][7];
            wa[j] += opacCH4[j1][i]*xx[j1][21];
            wa[j] += opacNH3[j1][i]*xx[j1][9];
            /* wa[j] += opacOH[j1][i]*xx[j1][4];
            wa[j] += opacH2CO[j1][i]*xx[j1][22];
            wa[j] += opacH2O2[j1][i]*xx[j1][6];
            wa[j] += opacHO2[j1][i]*xx[j1][5];
            wa[j] += opacCO[j1][i]*xx[j1][20];
            wa[j] += opacO3[j1][i]*xx[j1][2];
            wa[j] += opacC2H2[j1][i]*xx[j1][27];
            wa[j] += opacC2H4[j1][i]*xx[j1][29];
            wa[j] += opacC2H6[j1][i]*xx[j1][31];
            wa[j] += opacCH2O2[j1][i]*xx[j1][23];
            wa[j] += opacHCN[j1][i]*xx[j1][37];
            wa[j] += opacN2O[j1][i]*xx[j1][11];
            wa[j] += opacNO[j1][i]*xx[j1][12];
            wa[j] += opacNO2[j1][i]*xx[j1][13];
            wa[j] += opacHNO3[j1][i]*xx[j1][18];
            wa[j] += opacH2S[j1][i]*xx[j1][45];
            wa[j] += opacSO2[j1][i]*xx[j1][43];
            wa[j] += opacOCS[j1][i]*xx[j1][49];*/
            
            wa[j] += H2H2CIA[j1][i]*xx[j1][53]*xx[j1][53]; /* H2H2CIA */
            wa[j] += H2HeCIA[j1][i]*xx[j1][53]*heliumnumber[j1]; /* H2HeCIA */
            /* wa[j] += H2HCIA[j1][i]*xx[j1][53]*xx[j1][3]; /* H2HCIA */
            /* wa[j] += N2N2CIA[j1][i]*xx[j1][55]*xx[j1][55]; /* N2N2CIA */
            /* wa[j] += N2H2CIA[j1][i]*xx[j1][55]*xx[j1][53]; /* N2H2CIA */
            /* wa[j] += CO2CO2CIA[j1][i]*xx[j1][52]*xx[j1][52]; /* CO2CO2CIA */
            
            ws[j] += crossr[i]*MM[j1];
            
            /* wa[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*(1.0-sinab[1][i]);
            wa[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*(1.0-sinab[2][i]);
            ws[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*sinab[1][i];
            ws[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*sinab[2][i]; */
            
            /* wa[j] += cH2O[j1][i]*(1.0-aH2O[j1][i])/MM[j1];
            wa[j] += cNH3[j1][i]*(1.0-aNH3[j1][i])/MM[j1];
            ws[j] += cH2O[j1][i]*aH2O[j1][i]/MM[j1];
            ws[j] += cNH3[j1][i]*aNH3[j1][i]/MM[j1]; */
            
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
			/* g[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*sinab[1][i]*asym[1][i];
			g[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*sinab[2][i]*asym[2][i]; */
            /*g[j] += cH2O[j1][i]*aH2O[j1][i]*gH2O[j1][i]/MM[j1];
			g[j] += cNH3[j1][i]*aNH3[j1][i]*gNH3[j1][i]/MM[j1];*/
			if (ws[j] > 0.0) {
				g[j] = g[j]/ws[j];
			} else {
				g[j] = 0.0;
			}
		}
		
		/* Optical Depth of Each Layer */
		for (j=1; j <=zbin; j++) {
            j1   = zbin+1-j;
            tau[j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-P[j1])*1.0E-4;
			/* printf("%f %e %f %f\n", wavelength[i], tau[j], w[j], g[j]); */
		}
		
        /* Delta Adjustment */
        if (DELADJUST == 1) {
            for (j=1; j <= zbin; j++) {
                tau[j] *= (1-w[j]*g[j]*g[j]);
                w[j]   =  (1-g[j]*g[j])*w[j]/(1-w[j]*g[j]*g[j]);
                g[j]   =  g[j]/(1+g[j]);
            }
        }
        
		/* Cumulative Optical Depth */
		for (j=0; j <= zbin; j++) {
			tc[j] = 0.0;
		}
		for (j=1; j <= zbin; j++) {
			for (jj=1; jj <= j; jj++) {
				tc[j] += tau[jj];
			}
		}
		
		/* Determine if multi-layer computation is needed */
		/* if (tc[zbin] < TAUTHRESHOLD) {
			
			B0      = 2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/temp[zbin])-1.0)*1.0E-9;
			flux[i] = ((1-PSURFAB)*PI*B0 + PSURFAB*cos(THETAREF)*exp(-tc[zbin]/cos(THETAREF))*solar[i]*exp(-tc[zbin]/miu1);
			ratio[i] = flux[i]/solar[i];
			
		} else { */
			
			/* Determine the optical depth scale for computation */
		
		if (wavelength[i]>3000.0) {
			surfaceref = 1.0-PSURFEM;
		} else {
			surfaceref = PSURFAB;
		}
        
        j=1;
        while (j<zbin) { /* calculate all levels */
            j++;
        }
        ntau=j;
        nvector = 2*ntau;
        
        /*if (i==4600) {
         printf("%f %s %d %f %f\n", wavelength[i], "ntau", ntau, tc[zbin], tc[0]);
         }*/
        
        tc1  = dvector(0,ntau);
        temp1= dvector(0,ntau);
        tau1 = dvector(1,ntau);
        w1   = dvector(1,ntau);
        g1   = dvector(1,ntau);
        
        for (j=1; j<=ntau; j++) {
            tau1[j]=tau[j];
            /*markus2021:testing for small tau*/
            //if(tau1[j]<1.e-6) tau1[j]=1.e-6;
            w1[j]=w[j];
            g1[j]=g[j];
        }
        for (j=0; j<=ntau; j++) {
            tc1[j]=tc[j];
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
            CPh    = dvector(1,ntau);
            CMh    = dvector(1,ntau);
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
			
		/* Two stream parameters */
		if (wavelength[i] > 3000.0) {
                    //Hemispheric Mean:
			for (j=1; j <= ntau; j++) {
				gamma1[j] = 2.0-w1[j]*(1+g1[j]);
				gamma2[j] = w1[j]*(1.0-g1[j]);
				gamma3[j] = 0.5-g1[j]*cos(THETAREF);
				gamma4[j] = 1.0-gamma3[j];
			}
		} else {
			/* Eddington */
			for (j=1; j <= ntau; j++) {
				gamma1[j] = (7.0-(4.0+3.0*g1[j])*w1[j])/4.0;
				gamma2[j] = -(1.0-(4.0-3.0*g1[j])*w1[j])/4.0;
				gamma3[j] = (2.0-3.0*g1[j]*cos(THETAREF))/4.0;
				gamma4[j] = 1.0-gamma3[j];
			}
		}

			/* Two-stream computation scheme */
			/* some expressions */
			for (j=1; j <= ntau; j++) {
				lambda[j] = sqrt(gamma1[j]*gamma1[j]-gamma2[j]*gamma2[j]);
				gamma[j]  = gamma2[j]/(gamma1[j]+lambda[j]);
				e1[j] = 1.0+gamma[j]*exp(-lambda[j]*tau1[j]);
				e2[j] = 1.0-gamma[j]*exp(-lambda[j]*tau1[j]);
				e3[j] = gamma[j]+exp(-lambda[j]*tau1[j]);
				e4[j] = gamma[j]-exp(-lambda[j]*tau1[j]);
				/* printf("%d %f %f %f %f\n", j, e1[j], e2[j], e3[j], e4[j]); */
			}
        
        /*fprintf(fp, "%e\t%e\t", wavelength[i], solar[i]); */
        
        /* Loop the temperature variation */
        for (kk=0; kk<=nrl+1; kk++) {
            
            /* Get temperature */
            for (j=0; j<=ntau; j++) {
                temp1[j] = Tvar[j][kk];
            }
            
            /* source functions */
			for (j=1; j <= ntau; j++) {
				CP0[j] = 0.0;
				CP1[j] = 0.0;
				CM0[j] = 0.0;
				CM1[j] = 0.0;
                CPh[j] = 0.0;
				CMh[j] = 0.0;
                
                if (wavelength[i] > 3000.0) {
				B0     = 2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/temp1[j-1])-1.0)*1.0E-9;
		  		B1     = ( 2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/temp1[j])-1.0)*1.0E-9 - B0 )/tau1[j];
				CP0[j] += 2*PI*miu1*(B0+B1/(gamma1[j]+gamma2[j]));
				CP1[j] += 2*PI*miu1*(B0+B1*(tau1[j]+1.0/(gamma1[j]+gamma2[j])));
                CPh[j] += 2*PI*miu1*(B0+B1*(tau1[j]*0.5+1.0/(gamma1[j]+gamma2[j])));
				CM0[j] += 2*PI*miu1*(B0-B1/(gamma1[j]+gamma2[j]));
				CM1[j] += 2*PI*miu1*(B0+B1*(tau1[j]-1.0/(gamma1[j]+gamma2[j])));
                CMh[j] += 2*PI*miu1*(B0+B1*(tau1[j]*0.5-1.0/(gamma1[j]+gamma2[j])));
                } else {
                
				CP0[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
				CP1[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                CPh[j] += solarfraction*solar[i]*w1[j]*exp(-(tc1[j-1]+0.5*tau1[j])/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
				CM0[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
				CM1[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                CMh[j] += solarfraction*solar[i]*w1[j]*exp(-(tc1[j-1]+0.5*tau1[j])/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                    
                }
                
				/* printf("%s %f %d %e\n", "CP0", wavelength[i], j, CP0[j]);
                 printf("%s %f %d %e\n", "CP1", wavelength[i], j, CP1[j]);
                 printf("%s %f %d %e\n", "CM0", wavelength[i], j, CM0[j]);
                 printf("%s %f %d %e\n", "CM1", wavelength[i], j, CM1[j]); */
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
                                //markus2021:comparison with Toon89: -/+ sign changed:
				E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) - e4[j+1]*(CM0[j+1]-CM1[j]);
				//E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) + e4[j+1]*(CM0[j+1]-CM1[j]);
			}
			l = nvector;
            surfaceref = 1.0;
			A[l] = e1[ntau] - surfaceref*e3[ntau];
			B[l] = e2[ntau] - surfaceref*e4[ntau];
			D[l] = 0;
			E[l] = -CP1[ntau] + surfaceref*CM1[ntau];
			E[l] += (1-surfaceref)*PI*2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/temp1[ntau])-1.0)*1.0E-9;
            E[l] += PI*2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/Tint)-1.0)*1.0E-9;
			/* printf("%s %f %e\n","surf temp",temp1[ntau-1],tc1[ntau]); */
			E[l] += surfaceref*cos(THETAREF)*solar[i]*solarfraction*exp(-tc1[ntau]/cos(THETAREF));
			/* for (l=1; l<=nvector; l++) {
             printf("%s %f %d %f %f %f %f\n", "ABDE", wavelength[i], l, A[l], B[l], D[l], E[l]);
             } */
            
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
            
            /* obtain net flux at zero and midpoint */
            /*for (j=1; j<=ntau; j++) {
                aflux[j] = Y1[j]*(e1[j]-e3[j])+Y2[j]*(e2[j]-e4[j])+CP1[j]-CM1[j];
            }*/
            dflux[0] = Y1[1]*e3[1] - Y2[1]*e4[1] + CP0[1];
            for (j=1; j<=ntau; j++) {
                dflux[j] =  2.0*Y2[j]*(1.0-gamma[j])*exp(-lambda[j]*tau1[j]*0.5) + CPh[j] - CMh[j];
            }
            
            /* fprintf(fp, "%e\t", dflux[0]);
            
            if (kk==nrl+1) {
                fprintf(fp, "\n");
            } */
            
            /* add direct sunlight */
            dflux[0] -= cos(THETAREF)*solar[i]*solarfraction;
            for (j=1; j<=ntau; j++) {
                dflux[j] -= cos(THETAREF)*solar[i]*solarfraction*exp(-(tc1[j-1]+tau1[j]*0.5)/cos(THETAREF));
            }

            
            NetFlux[0][kk] += dflux[0]*(wavelength[i+1]-wavelength[i]);
            radcount = 0;
            for (k=1; k<=zbin; k++) {
                if (isconv1[k] == 0) {
                    radcount += 1;
                    NetFlux[radcount][kk] += dflux[k]*(wavelength[i+1]-wavelength[i]);
                }
            }
            
            if (kk==0) {
                radiationI1 += cos(THETAREF)*solar[i]*solarfraction*exp(-tc[zbin]/cos(THETAREF))*(wavelength[i+1]-wavelength[i]); /* BOA incoming flux */
                radiationI0 += cos(THETAREF)*solar[i]*solarfraction*(wavelength[i+1]-wavelength[i]); /* TOA incoming flux */
                radiationO += dflux[0]*(wavelength[i+1]-wavelength[i]); /* TOA net outgoing flux */
            }
            
        }
        
		
			free_dvector(tc1,0,ntau);
			free_dvector(temp1,0,ntau);
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
            free_dvector(CPh,1,ntau);
            free_dvector(CMh,1,ntau);
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
        
    }
    
    fclose(fp);
    fclose(fopa);

    /* Obtain the Flux and the Flux Jacobian, from bottom to up */
    for (j=0; j<=nrl; j++) {
        Rflux[j] = NetFlux[j][0] - SIGMA*pow(Tint,4.0);
    }
    
    /*printf("%e\t%f\n",P[zbin],Rflux[0]);
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k] == 0) {
            radcount += 1;
            printf("%e\t%f\n",pl[zbin-k+1],Rflux[radcount]);
        }
    }*/
    
    double **jmax;
    jmax=dmatrix(1,nrl+1,1,nrl+1);
    for (j=0; j<=nrl; j++) {
        for (kk=0; kk<=nrl; kk++) {
            jmax[j+1][kk+1] = (NetFlux[j][kk+1] - NetFlux[j][0])/tempvar;
        }
    }
    
    /* for (j=0; j<=nrl; j++) {
        for (k=0; k<=5; k++) {
            printf("%f\t",NetFlux[j][k]);
        }
        printf("\n");
    } */
    
    /* for (j=1; j<=nrl+1; j++) {
        for (k=1; k<=6; k++) {
            printf("%e\t",jmax[j][k]);
        }
        printf("\n");
    } */
    
    free_dmatrix(NetFlux,0,nrl,0,nrl+1);
    
    /* Solve Equation */
    int *indx;
    double ddd;
    double *deltaT;
    double relaxationfactor;
    indx=ivector(1,nrl+1);
    deltaT=dvector(1,nrl+1);
    for (j=0; j<=nrl; j++) {
        deltaT[j+1] = -Rflux[j];
    }
    
    ludcmp(jmax,nrl+1,indx,&ddd);
    lubksb(jmax,nrl+1,indx,deltaT);
    /* for (j=0; j<=nrl; j++) {
        printf("%f\n",deltaT[j+1]);
    } */
    
    /* calculate the new temperature profile */
    radcount = 0;
    relaxationfactor = R_RELAX;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer*/
            radcount = radcount+1;
            relaxationfactor=fmin(relaxationfactor,0.1*Tvar[k][0]/fabs(deltaT[radcount+1])); /* do not allow temperature to change more than 10% */
            printf("%s\t%f\n", "Relaxation factor is", relaxationfactor);
        }
    }
    for (j=1; j<=nrl+1; j++) {
        deltaT[j] *= relaxationfactor; /* relaxation */
    }
    Tvar[0][0] = Tvar[0][0] + deltaT[1];
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            tempvarc = deltaT[radcount+1];
            Tvar[k][0] = Tvar[k][0] + tempvarc;
            /* printf("%s %e %f\n", "update p=", P[zbin-k], tempvarc); */
            kk=1;
            while (k+kk<=zbin && isconv1[k+kk] == 1) { /* also perturb all subsequent convective layer */
                tempvarc  = tempvarc*pow(P[zbin-k-kk]/P[zbin-k-kk+1],lapse1[k+kk]);
                Tvar[k+kk][0] = Tvar[k+kk][0] + tempvarc;
                /* printf("%s %e %f\n", "convective update p=", P[zbin-k-kk], tempvarc); */
                kk = kk+1;
            }
        }
    }
    
    for (j=0; j<=zbin; j++) {
        tempbnew[j] = fmax(Tvar[zbin-j][0],1.0);
    }
    
    free_dmatrix(Tvar,0,zbin,0,nrl+1);
    free_dmatrix(jmax,1,nrl+1,1,nrl+1);
    free_dvector(deltaT,1,nrl+1);
    free_ivector(indx,1,nrl+1);
    
    /*printf("%s\t%f\n", "Top-of-Atmosphere incoming radiation flux is", radiationI0);
    printf("%s\t%f\n", "Bottom-of-Atmospehre incoming radiation flux is", radiationI1);*/
    printf("%s\t%f\n", "Tof-of-Atmosphere outgoing net radiation flux is", radiationO);
	
}
