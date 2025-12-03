/*Function to calculate the reflected solar light */
/*
 * Original Author: Renyu Hu (renyu.hu@jpl.nasa.gov)
 * Version: 2021, 2-Stream update to Heng+2019 flux equations including direct stellar beam
 * Editor v2023: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * ms tc -Update 2023: Tvar now center layer tl, rather than boundary T (fluxes at boundaries informed by T @ centers)
 * also now always all layers used in RT
 * compare against ms_radflux_temp.c for changes
 * Modifier symbols: 
    * 'ms_' ... in file names
    * 'markus<date>' ... in description lines/blocks 
    * 'ms:' ... in-line comments
    * 'OJO' ... warnings to be addressed before public release
*/


#include <math.h>
#include "constant.h"
//#include "ms_2stream.c" //ms: solving for radiative fluxes
//#include "ms_d2stream.c" //ms: solving for radiative fluxes
#include "ms_tc2stream.c" //ms: solving for radiative fluxes

void ms_RadFlux(double Rflux[], double tempcnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[]);

void ms_RadFlux(double Rflux[], double tempcnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[])
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
//ms: for flux jacobian to solve for steady-state temperatures
    int nrl, jTvar;
    nrl = zbin-ncl;
    double **Tvar;
    double tempvar=RJACOB_TEMPVAR;
    double tempvarc;
    int isconv1[zbin+1],radcount;
    double lapse1[zbin+1];
    for (k=1; k<=zbin; k++) {
        isconv1[k] = isconv[zbin+1-k]; /* reverse, from top to bottom */
        lapse1[k] = lapse[zbin+1-k]; // OJO: lapse is defined in one layer, but used now for tl of layer below
    }
    radcount = 0;

    jTvar = zbin+1; //RT always has to run for all layers 
//    if (isconv[0]==0) jTvar=jTvar+1; //if BOA is radiative, then also surface in Tvar

    if (TIME_STEPPING) jTvar = 0;
    Tvar = dmatrix(0,zbin+1,0,jTvar);
    for (j=1; j<=zbin; j++) {
        Tvar[j][0] = tl[zbin+1-j];
    }
    Tvar[0][0] = T[zbin];
    Tvar[zbin+1][0] = T[0]; //surface
    if (!TIME_STEPPING)//otherwise Tvar is effectively a vector
    {
//    for (j=1; j<=zbin; j++) {
//        Tvar[j][1] = tl[zbin+1-j];
//    }
//    Tvar[zbin+1][0] = T[0];
//    Tvar[0][1] = Tvar[0][1]+tempvar; //here no testing for convection below...
    
//    Tvar[0][1] = Tvar[0][1]*(1+0.01*tempvar); //ms22: testing

    for (k=1; k<=zbin+1; k++) {
            for (j=1; j<=zbin; j++) {
                Tvar[j][k] = tl[zbin+1-j]; //ms2021: waste of computation within k loop
            }
            Tvar[0][k] = T[zbin];
            Tvar[zbin+1][k] = T[0];
            Tvar[k][k] = Tvar[k][k]+tempvar;
//            Tvar[k][radcount] = Tvar[k][radcount]*(1+0.01*tempvar);
//ms tc            kk=1; //ms: longer than just writing +1 later?
//ms tc            tempvarc=tempvar;
            /* printf("%s %e %f\n", "perturb p=", P[zbin-k], tempvarc); */
//ms tc            while (k+kk<=zbin && isconv1[k+kk] == 1) { /* also perturb all subsequent convective layer */
//ms tc                tempvarc  = tempvarc*pow(pl[zbin-k-kk]/pl[zbin-k-kk+1],lapse1[k+kk]);
                /* printf("%s %e %f\n", "convective perturb p=", P[zbin-k-kk], tempvarc); */
//ms tc                Tvar[k+kk][radcount] = Tvar[k+kk][radcount]+tempvarc;
//ms tc                kk = kk+1;
//ms tc            }
    }
    Tvar[0][1] = T[zbin]+0.5*tempvar; //TOA boundary changes half of TOA layer change
    /*for (j=0; j<=zbin; j++) {
        for (k=0; k<=nrl+1; k++) {
            printf("%f\t",Tvar[j][k]);
        }
        printf("\n");
    }*/

    } //END OF !TIME_STEPPING

//markus2021: Here starts radiative transfer:
	double w[zbin+1], wa[zbin+1], ws[zbin+1], g[zbin+1], tau[zbin+1], tc[zbin+1], gc[zbin+1], wc[zbin+1];
	double *gamma1, *gamma2, *gamma3, *gamma4, *lambda, *gamma, *e1, *e2, *e3, *e4, *CP0, *CP1, *CM0, *CM1, *CPh, *CMh;
	double B0, B1, miu1=0.5;
	int    ntau, nvector;
	double *A, *AS, *B, *D, *DS, *E, *ES, *Y, *Y1, *Y2, dflux[zbin+1];
	double *tc1, *gc1, *wc1, *tau1, *w1, *g1;
	double temp[zbin+1], *temp1;
	double surfaceref;
        double dt[zbin+1]; //time step if TIME_STEPPING 
        double sumBint; // check against sigma*tint^4
    
    double solarfraction;
    solarfraction = FADV/cos(THETAREF); //ms2021: global average, I guess
	
/*MS: DEBUGGING:

printf("%s\n","\n==== DEBUGGING ===="); //template
printf("%s\t%f\n","THETAREF = ",THETAREF); 
printf("%s\t%f\n","cos(THETAREF) = ",cos(THETAREF)); 
//atexit(pexit);exit(0); //ms debugging mode
*/
	/* process temperature for UV cross sections */
//markus: we only have uv sigmas between 200-300K? (see input file)
	double temperature[zbin+1], crossl;
	for (j=1; j<=zbin; j++) {
		//temperature[j] = tl[j];
		temperature[j] = 0.5*(T[j]+T[j-1]); //ms2023 to be consistent through iterations
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
    NetFlux = dmatrix(0,zbin,0,jTvar);
    
	for (j=0; j<=zbin; j++) {
        for (k=0; k<=jTvar; k++) {
            NetFlux[j][k] = 0.0;
        }
	}
    
    //Up/Dn fluxes for diagnostics
    double Fup[zbin+1],Fdn[zbin+1];
	for (j=0; j<=zbin; j++) 
        {
            Fup[j] = 0.0;
            Fdn[j] = 0.0;

            dt[j] = 1.0e-1; //for TIME_STEPPING: init in some random unit to ...
        }
    
    double radiationI0, radiationI1, radiationO;
    radiationI0=0;
    radiationI1=0;
    radiationO=0;
    
//    FILE *fp, *fopa;
//    fp=fopen("AuxillaryOut/checkemission.dat","w");
//    fopa=fopen("AuxillaryOut/checkopacity.dat","w");


//========================================================    
//====== Wavelength loop =================================	
        for (i=0; i<(NLAMBDA-1); i++) {
//	for (i=0; i<1; i++) {
//========================================================    
//========================================================    


		/* printf("%s %f %f %f %f\n", "crossa", crossa[1][i], cross[2][i], sinab[1][i], sinab[2][i]); */
		//printf("%s\t%f\n","wavelength",wavelength[i]);
		//printf("%s%d\t%s%d\n","lambda bin ",i," of ",NLAMBDA-1);
		
		/* Single Scattering Albedo */	
		for (j=1; j <= zbin; j++) {
            wa[j] = 0;
            ws[j] = 0;
            w[j]  = 0;
            j1    = zbin+1-j;

/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 *The following assumes "perfectly correlated" absorbtion features!
 *See Scheucher+2020 (doi:10.3847/1538-4357/ab9084) for a more generalized aproach assuming "random overlap"
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/

            /*for (k=1; k<=nump; k++) {
                crossl = cross[k][i] + crosst[k][i]* ( temperature[j1] - 295.0 ) ;
                wa[j] += xx[j1][stdn[k]]*crossl/qysum[k];
            }*/
            wa[j] += opacH2O[j1][i]*xx[j1][7];      //
            wa[j] += opacCH4[j1][i]*xx[j1][21];     //
            wa[j] += opacNH3[j1][i]*xx[j1][9];      //
            wa[j] += opacCO2[j1][i]*xx[j1][52];     //
            //wa[j] += opacO2[j1][i]*xx[j1][54];    //
            //wa[j] += opacOH[j1][i]*xx[j1][4];     //
            //wa[j] += opacH2CO[j1][i]*xx[j1][22];  //
            //wa[j] += opacH2O2[j1][i]*xx[j1][6];   //
            //wa[j] += opacHO2[j1][i]*xx[j1][5];    //
            wa[j] += opacCO[j1][i]*xx[j1][20];      //
            //wa[j] += opacO3[j1][i]*xx[j1][2];     //
            //wa[j] += opacC2H2[j1][i]*xx[j1][27];  //
            //wa[j] += opacC2H4[j1][i]*xx[j1][29];  //
            //wa[j] += opacC2H6[j1][i]*xx[j1][31];  //
            //wa[j] += opacCH2O2[j1][i]*xx[j1][23]; //
            //wa[j] += opacHCN[j1][i]*xx[j1][37];   //
            //wa[j] += opacN2O[j1][i]*xx[j1][11];   //
            //wa[j] += opacN2[j1][i]*xx[j1][55];    //
            //wa[j] += opacNO[j1][i]*xx[j1][12];    //
            //wa[j] += opacNO2[j1][i]*xx[j1][13];   //
            //wa[j] += opacHNO3[j1][i]*xx[j1][18];  //
            //wa[j] += opacH2S[j1][i]*xx[j1][45];   //
            //wa[j] += opacSO2[j1][i]*xx[j1][43];   //
            //wa[j] += opacOCS[j1][i]*xx[j1][49];   //
            
            wa[j] += H2H2CIA[j1][i]*xx[j1][53]*xx[j1][53]; /* H2H2CIA */
            wa[j] += H2HeCIA[j1][i]*xx[j1][53]*heliumnumber[j1]; /* H2HeCIA */
            //wa[j] += H2HCIA[j1][i]*xx[j1][53]*xx[j1][3]; /* H2HCIA */
            //wa[j] += N2N2CIA[j1][i]*xx[j1][55]*xx[j1][55]; /* N2N2CIA */
            //wa[j] += N2H2CIA[j1][i]*xx[j1][55]*xx[j1][53]; /* N2H2CIA */
            //wa[j] += CO2CO2CIA[j1][i]*xx[j1][52]*xx[j1][52]; /* CO2CO2CIA */
            
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

            //ms:HELIOS GJ1214 Comparison:
            //w[j] = 0.5;
            //w[j] = 0.0; //testing pure absorption

        }
		
/*printf("%s\n","\n==== w_0 ===="); //template
for (i=0;i<=zbin;i++)
{
    printf("%s %d\t%.3e\n","w_0 ",i,w[i]);
}*/
		/* Asymmetry Factor */
/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 * ms_two_str can deal with scattering of larger molecules. 
 * g0 needs to be calculated consistently (so does ws above)
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/
                tau[0] = 0.0;    //ms2022: just to make sure
                TAUdoub[0] = 0.0; //ms2023: double grid
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
                        //g[j] = 1.0; //testing pure absorption
		}
		
		/* Optical Depth of Each Layer */
		for (j=1; j <=zbin; j++) {
                    j1   = zbin+1-j;
                    tau[j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-P[j1])*1.0E-4;
                    TAUdoub[2*j-1] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(pl[j1]-P[j1])*1.0E-4; //ms2023: double grid
                    TAUdoub[2*j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-pl[j1])*1.0E-4; //ms2023: double grid
			/* printf("%f %e %f %f\n", wavelength[i], tau[j], w[j], g[j]); */
		}
/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 *if any tau[j] >> 1.0, radiative fluxes may become diffusion limited
 *in such cases a control run with modified grid may be performed.
 *
 *At this point it is worth checking also for Photon deposition depth (region absorbing a majority of stellar radiation)
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/
		
//========================================================    
//====== 2-STREAM switch =================================	
    if (TWO_STR_SOLVER == 1 ) ms_two_str_solver( i, w, g, tau, nrl, isconv, Tvar,jTvar, P, NetFlux,Fup, Fdn, Tint, sumBint);
    //else toon_two_str_solver( i, w, g, tau, nrl, isconv, Tvar, jTvar, P, NetFlux,Fup, Fdn, Tint, sumBint);
    else {

        /* Delta Adjustment */
    //ms22: this follows Joseph76, eqs. 5b,13&14
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
		
/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 *the following should depend on stellar spectrum? 
 *and assumed surface albedos f[nu]?
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/
		if (wavelength[i]>3000.0) {
			surfaceref = 1.0-PSURFEM;
		} else {
			surfaceref = PSURFAB;
		}
//----------------------------------------------------                
//markus: I guess here starts the 2-stream with tri-diagonal solver:
//----------------------------------------------------                
        j=1;
        while (j<zbin) { /* calculate all levels */
            j++; //ms: just some slow version of setting j=zbin ...
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
						
			gamma1 = dvector(1,ntau); //ms: closures
			gamma2 = dvector(1,ntau); //ms: closures
			gamma3 = dvector(1,ntau); //ms: closures
			gamma4 = dvector(1,ntau); //ms: closures
			gamma  = dvector(1,ntau); //ms: for solving Toon89 eqs. 19&20 / 31&32
			lambda = dvector(1,ntau); //ms: for solving Toon89 eqs. 19&20 / 31&32
			e1     = dvector(1,ntau); //ms: to populate A,B,D,E
			e2     = dvector(1,ntau); //ms: to populate A,B,D,E
			e3     = dvector(1,ntau); //ms: to populate A,B,D,E
			e4     = dvector(1,ntau); //ms: to populate A,B,D,E
			CP0    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
			CP1    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
			CM0    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
			CM1    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
                        CPh    = dvector(1,ntau); //ms: Sources in eqs. 23&27
                        CMh    = dvector(1,ntau); //ms: Sources in eqs. 24&27
//ms: For 2N eqs. in Tridiagonal matrix: A_l Y_(l-1) + B_l Y_l +D_l Y_(l+1) = E_l
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
				gamma3[j] = 0.5-g1[j]*cos(THETAREF); //ms: haven't seen this one before.
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
//----------------------------------------------------                
//ms: All after Toon+ 1989:
//----------------------------------------------------                
	for (j=1; j <= ntau; j++) {
		lambda[j] = sqrt(gamma1[j]*gamma1[j]-gamma2[j]*gamma2[j]); //ms: eq.21
		gamma[j]  = gamma2[j]/(gamma1[j]+lambda[j]); //ms: eq. 22
		e1[j] = 1.0+gamma[j]*exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		e2[j] = 1.0-gamma[j]*exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		e3[j] = gamma[j]+exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		e4[j] = gamma[j]-exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		/* printf("%d %f %f %f %f\n", j, e1[j], e2[j], e3[j], e4[j]); */
	}
        
        /*fprintf(fp, "%e\t%e\t", wavelength[i], solar[i]); */
        
        /* Loop the temperature variation */
        for (kk=0; kk<=jTvar; kk++) {
            
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

                    //ms: why += instead of =?
		    CP0[j] += 2*PI*miu1*(B0+B1/(gamma1[j]+gamma2[j]));
		    CP1[j] += 2*PI*miu1*(B0+B1*(tau1[j]+1.0/(gamma1[j]+gamma2[j])));
                
                    CPh[j] += 2*PI*miu1*(B0+B1*(tau1[j]*0.5+1.0/(gamma1[j]+gamma2[j]))); //ms: eq. 27

		    CM0[j] += 2*PI*miu1*(B0-B1/(gamma1[j]+gamma2[j]));
		    CM1[j] += 2*PI*miu1*(B0+B1*(tau1[j]-1.0/(gamma1[j]+gamma2[j])));
                
                    CMh[j] += 2*PI*miu1*(B0+B1*(tau1[j]*0.5-1.0/(gamma1[j]+gamma2[j]))); //ms: eq. 27
                } else {
               
		    CP0[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
		    CP1[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                        
                    CPh[j] += solarfraction*solar[i]*w1[j]*exp(-(tc1[j-1]+0.5*tau1[j])/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF)); //ms:eq. 23
				
                    CM0[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
		    CM1[j] += solarfraction*solar[i]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                                
                    CMh[j] += solarfraction*solar[i]*w1[j]*exp(-(tc1[j-1]+0.5*tau1[j])/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF)); //ms:eq.24
                    
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
                                //markus2021: possible error in Toon89: -/+ sign changed:
				E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) - e4[j+1]*(CM0[j+1]-CM1[j]);
				//E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) + e4[j+1]*(CM0[j+1]-CM1[j]);
			}
			l = nvector;
            surfaceref = 1.0;
			A[l] = e1[ntau] - surfaceref*e3[ntau];
			B[l] = e2[ntau] - surfaceref*e4[ntau];
			D[l] = 0;
			E[l] = -CP1[ntau] + surfaceref*CM1[ntau];
                        //ms: not sure if we can add Planck (T_Surf) and Direct Stellar beam together in these equations.
			E[l] += (1-surfaceref)*PI*2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/temp1[ntau])-1.0)*1.0E-9;
            E[l] += PI*2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/Tint)-1.0)*1.0E-9;
			/* printf("%s %f %e\n","surf temp",temp1[ntau-1],tc1[ntau]); */
			E[l] += surfaceref*cos(THETAREF)*solar[i]*solarfraction*exp(-tc1[ntau]/cos(THETAREF));
			/* for (l=1; l<=nvector; l++) {
             printf("%s %f %d %f %f %f %f\n", "ABDE", wavelength[i], l, A[l], B[l], D[l], E[l]);
             } */
            
            /* solve equations */
//ms: solving the tridiagonal matrix
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
            
//ms: What's wrong with aflux below + direct flux, as eqs. 48-51 in Toon+1989?
//
            /* obtain net flux at zero and midpoint */
            /*for (j=1; j<=ntau; j++) {
                aflux[j] = Y1[j]*(e1[j]-e3[j])+Y2[j]*(e2[j]-e4[j])+CP1[j]-CM1[j];
            }*/
            dflux[0] = Y1[1]*e3[1] - Y2[1]*e4[1] + CP0[1]; //ms:should be ok, but let's write it out:
            //dflux[0] = Y1[1]*(e1[1]-e3[1])+Y2[1]*(e2[1]-e4[1])+CP0[1]-CM0[1];
            for (j=1; j<=ntau; j++) {
                dflux[j] =  2.0*Y2[j]*(1.0-gamma[j])*exp(-lambda[j]*tau1[j]*0.5) + CPh[j] - CMh[j]; //ms:should be ok, but let's write it out:
                //dflux[j] = Y1[j]*(e1[j]-e3[j])+Y2[j]*(e2[j]-e4[j])+CP1[j]-CM1[j];
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

//----------------------------------------------------                
//markus: I guess here ends the 2-stream flux calculation
//----------------------------------------------------                

            
            NetFlux[0][kk] += dflux[0]*(wavelength[i+1]-wavelength[i]); //ms: could be within loop...
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
//moved downward                radiationO += dflux[0]*(wavelength[i+1]-wavelength[i]); /* TOA net outgoing flux */
                //ms: actually sth like a version of outgoing, see l:494
            }
/*            
printf("%s\n","\n==== CHECKING RENYU's 2-STREAM ===="); //template
printf("%s\t%d\n","kk = ",kk); //template
printf("%s\n","\n==== (i) X Y Z S ===="); //template
for (i=1;i<=2*zbin;i++)
{
    printf("%d\t%.3e\t%.3e\t%.3e\t%.3e\n",i,A[i],B[i],D[i],E[i]);
}
printf("%s\n","\n==== Up/Dn FLUXES ===="); //template
for (i=0;i<zbin;i++)
{
    printf("%s %d\t%.3e\n","Fup - Fdn ",i,dflux[i]);
}
printf("%s\n","\n==== NET FLUXES ===="); //template
for (i=1;i<=nrl;i++)
{
    printf("%s %d\t%.3e\n","NetFlux ",i,NetFlux[i][kk]);
}
*/
			
        } // END Tvar loop
        
            /*markus2021 */
            //if (i==4530 || i==4531) {
            /*printf("%s%i\t","i= ",i);
            printf("%s\t", "tau1:");
            printf("%f\n",tau1[ntau]);
        if (i==15998) {
            printf("%s%i\n","i= ",i);
            printf("%s\t", "tau1:");

            //printf("%f\n",tau1[ntau]);
            for (j=1; j<=ntau; j++) {
                printf("%f\t", tau1[j]);
            }
            printf("\n"); 

            printf("%s%f\n", "B1: ",B1);
            
            printf("%s\t", "E:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", E[l]);
            }
            printf("\n"); 

            printf("%s\t", "AS:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", AS[l]);
            }
            printf("\n"); 

            printf("%s\t", "DS:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", DS[l]);
            }
            printf("\n"); 

            printf("%s\t", "Y:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", Y[l]);
            }
            printf("\n"); 

            printf("%s\t", "dflux:");
            //printf("%f\n", dflux[ntau]);
            for (j=0; j<=ntau; j++) {
                printf("%f\t", dflux[j]);
            }
            printf("\n"); 

        }*/
            

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
        
//========================================================    
//====== END of 2-STREAM switch ==========================	
    }
//========================================================    
//========================================================   
    
//atexit(pexit);exit(0); //ms debugging mode

//========================================================    
//====== END of Wavelength loop ==========================	
    } 
//========================================================    
//========================================================   

//    fclose(fp);
//    fclose(fopa);


printf("%s\n","Wavelength-loop DONE");

radiationO = NetFlux[0][0]; /* TOA net outgoing flux */

//MS: DEBUGGING:

//printf("%s\n","\n==== DEBUGGING ===="); //template
/*printf("%s\t%f\n","wavelength <= ",wavelength[i]);
printf("%s\n","\n==== NET FLUXES ===="); //template
printf("%s\t%s\n","j   ","NetFlux[:][j] ");
for (j=0;j<=nrl;j++)
{   
    printf("%d\t",j);
    for (k=0;k<=nrl;k++)
    {
        printf("%.2e\t",NetFlux[k][j]);
    }
    printf("\n");

}
printf("%s\n","\n==== END FLUXES ===="); //template
*/
//printf("%s\n","\n==== Fup/Fdn integrated ===="); //template
//for (i=0;i<=zbin;i++)
//{
//    printf("%s %d\t%.3e\t%.3e\n","Fup,Fdn ",i,Fup[i],Fdn[i]);
//}
//printf("%s\n","==== END FLUXES ===="); //template

//atexit(pexit);exit(0); //ms debugging mode
//
//check Tint ve Bint:
//    printf("%s %e %s %e %s %e\n","Energy from Interior:",SIGMA*pow(Tint,4.0)," SUM(pi*Bint):",PI*sumBint, " DELTA:",SIGMA*pow(Tint,4.0) - PI*sumBint);
//    printf("%s %e\n", "For comparison, NetFlux[Surf] =",NetFlux[nrl][0]);

    /* Obtain the Flux and the Flux Jacobian, from bottom to up */
    for (j=0; j<=zbin; j++) {
        Rflux[j] = NetFlux[j][0] - SIGMA*pow(Tint,4.0); 
    }
    
printf("%s\n","\n==== RFLUX ===="); //template
    printf("%e\t%f\n",P[zbin],Rflux[0]);
    for (k=1; k<=zbin; k++) {
            printf("%e\t%f\n",P[zbin-k],Rflux[k]);
    }
printf("%s\n","==== END RFLUX ===="); //template

//=====================================
//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(!TIME_STEPPING) 
{//ms22tss
//=====================================
//=====================================

    double **jmax;
    jmax=dmatrix(1,zbin+1,1,zbin+1);
    for (j=0; j<=zbin; j++) {
        for (kk=0; kk<jTvar; kk++) {
            jmax[j+1][kk+1] = (NetFlux[j][kk+1] - NetFlux[j][0])/tempvar;
        }
    }
    
    /*for (j=0; j<=nrl; j++) {
        for (k=0; k<=5; k++) {
            printf("%f\t",NetFlux[j][k]);
        }
        printf("\n");
    } */
    
    /*for (j=1; j<=nrl+1; j++) {
        for (k=1; k<=nrl; k++) {
            printf("%e\t",jmax[j][k]);
        }
        printf("\n");
    }*/ 
    
    /* Solve Equation */
    int *indx;
    double ddd;
    double *deltaT;
    double relaxationfactor; //Renyu2021
    indx=ivector(1,zbin+1);
    deltaT=dvector(1,zbin+1);

    for (j=0; j<=zbin; j++) {
        deltaT[j+1] = -Rflux[j];
    }

    //markus2021:
/*    printf("%s\t","Jmax[] = ");
    for (j=1; j<=jTvar; j++) {
        for (k=1; k<=jTvar; k++) {
            printf("%f\t",jmax[j][k]);
        }
    }
    printf("\n");
*/
//----------------------------------------------------                
//markus: this solves for deltaT from flux jacobian: 
//----------------------------------------------------                
//ms:rfluxconv    
    ludcmp(jmax,zbin+1,indx,&ddd); 
    lubksb(jmax,zbin+1,indx,deltaT);
    
//atexit(pexit);exit(0); //ms debugging mode
    for (j=0; j<=zbin; j++) {
        printf("%s %d %f\n","deltaT ",j+1,deltaT[j+1]);
    } 
    
    /* calculate the new temperature profile */
//Renyu2021:    
    radcount = 0;
    relaxationfactor = R_RELAX;
    for (k=1; k<=zbin+1; k++) {
            relaxationfactor=fmin(relaxationfactor,1e-2*Tvar[k][0]/fabs(deltaT[k])); /* do not allow temperature to change more than 10% */
            //printf("%s\t%f\n", "Relaxation factor is", relaxationfactor);
    }
    for (j=1; j<=zbin+1; j++) {
//Renyu2021:        deltaT[j] *= R_RELAX; /* relaxation */
//ms:rfluxconv
        //relaxationfactor=fmin(1.0,R_RELAX*Tvar[j][0]/fabs(deltaT[j])); // for each layer separately
        deltaT[j] *= relaxationfactor; /* relaxation */
        Tvar[j][0] = Tvar[j][0] + deltaT[j];
    }
//    Tvar[0][0] = Tvar[0][0] + deltaT[1];
//ms tc    for (k=1; k<=zbin+1; k++) {
//ms tc            tempvarc = deltaT[k];
//ms: correct for temperature spikes: 
/*            if(fabs(Tvar[k][0]+tempvarc-Tvar[k-1][0]) > ) 
            {
                tempvarc = ;
            }*/
//ms: end correction
            //if (k<zbin) Tvar[k][0] = Tvar[k][0] + tempvarc; //ms22 treat surface separately below
            //if (k == zbin && isconv1[k] == 0) Tvar[k][0] = Tvar[k-1][0]; //ms22 for deep atmospheres BOA temp can be set to last atmo layer temp before convection
//ms tc            Tvar[k][0] = Tvar[k][0] + tempvarc;
//end
            /* printf("%s %e %f\n", "update p=", P[zbin-k], tempvarc); */
//            kk=1;
//            while (k+kk<=zbin && isconv1[k+kk] == 1) { /* also perturb all subsequent convective layer */
//                tempvarc  = tempvarc*pow(P[zbin-k-kk]/P[zbin-k-kk+1],lapse1[k+kk]);
//                Tvar[k+kk][0] = Tvar[k+kk][0] + tempvarc;
                /* printf("%s %e %f\n", "convective update p=", P[zbin-k-kk], tempvarc); */
//                kk = kk+1;
//            }
//        }
//ms tc    }

    free_dmatrix(jmax,1,zbin+1,1,zbin+1);
    free_dvector(deltaT,1,zbin+1);
    free_ivector(indx,1,zbin+1);
//=====================================
//=====================================
} //eo ms22tss: no time stepping
    free_dmatrix(NetFlux,0,zbin,0,jTvar);
//=====================================
//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(TIME_STEPPING)//OJO: NOT UPDATED after TVAR is tl update
{//ms22tss
//This is experimental at the moment (August 2022). Need to find optimal dt or dTemp ratio from TOA to BOA so that
// (A) upper atmosphere converges quickly enough
// (B) lower atmosphere converges quickly enough
// (C) no spurious oscillations occur numerically
// (D) convergence of whole atmosphere occurs in a reasonable runtime frame
// (E) it works for thin terrestrial (e.g. Mars-like) as well as deep gas giant (e.g. Jupiter-like) atmospheres
// (F) a p[BOA] or T[BOA] based (automatic)switch between (E) cases is conveivable

/*printf("%s\n","\n==== TVAR ===="); //template
printf("%s\t%s\n","jTvar   ","TVAR[:][jTvar] ");
for (j=0;j<=jTvar;j++)
{   
    printf("%d\t",j);
    for (k=0;k<=nrl;k++)
    {
        printf("%.2e\t",NetFlux[k][j]);
    }
    printf("\n");

}
printf("%s\n","\n==== TVAR ===="); //template
*/

    double drfluxmax, drflux[zbin+1], dTemp, f0, f1, f2;
    drfluxmax = 0.0;
/*printf("%s\n","\n==== FLUXES ===="); //template
for (i=0;i<=zbin;i++)
{
    //printf("%s %d\t%.3e\n","Netflux ",i,Fup[i]-Fdn[i]);
    if(i>0) printf("%s %d\t%.3e\n","dNetflux ",i,Fup[i]-Fdn[i]-Fup[i-1]+Fdn[i-1]);
    //printf("%s %d\t%.3e\n","Rflux ",i,Fup[i]-Fdn[i]-SIGMA*pow(Tint,4.0));
}
printf("%s\n","==== END FLUXES ===="); //template
*/
    for (j=1; j<=zbin; j++)
    {
        drflux[j] = Fup[j-1]-Fdn[j-1]-Fup[j]+Fdn[j];
        drfluxmax = fmax(drfluxmax,fabs(drflux[j]) );
    }
    printf("%s %f\n","dRFLUX_max= ",drfluxmax);
    if (rtstepcount==1) rt_drfluxmax_init=drfluxmax;//ms: let's scale things as a ratio to initial fluxes

//    for (j=1; j<zbin; j++)//ms22 treat surface separately below
    for (j=1; j<zbin; j++)//surface treatment included

    {
        dt[j] = cp[zbin+1-j]/GA/SIGMA/pow(Tvar[j][0],3.0)*P[zbin-j];
        dt[j] *= 10.E0/pow(abs(drflux[j]),0.9);
//        dt[j] = 9E+2; //as in old 1D-TERRA
//        if (!isconv1[j])
        //dT = dt * g/cp * dFnet/dP
        //since Tvar are boundary temperatures, both boundaries are affected:
        //dTemp = fmin(1e-2*Tvar[j-1][0]*drflux[j]/fabs(drflux[j]), fabs(dt[j] * (GA/cp[zbin+1-j]) * drflux[j] / (P[zbin-j]-P[zbin+1-j])));
        //Tvar[j-1][0] = Tvar[j-1][0] - drflux[j]/fabs(drflux[j])*dTemp ;
        //Tvar[j][0] = Tvar[j][0] - drflux[j]/fabs(drflux[j])*dTemp ;
        
        //Tvar[j-1][0] = Tvar[j-1][0] - drflux[j]/fabs(drflux[j]) * fmin(fabs(dt[j]*drflux[j]/drfluxmax/pow(P[zbin-j]-P[zbin+1-j],0.5)),1E-2*Tvar[j-1][0]);
        //Tvar[j][0] = Tvar[j][0] - drflux[j]/fabs(drflux[j]) * fmin(fabs(dt[j]*drflux[j]/drfluxmax/pow(P[zbin-j]-P[zbin+1-j],0.5)),1E-2*Tvar[j][0]);
        
//testing something unphysical for stability reasons:  
//-----------------------------------------------------------
        dt[j] = 1.0; //just a prefactor, has little to do with actual time
        f0 = drflux[j]/fabs(drflux[j]);
        //f1 = drflux[j]/drfluxmax;
        f1 = drflux[j];
        //f1 = pow(drflux[j],0.1);
        //f2 = fabs(drfluxmax/rt_drfluxmax_init);
        //f2 = pow(drfluxmax/rt_drfluxmax_init,4.0);
        f2 = 0.75;
        //f2=1.0;
        //f2=0.0;
        Tvar[j-1][0] = Tvar[j-1][0] - f0 * fmin(fabs(dt[j]* f1 /pow(P[zbin-j]-P[zbin+1-j],1.0 - f2 )),1E-2*Tvar[j-1][0]);
        Tvar[j][0] = Tvar[j][0] - f0 * fmin(fabs(dt[j]* f1 /pow(P[zbin-j]-P[zbin+1-j],1.0- f2 )),1E-2*Tvar[j][0]);
//-----------------------------------------------------------
        
//testing HELIOS implementation:  
//-----------------------------------------------------------
        //dt[j] = 1.0; //should be 0.5 * 1e5
        //Tvar[j-1][0] = Tvar[j-1][0] + pow(drflux[j],0.1)/pow(Tvar[j-1][0],3.0)*dt[j]/SIGMA/log(P[zbin+1-j]/P[zbin-j]);
        //Tvar[j][0] = Tvar[j][0] + pow(drflux[j],0.1)/pow(Tvar[j][0],3.0)*dt[j]/SIGMA/log(P[zbin+1-j]/P[zbin-j]) ;
       
    }
//    Tvar[zbin][0] = Tvar[zbin-1][0];    //ms22 for deep atmospheres BOA temp can be set to last atmo layer temp before convection
}//eo ms22tss

    
    for (j=0; j<=zbin; j++) {
        tempcnew[j] = fmin(10000.0,fmax(Tvar[zbin+1-j][0],1.0)); //ms22: also include max temp
    }
    
    free_dmatrix(Tvar,0,zbin+1,0,jTvar);
    
    /*printf("%s\t%f\n", "Top-of-Atmosphere incoming radiation flux is", radiationI0);
    printf("%s\t%f\n", "Bottom-of-Atmospehre incoming radiation flux is", radiationI1);*/
    printf("%s\t%f\n", "TOA outgoing net radiation flux is", radiationO);
	
//atexit(pexit);exit(0); //ms debugging mode
}
