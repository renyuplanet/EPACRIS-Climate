/*Function to calculate the reflected solar light */
/*
 * Original Author: Renyu Hu (renyu.hu@jpl.nasa.gov)
 * Version: 2021, 2-Stream update to Heng+2019 flux equations including direct stellar beam
 * Editor v2021: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Modifier symbols: 
    * 'ms_' ... in file names
    * 'markus<date>' ... in description lines/blocks 
    * 'ms:' ... in-line comments
    * 'OJO' ... warnings to be addressed before public release
*/

// Include header files
#include <math.h>
#include "constant.h"
#include <stdio.h>
#include "config.h"
#include "global_temp.h"
#include "nrutil.h"
#include "lu_decomp.h"  // LU decomposition functions (ludcmp, lubksb)

// Include C files
#include "ms_2stream_test.c" //ms: solving for radiative fluxes using Heng+2018; Deitrick+2021
#include "toon_2stream.c"    //ms: solving for radiative fluxes using Toon+89
#include "mm_2stream.c" //ms: solving for radiative fluxes using Heng+2018; Deitrick+2021

void ms_RadTrans(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[], double dt[], int isequil[], double *radiationI0_out, double *radiationI1_out, double *radiationO_out);

void ms_RadTrans(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[], double dt[], int isequil[], double *radiationI0_out, double *radiationI1_out, double *radiationO_out)
{

    int i=0, j, k, jj, l, j1, kk;
    static int iter = 0;


    // GA is global (declared in global_temp.h, defined in epacris_main.c)
    // No local declaration needed - use the global GA directly
	double mole2dust, planck, xlll;
    // Calculate mole2dust???? what is this? particle size and density defined in config file
	mole2dust = PI*pow(AERSIZE,3)*AERDEN/6.0/AMU*2.0558;
	
    // convert to meters for toon solver? (#define TWO_STR_SOLVER  0)
	double lam1[NLAMBDA];
	for (i=0; i<NLAMBDA; i++) {
		lam1[i] = wavelength[i]*1.0E-9; /* convert to m */
	}
    
    /* Temperature variation */
    //for flux jacobian to solve for steady-state temperatures
    int nrl, jTvar; //number of radiative layers, number of temperature variables
    nrl = zbin-ncl; //number of radiative layers
    double **Tvar; // 2D array to store temperature variation? i think
    double tempvar[zbin+1]; //temperature variation
    double tempvarc; //temperature variation
    int isconv1[zbin+1],radcount; //convective layers
    double lapse1[zbin+1]; //convective layers

    // reverse the convective layers
    for (k=1; k<=zbin; k++) {
        isconv1[k] = isconv[zbin+1-k]; /* reverse, from top to bottom */
        lapse1[k] = lapse[zbin+1-k];
    }
    radcount = 0; //number of radiative layers

    jTvar = nrl+1;

    // For time stepping solver 
    if (TIME_STEPPING) jTvar = 0;
    Tvar = dmatrix(0,zbin,0,jTvar);
    for (j=0; j<=zbin; j++) {
        Tvar[j][0] = T[zbin-j];
        tempvar[j] = RJACOB_TEMPVAR;//*T[zbin-j]; //ms23
    }

    // For Jacobian solver
    if (!TIME_STEPPING)
    {
        for (j=0; j<=zbin; j++) {
            Tvar[j][1] = T[zbin-j];
        }
    Tvar[0][1] = Tvar[0][1]+tempvar[0];
    
//    Tvar[0][1] = Tvar[0][1]*(1+tempvar); //ms22: testing

    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            for (j=0; j<=zbin; j++) {
                Tvar[j][radcount+1] = T[zbin-j]; //ms2021: waste of computation within k loop
            }
            Tvar[k][radcount+1] = Tvar[k][radcount+1]+tempvar[k];
//            Tvar[k][radcount+1] = Tvar[k][radcount+1]*(1+tempvar);
            kk=1; //ms: longer than just writing +1 later?
            tempvarc=tempvar[k];
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

    } //END OF !TIME_STEPPING

//markus2021: Here starts radiative transfer:
	double w[zbin+1], wa[zbin+1], ws[zbin+1], g[zbin+1], tau[zbin+1], tc[zbin+1], gc[zbin+1], wc[zbin+1];
	double *gamma1, *gamma2, *gamma3, *gamma4, *lambda, *gamma, *e1, *e2, *e3, *e4, *CP0, *CP1, *CM0, *CM1, *CPh, *CMh;
	double B0, B1, miu1=0.5;
	int    ntau, nvector;
	double *A, *AS, *B, *D, *DS, *E, *ES, *Y, *Y1, *Y2, dflux[zbin+1];
	double *tc1, *gc1, *wc1, *tau1, *w1, *g1, *temp1;
	double temp[zbin+1];
	double surfaceref;
        double sumBint; // check against sigma*tint^4
        double dFnetC[zbin+1]; //ms23: taking time stepping approach to jacobian solver
    

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
    NetFlux = dmatrix(0,nrl,0,jTvar);
    
	for (j=0; j<=nrl; j++) {
        for (k=0; k<=jTvar; k++) {
            NetFlux[j][k] = 0.0;
        }
	}
    
    //Up/Dn fluxes for diagnostics
//    double Fup[zbin+1],Fdn[zbin+1],Fcup[zbin+1],Fcdn[zbin+1];
    double *Fup,*Fdn,*Fcup,*Fcdn;
    Fup = dvector(0,zbin);
    Fdn = dvector(0,zbin);
    Fcup = dvector(0,zbin);
    Fcdn = dvector(0,zbin);
	for (j=0; j<=zbin; j++) 
        {
            Fup[j] = 0.0;
            Fdn[j] = 0.0;
            Fcup[j] = 0.0; //zero element unused
            Fcdn[j] = 0.0; //zero element unused
        }
    
    double radiationI0, radiationI1, radiationO;
    radiationI0=0; //TOA incoming flux [W/m²]
    radiationI1=0; //BOA incoming flux [W/m²]
    radiationO=0; //TOA NET flux (upward-downward) [W/m²]
    
//    FILE *fp, *fopa;
//    fp=fopen("AuxillaryOut/checkemission.dat","w");
//    fopa=fopen("AuxillaryOut/checkopacity.dat","w");


//========================================================    
//====== Wavelength loop =================================	
	//DEBUGGprintf("%s\t%d\n","NLAMBDA ",NLAMBDA);
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

            // --- UV-VIS cross sections used in PhotoChem ---
            for (k=1; k<=nump; k++) {
                // crossl = cross[k][i] + crosst[k][i]* ( temperature[j1] - 295.0 ) ;
                crossl = cross[k][i];
                //wa[j] += xx[j1][stdn[k]]*crossl/qysum[k]; stdn not found anywhere, changed out with stdcross?
                wa[j] += xx[j1][stdcross[k]]*crossl/qysum[k];
            } 

            // --- IR Opacities ---
            wa[j] += opacH2O[j1][i]*xx[j1][7];      //
            wa[j] += opacCH4[j1][i]*xx[j1][21];     //
            wa[j] += opacNH3[j1][i]*xx[j1][9];      //
            wa[j] += opacCO2[j1][i]*xx[j1][52];     //
            wa[j] += opacO2[j1][i]*xx[j1][54];    //
            wa[j] += opacOH[j1][i]*xx[j1][4];     //
            wa[j] += opacH2CO[j1][i]*xx[j1][22];  //
            wa[j] += opacH2O2[j1][i]*xx[j1][6];   //
            wa[j] += opacHO2[j1][i]*xx[j1][5];    //
            wa[j] += opacCO[j1][i]*xx[j1][20];      //
            wa[j] += opacO3[j1][i]*xx[j1][2];     //
            wa[j] += opacC2H2[j1][i]*xx[j1][27];  //
            wa[j] += opacC2H4[j1][i]*xx[j1][29];  //
            wa[j] += opacC2H6[j1][i]*xx[j1][31];  //
            wa[j] += opacCH2O2[j1][i]*xx[j1][23]; //
            wa[j] += opacHCN[j1][i]*xx[j1][37];   //
            wa[j] += opacN2O[j1][i]*xx[j1][11];   //
//            wa[j] += opacN2[j1][i]*xx[j1][55];    //
            wa[j] += opacNO[j1][i]*xx[j1][12];    //
            wa[j] += opacNO2[j1][i]*xx[j1][13];   //
            wa[j] += opacHNO3[j1][i]*xx[j1][18];  //
            wa[j] += opacH2S[j1][i]*xx[j1][45];   //
            wa[j] += opacSO2[j1][i]*xx[j1][43];   //
            wa[j] += opacOCS[j1][i]*xx[j1][49];   //
            
            // --- CIA ---
            wa[j] += H2H2CIA[j1][i]*xx[j1][53]*xx[j1][53]; /* H2H2CIA */
            // wa[j] += H2HeCIA[j1][i]*xx[j1][53]*heliumnumber[j1]; /* H2HeCIA */
            wa[j] += H2HCIA[j1][i]*xx[j1][53]*xx[j1][3]; /* H2HCIA */
            wa[j] += N2N2CIA[j1][i]*xx[j1][55]*xx[j1][55]; /* N2N2CIA */
            wa[j] += N2H2CIA[j1][i]*xx[j1][55]*xx[j1][53]; /* N2H2CIA */
            wa[j] += CO2CO2CIA[j1][i]*xx[j1][52]*xx[j1][52]; /* CO2CO2CIA */
            
            // --- Rayleigh scattering cross sections ---
            ws[j] += crossr[i]*MM[j1];
            
            // --- Aerosol cross sections ---
            //wa[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*(1.0-sinab[1][i]);
            //wa[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*(1.0-sinab[2][i]);
            //ws[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*sinab[1][i];
            //ws[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*sinab[2][i]; */
           
            // --- Cloud opacities and albedos ---     
            // cH2O is already multiplied by number density (particles/m³)
            // and is in units of cm^-1 (extinction coefficient), so no need to adjust
            // *** Assymetry factor (g) is below 30 lines below ***

            // ADD SPECIES HERE
            // Absorption: extinction × (1 - albedo)
            // H2O
            wa[j] += cH2O[j1][i]*(1.0-aH2O[j1][i]);
            // NH3
            wa[j] += cNH3[j1][i]*(1.0-aNH3[j1][i]); //portion is reflected with albedo

            // Scattering: extinction × albedo
            // H2O
            ws[j] += cH2O[j1][i]*aH2O[j1][i];
            // NH3
            ws[j] += cNH3[j1][i]*aNH3[j1][i];


            // Calculate single scattering albedo w
            if (ws[j] > 0.0) {
                w[j]  = ws[j]/(wa[j]+ws[j]);
            } else {
                w[j]  = 0.0;
            }
            
            // Limit w to avoid numerical instability
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
		
        // Old Marcus Debugging:
        /*printf("%s\n","\n==== w_0 ===="); //template
        for (i=0;i<=zbin;i++)
        {
            printf("%s %d\t%.3e\n","w_0 ",i,w[i]);
        }*/


        // Asymmetry Factor
        /* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
        * ms_two_str can deal with scattering of larger molecules. 
        *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/

        tau[0] = 0.0;    //ms2022: just to make sure
        TAUdoub[0] = 0.0; //ms2023: double grid
		for (j=1; j <= zbin; j++) {
			g[j] = 0.0; 
			j1   = zbin+1-j;

            // Old Marcus dust asymmetry parameter?
			/* g[j] += crossa[1][i]*xx[j1][78]*98.0/mole2dust*sinab[1][i]*asym[1][i];
			g[j] += crossa[2][i]*xx[j1][111]*256.0/mole2dust*sinab[2][i]*asym[2][i]; */

            // Cloud asymmetry g
            // Weighted mean: g = sum(scattering_opacity × asymmetry) / sum(scattering_opacity)
            // cH2O*aH2O is scattering opacity (cm^-1), gH2O is asymmetry (dimensionless)


            // ADD SPECIES HERE
            // H2O
            g[j] += cH2O[j1][i]*aH2O[j1][i]*gH2O[j1][i];
            // NH3
            g[j] += cNH3[j1][i]*aNH3[j1][i]*gNH3[j1][i];
            
			if (ws[j] > 0.0) {
				g[j] = g[j]/ws[j];
			} else {
				g[j] = 0.0;
			}


            //Print cH2O, aH2O, gH2O only if g not zero
            // if (aH2O[j1][i] != 1.0) {
            //     printf("cH2O[%d][%d] = %f, aH2O[%d][%d] = %f, gH2O[%d][%d] = %f\n", j1, i, cH2O[j1][i], j1, i, aH2O[j1][i], j1, i, gH2O[j1][i]);
            // }



                        //g[j] = 1.0; //testing pure absorption
		}

        //Print max and minimum values of g
        //printf("Max value of g: %f, Min value of g: %f\n", fmax(g[1], g[zbin]), fmin(g[1], g[zbin]));
		/* Optical Depth of Each Layer */
		//printf("%s\t%s\t%s\t%s\t%s\t%s\n","j", "wavelength[i]", "tau[j]","TAUdoub[2*j-1]","TAUdoub[2*j]","TAUdoub[sum]"); 
		//if(i%100==0) printf("%s\t%12s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\t%16s\n","j", "wavelength[i]", "tau[j]","wa[j]","ws[j]","w0[j]","g0[j]","MM[j]","AMU","Meanmol.","GA","P[-1]","P[j]","Tvar[j][0]","solar"); 
		//if(i>=6000 && i<8000) printf("%s\t%12s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\t%26s\n","j", "wavelength[i]", "tau[j]","wa[j]","ws[j]","w0[j]","g0[j]","MM[j]","AMU","Meanmol.","GA","P[-1]","P[j]","Tvar[j][0]","solar"); 
		//if(i%100==0) printf("%s\t%s\t%16s\t%16s\t%16s\t%16s\t%16s\n","j", "wavelength[i]", "tau[j]","wa[j]","ws[j]","w0[j]","g0[j]"); 
		for (j=1; j <=zbin; j++) {
                    j1   = zbin+1-j;
                    tau[j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-P[j1])*1.0E-4;
                    TAUdoub[2*j-1] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(pl[j1]-P[j1])*1.0E-4; //ms2023: double grid
                    TAUdoub[2*j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-pl[j1])*1.0E-4; //ms2023: double grid
		    //printf("%d\t%f\t%e\t%e\t%e\t%e\n",j, wavelength[i], tau[j],TAUdoub[2*j-1],TAUdoub[2*j],TAUdoub[2*j-1]+TAUdoub[2*j]); 
		    //if(i%100==0) printf("%d\t%f\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",j, wavelength[i], tau[j], wa[j], ws[j], w[j], g[j],MM[j1],AMU,meanmolecular[j1],GA,P[j1-1],P[j1],Tvar[j][0],solar[i]); 
		    //if(i>=6000 && i<8000) printf("%d\t%f\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\t%.22e\n",j, wavelength[i], tau[j], wa[j], ws[j], w[j], g[j],MM[j1],AMU,meanmolecular[j1],GA,P[j1-1],P[j1],Tvar[j][0],solar[i]); 
		    //if(i%100==0) printf("%d\t%f %.12e %.12e %.12e %.12e %.12e\n",j, wavelength[i], tau[j], wa[j], ws[j], w[j], g[j]); 
		}
/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 *if any tau[j] >> 1.0, radiative fluxes may become diffusion limited
 *in such cases a control run with modified grid may be performed.
 *
 *At this point it is worth checking also for Photon deposition depth (region absorbing a majority of stellar radiation)
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/
		
//========================================================    
//====== 2-STREAM switch =================================	

    if (TWO_STR_SOLVER == 2 ) mm_two_str_solver( i, w, g, tau, nrl, isconv, Tvar,jTvar, P, NetFlux,Fup,Fdn,Fcup,Fcdn, Tint);
    if (TWO_STR_SOLVER == 1 ) ms_two_str_solver( i, w, g, tau, nrl, isconv, Tvar,jTvar, P, NetFlux,Fup,Fdn,Fcup,Fcdn, Tint);
    if (TWO_STR_SOLVER == 0 ) toon_two_str_solver( i, w, g, tau, nrl, isconv, Tvar, jTvar, P, NetFlux,Fup, Fdn,Fcup,Fcdn, Tint);

//========================================================    
//====== END of Wavelength loop ==========================	
    } 
//========================================================    
//========================================================   
//DEBUGG    printf("%s\n","Wavelength-loop DONE");

    // Calculate TOA incoming stellar flux (matches 2-stream solver: muD[1]*ffrac[1]*solar = FADV*solar)
    for (i=0; i<(NLAMBDA-1); i++) {
        double lambda1 = wavelength[i];
        double lambda2 = wavelength[i+1];
        double dlambda = lambda2 - lambda1;
        double solar_flux_avg = 0.5 * (solar[i] + solar[i+1]);
        radiationI0 += solar_flux_avg * dlambda * FADV;
    }

// Energy balance: F_out = F_absorbed_stellar + F_internal
// Stellar absorption controlled by FaintSun parameter in main code

radiationO = NetFlux[0][0]; /* TOA net outgoing flux */
    /* Obtain the Flux and the Flux Jacobian, from bottom to up */

    //Heng eqs.:
    if (TWO_STR_SOLVER >= 1 ) 
    {
    // Populate Rflux with netflux differences
    for (j=0; j<=nrl; j++) {
        if(j<nrl || isconv1[zbin] == 1) Rflux[j] = NetFlux[j][0] - NetFlux[j+1][0]; //now as delta netflux
        // if(j<nrl || isconv1[zbin] == 1) Rflux[j] = NetFlux[j][0] - SIGMA*pow(Tint,4.0); //now as delta netflux
        if(isnan(Rflux[j])) {printf("%s%d%s %f %s","isNaN Rflux[",j,"] ",Rflux[j]," --- set to 0.0");Rflux[j] = 0.0;}
        if (j<zbin && j>0) dFnetC[j] = Fcup[j]-Fcdn[j]-Fcup[j+1]+Fcdn[j+1]; //dFnet from fluxes at layer centers, direclty influenced by boundary temperatures
    }
    // Surface heating
    if(isconv1[zbin]== 0) {
        Rflux[nrl] = NetFlux[nrl][0] - SIGMA*pow(Tint,4.0);
        if(isnan(Rflux[nrl])) {printf("%s %f %s","isNaN Rflux[nrl] ",Rflux[nrl]," --- set to 0.0");Rflux[nrl] = 0.0;}
        dFnetC[nrl] = Fcup[nrl]-Fcdn[nrl] - SIGMA*pow(Tint,4.0); //BOA condition
    }

    }
    // Toon eqs.:
    if (TWO_STR_SOLVER == 0 ) for (j=0; j<=nrl; j++) Rflux[j] = NetFlux[j][0] - SIGMA*pow(Tint,4.0);

/*printf("%s\n","\n==== NetFLUX ===="); //template
    for (k=0; k<=nrl; k++) printf("%d\t%f\n",k,NetFlux[k][0]);
printf("%s\n","==== END NetFLUX ===="); //template
*/
//DEBUGGprintf("%s\n","\n==== RFLUX ===="); //template
    //DEBUGGprintf("%d\t%e\t%e\t%e\n",zbin,P[zbin],Rflux[0],NetFlux[0][0]);
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k] == 0) {
            radcount += 1;
            //DEBUGGprintf("%d\t%e\t%e\t%e\n",zbin-k,P[zbin-k],Rflux[radcount],NetFlux[radcount][0]);
        }
    }
//DEBUGGprintf("%s\n","==== END RFLUX ===="); //template

//=====================================
//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(!TIME_STEPPING)
{//ms22tss
//=====================================
//=====================================

    double **jmax;
    jmax=dmatrix(1,nrl+1,1,nrl+1);
    j=0;
    for (k=0; k<=zbin; k++) {
        if ((k==0 && isconv1[1] == 0) || (k>0 && isconv1[k] == 0)) {
        for (kk=0; kk<jTvar; kk++) {
            if (TWO_STR_SOLVER == 0 ) jmax[j+1][kk+1] = (NetFlux[j][kk+1] - NetFlux[j][0])/tempvar[k]; //here as netfluxes
            if (TWO_STR_SOLVER >= 1 ) 
            {
            if(j<nrl || isconv1[zbin] == 1) jmax[j+1][kk+1] = (NetFlux[j][kk+1] - NetFlux[j+1][kk+1] - NetFlux[j][0] + NetFlux[j+1][0])/tempvar[k]; //now as delta-Netfluxes through boundary
            if(j==nrl && isconv1[zbin] == 0) jmax[j+1][kk+1] = (NetFlux[nrl][kk+1] - NetFlux[nrl][0])/tempvar[k]; //now as delta-Netfluxes through boundary
            }
        }
        j+=1;
        }
    }
    
    /* Solve Equation */
    int *indx;
    double ddd;
    double *deltaT;
    double relaxationfactor; //Renyu2021
    double relaxfvector[nrl+2]; //so I don't have to define as dvector, and free later
    indx=ivector(1,nrl+1);
    deltaT=dvector(1,nrl+1);

    for (j=0; j<=nrl; j++) {
    //for (j=0; j<nrl; j++) {
        deltaT[j+1] = -Rflux[j];
    }

//----------------------------------------------------                
//markus: this solves for deltaT from flux jacobian: 
//----------------------------------------------------                
//ms:rfluxconv    
    ludcmp(jmax,nrl+1,indx,&ddd); 
    lubksb(jmax,nrl+1,indx,deltaT);
    
    // Print deltaT for jacobian solver
    // for (j=0; j<=nrl; j++) {
    //     printf("%s %d %.3e\n","deltaT ",j+1,deltaT[j+1]);
    // } 
    
    /* calculate the new temperature profile */
    radcount = 0;
    relaxationfactor = R_RELAX;
    relaxationfactor=fmin(relaxationfactor,DT_MAX*Tvar[0][0]/fabs(deltaT[1])); //ms2023: include TOA
    relaxfvector[1]=fmin(R_RELAX,DT_MAX*Tvar[0][0]/fabs(deltaT[1])); /* do not allow temperature to change more than 10% */
    for (k=1; k<=zbin; k++) { 
        if (isconv1[k]==0) { /* radiative layer*/
            radcount = radcount+1;
            relaxationfactor=fmin(relaxationfactor,DT_MAX*Tvar[k][0]/fabs(deltaT[radcount+1])); /* do not allow temperature to change more than xx% */
            relaxfvector[radcount+1]=fmin(R_RELAX,DT_MAX*Tvar[k][0]/fabs(deltaT[radcount+1])); /* do not allow temperature to change more than xx% */
            //printf("%s\t%f\n", "Relaxation factor is", relaxationfactor);
        }
    }
    for (j=1; j<=nrl+1; j++) {
        deltaT[j] *= relaxationfactor; /* relaxation */
        //deltaT[j] *= relaxfvector[j]; /* relaxation */
    }
    Tvar[0][0] = Tvar[0][0] + deltaT[1];
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            tempvarc = deltaT[radcount+1];
            Tvar[k][0] = Tvar[k][0] + tempvarc; //ms2023
            kk=1;
            while (k+kk<=zbin && isconv1[k+kk] == 1) { /* also perturb all subsequent convective layer */
                tempvarc  = tempvarc*pow(P[zbin-k-kk]/P[zbin-k-kk+1],lapse1[k+kk]);
                Tvar[k+kk][0] = Tvar[k+kk][0] + tempvarc;
                /* printf("%s %e %f\n", "convective update p=", P[zbin-k-kk], tempvarc); */
                kk = kk+1;
            }
        }
    }
    //Tvar[zbin][0] = Tvar[zbin-1][0] + (Tvar[zbin-1][0]-Tvar[zbin-2][0]); // just avoid weird surface behavior

    free_dmatrix(jmax,1,nrl+1,1,nrl+1);
    free_dvector(deltaT,1,nrl+1);
    free_ivector(indx,1,nrl+1);
//=====================================
//=====================================
} //eo ms22tss: no time stepping
//=====================================
//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(TIME_STEPPING)
{//ms22tss
//Challenge: find optimal dt or dTemp ratio from TOA to BOA so that
// (A) upper atmosphere converges quickly enough
// (B) lower atmosphere converges quickly enough
// (C) no spurious oscillations occur numerically
// (D) convergence of whole atmosphere occurs in a reasonable runtime frame
// (E) it works for thin terrestrial (e.g. Mars-like) as well as deep gas giant (e.g. Jupiter-like) atmospheres
// (F) a p[BOA] or T[BOA] based (automatic)switch between (E) cases is conveivable

    double drfluxmax, drflux[zbin+2], dTemp, f0, f1, f2, deltaT[zbin+2], fratio, change;
    double dcfluxmax, dcflux[zbin+1], dcTemp, deltaTc[zbin+1]; //add center fluxes and dTemps
    drfluxmax = 0.0;

//TS_SCHEME switch
if (TS_SCHEME == 0){ 
    for (j=1; j<=zbin; j++)
    {
        drflux[j] = Fup[j-1]-Fdn[j-1]-Fup[j]+Fdn[j];
        drfluxmax = fmax(drfluxmax,fabs(drflux[j]) );
    }
    //if (RTstepcount % PRINT_ITER == 0) printf("%s %.3e\n","dRFLUX_max= ",drfluxmax);
    if (RTstepcount==1) rt_drfluxmax_init=drfluxmax;//ms: let's scale things as a ratio to initial fluxes

    for (j=1; j<zbin; j++)//surface treatment included
    {
        //dt[j] = cp[zbin+1-j]/GA/SIGMA/pow(Tvar[j][0],3.0)*P[zbin-j];
        //dt[j] *= 10.E0/pow(abs(drflux[j]),0.9);
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
    }
}// eo TS_SCHEME=0

if (TS_SCHEME == 1){
    for (j=1; j<zbin; j++)
    {
        drflux[j] = NetFlux[j-1][0]-NetFlux[j+1][0];
        drfluxmax = fmax(drfluxmax,fabs(drflux[j]) );
        if (TWO_STR_SOLVER == 0 ) printf("%s %d %f\n","drflux ",zbin-j,drflux[j]);

        if(j<zbin){
        dcflux[j] = Fcup[j]-Fcdn[j]-Fcup[j+1]+Fcdn[j+1];
        dcfluxmax = fmax(dcfluxmax,fabs(dcflux[j]) );
        //DEBUGGif (TWO_STR_SOLVER >= 1 ) printf("%s %d %.3e\n","dcflux ",zbin-j,dcflux[j]);
        }
    }
    drflux[zbin] = NetFlux[zbin-1][0] - SIGMA*pow(Tint,4.0) ; //Surface "layer" separately
    drfluxmax = fmax(drfluxmax,fabs(drflux[zbin]));
    //if (RTstepcount % PRINT_ITER == 0) printf("%s %.3e\n","dRFLUX_max= ",drfluxmax);
    dcflux[zbin] = Fcup[zbin]-Fcdn[zbin] - SIGMA*pow(Tint,4.0) ; //Surface "layer" separately
    dcfluxmax = fmax(dcfluxmax,fabs(dcflux[zbin]));
    //if (RTstepcount % PRINT_ITER == 0) printf("%s %.3e\n","dcFLUX_max= ",dcfluxmax);
    if (RTstepcount==1) rt_drfluxmax_init=drfluxmax;//ms: let's store in case we need it

    for (j=1; j<=zbin; j++)
    {
        deltaT[j] = - dt[j] * pl[1]/(P[0]-P[1]) * drflux[j] / pow(fabs(drflux[j]),0.9); //HELIOS code
        if(fabs(drflux[j])<1.0e-20) deltaT[j]=0.0; //just avoid NAN from dFnet=0
        if (j<=zbin) {
            deltaTc[j] = - dt[j] * pl[1]/(P[0]-P[1]) * dcflux[j] / pow(fabs(dcflux[j]),0.9); //HELIOS code
            if(fabs(dcflux[j])<1.0e-20) deltaTc[j]=0.0; //just avoid NAN from dFnet=0
        }
        //check if |drflux|/sigma*T < epsval -> layer in equilibrium:
        fratio = fabs(drflux[j]) / SIGMA / pow(Tvar[j-1][0],4.0);
        if (fratio <= Tol_FRATIO) {isequil[zbin+1-j] = 1;} else {isequil[zbin+1-j] = 0;}
        if (TWO_STR_SOLVER == 0 )printf("%s %d %f \t%s %.2e\n","deltaT ",zbin-j,deltaT[j], "dF/sigma*T^4 ", fratio);
        if (j<zbin+1) fratio = fabs(dcflux[j]) / SIGMA / pow(Tvar[j][0],4.0);
        if (fratio <= Tol_FRATIO) {isequil[zbin+1-j] = 1;} else {isequil[zbin+1-j] = 0;} //testing
        //DEBUGGif (TWO_STR_SOLVER >= 1 ) if (j<zbin+1) printf("%s %d %.2e \t%s %.2e\n","deltaTc ",zbin-j,deltaTc[j], "dF/sigma*T^4 ", fratio);
    }
//atexit(pexit);exit(0); //ms debugging mode
    if (TWO_STR_SOLVER >= 1 )
    {
    Tvar[zbin][0] = Tvar[zbin][0] + deltaTc[zbin];
    for (j=zbin-1; j>0; j--)Tvar[j][0] = Tvar[j][0] + deltaTc[j];
    Tvar[0][0] = Tvar[0][0] + 2.0/3.0*deltaTc[1] + 1.0/3.0*deltaTc[2];

    /*    //advance temperatures for non-equilibrium layers
    if(isequil[0]==0) Tvar[zbin][0] = Tvar[zbin][0] + deltaTc[zbin];
    for (j=zbin-1; j>0; j--) if(isequil[zbin-j]==0) Tvar[j][0] = Tvar[j][0] + deltaTc[j];
    if(isequil[zbin]==0) Tvar[0][0] = Tvar[0][0] + 2.0/3.0*deltaTc[1] + 1.0/3.0*deltaTc[2];
*/
    }
    if (TWO_STR_SOLVER == 0 )
    {
    Tvar[zbin][0] = Tvar[zbin][0] + deltaT[zbin];
    for (j=zbin-1; j>0; j--)Tvar[j][0] = Tvar[j][0] + deltaT[j];
    Tvar[0][0] = Tvar[0][0] + 2.0/3.0*deltaT[1] + 1.0/3.0*deltaT[2];
    }


}//eo TS_SCHEME=1

}//eo ms22tss

    
    for (j=0; j<=zbin; j++) {
        tempbnew[j] = fmin(10000.0,fmax(Tvar[zbin-j][0],1.0)); //ms22: also include max temp
        //tempbnew[j] = fmax(Tvar[zbin-j][0],1.0); //ms22: also include max temp
    }
    
    free_dmatrix(Tvar,0,zbin,0,jTvar);
    free_dmatrix(NetFlux,0,nrl,0,jTvar);
    free_dvector(Fup,0,zbin);
    free_dvector(Fdn,0,zbin);
    free_dvector(Fcup,0,zbin);
    free_dvector(Fcdn,0,zbin);

    iter++;
    
    /*printf("%s\t%f\n", "Top-of-Atmosphere incoming radiation flux is", radiationI0);
    printf("%s\t%f\n", "Bottom-of-Atmospehre incoming radiation flux is", radiationI1);*/
    
    // Energy balance diagnostics
    if (iter % PRINT_ITER == 0 && TWO_STR_SOLVER >= 1) {
        double absorbed_stellar = radiationI0;  // Already scaled by FaintSun in main code
        double internal_flux = SIGMA * pow(Tint, 4.0);
        double energy_balance = radiationO - absorbed_stellar - internal_flux;
        
        printf("=======================================================\n");
        printf("Diagnostic from the ms_radtrans_test.c file:\n");
        printf("Energy Balance Check:\n");
        printf("  Absorbed stellar:     %.2e W/m2\n", absorbed_stellar);
        printf("  Internal heat:        %.2e W/m2\n", internal_flux);
        printf("  TOA net outgoing:     %.2e W/m2\n", radiationO);
        printf("  Energy imbalance:     %.2e W/m2\n", energy_balance);
        printf("=======================================================\n");
    }
    
    // Return radiation flux values for diagnostics
    *radiationI0_out = radiationI0;
    *radiationI1_out = radiationI1;
    *radiationO_out = radiationO;
	
//atexit(pexit);exit(0); //ms debugging mode
}
