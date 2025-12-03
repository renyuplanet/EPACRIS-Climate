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


#include <math.h>
#include "constant.h"
//#include "ms_2stream.c" //ms: solving for radiative fluxes using Heng+2018
#include "ms_d2stream.c" //ms: solving for radiative fluxes using Malik+2019; Deitrick+2021
#include "toon_2stream.c" //ms: solving for radiative fluxes using Toon+89

void ms_RadFlux(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[]);

void ms_RadFlux(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[])
{

    int i=0, j, k, jj, l, j1, kk;
    
    double GA;
    GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET;
	double mole2dust;
	mole2dust = PI*pow(AERSIZE,3)*AERDEN/6.0/AMU*2.0558;
	
    
    /* Temperature variation */
//ms: for flux jacobian to solve for steady-state temperatures
    int nrl, jTvar;
    nrl = zbin-ncl;
    double **Tvar,Told[zbin+1];
    double tempvar=RJACOB_TEMPVAR;
    double tempvarc;
    int isconv1[zbin+1],radcount;
    double lapse1[zbin+1];
    for (k=1; k<=zbin; k++) {
        isconv1[k] = isconv[zbin+1-k]; /* reverse, from top to bottom */
        lapse1[k] = lapse[zbin+1-k];
    }
    radcount = 0;

    jTvar = nrl+1;
    if (TIME_STEPPING) jTvar = 0;
    Tvar = dmatrix(0,zbin,0,jTvar);
    for (j=0; j<=zbin; j++) {
        Tvar[j][0] = T[zbin-j];
        Told[j] = T[zbin-j]; //conserve for implicit loop
    }
    if (!TIME_STEPPING)//otherwise Tvar is effectively a vector
    {
//ORIGINAL:
 // This changes single boundary layer temperature, and compares it's effect on boundary fluxes
 // Not really physically correct, to change one single boundary. (see Test approach below)
    for (j=0; j<=zbin; j++) {
        Tvar[j][1] = T[zbin-j];
    }
    Tvar[0][1] = Tvar[0][1]*(1.0+tempvar);
    
//    Tvar[0][1] = Tvar[0][1]*(1+0.01*tempvar); //ms22: testing

    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { // radiative layer, perturb temperature
            radcount = radcount+1;
            for (j=0; j<=zbin; j++) {
                Tvar[j][radcount+1] = T[zbin-j]; //ms2021: waste of computation within k loop
            }
            Tvar[k][radcount+1] = Tvar[k][radcount+1]*(1.0+tempvar);
//            Tvar[k][radcount+1] = Tvar[k][radcount+1]*(1+0.01*tempvar);
            kk=1; //ms: longer than just writing +1 later?
            tempvarc=tempvar*Tvar[k][radcount+1]/(1.0+tempvar);
            // printf("%s %e %f\n", "perturb p=", P[zbin-k], tempvarc);
            while (k+kk<=zbin && isconv1[k+kk] == 1) { // also perturb all subsequent convective layer
                tempvarc  = tempvarc*pow(P[zbin-k-kk]/P[zbin-k-kk+1],lapse1[k+kk]);
                // printf("%s %e %f\n", "convective perturb p=", P[zbin-k-kk], tempvarc); 
                Tvar[k+kk][radcount+1] = Tvar[k+kk][radcount+1]+tempvarc;
                kk = kk+1;
            }
        }
    }
    //for (j=0; j<=zbin; j++) {
    //    for (k=0; k<=nrl+1; k++) {
    //        printf("%f\t",Tvar[j][k]);
    //    }
    //    printf("\n");
    //}
//END ORIGINAL 

//MS2022: NOW let's test a different approach to Temperature Variation: 
// If a layer temperature is changed, then there is a delta (dT) and pairwise different boundary tempertures Tvar
// Let's now implement this thought and see if this is more stable than the (numerical?) oscillations in the deep 
// atmosphere when solving with jacobian.
// Also while implementing this, the direction of convective treatment is changed to bottom up for less change in internal energy between time steps.

/*    for (j=zbin; j>=0; j--){ for(k=jTvar;k>=0;k--){ Tvar[j][k] = T[zbin-j];}} //fill the matrix first
    for (k=zbin;k>=1;k--){ 
        if(isconv1[k]==0) radcount=radcount+1; //if radiative layer
        //the following represents a dT in layer k of value tempvar:
        Tvar[k][jTvar+1-radcount] = Tvar[k][jTvar+1-radcount] + tempvar/2.0; //lower boundary
        Tvar[k-1][jTvar+1-radcount] = Tvar[k-1][jTvar+1-radcount] + tempvar/2.0; //upper boundary
        kk=1; //cycle through convective layers above rad layers
        tempvarc=tempvar/2.0; //for convective layers
        while (k-kk>=0 && isconv1[k-kk] ==1) //cycle through convective layers above rad layers
        {
            tempvarc = tempvarc * pow(P[zbin-k-kk+1]/P[zbin-k-kk],lapse1[k-kk]);
            Tvar[k-kk][jTvar+1-radcount] = Tvar[k-kk][jTvar+1-radcount] + tempvarc;
            kk = kk+1;
        }//e.o.while loop
    }//e.o. "k" loop
*/
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
        double sumBint=0.0; // check against sigma*tint^4
    
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
		//temperature[j] = tl[j]; //tl not updated before RT iterations have finished
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
        dt[j] = 1.0e-1; //for TIME_STEPPING: init in some random unit to ...
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
    
    //Up/Dn fluxes for diagnostics
    double Fup[zbin+1],Fdn[zbin+1];
    
    double radiationI0, radiationI1, radiationO;
    radiationI0=0;
    radiationI1=0;
    radiationO=0;
    
//    FILE *fp, *fopa;
//    fp=fopen("AuxillaryOut/checkemission.dat","w");
//    fopa=fopen("AuxillaryOut/checkopacity.dat","w");


//========================================================    
//========================================================    
//----------------------------------------------------                
// Here starts the implicit "time-stepping" loop 
//----------------------------------------------------
//First time it calculates initial fluxes and temperature changes
//Subsequent times it re-runs RT with suggested temperature update and compares fluxes to initial state
//If fluxes increased, then lower dTemp 
//This loop breaks if fluxes gotten better and exit to ms_rad_conv.c
//========================================================    
//========================================================    

    double sumflux, sumflux_init=0.0; 
    int iTvar = jTvar;//subsequent RT calls require Tvar vector only, but allocated memory for arrays has to be fully freed later
    int Tbetter=0;//once fluxes gotten better, Tbetter = 1 and test loop exits
    int jtest=0; //count test loops
    double *deltaT,dtmul;
    deltaT=dvector(1,nrl+1);

    while (Tbetter==0)
    {
        jtest += 1;
        printf("%s %d\t %s %d\n","jtest = ",jtest,"Tbetter = ",Tbetter); 
        sumflux=0.0; //re-calc every step
        if (jtest > 1) jTvar = 0;

	for (j=0; j<=nrl; j++) {for (k=0; k<=jTvar; k++) {NetFlux[j][k] = 0.0;}}//reset fluxes
	for (j=0; j<=zbin; j++) {Fup[j] = 0.0;Fdn[j] = 0.0;}
    

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
            //wa[j] += opacNH3[j1][i]*xx[j1][9];      //
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
		
                //printf("%s %s %s %s %s\n","Lambda", "tau", "TAUdoub[center]", "TAUdoub[bound]", "TAU");
		/* Optical Depth of Each Layer */
		for (j=1; j <=zbin; j++) {
                    j1   = zbin+1-j;
                    tau[j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-P[j1])*1.0E-4;
                    TAUdoub[2*j-1] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(pl[j1]-P[j1])*1.0E-4; //ms2023: double grid
                    TAUdoub[2*j] = (wa[j] + ws[j])/MM[j1]/AMU/meanmolecular[j1]/GA*(P[j1-1]-pl[j1])*1.0E-4; //ms2023: double grid
	            //printf("%f %e %f %f\n", wavelength[i], tau[j], w[j], g[j]);
	            //printf("%f %e %e %e %e\n", wavelength[i], tau[j], TAUdoub[2*j-1], TAUdoub[2*j], TAUdoub[2*j-1]+TAUdoub[2*j]);
		}
//atexit(pexit);exit(0); //ms debugging mode
/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 *if any tau[j] >> 1.0, radiative fluxes may become diffusion limited
 *in such cases a control run with modified grid may be performed.
 *
 *At this point it is worth checking also for Photon deposition depth (region absorbing a majority of stellar radiation)
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/

//========================================================    
//====== 2-STREAM switch =================================	
    if (TWO_STR_SOLVER == 1 ) ms_two_str_solver( i, w, g, tau, nrl, isconv, Tvar, jTvar, P, NetFlux,Fup, Fdn, Tint, sumBint);
    else toon_two_str_solver( i, w, g, tau, nrl, isconv, Tvar, jTvar, P, NetFlux,Fup, Fdn, Tint, sumBint);
//========================================================    
//====== END of 2-STREAM switch ==========================	
    
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
/*printf("%s\n","\n==== Fup/Fdn integrated ===="); //template
for (i=0;i<=zbin;i++)
{
    printf("%s %d\t%.3e\t%.3e\n","Fup,Fdn ",i,Fup[i],Fdn[i]);
}
printf("%s\n","==== END FLUXES ===="); //template
*/
//atexit(pexit);exit(0); //ms debugging mode
//
//check Tint ve Bint:
//printf("%s %e %s %e %s %e\n","Energy from Interior:",SIGMA*pow(Tint,4.0)," SUM(pi*Bint):",PI*sumBint, " DELTA:",SIGMA*pow(Tint,4.0) - PI*sumBint);
//    printf("%s %e\n", "For comparison, NetFlux[Surf] =",NetFlux[nrl][0]);

    double drfluxmax, drflux[zbin+1];    
    drfluxmax = 0.0;
    for (j=1; j<zbin; j++) // exclude BOA for stability
    {
        drflux[j] = Fup[j-1]-Fdn[j-1]-Fup[j]+Fdn[j];
        drfluxmax = fmax(drfluxmax,fabs(drflux[j]) );
    }
    printf("%s %f\n","dRFLUX_max= ",drfluxmax);
    /*printf("%s\n","\n==== dRFLUX ===="); //template
    for (i=0;i<=zbin;i++)
    {
    //printf("%s %d\t%.3e\n","Netflux ",i,Fup[i]-Fdn[i]);
    if(i>0) printf("%s %d\t%.3e\n","dNetflux ",i,Fup[i]-Fdn[i]-Fup[i-1]+Fdn[i-1]);
    //printf("%s %d\t%.3e\n","Rflux ",i,Fup[i]-Fdn[i]-SIGMA*pow(Tint,4.0));
    }
    printf("%s\n","==== END dRFLUX ===="); //template
    */

    /* Obtain the Flux and the Flux Jacobian, from bottom to up */
    for (j=0; j<=nrl; j++) {
        if (j < zbin) { // exclude BOA for stability
            Rflux[j] = NetFlux[j][0] - SIGMA*pow(Tint,4.0); 
            if(jtest==1) sumflux_init += fabs(Rflux[j]);
            else sumflux += fabs(Rflux[j]);
        }
    }
    printf("%s %e\t %s %e\n","sumflux_init= ",sumflux_init,"sumflux= ",sumflux);
    
    printf("%s\n","\n==== RFLUX ===="); //template
    printf("%e\t%f\n",P[zbin],Rflux[0]);
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k] == 0) {
            radcount += 1;
            printf("%e\t%f\n",P[zbin-k],Rflux[radcount]);
        }
    }
    printf("%s\n","==== END RFLUX ===="); //template

//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(!TIME_STEPPING)
{//ms22tss
//=====================================
    if(jtest==1) //jacobian and dTemp vector only at first iteration
    {
//=====================================
    double **jmax;
    jmax=dmatrix(1,nrl+1,1,nrl+1);
    for (j=0; j<=nrl; j++) {
        for (kk=0; kk<jTvar; kk++) {
            jmax[j+1][kk+1] = (NetFlux[j][kk+1] - NetFlux[j][0])/tempvar; //ms22: this should stay tempvar rather than tempvar/2.0, since 1 layer was changed by tmepvar - leading to two boundaries changing with tempvar/2.0?
        }
    }
/*    
printf("%s\n","\n==== NetFlux ===="); //template
    for (j=0; j<jTvar; j++) {
        printf("j=%d: ",j);
        for (k=0; k<jTvar; k++) {
            printf("%e\t",NetFlux[j][k]);
        }
        printf("\n");
    } 
    
printf("%s\n","\n==== Jmax ===="); //template
    for (j=1; j<=jTvar; j++) {
        printf("j=%d: ",j);
        for (k=1; k<=jTvar; k++) {
            printf("%e\t",jmax[j][k]);
        }
        printf("\n");
    } 
*/    

    /* Solve Equation */
    int *indx;
    double ddd;
    double *deltaT;
    double relaxationfactor; //Renyu2021
    indx=ivector(1,nrl+1);
    deltaT=dvector(1,nrl+1);
    double dtmax=0.0; //MS-dTlayer: for stepping towards equilibrium

    for (j=0; j<=nrl; j++) {
        deltaT[j+1] = -Rflux[j];
    }

//----------------------------------------------------                
//markus: this solves for deltaT from flux jacobian: 
//----------------------------------------------------                
//ms:rfluxconv    
    ludcmp(jmax,nrl+1,indx,&ddd); 
    lubksb(jmax,nrl+1,indx,deltaT);

    free_ivector(indx,1,nrl+1);
    free_dmatrix(jmax,1,nrl+1,1,nrl+1);
    //MS-dTlayer: ok so here now deltaT should be the change for "layers"(incl. one above model grid) - not boundaries?
    //let's test this hypothesis:

//atexit(pexit);exit(0); //ms debugging mode
    for (j=1; j<=jTvar; j++) {
        printf("%s %d %f\n","deltaT ",j,deltaT[j]);
    } 

    //limit maximum change to 
    
    /* calculate the new temperature profile */
//Renyu2021:    
    radcount = 0;
    //MS-dTlayer: relaxationfactor = R_RELAX;
    if (isconv1[0]==0) dtmax = fabs(deltaT[1])/Told[0];
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer*/
            radcount = radcount+1;
            dtmax = fmax(dtmax,fabs(deltaT[radcount+1])/Told[k]); //MS-dTbound: find maximum suggested relative change in dTemp
            //MS-dtbound: relaxationfactor=fmin(relaxationfactor,1e-2*Tvar[k][0]/fabs(deltaT[radcount+1])); /* do not allow temperature to change more than 10% */
            //printf("%s\t%f\n", "Relaxation factor is", relaxationfactor);
        }
    }
    printf("%s %f\t%s %f\n","dT_max = ",dtmax," R_Relax = ",R_RELAX);
    if(dtmax>R_RELAX) //only if we exceed the maximum allowed change per step:
    {
        printf("%s\n","deltaT adjusted so dtmax < R_RELAX");
        for (j=1; j<=nrl+1; j++) {
//Renyu2021:        deltaT[j] *= R_RELAX; /* relaxation */
//ms:rfluxconv
//        deltaT[j] *= relaxationfactor; /* relaxation */
          deltaT[j] *= R_RELAX / dtmax; //MS-dTlayer: scale dT  vector againstmaximum allowed % change in a layer per step
          printf("%s %d %f\n","deltaT_adj ",j,deltaT[j]);
        }
    }

    dtmul = 1.0; //dT multiplicator
//    dtmul = 10.0; //dT multiplicator
    } //e.o.jtest==1

    else //jtest >1
    {
        if(sumflux == sumflux_init) dtmul *= 2.0; //theoretical, but : maybe change too small?
//        if(sumflux > sumflux_init) dtmul *= fmax(0.1*dtmul, fmin(pow(sumflux_init / sumflux,1.0), dtmul)); // change dt-multiplicator between 0.1 - 1.0 via (flxs/flxs_init)^2
        if(sumflux > sumflux_init) dtmul *= 0.5; //fmax(0.1, fmin(pow(sumflux_init / sumflux,1.0), 1.0)); // change dt-multiplicator between 0.1 - 1.0 via (flxs/flxs_init)^2
    }

    printf("%s %e\n","dTmul = ",dtmul);
    if (isconv1[0]==0) Tvar[0][0] = Told[0] + dtmul*deltaT[1];//MS-dTbound
    //Tvar[0][0] = Told[0] + deltaT[1]/2.0;//MS-dTlayer
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            tempvarc = dtmul*deltaT[radcount+1];
//ms: correct for temperature spikes: 
/*            if(fabs(Tvar[k][0]+tempvarc-Tvar[k-1][0]) > ) 
            {
                tempvarc = ;
            }*/
//ms: end correction
            //if (k<zbin) Tvar[k][0] = Tvar[k][0] + tempvarc; //ms22 treat surface separately below
            //if (k == zbin && isconv1[k] == 0) Tvar[k][0] = Tvar[k-1][0]; //ms22 for deep atmospheres BOA temp can be set to last atmo layer temp before convection
            Tvar[k][0] = Told[k] + tempvarc;
//end
            /* printf("%s %e %f\n", "update p=", P[zbin-k], tempvarc); */
            kk=1;
            while (k+kk<=zbin && isconv1[k+kk] == 1) { /* also perturb all subsequent convective layer */
                tempvarc  = tempvarc*pow(P[zbin-k-kk]/P[zbin-k-kk+1],lapse1[k+kk]);
                Tvar[k+kk][0] = Told[k+kk] + tempvarc;
                /* printf("%s %e %f\n", "convective update p=", P[zbin-k-kk], tempvarc); */
                kk = kk+1;
            }
        }
    }

//=====================================
//=====================================
} //eo ms22tss: no time stepping
//=====================================
//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(TIME_STEPPING)
{//ms22tss
//This is experimental at the moment (August 2022). Need to find optimal dt or dTemp ratio from TOA to BOA so that
// (A) upper atmosphere converges quickly enough
// (B) lower atmosphere converges quickly enough
// (C) no spurious oscillations occur numerically
// (D) convergence of whole atmosphere occurs in a reasonable runtime frame
// (E) it works for thin terrestrial (e.g. Mars-like) as well as deep gas giant (e.g. Jupiter-like) atmospheres
// (F) a p[BOA] or T[BOA] based (automatic)switch between (E) cases is conceivable

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

    double  dTemp, f0, f1, f2;
    if (rtstepcount==1) rt_drfluxmax_init=drfluxmax;//ms: let's scale things as a ratio to initial fluxes

//    for (j=1; j<zbin; j++)//ms22 treat surface separately below
    for (j=1; j<zbin; j++)//surface treatment included

    {
        dt[j] = cp[zbin+1-j]/GA/SIGMA/pow(Told[j],3.0)*P[zbin-j];
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
        Tvar[j-1][0] = Told[j-1] - f0 * fmin(fabs(dt[j]* f1 /pow(P[zbin-j]-P[zbin+1-j],1.0 - f2 )),1E-2*Told[j-1]);
        Tvar[j][0] = Told[j] - f0 * fmin(fabs(dt[j]* f1 /pow(P[zbin-j]-P[zbin+1-j],1.0- f2 )),1E-2*Told[j]);
//-----------------------------------------------------------
        
//testing HELIOS implementation:  
//-----------------------------------------------------------
        //dt[j] = 1.0; //should be 0.5 * 1e5
        //Tvar[j-1][0] = Tvar[j-1][0] + pow(drflux[j],0.1)/pow(Tvar[j-1][0],3.0)*dt[j]/SIGMA/log(P[zbin+1-j]/P[zbin-j]);
        //Tvar[j][0] = Tvar[j][0] + pow(drflux[j],0.1)/pow(Tvar[j][0],3.0)*dt[j]/SIGMA/log(P[zbin+1-j]/P[zbin-j]) ;
       
    }
//    Tvar[zbin][0] = Tvar[zbin-1][0];    //ms22 for deep atmospheres BOA temp can be set to last atmo layer temp before convection
}//eo ms22tss

if(jtest >1 && sumflux < sumflux_init) Tbetter = 1; //END the test loop
if(jtest>100){ atexit(pexit);exit(0);} //ms debugging mode
//=====================================
if(jtest=1) Tbetter = 1; //ms quick fix for now
//=====================================
//-----------------------------------------------------------
} //e.o. while loop: implicit time-stepping 
//-----------------------------------------------------------
//=====================================
//=====================================
    printf("%s\t%d\n", "test steps needed: ", jtest-1);
    
    for (j=0; j<=zbin; j++) {
        tempbnew[j] = fmin(10000.0,fmax(Tvar[zbin-j][0],1.0)); //ms22: also include max temp
    }
    
    free_dmatrix(Tvar,0,zbin,0,iTvar);
    free_dmatrix(NetFlux,0,nrl,0,iTvar);
    free_dvector(deltaT,1,nrl+1);
    
    /*printf("%s\t%f\n", "Top-of-Atmosphere incoming radiation flux is", radiationI0);
    printf("%s\t%f\n", "Bottom-of-Atmospehre incoming radiation flux is", radiationI1);*/
    printf("%s\t%f\n", "TOA outgoing net radiation flux is", radiationO);
	
//atexit(pexit);exit(0); //ms debugging mode
}
