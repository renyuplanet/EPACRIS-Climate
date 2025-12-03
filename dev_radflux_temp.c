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
#include "ms_d2stream.c" //ms: solving for radiative fluxes using Heng+2018; Deitrick+2021
#include "toon_2stream.c" //ms: solving for radiative fluxes using Toon+89

void ms_RadFlux(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[], double dt[], int isequil[]);

void ms_RadFlux(double Rflux[], double tempbnew[], double P[], int ncl, int isconv[], double lapse[], double T[], double Tint, double cp[], double dt[], int isequil[])
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
        lapse1[k] = lapse[zbin+1-k];
    }
    radcount = 0;

    jTvar = nrl+1;
    if (TIME_STEPPING) jTvar = 0;
    Tvar = dmatrix(0,zbin,0,jTvar);
    for (j=0; j<=zbin; j++) {
        Tvar[j][0] = T[zbin-j];
    }
    if (!TIME_STEPPING)//otherwise Tvar is effectively a vector
    {
    for (j=0; j<=zbin; j++) {
        Tvar[j][1] = T[zbin-j];
    }
    Tvar[0][1] = Tvar[0][1]+tempvar;
    
//    Tvar[0][1] = Tvar[0][1]*(1+0.01*tempvar); //ms22: testing

    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            for (j=0; j<=zbin; j++) {
                Tvar[j][radcount+1] = T[zbin-j]; //ms2021: waste of computation within k loop
            }
            Tvar[k][radcount+1] = Tvar[k][radcount+1]+tempvar;
//            Tvar[k][radcount+1] = Tvar[k][radcount+1]*(1+0.01*tempvar);
            kk=1; //ms: longer than just writing +1 later?
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
    NetFlux = dmatrix(0,nrl,0,jTvar);
    
	for (j=0; j<=nrl; j++) {
        for (k=0; k<=jTvar; k++) {
            NetFlux[j][k] = 0.0;
        }
	}
    
    //Up/Dn fluxes for diagnostics
    double Fup[zbin+1],Fdn[zbin+1],Fcup[zbin+1],Fcdn[zbin+1];
	for (j=0; j<=zbin; j++) 
        {
            Fup[j] = 0.0;
            Fdn[j] = 0.0;
            Fcup[j] = 0.0; //zero element unused
            Fcdn[j] = 0.0; //zero element unused
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
            
            //ws[j] += crossr[i]*MM[j1];
            
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
    if (TWO_STR_SOLVER == 1 ) ms_two_str_solver( i, w, g, tau, nrl, isconv, Tvar,jTvar, P, NetFlux,Fup, Fdn,Fcup,Fcdn, Tint, sumBint);
    if (TWO_STR_SOLVER == 0 ) toon_two_str_solver( i, w, g, tau, nrl, isconv, Tvar, jTvar, P, NetFlux,Fup, Fdn,Fcup,Fcdn, Tint, sumBint);
    
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
*/
/*printf("%s\n","\n==== TVAR ===="); //template
for (j=0;j<=nrl;j++)
{   
    printf("%d\t",j);
    printf("%.2e\t",Tvar[j][0]);
    printf("\n");
}
printf("%s\n","\n==== END TVAR ===="); //template
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
    for (j=0; j<=nrl; j++) {
        Rflux[j] = NetFlux[j][0] - SIGMA*pow(Tint,4.0); 
    }
    
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
//=====================================
//ms22: testing a TIME_STEPPING switch: (commented as ms22tss)
if(!TIME_STEPPING)
{//ms22tss
//=====================================
//=====================================

    double **jmax;
    //nosurf jmax=dmatrix(1,nrl+1,1,nrl+1);
    jmax=dmatrix(1,nrl,1,nrl);
    //nosurf for (j=0; j<=nrl; j++) {
    for (j=0; j<nrl; j++) {//excl.surf
        //nosurf for (kk=0; kk<jTvar; kk++) {
        for (kk=0; kk<jTvar-1; kk++) {
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
    double relaxfvector[nrl+2]; //so I don't have to define as dvector, and free later
    //nosurf indx=ivector(1,nrl+1);
    indx=ivector(1,nrl);
    //nosurf deltaT=dvector(1,nrl+1);
    deltaT=dvector(1,nrl);

    //for (j=0; j<=nrl; j++) {
    for (j=0; j<nrl; j++) {
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
    //nosurf ludcmp(jmax,nrl+1,indx,&ddd); 
    //nosurf lubksb(jmax,nrl+1,indx,deltaT);
    ludcmp(jmax,nrl,indx,&ddd); 
    lubksb(jmax,nrl,indx,deltaT);
    
//atexit(pexit);exit(0); //ms debugging mode
    //nosurf for (j=0; j<=nrl; j++) {
    for (j=0; j<nrl; j++) {
        printf("%s %d %f\n","deltaT ",j+1,deltaT[j+1]);
    } 
    
    /* calculate the new temperature profile */
//Renyu2021:    
    radcount = 0;
    relaxationfactor = R_RELAX;
    relaxationfactor=fmin(relaxationfactor,1e-2*Tvar[0][0]/fabs(deltaT[1])); //ms2023: include TOA
    //relaxationfactor=fmin(relaxationfactor,1e-2*Tvar[0][0]/fabs(deltaT[1]*Rflux[1])); //ms2023: include TOA
    relaxfvector[1]=fmin(R_RELAX,1e-2*Tvar[0][0]/fabs(deltaT[1]*Rflux[1])); /* do not allow temperature to change more than 10% */
    //nosurf for (k=1; k<=zbin; k++) { 
    for (k=1; k<zbin; k++) { 
        if (isconv1[k]==0) { /* radiative layer*/
            radcount = radcount+1;
            relaxationfactor=fmin(relaxationfactor,1e-2*Tvar[k][0]/fabs(deltaT[radcount+1])); /* do not allow temperature to change more than 10% */
            //relaxationfactor=fmin(relaxationfactor,1e-2*Tvar[k][0]/fabs(deltaT[radcount+1]*Rflux[radcount+1])); /* do not allow temperature to change more than 10% */
            relaxfvector[radcount+1]=fmin(R_RELAX,1e-2*Tvar[k][0]/fabs(deltaT[radcount+1]*Rflux[radcount+1])); /* do not allow temperature to change more than 10% */
            //printf("%s\t%f\n", "Relaxation factor is", relaxationfactor);
        }
    }
    //nosurf for (j=1; j<=nrl+1; j++) {
    for (j=1; j<nrl+1; j++) {
//Renyu2021:        deltaT[j] *= R_RELAX; /* relaxation */
//ms:rfluxconv
        deltaT[j] *= relaxationfactor; /* relaxation */
//        deltaT[j] *= relaxfvector[j]; /* relaxation */
    }
    Tvar[0][0] = Tvar[0][0] + deltaT[1];
    radcount = 0;
    for (k=1; k<=zbin; k++) {
        if (isconv1[k]==0) { /* radiative layer, perturb temperature */
            radcount = radcount+1;
            //nosurf tempvarc = deltaT[radcount+1];
            if (k<zbin) tempvarc = deltaT[radcount+1];
//ms: correct for temperature spikes: 
/*            if(fabs(Tvar[k][0]+tempvarc-Tvar[k-1][0]) > ) 
            {
                tempvarc = ;
            }*/
//ms: end correction
            if (k<zbin) Tvar[k][0] = Tvar[k][0] + tempvarc; //ms22 treat surface separately below
            if (k == zbin && isconv1[k] == 0) Tvar[k][0] = Tvar[k-1][0]; //ms22 for deep atmospheres BOA temp can be set to last atmo layer temp before convection
            //Tvar[k][0] = Tvar[k][0] + tempvarc; //ms2023
//end
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

    //nosurf free_dmatrix(jmax,1,nrl+1,1,nrl+1);
    //nosurf free_dvector(deltaT,1,nrl+1);
    //nosurf free_ivector(indx,1,nrl+1);
    free_dmatrix(jmax,1,nrl,1,nrl);
    free_dvector(deltaT,1,nrl);
    free_ivector(indx,1,nrl);
//=====================================
//=====================================
} //eo ms22tss: no time stepping
    free_dmatrix(NetFlux,0,nrl,0,jTvar);
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

    double drfluxmax, drflux[zbin+2], dTemp, f0, f1, f2, deltaT[zbin+2], fratio, change;
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

//TS_SCHEME switch
if (TS_SCHEME == 0){ 
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

}// eo TS_SCHEME=0
if (TS_SCHEME == 1){
    for (j=1; j<=zbin; j++)
    {
        drflux[j] = Fup[j-1]-Fdn[j-1]-Fup[j]+Fdn[j];
        drfluxmax = fmax(drfluxmax,fabs(drflux[j]) );
        printf("%s %d %f\n","drflux ",j,drflux[j]);
    }
    drflux[zbin+1] = Fup[zbin]-Fdn[zbin] - SIGMA*pow(Tint,4.0) ; //Surface "layer" separately
    drfluxmax = fmax(drfluxmax,fabs(drflux[zbin+1]));
    printf("%s %f\n","dRFLUX_max= ",drfluxmax);
    if (rtstepcount==1) rt_drfluxmax_init=drfluxmax;//ms: let's store in case we need it

    for (j=1; j<=zbin+1; j++)
    {
        deltaT[j] = - dt[j] * pl[1]/(P[0]-P[1]) * drflux[j] / pow(fabs(drflux[j]),0.9); //HELIOS code
        //check if |drflux|/sigma*T < epsval -> layer in equilibrium:
        fratio = fabs(drflux[j]) / SIGMA / pow(Tvar[j-1][0],4.0);
        if (fratio <= Tol_FRATIO) {isequil[zbin+1-j] = 1;} else {isequil[zbin+1-j] = 0;}
        printf("%s %d %f \t%s %.2e\n","deltaT ",j,deltaT[j], "dF/sigma*T^4 ", fratio);
    }
//atexit(pexit);exit(0); //ms debugging mode
    //add smoothing:
/*    double smoothexp=4.0;
    for (j=2; j<=zbin-1; j++)
    {
        deltaT[j] += pow( - (Tvar[j-2][0]-Tvar[j-1][0]-Tvar[j][0]+Tvar[j+1][0])/4.0 ,smoothexp);
        printf("%s %d %f\n","deltaT smoothed",j,deltaT[j]);
    }
    deltaT[1] += pow(0.25 * (Tvar[2][0] - Tvar[1][0]),smoothexp); //TOA: above isothermal
    deltaT[zbin] += pow(0.25 * (Tvar[zbin-2][0] - Tvar[zbin-1][0]),smoothexp); //BOA: below isothermal
    deltaT[zbin] += pow(0.5 * (Tvar[zbin-1][0] - Tvar[zbin][0]),smoothexp); //Surf: below isothermal
    printf("%s %f\n","deltaT smoothed 1",deltaT[1]);
    printf("%s %f\n","deltaT smoothed zbin",deltaT[zbin]);
    printf("%s %f\n","deltaT smoothed zbin+1",deltaT[zbin+1]);
*/
/*    Tvar[zbin][0] = Tvar[zbin][0] + deltaT[zbin+1];
    for (j=zbin-1; j>0; j--)
    {
        Tvar[j][0] = Tvar[j][0] + (deltaT[j] + deltaT[j+1]) / 2.0;
    }
    Tvar[0][0] = Tvar[0][0] + deltaT[1] + (deltaT[1]-deltaT[2]) / 2.0;
*/
    Tvar[zbin][0] = Tvar[zbin][0] + deltaT[zbin+1];
    for (j=zbin-1; j>0; j--)
    {
        Tvar[j][0] = Tvar[j][0] + (deltaT[j] + deltaT[j+1]) / 2.0;
    }
    Tvar[0][0] = Tvar[0][0] + deltaT[1] + (deltaT[1]-deltaT[2]) / 2.0;

    //test Renyu's version:
/*    change = deltaT[zbin+1];
    for (j=zbin; j>=0; j--)
    {
        Tvar[j][0] = Tvar[j][0] + change; 
        change = 2.0 * ( deltaT[j] - 0.5*change );
    }
*/    
    //test reverse (TOA-BOA) version:
/*    change = deltaT[1];
    for (j=0; j<=zbin; j++)
    {
        Tvar[j][0] = Tvar[j][0] + change; 
        change = 2.0 * ( deltaT[j+1] - 0.5*change );
    }
*/
}//eo TS_SCHEME=1

}//eo ms22tss

    
    for (j=0; j<=zbin; j++) {
        tempbnew[j] = fmin(10000.0,fmax(Tvar[zbin-j][0],1.0)); //ms22: also include max temp
    }
    
    free_dmatrix(Tvar,0,zbin,0,jTvar);
    
    /*printf("%s\t%f\n", "Top-of-Atmosphere incoming radiation flux is", radiationI0);
    printf("%s\t%f\n", "Bottom-of-Atmospehre incoming radiation flux is", radiationI1);*/
    printf("%s\t%f\n", "TOA outgoing net radiation flux is", radiationO);
	
//atexit(pexit);exit(0); //ms debugging mode
}
