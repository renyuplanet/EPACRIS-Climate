/*Function to provide solar radiation at a given height*/
/*
 * Original Author: Renyu Hu (renyu.hu@jpl.nasa.gov)
 * Version: 2021, 2-Stream update to Heng+2019 flux equations including direct stellar beam
 * Editor v2021: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Modifier symbols: 
    * 'ms_' ... in file names
    * 'markus<date>' ... in description lines/blocks 
    * 'ms:' ... in-line comments
*/

#include <math.h>
#include "constant.h"
#include "ludcmp.c"
#include "lubksb.c"
#include "ms_ludcmp.c"
#include "ms_lusolve.c"
#include "BTridiagonal.c"
//#include "ms_radflux_temp.c"
#include "ms_dradflux_temp.c"
//#include "ms_radflux_test.c"
#include "ms_conv_funcs.c"

void ms_RadConv(double tempeq[], double P[], double T[], double Tint, char outnewtemp[]);

void ms_RadConv(double tempeq[], double P[], double T[], double Tint, char outnewtemp[])
{
    int i,j,k,radcount;
    
    /* Initial profile */
    double tempb[zbin+1], tempbnew[zbin+1], tempc[zbin+1];
    for (j=0; j<=zbin; j++) {
		tempb[j] = T[j]; //ms23: double grid
		//tempb[j] = Tdoub[2*j]; //bottom-up!
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
    
    /* Define variables for iteration */
    double znew[zbin+1], scaleheight, GA;
    GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET;
    double lapse[zbin+1], gasheat, cp[zbin+1];
    double Rflux[zbin+1];
    double Rfluxmax,dRfluxmax;
    double drylapse[zbin+1];
    
    int isconv[zbin+1], ncl, nrl; /* if a layer is convective */
    int isconv_old[zbin+1], prevconv[zbin+1];
    int isconv_sum, deltaconv;
    double t_old[zbin+1],tl_old[zbin+1]; //for convective adj test
    
    for (j=0; j<=zbin; j++) {
        isconv[j] = 0;
    }
    
    /* gas heat capacity and laspe rate */
    for (j=1; j<=zbin; j++) {
        gasheat=(CO2Heat(tl[j])*xx[j][52]+HeHeat(tl[j])*heliumnumber[j]+N2Heat(tl[j])*xx[j][55]+NH3Heat(tl[j])*xx[j][9]+CH4Heat(tl[j])*xx[j][21]+H2Heat(tl[j])*xx[j][53]+O2Heat(tl[j])*xx[j][54]+COHeat(tl[j])*xx[j][20]+H2OHeat(tl[j])*xx[j][7])/(xx[j][52]+heliumnumber[j]+xx[j][55]+xx[j][9]+xx[j][21]+xx[j][53]+xx[j][54]+xx[j][20]+xx[j][7])*1000.0/meanmolecular[j];
        //lapse[j] = KBOLTZMANN/meanmolecular[j]/AMU/gasheat ; /* dimensionless dlog(T)/dlog(P) */
        //printf("%s %d %s %e\n", "DRY ADIABAT as kb/MM/AMU/cp/1000*MM", j, " = ", lapse[j]);
        cp[j] = gasheat*meanmolecular[j]/1000.0; //ms22: for convective adjustment update
        lapse[j] = R_GAS / cp[j]; //for numerical stability - avoiding div by e-27
        //printf("%s %d %s %e\n", "DRY ADIABAT: R/cp ", j, " = ", lapse[j]);
        drylapse[j] = lapse[j]; //testing
//comment out for testing        ms_adiabat(j,lapse,heliumnumber[j]); //calculate appropriate lapse rate depending on dry, moist, and condensate conditions
        //printf("%s %d %s %e\n", "MOIST ADIABAT ", j, " = ", lapse[j]);
        //printf("%s %d %s %e\n", "Dry - Moist ", j, " = ", drylapse[j] - lapse[j]);
        //printf("%s\n\n",fillmi);
    }

//atexit(pexit);exit(0); //ms debugging mode
//*************************************************************
// main iterative loop:
//*************************************************************
    for (i=1; i<=NMAX_RC; i++) {
        
        printf("%s %d\n", "Rad-Conv loop i = ", i);
        /* store last RC boudary */
        for (j=0; j<=zbin; j++) t_old[j] = tempb[j];
        for (j=1; j<=zbin; j++) {
            //tl_old[j] = tl[j]; //ms23: seems weird to overwrite here
            isconv_old[j] = isconv[j];
            prevconv[j] = isconv[j]; //for convective adjustment iteration
            isconv[j] = 0;
        }

        if(CONVEC_ADJUST==2)
        //if(CONVEC_ADJUST==2 && i!=1) //sm22: always start all radiative to avoid input influenced solutions
        {
            ncl = 0;
            nrl = 0;
            double isconvTB[zbin+1], isconvBT[zbin+1]; //check for convection Top->Bottom and reverse
            double enth, denth; //delta enthalpy
            double tempTB[zbin+1], tempBT[zbin+1];
            for (j=0; j<=zbin; j++) tempTB[j] = tempb[j];
            for (j=0; j<=zbin; j++) tempBT[j] = tempb[j];

            //first TOA -> BOA:
            for (j=zbin; j>=1; j--) {
                if ( tempTB[j-1] > (tempTB[j] * pow(P[j-1]/P[j],lapse[j]) )-1e-2 ) {
                    tempTB[j-1] = (tempTB[j] * pow(P[j-1]/P[j],lapse[j]) );
                    if(i!=1) isconvTB[j] = 1;
                }
            }
            for (j=0; j<zbin; j++) {
                if ( tempBT[j+1] < (tempBT[j] * pow(P[j+1]/P[j],lapse[j+1]) )+1e-2 ) {
                    tempBT[j+1] = (tempBT[j] * pow(P[j+1]/P[j],lapse[j+1]) );
                    if(i!=1) isconvBT[j+1] = 1;
                }
            }
            printf("%s\t%s\t%s\t%s\n","Layer","isconvTB","isconvBT","tempb");
            for (j=zbin; j>=0; j--) {printf("%d\t%f\t%f\t%f\n",j,isconvTB[j],isconvBT[j], tempb[j]);}
            
        }
        //break;
        if(CONVEC_ADJUST==1 || CONVEC_ADJUST==0)
        //if((CONVEC_ADJUST==1 || CONVEC_ADJUST==0) && i!=1) //sm22: always start all radiative to avoid input influenced solutions
        {
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//ms22: convective adjustment iteration:
//needed to appropriately treat upper/lower tails of convective regions 
//to avoid temperature jumps where radiative meet convective regions:
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        deltaconv = zbin; //for convective adjustment iteration
        while (deltaconv > 0) 
        {
            /* new RC boundary */
            ncl = 0;
            nrl = 0;
            /* determine convection, and record convective layers */

            if(CONVEC_ADJUST==0)
            {
                for (j=zbin; j>=1; j--) {
                    //if ( tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.999999) ) {
                    if ( i!=1 && (tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.98)) ) {
                    //if ( tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.95) ) {
                        tempb[j-1] = tempb[j] * pow(P[j-1]/P[j],lapse[j]); //ms2022: moved down to include manually modded layers
                        isconv[j] = 1;
                        ncl = ncl+1;
                    } else {
                        isconv[j] = 0;
                        nrl = nrl+1;
                    }
                }
            } //END: CONVEC_ADJUST=0

            if(CONVEC_ADJUST==1) ms_conv_check(tempb, P, lapse, isconv, &ncl, &nrl); //ms22: check convective layers

//          printf("%s %d %s %d\n", "ncl", ncl, "nrl", nrl);

        // remove single radiative layer between convective layers
/*          for (j=2; j<zbin; j++) { //rh->ms 2021
                if (isconv[j]==0 && (isconv[j+1]==1 && isconv[j-1]==1)) {
                    isconv[j]=1;
                    ncl = ncl+1;
                    nrl = nrl-1;
                }
            }
*/
        
        // remove single convective layer between radiative layers
/*          for (j=2; j<zbin; j++) { //rh->ms 2021
                if (isconv[j]==1 && (isconv[j+1]==0 && isconv[j-1]==0)) {
                    isconv[j]=0;
                    ncl = ncl-1;
                    nrl = nrl+1;
                }
            }
*/        
        /* force all radiative layers initially */
/*            if (i==1) {//rh->ms 2021
                for (j=1; j<=zbin; j++) {
                    isconv[j] = 0; //ms2021: duplicate from above
                }
                ncl = 0;
                nrl = zbin;
            } 
*/        
            if(CONVEC_ADJUST==0 && i!=1)
            {
// if no sinlge layer stitching above is performed, this is done already above                
//                for (j=zbin; j>=1; j--) {
//                    if (isconv[j]) tempb[j-1] = tempb[j] * pow(P[j-1]/P[j],lapse[j]);
//                }
                printf("%s\t%s\t\t%s\t\t%s\t%s\t%s\n","Layer","cp","lapse","isconv","T_old","T_new");
                for (j=zbin; j>=0; j--) {printf("%d\t%f\t%f\t%d\t%f\t%f\n",j,cp[j],lapse[j],isconv[j], t_old[j], tempb[j]);}
//atexit(pexit);exit(0); //ms debugging mode
            }

            if(CONVEC_ADJUST==1 && ncl>=0) ms_temp_adj(tempb, P, lapse, isconv, cp); //ms22: assign convective regimes, adjust temperatures

            //ms22: check if convective regions were fully addressed:
            deltaconv = 0;
            for (j=1;j<=zbin;j++)
            {
                if (abs(prevconv[j]-isconv[j])>0) deltaconv += 1;
                prevconv[j] = isconv[j];
            }
            printf("%s %d\n","deltaconv=",deltaconv);
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        } //ms22: END of convective adjustment iteration
        } //ms22: end convec_adjust 1 or 2 iteration
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        /* hydrostatic equilibrium */
        znew[0]=0.0;
        for (j=1; j<=zbin; j++) {
            tempc[j] = (tempb[j]+tempb[j-1])/2.0;
            scaleheight = KBOLTZMANN * tempc[j] / meanmolecular[j] / AMU / GA / 1000 ;
            znew[j] = znew[j-1] - scaleheight*log(P[j]/P[j-1]);
            //ms23: double grid
            Tdoub[2*j] = tempb[j];
            Tdoub[2*j-1] = tempc[j];
        }
        Tdoub[0] = tempb[0];
        /* print out new RC boundary */
/*        for (j=1; j<=zbin; j++) {
            printf("%d\n",isconv[zbin+1-j]);
        } */
        printf("%s %d %s %d\n", "conv layer", ncl, "rad layer", nrl);
        /* compare old and new RC boundary */
        isconv_sum=0;
        for (j=1; j<=zbin; j++) {
            isconv_sum += abs(isconv_old[j]-isconv[j]);
        }
        printf("%s %d\n", "change in convective layer is", isconv_sum);
        /* determine whether converged */
        if (isconv_sum == 0 && i!=1) {
            printf("Climate converged\n");
            //ms22: smoothing temperature profile
//            tempeq[0] = tempb[0]/2. + tempb[1]/2.; //surf
//          tempeq[zbin] = tempb[zbin-1]/2. + tempb[zbin]/2.; //TOA
            for (j=1; j<zbin; j++)
            {
                tempeq[j] = tempb[j-1]/3. + tempb[j]/3. + tempb[j+1]/3.; //ms: smoothing
            }
          
            //ms23: different smoothing for TOA and BOA
            tempeq[0] = tempeq[1] + 1.0/3.0 * (tempeq[1] - tempeq[4]); //trend of smoothed T 1 through T 4
            tempeq[zbin] = tempeq[zbin-1] + 1.0/3.0 * (tempeq[zbin-1] - tempeq[zbin-4]); //trend of smoothed T N-1 through T N-4

            /*for (j=0; j<=zbin; j++) {
                tempeq[j] = tempb[j];
            }*/
            break;
        }
    
        /* Radiative-equilibrium iteration, report flux */
        //Rflux = dvector(0,nrl);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Radiative Transfer iteration
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //for (rstepcount=1; rstepcount<=NMAX_RC; rstepcount++) { /* roughly balance radiation between each convective update */
//ms: let's try without sub iteration:
        //for (rstepcount=1; rstepcount<=1; rstepcount++) { 
//ms: let's enforce rt equilibrium:
        int pevery=20, pcount=0;//Added for printing every few steps
        double dt[zbin+2]; //time step if TIME_STEPPING
        double Tsaved[zbin+1];
        int saveevery=20; //save Told every n steps for dt progression
        int isequil[zbin+1], sumisequil;

        for (j=1;j<=zbin+1;j++)  dt[j] = 1.0e+0; //for TIME_STEPPING: init in some random unit to ...

        rtstepcount=0;
        for (j=0;j<=zbin;j++) Tsaved[j] = tempb[j]; //store T for dt progression
        for (j=0;j<=zbin;j++) isequil[j] = 0; //store T for dt progression

        while ( (((Rfluxmax >= Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax >= Tol_RC)) || rtstepcount==0 || sumisequil < zbin) && rtstepcount<NMAX_RT ) { 
        //while ( (((Rfluxmax >= Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax >= Tol_RC)) || rtstepcount==0 ) && rtstepcount<1000 ) { 
        //while (rcount<=NMAX_RC ) { 
            rtstepcount +=1;
            pcount +=1;
            sumisequil=0;

//ms: end rt equil
            printf("%s %d\n", "Radiative loop RTstepCount = ", rtstepcount);
//                printf("%s\t%s\n","Layer","isconv");
//                for (j=zbin; j>=0; j--) {printf("%d\t%d\n",j,isconv[j]);}
            ms_RadFlux(Rflux, tempbnew, P, ncl, isconv, lapse, tempb, Tint, cp, dt, isequil);

//...............................................
//ms: sanity checks for temperature profile:
/*
//1st limit spikes:
            for (j=1; j<zbin; j++) { //check for spikes:
                if( abs(tempbnew[j] -  tempbnew[j+1]) > 0.2* fmin(tempbnew[j],tempbnew[j+1]) && abs(tempbnew[j-1] -  tempbnew[j]) > 0.2* fmin(tempbnew[j-1],tempbnew[j]) && abs(tempbnew[j-1] -  tempbnew[j+1]) < 0.2* fmin(tempbnew[j-1],tempbnew[j+1]) ) 
                {
                    printf("%s %d %s %f %s %f\n","T_Rad spike eliminated: layer ",j, "T_Rad ", tempbnew[j], "set to ", 0.5*(tempbnew[j+1]+tempbnew[j-1]));
                    tempbnew[j]=0.5*(tempbnew[j+1]+tempbnew[j-1]);
                }
            }
//2nd near surface T>= layer above:
            for (j=zbin; j>1; j--) { //modify layers above 10 bars pressure (stellar rad should mostly be absorbed at this point)
                if(P[j] > 1.0E+6 && tempbnew[j] > tempbnew[j-1]) 
                {
                    printf("%s %d %s %f %s %f\n","T_BOA unphysical: layer ",j, "T_Rad ", tempbnew[j-1], "set to ", 0.5*(tempbnew[j]+tempbnew[j-2]));
                    tempbnew[j-1]=0.5*(tempbnew[j]+tempbnew[j-2]);
                }
            }
            if(tempbnew[0] < tempbnew[1])
            {
                    printf("%s %d %s %f %s %f\n","T_BOA unphysical: layer ","0", "T_Rad ", tempbnew[0], "set to ", 2.0*tempbnew[1]-tempbnew[2]);
                    tempbnew[0]= 2.0*tempbnew[1] - tempbnew[2];

            }
*/            
//...............................................

            printf("%s\t%s\t%s\t%s\t%s\n","==== P[:]","T_old[:]","T_new[:]","Layer_equil","dt[:] ====");

            for (j=0; j<=zbin; j++) {
                printf("%e\t%f\t%f\t%d\t%f\n",P[j],tempb[j],tempbnew[j],isequil[j],dt[zbin+1-j]);
                sumisequil += isequil[j];
            }
            
            if (rtstepcount%saveevery==1) {//store T every n steps for dt progression
                for (j=0;j<=zbin;j++) {
                    if (fabs(tempbnew[j] - Tsaved[j]) < 0.5*saveevery*(fabs(tempbnew[j]-tempb[j]))) {
                        dt[zbin+1-j] /=1.5; //as in HELIOS
                    } else {
                        dt[zbin+1-j] *= 1.1; //as in HELIOS
                    }
                    Tsaved[j] = tempb[j];
                }
            }
            for (j=0; j<=zbin; j++) {
                tempb[j]=tempbnew[j];
            }
            //atexit(pexit);exit(0); //ms debugging mode
            Rfluxmax=0.0;
            dRfluxmax=0.0;
            radcount = 0;
            for (k=1; k<=zbin; k++) {
                if (isconv[k] == 0) {
                    radcount += 1;
                    Rfluxmax = fmax(Rfluxmax,fabs(Rflux[radcount]));
                    dRfluxmax = fmax(dRfluxmax,fabs((Rflux[radcount-1] - Rflux[radcount])));
                }
            }
            //for (j=0; j<=nrl; j++) { 
                //if (Rfluxmax<fabs(Rflux[j])) {
                //    Rfluxmax = fabs(Rflux[j]);
                //}
            //}
            printf("%s\t%s\t%s\t%s\n","Conv_loop","Rad_loop","Residual flux max","dFnet max");
            printf("%d\t\t%d\t\t%f\t%f\n", i, rtstepcount, Rfluxmax, dRfluxmax);
            printf("%s\n",filleq);
            if ( (Rfluxmax < Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax < Tol_RC) || (dRfluxmax < 0.2*Tol_RC)  ) {
                printf("Radiative transfer done\n");
                break;
            }
            if (pcount = pevery)
            {
            /* Print-out new TP profile*/
            FILE *fp;
            fp=fopen(outnewtemp,"w");
            for (j=0; j<=zbin; j++) fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempb[j]);
            fclose(fp);
            pcount = 0;
            }
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        } //ms22: END of radiative transfer iteration
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //ms22: smoothing temperature profile
        tempeq[0] = tempb[0]/2. + tempb[1]/2.; //surf
        tempeq[zbin] = tempb[zbin-1]/2. + tempb[zbin]/2.; //TOA
        for (j=1; j<zbin; j++)
        {
            tempeq[j] = tempb[j-1]/3. + tempb[j]/3. + tempb[j+1]/3.;
        }
        
        for (j=1; j<=zbin; j++)//update tl before convection
        {
            tl[j] = 0.5* (tempb[j]+tempb[j-1]);
        }

//atexit(pexit);exit(0); //ms debugging mode
    
        //free_dvector(Rflux,0,nrl);
    }
//*************************************************************
    
    for (j=0; j<=zbin; j++) {
        tempeq[j] = tempb[j];
    }

//ms22: smoothing temperature profile
    tempeq[0] = tempb[0]/2. + tempb[1]/2.; //surf
    tempeq[zbin] = tempb[zbin-1]/2. + tempb[zbin]/2.; //TOA
    for (j=1; j<zbin; j++)
    {
        tempeq[j] = tempb[j-1]/3. + tempb[j]/3. + tempb[j+1]/3.;
    }
    
    /* Print-out new pressure */
    FILE *fp;
    fp=fopen(outnewtemp,"w");
    for (j=0; j<=zbin; j++) {
        fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempeq[j]);
    }
    fclose(fp);
    
}
