/*Function to provide solar radiation at a given height*/
#include <math.h>
#include "constant.h"
#include "gasheat.c"
#include "RadFlux.c"

void RadConv(double tempeq[], double P[], double Tint, double **xxf, int stdn[], double qysum[], double **cross, double **crosst);

void RadConv(double tempeq[], double P[], double Tint, double **xxf, int stdn[], double qysum[], double **cross, double **crosst)
{
	int i,j,k;
    
    /* Initial Guess */
    double tempb[zbin+1], tempbnew[zbin+1], tempc[zbin+1], presb[zbin+1];
    for (j=1; j<zbin; j++) {
		tempb[j] = (tl[j] + tl[j+1])/2.0;
	}
	tempb[0] = 1.5*tl[1] - 0.5*tl[2];
	tempb[zbin] = 1.5*tl[zbin] - 0.5*tl[zbin-1];
    for (j=1; j<=zbin; j++) {
		tempc[j] = tl[j];
	}
    
    double heliumnumber[zbin+1], meanmolecular[zbin+1];
    for (j=1; j<=zbin; j++) {
        heliumnumber[j] = MM[j];
        for (i=1; i<=NSP; i++) {
            heliumnumber[j] -= xx[j][i];
        }
        if (heliumnumber[j]<0.0) {
            heliumnumber[j] = 0.0;
        }
    }
    for (j=1; j<=zbin; j++) {
        meanmolecular[j] = (2.0*xx[j][53]+4.0*heliumnumber[j])/(xx[j][53]+heliumnumber[j]);
        /*printf("%s %f\n","Mean Molecular",meanmolecular[j]);*/
    }
    
    /* Define variables for iteration */
    double znew[zbin+1], scaleheight;
    int isconv[zbin+1], ncl, nrl; /* if a layer is convective */
    double lapse[zbin+1], dtdz, gasheat, gasheat1, gasheat2;
    double *Rflux;
    double Rfluxmax;
    double romega=0.1;
    double imax=100;
    
    /* compute hydrostatitc znew */
	/* znew[0]=0.0;
	for (j=1; j<=zbin; j++) {
		scaleheight = KBOLTZMANN * tempc[j] / AIRM / AMU / GA / 1000 ;
		znew[j] = znew[j-1] - scaleheight*log(P[j]/P[j-1]);
	} */
    /* for (j=0; j<=zbin; j++) {
     printf("%f %f\n",znew[j],tempb[j]);
     } */
    
    /* Determine the layer of convection */
    int IFRC=1;
    
    if (IFRC==0) {
    ncl = 0;
    nrl = 0;
    for (j=1; j<=zbin; j++) {
        if (RefIdxType == 0) { gasheat=AirHeat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 1) { gasheat=CO2Heat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 2) { gasheat=HeHeat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 3) { gasheat=N2Heat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 4) { gasheat=NH3Heat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 5) { gasheat=CH4Heat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 6) { gasheat=H2Heat(tempc[j])*1000.0/AIRM;}
        if (RefIdxType == 7) { gasheat=O2Heat(tempc[j])*1000.0/AIRM;} /* J/K/kg */
        if (RefIdxType == 8) { gasheat=COHeat(tempc[j])*1000.0/AIRM;} /* J/K/kg */
        if (RefIdxType == 9) { gasheat=H2OHeat(tempc[j])*1000.0/AIRM;} /* J/K/kg */
        gasheat1 = H2Heat(tempc[j]);
        gasheat2 = HeHeat(tempc[j]);
        gasheat  = (gasheat1*xx[j][53]+gasheat2*heliumnumber[j])/(xx[j][53]+heliumnumber[j])*1000.0/meanmolecular[j];
        lapse[j] = KBOLTZMANN/meanmolecular[j]/AMU/gasheat ; /* dimension less dlog(T)/dlog(P) */
    }
    for (j=zbin; j>=1; j--) {
        if ( (log(tempb[j-1]/tempb[j])/log(P[j-1]/P[j]) > lapse[j] || j<0 ) && j<zbin ) {
            tempb[j-1] = tempb[j] * pow(P[j-1]/P[j],lapse[j]);
            isconv[j] = 1;
            ncl = ncl+1;
        } else {
            isconv[j] = 0;
            nrl = nrl+1;
        }
    }
    for (j=1; j<=zbin; j++) {
        printf("%d\n",isconv[zbin+1-j]);
    }
    printf("%s %d %s %d\n", "conv layer", ncl, "rad layer", nrl);
    }
    
    for (i=1; i<=imax; i++) {
        
        if (IFRC==1) {

        /* Determine the layer of convection */
        ncl = 0;
        nrl = 0;
        for (j=1; j<=zbin; j++) {
            if (RefIdxType == 0) { gasheat=AirHeat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 1) { gasheat=CO2Heat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 2) { gasheat=HeHeat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 3) { gasheat=N2Heat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 4) { gasheat=NH3Heat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 5) { gasheat=CH4Heat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 6) { gasheat=H2Heat(tempc[j])*1000.0/AIRM;}
            if (RefIdxType == 7) { gasheat=O2Heat(tempc[j])*1000.0/AIRM;} /* J/K/kg */
            if (RefIdxType == 8) { gasheat=COHeat(tempc[j])*1000.0/AIRM;} /* J/K/kg */
            if (RefIdxType == 9) { gasheat=H2OHeat(tempc[j])*1000.0/AIRM;} /* J/K/kg */
            gasheat1 = H2Heat(tempc[j]);
            gasheat2 = HeHeat(tempc[j]);
            gasheat  = (gasheat1*xx[j][53]+gasheat2*heliumnumber[j])/(xx[j][53]+heliumnumber[j])*1000.0/meanmolecular[j];
            /* printf("%f\t%f\t%f\t%f\n",gasheat1,gasheat2,gasheat,gasheat1*1000.0/meanmolecular[j]); */
            lapse[j] = KBOLTZMANN/meanmolecular[j]/AMU/gasheat ; /* dimension less dlog(T)/dlog(P) */
        }
            for (j=zbin; j>=1; j--) {
                printf("%f %f\n",log(tempb[j-1]/tempb[j])/log(P[j-1]/P[j]),lapse[j]);
            }
        for (j=zbin; j>=1; j--) {
            if ( (log(tempb[j-1]/tempb[j])/log(P[j-1]/P[j]) > lapse[j]*0.95 || j<0 ) && j<zbin ) {
                tempb[j-1] = tempb[j] * pow(P[j-1]/P[j],lapse[j]);
                isconv[j] = 1;
                ncl = ncl+1;
            } else {
                isconv[j] = 0;
                nrl = nrl+1;
            }
        }
        for (j=1; j<=zbin; j++) {
            printf("%d\n",isconv[zbin+1-j]);
        }
        printf("%s %d %s %d\n", "conv layer", ncl, "rad layer", nrl);
            
        }
    
    /* Calculate midpoint fluxes and Jacobian */
    Rflux = dvector(0,nrl);
    RadFlux(Rflux, tempbnew, P, ncl, isconv, lapse, xxf, tempb, Tint, stdn, qysum, cross, crosst, romega);
    for (j=0; j<=nrl; j++) {
        printf("%f\t",Rflux[j]);
        /*for (k=0; k<=5; k++) {
            printf("%f\t",FJacob[j][k]);
        }*/
        printf("\n");
    }
    for (j=0; j<=zbin; j++) {
        printf("%e\t%f\t%f\n",P[j],tempb[j],tempbnew[j]);
    }
    
    /* Determine whether converged */
    Rfluxmax=0;
    for (j=0; j<=nrl; j++) {
        if (Rfluxmax<fabs(Rflux[j])) {
            Rfluxmax = fabs(Rflux[j]);
        }
    }
    if ( (Rfluxmax < 1.0E-3 * SIGMA*pow(Tint,4.0))  || (Rfluxmax < 1.0E-2)  ) {
        printf("%s %d %f\n", "Residual flux max", i, Rfluxmax);
        printf("converged\n");
        for (j=0; j<=zbin; j++) {
            tempeq[j] = tempb[j];
        }
        break;
    } else {
        printf("%s %d %f\n", "Residual flux max", i, Rfluxmax);
        for (j=0; j<=zbin; j++) {
            tempb[j] = tempbnew[j];
        }
        for (j=1; j<=zbin; j++) {
            tempc[j] = 0.5*tempb[j]+0.5*tempb[j-1];
        }
    }
    
    free_dvector(Rflux,0,nrl);
        
    }
    
    for (j=0; j<=zbin; j++) {
        tempeq[j] = tempb[j];
    }
    
}