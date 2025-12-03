/*Function to provide solar radiation at a given height*/
#include <math.h>
#include "constant.h"
#include "ludcmp.c"
#include "lubksb.c"
#include "RadFlux_Temp2.c"

void RadConv(double tempeq[], double P[], double T[], double Tint, char outnewtemp[]);

void RadConv(double tempeq[], double P[], double T[], double Tint, char outnewtemp[])
{
	int i,j,k,rcount;
    
    /* Initial profile */
    double tempb[zbin+1], tempbnew[zbin+1], tempc[zbin+1];
    for (j=0; j<=zbin; j++) {
		tempb[j] = T[j];
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
    double lapse[zbin+1], gasheat;
    double *Rflux;
    double Rfluxmax;
    
    int isconv[zbin+1], ncl, nrl; /* if a layer is convective */
    int isconv_old[zbin+1];
    int isconv_sum;
    
    for (j=1; j<=zbin; j++) {
        isconv[j] = 0;
    }
    
    /* gas heat capacity and laspe rate */
    for (j=1; j<=zbin; j++) {
        gasheat=(CO2Heat(tl[j])*xx[j][52]+HeHeat(tl[j])*heliumnumber[j]+N2Heat(tl[j])*xx[j][55]+NH3Heat(tl[j])*xx[j][9]+CH4Heat(tl[j])*xx[j][21]+H2Heat(tl[j])*xx[j][53]+O2Heat(tl[j])*xx[j][54]+COHeat(tl[j])*xx[j][20]+H2OHeat(tl[j])*xx[j][7])/(xx[j][52]+heliumnumber[j]+xx[j][55]+xx[j][9]+xx[j][21]+xx[j][53]+xx[j][54]+xx[j][20]+xx[j][7])*1000.0/meanmolecular[j];
        lapse[j] = KBOLTZMANN/meanmolecular[j]/AMU/gasheat ; /* dimension less dlog(T)/dlog(P) */
    }
    
    for (i=1; i<=NMAX_RC; i++) {
            
        /* store last RC boudary */
        for (j=1; j<=zbin; j++) {
            isconv_old[j] = isconv[j];
        }
        /* new RC boundary */
        ncl = 0;
        nrl = 0;
        /* determine convection, and record convective layer */
        for (j=zbin; j>=1; j--) {
            if ( tempb[j-1] >= (tempb[j] * pow(P[j-1]/P[j],lapse[j]) * 0.999999) ) {
                isconv[j] = 1;
                ncl = ncl+1;
            } else {
                isconv[j] = 0;
                nrl = nrl+1;
            }
        }
        /* process isconv to remove single radiative layer between convective layers */
        for (j=2; j<zbin; j++) {
            if (isconv[j]==0 && (isconv[j+1]==1 && isconv[j-1]==1)) {
                isconv[j]=1;
                ncl = ncl+1;
                nrl = nrl-1;
            }
        }
        
        /*test*/
        if (i==1) {
            for (j=1; j<=60; j++) {
                isconv[j] = 0;
            }
            ncl = 0;
            nrl = 60;
        }
        
        
        /* adjust temperature */
        for (j=zbin; j>=1; j--) {
            if ( isconv[j]==1 ) {
                tempb[j-1] = tempb[j] * pow(P[j-1]/P[j],lapse[j]);
            }
        }
        /* hydrostatic equilibrium */
        znew[0]=0.0;
        for (j=1; j<=zbin; j++) {
            tempc[j] = (tempb[j]+tempb[j-1])/2.0;
            scaleheight = KBOLTZMANN * tempc[j] / meanmolecular[j] / AMU / GA / 1000 ;
            znew[j] = znew[j-1] - scaleheight*log(P[j]/P[j-1]);
        }
        /* print out new RC boundary */
        for (j=1; j<=zbin; j++) {
            printf("%d\n",isconv[zbin+1-j]);
        }
        printf("%s %d %s %d\n", "conv layer", ncl, "rad layer", nrl);
        /* compare old and new RC boundary */
        isconv_sum=0;
        for (j=1; j<=zbin; j++) {
            isconv_sum += abs(isconv_old[j]-isconv[j]);
        }
        printf("%s %d\n", "change in convective layer is", isconv_sum);
        /* determine whether converged */
        if (isconv_sum == 0 && i!=1) {
            printf("converged\n");
            for (j=0; j<=zbin; j++) {
                tempeq[j] = tempb[j];
            }
            break;
        }
        /* Radiative-equilibrium iteration, report flux */
        Rflux = dvector(0,nrl);
        for (rcount=1; rcount<=100; rcount++) { /* roughlt balance radiation between each convective update */
            RadFlux(Rflux, tempbnew, P, ncl, isconv, lapse, tempb, Tint);
            Rfluxmax=0;
            for (j=0; j<=nrl; j++) {
                if (Rfluxmax<fabs(Rflux[j])) {
                    Rfluxmax = fabs(Rflux[j]);
                }
            }
            printf("%s %d %d %f\n", "Residual flux max", i, rcount, Rfluxmax);
            if ( (Rfluxmax < Tol_RC_R * SIGMA*pow(Tint,4.0))  || (Rfluxmax < Tol_RC)  ) {
                for (j=0; j<=zbin; j++) {
                    printf("%e\t%f\t%f\n",P[j],tempb[j],tempbnew[j]);
                    tempb[j]=tempbnew[j];
                }
                printf("radiative equilibrium reached\n");
                break;
            } else {
                for (j=0; j<=zbin; j++) {
                    printf("%e\t%f\t%f\n",P[j],tempb[j],tempbnew[j]);
                    tempb[j]=tempbnew[j];
                }
            }
        }
        free_dvector(Rflux,0,nrl);
    }
    
    for (j=0; j<=zbin; j++) {
        tempeq[j] = tempb[j];
    }
    
    /* Print-out new pressure */
    FILE *fp;
    fp=fopen(outnewtemp,"w");
    for (j=0; j<=zbin; j++) {
        fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), tempeq[j]);
    }
    fclose(fp);
    
}