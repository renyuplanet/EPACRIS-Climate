
/*Function to calculate temperature at a given height with grey atmosphere assumption */

#include <math.h>
//#include "gasheat.c"

void GreyTemp(double P[], char outnewtemp[], double tempint);

void GreyTemp(double P[], char outnewtemp[], double tempint)
{

	int j, i;
    
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
	
	double Tirr, Tnew4, znew[zbin+1];
	Tirr = STAR_TEMP/pow(ORBIT,0.5);
	
	double Tint;
	Tint=TINTSET;
    Tint=tempint;
	printf("%s %f\n", "Tint is", Tint);
	
	double miu, gamma;
	miu=cos(THETAREF/2.0);
	
	double GA, scaleheight;
	GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET;
	
	double tau[zbin+1], kth[zbin+1], ktv[zbin+1];
	
	for (j=1; j<=zbin; j++) {
        
		kth[j]=0.0;
		kth[j] += MeanCO2[j]*xx[j][52]/MM[j];  /* CO2  */
		kth[j] += MeanO2[j]*xx[j][54]/MM[j];   /* O2   */
		kth[j] += MeanH2O[j]*xx[j][7]/MM[j];   /* H2O  */
		kth[j] += MeanOH[j]*xx[j][4]/MM[j];    /* OH   */
		kth[j] += MeanH2CO[j]*xx[j][22]/MM[j]; /* H2CO */
		kth[j] += MeanH2O2[j]*xx[j][6]/MM[j];  /* H2O2 */
		kth[j] += MeanHO2[j]*xx[j][5]/MM[j];   /* HO2  */
		kth[j] += MeanCO[j]*xx[j][20]/MM[j];   /* CO   */
		kth[j] += MeanO3[j]*xx[j][2]/MM[j];    /* O3   */
		kth[j] += MeanCH4[j]*xx[j][21]/MM[j];   /* CH4   */
		kth[j] += MeanC2H2[j]*xx[j][27]/MM[j];
		kth[j] += MeanC2H4[j]*xx[j][29]/MM[j];
		kth[j] += MeanC2H6[j]*xx[j][31]/MM[j];
		kth[j] += MeanCH2O2[j]*xx[j][23]/MM[j];
        
        kth[j] += MeanHCN[j]*xx[j][37]/MM[j];
        kth[j] += MeanNH3[j]*xx[j][9]/MM[j];
        kth[j] += MeanN2O[j]*xx[j][11]/MM[j];
        kth[j] += MeanNO[j]*xx[j][12]/MM[j];
        kth[j] += MeanNO2[j]*xx[j][13]/MM[j];
        kth[j] += MeanHNO3[j]*xx[j][18]/MM[j];
        kth[j] += MeanH2S[j]*xx[j][45]/MM[j];
        kth[j] += MeanSO2[j]*xx[j][43]/MM[j];
        kth[j] += MeanOCS[j]*xx[j][49]/MM[j];
        
        //printf("%s %d %e\n", "MM Lay", j, MM[j]);
        /*printf("%s %f %e %e\n", "CO2", zl[j], MeanCO2[j], xx[j][52]);
        printf("%s %f %e %e\n", "O2", zl[j], MeanO2[j], xx[j][54]);
        printf("%s %f %e %e\n", "H2O", zl[j], MeanH2O[j], xx[j][7]);
        printf("%s %f %e %e\n", "OH", zl[j], MeanOH[j], xx[j][4]);
        printf("%s %f %e %e\n", "H2CO", zl[j], MeanH2CO[j], xx[j][22]);
        printf("%s %f %e %e\n", "H2O2", zl[j], MeanH2O2[j], xx[j][6]);
        printf("%s %f %e %e\n", "HO2", zl[j], MeanHO2[j], xx[j][5]);
        printf("%s %f %e %e\n", "CO", zl[j], MeanCO[j], xx[j][20]);
        printf("%s %f %e %e\n", "O3", zl[j], MeanO3[j], xx[j][2]);
        printf("%s %f %e %e\n", "CH4", zl[j], MeanCH4[j], xx[j][21]);
        printf("%s %f %e %e\n", "C2H2", zl[j], MeanC2H2[j], xx[j][27]);
        printf("%s %f %e %e\n", "C2H4", zl[j], MeanC2H4[j], xx[j][29]);
        printf("%s %f %e %e\n", "C2H6", zl[j], MeanC2H6[j], xx[j][31]);
        printf("%s %f %e %e\n", "CH2O2", zl[j], MeanCH2O2[j], xx[j][23]);
        printf("%s %f %e %e\n", "HCN", zl[j], MeanHCN[j], xx[j][37]);
        printf("%s %f %e %e\n", "NH3", zl[j], MeanNH3[j], xx[j][9]);
        printf("%s %f %e %e\n", "N2O", zl[j], MeanN2O[j], xx[j][11]);
        printf("%s %f %e %e\n", "NO", zl[j], MeanNO[j], xx[j][12]);
        printf("%s %f %e %e\n", "NO2", zl[j], MeanNO2[j], xx[j][13]);
        printf("%s %f %e %e\n", "HNO3", zl[j], MeanHNO3[j], xx[j][18]);
        printf("%s %f %e %e\n", "H2S", zl[j], MeanH2S[j], xx[j][45]);
        printf("%s %f %e %e\n", "SO2", zl[j], MeanSO2[j], xx[j][43]);
        printf("%s %f %e %e\n", "OCS", zl[j], MeanOCS[j], xx[j][49]); */
        
	kth[j] += MeanH2H2CIA[j]*xx[j][53]*xx[j][53]/MM[j]; /* H2H2CIA */
        kth[j] += MeanH2HeCIA[j]*xx[j][53]*heliumnumber[j]/MM[j]; /* H2HeCIA */
        kth[j] += MeanH2HCIA[j]*xx[j][53]*xx[j][3]/MM[j]; /* H2HCIA */
        kth[j] += MeanN2N2CIA[j]*xx[j][55]*xx[j][55]/MM[j]; /* N2N2CIA */
        kth[j] += MeanN2H2CIA[j]*xx[j][55]*xx[j][53]/MM[j]; /* N2H2CIA */
        kth[j] += MeanCO2CO2CIA[j]*xx[j][52]*xx[j][52]/MM[j]; /* CO2CO2CIA */
        
        /*printf("%s %f %e %e\n", "H2H2CIA", zl[j], MeanH2H2CIA[j], xx[j][53]);
        printf("%s %f %e %e\n", "H2HeCIA", zl[j], MeanH2HeCIA[j], heliumnumber[j]);
        printf("%s %f %e %e\n", "H2HCIA", zl[j], MeanH2HCIA[j], xx[j][3]);
        printf("%s %f %e %e\n", "N2N2CIA", zl[j], MeanN2N2CIA[j], xx[j][55]);
        printf("%s %f %e %e\n", "N2H2CIA", zl[j], MeanN2H2CIA[j], xx[j][55]);
        printf("%s %f %e %e\n", "CO2CO2CIA", zl[j], MeanCO2CO2CIA[j], xx[j][52]); */
        
        /*kth[j] += MeanN2N2CIA[j]*xx[j][54]*xx[j][54]/MM[j]; */
		
	}
    
    for (j=1; j<=zbin; j++) {
        
        ktv[j]=0.0;
        ktv[j] += SMeanCO2[j]*xx[j][52]/MM[j];  /* CO2  */
        ktv[j] += SMeanO2[j]*xx[j][54]/MM[j];   /* O2   */
        ktv[j] += SMeanH2O[j]*xx[j][7]/MM[j];   /* H2O  */
        ktv[j] += SMeanOH[j]*xx[j][4]/MM[j];    /* OH   */
        ktv[j] += SMeanH2CO[j]*xx[j][22]/MM[j]; /* H2CO */
        ktv[j] += SMeanH2O2[j]*xx[j][6]/MM[j];  /* H2O2 */
        ktv[j] += SMeanHO2[j]*xx[j][5]/MM[j];   /* HO2  */
        ktv[j] += SMeanCO[j]*xx[j][20]/MM[j];   /* CO   */
        ktv[j] += SMeanO3[j]*xx[j][2]/MM[j];    /* O3   */
        ktv[j] += SMeanCH4[j]*xx[j][21]/MM[j];   /* CH4   */
        ktv[j] += SMeanC2H2[j]*xx[j][27]/MM[j];
        ktv[j] += SMeanC2H4[j]*xx[j][29]/MM[j];
        ktv[j] += SMeanC2H6[j]*xx[j][31]/MM[j];
        ktv[j] += SMeanCH2O2[j]*xx[j][23]/MM[j];
        
        ktv[j] += SMeanHCN[j]*xx[j][37]/MM[j];
        ktv[j] += SMeanNH3[j]*xx[j][9]/MM[j];
        ktv[j] += SMeanN2O[j]*xx[j][11]/MM[j];
        ktv[j] += SMeanNO[j]*xx[j][12]/MM[j];
        ktv[j] += SMeanNO2[j]*xx[j][13]/MM[j];
        ktv[j] += SMeanHNO3[j]*xx[j][18]/MM[j];
        ktv[j] += SMeanH2S[j]*xx[j][45]/MM[j];
        ktv[j] += SMeanSO2[j]*xx[j][43]/MM[j];
        ktv[j] += SMeanOCS[j]*xx[j][49]/MM[j];
        
        ktv[j] += SMeanH2H2CIA[j]*xx[j][53]*xx[j][53]/MM[j]; /* H2H2CIA */
        ktv[j] += SMeanH2HeCIA[j]*xx[j][53]*heliumnumber[j]/MM[j]; /* H2HeCIA */
        ktv[j] += SMeanH2HCIA[j]*xx[j][53]*xx[j][3]/MM[j]; /* H2HCIA */
        ktv[j] += SMeanN2N2CIA[j]*xx[j][55]*xx[j][55]/MM[j]; /* N2N2CIA */
        ktv[j] += SMeanN2H2CIA[j]*xx[j][55]*xx[j][53]/MM[j]; /* N2H2CIA */
        ktv[j] += SMeanCO2CO2CIA[j]*xx[j][52]*xx[j][52]/MM[j]; /* CO2CO2CIA */
        
        /*ktv[j] += SMeanN2N2CIA[j]*xx[j][54]*xx[j][54]/MM[j]; */
        
    }
	
	for (j=1; j<=zbin; j++) {
		printf("%s %f %s %e %e %f\n", "The mean opacity at altitude", zl[j], "km is", ktv[j], kth[j], ktv[j]/kth[j]);
	}
	
	double kvtotal, kthtotal, tauint;
	kvtotal=0.0;
	kthtotal=0.0;
	tauint=0.0;
	for (j=zbin; j>=1; j--) {
		if (tauint<100.0) {
			kvtotal += ktv[j]*NAVOGADRO/(meanmolecular[j]*0.001)/GA*(P[j-1]-P[j])*1.0E-4;
			kthtotal+= kth[j]*NAVOGADRO/(meanmolecular[j]*0.001)/GA*(P[j-1]-P[j])*1.0E-4;
			tauint += kth[j]*NAVOGADRO/(meanmolecular[j]*0.001)/GA*(P[j-1]-P[j])*1.0E-4;
		}
	}
	gamma = kvtotal/kthtotal*0.35;
	printf("%s %f\n","The gamma factor is", gamma);
    
    /* correction for Rossland mean opacity */
    double logpp;
    for (j=1; j<=zbin; j++) {
        logpp = log10(pl[j]);
        if (logpp<1.5) {
            logpp=1.5;
        }
        if (logpp>8.0) {
            logpp=8.0;
        }
        kth[j] = kth[j]* pow(10.0, -0.049*logpp*logpp + 0.97*logpp - 5.6);
    }
	
	tau[zbin]=0;
	j=zbin-1;
	while (j>=0) {
		/* tau[j] = tau[j+1]+kth[j+1]*thickl*MM[j+1]; */
		tau[j] = tau[j+1]+kth[j+1]*NAVOGADRO/(meanmolecular[j+1]*0.001)/GA*(P[j]-P[j+1])*1.0E-4;
        printf("%s %e %e\n","The optical depth is", P[j],tau[j]);
		j=j-1;
	}
	
	FILE *fin;
	int s;
	fin=fopen("Library/E2/E2.dat","r");
	s=LineNumber(fin, 1000);
	fclose(fin);
	fin=fopen("Library/E2/E2.dat","r");
	double zx[s];
	double E2y[s];
	double E2[1];
	double gt[1];
	GetData(fin, 1000, s, zx, E2y);
	fclose(fin);
	
	for (j=0; j<=zbin; j++) {
		
		Tnew4=3.0*pow(Tirr,4.0)/4.0*miu*(2.0/3.0+miu/gamma+(gamma/3.0/miu-miu/gamma)*exp(-gamma*tau[j]/miu))*FADV*(1.0-PAB);
		
		if (gamma*tau[j]<=zx[0]) {
			E2[0]=1.0;
		} else if (gamma*tau[j]>=zx[s-1]) {
			E2[0]=0.0;
		} else {
			gt[0]=gamma*tau[j];
			Interpolation(gt, 1, E2, zx, E2y, s, 0);
		}

		Tnew4=3.0*pow(Tirr,4.0)/4.0*FADV*(1.0-PAB)*(2.0/3.0 + 2.0/3.0/gamma*(1+(gamma*tau[j]/2.0-1.0)*exp(-gamma*tau[j])) + 2.0*gamma/3.0*(1-tau[j]*tau[j]/2.0)*E2[0] );
		
		Tnew4 += 3.0*pow(Tint,4.0)/4.0*(2.0/3.0+tau[j]);
		Tnew[j]=pow(Tnew4,0.25);
		
		/* printf("%s %d %e %f\n", "Tnew", j, P[j], Tnew[j]); */
		
	}
	
	/* check adiabats */
	double lapse[zbin+1], gasheat;
    for (j=1; j<=zbin; j++) {
        gasheat=(CO2Heat(tl[j])*xx[j][52]+HeHeat(tl[j])*heliumnumber[j]+N2Heat(tl[j])*xx[j][55]+NH3Heat(tl[j])*xx[j][9]+CH4Heat(tl[j])*xx[j][21]+H2Heat(tl[j])*xx[j][53]+O2Heat(tl[j])*xx[j][54]+COHeat(tl[j])*xx[j][20]+H2OHeat(tl[j])*xx[j][7])/(xx[j][52]+heliumnumber[j]+xx[j][55]+xx[j][9]+xx[j][21]+xx[j][53]+xx[j][54]+xx[j][20]+xx[j][7])*1000.0/meanmolecular[j];
        lapse[j] = KBOLTZMANN/meanmolecular[j]/AMU/gasheat ; /* dimension less dlog(T)/dlog(P) */
    }
    /* determine convection, adjust */
    for (j=zbin; j>=1; j--) {
        if ( Tnew[j-1] >= (Tnew[j] * pow(P[j-1]/P[j],lapse[j])) ) {
            Tnew[j-1] = Tnew[j] * pow(P[j-1]/P[j],lapse[j]);
        }
    }
    /* compute hydrostatitc znew */
    znew[0]=0.0;
    for (j=1; j<=zbin; j++) {
        scaleheight = KBOLTZMANN * (Tnew[j-1]+Tnew[j]) /2.0 / meanmolecular[j] / AMU / GA /1000.0 ;
        znew[j] = znew[j-1] - scaleheight*log(P[j]/P[j-1]);
    }
	
	/* Print-out new pressure */
	FILE *fp;
	fp=fopen(outnewtemp,"w");
	for (j=0; j<=zbin; j++) {
		fprintf(fp, "%f\t%f\t%f\n", znew[j], log10(P[j]), Tnew[j]);
	}
	fclose(fp);
	
    printf("%s\n","GreyTemp DONE");
	
}
