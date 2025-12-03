#include <math.h>
#include "constant.h"

void GlobalBalance(double x[], int lx[], int laer[], int zone_r[], int zone_m[], int zone_t[], double **JJ, int zone_p[],
				   double Upflux[], double Loflux[], double Depo[], char outbalance[]);

void GlobalBalance(double x[], int lx[], int laer[], int zone_r[], int zone_m[], int zone_t[], double **JJ, int zone_p[],
				   double Upflux[], double Loflux[], double Depo[], char outbalance[])
{
	
	int i, j, jj, std;
	double rate, rater;
	double **P, **L;
	P=dmatrix(1, zbin, 1, NSP);
	L=dmatrix(1, zbin, 1, NSP);
	/* Initialization of matrix */
	for (i=1; i<=zbin; i++) {
		for (j=1; j<=NSP; j++) {
			P[i][j]=0;
			L[i][j]=0;
		}
	}
	
	/*Chemical Terms*/	
	for (i=1; i<=numr; i++) {
        std=zone_r[i]; /* Get the standard reaction number from list */
	    for (j=1; j<=zbin; j++) {
			rate=kk[j][std]*xx[j][ReactionR[std][1]]*xx[j][ReactionR[std][2]];
			L[j][ReactionR[std][1]] += kk[j][std]*xx[j][ReactionR[std][2]];
			L[j][ReactionR[std][2]] += kk[j][std]*xx[j][ReactionR[std][1]];
			P[j][ReactionR[std][4]] += rate;
			if (ReactionR[std][3]>1) {
				P[j][ReactionR[std][5]] += rate;
			}
			if (ReactionR[std][3]>2) {
				P[j][ReactionR[std][6]] += rate;
			}
			if (ReactionR[std][3]==1) {
				rater=Rkk[j][std]*xx[j][ReactionR[std][4]];
				P[j][ReactionR[std][1]] += rater;
				P[j][ReactionR[std][2]] += rater;
				L[j][ReactionR[std][4]] += Rkk[j][std];
			}
			if (ReactionR[std][3]==2) {
				rater=Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][5]];
				P[j][ReactionR[std][1]] += rater;
				P[j][ReactionR[std][2]] += rater;
				L[j][ReactionR[std][4]] += Rkk[j][std]*xx[j][ReactionR[std][5]];
				L[j][ReactionR[std][5]] += Rkk[j][std]*xx[j][ReactionR[std][4]];
			}
			if (ReactionR[std][3]==3) {
				rater=Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][5]]*xx[j][ReactionR[std][6]];
				P[j][ReactionR[std][1]] += rater;
				P[j][ReactionR[std][2]] += rater;
				L[j][ReactionR[std][4]] += Rkk[j][std]*xx[j][ReactionR[std][5]]*xx[j][ReactionR[std][6]];
				L[j][ReactionR[std][5]] += Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][6]];
				L[j][ReactionR[std][6]] += Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][5]];
			}
		}
	}
	
	for (i=1; i<=numm; i++) {
        std=zone_m[i];
		for (j=1; j<=zbin; j++) {
			rate=kkM[j][std]*xx[j][ReactionM[std][1]]*xx[j][ReactionM[std][2]];
			L[j][ReactionM[std][1]] += kkM[j][std]*xx[j][ReactionM[std][2]];
			L[j][ReactionM[std][2]] += kkM[j][std]*xx[j][ReactionM[std][1]];
			P[j][ReactionM[std][3]] += rate;
		    if (ReactionM[std][4]>0) {
				P[j][ReactionM[std][4]] += rate;
		    }
			if (ReactionM[std][4]==0) {
				rater=RkkM[j][std]*xx[j][ReactionR[std][3]];
				P[j][ReactionM[std][1]] += rater;
				P[j][ReactionM[std][2]] += rater;
				L[j][ReactionM[std][3]] += RkkM[j][std];
			}
			if (ReactionM[std][4]>0) {
				rater=RkkM[j][std]*xx[j][ReactionM[std][3]]*xx[j][ReactionM[std][4]];
				P[j][ReactionM[std][1]] += rater;
				P[j][ReactionM[std][2]] += rater;
				L[j][ReactionM[std][3]] += RkkM[j][std]*xx[j][ReactionM[std][4]];
				L[j][ReactionM[std][4]] += RkkM[j][std]*xx[j][ReactionM[std][3]];
			}
		}
	}
	
	for (i=1; i<=nump; i++) {
        std=zone_p[i];
		for (j=1; j<=zbin; j++) {
			rate=JJ[j][i]*xx[j][ReactionP[std][1]];
			L[j][ReactionP[std][1]] += JJ[j][i];
			P[j][ReactionP[std][2]] += rate;
			if (ReactionP[std][3]>0) {
				P[j][ReactionP[std][3]] += rate;
			}
			if (ReactionP[std][4]>0) {
				P[j][ReactionP[std][4]] += rate;
			}
		}
	}
	
	for (i=1; i<=numt; i++) {
        std=zone_t[i];
		for (j=1; j<=zbin; j++) {
			rate=kkT[j][std]*xx[j][ReactionT[std][1]];
			L[j][ReactionT[std][1]] += kkT[j][std];
			P[j][ReactionT[std][2]] += rate;
			if (ReactionT[std][3]>0) {
				P[j][ReactionT[std][3]] += rate;
			}
			if (ReactionT[std][3]==0) {
				rater=RkkT[j][std]*xx[j][ReactionT[std][2]];
				P[j][ReactionT[std][1]] += rater;
				L[j][ReactionT[std][2]] += RkkT[j][std];
			}
			if (ReactionT[std][3]>0) {
				rater=RkkT[j][std]*xx[j][ReactionT[std][2]]*xx[j][ReactionT[std][3]];
				P[j][ReactionT[std][1]] += rater;
				L[j][ReactionT[std][2]] += RkkT[j][std]*xx[j][ReactionT[std][3]];
				L[j][ReactionT[std][3]] += RkkT[j][std]*xx[j][ReactionT[std][2]];
			}
		}
	}
	
	/* Definition of fluxes (all in molecule cm-2 s-1) */
	double Outgassing[numx+1];
	double ChemProd[numx+1];
	double ChemLoss[numx+1];
	double DryDepo[numx+1];
	double WetDepo[numx+1];
	double Escape[numx+1];
	double Condensation[numx+1];
	double Net[numx+1];
	double GlobalRate[numx+1];
	
	/* Outgassing */
	for (i=1; i<=numx; i++) {
		Outgassing[i] = Loflux[i]*thickl;
	}
	
	/* Chemical Production and Loss */
	for (i=1; i<=numx; i++) {
		ChemProd[i]=0.0;
		ChemLoss[i]=0.0;
		for (j=1; j<=zbin; j++) {
			ChemProd[i] += P[j][lx[i]]*thickl;
			ChemLoss[i] -= L[j][lx[i]]*xx[j][lx[i]]*thickl;
		}
	}
	
	/* Dry Deposition including the falling of aerosols in the boundary layers */
	for (i=1; i<=numx; i++) {
		DryDepo[i] = -xx[1][lx[i]]*Depo[i]*thickl;
	}
	if (numa>0) {
		for (j=1; j<=numa; j++) {
			DryDepo[laer[j]] -= x[laer[j]]*VFall[1];
		}
	}
	
	/* Wet Deposition or rainout */
	for (i=1; i<=numx; i++) {
		WetDepo[i] = 0.0;
		for (j=1; j<=zbin; j++) {
			if (lx[i]!=7) {
				WetDepo[i] -= xx[j][lx[i]]*rainoutrate[j][lx[i]]*xx[j][7]*thickl;
			}
		}
	}
	
	/* Escape and Upper boundary conditions */
	for (i=1; i<=numx; i++) {
		Escape[i] = -(xx[zbin][lx[i]]*Vesc[lx[i]]+Upflux[i])*thickl;
	}
	
	/* Condensation */
	for (i=1; i<=numx; i++) {Condensation[i] = 0.0;}
	for (i=1; i<=numx; i++) {
		if (lx[i] == 7) {
			for (j=1; j<=zbin; j++) {
				if (xx[j][7] > nsH2O[j]) {
					Condensation[i] -= tcondfH2O[j]*(xx[j][7]-nsH2O[j])*xx[j][7]*thickl;
				}
			}
		}
	}
	if (numa>0) {
		for (j=1; j<=numa; j++) {
			/* H2SO4 */
			if (lx[laer[j]] == 78) {
				for (jj=1; jj<=numx; jj++) {
					if (lx[jj] == 73) {
						for (i=1; i<=zbin; i++) {
							if (xx[i][73] > nsH2SO4[i]) { /* condensation */
								Condensation[jj] -= tcondfH2SO4[i]*(xx[i][73]-nsH2SO4[i])*xx[i][73]*thickl;
								Condensation[laer[j]] += tcondfH2SO4[i]*(xx[i][73]-nsH2SO4[i])*xx[i][73]*thickl;
							} else { /* evaporation */
								Condensation[jj] += tcondfH2SO4[i]*(nsH2SO4[i]-xx[i][73])*xx[i][78]*thickl;
								Condensation[laer[j]] -= tcondfH2SO4[i]*(nsH2SO4[i]-xx[i][73])*xx[i][78]*thickl;
							}

						}
					}
				}
			}
			/* S8 */
			if (lx[laer[j]] == 111) {
				for (jj=1; jj<=numx; jj++) {
					if (lx[jj] == 79) {
						for (i=1; i<=zbin; i++) {
							if (xx[i][79] > nsS8[i]) { /* condensation */
								Condensation[jj] -= tcondfS8[i]*(xx[i][79]-nsS8[i])*xx[i][79]*thickl;
								Condensation[laer[j]] += tcondfS8[i]*(xx[i][79]-nsS8[i])*xx[i][79]*thickl;
							} else { /* evaporation */
								Condensation[jj] += tcondfS8[i]*(nsS8[i]-xx[i][79])*xx[i][111]*thickl;
								Condensation[laer[j]] -= tcondfS8[i]*(nsS8[i]-xx[i][79])*xx[i][111]*thickl;
							}
							
						}
					}
				}
			}
		}
	}
	
	
	/* Net flux */
	for (i=1; i<=numx; i++) {
		Net[i] = Outgassing[i]+ChemProd[i]+ChemLoss[i]+DryDepo[i]+WetDepo[i]+Escape[i]+Condensation[i];
	}
	
	/* Integration over the column and get the column density */
	double ntotal;
	for (i=1; i<=numx; i++) {
		ntotal = 0.0;
		for (j=1; j<=zbin; j++) {
			ntotal = ntotal + xx[j][lx[i]]*thickl;
		}
		GlobalRate[i] = Net[i]/ntotal;
	}
	
	/* Output */
	FILE *fp;
	fp=fopen(outbalance, "w");
	fprintf(fp, "%s\t%s\t%s\t%s\t%s\t\t%s\t\t%s\t\t%s\t%s\t%s\n", "Std No.", "Outgassing", "ChemProd", "ChemLoss", "DryDepo", "WetDepo", "Escape", "Condensation", "Net Flux", "Global Change");
	for (i=1; i<=numx; i++) {
		fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", lx[i], Outgassing[i],ChemProd[i],ChemLoss[i],DryDepo[i],WetDepo[i],Escape[i],Condensation[i],Net[i],GlobalRate[i]);
	}
	fclose(fp);

	
	free_dmatrix(P, 1, zbin, 1, NSP);
	free_dmatrix(L, 1, zbin, 1, NSP);
	
}