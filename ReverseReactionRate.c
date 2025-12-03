/* Function that computes reverse reaction rates */

#include <math.h>
#include "constant.h"

void ReverseReactionRate(int zone_r[], int zone_m[], int zone_t[]);

void ReverseReactionRate(int zone_r[], int zone_m[], int zone_t[])
{
	int i,j;
	int nn;
	int nr;
	int np;
	double deltaG;
	double kbmodi=1.38065E-22;
	
	/* Initialization */
	for (j=1; j<=zbin; j++) {
		for (i=1; i<=NKin; i++) {
			Rkk[j][i] = 0.0;
		}
		for (i=1; i<=NKinM; i++) {
			RkkM[j][i] = 0.0;
		}
		for (i=1; i<=NKinT; i++) {
			RkkT[j][i] = 0.0;
		}
	}
	
	/* Compute reverse rates for 2-molecule forward reactions */
	for (i=1; i<=numr; i++) {
		nn=zone_r[i];
		nr=2;
		np=ReactionR[nn][3];
		if (ReactionR[nn][1]!=74 && ReactionR[nn][1]!=75 && ReactionR[nn][1]!=110 && ReactionR[nn][1]!=70 && ReactionR[nn][1]!=105) {
			if (ReactionR[nn][2]!=74 && ReactionR[nn][2]!=75 && ReactionR[nn][2]!=110 && ReactionR[nn][2]!=70 && ReactionR[nn][2]!=105) {
				if (ReactionR[nn][4]!=74 && ReactionR[nn][4]!=75 && ReactionR[nn][4]!=110 && ReactionR[nn][4]!=70 && ReactionR[nn][4]!=105) {
					if (ReactionR[nn][5]!=74 && ReactionR[nn][5]!=75 && ReactionR[nn][5]!=110 && ReactionR[nn][5]!=70 && ReactionR[nn][5]!=105) {
						if (ReactionR[nn][6]!=74 && ReactionR[nn][6]!=75 && ReactionR[nn][6]!=110 && ReactionR[nn][6]!=70 && ReactionR[nn][6]!=105) {
							for (j=1; j<=zbin; j++) {
								deltaG = 0.0;
								deltaG -= (GibbsForm[ReactionR[nn][1]][j] + GibbsForm[ReactionR[nn][2]][j]);
								deltaG += GibbsForm[ReactionR[nn][4]][j];
								if (np>1) {
									deltaG += GibbsForm[ReactionR[nn][5]][j];
								}
								if (np>2) {
									deltaG += GibbsForm[ReactionR[nn][6]][j];
								}
								Rkk[j][nn] = kk[j][nn]*exp(deltaG/R_GAS/tl[j])*pow(kbmodi*tl[j],np-nr);
							}
						}
					}
				}
			}
		}
	}
	

	/* Compute reverse rates for 3-molecule forward reactions */
	for (i=1; i<=numm; i++) {
		nn=zone_m[i];
		nr=2;
		if (ReactionM[nn][4]>0) {
			np=2;
		}
		else {
			np=1;
		}
		if (ReactionM[nn][1]!=74 && ReactionM[nn][1]!=75 && ReactionM[nn][1]!=110 && ReactionM[nn][1]!=70 && ReactionM[nn][1]!=105) {
			if (ReactionM[nn][2]!=74 && ReactionM[nn][2]!=75 && ReactionM[nn][2]!=110 && ReactionM[nn][2]!=70 && ReactionM[nn][2]!=105) {
				if (ReactionM[nn][4]!=74 && ReactionM[nn][4]!=75 && ReactionM[nn][4]!=110 && ReactionM[nn][4]!=70 && ReactionM[nn][4]!=105) {
					if (ReactionM[nn][3]!=74 && ReactionM[nn][3]!=75 && ReactionM[nn][3]!=110 && ReactionM[nn][3]!=70 && ReactionM[nn][3]!=105) {
						for (j=1; j<=zbin; j++) {
							deltaG = 0.0;
							deltaG -= (GibbsForm[ReactionM[nn][1]][j] + GibbsForm[ReactionM[nn][2]][j]);
							deltaG += GibbsForm[ReactionM[nn][3]][j];
							if (np>1) {
								deltaG += GibbsForm[ReactionM[nn][4]][j];
							}
							RkkM[j][nn] = kkM[j][nn]*exp(deltaG/R_GAS/tl[j])*pow(kbmodi*tl[j],np-nr);
						}
					}
				}
			}
		}
	}
	
	
	/* Compute reverse rates for dissociative forward reactions */
	for (i=1; i<=numt; i++) {
		nn=zone_t[i];
		nr=1;
		if (ReactionT[nn][3]>0) {
			np=2;
		}
		else {
			np=1;
		}
		if (ReactionT[nn][1]!=74 && ReactionT[nn][1]!=75 && ReactionT[nn][1]!=110 && ReactionT[nn][1]!=70 && ReactionT[nn][1]!=105) {
			if (ReactionT[nn][2]!=74 && ReactionT[nn][2]!=75 && ReactionT[nn][2]!=110 && ReactionT[nn][2]!=70 && ReactionT[nn][2]!=105) {
					if (ReactionT[nn][3]!=74 && ReactionT[nn][3]!=75 && ReactionT[nn][3]!=110 && ReactionT[nn][3]!=70 && ReactionT[nn][3]!=105) {
						for (j=1; j<=zbin; j++) {
							deltaG = 0.0;
							deltaG -= GibbsForm[ReactionT[nn][1]][j];
							deltaG += GibbsForm[ReactionT[nn][2]][j];
							if (np>1) {
								deltaG += GibbsForm[ReactionT[nn][3]][j];
							}
							RkkT[j][nn] = kkT[j][nn]*exp(deltaG/R_GAS/tl[j])*pow(kbmodi*tl[j],np-nr);
						}
					}
			}
			
		}
	}
	
	
	
}