#include <math.h>
#include "constant.h"

void ChemEqu(double x[], double xc[], double xf[], double fvec[], double **Jvalue, int lx[], int lc[], int lf[], int laer[],
			   int zone_r[], int zone_m[], int zone_t[], double **JJ, int zone_p[], double Upflux[], double Loflux[], double Depo[], int listFix[])
{
	int i, j, jj, temp, tempc, tempf, tempx, std;
	double rate, rate1, rate2, rater, rater1, rater2, rater3;
	int ntotal=zbin*NSP;
	int nn=zbin*numx; /* dimension of Jvalue */
	double **zz, **P, **L;
	P=dmatrix(1, zbin, 1, NSP);
	L=dmatrix(1, zbin, 1, NSP);
	/* Initialization of matrix */
	for (i=1; i<=zbin; i++) {
		for (j=1; j<=NSP; j++) {
			P[i][j]=0;
			L[i][j]=0;
		}
	}
	zz=dmatrix(1, ntotal, 1, ntotal);
	for (i=1; i<=ntotal; i++) {
		for (j=1; j<=ntotal; j++) {
			zz[i][j]=0;
		}
	}
	
	/**************************************************************************************/
	/* Calculate Jvalue */
	
	/* clean-up Jvalue */
	for (i=1; i<=nn; i++) {
		for (j=1; j<=nn; j++) {
			Jvalue[i][j]=0;
		}
	}
	
	/*Difussion Part*/
	for (j=1; j<=numx; j++) {
		Jvalue[j][j] = -(KE[1]+DM[1][j])/MM[1];
		Jvalue[j][j+numx] = (KE[1]+DM[1][j])/MM[2];
	for (i=2; i<zbin; i++) {
		temp=numx*(i-1);
		Jvalue[temp+j][temp+j+numx] = (KE[i]+DM[i][j])/MM[i+1];
		Jvalue[temp+j][temp+j] = -((KE[i]+DM[i][j])/MM[i]+(KE[i-1]+DM[i-1][j])/MM[i]);
		Jvalue[temp+j][temp+j-numx] = (KE[i-1]+DM[i-1][j])/MM[i-1];
	}
	temp=numx*(zbin-1);
	i=zbin;
	Jvalue[temp+j][temp+j] = -(KE[i-1]+DM[i-1][j])/MM[i];
	Jvalue[temp+j][temp+j-numx] = (KE[i-1]+DM[i-1][j])/MM[i-1];
	}
	for (j=1; j<=numx; j++) {
		Jvalue[j][j] -= dl[1][j]/MM[1];
		Jvalue[j][j+numx] -= dl[1][j]/MM[2];
		for (i=2; i<zbin; i++) {
			temp=numx*(i-1);
			Jvalue[temp+j][temp+j+numx] -= dl[i][j]/MM[i+1];
			Jvalue[temp+j][temp+j] -= (dl[i][j]/MM[i]-dl[i-1][j]/MM[i]);
			Jvalue[temp+j][temp+j-numx] += dl[i-1][j]/MM[i-1];
		}
		temp=numx*(zbin-1);
		i=zbin;
		Jvalue[temp+j][temp+j] += dl[i-1][j]/MM[i];
		Jvalue[temp+j][temp+j-numx] += dl[i-1][j]/MM[i-1];
	}
	/* Reset the diffusion terms for aerosol particles */
	/* if (numa>0) {
		for (j=1; j<=numa; j++) {
			for (i=1; i<=zbin; i++) {
				temp=numx*(i-1);
				for (jj=1; jj<=nn; jj++) {
					Jvalue[temp+laer[j]][jj] = 0.0;
				}
			}
		}
	} */
	
	/* Aerosols */
	if (numa>0) {
		for (j=1; j<=numa; j++) {
			/* Falling of aerosols */
			for (i=1; i<zbin; i++) {
				temp=numx*(i-1)+laer[j];
				Jvalue[temp][temp+numx] += VFall[i+1]/thickl;
				Jvalue[temp][temp] -= VFall[i]/thickl;
			}
			i=zbin;
			temp=numx*(i-1)+laer[j];
			Jvalue[temp][temp] -= VFall[i]/thickl;
			/* Condensation/Evaporation of H2SO4 */
			if (lx[laer[j]] == 78) {
				for (jj=1; jj<=numx; jj++) {
					if (lx[jj] == 73) {
						for (i=1; i<=zbin; i++) {
							if (xx[i][73] > nsH2SO4[i]) { /* condensation */
								temp=numx*(i-1);
								Jvalue[temp+laer[j]][temp+jj] += tcondfH2SO4[i]*(2*xx[i][73]-nsH2SO4[i]);
								Jvalue[temp+jj][temp+jj]      -= tcondfH2SO4[i]*(2*xx[i][73]-nsH2SO4[i]);
							} else { /* evaporation */
								temp=numx*(i-1);
								Jvalue[temp+laer[j]][temp+laer[j]] -= tcondfH2SO4[i]*(nsH2SO4[i]-xx[i][73]);
								Jvalue[temp+laer[j]][temp+jj]      += tcondfH2SO4[i]*xx[i][78];
								Jvalue[temp+jj][temp+laer[j]]      += tcondfH2SO4[i]*(nsH2SO4[i]-xx[i][73]);
								Jvalue[temp+jj][temp+jj]           -= tcondfH2SO4[i]*xx[i][78];
							}
						}
					}
				}
				
			}
			/* Condensation/Evaporation of S8 */
			if (lx[laer[j]] == 111) {
				for (jj=1; jj<=numx; jj++) {
					if (lx[jj] == 79) {
						for (i=1; i<=zbin; i++) {
							if (xx[i][79] > nsS8[i]) { /* condensation */
								temp=numx*(i-1);
								Jvalue[temp+laer[j]][temp+jj] += tcondfS8[i]*(2*xx[i][79]-nsS8[i]);
								Jvalue[temp+jj][temp+jj]      -= tcondfS8[i]*(2*xx[i][79]-nsS8[i]);
							} else { /* evaporation */
								temp=numx*(i-1);
								Jvalue[temp+laer[j]][temp+laer[j]] -= tcondfS8[i]*(nsS8[i]-xx[i][79]);
								Jvalue[temp+laer[j]][temp+jj]      += tcondfS8[i]*xx[i][111];
								Jvalue[temp+jj][temp+laer[j]]      += tcondfS8[i]*(nsS8[i]-xx[i][79]);
								Jvalue[temp+jj][temp+jj]           -= tcondfS8[i]*xx[i][111];
							}
						}
					}
				}
				
			}
		}
	}
	
	/* Condensation of H2O */
	for (i=1; i<=numx; i++) {
		if (lx[i] == 7) {
			for (j=1; j<=zbin; j++) {
				if (xx[j][7] > nsH2O[j]) {
					temp = numx*(j-1);
					Jvalue[temp+i][temp+i] -= tcondfH2O[j]*(2*xx[j][7]-nsH2O[j]);
				}
			}
		}
	}
	
	/* Chemistry Part*/
	for (i=1; i<=numr; i++) {
        std=zone_r[i]; /* Get the standard reaction number from list */
	    for (j=1; j<=zbin; j++) {
			temp=NSP*(j-1);
			rate1=kk[j][std]*xx[j][ReactionR[std][2]];
			rate2=kk[j][std]*xx[j][ReactionR[std][1]];
            zz[temp+ReactionR[std][1]][temp+ReactionR[std][1]] -= rate1;
			zz[temp+ReactionR[std][1]][temp+ReactionR[std][2]] -= rate2;
			zz[temp+ReactionR[std][2]][temp+ReactionR[std][1]] -= rate1;
			zz[temp+ReactionR[std][2]][temp+ReactionR[std][2]] -= rate2;
			zz[temp+ReactionR[std][4]][temp+ReactionR[std][1]] += rate1;
			zz[temp+ReactionR[std][4]][temp+ReactionR[std][2]] += rate2;
			if (ReactionR[std][3]>1) {
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][1]] += rate1;
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][2]] += rate2;
			}
			if (ReactionR[std][3]>2) {
				zz[temp+ReactionR[std][6]][temp+ReactionR[std][1]] += rate1;
				zz[temp+ReactionR[std][6]][temp+ReactionR[std][2]] += rate2;
			}
			if (ReactionR[std][3]==1) {
				rater=Rkk[j][std];
				zz[temp+ReactionR[std][1]][temp+ReactionR[std][4]] += rater;
				zz[temp+ReactionR[std][2]][temp+ReactionR[std][4]] += rater;
				zz[temp+ReactionR[std][4]][temp+ReactionR[std][4]] -= rater;
			}
			if (ReactionR[std][3]==2) {
				rater1=Rkk[j][std]*xx[j][ReactionR[std][5]];
				rater2=Rkk[j][std]*xx[j][ReactionR[std][4]];
				zz[temp+ReactionR[std][1]][temp+ReactionR[std][4]] += rater1;
				zz[temp+ReactionR[std][1]][temp+ReactionR[std][5]] += rater2;
				zz[temp+ReactionR[std][2]][temp+ReactionR[std][4]] += rater1;
				zz[temp+ReactionR[std][2]][temp+ReactionR[std][5]] += rater2;
				zz[temp+ReactionR[std][4]][temp+ReactionR[std][4]] -= rater1;
				zz[temp+ReactionR[std][4]][temp+ReactionR[std][5]] -= rater2;
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][4]] -= rater1;
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][5]] -= rater2;
			}
			if (ReactionR[std][3]==3) {
				rater1=Rkk[j][std]*xx[j][ReactionR[std][5]]*xx[j][ReactionR[std][6]];
				rater2=Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][6]];
				rater3=Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][5]];
				zz[temp+ReactionR[std][1]][temp+ReactionR[std][4]] += rater1;
				zz[temp+ReactionR[std][1]][temp+ReactionR[std][5]] += rater2;
				zz[temp+ReactionR[std][1]][temp+ReactionR[std][6]] += rater3;
				zz[temp+ReactionR[std][2]][temp+ReactionR[std][4]] += rater1;
				zz[temp+ReactionR[std][2]][temp+ReactionR[std][5]] += rater2;
				zz[temp+ReactionR[std][2]][temp+ReactionR[std][6]] += rater3;
				zz[temp+ReactionR[std][4]][temp+ReactionR[std][4]] -= rater1;
				zz[temp+ReactionR[std][4]][temp+ReactionR[std][5]] -= rater2;
				zz[temp+ReactionR[std][4]][temp+ReactionR[std][6]] -= rater3;
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][4]] -= rater1;
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][5]] -= rater2;
				zz[temp+ReactionR[std][5]][temp+ReactionR[std][6]] -= rater3;
				zz[temp+ReactionR[std][6]][temp+ReactionR[std][4]] -= rater1;
				zz[temp+ReactionR[std][6]][temp+ReactionR[std][5]] -= rater2;
				zz[temp+ReactionR[std][6]][temp+ReactionR[std][6]] -= rater3;
			}
		}
	}
	
	for (i=1; i<=numm; i++) {
        std=zone_m[i];
		for (j=1; j<=zbin; j++) {
			temp=NSP*(j-1);
			rate1=kkM[j][std]*xx[j][ReactionM[std][2]];
			rate2=kkM[j][std]*xx[j][ReactionM[std][1]];
            zz[temp+ReactionM[std][1]][temp+ReactionM[std][1]] -= rate1;
			zz[temp+ReactionM[std][1]][temp+ReactionM[std][2]] -= rate2;
			zz[temp+ReactionM[std][2]][temp+ReactionM[std][1]] -= rate1;
			zz[temp+ReactionM[std][2]][temp+ReactionM[std][2]] -= rate2;
			zz[temp+ReactionM[std][3]][temp+ReactionM[std][1]] += rate1;
			zz[temp+ReactionM[std][3]][temp+ReactionM[std][2]] += rate2;
		    if (ReactionM[std][4]>0) {
				zz[temp+ReactionM[std][4]][temp+ReactionM[std][1]] += rate1;
				zz[temp+ReactionM[std][4]][temp+ReactionM[std][2]] += rate2;
		    }
			if (ReactionM[std][4]==0) {
				rater=RkkM[j][std];
				zz[temp+ReactionM[std][1]][temp+ReactionM[std][3]] += rater;
				zz[temp+ReactionM[std][2]][temp+ReactionM[std][3]] += rater;
				zz[temp+ReactionM[std][3]][temp+ReactionM[std][3]] -= rater;
			}
			if (ReactionM[std][4]>0) {
				rater1=RkkM[j][std]*xx[j][ReactionM[std][4]];
				rater2=RkkM[j][std]*xx[j][ReactionM[std][3]];
				zz[temp+ReactionM[std][1]][temp+ReactionM[std][3]] += rater1;
				zz[temp+ReactionM[std][1]][temp+ReactionM[std][4]] += rater2;
				zz[temp+ReactionM[std][2]][temp+ReactionM[std][3]] += rater1;
				zz[temp+ReactionM[std][2]][temp+ReactionM[std][4]] += rater2;
				zz[temp+ReactionM[std][3]][temp+ReactionM[std][3]] -= rater1;
				zz[temp+ReactionM[std][3]][temp+ReactionM[std][4]] -= rater2;
				zz[temp+ReactionM[std][4]][temp+ReactionM[std][3]] -= rater1;
				zz[temp+ReactionM[std][4]][temp+ReactionM[std][4]] -= rater2;
			}
		}
	}
	
	for (i=1; i<=nump; i++) {
        std=zone_p[i];
		for (j=1; j<=zbin; j++) {
			temp=NSP*(j-1);
			rate=JJ[j][i];
			zz[temp+ReactionP[std][1]][temp+ReactionP[std][1]] -= rate;
			zz[temp+ReactionP[std][2]][temp+ReactionP[std][1]] += rate;
			if (ReactionP[std][3]>0) {
				zz[temp+ReactionP[std][3]][temp+ReactionP[std][1]] += rate;
			}
			if (ReactionP[std][4]>0) {
				zz[temp+ReactionP[std][4]][temp+ReactionP[std][1]] += rate;
			}
		}
	}
	
	 for (i=1; i<=numt; i++) {
        std=zone_t[i];
		for (j=1; j<=zbin; j++) {
			temp=NSP*(j-1);
			rate=kkT[j][std];
			zz[temp+ReactionT[std][1]][temp+ReactionT[std][1]] -= rate;
			zz[temp+ReactionT[std][2]][temp+ReactionT[std][1]] += rate;
			if (ReactionT[std][3]>0) {
				zz[temp+ReactionT[std][3]][temp+ReactionT[std][1]] += rate;
			}
			if (ReactionT[std][3]==0) {
				rater=RkkT[j][std];
				zz[temp+ReactionT[std][1]][temp+ReactionT[std][2]] += rater;
				zz[temp+ReactionT[std][2]][temp+ReactionT[std][2]] -= rater;
			}
			if (ReactionT[std][3]>0) {
				rater1=RkkT[j][std]*xx[j][ReactionT[std][3]];
				rater2=RkkT[j][std]*xx[j][ReactionT[std][2]];
				zz[temp+ReactionT[std][1]][temp+ReactionT[std][2]] += rater1;
				zz[temp+ReactionT[std][1]][temp+ReactionT[std][3]] += rater2;
				zz[temp+ReactionT[std][2]][temp+ReactionT[std][2]] -= rater1;
				zz[temp+ReactionT[std][2]][temp+ReactionT[std][3]] -= rater2;
				zz[temp+ReactionT[std][3]][temp+ReactionT[std][2]] -= rater1;
				zz[temp+ReactionT[std][3]][temp+ReactionT[std][3]] -= rater2;
			}
		}
	}
	
	for (i=1; i<=zbin; i++) {
		temp=NSP*(i-1);
		tempx=numx*(i-1);
		for (j=1; j<=numx; j++) {
			for (jj=1; jj<=numx; jj++) {
				Jvalue[tempx+j][tempx+jj] += zz[temp+lx[j]][temp+lx[jj]];
			}
			/* rainout */
			if (lx[j]!=7) {
				Jvalue[tempx+j][tempx+j] -= rainoutrate[i][lx[j]]*xx[i][7];
				if (waterx==1) {
					Jvalue[tempx+j][tempx+waternum] -= rainoutrate[i][lx[j]]*xx[i][lx[j]];
				}
			}
		}
	}
	/* Diffusion-limited Escape */
	i=zbin;
	tempx=numx*(i-1);
	for (j=1; j<=numx; j++) {
		Jvalue[tempx+j][tempx+j] -= Vesc[lx[j]];
	}
	/* Dry Deposition */
	i=1;
	for (j=1; j<=numx; j++) {
		Jvalue[j][j] -= Depo[j];
	}
	/* Fix species at the lower boundary */
	for (j=1; j<=numx; j++) {
		if (listFix[j] == 1) {
			for (jj=1; jj<=numx+numx; jj++) {
				Jvalue[j][jj] = 0.0;
				Jvalue[jj][j] = 0.0;
			}
		}
	}
	
	/**************************************************************************************/
	/* Calculate fvec */
	
	/* Diffusion terms */
	for (j=1; j<=numx; j++) {
		/* Lower Boundary Condition */
		fvec[j]=x[j+numx]*(KE[1]+DM[1][j])/MM[2] - x[j]*(KE[1]+DM[1][j])/MM[1];
		/* Diffusion-Continuity Equatipn */
		for (i=2; i<zbin; i++) {
			temp=numx*(i-1);
			fvec[temp+j]=x[temp+j+numx]*(KE[i]+DM[i][j])/MM[i+1] - x[temp+j]*((KE[i]+DM[i][j])/MM[i]+(KE[i-1]+DM[i-1][j])/MM[i]) + x[temp+j-numx]*(KE[i-1]+DM[i-1][j])/MM[i-1];
		}
		/* Upper Boundary Condition */
		i=zbin;
		temp=numx*(zbin-1);
		fvec[temp+j]=-x[temp+j]*(KE[i-1]+DM[i-1][j])/MM[i] + x[temp+j-numx]*(KE[i-1]+DM[i-1][j])/MM[i-1];
	} 
	for (j=1; j<=numx; j++) {
		fvec[j] +=  (-x[j+numx]*dl[1][j]/MM[2] - x[j]*dl[1][j]/MM[1]);
		for (i=2; i<zbin; i++) {
			temp=numx*(i-1);
			fvec[temp+j] += (-x[temp+j+numx]*dl[i][j]/MM[i+1] - x[temp+j]*(dl[i][j]/MM[i]-dl[i-1][j]/MM[i]) + x[temp+j-numx]*dl[i-1][j]/MM[i-1]);
		}
		i=zbin;
		temp=numx*(zbin-1);
		fvec[temp+j] += (x[temp+j]*dl[i-1][j]/MM[i] + x[temp+j-numx]*dl[i-1][j]/MM[i-1]);
	}
	
	/* Reset the diffusion terms for aerosol particles */
	/* if (numa>0) {
		for (j=1; j<=numa; j++) {
			for (i=1; i<=zbin; i++) {
				temp=numx*(i-1);
				fvec[temp+laer[j]] = 0.0;
			}
		}
	} */
	
	if (numa>0) {
	/* Aerosols */
	for (j=1; j<=numa; j++) {
		/* Falling of Aerosols */
		for (i=1; i<zbin; i++) {
			temp=numx*(i-1)+laer[j];
			fvec[temp] -= VFall[i]*x[temp]/thickl;
			fvec[temp] += VFall[i+1]*x[temp+numx]/thickl;
		}
		i=zbin;
		temp=numx*(i-1)+laer[j];
		fvec[temp] -= VFall[i]*x[temp]/thickl;
		/* Condensation/evaporation of H2SO4 */
		if (lx[laer[j]] == 78) {
			for (jj=1; jj<=numx; jj++) {
				if (lx[jj] == 73) {
						for (i=1; i<=zbin; i++) {
							if (xx[i][73] > nsH2SO4[i]) {  /* condensation */
								temp=numx*(i-1);
								fvec[temp+laer[j]] += tcondfH2SO4[i]*(xx[i][73]-nsH2SO4[i])*xx[i][73];
								fvec[temp+jj]  -= tcondfH2SO4[i]*(xx[i][73]-nsH2SO4[i])*xx[i][73];
							} else { /* evaporation */
								temp=numx*(i-1);
								fvec[temp+laer[j]] -= tcondfH2SO4[i]*(nsH2SO4[i]-xx[i][73])*xx[i][78];
								fvec[temp+jj]  += tcondfH2SO4[i]*(nsH2SO4[i]-xx[i][73])*xx[i][78];
							}
						}
				}
			}

		}
		/* Condensation/evaporation of S8 */
		if (lx[laer[j]] == 111) {
			for (jj=1; jj<=numx; jj++) {
				if (lx[jj] == 79) {
					for (i=1; i<=zbin; i++) {
						if (xx[i][79] > nsS8[i]) {  /* condensation */
							temp=numx*(i-1);
							fvec[temp+laer[j]] += tcondfS8[i]*(xx[i][79]-nsS8[i])*xx[i][79];
							fvec[temp+jj]  -= tcondfS8[i]*(xx[i][79]-nsS8[i])*xx[i][79];
						} else { /* evaporation */
							temp=numx*(i-1);
							fvec[temp+laer[j]] -= tcondfS8[i]*(nsS8[i]-xx[i][79])*xx[i][111];
							fvec[temp+jj]  += tcondfS8[i]*(nsS8[i]-xx[i][79])*xx[i][111];
						}
					}
				}
			}
			
		}
	}
	}
	
	/* Condensation of H2O */
	for (i=1; i<=numx; i++) {
		if (lx[i] == 7) {
			for (j=1; j<=zbin; j++) {
				if (xx[j][7] > nsH2O[j]) {
					temp = numx*(j-1);
					fvec[temp+i] -= tcondfH2O[j]*(xx[j][7]-nsH2O[j])*xx[j][7];
				}
			}
		}
	}
		
	/*Chemical Terms*/	
	for (i=1; i<=numr; i++) {
        std=zone_r[i];
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
				rater=RkkM[j][std]*xx[j][ReactionM[std][3]];
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
	
	/* for (i=1; i<=zbin; i++) {
		printf("%d %d %d\n", numx, numf, numa);
		printf("%s %f %f %f\n", "chemical flux at altitude", zl[i], P[i][43], L[i][43]*xx[i][43]);
	} */
	
	/* check sulfur chemistry */
	/* for (i=1; i<=zbin; i++) {
	printf("%s %f %e\n", "total sulfur chem production at altitude", zl[i], 
		   P[i][40]+2*P[i][41]+P[i][42]+P[i][43]+P[i][44]+P[i][45]+P[i][46]+P[i][47]+P[i][48]+P[i][49]
		   +P[i][50]+P[i][51]+P[i][68]+P[i][69]+P[i][73]+P[i][74]+P[i][75]+3*P[i][76]+4*P[i][77]+P[i][78]+8*P[i][79]);
	printf("%s %f %e\n", "total sulfur chem loss at altitude", zl[i], 
		   L[i][40]*xx[i][40]+2*L[i][41]*xx[i][41]+L[i][42]*xx[i][42]+L[i][43]*xx[i][43]+L[i][44]*xx[i][44]
		   +L[i][45]*xx[i][45]+L[i][46]*xx[i][46]+L[i][47]*xx[i][47]+L[i][48]*xx[i][48]+L[i][49]*xx[i][49]
		   +L[i][50]*xx[i][50]+L[i][51]*xx[i][51]+L[i][68]*xx[i][68]+L[i][69]*xx[i][69]+L[i][73]*xx[i][73]
		   +L[i][74]*xx[i][74]+L[i][75]*xx[i][75]+3*L[i][76]*xx[i][76]+4*L[i][77]*xx[i][77]+L[i][78]*xx[i][78]+8*L[i][79]*xx[i][79]);
	} */
	
	for (i=1; i<=zbin; i++) {
		temp=numx*(i-1);
		for (j=1; j<=numx; j++) {
			fvec[temp+j] += P[i][lx[j]];
			fvec[temp+j] -= L[i][lx[j]]*xx[i][lx[j]];
			/* rainout */
			if (lx[j]!=7) {
				fvec[temp+j] -= xx[i][lx[j]]*rainoutrate[i][lx[j]]*xx[i][7];
			}
		}
	}

	/* Diffusion-limited Escape & Additional Upper Boundary Flux */
	i=zbin;
	temp=numx*(i-1);
	for (j=1; j<=numx; j++) {
		fvec[temp+j] -= xx[i][lx[j]]*Vesc[lx[j]];
		fvec[temp+j] -= Upflux[j];
	}
	/* Dry Deposition & Additional Lower Boundary Flux */
	i=1;
	for (j=1; j<=numx; j++) {
		fvec[j] -= xx[1][lx[j]]*Depo[j];
		fvec[j] += Loflux[j];
	}
	/* PhotoChemical Equilibrium Model for Fast */
	for (i=1; i<=zbin; i++) {
		temp=numf*(i-1);
		for (j=1; j<=numf; j++) {
			xf[temp+j]=P[i][lf[j]]/(L[i][lf[j]]+1.0E-24); /* Photochemical equilibrium*/
		}
	}
	/* Fix species at the lower boundary */
	for (j=1; j<=numx; j++) {
		if (listFix[j] == 1) {
			fvec[j] = 0.0;
		}
	}

	free_dmatrix(zz, 1, ntotal, 1, ntotal);
	free_dmatrix(P, 1, zbin, 1, NSP);
	free_dmatrix(L, 1, zbin, 1, NSP);
}

