#include <math.h>
#include "constant.h"

/* Calculating the corresponding z profile from T-P Profile for Hot Jupiter
 
 Output: z: altitude in km
 Input:  P: Pressure profile in SI
         T: Temperature profile
         NL: number of layer */

double TPScale(double P[], double T[], int NL, double z[])
{
	int i, j;
	double GA;
	
	GA=GRAVITY*MASS_PLANET/RADIUS_PLANET/RADIUS_PLANET; /* Planet Surface Gravity Acceleration, in SI */
	z[NL-1]=0;
	for (i=NL-2; i>=0; i--) {
		z[i]=z[i+1]+(log(P[i+1])-log(P[i]))*KBOLTZMANN*(T[i+1]+T[i])/2/AIRM/AMU/GA/1000;
	}
	
	Reverse(P,zbin+1);
	Reverse(T,zbin+1);
	Reverse(z,zbin+1);
	
	double pn[zbin+1], tn[zbin+1], zn[zbin+1], zmax1, zmin1, zd;
	
	zmax1 = z[1];
	zmin1 = z[1];
	for (i=0; i<=zbin; i++) {
		if (z[i]>zmax1) {
			zmax1 = z[i];
		}
		if (z[i]<zmin1) {
			zmin1 = z[i];
		}
	}
	zd = (zmax1-zmin1)/zbin;
	for (j=0; j<=zbin; j++) { zn[j] = zmin1+j*zd; }
	Interpolation(zn, zbin+1, pn, z, P, NL, 0);
	Interpolation(zn, zbin+1, tn, z, T, NL, 0);
	
	for (j=0; j<=zbin; j++) {
		z[j] = zn[j];
		P[j] = pn[j];
		T[j] = tn[j];
	}
	
	return zd;
	
}
