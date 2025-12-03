/*Function to calculate water saturation pressure from temperature*/

/* Formulation from Murphy & Koop (2005)*/
/* SI */


#include <math.h>

void waterpressure(double P[]);

void waterpressure(double P[])
{
	int i;
	double temp, a;

	for (i=1; i<=zbin; i++) {
		if (tl[i]<273.16) {
			/* over ice */
			/* Formulation from Murphy & Koop (2005)*/
			P[i] = exp(9.550426-5723.265/tl[i]+3.53068*log(tl[i])-0.00728332*tl[i]); /* Pa */
		}else {
			/* over water */
			/* Formulation from Seinfeld & Pandis (2006) */
			a = 1-373.15/tl[i];
			P[i] = 101325*exp(13.3185*a-1.97*a*a-0.6445*a*a*a-0.1229*a*a*a*a); /* Pa */
		}
	}
	
}

void sulfuridpressure(double P[]);

void sulfuridpressure(double P[])
{
	int i;
	double temp;
	
	for (i=1; i<=zbin; i++) {
		/* Expression of Seinfeld and Pandis (2006), Ayers et al. (1980) measurement and corrected by Kulmala & Laaksonen (1990), H2SO4 weight 98.01% */
		temp = -11.9521+10156*(-1.0/tl[i]+1.0/360.0+0.38/545.0*(1+log(360.0/tl[i])-360.0/tl[i]));
		/* correct for H2SO4 weight 75%, using chemical potential data from Gianque et al. (1960) at 25 oC */
		temp=temp - 3567.0/1.98726/tl[i];
		P[i]=exp(temp)*101325; /* Pa */
	}
}

void sulfurpressure(double P[]);

void sulfurpressure(double P[])
{
	int i;
	double temp;
	
	for (i=1; i<=zbin; i++) {
		if (tl[i]<=392.0) { /* total sulfur over monoclinic sulfur by Lyons (2008) */
			temp = 8.4832 - 5082.0/tl[i];
			P[i] = pow(10,temp)*101325; /* Pa */
		} else { /* total sulfur over liquid sulfur by Lyons (2008) */
			temp = 7.0240 - 6091.2/tl[i];
			P[i] =pow(10,temp)*101325; /* Pa */
			temp = 6.3428 - 6202.2/tl[i];
			P[i] += pow(10,temp)*101325; /* Pa */
			temp = 6.0028 - 6047.5/tl[i];
			P[i] += pow(10,temp)*101325; /* Pa */
			temp = 5.1609 - 4714.8/tl[i];
			P[i] += pow(10,temp)*101325; /* Pa */
			temp = 4.8039 - 3814.1/tl[i];
			P[i] += pow(10,temp)*101325; /* Pa */
			temp = 5.2127 - 4113.6/tl[i];
			P[i] += pow(10,temp)*101325; /* Pa */
			temp = 4.1879 - 3269.1/tl[i];
			P[i] += pow(10,temp)*101325; /* Pa */
		}

	}
}

