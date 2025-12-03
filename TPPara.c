#include <math.h>

/* Calculating T-P Profile for Hot Jupiter
 
 Output: P: Pressure profile in SI
         T: Temperature profile
 Input:  NL: number of layer
         P0: Pressure at the top of atmosphere in bar
		 T0: Temperature at P0
         P1: Pressure at the top of troposphere or stratosphere (if inversion)
         T1: Temperature at P1
         TINV: 1 for inversion, 0 for no inversion
         P2: Pressure at the top of troposphere (only used if inversion)
         T2: Temperature at T2
         P3: Pressure at the bottom of troposphere
         T3: Temperature at T3
         P4: Pressure at the bottom of atmosphere */

void TPPara(double P[], double T[], int Tinv, int NL, 
			double P0, double T0, double P1, double T1, double P2, double T2, double P3, double T3, double P4)
{
	double grid[NL];
	double bar, gridi, a1, a2, a3;
	int i;
	
	bar=pow(10,5);
	grid[0]=log(P0);
	T[0]=T0;
	P[0]=P0*bar;
	gridi=(log(P4)-grid[0])/(NL-1);
	if (Tinv==1) {
		i=1;
		a1=(log(P1)-grid[0])/sqrt(T1-T0);
		while (P[i-1]*exp(gridi)<P1*bar) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-grid[0])/a1,2)+T0;
			i++;
		}
		a2=(log(P1)-log(P2))/sqrt(T1-T2);
		while (P[i-1]*exp(gridi)<P2*bar) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-log(P2))/a2,2)+T2;
			i++;
		}
        a3=(log(P3)-log(P2))/sqrt(T3-T2);
        while (P[i-1]*exp(gridi)<P3*bar) {
            grid[i]=grid[i-1]+gridi;
            P[i]=exp(grid[i])*bar;
            T[i]=pow((grid[i]-log(P2))/a3,2)+T2;
            i++;
        }
        while (i<NL) {
            grid[i]=grid[i-1]+gridi;
            P[i]=exp(grid[i])*bar;
            T[i]=T3;
            i++;
        }
	} else {
		i=1;
		a1=(log(P1)-grid[0])/sqrt(T1-T0);
		while (P[i-1]*exp(gridi)<P1*bar) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-grid[0])/a1,2)+T0;
			i++;
		}
		a2=(log(P3)-log(P1))/sqrt(T3-T1);
		while (P[i-1]*exp(gridi)<P3*bar) {
			grid[i]=grid[i-1]+gridi;
			P[i]=exp(grid[i])*bar;
			T[i]=pow((grid[i]-log(P1))/a2,2)+T1;
			i++;
		}
        while (i<NL) {
            grid[i]=grid[i-1]+gridi;
            P[i]=exp(grid[i])*bar;
            T[i]=T3;
            i++;
        }
	}

}
