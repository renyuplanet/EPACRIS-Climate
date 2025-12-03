#include <math.h>
#include "nrutil.h"
void ms_LUdcmp(double **lu, int n, int *indx);

void ms_LUdcmp(double **lu, int n, int *indx)
{
    const double TINY=1e-40;
    int i,imax,j,k,d;
    double big,temp;
    double vv[n];

    d=1;
    for (i=0;i<n;i++) 
    {
	big=0.0;
	for (j=0;j<n;j++) if ((temp=fabs(lu[i][j])) > big) big=temp;
	if (big == 0.0) nrerror("Singular matrix in ms_LUdcmp");
	vv[i]=1.0/big;
    }
    for (k=0;k<n;k++) 
    {
	big=0.0;
	for (i=k;i<n;i++) 
        {
	    temp=vv[i]*fabs(lu[i][k]);
            if (temp > big) 
            {
		big=temp;
		imax=i;
	    }
	}
	if (k != imax) {
	    for (j=0;j<n;j++) {
		temp=lu[imax][k];
		lu[imax][j]=lu[k][j];
		lu[k][j]=temp;
	    }
	    d = -d;
	    vv[imax]=vv[k];
	}
	indx[k]=imax;
	if (lu[k][k] == 0.0) lu[k][k]=TINY;
	for (i=k+1;i<n;i++) 
        {
            temp = lu[i][k] / lu[k][k];
            for (j=k+1;j<n;j++) lu[i][j] -= temp*lu[k][j];
	}
    }
}
