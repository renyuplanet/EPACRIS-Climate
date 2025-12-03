#include <math.h>
#include "Trapz.c"

void Photorate(double **Radiation, double **cross, double **crosst, double **qy, double **qyt, int zone_p[], int nump, double **JJ)
{
     
	int imin, imax, i, j, k, zz, std, iref;
	double *x1, *y1, JJ1[zbin+1];
	double crossl, qyl; /* local cross and qy */
	
	/* index for reference wavelength of 300 nm */
    iref=0;
    for (i=0; i<NLAMBDA; i++) {
        if (wavelength[i]<=300.0 && wavelength[i]<=WaveMax) {iref=i;}
    }
    /* printf("%s %f %d\n", "photorate", wavelength[iref], iref); */
	
	/* process temperature for cross sections */
	double temperature[zbin+1];
	temperature[0] = tl[1];
	temperature[zbin] = tl[zbin];
	for (zz=1; zz<zbin; zz++) {
		temperature[zz] = (tl[zz] + tl[zz+1])/2.0;
	}
	for (zz=0; zz<=zbin; zz++) {
		if (temperature[zz] > TDEPMAX) {
			temperature[zz] = TDEPMAX; 
		}
		if (temperature[zz] < TDEPMIN) {
			temperature[zz] = TDEPMIN; 
		}
	}
	
	for (k=1;k<=nump;k++)
	{
         std=zone_p[k];
         if (std!= 52 ) {
             imin=0;
             imax=0;
             for (i=0; i<NLAMBDA; i++) {
                 if (wavelength[i]<=ReactionP[std][5] && wavelength[i]<=WaveMax) {imin=i;}
                 if (wavelength[i]<=ReactionP[std][6] && wavelength[i]<=WaveMax) {imax=i;}
             }
             /* printf("%s %d %d %f %f\n", "photorate", imin, imax, wavelength[imin], wavelength[imax]); */
         x1=dvector(0,imax-imin-1);
         y1=dvector(0,imax-imin-1);
       	 for(zz=0;zz<=zbin;zz++) {
           i=imin;
	       j=0;
           while (i<imax)
		   {
			   x1[j]=wavelength[i];
			   crossl = cross[k][i] + crosst[k][i] * ( temperature[zz] - 295.0 ) ;
			   qyl = qy[k][i] + qyt[k][i] * ( temperature[zz] - 295.0 ) ;
			   y1[j]=crossl*qyl*Radiation[i][zz]*(wavelength[i]*pow(10, -13))/HPLANCK/CLIGHT;
			   i++;
			   j++;
			}
           if (TIDELOCK!=1) {
			   JJ1[zz]=Trapz(x1, y1, imax-imin)/2.0;  /* Diurnal variation */
		   } else {
			   JJ1[zz]=Trapz(x1, y1, imax-imin);
		   }

         }
         free_dvector(x1,0,imax-imin-1);
         free_dvector(y1,0,imax-imin-1);
         }

		 if (std == 52) {
			 for (zz=1; zz<=zbin; zz++) {
				 JJ1[zz]=4.0E-3/solar[iref]*Radiation[iref][zz]/ORBIT/ORBIT; } /* Special Treatment for S2 */
         }
		
		for (zz=1; zz<=zbin; zz++) {
			JJ[zz][k] = (JJ1[zz-1] + JJ1[zz])/2.0;
		}
		
		/* printf("%s %d %e\n", "P", k, JJ1[zbin]); */
		
		
	}

}
