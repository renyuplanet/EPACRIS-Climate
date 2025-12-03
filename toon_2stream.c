/*******************************************************************
 * File: toon_2stream.c
 * Function: Solve for equilibrium 2-stream fluxes using Toon+89 Eqs. 
 * Original Author: Renyu Hu (renyu.hu@jpl.nasa.gov)
 * Version: v01, 2022: Original upload by Markus Scheucher
 *******************************************************************
*/                  

#include <math.h>
#include "constant.h"

//=========================================================
//=== Function identifiers=================================
void toon_two_str_solver(int lbin, double w[], double g[], double tau[], int nrl, int isconv[], double **Tvar, int jTvar, double P[], double **NetFlux, double Fup[], double Fdn[], double Fcup[], double Fcdn[], double Tint);

//=========================================================
//=== Functions ===========================================
void toon_two_str_solver(int lbin, double w[], double g[], double tau[], int nrl, int isconv[], double **Tvar, int jTvar, double P[], double **NetFlux, double Fup[], double Fdn[], double Fcup[], double Fcdn[], double Tint)
{
/* This function solves for Fnet using inversion of a tri-diagonal matrix. It includes a loop (tv) for temperature variation to populate a jacobian matrix later, to solve for equilibrium.
 *Equations setup for tri-diagonal matrix (see Toon+89): X(i) F(i-1) + Y(i) F(i) + Z(i) F(i+1) = S(i)
 *
 */ 
//=========================================================
//printf("%s\n", "START TOON 2 STREAM"); //template
//=== DECLARE variables ===================================

    int i,j,jj,k,kk,l,radcount;
    double tc[zbin+1], gc[zbin+1], wc[zbin+1];
    double *gamma1, *gamma2, *gamma3, *gamma4, *lambda, *gamma, *e1, *e2, *e3, *e4, *CP0, *CP1, *CM0, *CM1, *CPh, *CMh;
    double B0, B1, miu1=0.5;
    int    ntau, nvector;
    double *A, *AS, *B, *D, *DS, *E, *ES, *Y, *Y1, *Y2, dflux[zbin+1];
    double *tc1, *gc1, *wc1, *tau1, *w1, *g1;
    double temp[zbin+1], *temp1;
    double surfaceref;
    double solarfraction;
    solarfraction = FADV/cos(THETAREF); //ms2021: global average, I guess
//    jTvar = nrl+1;
//    if (TIME_STEPPING) jTvar = 0;
    double lam1[NLAMBDA];
    for (l=0; l<NLAMBDA; l++) lam1[l] = wavelength[l]*1.0E-9; /* convert to m */
    double xlll;
//=========================================================
//=== START Two-Stream solver =============================
//
        /* Delta Adjustment */
    //ms22: this follows Joseph76, eqs. 5b,13&14
        if (DELADJUST == 1) {
            for (j=1; j <= zbin; j++) {
                tau[j] *= (1-w[j]*g[j]*g[j]);
                w[j]   =  (1-g[j]*g[j])*w[j]/(1-w[j]*g[j]*g[j]);
                g[j]   =  g[j]/(1+g[j]);
            }
        }
        
		/* Cumulative Optical Depth */
		for (j=0; j <= zbin; j++) {
			tc[j] = 0.0;
		}
		for (j=1; j <= zbin; j++) {
			for (jj=1; jj <= j; jj++) {
				tc[j] += tau[jj];
			}
		}
		
		/* Determine if multi-layer computation is needed */
		/* if (tc[zbin] < TAUTHRESHOLD) {
			
			B0      = 2*HPLANCK*pow(CLIGHT,2)/pow(lam1[i],5)/(exp(HPLANCK*CLIGHT/lam1[i]/KBOLTZMANN/temp[zbin])-1.0)*1.0E-9;
			flux[i] = ((1-PSURFAB)*PI*B0 + PSURFAB*cos(THETAREF)*exp(-tc[zbin]/cos(THETAREF))*solar[i]*exp(-tc[zbin]/miu1);
			ratio[i] = flux[i]/solar[i];
			
		} else { */
			
			/* Determine the optical depth scale for computation */
		
/* +++ OJO +++++++++++++++++++++++++++++++++++++++++++++++
 *the following should depend on stellar spectrum? 
 *and assumed surface albedos f[nu]?
 *+++ To be addressed ++++++++++++++++++++++++++++++++++++*/
		if (wavelength[lbin]>3000.0) {
			surfaceref = 1.0-PSURFEM;
		} else {
			surfaceref = PSURFAB;
		}
//----------------------------------------------------                
//markus: here starts the 2-stream with tri-diagonal solver:
//----------------------------------------------------                
        ntau=zbin;
        nvector = 2*ntau;

        /*if (lbin==4600) {
         printf("%f %s %d %f %f\n", wavelength[lbin], "ntau", ntau, tc[zbin], tc[0]);
         }*/
        
        tc1  = dvector(0,ntau);
        temp1= dvector(0,ntau);
        tau1 = dvector(1,ntau);
        w1   = dvector(1,ntau);
        g1   = dvector(1,ntau);
        
        for (j=1; j<=ntau; j++) {
            tau1[j]=tau[j];
            /*markus2021:testing for small tau*/
            //if(tau1[j]<1.e-6) tau1[j]=1.e-6;
            w1[j]=w[j];
            g1[j]=g[j];
        }
        for (j=0; j<=ntau; j++) {
            tc1[j]=tc[j];
        }
						
			gamma1 = dvector(1,ntau); //ms: closures
      			gamma2 = dvector(1,ntau); //ms: closures
			gamma3 = dvector(1,ntau); //ms: closures
			gamma4 = dvector(1,ntau); //ms: closures
			gamma  = dvector(1,ntau); //ms: for solving Toon89 eqs. 19&20 / 31&32
			lambda = dvector(1,ntau); //ms: for solving Toon89 eqs. 19&20 / 31&32
			e1     = dvector(1,ntau); //ms: to populate A,B,D,E
			e2     = dvector(1,ntau); //ms: to populate A,B,D,E
			e3     = dvector(1,ntau); //ms: to populate A,B,D,E
			e4     = dvector(1,ntau); //ms: to populate A,B,D,E
			CP0    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
			CP1    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
			CM0    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
			CM1    = dvector(1,ntau); //ms: Sources for use in eqs. 41&42
                        CPh    = dvector(1,ntau); //ms: Sources in eqs. 23&27
                        CMh    = dvector(1,ntau); //ms: Sources in eqs. 24&27
//ms: For 2N eqs. in Tridiagonal matrix: A_l Y_(l-1) + B_l Y_l +D_l Y_(l+1) = E_l
			A      = dvector(1,nvector);
			B      = dvector(1,nvector);
			D      = dvector(1,nvector);
			E      = dvector(1,nvector);
			AS     = dvector(1,nvector);
			DS     = dvector(1,nvector);
			ES     = dvector(1,nvector);
			Y      = dvector(1,nvector);
			Y1     = dvector(1,ntau);
			Y2     = dvector(1,ntau);
			
		/* Two stream parameters */
		if (wavelength[lbin] > 3000.0) {
                    //Hemispheric Mean:
			for (j=1; j <= ntau; j++) {
				gamma1[j] = 2.0-w1[j]*(1+g1[j]);
				gamma2[j] = w1[j]*(1.0-g1[j]);
				gamma3[j] = 0.5-g1[j]*cos(THETAREF); //ms: haven't seen this one before.
				gamma4[j] = 1.0-gamma3[j];
			}
		} else {
			/* Eddington */
			for (j=1; j <= ntau; j++) {
				gamma1[j] = (7.0-(4.0+3.0*g1[j])*w1[j])/4.0;
				gamma2[j] = -(1.0-(4.0-3.0*g1[j])*w1[j])/4.0;
				gamma3[j] = (2.0-3.0*g1[j]*cos(THETAREF))/4.0;
				gamma4[j] = 1.0-gamma3[j];
			}
		}

			/* Two-stream computation scheme */
			/* some expressions */
//----------------------------------------------------                
//ms: All after Toon+ 1989: 
//----------------------------------------------------                
	for (j=1; j <= ntau; j++) {
		lambda[j] = sqrt(gamma1[j]*gamma1[j]-gamma2[j]*gamma2[j]); //ms: eq.21
		gamma[j]  = gamma2[j]/(gamma1[j]+lambda[j]); //ms: eq. 22
		e1[j] = 1.0+gamma[j]*exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		e2[j] = 1.0-gamma[j]*exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		e3[j] = gamma[j]+exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		e4[j] = gamma[j]-exp(-lambda[j]*tau1[j]); //ms: eqs. 44
		/* printf("%d %f %f %f %f\n", j, e1[j], e2[j], e3[j], e4[j]); */
	}
        
        /*fprintf(fp, "%e\t%e\t", wavelength[lbin], solar[lbin]); */
        
        /* Loop the temperature variation */
        for (kk=0; kk<=jTvar; kk++) {
            
            /* Get temperature */
            for (j=0; j<=ntau; j++) {
                temp1[j] = Tvar[j][kk];
            }

            /* source functions */
	    for (j=1; j <= ntau; j++) {
		CP0[j] = 0.0;
		CP1[j] = 0.0;
		CM0[j] = 0.0;
		CM1[j] = 0.0;
                CPh[j] = 0.0;
		CMh[j] = 0.0;
                //printf("%s%d%s %e\n", "CP0[",j,"]A",CP0[j]);
                
                if (wavelength[lbin] > 3000.0) {
		    B0 = 2*HPLANCK*pow(CLIGHT,2)/pow(lam1[lbin],5)/(exp(HPLANCK*CLIGHT/lam1[lbin]/KBOLTZMANN/temp1[j-1])-1.0)*1.0E-9;
                    B1 = 0.0;
		    //if (tau1[j] > 1e-4){ //ms: isothermal if tau <<
                    B1 = ( 2*HPLANCK*pow(CLIGHT,2)/pow(lam1[lbin],5)/(exp(HPLANCK*CLIGHT/lam1[lbin]/KBOLTZMANN/temp1[j])-1.0)*1.0E-9 - B0 )/tau1[j];
                    //}

                    //ms: why += instead of =?
		    CP0[j] += 2*PI*miu1*(B0+B1/(gamma1[j]+gamma2[j]));
		    CP1[j] += 2*PI*miu1*(B0+B1*(tau1[j]+1.0/(gamma1[j]+gamma2[j])));
                
                    CPh[j] += 2*PI*miu1*(B0+B1*(tau1[j]*0.5+1.0/(gamma1[j]+gamma2[j]))); //ms: eq. 27

		    CM0[j] += 2*PI*miu1*(B0-B1/(gamma1[j]+gamma2[j]));
		    CM1[j] += 2*PI*miu1*(B0+B1*(tau1[j]-1.0/(gamma1[j]+gamma2[j])));
                
                    CMh[j] += 2*PI*miu1*(B0+B1*(tau1[j]*0.5-1.0/(gamma1[j]+gamma2[j]))); //ms: eq. 27
                } else {
               
                //printf("%s%d%s %e\n", "CP0[",j,"]B",CP0[j]);
		    CP0[j] += solarfraction*solar[lbin]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                //printf("%s%d%s %e\n", "CP0[",j,"]C",CP0[j]);
		    CP1[j] += solarfraction*solar[lbin]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                        
                    CPh[j] += solarfraction*solar[lbin]*w1[j]*exp(-(tc1[j-1]+0.5*tau1[j])/cos(THETAREF))*((gamma1[j]-1.0/cos(THETAREF))*gamma3[j]+gamma2[j]*gamma4[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF)); //ms:eq. 23
				
                    CM0[j] += solarfraction*solar[lbin]*w1[j]*exp(-tc1[j-1]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
		    CM1[j] += solarfraction*solar[lbin]*w1[j]*exp(-tc1[j]/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF));
                                
                    CMh[j] += solarfraction*solar[lbin]*w1[j]*exp(-(tc1[j-1]+0.5*tau1[j])/cos(THETAREF))*((gamma1[j]+1.0/cos(THETAREF))*gamma4[j]+gamma2[j]*gamma3[j])/(lambda[j]*lambda[j]-1.0/cos(THETAREF)/cos(THETAREF)); //ms:eq.24
                    
                }
                /*printf("%s%d%s%d%s %e\n", "solarfraction*solar[",i,"]*w1[",j,"]",solarfraction*solar[i]*w1[j]);
                printf("%s%d%s %e\n", "tc1[",j-1,"]",tc1[j-1]);
                printf("%s%d%s %e\n", "tc1[",j,"]",tc1[j]);
                printf("%s %e\n", "cos(thetaref)",cos(THETAREF));
                printf("%s%d%s %e\n", "gamma1[",j,"]",gamma1[j]);
                printf("%s%d%s %e\n", "gamma2[",j,"]",gamma2[j]);
                printf("%s%d%s %e\n", "gamma3[",j,"]",gamma3[j]);
                printf("%s%d%s %e\n", "gamma4[",j,"]",gamma4[j]);
                printf("%s%d%s %e\n", "lambda[",j,"]",lambda[j]);*/
                 /*printf("%s%d%s %e\n", "CP0[",j,"]",CP0[j]);
                 printf("%s%d%s %e\n", "CP1[",j,"]",CP1[j]);
                 printf("%s%d%s %e\n", "CM0[",j,"]",CM0[j]);
                 printf("%s%d%s %e\n", "CM1[",j,"]",CM1[j]);*/
                
				/* printf("%s %f %d %e\n", "CP0", wavelength[lbin], j, CP0[j]);
                 printf("%s %f %d %e\n", "CP1", wavelength[lbin], j, CP1[j]);
                 printf("%s %f %d %e\n", "CM0", wavelength[lbin], j, CM0[j]);
                 printf("%s %f %d %e\n", "CM1", wavelength[lbin], j, CM1[j]); */
            }
            
            /* equations */
			for (l=1; l<=nvector; l++) {
				A[l] = 0.0;
				B[l] = 0.0;
				D[l] = 0.0;
				E[l] = 0.0;
				Y[l] = 0.0;
				AS[l] = 0.0;
				DS[l] = 0.0;
				ES[l] = 0.0;
			}
			A[1] = 0;
			B[1] = e1[1];
			D[1] = -e2[1];
			E[1] = -CM0[1];
			for (j=1; j < ntau; j++) {
				l = 2*j+1;
				A[l] = e2[j]*e3[j] - e4[j]*e1[j];
				B[l] = e1[j]*e1[j+1] - e3[j]*e3[j+1];
				D[l] = e3[j]*e4[j+1] - e1[j]*e2[j+1];
				E[l] = e3[j]*(CP0[j+1]-CP1[j]) - e1[j]*(CM0[j+1]-CM1[j]);
			}
			for (j=1; j < ntau; j++) {
				l = 2*j;
				A[l] = e2[j+1]*e1[j] - e3[j]*e4[j+1];
				B[l] = e2[j]*e2[j+1] - e4[j]*e4[j+1];
				D[l] = e1[j+1]*e4[j+1] - e2[j+1]*e3[j+1];
                                //markus2021: possible error in Toon89: -/+ sign changed:
				E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) - e4[j+1]*(CM0[j+1]-CM1[j]);
				//E[l] = e2[j+1]*(CP0[j+1]-CP1[j]) + e4[j+1]*(CM0[j+1]-CM1[j]);
			}
			l = nvector;
            surfaceref = 1.0;
			A[l] = e1[ntau] - surfaceref*e3[ntau];
			B[l] = e2[ntau] - surfaceref*e4[ntau];
			D[l] = 0;
			E[l] = -CP1[ntau] + surfaceref*CM1[ntau];
                        //ms: not sure if we can add Planck (T_Surf) and Direct Stellar beam together in these equations.
			// E[l] += (1-surfaceref)*PI*2*HPLANCK*pow(CLIGHT,2)/pow(lam1[lbin],5)/(exp(HPLANCK*CLIGHT/lam1[lbin]/KBOLTZMANN/temp1[ntau])-1.0)*1.0E-9;
            E[l] += PI*2*HPLANCK*pow(CLIGHT,2)/pow(lam1[lbin],5)/(exp(HPLANCK*CLIGHT/lam1[lbin]/KBOLTZMANN/Tint)-1.0)*1.0E-9;
			/* printf("%s %f %e\n","surf temp",temp1[ntau-1],tc1[ntau]); */
			E[l] += surfaceref*cos(THETAREF)*solar[lbin]*solarfraction*exp(-tc1[ntau]/cos(THETAREF));
			/* for (l=1; l<=nvector; l++) {
             printf("%s %f %d %f %f %f %f\n", "ABDE", wavelength[lbin], l, A[l], B[l], D[l], E[l]);
             } */
            
            /* solve equations */
//ms: solving the tridiagonal matrix
			l=nvector;
			AS[l] = A[l]/B[l];
			DS[l] = E[l]/B[l];
			for (l=nvector-1; l>0; l--) {
				xlll = 1.0/(B[l] - D[l]*AS[l+1]);
				AS[l] = A[l]*xlll;
				DS[l] = (E[l]-D[l]*DS[l+1])*xlll;
			}
			Y[1] = DS[1];
			for (l=2; l <= nvector; l++) {
				Y[l] = DS[l] - AS[l]*Y[l-1];
			}
            
            for (j=1; j <= ntau; j++) {
                Y1[j] = Y[2*j-1];
                Y2[j] = Y[2*j];
            }
            
//ms: What's wrong with aflux below + direct flux, as eqs. 48-51 in Toon+1989?
//
            /* obtain net flux at zero and midpoint */
            /*for (j=1; j<=ntau; j++) {
                aflux[j] = Y1[j]*(e1[j]-e3[j])+Y2[j]*(e2[j]-e4[j])+CP1[j]-CM1[j];
            }*/
            dflux[0] = Y1[1]*e3[1] - Y2[1]*e4[1] + CP0[1]; //ms:should be ok, but let's write it out:
            //dflux[0] = Y1[1]*(e1[1]-e3[1])+Y2[1]*(e2[1]-e4[1])+CP0[1]-CM0[1];
            for (j=1; j<=ntau; j++) {
                dflux[j] =  2.0*Y2[j]*(1.0-gamma[j])*exp(-lambda[j]*tau1[j]*0.5) + CPh[j] - CMh[j]; //ms:should be ok, but let's write it out:
                //dflux[j] = Y1[j]*(e1[j]-e3[j])+Y2[j]*(e2[j]-e4[j])+CP1[j]-CM1[j];
            }
            
            /* fprintf(fp, "%e\t", dflux[0]);
            
            if (kk==nrl+1) {
                fprintf(fp, "\n");
            } */
            
            /* add direct sunlight */
            dflux[0] -= cos(THETAREF)*solar[lbin]*solarfraction;
            for (j=1; j<=ntau; j++) {
                dflux[j] -= cos(THETAREF)*solar[lbin]*solarfraction*exp(-(tc1[j-1]+tau1[j]*0.5)/cos(THETAREF));
            }

//----------------------------------------------------                
//markus: I guess here ends the 2-stream flux calculation
//----------------------------------------------------                

            
            NetFlux[0][kk] += dflux[0]*(wavelength[lbin+1]-wavelength[lbin]); //ms: could be within loop...
            radcount = 0;
            for (k=1; k<=zbin; k++) {
                if (isconv[zbin+1-k] == 0) {
                    radcount += 1;
                    NetFlux[radcount][kk] += dflux[k]*(wavelength[lbin+1]-wavelength[lbin]);
                }
            }
            
/*            if (kk==0) {
                radiationI1 += cos(THETAREF)*solar[lbin]*solarfraction*exp(-tc[zbin]/cos(THETAREF))*(wavelength[lbin+1]-wavelength[lbin]); // BOA incoming flux //
                radiationI0 += cos(THETAREF)*solar[lbin]*solarfraction*(wavelength[lbin+1]-wavelength[lbin]); // TOA incoming flux //
//moved downward                radiationO += dflux[0]*(wavelength[lbin+1]-wavelength[lbin]); // TOA net outgoing flux //
                //ms: actually sth like a version of outgoing, see l:494
            }*/
//if(lbin%100==0 && kk==0){            
/*if(lbin>=6000 && lbin<8000 && kk==0){            
printf("%s\n","\n==== CHECKING 2-STREAM ===="); //template
printf("%s\t%d\n","wavelength bin i = ",lbin); //template
printf("%s\n","\n==== (i) X Y Z S Lam Tau ===="); //template
for (i=1;i<=nvector;i++)
{
    printf("%d\t%.3e\t%.3e\t%.3e\t%.3e\t%.50e\t%.50e\n",i,A[i],B[i],D[i],E[i],lam1[lbin],tau1[i]);
}*/
/*printf("%s\n","\n==== Up/Dn FLUXES ===="); //template
for (i=0;i<=ntau;i++)
{
    printf("%s %d\t%.3e\n","Fup - Fdn ",i,dflux[i]);
}
printf("%s\n","\n==== NET FLUXES ===="); //template
for (i=1;i<=nrl;i++)
{
    printf("%s %d\t%.3e\n","NetFlux ",i,NetFlux[i][kk]);
}
}*/
        } // END Tvar loop
        

			
            /*markus2021 */
            //if (i==4530 || i==4531) {
            /*printf("%s%i\t","i= ",i);                        
            printf("%s\t", "tau1:");
            printf("%f\n",tau1[ntau]);
        if (i==15998) {
            printf("%s%i\n","i= ",i);
            printf("%s\t", "tau1:");

            //printf("%f\n",tau1[ntau]);
            for (j=1; j<=ntau; j++) {
                printf("%f\t", tau1[j]);
            }
            printf("\n"); 

            printf("%s%f\n", "B1: ",B1);
            
            printf("%s\t", "E:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", E[l]);
            }
            printf("\n"); 

            printf("%s\t", "AS:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", AS[l]);
            }
            printf("\n"); 

            printf("%s\t", "DS:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", DS[l]);
            }
            printf("\n"); 

            printf("%s\t", "Y:");
	    for (l=1; l<=nvector; l++) {
                printf("%f\t", Y[l]);
            }
            printf("\n"); 

            printf("%s\t", "dflux:");
            //printf("%f\n", dflux[ntau]);
            for (j=0; j<=ntau; j++) {
                printf("%f\t", dflux[j]);
            }
            printf("\n"); 

        }*/
            

			free_dvector(tc1,0,ntau);
			free_dvector(temp1,0,ntau);    
			free_dvector(w1,1,ntau);
			free_dvector(g1,1,ntau);
			free_dvector(tau1,1,ntau);
			free_dvector(gamma1,1,ntau);
			free_dvector(gamma2,1,ntau);
			free_dvector(gamma3,1,ntau);
			free_dvector(gamma4,1,ntau);
			free_dvector(gamma,1,ntau);
			free_dvector(lambda,1,ntau);
			free_dvector(e1,1,ntau);
			free_dvector(e2,1,ntau);
			free_dvector(e3,1,ntau);
			free_dvector(e4,1,ntau);
			free_dvector(CP0,1,ntau);
			free_dvector(CP1,1,ntau);
			free_dvector(CM0,1,ntau);
			free_dvector(CM1,1,ntau);
                        free_dvector(CPh,1,ntau);
                        free_dvector(CMh,1,ntau);
			free_dvector(Y1,1,ntau);
			free_dvector(Y2,1,ntau);
			free_dvector(Y,1,nvector);
			free_dvector(A,1,nvector);
			free_dvector(B,1,nvector);
			free_dvector(D,1,nvector);
			free_dvector(E,1,nvector);
			free_dvector(AS,1,nvector);
			free_dvector(DS,1,nvector);
			free_dvector(ES,1,nvector);
        



//atexit(pexit);exit(0); //ms debugging mode

//clock_gettime(CLOCK_REALTIME, &end);
//t2_passed = ((double)end.tv_sec*1e9 + end.tv_nsec) - ((double)start.tv_sec*1e9 + start.tv_nsec);
//printf("\n%s %f\n","ms_2stream took [ms]:", t2_passed/1e6); 
//atexit(pexit);exit(0); //ms debugging mode
//=== HERE ends the temperature variation loop ============
//printf("%s\n", "END NEW 2 STREAM"); //template
//=========================================================

/*printf("%s\n","\n==== DEBUGGING ===="); //template
printf("%s %d\n","\n==== Up/Dn FLUXES",lbin," ===="); //template
for (i=0;i<zbin;i++)
{
    printf("%s %d\t%.3e\t%.3e\n","Fup,Fdn ",i,Fmat[2*i],Fmat[2*i+1]);
}
//    printf("%s %d\t%.3e\t%.3e\n","Fup,Fdn ",zbin-1,Fmat[2*(zbin-1)],Fmat[2*(zbin-1)+1]);
//    printf("%s %d\t%.3e\t%.3e\n","Fup,Fdn ",zbin,Fmat[2*(zbin)],Fmat[2*(zbin)+1]);
printf("%s\n","==== END FLUXES ===="); //template
*/

//printf("%s\t%e\n","Mu = cos(THETAREF) = ",mu[0]); //check
}
// END: void ms_two_str_solver()

//atexit(pexit);exit(0); //ms debugging mode
//printf("%s\n", "HERE"); //template
