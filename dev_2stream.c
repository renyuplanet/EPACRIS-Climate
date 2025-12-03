/*******************************************************************
 * File: ms_2stream.c
 * Function: Solve for equilibrium 2-stream fluxes from Heng+2018, Eqs. 25&26 (doi:10.3847/1538-4365/aad199)
 * Original Author: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Version: v01, 2022: Original upload
 *******************************************************************
 * (All (Eq.xx) refer to Heng+2018 unless stated otherwise)
 * Eq.25 & 26 represent Fup & Fdn for a pair of atmospheric layers i and i-1(above).
 * (Eq.25) Fup(i) = 1/Xeq { -Yeq Fup(i+1) - Zeq Fdn(i) + PIt [ Xeq B+(i) + Yeq B+(i+1) + Zeq B-(i) ]
 *                  + Yeq C-TsTc + Tc [ zeta+ ( C+zeta- - C-zeta+ ) + zeta-T^2 ( C-zeta- - C+zeta+  )  ] }
 * (Eq.26) Fdn(i) = 1/Xeq  { -Yeq Fdn(i-1) - Zeq Fup(i) + PIt [ Xeq B-(i) + Yeq B-(i-1) + Zeq B+(i) ]
 *                  + Yeq C+Tc + TsTc [ zeta+ ( C-zeta- - C+zeta+ ) + zeta-T^2 ( C+zeta- - C-zeta+  )  ] }
*/                  

#include <math.h>
#include "constant.h"

//=========================================================
//=== Function identifiers=================================
void ms_two_str_solver(int lbin, double w[], double g[], double tau[], int nrl, int isconv[], double **Tvar, int jTvar, double P[], double **NetFlux, double Fup[], double Fdn[], double Tint, double sumBint);

//=========================================================
//=== Functions ===========================================
void ms_two_str_solver(int lbin, double w[], double g[], double tau[], int nrl, int isconv[], double **Tvar, int jTvar, double P[], double **NetFlux, double Fup[], double Fdn[], double Tint, double sumBint)
{
/* This function solves for Fnet using inversion of a tri-diagonal matrix. It includes a loop (tv) for temperature variation to populate a jacobian matrix later, to solve for equilibrium.
 *Equations setup for tri-diagonal matrix (see Toon+89): X(i) F(i-1) + Y(i) F(i) + Z(i) F(i+1) = S(i)
 *
 * => (Eq.25) 
 * Xeq Fup(i) + Zeq Fdn(i) + Yeq Fup(i+1) = 
 * PIt [ Xeq B+(i) + Yeq B+(i+1) + Zeq B-(i) ] + Yeq C-TsTc + Tc [ zeta+ ( C+zeta- - C-zeta+ ) + zeta-T^2 ( C-zeta- - C+zeta+  )  ]
 * 
 * & (Eq.26)
 * Yeq Fdn(i-1) + Zeq Fup(i) + Xeq Fdn(i) =
 * PIt [ Xeq B-(i) + Yeq B-(i-1) + Zeq B+(i) ] + Yeq C+Tc + TsTc [ zeta+ ( C-zeta- - C+zeta+ ) + zeta-T^2 ( C+zeta- - C-zeta+  )  ]
*/ 
//=========================================================
//printf("%s\n", "START NEW 2 STREAM"); //template
//=== DECLARE variables ===================================

    double time,t_passed,t2_passed;
    struct timespec start,end;
//clock_gettime(CLOCK_REALTIME, &start);

    int i, j, ii, tv, rl;   //...loops: usually arrays are filled from 1 .. N (not from 0)
    double E[zbin+1];
    double mu[zbin+1],muD[zbin+1],ffrac[zbin+1];
    double zetap[zbin+1],zetam[zbin+1];
    double eps2[zbin+1],gammaA[zbin+1],gammaS[zbin+1],gammaD[zbin+1];
    double PIt[zbin+1];
    double CD[zbin+1],CP[zbin+1],CM[zbin+1];
    double dtau[zbin+1],Trans[zbin+1],TransD[zbin+1],TransC[zbin+1],tauC[zbin+1];
    double B[zbin+1],dB[zbin+1],BPup[zbin+1],BMup[zbin+1],BPlo[zbin+1],BMlo[zbin+1], Bint; //helper terms
    double Xeq[zbin+1],Yeq[zbin+1],Zeq[zbin+1]; //X,Y,Z in Eq25&26 description above
    double X[2*zbin+2],Y[2*zbin+2],Z[2*zbin+2],S[2*zbin+2]; //A,B,D,E in Toon+89
    double Xdash[2*zbin+2],Sdash[2*zbin+2],Zdash[2*zbin+2];   //AS and DS in Toon+89
    double Fmat[2*zbin+2]; //matrix solution vector Y in Toon+89
    
//- FLux Matrix ------------------------------------------
    double **fcoef, *Scoef;
    fcoef = dmatrix(1,2*zbin+2,1,2*zbin+2); //Btridiag
    Scoef = dvector(1,2*zbin+2); //Btridiag&penta
    double *AA, *BB, *CC, *DD, *EE;
    AA = dvector(1,2*zbin+2); //penta
    BB = dvector(1,2*zbin+2); //penta
    CC = dvector(1,2*zbin+2); //penta
    DD = dvector(1,2*zbin+2); //penta
    EE = dvector(1,2*zbin+2); //penta

    double Emat[2*zbin+2], Cmat[2*zbin+2], Dmat[2*zbin+2], Amat[2*zbin+2], Bmat[2*zbin+2], Smat[2*zbin+2]; //for pentadiagonal solvers
    double alphaPen[2*zbin+2], betaPen[2*zbin+2], gammaPen[2*zbin+2], muPen[2*zbin+2], zPen[2*zbin+2]; //for PTRANS-I
    double fSog[2*zbin+2], gSog[2*zbin+2], cSog[2*zbin+2], eSog[2*zbin+2], ytilSog[2*zbin+2]; //for Sogabe solver 
    double cKPen[2*zbin+2], fKPen[2*zbin+2], eKPen[2*zbin+2], vKPen[2*zbin+2], hKPen[2*zbin+2], zKPen[2*zbin+2]; //for KPENTA solver

//=========================================================
//=== INITIALIZE variables ================================

/*    for (i=0; i<=zbin; i++)
    {
        E[i] = 0.0;
    }
*/
//=========================================================
//printf("%s\n", "POPULATE helper terms"); //template
//=== POPULATE helper terms (mostly Heng+2018) ============
    tauC[0] = 0.0;  //initializing 'TOA'
    for (i=1;i<=zbin;i++)   //Layloop (through layers)
    {
//printf("%s %d\n", "loop i = ", i); //template
        if(w[i] == 1.0) printf("%s %d\n", "WARNING: w[",i,"] = 1.0 at lambda: ", wavelength[lbin]); //template
    //E-factor (Eq.31)
//        printf("%s %f\n", "g[i] = ",g[i]);
//        printf("%s %f\n", "w[i] = ",w[i]);
        E[i] = 1.225 - 0.1582*g[i] - 0.1777*w[i] - 0.07465*g[i]*g[i] + 0.2351*w[i]*g[i] - 0.05582*w[i]*w[i];
        //E[i] = 1.0; //testing
//        printf("%s %f\n", "E[i] = ",E[i]);
        if(E[i]<1.0)
        {
            //if(g[i]>0.0) printf("%s %f %s\n","WARNING: E-factor",E[i],"<1.0, adjusted!");
            //printf("%s %f %s\n","WARNING: E-factor",E[i],"<1.0, adjusted!");
            E[i] = 1.0;
        }
//printf("%s\n", "E done"); //template
    
    //angular parameters
        mu[i] = cos(THETAREF);  //...for now constant
        muD[i] = cos(THETAREF); //...for now constant and the same
        //ffrac[i] = FADV/muD[i]; //...depends on the definition of FADV, since FADV*muD = 1/4 for global energy balance
        ffrac[i] = FADV; //ms05/22 testing
//printf("%s\t%e\n","Mu = cos(THETAREF) = ",mu[1]); //check
//atexit(pexit);exit(0); //ms debugging mode
        //for numerical stability
        if(muD[i] == 0.5*sqrt(E[i]*(E[i]-w[i])*(1-w[i]*g[i])) ) 
        {    
            muD[i] *= 0.99;
            printf("%s\n","WARNING: muD adjustment actually happened!");
        }
//printf("%s\n", "mu done"); //template

    //Coupling Coefficients (Eq.24)
        zetap[i] = 0.5*(1+sqrt( (E[i]-w[i]) / (E[i]*(1-w[i]*g[i])) ));
        zetam[i] = 0.5*(1-sqrt( (E[i]-w[i]) / (E[i]*(1-w[i]*g[i])) ));
//printf("%s\n", "zeta done"); //template

    //2-stream closures
        //eps2[i] = 0.5*(1-sqrt(3)*g[i]*w[i]);    //...2nd Eddington coefficient (Toon+89, Tbl.1 quadrature closure)
        //eps2[i] = 1.0/sqrt(3);    //... quadrature closure (Eq.21)
        eps2[i] = 2.0/3.0;    //... Eddington closure (Eq.21)
        gammaA[i] = 2*E[i] - w[i]*(1+E[i]*g[i]); //...absorbtion (Eq.20a)
        gammaS[i] = w[i]*(1-E[i]*g[i]);          //...scattering (Eq.20b)
        gammaD[i] = w[i]*g[i]*muD[i]*ffrac[i]*solar[lbin] / eps2[i]; //...direct beam (Eq.18)
//        printf("%s %f\n", "gammaA[i] = ",gammaA[i]);
//        printf("%s %f\n", "gammaS[i] = ",gammaS[i]);
//printf("%s\n", "gamma done"); //template

    //PI term (Eq.27a)
        PIt[i] = PI*( (1-w[i]) / (E[i]-w[i]) );  
//printf("%s\n", "PIt done"); //template

    //Direct beam coefficients
        CD[i] = w[i]*ffrac[i]*solar[lbin]*( 2*E[i]*(1-w[i]*g[i]) + g[i]/eps2[i] );  //... see C* (Eq.23)
        CD[i] /= 4*E[i]*( E[i] - w[i]  )*( 1 - w[i]*g[i] ) - 1/( muD[i]*muD[i] );
        CP[i] = 0.5* ( CD[i] + ( CD[i]/muD[i] + gammaD[i] ) / ( 2*E[i]*(1 - w[i]*g[i]) ) ); //... (Eq.27b)
        CM[i] = 0.5* ( CD[i] - ( CD[i]/muD[i] + gammaD[i] ) / ( 2*E[i]*(1 - w[i]*g[i]) ) ); //... (Eq.27b)
//printf("%s\n", "CM/CP done"); //template

    //Transmission functions (Eq.27c-f)
//        printf("%s %e\n", "tau[i] = ",tau[i]); //OJO: tau defined per layer! not cumulative! i.e.: tau[i] = dtau[i]
//        dtau[i] = tau[i] - tau[i-1]; //... TOA: tau[0] should be 0, this puts fluxes at cell boundaries!
//        printf("%s %f\n", "dtau[i] = ",dtau[i]);
        tauC[i] = tauC[i-1] + tau[i];   //... tau above (cumulative)
//        printf("%s %e\n", "tauC[i] = ",tauC[i]);
        Trans[i] = exp( -tau[i]*sqrt( (gammaA[i]+gammaS[i]) * (gammaA[i]-gammaS[i]) ) );     //... T
//        printf("%s %e\n", "trans[i] = ",Trans[i]);
        TransD[i] = exp( -tau[i] / muD[i]  );    //.. .T*
//        printf("%s %e\n", "transD[i] = ",TransD[i]);
        TransC[i] = exp( -tauC[i-1] / muD[i]  );    //... Tabove
//        printf("%s %e\n", "transC[i] = ",TransC[i]);
//printf("%s\n", "Trans done"); //template

    //Abbreviations in Eq25&26
        Xeq[i] = ( zetam[i]*zetam[i]*Trans[i]*Trans[i] - zetap[i]*zetap[i] );
        Yeq[i] = Trans[i] * ( zetap[i]*zetap[i] - zetam[i]*zetam[i] );
        Zeq[i] = zetam[i]*zetap[i] * ( 1 - Trans[i]*Trans[i] );
//printf("%s\n", "XYZ done"); //template
        if(isnan(Xeq[i])) {printf("%s %d","isNaN Xeq ",i);atexit(pexit);exit(0);} 
        if(isnan(Yeq[i])) {printf("%s %d","isNaN Yeq ",i);atexit(pexit);exit(0);} 
        if(isnan(Zeq[i])) {printf("%s %d","isNaN Zeq ",i);atexit(pexit);exit(0);} 
    }//END: Layloop

//=========================================================
//printf("%s\n", "TVar loop"); //template
//=== HERE starts the temperature variation loop ==========

//    jTvar = nrl+1;
//    if (TIME_STEPPING) jTvar = 0; //Only fill Tvar[:][0] elements to avoid excessive runtime
    for(j=0; j<=jTvar; j++) //Tvar loop
    {  

//=========================================================
//printf("%s\n", "CALCULATE Planck"); //template
//=== Calculate Planck functions and derivatives ==========
        for(i=0; i<=zbin; i++) //Layloop 2; Thermal atmospheric radiation
        {
//          printf("%s %e\n", "wavelength[lbin] = ",wavelength[lbin]);
            B[i] = 2*HPLANCK*CLIGHT*CLIGHT / pow(wavelength[lbin]*1e-9,5) * 1e-9 ;
            B[i] /= exp(HPLANCK*CLIGHT / wavelength[lbin]/1e-9/KBOLTZMANN/Tvar[i][j]) - 1.0 ;
            if(i>0)
            {
                if(tau[i] > 1.0e-4) dB[i] = (B[i] - B[i-1]) / tau[i]; //assumes opacities const. over Tvar range... 
                else 
                {
                    dB[i] = 0.0; //isothermal approach where no absorption
                    //printf("%s %d %s %d %s %e\n", "Layer",i,"set to isothermal, tau[i] < 1e-4 in lbin",lbin,"at lambda ",wavelength[lbin]); //template
                }
            }
        }//END: Layloop 2
        //Internal Temperature
            Bint = 2*HPLANCK*CLIGHT*CLIGHT / pow(wavelength[lbin]*1e-9,5) * 1e-9 ;
            Bint /= exp(HPLANCK*CLIGHT / wavelength[lbin]/1e-9/KBOLTZMANN/Tint) - 1.0 ;
            
            if(lbin==0) sumBint = 0.0;
            if(j==0) sumBint += Bint;

        for(i=0; i<=zbin; i++) //Layloop 3; setting up planck related functions (Eq.27g)
        {
            if(i>0)
            {
                BPup[i] = B[i] + dB[i] / ( 2*E[i]*(1 - w[i]*g[i])  );
                BMup[i] = B[i] - dB[i] / ( 2*E[i]*(1 - w[i]*g[i])  );
            }
            if(i<zbin)
            {
                BPlo[i] = B[i] + dB[i+1] / ( 2*E[i+1]*(1 - w[i+1]*g[i+1])  );
                BMlo[i] = B[i] - dB[i+1] / ( 2*E[i+1]*(1 - w[i+1]*g[i+1])  );
            }
        }//END: Layloop 3

//printf("%s\n", "POPULATE matrix"); //template

//==================================================================================================
//=== POPULATE flux matrix coefficients ============================================================
//==================================================================================================
    if(RT_FLUX_SOLVER < 3)
    {
        //This correlates Xeq,Yeq,Zeq to the matrix coefficients for X(ii) Fmat(ii-1) + Y(ii) Fmat(ii) + Z(ii) Fmat(ii+1) = S(ii)
        // Matrix sequence: Eq.26 (1), Eq.25 (1), Eq.26 (2), Eq.25 (2), ... Eq.25(2N-1), Eq.26 (2N), Boundary Eq. (2N)

        //TOA:
            Fmat[1] = muD[1]*ffrac[1]*solar[lbin]; //incoming TOA flux
            //LHS: Flux terms   
            X[1] = 0.0;
            Y[1] = Zeq[1];
            Z[1] = Xeq[1];
            //RHS: Source terms
            S[1] = PIt[1] * ( BMup[1]*Xeq[1] + BMlo[0]*Yeq[1] + BPup[1]*Zeq[1] ) + CP[1]*TransC[1]*Yeq[1];
            S[1] += TransD[1]*TransC[1] * ( zetap[1]*( CM[1]*zetam[1]-CP[1]*zetap[1] ) + zetam[1]*Trans[1]*Trans[1]*( CP[1]*zetam[1]-CM[1]*zetap[1] ) );
            S[1] -= Yeq[1] * Fmat[1]; //Since we use Eq.26 here without Fdn0, we need to shift it to RHS as Source
            //For Fup[0] i.e. Fmat[0] using Eq.25 after matrix is solved:
            S[0] = PIt[1] * ( BPlo[0]*Xeq[1] + BPup[1]*Yeq[1] + BMlo[0]*Zeq[1] ) + CM[1]*TransD[1]*TransC[1]*Yeq[1]; 
            S[0] += TransC[1] * ( zetap[1]*( CP[1]*zetam[1]-CM[1]*zetap[1] ) + zetam[1]*Trans[1]*Trans[1]*( CM[1]*zetam[1]-CP[1]*zetap[1] ) );

        //BOA:
            //LHS: Flux terms
            X[2*zbin] = 1.0;
            //Y[2*zbin] = -PSURFAB;
            Y[2*zbin] = -1.0; //testing
            Z[2*zbin] = 0.0;
            //RHS: Source terms
            //S[2*zbin] = PSURFAB * ( muD[zbin]*ffrac[zbin]*solar[lbin]*TransC[zbin] ) + PSURFEM * PIt[zbin] * B[zbin]; // Reflected Fdn [zbin] + emmission of surface (less the scattered part, hence PIt[]) 
            //S[2*zbin] = PSURFAB * ( muD[zbin]*ffrac[zbin]*solar[lbin]*TransC[zbin] ) + PI*Bint; //testing
            //S[2*zbin] =  PI* (Bint + B[zbin]); //testing
            S[2*zbin] =  PI* Bint; //testing
            //S[2*zbin] =  0.0; //testing

        //remainig matrix:
        for(ii=2; ii<=zbin; ii++) //POPUloop
        {
        //LHS: Flux terms
            //Eq.25
            X[2*ii-2] = Xeq[ii];
            Y[2*ii-2] = Zeq[ii];
            Z[2*ii-2] = Yeq[ii];
            //Eq.26
            X[2*ii-1] = Yeq[ii];
            Y[2*ii-1] = Zeq[ii];
            Z[2*ii-1] = Xeq[ii];
        //RHS: Source terms
            //Eq.25
            S[2*ii-2] = PIt[ii] * ( BPlo[ii-1]*Xeq[ii] + BPup[ii]*Yeq[ii] + BMlo[ii-1]*Zeq[ii] ) + CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            S[2*ii-2] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
            //Eq.26
            S[2*ii-1] = PIt[ii] * ( BMup[ii]*Xeq[ii] + BMlo[ii-1]*Yeq[ii] + BPup[ii]*Zeq[ii] ) + CP[ii]*TransC[ii]*Yeq[ii];
            S[2*ii-1] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );

        if(isnan(X[2*ii-2])) {printf("%s %d","isNaN X ",2*ii-2);atexit(pexit);exit(0);} 
        if(isnan(X[2*ii-1])) {printf("%s %d","isNaN X ",2*ii-1);atexit(pexit);exit(0);} 
        if(isnan(Y[2*ii-2])) {printf("%s %d","isNaN Y ",2*ii-2);atexit(pexit);exit(0);} 
        if(isnan(Y[2*ii-1])) {printf("%s %d","isNaN Y ",2*ii-1);atexit(pexit);exit(0);} 
        if(isnan(Z[2*ii-2])) {printf("%s %d","isNaN Z ",2*ii-2);atexit(pexit);exit(0);} 
        if(isnan(Z[2*ii-1])) {printf("%s %d","isNaN Z ",2*ii-1);atexit(pexit);exit(0);} 
        if(isnan(S[2*ii-2])) {printf("%s %d","isNaN S ",2*ii-2);atexit(pexit);exit(0);} 
        if(isnan(S[2*ii-1])) {printf("%s %d","isNaN S ",2*ii-1);atexit(pexit);exit(0);} 
        }//END: POPUloop
    }
//--------------------------------------------------------
    if(RT_FLUX_SOLVER == 3)
    {
//- RT flux matrix ---------------------------------------
//compared to above methods, here it is more efficient to cast coefficients directly into matrix
        //This casts flux coeffiecients Xeq,Yeq,Zeq into the matrix
        // Matrix sequence: TOA Eq., Eq.25(0), Eq.26(1), Eq.25(1), Eq.26(2), Eq.25(2), ... Eq.25(2N-1), Eq.26(2N), Surf.Eq.(2N)
        // where the index (i) refers to index in Heng+2018 Eq.25&26: F_up/dn(i) = f(F_up/dn,B,w_0,g_0,etc.)

//printf("%s\n", "Filling RT flux matrix"); //template
        //TOA:
            //LHS: Flux coef:
            fcoef[1][1] = 1.0;//btridiag
            //RHS: Source term
            Scoef[1] = muD[1]*ffrac[1]*solar[lbin]; //incoming TOA flux
        //BOA:
            //LHS: Flux terms
            fcoef[2*zbin+2][2*zbin+1] = -PSURFAB;
            fcoef[2*zbin+2][2*zbin+2] = 1.0;
            //RHS: Source terms
            Scoef[2*zbin+2] =  PI* Bint; //testing

        //remainig matrix:
        for(ii=1; ii<=zbin; ii++) //POPUloop
        {
        //Eq.25
            //LHS: Flux terms
            fcoef[2*ii][2*ii-1] = Zeq[ii];
            fcoef[2*ii][2*ii] = Xeq[ii];
            fcoef[2*ii][2*ii+2] = Yeq[ii];
            //RHS: Source terms
            Scoef[2*ii] = PIt[ii] * ( BPlo[ii-1]*Xeq[ii] + BPup[ii]*Yeq[ii] + BMlo[ii-1]*Zeq[ii] ) + CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            Scoef[2*ii] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
        //Eq.26
            //LHS: Flux terms
            fcoef[2*ii+1][2*ii-1] = Yeq[ii];
            fcoef[2*ii+1][2*ii+1] = Xeq[ii];
            fcoef[2*ii+1][2*ii+2] = Zeq[ii];
            //RHS: Source terms
            Scoef[2*ii+1] = PIt[ii] * ( BMup[ii]*Xeq[ii] + BMlo[ii-1]*Yeq[ii] + BPup[ii]*Zeq[ii] ) + CP[ii]*TransC[ii]*Yeq[ii];
            Scoef[2*ii+1] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );

        //if(isnan(X[2*ii-2])) {printf("%s %d","isNaN X ",2*ii-2);atexit(pexit);exit(0);} 
        }//END: POPUloop
    }
//--------------------------------------------------------
    if(RT_FLUX_SOLVER == 14)
    {
//- RT flux matrix ---------------------------------------
//compared to above methods, here it is more efficient to cast coefficients directly into matrix
        //This casts flux coeffiecients Xeq,Yeq,Zeq into the matrix
        // Matrix sequence: TOA Eq., Eq.25(0), Eq.26(1), Eq.25(1), Eq.26(2), Eq.25(2), ... Eq.25(2N-1), Eq.26(2N), Surf.Eq.(2N)
        // where the index (i) refers to index in Heng+2018 Eq.25&26: F_up/dn(i) = f(F_up/dn,B,w_0,g_0,etc.)

//printf("%s\n", "Filling RT flux matrix"); //template
        //TOA:
            //LHS: Flux coef:
            AA[1] = 0.0;//penta
            BB[1] = 0.0;//penta
            CC[1] = 1.0;//penta
            DD[1] = 0.0;//penta
            EE[1] = 0.0;//penta
            //RHS: Source term
            Scoef[1] = muD[1]*ffrac[1]*solar[lbin]; //incoming TOA flux
        //BOA:
            //LHS: Flux terms
            AA[2*zbin+2] = 0.0;
            BB[2*zbin+2] = -PSURFAB;
            CC[2*zbin+2] = 1.0;
            DD[2*zbin+2] = 0.0;
            EE[2*zbin+2] = 0.0;
            //RHS: Source terms
            Scoef[2*zbin+2] =  PI* Bint; //testing

        //remainig matrix:
        for(ii=1; ii<=zbin; ii++) //POPUloop
        {
        //Eq.25
            //LHS: Flux terms
            AA[2*ii] = 0.0;
            BB[2*ii] = Zeq[ii];
            CC[2*ii] = Xeq[ii];
            DD[2*ii] = 0.0;
            EE[2*ii] = Yeq[ii];
            //RHS: Source terms
            Scoef[2*ii] = PIt[ii] * ( BPlo[ii-1]*Xeq[ii] + BPup[ii]*Yeq[ii] + BMlo[ii-1]*Zeq[ii] ) + CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            Scoef[2*ii] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
        //Eq.26
            //LHS: Flux terms
            AA[2*ii+1] = Yeq[ii];
            BB[2*ii+1] = 0.0;
            CC[2*ii+1] = Xeq[ii];
            DD[2*ii+1] = Zeq[ii];
            EE[2*ii+1] = 0.0;
            //RHS: Source terms
            Scoef[2*ii+1] = PIt[ii] * ( BMup[ii]*Xeq[ii] + BMlo[ii-1]*Yeq[ii] + BPup[ii]*Zeq[ii] ) + CP[ii]*TransC[ii]*Yeq[ii];
            Scoef[2*ii+1] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );

        //if(isnan(X[2*ii-2])) {printf("%s %d","isNaN X ",2*ii-2);atexit(pexit);exit(0);} 
        }//END: POPUloop
    }
//--------------------------------------------------------
    if(RT_FLUX_SOLVER >= 4)
    {
//- RT flux matrix ---------------------------------------
//compared to above methods, here it is more efficient to cast coefficients directly into matrix
        //This casts flux coeffiecients Xeq,Yeq,Zeq into the matrix
        // Matrix sequence: TOA Eq., Eq.25(0), Eq.26(1), Eq.25(1), Eq.26(2), Eq.25(2), ... Eq.25(2N-1), Eq.26(2N), Surf.Eq.(2N)
        // where the index (i) refers to index in Heng+2018 Eq.25&26: F_up/dn(i) = f(F_up/dn,B,w_0,g_0,etc.)

//printf("%s\n", "Filling RT flux matrix"); //template
        //TOA:
            //LHS: Flux coef:
            Emat[0] = 0.0;//penta
            Cmat[0] = 0.0;//penta
            Dmat[0] = 1.0;//penta
            Amat[0] = 0.0;//penta
            Bmat[0] = 0.0;//penta
            //RHS: Source term
            Smat[0] = muD[1]*ffrac[1]*solar[lbin]; //incoming TOA flux
        //BOA:
            //LHS: Flux terms
            Emat[2*zbin+1] = 0.0;
            //Cmat[2*zbin+1] = -PSURFAB;
            Cmat[2*zbin+1] = -1.0;
            Dmat[2*zbin+1] = 1.0;
            Amat[2*zbin+1] = 0.0;
            Bmat[2*zbin+1] = 0.0;
            //RHS: Source terms
            Smat[2*zbin+1] =  PI* Bint; //testing

        //remainig matrix:
        for(ii=1; ii<=zbin; ii++) //POPUloop
        {
        //Eq.25
            //LHS: Flux terms
            Emat[2*ii-1] = 0.0;
            Cmat[2*ii-1] = Zeq[ii];
            Dmat[2*ii-1] = Xeq[ii];
            Amat[2*ii-1] = 0.0;
            Bmat[2*ii-1] = Yeq[ii];
            //RHS: Source terms
            Smat[2*ii-1] = PIt[ii] * ( BPlo[ii-1]*Xeq[ii] + BPup[ii]*Yeq[ii] + BMlo[ii-1]*Zeq[ii] ) + CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            Smat[2*ii-1] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
        //Eq.26
            //LHS: Flux terms
            Emat[2*ii] = Yeq[ii];
            Cmat[2*ii] = 0.0;
            Dmat[2*ii] = Xeq[ii];
            Amat[2*ii] = Zeq[ii];
            Bmat[2*ii] = 0.0;
            //RHS: Source terms
            Smat[2*ii] = PIt[ii] * ( BMup[ii]*Xeq[ii] + BMlo[ii-1]*Yeq[ii] + BPup[ii]*Zeq[ii] ) + CP[ii]*TransC[ii]*Yeq[ii];
            Smat[2*ii] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );

        //if(isnan(X[2*ii-2])) {printf("%s %d","isNaN X ",2*ii-2);atexit(pexit);exit(0);} 
        }//END: POPUloop
    }



//==================================================================================================
//=== SOLVING RT Fluxes: various solvers available below ===========================================
//==================================================================================================
    if(RT_FLUX_SOLVER == 0)
    {
//=== SOLVING tri-diagonal matrix "Thomas Algorithm" (e.g. Toon+1989) ==============================
// Can only be used on Eq. sequence 26,25,26,25,...
//time = clock();
//clock_gettime(CLOCK_REALTIME, &start);
//printf("%s\n", "Solving matrix"); //template
        //--- Elimination ---
        Xdash[2*zbin] = X[2*zbin] / Y[2*zbin];
        Sdash[2*zbin] = S[2*zbin] / Y[2*zbin];
//            printf("%s%d%s%e\t","Xdash[",2*zbin,"] = ",Xdash[2*zbin]);
//            printf("%s%d%s%e\t","Sdash[",2*zbin,"] = ",Sdash[2*zbin]);
//            printf("%s%d%s%e\t","X[",2*zbin,"] = ",X[2*zbin]);
//            printf("%s%d%s%e\t","Z[",2*zbin,"] = ",Z[2*zbin]);
//            printf("%s%d%s%e\t","S[",2*zbin,"] = ",S[2*zbin]);
//            printf("%s%d%s%e\n","Y[",2*zbin,"] = ",Y[2*zbin]);
        if(isnan(Sdash[2*zbin])) {
            printf("%s %d","isNaN Sdash ",2*zbin);
            printf("%s%d%s%e","S[",2*zbin,"] = ",S[2*zbin]);
            printf("%s%d%s%e","Y[",2*zbin,"] = ",Y[2*zbin]);
            atexit(pexit);exit(0);} 
        if(isnan(Xdash[2*zbin])) {
            printf("%s %d","isNaN Xdash ",2*zbin);
            printf("%s%d%s%e","X[",2*zbin,"] = ",X[2*zbin]);
            printf("%s%d%s%e","Y[",2*zbin,"] = ",Y[2*zbin]);
            atexit(pexit);exit(0);} 
        for(ii=2*zbin-1; ii>=1; ii--) //Matrixloop BOA -> TOA
        {
            Xdash[ii] = X[ii] / ( Y[ii] - Z[ii]*Xdash[ii+1] );
            Sdash[ii] = ( S[ii] - Z[ii]*Sdash[ii+1] ) / ( Y[ii] - Z[ii]*Xdash[ii+1] );
//            printf("%s%d%s%e\t","Xdash[",ii,"] = ",Xdash[ii]);
//            printf("%s%d%s%e\t","Sdash[",ii,"] = ",Sdash[ii]);
//            printf("%s%d%s%e\t","X[",ii,"] = ",X[ii]);
//            printf("%s%d%s%e\t","Z[",ii,"] = ",Z[ii]);
//            printf("%s%d%s%e\t","S[",ii,"] = ",S[ii]);
//            printf("%s%d%s%e\n","Y[",ii,"] = ",Y[ii]);
        if(isnan(Xdash[ii])) {printf("%s %d","isNaN Xdash ",ii);atexit(pexit);exit(0);} 
        if(isnan(Sdash[ii])) {printf("%s %d","isNaN Sdash ",ii);atexit(pexit);exit(0);} 
        }//END: Matrixloop BOA -> TOA

        //--- Backward Substitution ---
        //Fmat is like Yl in Toon+89, but already as Fup and Fdn, in sequence:
        //Fmat[2n] = Fup[n] & Fmat[2n+1] = Fdn[n]
        Fmat[2] = Sdash[1];
        for(ii=2; ii<=2*zbin; ii++) //Matrixloop TOA -> BOA
        {
            Fmat[ii+1] = Sdash[ii] - Xdash[ii]*Fmat[ii];
        }//END: Matrixloop TOA -> BOA

    }
//--- HERE ends tri-diagonal matrix (Thomas algorithm) ---------------------------------------------
    if(RT_FLUX_SOLVER == 1)
    {
//=== TOA-BOA and back substitution BOA to TOA =====================================================
        //--- Elimination ---
        Zdash[1] = Z[1] / Y[1];
        Sdash[1] = S[1] / Y[1];
//            printf("%s%d%s%e\t","Zdash[",2*zbin,"] = ",Zdash[2*zbin]);
//            printf("%s%d%s%e\t","Sdash[",2*zbin,"] = ",Sdash[2*zbin]);
//            printf("%s%d%s%e\t","X[",2*zbin,"] = ",X[2*zbin]);
//            printf("%s%d%s%e\t","Z[",2*zbin,"] = ",Z[2*zbin]);
//            printf("%s%d%s%e\t","S[",2*zbin,"] = ",S[2*zbin]);
//            printf("%s%d%s%e\n","Y[",2*zbin,"] = ",Y[2*zbin]);
        if(isnan(Sdash[1])) {
            printf("%s %d","isNaN Sdash ",1);
            printf("%s%d%s%e","S[",1,"] = ",S[1]);
            printf("%s%d%s%e","Y[",1,"] = ",Y[1]);
            atexit(pexit);exit(0);} 
        if(isnan(Zdash[1])) {
            printf("%s %d","isNaN Zdash ",1);
            printf("%s%d%s%e","Z[",1,"] = ",Z[1]);
            printf("%s%d%s%e","Y[",1,"] = ",Y[1]);
            atexit(pexit);exit(0);} 
        for(ii=2; ii<=2*zbin; ii++) //Matrixloop BOA -> TOA
        {
            Zdash[ii] = Z[ii] / ( Y[ii] - X[ii]*Zdash[ii-1] );
            Sdash[ii] = ( S[ii] - X[ii]*Sdash[ii-1] ) / ( Y[ii] - X[ii]*Zdash[ii-1] );
//            printf("%s%d%s%e\t","Zdash[",ii,"] = ",Zdash[ii]);
//            printf("%s%d%s%e\t","Sdash[",ii,"] = ",Sdash[ii]);
//            printf("%s%d%s%e\t","X[",ii,"] = ",X[ii]);
//            printf("%s%d%s%e\t","Z[",ii,"] = ",Z[ii]);
//            printf("%s%d%s%e\t","S[",ii,"] = ",S[ii]);
//            printf("%s%d%s%e\n","Y[",ii,"] = ",Y[ii]);
        if(isnan(Zdash[ii])) {printf("%s %d","isNaN Zdash ",ii);atexit(pexit);exit(0);} 
        if(isnan(Sdash[ii])) {printf("%s %d","isNaN Sdash ",ii);atexit(pexit);exit(0);} 
        }//END: Matrixloop BOA -> TOA

        //--- Backward Substitution ---
        //Fmat is like Yl in Toon+89, but already as Fup and Fdn, in sequence:
        //Fmat[2n] = Fup[n] & Fmat[2n+1] = Fdn[n]
        Fmat[2*zbin+1] = Sdash[2*zbin];
        for(ii=2*zbin; ii>=2; ii--) //Matrixloop TOA -> BOA
        {
            Fmat[ii] = Sdash[ii-1] - Zdash[ii-1]*Fmat[ii+1];
        }//END: Matrixloop TOA -> BOA

    }
//--- HERE ends tri-diagonal matrix (Thomas algorithm) ---------------------------------------------
    if(RT_FLUX_SOLVER == 2)
    {
//=== SOLVING fluxes using LU decomposition (e.g. numerical recipes 2nd edition) ===================
//clock_gettime(CLOCK_REALTIME, &start);
    double **LU;
    int *indx;
    double bonus;
    double *LUsol;
    LU=dmatrix(1,2*zbin,1,2*zbin);
    indx=ivector(1,2*zbin);
    LUsol=dvector(1,2*zbin); //Input: S[1:2*zbin]; Output: Fmat[2:2*zbin+1] !!
// first put everything into LU matrix:
    //printf("%s%d%s%d\n","j=",j,"lbin = ",lbin);
    for (i=1;i<=2*zbin;i++) // rows
    {
        //if(j==0 && lbin==0) printf("%s%d%s","LU[",i,"] = ");
        LUsol[i] = S[i]; // fill first with RHS
        for (ii=1;ii<=2*zbin;ii++) //columns
        {
            if (ii == i-1) LU[i][ii] = X[i]; //lower diagonal
            else if (ii == i) LU[i][ii] = Y[i]; //main diagonal
            else if (ii == i+1) LU[i][ii] = Z[i]; //upper diagonal
            else LU[i][ii] = 0.0;
            //if(j==0 && lbin==0) printf("%e\t",LU[i][ii]);
        }
        //if(j==0 && lbin==0) printf("\n");
    }

    ludcmp(LU,2*zbin,indx,&bonus);
    if(j==-1 && lbin==0){
        for(i=1;i<=2*zbin;i++){ 
            printf("%s%d%s","LU[",i,"] = ");
            for(j=1;j<=2*zbin;j++) printf("%e\t",LU[i][j]);
            printf("\n");
        }
    }
    lubksb(LU,2*zbin,indx,LUsol); //LUsol overwritten with Fupdn solution

    for (i=1;i<=2*zbin;i++)
    {
        Fmat[i+1] = LUsol[i]; //copy solutions to appropriate index
        //if(j==0 && lbin==0)  printf("%s%d%s%e\n","Fmat[",i+1,"] = ",Fmat[i+1]);
    }

    free_dmatrix(LU,1,2*zbin,1,2*zbin);
    free_ivector(indx,1,2*zbin);
    free_dvector(LUsol,1,2*zbin);

    }
//--- HERE ends LU decomposition -------------------------------------------------------------------
    if(RT_FLUX_SOLVER == 12) //supposedly faster LU decomposition from Numerical Recipes 3rd edition
    {
//=== SOLVING fluxes using LU decomposition (e.g. numerical recipes 3rd edition) ===================
//!!not working yet: needs to be set up with struct LUdcmp !!!!!!!!!!!!!!!!!!!!!!!!
//clock_gettime(CLOCK_REALTIME, &start);
    double **LU;
    LU = dmatrix(0,2*zbin-1,0,2*zbin-1);
    int *indx;
    indx = ivector(0,2*zbin-1);
    double *LUsol;
    LUsol = dvector(0,2*zbin-1);
// first put everything into LU matrix:
    //printf("%s%d%s%d\n","j=",j,"lbin = ",lbin);
    for (i=0;i<2*zbin;i++) // rows
    {
        //if(j==0 && lbin==0) printf("%s%d%s","LU[",i,"] = ");
        LUsol[i] = S[i+1]; // fill first with RHS
        for (ii=0;ii<2*zbin;ii++) //columns
        {
            if (ii == i-1) LU[i][ii] = X[i+1]; //lower diagonal
            else if (ii == i) LU[i][ii] = Y[i+1]; //main diagonal
            else if (ii == i+1) LU[i][ii] = Z[i+1]; //upper diagonal
            else LU[i][ii] = 0.0;
            //if(j==0 && lbin==0) printf("%e\t",LU[i][ii]);
        }
        //if(j==0 && lbin==0) printf("\n");
    }

//    if(j==0 && lbin==0){
//        for(i=0;i<2*zbin;i++){ 
//            printf("%s%d%s","input LU[",i+1,"] = ");
//            for(j=0;j<2*zbin;j++) printf("%e\t",LU[i][j]);
//            printf("\n");
//        }
//    }
    ms_LUdcmp(LU,2*zbin,indx);
    if(j==0 && lbin==0){
        for(i=0;i<2*zbin;i++){ 
            printf("%s%d%s","output LU[",i+1,"] = ");
            for(j=0;j<2*zbin;j++) printf("%e\t",LU[i][j]);
            printf("\n");
        }
    }
    ms_LUsolve(LU,2*zbin,indx,LUsol); //LUsol overwritten with Fupdn solution

    for (i=1;i<=2*zbin;i++)
    {
        Fmat[i+1] = LUsol[i-1]; //copy solutions to appropriate index
        //if(j==0 && lbin==0)  printf("%s%d%s%e\n","Fmat[",i+1,"] = ",Fmat[i+1]);
    }

    free_dmatrix(LU,0,2*zbin-1,0,2*zbin-1);
    free_ivector(indx,0,2*zbin-1);
    free_dvector(LUsol,0,2*zbin-1);

    }
//--- HERE ends ms_flux solver ---------------------------------------------------------------------
    if(RT_FLUX_SOLVER == 3) 
    {
//--- SOLVING fluxes in block tridiagonal setup - Eq. sequence 25,26,25,26,... ---------------------
//printf("%s\n", "Solving RT flux matrix"); //template
//clock_gettime(CLOCK_REALTIME, &start);
    BTridiagonal(fcoef,Scoef,zbin+1,2);
//printf("%s\n", "Filling Fmat"); //template
    for (i=0;i<=2*zbin+1;i++)
    {
        Fmat[i] = Scoef[i+1]; //copy solutions to appropriate index
        //if(j==0 && lbin==0)  printf("%s%d%s%e\n","Fmat[",i,"] = ",Fmat[i]);
    }

    }
//--- HERE ends Block Tridiagonal solver -----------------------------------------------------------
    if(RT_FLUX_SOLVER == 14) 
    {
//--- SOLVING fluxes in pentadiagonal setup: Eq. sequence 25,26,25,26,... --------------------------
//printf("%s\n", "Solving RT flux matrix"); //template
//clock_gettime(CLOCK_REALTIME, &start);
//printf("%s\n", "Filling Fmat"); //template
    double AdivD;

    for (i=2;i<=2*zbin+1;i++)
    {
        AdivD = BB[i-1] / CC[i-1];
        CC[i] = CC[i] - AdivD * DD[i-1];
        DD[i] = DD[i] - AdivD * EE[i-1];
        Scoef[i] = Scoef[i] - AdivD * Scoef[i-1];
        
        AdivD = AA[i-1] / CC[i-1];
        BB[i] = BB[i] - AdivD * DD[i-1];
        CC[i+1] = CC[i+1] - AdivD * EE[i-1];
        Scoef[i+1] = Scoef[i+1] - AdivD * Scoef[i-1];
    }
        
    AdivD = BB[2*zbin+1] / CC[2*zbin+1];
    CC[2*zbin+2] = CC[2*zbin+2] - AdivD * DD[zbin+1];
    Fmat[2*zbin+1] = (Scoef[2*zbin+2] - AdivD * Scoef[2*zbin+1]) / CC[2*zbin+2];
    Fmat[2*zbin] = (Scoef[2*zbin+1] - DD[2*zbin+1] * Fmat[2*zbin+1]) / CC[2*zbin+1];
        
    for (i=2*zbin;i>=1;i--)
    {
        Fmat[i-1] = (Scoef[i] - EE[i] * Fmat[i+1] - DD[i] * Fmat[i]) / CC[i]; //copy solutions to appropriate index
        //if(j==0 && lbin==0)  printf("%s%d%s%e\n","Fmat[",i,"] = ",Fmat[i]);
    }

    }
//--- HERE ends Pentadiagonal solver ---------------------------------------------------------------
    if(RT_FLUX_SOLVER == 4) 
    {
//--- SOLVING fluxes with PTRANS-I: Eq. sequence 25,26,25,26,... -----------------------------------
//printf("%s\n", "Solving RT flux matrix"); //template
//clock_gettime(CLOCK_REALTIME, &start);
//printf("%s\n", "Filling Fmat"); //template

    //TOA:
    muPen[0] = Dmat[0]; 
    alphaPen[0] = Amat[0]/muPen[0];
    betaPen[0] = Bmat[0]/muPen[0];
    zPen[0] = Smat[0]/muPen[0];

    gammaPen[1] = Cmat[1];
    muPen[1] = Dmat[1] - alphaPen[0]*gammaPen[1];
    alphaPen[1] = (Amat[1] - betaPen[0]*gammaPen[1])/muPen[1];
    betaPen[1] = Bmat[1]/muPen[1];
    zPen[1] = (Smat[1] - zPen[0]*gammaPen[1])/muPen[1];

    //iteration:
    for (i=2;i<=2*zbin-1;i++)
    {
        gammaPen[i] = Cmat[i] - alphaPen[i-2]*Emat[i];
        muPen[i] = Dmat[i] - betaPen[i-2]*Emat[i] - alphaPen[i-1]*gammaPen[i];
        alphaPen[i] = (Amat[i] - betaPen[i-1]*gammaPen[i]) / muPen[i];
        betaPen[i] = Bmat[i] / muPen[i];
        zPen[i] = (Smat[i] - zPen[i-2]*Emat[i] - zPen[i-1]*gammaPen[i]) / muPen[i];
    }

    //BOA:
    gammaPen[2*zbin] = Cmat[2*zbin] - alphaPen[2*zbin-2]*Emat[2*zbin];
    muPen[2*zbin] = Dmat[2*zbin] - betaPen[2*zbin-2]*Emat[2*zbin] - alphaPen[2*zbin-1]*gammaPen[2*zbin];
    alphaPen[2*zbin] = (Amat[2*zbin] - betaPen[2*zbin-1]*gammaPen[2*zbin]) / muPen[2*zbin];
    zPen[2*zbin] = (Smat[2*zbin] - zPen[2*zbin-1]*Emat[2*zbin] - zPen[2*zbin-1]*gammaPen[2*zbin]) / muPen[2*zbin];


    gammaPen[2*zbin+1] = Cmat[2*zbin+1] - alphaPen[2*zbin-1]*Emat[2*zbin+1];
    muPen[2*zbin+1] = Dmat[2*zbin+1] - betaPen[2*zbin-1]*Emat[2*zbin+1] - alphaPen[2*zbin]*gammaPen[2*zbin+1];
    zPen[2*zbin+1] = (Smat[2*zbin+1] - zPen[2*zbin]*Emat[2*zbin+1] - zPen[2*zbin]*gammaPen[2*zbin+1]) / muPen[2*zbin+1];

    //Solution vector back-substitution:
    Fmat[2*zbin+1] = zPen[2*zbin+1];
    Fmat[2*zbin] = zPen[2*zbin] - alphaPen[2*zbin]*Fmat[2*zbin+1];

    for (i=2*zbin-1;i>=0;i--)
    {
        Fmat[i] = zPen[i] - alphaPen[i]*Fmat[i+1] - betaPen[i]*Fmat[i+2];
    }

    }
//--- HERE ends PTRANS-I solver --------------------------------------------------------------------
    if(RT_FLUX_SOLVER == 5) 
    {
//--- SOLVING fluxes with Sogabe-3: Eq. sequence 25,26,25,26,... -----------------------------------
//for simplicity, the symbolic substitution with lambda is not implemented. this may increase instability but implementing a symbolic algorithm seems excessive for our use
//printf("%s\n", "Solving RT flux matrix"); //template
//clock_gettime(CLOCK_REALTIME, &start);
//printf("%s\n", "Filling Fmat"); //template

    //TOA:
    cSog[0] = Dmat[0];
    eSog[0] = Amat[0];
    fSog[1] = Cmat[1] / cSog[0];
    cSog[1] = Dmat[1] - fSog[1]*Amat[0];

    //iteration:
    for (i=2;i<=2*zbin+1;i++)
    {
        eSog[i-1] = Amat[i-1] - fSog[i-1]*Bmat[i-2];
        gSog[i] = Emat[i] / cSog[i-2];
        fSog[i] = (Cmat[i] - gSog[i]*eSog[i-2]) / cSog[i-1];
        cSog[i] = Dmat[i] - fSog[i]*eSog[i-1] - gSog[i]*Bmat[i-2];
    }

    //Backward substitution:
    ytilSog[0] = Smat[0];
    ytilSog[1] = Smat[1] - fSog[1]*ytilSog[0];
    for (i=2;i<=2*zbin+1;i++)
    {
        ytilSog[i] = Smat[i] - gSog[i]*ytilSog[i-2] - fSog[i]*ytilSog[i-1];
    }

    //Forward substitution:
    Fmat[2*zbin+1] = ytilSog[2*zbin+1] / cSog[2*zbin+1];
    Fmat[2*zbin] = (ytilSog[2*zbin] - eSog[2*zbin]*Fmat[2*zbin+1]) / cSog[2*zbin];
    for (i=2*zbin-1;i<=0;i--)
    {
        Fmat[i] = (ytilSog[i] - Bmat[i]*Fmat[i+2] - eSog[i]*Fmat[i+1]) / cSog[i];
    }


    }
//--- HERE ends Sogabe-3 solver --------------------------------------------------------------------
    if(RT_FLUX_SOLVER == 6) 
    {
//--- SOLVING fluxes with KPENTA: Eq. sequence 25,26,25,26,... -------------------------------------
//for simplicity, the symbolic substitution with lambda is not implemented. this may increase instability but implementing a symbolic algorithm seems excessive for our use
//printf("%s\n", "Solving RT flux matrix"); //template
//clock_gettime(CLOCK_REALTIME, &start);
//printf("%s\n", "Filling Fmat"); //template

    //TOA:
    cKPen[0] = Dmat[0];
    eKPen[0] = Amat[0];
    vKPen[0] = 0.0;
    hKPen[0] = 0.0;
    fKPen[1] = Cmat[1] / cKPen[0];
    cKPen[1] = Dmat[1] - fKPen[1]*Amat[0];
    vKPen[1] = -fKPen[1]*vKPen[0];
    hKPen[1] = -eKPen[0]*hKPen[0] / cKPen[1];

    //iteration:
    for (i=2;i<=2*zbin;i++)
    {
        eKPen[i-1] = Amat[i-1] - fKPen[i-1] * Bmat[i-2];
        fKPen[i] = Cmat[i] / cKPen[i-1] - (Emat[i] / cKPen[i-2]) * (eKPen[i-2] / cKPen[i-1]);
        cKPen[i] = Dmat[i] - fKPen[i] * eKPen[i-1] - Emat[i] * Bmat[i-2] / cKPen[i-2];
    }

    for (i=2;i<=2*zbin-2;i++)
    {
        vKPen[i] = Emat[i] * vKPen[i-2] / cKPen[i-2] - fKPen[i] * vKPen[i-1];
        hKPen[i] = - Bmat[i-2] * hKPen[i-2] / cKPen[i] - eKPen[i-1] * hKPen[i-1] / cKPen[i];
    }
    vKPen[2*zbin-1] = Bmat[2*zbin-1] - Emat[2*zbin-1] * vKPen[2*zbin-3] / cKPen[2*zbin-3] - fKPen[2*zbin-1] * vKPen[2*zbin-2];
    vKPen[2*zbin] = Amat[2*zbin] - Emat[2*zbin] * vKPen[2*zbin-2] / cKPen[2*zbin-2] - fKPen[2*zbin] * vKPen[2*zbin-1];
    hKPen[2*zbin-1] = Emat[2*zbin+1] / cKPen[2*zbin-1] - Bmat[2*zbin-3] * hKPen[2*zbin-3] / cKPen[2*zbin-1] - eKPen[2*zbin-2] * hKPen[2*zbin-2] / cKPen[2*zbin-1];
    hKPen[2*zbin] = Cmat[2*zbin+1] / cKPen[2*zbin] - Bmat[2*zbin-2] * hKPen[2*zbin-2] / cKPen[2*zbin] - eKPen[2*zbin-1] * hKPen[2*zbin-1] / cKPen[2*zbin];

    cKPen[2*zbin+1] = Dmat[2*zbin+1];
    for (i=0;i<=2*zbin;i++) cKPen[2*zbin+1] -= hKPen[i] * vKPen[i];

    //Backward substitution:
    zKPen[0] = Smat[0];
    zKPen[1] = Smat[1] - fKPen[1]*zKPen[0];
    for (i=2;i<=2*zbin;i++)
    {
        zKPen[i] = Smat[i] - Emat[i] * zKPen[i-2] / cKPen[i-2] - fKPen[i] * zKPen[i-1];
    }

    zKPen[2*zbin+1] = Smat[2*zbin+1];
    for (i=0;i<=2*zbin;i++) zKPen[2*zbin+1] -= hKPen[i] * zKPen[i];

    //Forward substitution:
    Fmat[2*zbin+1] = zKPen[2*zbin+1] / cKPen[2*zbin+1];
    Fmat[2*zbin] = (zKPen[2*zbin] - vKPen[2*zbin] * Fmat[2*zbin+1]) / cKPen[2*zbin];
    Fmat[2*zbin-1] = (zKPen[2*zbin-1] - vKPen[2*zbin-1] * Fmat[2*zbin+1] - eKPen[2*zbin-1] * Fmat[2*zbin]) / cKPen[2*zbin-1];
    for (i=2*zbin-2;i<=0;i--)
    {
        Fmat[i] = (zKPen[i] - vKPen[i] * Fmat[2*zbin+1] - eKPen[i] * Fmat[i+1] - Bmat[i] * Fmat[i+2]) / cKPen[i];
    }


    }
//--- HERE ends Sogabe-3 solver --------------------------------------------------------------------


       //Now use Eq.25 to calculate Fup[0] i.e. Fmat[0] from Fdn[0] & Fup[1]:
        //NOt for re-cast equations in solver 3 or 4
    if(RT_FLUX_SOLVER < 3)   Fmat[0] = 1/Xeq[1] * (S[0] - Yeq[1]*Fmat[2] - Zeq[1]*Fmat[1]);

//clock_gettime(CLOCK_REALTIME, &end);
//t2_passed = ((double)end.tv_sec*1e9 + end.tv_nsec) - ((double)start.tv_sec*1e9 + start.tv_nsec);
//printf("%s %f\n","\nflux solver took [ms]:", t2_passed/1e6); 

        for (ii=0;ii<=2*zbin+1;ii++) if(isnan(Fmat[ii])) {printf("%s %d %d","isNaN Fmat[lbin][ii] ",lbin,ii);atexit(pexit);exit(0);} 

//atexit(pexit);exit(0); //ms debugging mode

//==================================================================================================
//printf("%s\n", "Calculate NetFlux"); //template
//=== Calculate NetFlux for radiative layers =======================================================
//==================================================================================================
    
    if(RT_FLUX_SOLVER < 3) //re-written Hang+18 equations 
    {
    //Sequence is Fup[0], Fdn[0], Fup[1], Fdn[1], ..., Fup[zbin], Fdn[zbin]
    rl = 0;
    NetFlux[rl][j] += (Fmat[0] - Fmat[1]) * ( wavelength[lbin+1] - wavelength[lbin] );
    if(j==0) Fup[0] += Fmat[0] * ( wavelength[lbin+1] - wavelength[lbin] );
    if(j==0) Fdn[0] += Fmat[1] * ( wavelength[lbin+1] - wavelength[lbin] );
    for(i=1; i<=zbin; i++)
    {
        if(isconv[zbin+1-i] == 0)
        {
            rl += 1;
            NetFlux[rl][j] += ( Fmat[2*i]  - Fmat[2*i+1] ) * ( wavelength[lbin+1] - wavelength[lbin] );
        }
    if(j==0) Fup[i] += Fmat[2*i] * ( wavelength[lbin+1] - wavelength[lbin] ); //fill ALL levels for diagnostics
    if(j==0) Fdn[i] += Fmat[2*i+1] * ( wavelength[lbin+1] - wavelength[lbin] );
    }
    }

    if(RT_FLUX_SOLVER >= 3) //orig. Heng+18 equations
    {
    //Sequence is Fdn[0], Fup[0], Fdn[1], Fup[1], ..., Fdn[zbin], Fup[zbin]
    rl = 0;
    NetFlux[rl][j] += (Fmat[1] - Fmat[0]) * ( wavelength[lbin+1] - wavelength[lbin] );
    if(j==0) Fup[0] += Fmat[1] * ( wavelength[lbin+1] - wavelength[lbin] );
    if(j==0) Fdn[0] += Fmat[0] * ( wavelength[lbin+1] - wavelength[lbin] );
    for(i=1; i<=zbin; i++)
    {
        if(isconv[zbin+1-i] == 0)
        {
            rl += 1;
            NetFlux[rl][j] += ( Fmat[2*i+1]  - Fmat[2*i] ) * ( wavelength[lbin+1] - wavelength[lbin] );
        }
    if(j==0) Fup[i] += Fmat[2*i+1] * ( wavelength[lbin+1] - wavelength[lbin] ); //fill ALL levels for diagnostics
    if(j==0) Fdn[i] += Fmat[2*i] * ( wavelength[lbin+1] - wavelength[lbin] );
    }
    }


//atexit(pexit);exit(0); //ms debugging mode

//DEBUGGING:

//printf("%s\n","\n==== DEBUGGING ===="); //template
//printf("%s\t%d\n","j = ",j); //template
//printf("%s\n","\n==== (i) X Y Z S ===="); //template
//for (i=1;i<=2*zbin;i++)
//{
//    printf("%d\t%.3e\t%.3e\t%.3e\t%.3e\n",i,X[i],Y[i],Z[i],S[i]);
//}
//printf("%s\n","\n==== Up/Dn FLUXES ===="); //template
//for (i=0;i<zbin;i++)
//{
//    printf("%s %d\t%.3e\t%.3e\n","Fup,Fdn ",i,Fmat[2*i],Fmat[2*i+1]);
//}
//printf("%s\n","\n==== NET FLUXES ===="); //template
//for (i=0;i<=nrl;i++)
//{
//    printf("%s %d\t%.3e\n","NetFlux ",i,NetFlux[i][j]);
//}
//atexit(pexit);exit(0); //ms debugging mode

    }//END: Tvar loop
    
    free_dmatrix(fcoef,1,2*zbin+2,1,2*zbin+2);
    free_dvector(Scoef,1,2*zbin+2);
    free_dvector(AA,1,2*zbin+2);
    free_dvector(BB,1,2*zbin+2);
    free_dvector(CC,1,2*zbin+2);
    free_dvector(DD,1,2*zbin+2);
    free_dvector(EE,1,2*zbin+2);



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
