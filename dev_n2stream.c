/*******************************************************************
 * File: ms_2stream.c
 * Function: Solve for equilibrium 2-stream fluxes from Heng+2018, Eqs. 25&26 (doi:10.3847/1538-4365/aad199)
 * +++++ BUT in the form written out in Malik+2019 & Deitrick+2022 +++++
 * Original Author: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Version: v01, 2022: Original upload
 *******************************************************************
 * (All (Eq.xx) refer to Malik+2019 unless stated otherwise)
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
    double time,t_passed,t2_passed;
    struct timespec start,end;
    
//clock_gettime(CLOCK_REALTIME, &start);
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

//=========================================================
//printf("%s\n", "POPULATE helper terms"); //template
//=== POPULATE helper terms (mostly Heng+2018) ============
    tauC[0] = 0.0;  //initializing 'TOA'
    for (i=1;i<=zbin;i++)   //Layloop (through layers)
    {
        if(w[i] == 1.0) printf("%s %d\n", "WARNING: w[",i,"] = 1.0 at lambda: ", wavelength[lbin]); //template
    //E-factor (Eq.31)
        E[i] = 1.225 - 0.1582*g[i] - 0.1777*w[i] - 0.07465*g[i]*g[i] + 0.2351*w[i]*g[i] - 0.05582*w[i]*w[i];
        //E[i] = 1.0; //testing
        if(E[i]<1.0)
        {
            E[i] = 1.0;
        }
    
    //angular parameters
        mu[i] = cos(THETAREF);  //...for now constant
        muD[i] = cos(THETAREF); //...for now constant and the same
        //ffrac[i] = FADV/muD[i]; //...depends on the definition of FADV, since FADV*muD = 1/4 for global energy balance
        ffrac[i] = FADV; //ms05/22 testing
        //for numerical stability
        if(muD[i] == 0.5*sqrt(E[i]*(E[i]-w[i])*(1-w[i]*g[i])) ) 
        {    
            muD[i] *= 0.99;
            printf("%s\n","WARNING: muD adjustment actually happened!");
        }

    //Coupling Coefficients (Eq.24)
        zetap[i] = 0.5*(1+sqrt( (E[i]-w[i]) / (E[i]*(1-w[i]*g[i])) ));
        zetam[i] = 0.5*(1-sqrt( (E[i]-w[i]) / (E[i]*(1-w[i]*g[i])) ));

    //2-stream closures
       // eps2[i] = 0.5*(1-sqrt(3)*g[i]*w[i]);    //...2nd Eddington coefficient (Toon+89, Tbl.1 quadrature closure)
       // eps2[i] = 1.0/sqrt(3);    //... quadrature closure (Eq.21)
        eps2[i] = 2.0/3.0; // Eddington closure
        gammaA[i] = 2*E[i] - w[i]*(1+E[i]*g[i]); //...absorbtion (Eq.20a)
        gammaS[i] = w[i]*(1-E[i]*g[i]);          //...scattering (Eq.20b)
        gammaD[i] = w[i]*g[i]*muD[i]*ffrac[i]*solar[lbin] / eps2[i]; //...direct beam (Eq.18)

    //PI term (Eq.27a)
        PIt[i] = PI*( (1-w[i]) / (E[i]-w[i]) );  

    //Direct beam coefficients
        CD[i] = w[i]*ffrac[i]*solar[lbin]*( 2*E[i]*(1-w[i]*g[i]) + g[i]/eps2[i] );  //... see C* (Eq.23)
        CD[i] /= 4*E[i]*( E[i] - w[i]  )*( 1 - w[i]*g[i] ) - 1/( muD[i]*muD[i] );
        CP[i] = 0.5* ( CD[i] + ( CD[i]/muD[i] + gammaD[i] ) / ( 2*E[i]*(1 - w[i]*g[i]) ) ); //... (Eq.27b)
        CM[i] = 0.5* ( CD[i] - ( CD[i]/muD[i] + gammaD[i] ) / ( 2*E[i]*(1 - w[i]*g[i]) ) ); //... (Eq.27b)

    //Transmission functions (Eq.27c-f)
        tauC[i] = tauC[i-1] + tau[i];   //... tau above (cumulative)
        Trans[i] = exp( -tau[i]*sqrt( (gammaA[i]+gammaS[i]) * (gammaA[i]-gammaS[i]) ) );     //... T
        TransD[i] = exp( -tau[i] / muD[i]  );    //.. .T*
        TransC[i] = exp( -tauC[i-1] / muD[i]  );    //... Tabove

    //Abbreviations in Eq25&26
        Xeq[i] = ( zetam[i]*zetam[i]*Trans[i]*Trans[i] - zetap[i]*zetap[i] );
        Yeq[i] = Trans[i] * ( zetam[i]*zetam[i] - zetap[i]*zetap[i] );
        Zeq[i] = zetam[i]*zetap[i] * ( 1 - Trans[i]*Trans[i] );
    }//END: Layloop

//=========================================================
//printf("%s\n", "TVar loop"); //template
//=== HERE starts the temperature variation loop ==========

    for(j=0; j<=jTvar; j++) //Tvar loop
    {  

//=========================================================
//printf("%s\n", "CALCULATE Planck"); //template
//=== Calculate Planck functions and derivatives ==========
        for(i=0; i<=zbin; i++) //Layloop 2; Thermal atmospheric radiation
        {
            B[i] = 2*HPLANCK*CLIGHT*CLIGHT / pow(wavelength[lbin]*1e-9,5) * 1e-9 ; //convert B into per 1 nm
            B[i] /= exp(HPLANCK*CLIGHT / wavelength[lbin]/1e-9/KBOLTZMANN/Tvar[i][j]) - 1.0 ;
            if(i>0)
            {
                if(tau[i] > 1.0e-4) dB[i] = (B[i] - B[i-1]) / tau[i] / ( 2*E[i]*(1 - w[i]*g[i])  ); //assumes opacities const. over Tvar range... 
                else 
                {
                    dB[i] = 0.0; //isothermal approach where no absorption
                }
            }
        }//END: Layloop 2
        //Internal Temperature
            Bint = 2*HPLANCK*CLIGHT*CLIGHT / pow(wavelength[lbin]*1e-9,5) * 1e-9 ;
            Bint /= exp(HPLANCK*CLIGHT / wavelength[lbin]/1e-9/KBOLTZMANN/Tint) - 1.0 ;
            
            if(lbin==0) sumBint = 0.0;
            if(j==0) sumBint += Bint*(wavelength[lbin+1]-wavelength[lbin])*1e-9;

        for(i=0; i<=zbin; i++) //Layloop 3; setting up planck related functions (Eq.27g)
        {
            if(i>0)
            {
                BPup[i] = B[i] + dB[i];
                BMup[i] = B[i] - dB[i];
            }
            if(i<zbin)
            {
                BPlo[i] = B[i] + dB[i+1];
                BMlo[i] = B[i] - dB[i+1];
            }
        }//END: Layloop 3

//==================================================================================================
//=== POPULATE flux matrix coefficients ============================================================
//==================================================================================================
//printf("%s\n", "POPULATE matrix"); //template
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
            //Smat[2*zbin+1] = PSURFAB * ( muD[zbin]*ffrac[zbin]*solar[lbin]*TransC[zbin] ) + PSURFEM * PIt[zbin] * B[zbin]; //testing
            //Smat[2*zbin+1] =  PSURFEM * PI * B[zbin]; //testing
            //Smat[2*zbin+1] = PSURFAB * ( muD[zbin]*ffrac[zbin]*solar[lbin]*TransC[zbin] ) + PI * Bint; //testing
            Smat[2*zbin+1] =  PI * Bint; //testing

        //remainig matrix:
        for(ii=1; ii<=zbin; ii++) //POPUloop
        {
        //Eq.25: F_up [ii-1]
            //LHS: Flux terms
            Emat[2*ii-1] = 0.0;
            Cmat[2*ii-1] = Zeq[ii];
            Dmat[2*ii-1] = Xeq[ii];
            Amat[2*ii-1] = 0.0;
            Bmat[2*ii-1] = -Yeq[ii];
            //RHS: Source terms
            Smat[2*ii-1] = PIt[ii] * (B[ii-1] * (Xeq[ii]+Zeq[ii]) - B[ii] * Yeq[ii] + dB[ii] * (Xeq[ii] - Yeq[ii] - Zeq[ii]) ) - CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            Smat[2*ii-1] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
        //Eq.26: F_dn [ii]
            //LHS: Flux terms
            Emat[2*ii] = -Yeq[ii];
            Cmat[2*ii] = 0.0;
            Dmat[2*ii] = Xeq[ii];
            Amat[2*ii] = Zeq[ii];
            Bmat[2*ii] = 0.0;
            //RHS: Source terms
            Smat[2*ii] = PIt[ii] * (B[ii] * (Xeq[ii] + Zeq[ii]) - B[ii-1] * Yeq[ii] - dB[ii] * (Xeq[ii] - Yeq[ii] - Zeq[ii]) ) - CP[ii]*TransC[ii]*Yeq[ii];
            Smat[2*ii] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );

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


//printf("%s\t%e\n","Mu = cos(THETAREF) = ",mu[0]); //check
}
// END: void ms_two_str_solver()

//atexit(pexit);exit(0); //ms debugging mode
