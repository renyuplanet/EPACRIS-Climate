/*******************************************************************
 * File: ms_d2stream.c
 * Function: Solve for equilibrium 2-stream fluxes from Heng+2018, Eqs. 25&26 (doi:10.3847/1538-4365/aad199)
 * +++++ BUT with double grid for non-isothermal layers described in Malik+2019 & Deitrick+2022 +++++
 * +++++ uses only Penta solver (4) +++++
 * Original Author: Markus Scheucher (markus.scheucher@jpl.nasa.gov)
 * Version: v01, 2022: Original upload
 *******************************************************************
 * (All (Eq.xx) refer to Heng+2018 / Malik+2019 unless stated otherwise)
 * Eq.25 & 26 represent Fup & Fdn for a pair of atmospheric layers i and i-1(above).
 * (Eq.25) Fup(i) = 1/Xeq { -Yeq Fup(i+1) - Zeq Fdn(i) + PIt [ Xeq B+(i) + Yeq B+(i+1) + Zeq B-(i) ]
 *                  + Yeq C-TsTc + Tc [ zeta+ ( C+zeta- - C-zeta+ ) + zeta-T^2 ( C-zeta- - C+zeta+  )  ] }
 * (Eq.26) Fdn(i) = 1/Xeq  { -Yeq Fdn(i-1) - Zeq Fup(i) + PIt [ Xeq B-(i) + Yeq B-(i-1) + Zeq B+(i) ]
 *                  + Yeq C+Tc + TsTc [ zeta+ ( C-zeta- - C+zeta+ ) + zeta-T^2 ( C+zeta- - C-zeta+  )  ] }
*/                  

#include <math.h>
#include "constant.h"

void pexit(void) {printf("\n%s\n","=========================\n==== EXITING PROGRAM ====\n=========================\n");}

//=========================================================
//=== Function identifiers=================================
void ms_two_str_solver(int lbin, double w[], double g[], double tau[], int nrl, int isconv[], double **Tvar, int jTvar, double P[], double **NetFlux, double Fup[], double Fdn[], double Fcup[], double Fcdn[], double Tint);

//=========================================================
//=== Functions ===========================================
void ms_two_str_solver(int lbin, double w[], double g[], double tau[], int nrl, int isconv[], double **Tvar, int jTvar, double P[], double **NetFlux, double Fup[], double Fdn[], double Fcup[], double Fcdn[], double Tint)
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
//printf("%s\n", "ms_d2stream opened successfully");
//atexit(pexit);exit(0); //ms debugging mode
    double time,t_passed,t2_passed;
    struct timespec start,end;
//clock_gettime(CLOCK_REALTIME, &start);

    int i, j, ii, tv, rl;   //...loops: usually arrays are filled from 1 .. N (not from 0)
    double E[2*zbin+1];
    double mu[2*zbin+1],muD[2*zbin+1],ffrac[2*zbin+1];
    double zetap[2*zbin+1],zetam[2*zbin+1];
    double eps2[2*zbin+1],gammaA[2*zbin+1],gammaS[2*zbin+1],gammaD[2*zbin+1];
    double PIt[2*zbin+1];
    double CD[2*zbin+1],CP[2*zbin+1],CM[2*zbin+1];
    double dtau[2*zbin+1],Trans[2*zbin+1],TransD[2*zbin+1],TransC[2*zbin+1],tauC[2*zbin+1];
    double B[2*zbin+1],dB[2*zbin+1],Bint; //helper terms
    double Xeq[2*zbin+1],Yeq[2*zbin+1],Zeq[2*zbin+1]; //X,Y,Z in Eq25&26 description above
    double Fmat[4*zbin+2]; //matrix solution vector Y in Toon+89

//- interpolation for double grid:
    double wdoub[2*zbin+1],gdoub[2*zbin+1],Tdoub[2*zbin+1];
    
//- FLux Matrix ------------------------------------------
    double Emat[4*zbin+2], Cmat[4*zbin+2], Dmat[4*zbin+2], Amat[4*zbin+2], Bmat[4*zbin+2], Smat[4*zbin+2]; //for pentadiagonal solvers
    double alphaPen[4*zbin+2], betaPen[4*zbin+2], gammaPen[4*zbin+2], muPen[4*zbin+2], zPen[4*zbin+2]; //for PTRANS-I

//printf("%s\n", "ms_d2stream variables declared");
//=========================================================
//=== INITIALIZE variables ================================
    wdoub[1] = w[1];
    gdoub[1] = g[1];
    wdoub[2*zbin] = w[zbin];
    gdoub[2*zbin] = g[zbin];
    for (i=1;i<zbin;i++) //double grid interpolation
    {
        wdoub[2*i] = 2./3.*w[i] + 1./3.*w[i+1];
        wdoub[2*i+1] = 1./3.*w[i] + 2./3.*w[i+1];
        gdoub[2*i] = 2./3.*g[i] + 1./3.*g[i+1];
        gdoub[2*i+1] = 1./3.*g[i] + 2./3.*g[i+1];
    }
//printf("%s\n", "ms_d2stream wdoub and gdoub done");
//=========================================================
//printf("%s\n", "POPULATE helper terms"); //template
//=== POPULATE helper terms (mostly Heng+2018) ============
    tauC[0] = 0.0;  //initializing 'TOA'
    for (i=1;i<=2*zbin;i++)   //Layloop (through layers)
    {
        if(wdoub[i] >= 1.0) printf("WARNING: w[%d] >= 1.0 at lambda: %f\n", i, wavelength[lbin]); //template
    //E-factor (Eq.31)
        E[i] = 1.225 - 0.1582*gdoub[i] - 0.1777*wdoub[i] - 0.07465*gdoub[i]*gdoub[i] + 0.2351*wdoub[i]*gdoub[i] - 0.05582*wdoub[i]*wdoub[i];
        //E[i] = 1.0; //testing
        if(wdoub[i]<0.1 || E[i]<1.0)
        {
            E[i] = 1.0;
        }
        //E[i] = 1.0; //ms23 RC testing
    
    //angular parameters
        muD[i] = cos(THETAREF); //...for now constant and the same
        ffrac[i] = FADV/muD[i]; //...depends on the definition of FADV, since FADV*muD = 1/4 for global energy balance
        //ffrac[i] = FADV; //ms05/22 testing
        //for numerical stability
        if(muD[i] == 0.5*sqrt(E[i]*(E[i]-wdoub[i])*(1-wdoub[i]*gdoub[i])) ) 
        {    
            muD[i] *= 0.99;
            printf("%s\n","WARNING: muD adjustment actually happened!");
        }

    //Coupling Coefficients (Eq.24)
        zetap[i] = 0.5*(1+sqrt( (E[i]-wdoub[i]) / (E[i]*(1-wdoub[i]*gdoub[i])) ));
        zetam[i] = 0.5*(1-sqrt( (E[i]-wdoub[i]) / (E[i]*(1-wdoub[i]*gdoub[i])) ));

    //2-stream closures
       // eps2[i] = 0.5*(1-sqrt(3)*gdoub[i]*wdoub[i]);    //...2nd Eddington coefficient (Toon+89, Tbl.1 quadrature closure)
       // eps2[i] = 1.0/sqrt(3);    //... quadrature closure (Eq.21)
        eps2[i] = 2.0/3.0; // Eddington closure
       // eps2[i] = 0.5; // deep atmosphere - 100% absorption approximation
        gammaA[i] = 2*E[i] - wdoub[i]*(1+E[i]*gdoub[i]); //...absorbtion (Eq.20a)
        gammaS[i] = wdoub[i]*(1-E[i]*gdoub[i]);          //...scattering (Eq.20b)
        gammaD[i] = wdoub[i]*gdoub[i]*muD[i]*ffrac[i]*solar[lbin] / eps2[i]; //...direct beam (Eq.18)

    //PI term (Eq.27a)
        PIt[i] = PI*( (1-wdoub[i]) / (E[i]-wdoub[i]) );  

    //Direct beam coefficients
        CD[i] = wdoub[i]*ffrac[i]*solar[lbin]*( 2*E[i]*(1-wdoub[i]*gdoub[i]) + gdoub[i]/eps2[i] );  //... see C* (Eq.23)
        CD[i] /= 4*E[i]*( E[i] - wdoub[i]  )*( 1 - wdoub[i]*gdoub[i] ) - 1/( muD[i]*muD[i] );
        CP[i] = 0.5* ( CD[i] + ( CD[i]/muD[i] + gammaD[i] ) / ( 2*E[i]*(1 - wdoub[i]*gdoub[i]) ) ); //... (Eq.27b)
        CM[i] = 0.5* ( CD[i] - ( CD[i]/muD[i] + gammaD[i] ) / ( 2*E[i]*(1 - wdoub[i]*gdoub[i]) ) ); //... (Eq.27b)

    //Transmission functions (Eq.27c-f)
        tauC[i] = tauC[i-1] + TAUdoub[i];   //... tau above (cumulative)
        Trans[i] = exp( -TAUdoub[i]*sqrt( (gammaA[i]+gammaS[i]) * (gammaA[i]-gammaS[i]) ) );     //... T
        TransD[i] = exp( -TAUdoub[i] / muD[i]  );    //.. .T*
        TransC[i] = exp( -tauC[i-1] / muD[i]  );    //... Tabove

    //Abbreviations in Eq25&26
        Xeq[i] = ( zetam[i]*zetam[i]*Trans[i]*Trans[i] - zetap[i]*zetap[i] );
        Yeq[i] = Trans[i] * ( zetam[i]*zetam[i] - zetap[i]*zetap[i] );
        Zeq[i] = zetam[i]*zetap[i] * ( 1 - Trans[i]*Trans[i] );
    }//END: Layloop

//printf("%s\n", "ms_d2stream helper functions done");
//=========================================================
//printf("%s\n", "TVar loop"); //template
//=== HERE starts the temperature variation loop ==========

    for(j=0; j<=jTvar; j++) //Tvar loop
    {  

//    printf("%s %d\n", "ms_d2stream jTvar =",j);
        //Tdoub[0] = Tvar[0][j];
        for (i=0;i<=zbin;i++) //initialize temperature vector: cell interfaces
        {
            Tdoub[2*i] = Tvar[i][j];
        }
        for (i=1;i<2*zbin;i=i+2) //initialize temperature vector: cell centers
        {
            Tdoub[i] = (Tdoub[i-1] + Tdoub[i+1]) * 0.5;
        }
//printf("%s\n", "ms_d2stream Tdoub written");
//for (i=0;i<=2*zbin;i++) printf("%s%d%s %f\n", "Tdoub[",i,"]=",Tdoub[i]);
//atexit(pexit);exit(0); //ms debugging mode
//=========================================================
//printf("%s\n", "CALCULATE Planck"); //template
//=== Calculate Planck functions and derivatives ==========
        for(i=0; i<=2*zbin; i++) //Layloop 2; Thermal atmospheric radiation
        {
            B[i] = 2*HPLANCK*CLIGHT*CLIGHT / pow(wavelength[lbin]*1.0E-9,5) * 1.0E-9 ; //convert B into per 1 nm
            B[i] /= exp(HPLANCK*CLIGHT / wavelength[lbin]/1.0E-9/KBOLTZMANN/Tdoub[i]) - 1.0 ;
            if(i>0) dB[i] = (B[i] - B[i-1]) / TAUdoub[i] / ( 2*E[i]*(1 - wdoub[i]*gdoub[i])  ); //assumes opacities const. over Tvar range... 
        }//END: Layloop 2
        //Internal Temperature
            Bint = 2*HPLANCK*CLIGHT*CLIGHT / pow(wavelength[lbin]*1.0E-9,5) * 1.0E-9 ;
            Bint /= exp(HPLANCK*CLIGHT / wavelength[lbin]/1.0E-9/KBOLTZMANN/Tint) - 1.0 ;
            
//printf("%s\n", "ms_d2stream Planck functions done");
//==================================================================================================
//=== POPULATE flux matrix coefficients ============================================================
//==================================================================================================
//- RT flux matrix ---------------------------------------
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
            Emat[4*zbin+1] = 0.0;
            //Cmat[2*zbin+1] = 1.0 - (2*PSURFAB); //ms2023: well, but actually only the beam term, not the whole flux...
            // Cmat[4*zbin+1] = -1.0;
            Cmat[4*zbin+1] = 0.0; /* RYH testing */
            Dmat[4*zbin+1] = 1.0;
            Amat[4*zbin+1] = 0.0;
            Bmat[4*zbin+1] = 0.0;
            //RHS: Source terms
            //Smat[2*zbin+1] = PSURFAB * ( muD[zbin]*ffrac[zbin]*solar[lbin]*TransC[zbin] ) + PSURFEM * PIt[zbin] * B[zbin]; //testing
            //Smat[2*zbin+1] =  PSURFEM * PI * B[zbin]; //testing
            //Smat[2*zbin+1] = PSURFAB * ( muD[zbin]*ffrac[zbin]*solar[lbin]*TransC[zbin] ) + PI * Bint; //testing
            // Smat[4*zbin+1] =  PI * Bint + muD[2*zbin]*ffrac[1]*solar[lbin] * TransC[2*zbin]*TransD[2*zbin]; //testing
            Smat[4*zbin+1] =  PSURFEM * PI * B[2*zbin] + PSURFAB * muD[2*zbin]*ffrac[1]*solar[lbin] * TransC[2*zbin]*TransD[2*zbin]; //RYH testing

        //remainig matrix:
        for(ii=1; ii<=2*zbin; ii++) //POPUloop
        {
        //Heng Eq.25 / Deitrick Eq.4: F_up [ii-1]
            //LHS: Flux terms
            Emat[2*ii-1] = 0.0;
            Cmat[2*ii-1] = Zeq[ii];
            Dmat[2*ii-1] = Xeq[ii];
            Amat[2*ii-1] = 0.0;
            Bmat[2*ii-1] = -Yeq[ii];
            //RHS: Source terms Heng:
            if(TAUdoub[ii] >= 1.0e-4) // non-isothermal treatment
            {
            Smat[2*ii-1] = PIt[ii] * (B[ii-1] * (Xeq[ii]+Zeq[ii]) - B[ii] * Yeq[ii] + dB[ii] * (Xeq[ii] - Yeq[ii] - Zeq[ii]) ) - CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            Smat[2*ii-1] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
            } else { // isothermal treatment
            Smat[2*ii-1] = PIt[ii] * ((B[ii-1]+B[ii])/2.0 * (Xeq[ii]+Zeq[ii]-Yeq[ii])) - CM[ii]*TransD[ii]*TransC[ii]*Yeq[ii]; 
            Smat[2*ii-1] += TransC[ii] * ( zetap[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) );
            }
        //Heng Eq.26 / Deitrick Eq.6: F_dn [ii]
            //LHS: Flux terms
            Emat[2*ii] = -Yeq[ii];
            Cmat[2*ii] = 0.0;
            Dmat[2*ii] = Xeq[ii];
            Amat[2*ii] = Zeq[ii];
            Bmat[2*ii] = 0.0;
            //RHS: Source terms Heng:
            if(TAUdoub[ii] >= 1.0e-4) // non-isothermal treatment
            {
            Smat[2*ii] = PIt[ii] * (B[ii] * (Xeq[ii]+Zeq[ii]) - B[ii-1] * Yeq[ii] - dB[ii] * (Xeq[ii] - Yeq[ii] - Zeq[ii]) ) - CP[ii]*TransC[ii]*Yeq[ii];
            Smat[2*ii] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );
            } else { // isothermal treatment
            Smat[2*ii] = PIt[ii] * ((B[ii]+B[ii-1])/2.0 * (Xeq[ii]+Zeq[ii]-Yeq[ii])) - CP[ii]*TransC[ii]*Yeq[ii];
            Smat[2*ii] += TransD[ii]*TransC[ii] * ( zetap[ii]*( CM[ii]*zetam[ii]-CP[ii]*zetap[ii] ) + zetam[ii]*Trans[ii]*Trans[ii]*( CP[ii]*zetam[ii]-CM[ii]*zetap[ii] ) );
            }

        }//END: POPUloop


//printf("%s\n", "ms_d2stream Flux Matrix populated");

//==================================================================================================
//=== SOLVING RT Fluxes: various solvers available below ===========================================
//==================================================================================================
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
    for (i=2;i<=4*zbin-1;i++)
    {
        gammaPen[i] = Cmat[i] - alphaPen[i-2]*Emat[i];
        muPen[i] = Dmat[i] - betaPen[i-2]*Emat[i] - alphaPen[i-1]*gammaPen[i];
        alphaPen[i] = (Amat[i] - betaPen[i-1]*gammaPen[i]) / muPen[i];
        betaPen[i] = Bmat[i] / muPen[i];
        zPen[i] = (Smat[i] - zPen[i-2]*Emat[i] - zPen[i-1]*gammaPen[i]) / muPen[i];
    }

    //BOA:
    gammaPen[4*zbin] = Cmat[4*zbin] - alphaPen[4*zbin-2]*Emat[4*zbin];
    muPen[4*zbin] = Dmat[4*zbin] - betaPen[4*zbin-2]*Emat[4*zbin] - alphaPen[4*zbin-1]*gammaPen[4*zbin];
    alphaPen[4*zbin] = (Amat[4*zbin] - betaPen[4*zbin-1]*gammaPen[4*zbin]) / muPen[4*zbin];
    zPen[4*zbin] = (Smat[4*zbin] - zPen[4*zbin-1]*Emat[4*zbin] - zPen[4*zbin-1]*gammaPen[4*zbin]) / muPen[4*zbin];


    gammaPen[4*zbin+1] = Cmat[4*zbin+1] - alphaPen[4*zbin-1]*Emat[4*zbin+1];
    muPen[4*zbin+1] = Dmat[4*zbin+1] - betaPen[4*zbin-1]*Emat[4*zbin+1] - alphaPen[4*zbin]*gammaPen[4*zbin+1];
    zPen[4*zbin+1] = (Smat[4*zbin+1] - zPen[4*zbin]*Emat[4*zbin+1] - zPen[4*zbin]*gammaPen[4*zbin+1]) / muPen[4*zbin+1];

    //Solution vector back-substitution:
    Fmat[4*zbin+1] = zPen[4*zbin+1];
    Fmat[4*zbin] = zPen[4*zbin] - alphaPen[4*zbin]*Fmat[4*zbin+1];

    for (i=4*zbin-1;i>=0;i--)
    {
        Fmat[i] = zPen[i] - alphaPen[i]*Fmat[i+1] - betaPen[i]*Fmat[i+2];
    }

//printf("%s\n", "ms_d2stream Fluxes solved");
//--- HERE ends PTRANS-I solver --------------------------------------------------------------------


//clock_gettime(CLOCK_REALTIME, &end);
//t2_passed = ((double)end.tv_sec*1e9 + end.tv_nsec) - ((double)start.tv_sec*1e9 + start.tv_nsec);
//printf("%s %f\n","\nflux solver took [ms]:", t2_passed/1e6); 

        for (ii=0;ii<=4*zbin+1;ii++) if(isnan(Fmat[ii])) {printf("%s %d %d","isNaN Fmat[lbin][ii] ",lbin,ii);atexit(pexit);exit(0);} 

//atexit(pexit);exit(0); //ms debugging mode

//==================================================================================================
//printf("%s\n", "Calculate NetFlux"); //template
//=== Calculate NetFlux for radiative layers =======================================================
//==================================================================================================
    

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
            //NetFlux[rl][j] += ( Fmat[4*i+1]  - Fmat[4*i] ) * ( wavelength[lbin+1] - wavelength[lbin] );
            NetFlux[rl][j] += ( Fmat[4*i-1]  - Fmat[4*i-2] ) * ( wavelength[lbin+1] - wavelength[lbin] ); //let's use center fluxes instead
        }
    if(j==0) Fup[i] += Fmat[4*i+1] * ( wavelength[lbin+1] - wavelength[lbin] ); //fill ALL levels for diagnostics
    if(j==0) Fdn[i] += Fmat[4*i] * ( wavelength[lbin+1] - wavelength[lbin] );
    if(j==0) Fcup[i] += Fmat[4*i-1] * ( wavelength[lbin+1] - wavelength[lbin] ); //fill ALL levels for diagnostics
    if(j==0) Fcdn[i] += Fmat[4*i-2] * ( wavelength[lbin+1] - wavelength[lbin] );
    }

//printf("%s\n", "ms_d2stream Netflux solved");

//atexit(pexit);exit(0); //ms debugging mode


    }//END: Tvar loop
    
//printf("%s\n", "ms_d2stream e.o. Tvar loop");

//clock_gettime(CLOCK_REALTIME, &end);
//t2_passed = ((double)end.tv_sec*1e9 + end.tv_nsec) - ((double)start.tv_sec*1e9 + start.tv_nsec);
//printf("\n%s %f\n","ms_d2stream took [ms]:", t2_passed/1e6); 
//atexit(pexit);exit(0); //ms debugging mode


//printf("%s\t%e\n","MuD = cos(THETAREF) = ",muD[0]); //check
}
// END: void ms_two_str_solver()

//atexit(pexit);exit(0); //ms debugging mode
