/*----------------------- global.h --------------------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: July 20, 2011

--------------------------------------------------------------------- */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#include "config.h"

/*---- External Variables ------------------------------------------- */

extern double meanmolecular[];
extern double zl[];
extern double pl[];
extern double tl[];
extern double MM[], MMZ[];
extern double wavelength[];
extern double solar[];
extern double crossr[], crossa[3][NLAMBDA], sinab[3][NLAMBDA], asym[3][NLAMBDA];
extern double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO; 
extern double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4; 
extern double **opacNH3;
extern double **opacC2H2, **opacC2H4, **opacC2H6, **opacHCN, **opacCH2O2, **opacHNO3;
extern double **opacN2O, **opacN2, **opacNO, **opacNO2, **opacOCS;
extern double **opacHF, **opacHCl, **opacHBr, **opacHI, **opacClO, **opacHClO;
extern double **opacHBrO, **opacPH3, **opacCH3Cl, **opacCH3Br, **opacDMS, **opacCS2;
extern double MeanCO2[zbin+1], MeanO2[zbin+1], MeanSO2[zbin+1], MeanH2O[zbin+1], MeanOH[zbin+1], MeanH2CO[zbin+1];
extern double MeanH2O2[zbin+1], MeanHO2[zbin+1], MeanH2S[zbin+1], MeanCO[zbin+1], MeanO3[zbin+1], MeanCH4[zbin+1];
extern double MeanNH3[zbin+1];	
extern double MeanC2H2[zbin+1], MeanC2H4[zbin+1], MeanC2H6[zbin+1], MeanHCN[zbin+1], MeanCH2O2[zbin+1], MeanHNO3[zbin+1];
extern double MeanN2O[zbin+1], MeanN2[zbin+1], MeanNO[zbin+1], MeanNO2[zbin+1], MeanOCS[zbin+1];
extern double SMeanCO2[zbin+1], SMeanO2[zbin+1], SMeanSO2[zbin+1], SMeanH2O[zbin+1], SMeanOH[zbin+1], SMeanH2CO[zbin+1];
extern double SMeanH2O2[zbin+1], SMeanHO2[zbin+1], SMeanH2S[zbin+1], SMeanCO[zbin+1], SMeanO3[zbin+1], SMeanCH4[zbin+1];
extern double SMeanNH3[zbin+1];
extern double SMeanC2H2[zbin+1], SMeanC2H4[zbin+1], SMeanC2H6[zbin+1], SMeanCH2O2[zbin+1];
extern double SMeanHCN[zbin+1], SMeanN2O[zbin+1], SMeanNO[zbin+1], SMeanNO2[zbin+1], SMeanOCS[zbin+1], SMeanHNO3[zbin+1];
extern int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
extern int    numx, numc, numf, numa, waternum, waterx, numr, numm, numt, nump;
extern double xx[zbin+1][NSP+1];
extern double mkv[], Tnew[], Pnew[];
extern double clouds[zbin+1][NSP+1]; // Cloud abundances for condensible species (number density) molecules cm^-3

// Photochemistry arrays (PhotoCross)
extern double **cross, **crosst;  // PhotoCross cross sections
extern int *stdcross;             // PhotoCross species indices
extern double *qysum;              // PhotoCross quantum yields

// Cloud optical property arrays for radiative transfer [layer][wavelength]
extern double **cH2O, **aH2O, **gH2O;  // H2O cloud: cross-section (cm^-1), albedo, asymmetry g
extern double **cNH3, **aNH3, **gNH3;  // NH3 cloud: cross-section (cm^-1), albedo, asymmetry g
extern double H2H2CIA[zbin+1][NLAMBDA], H2HeCIA[zbin+1][NLAMBDA], H2HCIA[zbin+1][NLAMBDA], N2H2CIA[zbin+1][NLAMBDA], N2N2CIA[zbin+1][NLAMBDA], CO2CO2CIA[zbin+1][NLAMBDA];
extern double MeanH2H2CIA[], MeanH2HeCIA[], MeanH2HCIA[], MeanN2H2CIA[], MeanN2N2CIA[],MeanCO2CO2CIA[];
extern double SMeanH2H2CIA[], SMeanH2HeCIA[], SMeanH2HCIA[], SMeanN2H2CIA[], SMeanN2N2CIA[], SMeanCO2CO2CIA[];

// Dynamic condensibles management variables
extern int NCONDENSIBLES;
extern int CONDENSIBLES[];
// ALPHA_RAINOUT is now a single constant defined in AlphaAb.h

// Additional global variables
extern double Tdoub[];
extern double TAUdoub[];  // Double grid optical depth
extern double Pdoub[];    // Double grid pressure
extern double MMdoub[];   // Double grid number density
extern double zdoub[];    // Double grid altitude
extern double rt_drfluxmax_init; // Initial radiative flux maximum

extern int RTstepcount;
extern double GA; // Gravitational acceleration

// Enhanced cloud physics arrays
extern double particle_r2[zbin+1][MAX_CONDENSIBLES];
extern double particle_r1[zbin+1][MAX_CONDENSIBLES];  // Surface-area-weighted radius [μm] (for cloud optics)
extern double particle_r0[zbin+1][MAX_CONDENSIBLES];  // Mode radius (nucleation/monomer radius) [μm]
extern double particle_VP[zbin+1][MAX_CONDENSIBLES];  // Particle volume [cm³]
extern double particle_mass[zbin+1][MAX_CONDENSIBLES];  // Particle mass [kg]
extern double particle_number_density[zbin+1][MAX_CONDENSIBLES]; // Particle number density [particles/m³]
extern double fall_velocity_ms[zbin+1][MAX_CONDENSIBLES];
extern double cloud_retention[zbin+1][MAX_CONDENSIBLES];

// NOTE: Function declarations have been moved to module-specific headers:
// - Convection/condensation/cloud physics: conv_cond_funcs.h
// - Opacity reading/reinterpolation: readcross.h
// - CIA opacity: readcia.h
// - Cloud optical properties: cloud_optics.h

#endif /* !__GLOBAL_H__ */

/*---- end ------------------------ global.h ---------------------- */
