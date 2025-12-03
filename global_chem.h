/*----------------------- global.h --------------------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: July 20, 2011

--------------------------------------------------------------------- */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

/*---- External Variables ------------------------------------------- */

extern double meanmolecular[];
extern double thickl;
extern double zl[];
extern double pl[];
extern double tl[];
extern double MM[], MMZ[];
extern double wavelength[];
extern double longwavelength[];
extern double solar[];
extern double crossr[], crossa[3][NLAMBDA], sinab[3][NLAMBDA], asym[3][NLAMBDA];
extern double **opacCO2, **opacO2, **opacSO2, **opacH2O, **opacOH, **opacH2CO; 
extern double **opacH2O2, **opacHO2, **opacH2S, **opacCO, **opacO3, **opacCH4; 
extern double **opacNH3;
extern double **opacC2H2, **opacC2H4, **opacC2H6, **opacCH2O2;
extern double **opacHCN, **opacN2O, **opacNO, **opacNO2, **opacHNO3, **opacOCS;
extern double MeanCO2[], MeanO2[], MeanSO2[], MeanH2O[], MeanOH[], MeanH2CO[];
extern double MeanH2O2[], MeanHO2[], MeanH2S[], MeanCO[], MeanO3[], MeanCH4[];
extern double MeanNH3[];
extern double MeanC2H2[], MeanC2H4[], MeanC2H6[], MeanCH2O2[];
extern double MeanHCN[zbin+1], MeanN2O[zbin+1], MeanNO[zbin+1], MeanNO2[zbin+1], MeanOCS[zbin+1], MeanHNO3[zbin+1];
extern double SMeanCO2[], SMeanO2[], SMeanSO2[], SMeanH2O[], SMeanOH[], SMeanH2CO[];
extern double SMeanH2O2[], SMeanHO2[], SMeanH2S[], SMeanCO[], SMeanO3[], SMeanCH4[];
extern double SMeanNH3[];
extern double SMeanC2H2[], SMeanC2H4[], SMeanC2H6[], SMeanCH2O2[];
extern double SMeanHCN[zbin+1], SMeanN2O[zbin+1], SMeanNO[zbin+1], SMeanNO2[zbin+1], SMeanOCS[zbin+1], SMeanHNO3[zbin+1];
extern double H2H2CIA[zbin+1][NLAMBDA], H2HeCIA[zbin+1][NLAMBDA], H2HCIA[zbin+1][NLAMBDA], N2H2CIA[zbin+1][NLAMBDA], N2N2CIA[zbin+1][NLAMBDA], CO2CO2CIA[zbin+1][NLAMBDA];
extern double MeanH2H2CIA[], MeanH2HeCIA[], MeanH2HCIA[], MeanN2H2CIA[], MeanN2N2CIA[], MeanCO2CO2CIA[];
extern double SMeanH2H2CIA[], SMeanH2HeCIA[], SMeanH2HCIA[], SMeanN2H2CIA[], SMeanN2N2CIA[], SMeanCO2CO2CIA[];
extern double rainoutrate[zbin+1][NSP+1];
extern double Vesc[], VFall[];
extern double nsH2O[], nsH2SO4[], nsS8[], tcondfH2O[], tcondfH2SO4[], tcondfS8[];
extern double kk[zbin+1][NKin+1], kkM[zbin+1][NKinM+1], kkT[zbin+1][NKinT+1];
extern double Rkk[zbin+1][NKin+1], RkkM[zbin+1][NKinM+1], RkkT[zbin+1][NKinT+1];
extern int    ReactionR[NKin+1][7], ReactionM[NKinM+1][5], ReactionP[NPho+1][9], ReactionT[NKinT+1][4];
extern int    numr, numm, numt, nump, numx, numc, numf, numa, waternum, waterx;
extern double **DM, **dl, KE[];
extern double xx[zbin+1][NSP+1];
extern double mkv[], Tnew[], Pnew[];
extern double GibbsForm[NSP+1][zbin+1];

#endif /* !__GLOBAL_H__ */

/*---- end ------------------------ global.h ---------------------- */
