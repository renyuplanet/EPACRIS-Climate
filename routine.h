/*----------------------- routine.h --------------------------------

Author: Renyu Hu (hury@mit.edu)
Last modified: February 7, 2010

--------------------------------------------------------------------- */

#ifndef _ROUTINE_H_
#define _ROUTINE_H_

/*Chemical Equilibrium Calculations*/
void chemquil(double P[], double T[], int NL, int species[], int NS, double **mixequil, char atomfile[]);

/*TP profile*/
void TPPara(double P[], double T[], int Tinv, int NL, 
			double P0, double T0, double P1, double T1, double P2, double T2, double P3, double T3, double P4);
double TPScale(double P[], double T[], int NL, double z[]);

/*Chemical Kinetic Rates*/
void ReactionRate();
void ReactionRateM();
void ReactionRateT();
void Photorate(double **Radiation, double **cross, double **crosst, double **qy, double **qyt, int zone_p[], int nump, double **JJ);

/*Convertor*/
void Convert1(double Con[], double ConC[], double Conf[], int labelx[], int labelc[], int labelf[]);
void Convert2(double Con[], double ConC[], double Conf[], int labelx[], int labelc[], int labelf[]);

/*Get Data*/
int LineNumber(FILE *iop, int n);
/* Get absoprtion CrossSection data from a txt file and put data into two arrays*/      
void GetData3(FILE *iop, int n, int size, double *var1, double *var2, double *var3); 
void GetData4(FILE *iop, int n, int size, double *var1, double *var2, double *var3, double *var4);
void GetData(FILE *iop, int n, int size, double *wavelength, double *CrossSection); 
void GetReaction();

/*Integration*/
double Trapz(double x[], double y[], int size);

/*Equations*/
void ChemEqu(double x[], double xc[], double xf[], double fvec[], double **Jvalue, int lx[], int lc[], int lf[], int laer[],
			 int zone_r[], int zone_m[], int zone_t[], double **JJ, int zone_p[], double Upflux[], double Loflux[], double Depo[], int listFix[]);
			   
void BTridiagonal(double **m, double x[], int NB, int ND);
int BTridiagonalm(double **m, double x[], int NB, int ND);

/* Define structure */
struct Molecule {
	char name[20];
	char type[10];
	int num;
	int mass;
	double mix;
	double upper;
	int lowertype;
	double lower;
	double lower1;
};

struct Reaction {
       int dum;
       char type[10];
       int num;
};

#endif
