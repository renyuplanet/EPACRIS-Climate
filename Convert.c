/*Function to convert different representations of variables*/
#include <math.h>
#include "constant.h"


/*Inidivudal species concentrations are converted to xx array*/
void Convert1(double Con[], double ConC[], double Conf[], 
			  int labelx[], int labelc[], int labelf[])
{
	int i, j;
	for (i=1; i<=zbin; i++)
	{
		for (j=1; j<=NSP; j++) {xx[i][j]=0;} /* clean-up xx */
		for (j=1; j<=numx; j++) {xx[i][labelx[j]]=Con[(i-1)*numx+j];}
		for (j=1; j<=numc; j++) {xx[i][labelc[j]]=ConC[(i-1)*numc+j];}
		for (j=1; j<=numf; j++) {xx[i][labelf[j]]=Conf[(i-1)*numf+j];}
	}
}


// For debugging purposes, prints the original values of the species concentrations
// void Convert1(double Con[], double ConC[], double Conf[], 
// 			  int labelx[], int labelc[], int labelf[])
// {
// 	int i, j;
// 	for (i=1; i<=zbin; i++)
// 	{
// 		for (j=1; j<=NSP; j++) {xx[i][j]=0;} /* clean-up xx */
// 		for (j=1; j<=numx; j++) {
// 			if (labelx[j] == 7 || labelx[j] == 52 || labelx[j] == 20) {
// 				printf("Original Con value for X species %d (label %d): %e\n", 
// 					   j, labelx[j], Con[(i-1)*numx+j]);
// 			}
// 			xx[i][labelx[j]]=Con[(i-1)*numx+j];
// 			if (labelx[j] == 7 || labelx[j] == 52 || labelx[j] == 20) {
// 				printf("X species %d (label %d) at layer %d: %e\n", j, labelx[j], i, Con[(i-1)*numx+j]);
// 			}
// 		}
// 		for (j=1; j<=numc; j++) {
// 			xx[i][labelc[j]]=ConC[(i-1)*numc+j];
// 			if (labelc[j] == 7 || labelc[j] == 52 || labelc[j] == 20) {
// 				printf("C species %d (label %d) at layer %d: %e\n", j, labelc[j], i, ConC[(i-1)*numc+j]);
// 			}
// 		}
// 		for (j=1; j<=numf; j++) {
// 			xx[i][labelf[j]]=Conf[(i-1)*numf+j];
// 			if (labelf[j] == 7 || labelf[j] == 52 || labelf[j] == 20) {
// 				printf("F species %d (label %d) at layer %d: %e\n", j, labelf[j], i, Conf[(i-1)*numf+j]);
// 			}
// 		}
// 	}
// }

/*converts xx array to individual Con, ConC, Conf arrays*/
void Convert2(double Con[], double ConC[], double Conf[],
			  int labelx[], int labelc[], int labelf[])
{
	int i, j;
	for (i=1; i<=zbin; i++)
	{
		for (j=1; j<=numx; j++) {Con[(i-1)*numx+j]=xx[i][labelx[j]];}
		for (j=1; j<=numc; j++) {ConC[(i-1)*numc+j]=xx[i][labelc[j]];}
		for (j=1; j<=numf; j++) {Conf[(i-1)*numf+j]=xx[i][labelf[j]];}
	}
}

void Reverse(double P[], int NL);


void Reverse(double P[], int NL)
{
	int n, i;
	double temp;
	n=NL/2;
	
	for (i=0; i<n; i++) {
		temp=P[i];
		P[i]=P[NL-i-1];
		P[NL-i-1]=temp;
	}
}
