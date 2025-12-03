/*Function to read in Gibbs Free Energy of Formation data */
/* Save in J mole-1 */
/* Interpolation to specified T*/

#include <math.h>
#include "constant.h"

void GetGibbsForm();

void GetGibbsForm()
{
	int i,j;
	FILE *fim;
	double temp[31];
	double gibbs[31];
	double gs[zbin+1];
	
	/* open file */
	fim = fopen("Library/Gibbs/GibbsForm.dat","r");
	for (i=1; i<=30; i++) {
		fscanf(fim, "%lf", temp+i);
	}
	for (j=1; j<=NSP; j++) {
		for (i=1; i<=30; i++) {
			gibbs[i] = 0.0;
			fscanf(fim, "%lf", gibbs+i);
		}
		Interpolation(tl, zbin+1, gs, temp, gibbs, 31, 0);
		for (i=1; i<=zbin; i++) {
			GibbsForm[j][i] = gs[i];
			gs[i] = 0.0;
		}
	}
	fclose(fim);
	
}