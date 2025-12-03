void printoutrate(int zone_r[], int zone_m[], int zone_t[], double **JJ, int zone_p[], char outchemicalrate[]);

void printoutrate(int zone_r[], int zone_m[], int zone_t[], double **JJ, int zone_p[], char outchemicalrate[])
{
     int i,j,std;
	 double rate;
	 FILE *fout;
     
     fout=fopen(outchemicalrate,"w");
	 fprintf(fout, "%s\n", "Type	Std	Rate");
     
	for (i=1; i<=numr; i++) {
        std=zone_r[i]; /* Get the standard reaction number from list */
		fprintf(fout, "%s%d\t", "R", std);
	    for (j=1; j<=zbin; j++) {
			rate=kk[j][std]*xx[j][ReactionR[std][1]]*xx[j][ReactionR[std][2]];
			if (ReactionR[std][3]==1) { rate -= Rkk[j][std]*xx[j][ReactionR[std][4]]; }
			if (ReactionR[std][3]==2) { rate -= Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][5]]; }
			if (ReactionR[std][3]==3) { rate -= Rkk[j][std]*xx[j][ReactionR[std][4]]*xx[j][ReactionR[std][5]]*xx[j][ReactionR[std][6]]; }
			fprintf(fout, "%e\t", rate);
		}
		fprintf(fout, "\n");
	}
	
	for (i=1; i<=numm; i++) {
        std=zone_m[i]; /* Get the standard reaction number from list */
		fprintf(fout, "%s%d\t", "M", std);
	    for (j=1; j<=zbin; j++) {
			rate=kkM[j][std]*xx[j][ReactionM[std][1]]*xx[j][ReactionM[std][2]];
			if (ReactionM[std][4]==0) { rate -= RkkM[j][std]*xx[j][ReactionM[std][3]]; }
			if (ReactionM[std][4]>0) { rate -= RkkM[j][std]*xx[j][ReactionM[std][3]]*xx[j][ReactionM[std][4]]; }
			fprintf(fout, "%e\t", rate);
		}
		fprintf(fout, "\n");
	}
	
	for (i=1; i<=nump; i++) {
        std=zone_p[i]; /* Get the standard reaction number from list */
		fprintf(fout, "%s%d\t", "P", std);
	    for (j=1; j<zbin; j++) {
			rate=JJ[j][i]*xx[j][ReactionP[std][1]];
			fprintf(fout, "%e\t", rate);
		}
		j=zbin;
		rate=JJ[j][i]*xx[j][ReactionP[std][1]];
		fprintf(fout, "%e\n", rate);
	}
	
	for (i=1; i<=numt; i++) {
        std=zone_t[i]; /* Get the standard reaction number from list */
		fprintf(fout, "%s%d\t", "T", std);
	    for (j=1; j<=zbin; j++) {
			rate=kkT[j][std]*xx[j][ReactionT[std][1]];
			if (ReactionT[std][3]==0) { rate -= RkkT[j][std]*xx[j][ReactionT[std][2]]; }
			if (ReactionT[std][3]>0) { rate -= RkkT[j][std]*xx[j][ReactionT[std][2]]*xx[j][ReactionT[std][3]]; }
			fprintf(fout, "%e\t", rate);
		}
		fprintf(fout, "\n");
	}
	
     fclose(fout);
}
