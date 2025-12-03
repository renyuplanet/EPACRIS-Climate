void printout_std_t(double z[], char outstdt[]);

void printout_std_t(double z[], char outstdt[])
{
	int i, j;
	FILE *fp;
	fp=fopen(outstdt,"w");
	fprintf(fp, "%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t", "z", "z0", "z1", "T", "P");
	for (i=1; i<=NSP; i++) { fprintf(fp, "%d\t\t", i); }
	fprintf(fp, "\n");
	fprintf(fp, "%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t", "km", "km", "km", "K", "Pa");
	fprintf(fp, "\n");
	for (j=1; j<=zbin; j++) {
		fprintf(fp, "%f\t", zl[j]);
		fprintf(fp, "%f\t%f\t", z[j-1], z[j]);
		fprintf(fp, "%f\t%e\t", tl[j], pl[j]);
		for (i=1; i<=NSP; i++) {
			fprintf(fp, "%e\t", xx[j][i]/MM[j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}


void checkmixequil(int NM, double **mixequil);

void checkmixequil(int NM, double **mixequil)
{
	int i,j;
	int k1,k2,k;
	double f1,f2;
	
	for (i=1; i<=NM; i++) {
		j=1;
		while (j<=zbin) {
			if (isnan(mixequil[j][i])==0) {
				j=j+1;
			} else {
				if (j==1) {
					mixequil[j][i]=1.0E-30;
				} else {
					k1=j-1;
					k2=j;
					while (isnan(mixequil[k2][i])!=0 && k2<zbin) {
						k2=k2+1;
					}
					if (k2>=zbin) {
						j=k2+1;
						for (k=k1+1; k<=zbin; k++) {
							mixequil[k][i] = 1.0E-30;
						}
					} else {
						f1 = mixequil[k1][i];
						f2 = mixequil[k2][i];
						for (k=k1+1; k<k2; k++) {
							mixequil[k][i] = exp(log(f1) + (log(f2)-log(f1))*(k-k1)/(k2-k1));
						}
						j = k2;
					}
				}

				
			}

		}
	}
	
}