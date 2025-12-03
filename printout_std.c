void printout_std(double z[], char outstd[]);

void printout_std(double z[], char outstd[])
{
	int i, j;
	FILE *fp;
	fp=fopen(outstd,"w");
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
			fprintf(fp, "%e\t", xx[j][i]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
