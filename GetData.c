#include <stdio.h>

/* Return the number of data lines, to be used before GetData*/
int LineNumber(FILE *iop, int n)
{
      char dataline[n];
      int size=0;
      /*Count number of lines*/
      while (fgets(dataline, n, iop) != NULL )
      {
            size++;
      }
      return size;
}


/* Get absoprtion CrossSection data from a txt file and put data into two arrays*/      
void GetData(FILE *iop, int n, int size, double *wavelength, double *CrossSection) 
{    
      char dataline[n];
      int i=0;
      while (fgets(dataline, n, iop) != NULL )
      {
            sscanf(dataline, "%lf %lf", wavelength+i, CrossSection+i);
            i=i+1;
      }
}

/* Get absoprtion CrossSection data from a txt file and put data into two arrays*/      
void GetData3(FILE *iop, int n, int size, double *var1, double *var2, double *var3) 
{    
      char dataline[n];
      int i=0;
      while (fgets(dataline, n, iop) != NULL )
      {
            sscanf(dataline, "%lf %lf %lf", var1+i, var2+i, var3+i);
            i=i+1;
      }
}

/* Get absoprtion CrossSection data from a txt file and put data into two arrays*/      
void GetData4(FILE *iop, int n, int size, double *var1, double *var2, double *var3, double *var4) 
{    
      char dataline[n];
      int i=0;
      while (fgets(dataline, n, iop) != NULL )
      {
            sscanf(dataline, "%lf %lf %lf %lf", var1+i, var2+i, var3+i, var4+i);
            i=i+1;
      }
}

/* Get the reaction lists */
void GetReaction()
{
	int i,x1,x2,x3,x4,x5,x6;
	FILE *fp;
	char dataline[1000];
	
	fp=fopen("Library/ReactionList/Reaction_R.txt", "r");
	i=1;
	while (fgets(dataline, 1000, fp) != NULL )
	{
		sscanf(dataline, "%d %d %d %d %d %d", &x1, &x2, &x3, &x4, &x5, &x6);
		ReactionR[i][1]=x1;
		ReactionR[i][2]=x2;
		ReactionR[i][3]=x3;
		ReactionR[i][4]=x4;
		ReactionR[i][5]=x5;
		ReactionR[i][6]=x6;
		i=i+1;
	}
	fclose(fp);
	
	fp=fopen("Library/ReactionList/Reaction_M.txt", "r");
	i=1;
	while (fgets(dataline, 1000, fp) != NULL )
	{
		sscanf(dataline, "%d %d %d %d", &ReactionM[i][1], &ReactionM[i][2], &ReactionM[i][3], &ReactionM[i][4]);
		i=i+1;
	}
	fclose(fp);
	
	fp=fopen("Library/ReactionList/Reaction_P.txt", "r");
	i=1;
	while (fgets(dataline, 1000, fp) != NULL )
	{
		sscanf(dataline, "%d %d %d %d %d %d %d %d", &ReactionP[i][1], &ReactionP[i][2], &ReactionP[i][3], &ReactionP[i][4], &ReactionP[i][5], &ReactionP[i][6], &ReactionP[i][7], &ReactionP[i][8]);
		i=i+1;
	}
	fclose(fp);
	
	fp=fopen("Library/ReactionList/Reaction_T.txt", "r");
	i=1;
	while (fgets(dataline, 1000, fp) != NULL )
	{
		sscanf(dataline, "%d %d %d", &ReactionT[i][1], &ReactionT[i][2], &ReactionT[i][3]);
		i=i+1;
	}
	fclose(fp);
}

/* Get Concentration data from a txt file
void GetDataN(FILE *iop, int n, double Con[]);
void GetDataN(FILE *iop, int n, double Con[]) 
{    
	char dataline[n];
	int i=0;
	double a[2];
	while (fgets(dataline, n, iop) != NULL & i<=zbin)
	{
		sscanf(dataline, "%lf %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", 
			   a, &Con[i*NSolve+1], &Con[i*NSolve+2], &Con[i*NSolve+3], 
			   &Con[i*NSolve+4], &Con[i*NSolve+5], &Con[i*NSolve+6], &Con[i*NSolve+7], &Con[i*NSolve+8],
			   &Con[i*NSolve+9], &Con[i*NSolve+10], &Con[i*NSolve+11], &Con[i*NSolve+12], &Con[i*NSolve+13],
			   &Con[i*NSolve+14], &Con[i*NSolve+15], &Con[i*NSolve+16], &Con[i*NSolve+17], &Con[i*NSolve+18],
			   &Con[i*NSolve+19], &Con[i*NSolve+20], a+1);
		i=i+1;
	}
} */
