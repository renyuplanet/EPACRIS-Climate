/* Solve Equations for Block Tridiagonal Matrix */

#include <math.h>

void BTridiagonal(double **m, double x[], int NB, int ND)
/* m: the block tridiagonal matrix to solve
   x: the right hand side of equation, also the output result
   NB: number of blocks in diagonal
   ND: dimension of each blocks */
{
	double **a, **b, **c, *cc, ***r, **d, *dd, ddd;
	int *indx, i, j, k, z, temp, singularflag;
	indx=ivector(1,ND);
	cc=dvector(1,ND);
	d=dmatrix(1,ND,1,NB);
//printf("%s\n", "Start BTridiagonal.c"); //template
	for (i=1; i<=ND; i++) {
		for (j=1; j<=NB; j++) {
			d[i][j]=x[(j-1)*ND+i];
		}
	}
	dd=dvector(1,ND);
	r=f3tensor(1,NB-1,1,ND,1,ND);
	
	a=submatrix(m,1,ND,1,ND,1,1);
	c=submatrix(m,1,ND,ND+1,2*ND,1,1);
	ludcmp(a,ND,indx,&ddd);
	for (j=1; j<=ND; j++) {
		for (i=1; i<=ND; i++) {cc[i]=c[i][j];}
		lubksb(a,ND,indx,cc);
		for (i=1; i<=ND; i++) {r[1][i][j]=cc[i];}
	} /* Matrix inversion */
	for (i=1; i<=ND; i++) {dd[i]=d[i][1];}
	lubksb(a,ND,indx,dd);
	for (i=1; i<=ND; i++) {d[i][1]=dd[i];}
	
//printf("%s %d\n", "BTridiagonal.c, NB = ",NB); //template
    for (z=2; z<NB; z++) {
//printf("%s %d\n", "BTridiagonal.c, z = ",z); //template
		temp=(z-1)*ND;
		a=submatrix(m,temp+1,temp+ND,temp+1,temp+ND,1,1);
		c=submatrix(m,temp+1,temp+ND,temp+1+ND,temp+2*ND,1,1);
		b=submatrix(m,temp+1,temp+ND,temp+1-ND,temp,1,1);
		for (k=1; k<=ND; k++) {
			for (i=1; i<=ND; i++) {
				for (j=1; j<=ND; j++) {
					a[i][j] -= b[i][k]*r[z-1][k][j];
				}
			}
		} /* a=a-b*r */
		ludcmp(a,ND,indx,&ddd);
		for (j=1; j<=ND; j++) {
			for (i=1; i<=ND; i++) {cc[i]=c[i][j];}
			lubksb(a,ND,indx,cc);
			for (i=1; i<=ND; i++) {r[z][i][j]=cc[i];}
		} /* Matrix inversion */
		for (i=1; i<=ND; i++) {dd[i]=d[i][z];}
		for (j=1; j<=ND; j++) {
			for (i=1; i<=ND; i++) {
				dd[i] -= b[i][j]*d[j][z-1];
			}
		}
		lubksb(a,ND,indx,dd);
		for (i=1; i<=ND; i++) {d[i][z]=dd[i];}
	}
	
	z=NB;
	temp=(z-1)*ND;
	a=submatrix(m,temp+1,temp+ND,temp+1,temp+ND,1,1);
	b=submatrix(m,temp+1,temp+ND,temp+1-ND,temp,1,1);
	for (k=1; k<=ND; k++) {
		for (i=1; i<=ND; i++) {
			for (j=1; j<=ND; j++) {
				a[i][j] -= b[i][k]*r[z-1][k][j];
			}
		}
	}
	ludcmp(a,ND,indx,&ddd);
	for (i=1; i<=ND; i++) {dd[i]=d[i][z];}
	for (j=1; j<=ND; j++) {
		for (i=1; i<=ND; i++) {
			dd[i] -= b[i][j]*d[j][z-1];
		}
	}
	lubksb(a,ND,indx,dd);
	for (i=1; i<=ND; i++) {d[i][z]=dd[i];}
	
	for (i=1; i<=ND; i++) {
		x[temp+i]=d[i][z];
	}

	for (z=NB-1; z>=1; z--) {
		temp=(z-1)*ND;
		for (i=1; i<=ND; i++) {x[temp+i]=d[i][z];}
		for (j=1; j<=ND; j++) {
			for (i=1; i<=ND; i++) {
				x[temp+i] -= r[z][i][j]*x[temp+j+ND];
			}
		}
	}
//printf("%s\n", "BTridiagonal.c l:110"); //template
	
	free_f3tensor(r,1,NB-1,1,ND,1,ND);
	free_ivector(indx,1,ND);
	free_dvector(cc,1,ND);
	free_dmatrix(d,1,ND,1,ND);
	free_dvector(dd,1,ND);
	free_submatrix(a,1,ND,1,ND);
	free_submatrix(b,1,ND,1,ND);
	free_submatrix(c,1,ND,1,ND);
}
