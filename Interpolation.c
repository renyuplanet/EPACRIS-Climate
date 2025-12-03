/*Function to put the cross section data  into desired wavelength by interpolation*/

#include <math.h>
#include "locate.c"

void Interpolation(double wavelength[], int size, double Cross[], double xx[], double yy[], int s, int IFEX);

void Interpolation(double wavelength[], int size, double Cross[], double xx[], double yy[], int s, int IFEX)
{

      /*Convert into desired wavelength*/
      int i=0;
      unsigned long j=0;
      while (i<size)
      {
            if (wavelength[i]<xx[0]) {
				if (IFEX == 1) {
					Cross[i]=fmax(0, yy[0]+(wavelength[i]-xx[0])*(yy[1]-yy[0])/(xx[1]-xx[0]) ); /* linear extrapolation */
				} else if (IFEX==2) {
                    Cross[i]=yy[0];
                } else {
					Cross[i]=0;
				}
			}
            if (wavelength[i]>xx[s-1]) {
				if (IFEX == 1) {
					Cross[i]=fmax(0, yy[s-1]+(wavelength[i]-xx[s-1])*(yy[s-1]-yy[s-2])/(xx[s-1]-xx[s-2]) ); /* linear extrapolation */
				} else if (IFEX==2) {
                    Cross[i]=yy[s-1];
                } else {
					Cross[i]=0;
				}
			}
            if (wavelength[i]>=xx[0] && wavelength[i]<=xx[s-1])
            {
            j=locate(xx, s, wavelength[i]);
				if (xx[j+1] == xx[j]) {
					Cross[i]=yy[j+1];
				} else {
					Cross[i]=yy[j]+(wavelength[i]-xx[j])*(yy[j+1]-yy[j])/(xx[j+1]-xx[j]);
				}
            }
            i=i+1;
      }
}

void Interpolationr1(double wavelength[], int size, double Cross[], double xx[], double yy[], int s, int IFEX);

void Interpolationr1(double wavelength[], int size, double Cross[], double xx[], double yy[], int s, int IFEX)
{
    
    /*Convert into desired wavelength*/
    int i=0;
    unsigned long j=0;
    while (i<size)
    {
        if (wavelength[i]>xx[1]) {
            if (IFEX == 1) {
                Cross[i]=fmax(0, yy[1]+(wavelength[i]-xx[1])*(yy[2]-yy[1])/(xx[2]-xx[1]) ); /* extrapolation */
            }else {
                Cross[i]=0;
            }
        }
        if (wavelength[i]<xx[s-1]) {
            if (IFEX == 1) {
                Cross[i]=fmax(0, yy[s-1]+(wavelength[i]-xx[s-1])*(yy[s-1]-yy[s-2])/(xx[s-1]-xx[s-2]) ); /* extrapolation */
            }else {
                Cross[i]=0;
            }
        }
        if (wavelength[i]<=xx[1] && wavelength[i]>=xx[s-1])
        {
            j=locate(xx, s, wavelength[i]);
            if (xx[j+1] == xx[j]) {
                Cross[i]=yy[j+1];
            } else {
                Cross[i]=yy[j]+(wavelength[i]-xx[j])*(yy[j+1]-yy[j])/(xx[j+1]-xx[j]);
            }
        }
        i=i+1;
    }
}


void InterpolationLog(double wavelength[], int size, double Cross[], double xx[], double yy[], int s, int IFEX);

void InterpolationLog(double wavelength[], int size, double Cross[], double xx[], double yy[], int s, int IFEX)
{
	
	/* Linear Interpolation in log scale */
	int i=0;
	unsigned long j=0;
	
	double wavelengthlog[size];
	double xxlog[s];	
	for (i=0; i<size; i++) { wavelengthlog[i]=log(wavelength[i]); }
	for (i=0; i<s; i++) { xxlog[i]=log(xx[i]); }
	
	i=0;
	while (i<size)
	{
		if (wavelength[i]<xx[0]) {
			if (IFEX==1) {
				Cross[i]=fmax(0.0, yy[0]+(wavelengthlog[i]-xxlog[0])*(yy[1]-yy[0])/(xxlog[1]-xxlog[0]) );
			}else {
				Cross[i]=0;
			}
		}
		
		if (wavelength[i]>xx[s-1]) {
			if (IFEX==1) {
				Cross[i]=fmax(0.0, yy[s-1]+(wavelengthlog[i]-xxlog[s-1])*(yy[s-1]-yy[s-2])/(xxlog[s-1]-xxlog[s-2]) );
			}else {
				Cross[i]=0;
			}
		}
		
		if (wavelength[i]>=xx[0] && wavelength[i]<=xx[s-1])
		{
            j=locate(xx, s, wavelength[i]);
            Cross[i]=yy[j]+(wavelengthlog[i]-xxlog[j])*(yy[j+1]-yy[j])/(xxlog[j+1]-xxlog[j]);
		}
		i=i+1;
	}
}

double Interpolation2D(double x, double y, double xx[], int nx, double yy[], int ny, double **data);
/* x interpolated in linear, and y interpolated in log */
/* require xx and yy from small to large */
/* Index of xx and yy from 0 */
/* No extrapolation so far */

double Interpolation2D(double x, double y, double xx[], int nx, double yy[], int ny, double **data)
{
	unsigned long i, j;
	int flag=1;
	double logyy[ny];
	double logy;
	logy = log(y);
	for (i=0; i<ny; i++) { logyy[i] = log(yy[i]); }
	double f1, f2, f;

	i=0;
	j=0;
	if (x>=xx[0] && x<=xx[nx-1]) {
		i=locate(xx,nx,x);
	} else {
		flag=0;
		f=0.0;
	}
	if (y>=yy[0] && y<=yy[ny-1]) {
		j=locate(yy,ny,y);
	} else {
		flag=0;
		f=0.0;
	}
	
	if (flag==1) {
		f1 = data[i][j] + (x-xx[i])*(data[i+1][j]-data[i][j])/(xx[i+1]-xx[i]);
		f2 = data[i][j+1] + (x-xx[i])*(data[i+1][j+1]-data[i][j+1])/(xx[i+1]-xx[i]);
		f = f1 + (logy-logyy[j])*(f2-f1)/(logyy[j+1]-logyy[j]);
	} else {
		f = 0.0;
	}
	
	return f;
	
}