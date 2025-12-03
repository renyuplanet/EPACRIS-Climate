/*Trapz integration of two vectors*/
/*To get the size , use 
int size=sizeof(x)/sizeof(x[0]); */

double Trapz(double x[], double y[], int size)
{
      double sum=0;
      int i=0; 
      while (i+1<size)
      {
            sum=sum+(y[i]+y[i+1])*(x[i+1]-x[i])/2.0;
            i++;
      }
      return sum;
} 
