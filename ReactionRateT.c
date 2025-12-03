#include <math.h>

/* s-1 */
/* M in molecule cm-3 */

void ReactionRateT()
{
     int i;
	 double RH, LH, ind, FC;
	
     for(i=1;i<=zbin;i++)
     {
       kkT[i][1]=7.16E-10*exp(-11200.0/tl[i])*MM[i];
       kkT[i][2]=2.41E-8*pow(tl[i]/298.0, -1.18)*exp(-24415.0/tl[i])*MM[i];
       kkT[i][3]=2.01E-7*exp(-22852.0/tl[i])*MM[i];
	
	   LH=5.88E-10*exp(-28265.0/tl[i]);
	   RH=1.30E+11*exp(-30000.0/tl[i]);
	   ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	   kkT[i][4]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
       kkT[i][5]=1.88E-4*pow(tl[i]/298.0, -3.37)*exp(-37645.0/tl[i])*MM[i];
	
	   LH=2.5E-14*exp(-1230.0/tl[i]);
       RH=2.5E+6*exp(-6100.0/tl[i]);
	   ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	   kkT[i][6]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
       kkT[i][7]=1.0E-3*pow(tl[i]/298.0, -3.5)*exp(-11000.0/tl[i])*MM[i];
       kkT[i][8]=5.48E-7*pow(tl[i]/298.0, -1.24)*exp(-25312.0/tl[i])*MM[i];
		 
       LH=1.98E-3*pow(tl[i]/298.0, -3.8)*exp(-25257.0/tl[i]);
	   RH=1.09E+16*pow(tl[i]/298.0, -1.23)*exp(-25016.0/tl[i]);
	   ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	   kkT[i][9]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);	
		 
       LH=1.15E-6*exp(-23092.0/tl[i]);
	   RH=9.33E+15*exp(-24656.0/tl[i]);	 
	   ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	   kkT[i][10]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);	
		 
	   LH=8.0E-2*pow(tl[i]/298.0, -6.55)*exp(-26099.0/tl[i]);
	   RH=5.56E+17*pow(tl[i]/298.0, -2.27)*exp(-26340.0/tl[i]);
	   ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	   kkT[i][11]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
       kkT[i][12]=6.31E17*exp(-13110.0/tl[i]);
       kkT[i][13]=4.10E-5*exp(-10600.0/tl[i])*MM[i];
		 
		 LH=7.7E-9*exp(-33075.0/tl[i]);
		 RH=3.7E+13*exp(-36202.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][14]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
		 LH=2.81E-9*exp(-25738.0/tl[i]);
		 RH=4.46E+13*exp(-34398.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][15]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
		 LH=6.73E-9*exp(-26700.0/tl[i]);
		 RH=7.49E+14*exp(-34518.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][16]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 

       kkT[i][17]=5.0E4;
       kkT[i][18]=1.0E14*exp(-7500.0/tl[i]);
       kkT[i][19]=2.16E-8*exp(-33556.0/tl[i])*MM[i];
       kkT[i][20]=1.1E-7*exp(-33075.0/tl[i])*MM[i];
       kkT[i][21]=2.12E13*pow(tl[i]/298.0, 1.22)*exp(-43539.0/tl[i]);
		 
		 LH=1.484E+30*pow(tl[i],-10.2)*exp(-52454.0/tl[i])+1.223E+17*pow(tl[i],-6.577)*exp(-48007.0/tl[i]);
		 RH=3.121E+18*pow(tl[i],-1.017)*exp(-46156.0/tl[i]);
		 FC=0.9922*exp(-tl[i]/943.0)+0.0078*exp(-tl[i]/47310.0)+exp(-47110.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
		 kkT[i][22]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(FC,ind);
		 
       kkT[i][23]=6.0E14*exp(-21288.0/tl[i]);
       kkT[i][24]=4.15E13*exp(-22731.0/tl[i]);
       kkT[i][25]=0.19*pow(tl[i]/298.0, -7.5)*exp(-22852.0/tl[i])*MM[i];
       kkT[i][26]=5.8E-8*exp(-35961.0/tl[i])*MM[i];
       kkT[i][27]=1.69E-6*exp(-16838.0/tl[i])*MM[i];
       kkT[i][28]=8.11E17*pow(tl[i]/298.0, -1.23)*exp(-51356.0/tl[i]);
       kkT[i][29]=1.08E-8*exp(-29828.0/tl[i])*MM[i];
		 
		 LH=5.98E-9*exp(-29828.0/tl[i]);
		 RH=3.00E+14*exp(-35720.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][30]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);

       kkT[i][31]=2.2E12*exp(-15154.0/tl[i]);
       kkT[i][32]=6.82E-3*pow(tl[i]/298.0, -8.62)*exp(-11300.0/tl[i])*MM[i];
       kkT[i][33]=2.0E15*exp(-39810.0/tl[i]);
       kkT[i][34]=3.0E14*exp(-42216.0/tl[i]);
       kkT[i][35]=1.0E15*exp(-42817.0/tl[i]);
       kkT[i][36]=1.33E15*pow(tl[i]/298.0, -2.02)*exp(-12749.0/tl[i]);
       kkT[i][37]=6.19E11*exp(-11900.0/tl[i]);
       kkT[i][38]=1.07E14*pow(tl[i]/298.0, -0.69)*exp(-11187.0/tl[i]);
		 
		 LH=7.69E-4*pow(tl[i]/298.0, -3.10)*exp(-51356.0/tl[i]);
		 RH=6.00E+13*exp(-50154.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][39]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
       kkT[i][40]=1.69E-9*exp(-23453.0/tl[i])*MM[i];
       kkT[i][41]=9.47E-7*exp(-40051.0/tl[i])*MM[i];
       kkT[i][42]=1.4E-8*exp(-29467.0/tl[i])*MM[i];
		 
		 LH=4.16E-7*pow(tl[i]/298.0, -3.29)*exp(-9610.0/tl[i]);
		 RH=3.42E13*pow(tl[i]/298.0, 0.90 )*exp(-9240.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][43]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
		 LH=4.61E-4*pow(tl[i]/298.0, -5.87)*exp(-15635.0/tl[i]);
		 RH=4.53E16*pow(tl[i]/298.0, -1.07)*exp(-14312.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][44]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
		 LH=2.73E-8*pow(tl[i]/298.0, -2.82)*exp(3750.0/tl[i]);
		 RH=9.93E15*pow(tl[i]/298.0, -1.07)*exp(-3900.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][45]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
       kkT[i][46]=6.0E-11*exp(-7721.0/tl[i])*MM[i];
		 
		 LH=9.0E-11*pow(tl[i]/298.0, 0.00 )*exp(-6790.0/tl[i]);
		 RH=1.69E14*pow(tl[i]/298.0, -0.39)*exp(-13230.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][47]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);

		 LH=4.92E-29*pow(tl[i]/298.0, -2.40)*exp(-18882.0/tl[i]);
		 RH=29.85   *pow(tl[i]/298.0, 0.13 )*exp(-18401.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][48]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
		 LH=7.23E-31*pow(tl[i]/298.0, -3.15)*exp(-18642.0/tl[i]);
		 RH=125.0   *pow(tl[i]/298.0, 0.41 )*exp(-17800.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][49]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
       kkT[i][50]=1.32E-6*exp(-17199.0/tl[i])*MM[i];
       kkT[i][51]=3.98E13*exp(-19364.0/tl[i]);
       kkT[i][52]=1.0E13*exp(-16838.0/tl[i]);
      
		 LH=3.16E-4*pow(tl[i]/298.0, -4.72)*exp(-13590.0/tl[i]);
		 RH=4.34E14*pow(tl[i]/298.0, 0.0  )*exp(-12989.0/tl[i]);
		 ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
		 kkT[i][53]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
		 
		 kkT[i][54]=1.0E-12*MM[i] + 1.5E+3;
		 kkT[i][55]=1.0E-11*MM[i] + 2.2E+4;
		 kkT[i][56]=1.5E-13*MM[i] + 1.13E+3;
		 kkT[i][57]=1.93E-4*pow(tl[i]/298.0,-2.44)*exp(-62782.1/tl[i])*MM[i];
		 kkT[i][58]=3.16E-10*exp(-32954.6/tl[i])*MM[i];
		 kkT[i][59]=2.92E-8*exp(-33315.4/tl[i])*MM[i];
		 kkT[i][60]=4.07E-10*exp(-30910.0/tl[i])*MM[i];
		 kkT[i][61]=3.98E+12*exp(-18402.0/tl[i]);
		 kkT[i][62]=1.40E+13*exp(-30188.3/tl[i]);
		 kkT[i][63]=1.26E+13*exp(-16838.0/tl[i]);
		 kkT[i][64]=1.18E+18*pow(tl[i]/298.0, -1.2)*exp(-49191.0/tl[i]);
		 kkT[i][65]=2.50E+15*exp(-43659.0/tl[i]);
		 kkT[i][66]=4.00E+13*exp(-40291.0/tl[i]);
		 kkT[i][67]=1.09E+13*pow(tl[i]/298.0, 0.17)*exp(-17921.6/tl[i]);
		 kkT[i][68]=1.31E+13*pow(tl[i]/298.0, 0.87)*exp(-15274.6/tl[i]);
		 kkT[i][69]=1.58E+16*exp(-49071.0/tl[i]);
		 kkT[i][70]=1.07E+12*pow(tl[i]/298.0, -15.74)*exp(-49672.0/tl[i])*MM[i];
		 kkT[i][71]=3.38E+10;
		 kkT[i][72]=1.66E-8*exp(-29707.0/tl[i])*MM[i];
		 kkT[i][73]=6.3E+13*exp(-43779.0/tl[i]);
		 kkT[i][74]=1.15E+15*exp(-43539.0/tl[i]);
		 kkT[i][75]=1.6E+14*exp(-20807.0/tl[i]);
		 kkT[i][76]=3.57E-8*pow(tl[i]/298.0, 0.7)*exp(-21288.0/tl[i])*MM[i];
		 kkT[i][77]=7.7E+14*exp(-44260.0/tl[i]);
		 kkT[i][78]=5.0E+15*exp(-38126.0/tl[i]);
		 kkT[i][79]=2.0E+15*exp(-43178.0/tl[i]);
		 kkT[i][80]=2.0E+15*exp(-37645.0/tl[i]);
		 kkT[i][81]=1.58E+16*exp(-55445.0/tl[i]);
		 kkT[i][82]=1.8E+13*exp(-42817.0/tl[i]);
		 kkT[i][83]=1.0E+14*exp(-37645.0/tl[i]);
		 kkT[i][84]=1.0E+16*exp(-36683.0/tl[i]);
		 kkT[i][85]=1.1E+13*pow(tl[i]/298.0, 0.25)*exp(-17921.6/tl[i]);
		 kkT[i][86]=7.7E+13*pow(tl[i]/298.0, 0.77)*exp(-15394.8/tl[i]);
		 kkT[i][87]=1.58E+16*exp(-49312.0/tl[i]);
		 kkT[i][88]=1.0E+17*exp(-42576.0/tl[i]);
		 kkT[i][89]=7.94E+16*exp(-40411.0/tl[i]);
		 kkT[i][90]=6.31E+14*exp(-31151.0/tl[i]);
		 kkT[i][91]=6.61E-9*exp(-20567.0/tl[i])*MM[i];
		 kkT[i][92]=2.66E-11*exp(-24295.0/tl[i])*MM[i];
		 kkT[i][93]=1.36E-7*exp(-30669.0/tl[i])*MM[i];
     }
} 
