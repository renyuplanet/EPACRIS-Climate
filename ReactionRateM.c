#include <math.h>

/* cm3 molecule-1 s-1 */
/* M in molecule cm-3 */

void ReactionRateM()
{
  int i, j;
  double RH, LH, ind, K0, K2, K3, FC;
	
  for (i=1; i<=zbin; i++) {
    kkM[i][1]=5.46E-31*pow(tl[i]/298.0,-1.6)*MM[i];
	
	LH=7.00E-32;
	RH=2.06E-11*exp(-57.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][2]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);  
	  
    LH=1.68E-24*pow(tl[i]/298.0,-7.0)*exp(-1390.0/tl[i]);
	RH=1.12E-9*pow(tl[i],-0.5)*exp(-25.0/tl[i]);
	FC=0.62*exp(-tl[i]/1180.0)+0.38*exp(-tl[i]/73.0);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	kkM[i][3]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(FC,ind);  
    
    RH=1.2E-12*pow(tl[i]/300.0, 1.1);
    LH=4.0E-31*pow(tl[i]/300.0, -3.6);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][4]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    LH=4.1E-30*pow(tl[i]/298.0,-2.1);
	RH=3.0E-11*pow(tl[i]/298.0,-0.9);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][5]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);  
	  
	LH=3.0E-34*exp(-1910.0/tl[i]);
	RH=8.4E-13*exp(-3460.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][6]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);  
	  
	if (tl[i]<400.0) {
		LH=3.3E-30*exp(-740.0/tl[i]);
		RH=1.4E-11*exp(-1300.0/tl[i]);
	} else {
		LH=3.77E-8*pow(tl[i],-7.27)*exp(-3633.0/tl[i]);
		RH=7.5E-15*pow(tl[i],-1.50)*exp(-1300.0/tl[i]);
	}
	FC=0.44;
	kkM[i][7]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(FC,ind); 

	LH=1.38E-30;
	RH=2.66E-11*exp(-755.3/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][8]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind); 
	  
	  
    LH=7.69E-30*exp(-380.0/tl[i])*MM[i];
	RH=9.68E-12*pow(tl[i]/298.0,1.28)*exp(-650.0/tl[i]);
	FC=0.24*exp(-tl[i]/4.0)+0.76*exp(-tl[i]/1025.0);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	kkM[i][9]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(FC,ind);  
	  
    kkM[i][10]=3.21E-30*pow(tl[i]/298.0,-2.57)*exp(-215.0/tl[i])*MM[i];
	  
	LH=9.35E-30*pow(tl[i]/298.0,-2.0)*exp(-521.0/tl[i]);
	RH=1.73E-10*pow(tl[i]/298.0,-0.5);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][11]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    LH=5.29E-34*exp(-370.0/tl[i])*MM[i];
	RH=1.96E-13*exp(-1370.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][12]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    kkM[i][13]=6.04E-33*pow(tl[i]/298.0, -1.0)*MM[i];
	
	LH=5.33E-25*pow(tl[i],-2.0);
	RH=2.66E-11;
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][14]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=2.44E-10*pow(tl[i]/298.0, -0.41);
    LH=1.34E-31*pow(tl[i]/298.0, -1.32)*exp(-370.5/tl[i]);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2));
    kkM[i][15]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    LH=4.36E-32*pow(tl[i]/298.0,-1.0);
	RH=1.0E-11;
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2));
	kkM[i][16]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	      
    RH=7.5E-11*pow(tl[i]/300.0, 0.2);
    LH=4.4E-32*pow(tl[i]/300.0, -1.3);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][17]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    LH=6.87E-31*pow(tl[i]/298.0, -2.0);
	RH=1.58E-10*pow(tl[i]/298.0,0.23)*exp(57.7/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][18]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    LH=1.70E-33*exp(1000.0/tl[i]);
	RH=2.31E-13*exp(600.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][19]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    RH=2.9E-12*pow(tl[i]/300.0, -1.1);
    LH=2.0E-31*pow(tl[i]/300.0, -3.4);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][20]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    LH=9.18E-34*pow(tl[i]/298.0, -1.69);
	RH=2.01E-10*pow(tl[i]/298.0, 0.31);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][21]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    kkM[i][22]=9.4E-33*MM[i];
	  
    kkM[i][23]=5.0E-32*MM[i];
	  
	LH=8.3E-38;
	RH=1.94E-20;
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][24]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    kkM[i][25]=1.38E-33*exp(502.7/tl[i])*MM[i];
	  
    kkM[i][26]=5.46E-33*exp(155.0/tl[i])*MM[i];
    
    RH=3.8E-11*pow(tl[i]/300.0, -0.6);
    LH=2.3E-29*pow(tl[i]/300.0, -2.8);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][27]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.9E-11*pow(tl[i]/300.0, -1.8);
    LH=5.3E-29*pow(tl[i]/300.0, -4.4);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][28]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.4E-12*pow(tl[i]/300.0, -0.7);
    LH=2.0E-30*pow(tl[i]/300.0, -4.4);
    ind=1/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][29]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    kkM[i][30]=2.0E-34*MM[i];
	  
    LH=1.7E-33*exp(-1510.0/tl[i]);
	RH=2.66E-14*exp(-1459.0/tl[i]);
	ind=1/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][31]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind)/THREEBODY;
    
    RH=3.0E-11;
    LH=9.0E-32*pow(tl[i]/300, -1.5);
    ind=1.0/(1+pow(log10(LH*MM[i]/RH),2));
    kkM[i][32]=LH*MM[i]/(1+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=2.2E-11*pow(tl[i]/300.0, -0.7);
    LH=2.5E-31*pow(tl[i]/300.0, -1.8);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));		
    kkM[i][33]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    kkM[i][34]=5.21E-35*exp(900.0/tl[i])*MM[i];
	  
    kkM[i][35]=6.0E-34*pow(tl[i]/300.0, -2.4)*MM[i];
	  
    kkM[i][36]=2.8E-36*pow(tl[i]/300.0, -0.9)*MM[i];
    
    RH=8.3E-13*pow(tl[i]/300.0, 2.0);
    LH=5.5E-30;
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][37]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=7.5E-12*pow(tl[i]/300.0, -0.85);
    LH=1.0E-28*pow(tl[i]/300.0, -4.5);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2));
    kkM[i][38]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.1E-12*pow(tl[i]/300.0, 1.3);
    LH=5.9E-33*pow(tl[i]/300.0, -1.4);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][39]=LH*MM[i]/(1+LH*MM[i]/RH)*pow(0.6,ind);
    
    kkM[i][40]=1.8E-27*pow(tl[i]/298.0,-3.85)*MM[i];
	  
    K0=2.4E-14*exp(460.0/tl[i]);
    K2=2.7E-17*exp(2199.0/tl[i]);
    K3=6.5E-34*exp(1355.0/tl[i]);
    kkM[i][41]=K0+K3*MM[i]/(1.0+K3*MM[i]/K2);
    
    RH=3.6E-11*pow(tl[i]/300.0, -0.1);
    LH=7.0E-31*pow(tl[i]/300.0, -2.6);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][42]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=2.8E-11;
    LH=1.8E-30*pow(tl[i]/300.0, -3.0);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][43]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=2.6E-11;
    LH=6.9E-31*pow(tl[i]/300.0, -1.0);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][44]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    LH=1.98E-33*exp(206.0/tl[i]);
	RH=2.26E-14*exp(415.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][45]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	LH=4.82E-31*pow(tl[i]/298.0,-2.17);
	RH=5.31E-11;
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][46]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    LH=6.45E-29*pow(tl[i]/298.0,-3.48)*exp(-490.0/tl[i]);
	RH=4.47E-11*pow(tl[i]/298.0,0.50)*exp(200.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][47]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
    LH=5.61E-30*pow(tl[i]/298.0,-5.19)*exp(-2271.0/tl[i]);
	RH=7.58E-12*pow(tl[i]/298.0,1.59)*exp(1244.0/tl[i]);
	ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	kkM[i][48]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=4.2E-14*pow(tl[i]/300.0, 1.8);
    LH=1.8E-33*pow(tl[i]/300.0, 2.0);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][49]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
    
    RH=1.6E-12;
    LH=3.3E-31*pow(tl[i]/300.0, -4.3);
    ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
    kkM[i][50]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  RH=7.5E-11;
	  LH=5.7E-32*pow(tl[i]/300.0, -1.6);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][51]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][52]=1.1E-30*pow(tl[i]/300.0, -2.0)*MM[i];
	  kkM[i][53]=1.1E-30*pow(tl[i]/300.0, -2.0)*MM[i];
	  kkM[i][54]=4.0E-31*exp(900.0/tl[i])*MM[i];
	  kkM[i][55]=4.0E-31*exp(900.0/tl[i])*MM[i];
	  
	  RH=2.37E-12*exp(523.0/tl[i]);
	  LH=5.8E-30*exp(355.0/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][56]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  LH=6.2E-29*pow(tl[i]/298.0,-1.8);
	  RH=3.5E-10;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][57]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  RH=1.5E-10;
	  LH=5.5E-23*pow(tl[i],-2.0)*exp(-1040/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][58]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][59]=4.42E-10*pow(tl[i]/298.0,0.22)*exp(43.0/tl[i])/THREEBODY;
	  kkM[i][60]=1.9E-10*pow(tl[i]/298.0,0.1)*exp(16.0/tl[i])/THREEBODY;
	  kkM[i][61]=1.1E-10*pow(tl[i]/298.0,0.21)*exp(86.6/tl[i])/THREEBODY;
	  kkM[i][62]=9.49E-12*pow(tl[i]/298.0,1.74)*exp(3872.8/tl[i])/THREEBODY;
	  kkM[i][63]=9.61E-12*exp(-1560.0/tl[i])/THREEBODY;
	  kkM[i][64]=9.49E-12*pow(tl[i]/298.0,1.74)*exp(3872.8/tl[i])/THREEBODY;
	  kkM[i][65]=1.01E-11*pow(tl[i]/298.0,0.69)*exp(1509.0/tl[i])/THREEBODY;
	  kkM[i][66]=2.64E-10*pow(tl[i]/298.0,0.18)*exp(62.5/tl[i])/THREEBODY;
	  kkM[i][67]=7.59E-12*pow(tl[i]/298.0,0.51)*exp(-1318.2/tl[i])/THREEBODY;
	  kkM[i][68]=3.32E-11/THREEBODY;
	  kkM[i][69]=4.23E-12*pow(tl[i]/298.0,2.54)*exp(-3400/tl[i])/THREEBODY;
	  
	  RH=1.39E-10*exp(-1184/tl[i]);
	  LH=1.0E-28;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][70]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][71]=6.68E-11*pow(tl[i]/298.0,0.31)*exp(93.8/tl[i])/THREEBODY;
	  
	  RH=1.39E-10*exp(-1184/tl[i]);
	  LH=1.0E-28;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][72]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  RH=2.0E-10*pow(tl[i]/298.0,0.15);
	  LH=9.0E-31*exp(550.0/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][73]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][74]=1.2E-10/THREEBODY;
	  kkM[i][75]=5.6E-11/THREEBODY;
	  kkM[i][76]=8.3E-12/THREEBODY;
	  kkM[i][77]=8.3E-12/THREEBODY;
	  kkM[i][78]=2.72E-11*pow(tl[i]/298.0,-0.32)*exp(66.1/tl[i])/THREEBODY;
	  kkM[i][79]=1.27E-14*pow(tl[i]/298.0,2.67)*exp(-3447.0/tl[i])/THREEBODY;
	  kkM[i][80]=5.16E-11*pow(tl[i]/298.0,-0.32)/THREEBODY;
	  kkM[i][81]=2.09E-14*pow(tl[i]/298.0,1.9)*exp(-1060/tl[i])/THREEBODY;
	  
	  RH=8.2E-11;
	  LH=8.76E-6*pow(tl[i],-7.03)*exp(-1390/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][82]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  kkM[i][83]=6.5E-11/THREEBODY;
	  kkM[i][84]=1.9E-11/THREEBODY;
	  
	  LH=4.48E-14*pow(tl[i],-5.49)*exp(-1000.0/tl[i]);
	  RH=9.33E-10*pow(tl[i],-0.414)*exp(-33.0/tl[i]);
	  FC=0.31;
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	  kkM[i][85]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(FC,ind);
	  
	  LH=2.01E-18*pow(tl[i],-3.85);
	  RH=3.97E-12*pow(tl[i],0.42);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH),2.0));
	  kkM[i][86]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(0.6,ind);
	  
	  LH=1.932E+3*pow(tl[i],-9.88)*exp(-7544.0/tl[i])+5.109E-11*pow(tl[i],-6.25)*exp(-1433.0/tl[i]);
	  RH=1.031E-10*pow(tl[i],-0.018)*exp(16.74/tl[i]);
	  FC=0.1855*exp(-tl[i]/155.8)+0.8145*exp(-tl[i]/1675.0)+exp(-4531.0/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	  kkM[i][87]=LH*MM[i]/(1.0+LH*MM[i]/RH)*pow(FC,ind);
	  
	  LH=3.714E-8*pow(tl[i],-1.172)*exp(-132.0/tl[i])+1.535E-19*pow(tl[i],2.127)*exp(-547.8/tl[i]);
	  RH=7.976E-4*pow(tl[i],4.096)*exp(625/tl[i]);
	  FC=0.4863*exp(-tl[i]/321.4)+0.5137*exp(-tl[i]/30000.0)+exp(-2804.0/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	  kkM[i][88]=LH/(1.0+LH*MM[i]/RH)*pow(FC,ind);
	  
	  LH=1.092E-14*pow(tl[i],0.996)*exp(-1606.0/tl[i]);
	  RH=5.864E-6*pow(tl[i],5.009)*exp(-949.4/tl[i]);
	  FC=0.8622*exp(-tl[i]/9321.0)+0.1378*exp(-tl[i]/361.8)+exp(-3125.0/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	  kkM[i][89]=LH/(1.0+LH*MM[i]/RH)*pow(FC,ind);
	  
	  LH=2.359E-10*pow(tl[i],-1.234)*exp(19.57/tl[i])+3.141E+3*pow(tl[i],-4.484)*exp(-9188.0/tl[i]);
	  RH=5.876E-14*pow(tl[i],6.721)*exp(1521.0/tl[i]);
	  FC=0.5*exp(-tl[i]/22122.0)+0.5*exp(-tl[i]/174.9)+exp(-3047.0/tl[i]);
	  ind=1.0/(1.0+pow(log10(LH*MM[i]/RH)/(0.75-1.27*log10(FC)),2.0));
	  kkM[i][90]=LH/(1.0+LH*MM[i]/RH)*pow(FC,ind);
	  
	  for (j=1; j<=NKinM; j++) {
		  kkM[i][j] = kkM[i][j]*THREEBODY;
	  }
  }
}











































