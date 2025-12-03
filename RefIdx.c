/* Function to calculate the refractive index of air and other molecules according to wavelength */
/* condition at 273 K and 1 atm */
/* Unit of wavelength: nm */

#include <math.h>

double AirRefIdx(double ); /* Type 0 */
double CO2RefIdx(double ); /* Type 1 */
double HeRefIdx(double );  /* Type 2 */
double N2RefIdx(double );  /* Type 3 */
double NH3RefIdx(double ); /* Type 4 */
double CH4RefIdx(double ); /* Type 5 */
double H2RefIdx(double );  /* Type 6 */
double O2RefIdx(double );  /* Type 7 */

double AirRefIdx(double l)
{
	double n, x;
	x = l/1000.0;
    if (x<0.23) {x=0.23;}
    if (x>1.69) {x=1.69;}
	n = 1 + 5792105E-8/(238.0185-pow(x,-2.0)) + 167917E-8/(57.362-pow(x,-2.0));
	return n;
}

double CO2RefIdx(double l)
{
	double n, x, w;
	if (l<480.0) {l=480.0;}
	if (l>1800.0) {l=1800.0;}
	x = l*1.0E-3;
/*	n = 1 + 0.012055*(0.579925*pow(x,2)/(166.175*pow(x,2)-1) + 0.12005*pow(x,2)/(79.609*pow(x,2)-1) + 0.0053334*pow(x,2)/(56.3064*pow(x,2)-1) + 0.0043244*pow(x,2)/(46.0196*pow(x,2)-1) + 0.0001218145*pow(x,2)/(0.0584738*pow(x,2)-1)); */
	w = 1.0/x/x;
	n = 1+1.0E-5*(0.154489/(0.0584738-w)+8309.1927/(210.92417-w)+287.64190/(60.122959-w));
	return n;
}

double HeRefIdx(double l)
{
	double n, x;
	x = l/1000.0;
	n = 1 + 0.01470091/(423.98-pow(x,-2.0));
	return n;
}

double N2RefIdx(double l)
{
	double n, x;
	x = l/1000.0;
	if (x<0.145) {x=0.145;}
	if (x>2.06) {x=2.06;}
	if (x > 0.47) {
		n = 1 + 68.5520E-6 + 32431.57E-6*pow(x,2)/(144*pow(x,2)-1);
	} else {
		n = 1 + 1.9662731/(22086.66-pow(x,-2)) + 2.7450825E-2/(133.85688-pow(x,-2));
	}
	return n;
}

double NH3RefIdx(double l)
{
	double n, x;
	x = l/1000.0;
	n = 1.00038;
	return n;
}

double CH4RefIdx(double l)
{
	double n, x;
	x = l/1000.0;
	n = 1.00045;
	return n;
}

double H2RefIdx(double l)
{
	double n, x;
	x = l/1000.0;
	n = 1.00015;
	return n;
}

double O2RefIdx(double l)
{
	double n, x;
	x = l/1000.0;
	n = 1.0003;
	return n;
}



