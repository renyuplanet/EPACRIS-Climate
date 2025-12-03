/**
 * ms_gasheat.h
 * 
 * Header file for gas-phase heat capacity functions
 * Units: J/(mol·K)
 */

#ifndef MS_GASHEAT_H
#define MS_GASHEAT_H

// Gas-phase heat capacity function declarations
double AirHeat(double T);   // Type 0 - Air (0.8 N2 + 0.2 O2)
double CO2Heat(double T);   // Type 1 - Carbon dioxide
double HeHeat(double T);    // Type 2 - Helium
double N2Heat(double T);    // Type 3 - Nitrogen
double NH3Heat(double T);   // Type 4 - Ammonia
double CH4Heat(double T);   // Type 5 - Methane
double H2Heat(double T);    // Type 6 - Hydrogen
double O2Heat(double T);    // Type 7 - Oxygen
double COHeat(double T);    // Type 8 - Carbon monoxide
double H2OHeat(double T);   // Type 9 - Water vapor
double H2SHeat(double T);   // Type 10 - Hydrogen sulfide

#endif // MS_GASHEAT_H

