/*******************************************************************
 * File: condensed_heat.h
 * PURPOSE: Header file for condensed phase heat capacity functions
 * Units: J/(mol·K)
 *******************************************************************
*/

#ifndef CONDENSED_HEAT_H
#define CONDENSED_HEAT_H

// Function declarations for condensed-phase heat capacities

// Water (H2O) - Liquid and Ice using NIST Shomate Equation
double H2O_liquid_heat_capacity(double T);

// Ammonia (NH3) - Liquid and Solid
double NH3_liquid_heat_capacity(double T);

// Carbon Dioxide (CO2) - Solid (dry ice)
double CO2_solid_heat_capacity(double T);

// Carbon Monoxide (CO) - Liquid and Solid
double CO_condensed_heat_capacity(double T);

// Methane (CH4) - Liquid and Solid  
double CH4_condensed_heat_capacity(double T);

// Hydrogen (H2) - Liquid and Solid
double H2_condensed_heat_capacity(double T);

// Oxygen (O2) - Liquid and Solid
double O2_condensed_heat_capacity(double T);

// Nitrogen (N2) - Liquid and Solid
double N2_condensed_heat_capacity(double T);

// Hydrogen Sulfide (H2S) - Liquid and Solid
double H2S_condensed_heat_capacity(double T);

// Generic function for other condensed species
double generic_condensed_heat_capacity(int species, double T);

#endif // CONDENSED_HEAT_H 