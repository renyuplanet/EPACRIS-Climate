/**
 * conv_cond_funcs.h
 * 
 * Header file for convection, condensation, and cloud physics functions
 * Contains function declarations for adiabatic calculations, rainout, 
 * cloud physics, and condensation detection
 */

#ifndef CONV_COND_FUNCS_H
#define CONV_COND_FUNCS_H

#include "config.h"

// Convection function declarations
void condensation_and_lapse_rate(int lay, double lapse[], double xxHe, double* cp, double saturation_ratios[]); 
void simulate_rainout(int lay, double* mass_loss_ratio); 
void ms_conv_check(double tempb[],double P[],double lapse[],int isconv[],int* ncl,int* nrl);
void ms_temp_adj(double tempb[],double P[],double lapse[],int isconv[], double cp[], int ncreg[], double pot_temp[]);

// Cloud physics functions
void get_particle_properties(int species_id, double temperature, double *density, double *accommodation_coeff, double *molecular_mass);
void calculate_cloud_properties(double g, double T, double P, double mean_molecular_mass, int condensible_species_id, double Kzz, int layer, double *r0, double *r1, double *r2, double *VP, double *effective_settling_velocity, double *scale_height, double *mass_per_particle, double *n_density);
void exponential_cloud(double gravity, double P[], double **particle_r2);
void cloud_redistribution_none(double gravity, double P[]);

// Global alpha storage functions
void init_alpha_values();
void update_alpha_from_cold_trapping(int layer, int species_index, double alpha_reduction);
double get_global_alpha_value(int layer, int species_index);

// Dynamic condensation detection functions
void initialize_condensibles_mode();
int check_species_condensible(int species_id, double temp, double pressure, double partial_pressure);
void update_condensibles_list(int layer);
void detect_condensibles_atmosphere();
void report_condensibles_changes(int iteration);

// Helper function
const char* get_species_name(int species_id);

#endif // CONV_COND_FUNCS_H

